#include "Database.h"
#include "rsyn/io/parser/lef_def/DEFControlParser.h"


db::Database database;

namespace db {

void Database::initRsynService(){
    rsynService.init();
}

void Database::initPolicy(){
    vector<std ::string> filter_polices;

    boost::split(filter_polices, db::setting.runtimeImprvCRPPolicies, boost::is_any_of(","));

    for(auto policy : filter_polices){
        log() << "polices: " << policy << std::endl;
        policy_set.insert(policy);
    }
        
}
void Database::initIllegalPlacementBoxs(){
     const Rsyn::Session session;
    const Rsyn::PhysicalDesign& physicalDesign =
        static_cast<Rsyn::PhysicalService*>(session.getService("rsyn.physical"))->getPhysicalDesign();
    const DBU libDBU = physicalDesign.getDatabaseUnits(Rsyn::LIBRARY_DBU);


    for(int i = 0; i < lookup_tb.illegal_placement_boxs_strings.size(); i++){
        vector<std ::string> strs;
        boost::split(strs, lookup_tb.illegal_placement_boxs_strings[i],
                    boost::is_any_of(","));

        if(strs.size() != 4){
            // log() << "invalid illegal placement boxs!" << std::endl;
            // std::exit(1);
            continue;
        }

        illegal_placement_boxs.emplace_back(
            std::stof(strs[0])*libDBU,
            std::stof(strs[1])*libDBU,
            std::stof(strs[2])*libDBU,
            std::stof(strs[3])*libDBU
        );                                    

  
    }
    
}//initIllegalPlacementBoxs

void Database::init() {
    // init lookup tabels
    if(db::setting.refinePlacement){
        lookup_tb.run(db::setting.benchmarkName);
        initIllegalPlacementBoxs();
        origin_offset_die = lookup_tb.origin_offset_die;
    }
      
    // return;

    // some debugging metrics
    if(setting.debug){
        ss_legalizer << "inst_critical,area_lx,area_ly,area_hx,area_hy,total_cells,movable_cells" << std::endl;
    }

    // end debugging initialization
    initPolicy();

    total_cells_moved = 0;
    if(db::setting.debug){
        log() << std::endl;
        log() << "################################################################" << std::endl;
        log() << "Start initializing database" << std::endl;
        log() << std::endl;
    }
    
    
    auto dieBound = rsynService.physicalDesign.getPhysicalDie().getBounds();

    if(db::setting.debug){
        log() << "dieBound: " << dieBound << std::endl;
    }
    
    dieRegion = getBoxFromRsynBounds(dieBound);
    if (setting.dbVerbose >= +db::VerboseLevelT::MIDDLE) {
        log() << "Die region (in DBU): " << dieRegion << std::endl;
        log() << std::endl;
    }

    RouteGrid::init();


    // bool debug = false;
    // // movement loop
    // if (debug) {
    //     for (Rsyn::Instance instance : rsynService.module.allInstances()) {
    //         if (instance.getName() != "inst5638") continue;
    //         if (instance.getType() != Rsyn::CELL) continue;
    //         // phCell
    //         Rsyn::Cell cell = instance.asCell();
    //         Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cell);
    //         Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cell);
    //         rsynService.physicalDesign.placeCell(phCell,92500,78660,cell.getOrientation());
    //         //     log() << "after move move inst name: " << instance.getName() 
    //         //         << " ( " << instance.getX() << ", " << instance.getY() << ") "<< std::endl; 
    //         // libPin
    //     }//end movement loop
    // }
    
    NetList::init(rsynService);
    CellList::init(rsynService);
    initSites();
    
    markPinAndObsOccupancy();

    initMTSafeMargin();

    //add index to the vector
    for (auto& net : nets){
        netsToIdx[net.getName()] = net.idx;
        // log() << "net.idx: " << net.idx << ", net: " 
        //       << net.getName() << std::endl;
    } 
    for (auto& cell : cells){
        cellsToIdx[cell.getName()] =cell.idx;
        // log() << "cell.idx: " << cell.idx << ", cell: " 
        //       << cell.getName() << std::endl;
    } 

    for (auto& cell : cells){
        cell.initConnectedCells();
    } 

    initRtrees();
    

    cost_hist.resize(2);

    

    

    log() << "Finish initializing database" << std::endl;
    log() << "MEM: cur=" << utils::mem_use::get_current() << "MB, peak=" << utils::mem_use::get_peak() << "MB"
          << std::endl;
    log() << std::endl;
}

void Database::getNetCells(int net_idx,std::vector<int>& cells_idx){
    auto net = nets[net_idx];
    for (auto pin : net.rsynPins){
        std::string inst_name = pin.getInstanceName();
        int cell_idx = cellsToIdx[inst_name];
        cells_idx.push_back(cell_idx);
    }   
}//end getNetCells

void Database::getIntersectedSites(utils::BoxT<DBU>& box,std::vector<int>& sites){
    int site_lx = (box.lx()-sites_origin)/sites_step;
    int site_hx = (box.hx()-sites_origin)/sites_step;
    for(auto i = site_lx; i < site_hx; i++)
        sites.push_back(i);
}//end getIntersectedSites

void Database::initSites(){

    // DBU sites_origin = 0;
    // DBU sites_step = 0;
    // int sites_num = 0;
    for (Rsyn::PhysicalRow phRow : rsynService.physicalDesign.allPhysicalRows()) {
        sites_origin = phRow.getOrigin().x;
        sites_step = phRow.getStep(X);
        sites_num = phRow.getNumSites(X);
        rows_origin =   phRow.getOrigin().y;
        rows_num = rsynService.physicalDesign.getNumRows();
        rows_step = phRow.getHeight();
        break;
    }  // end for
    log() << "sites_origin: " << sites_origin << std::endl;
    log() << "sites_step: " << sites_step << std::endl;
    log() << "sites_num: " << sites_num << std::endl;
    log() << "rows_origin: " << rows_origin << std::endl;
    log() << "rows_num: " << rows_num << std::endl;
    log() << "rows_step: " << rows_step << std::endl;

}//end initsites
int Database::getSiteDBU(DBU x){
     int step = (x-sites_origin)/sites_step;
    // eq2: lx = step*sites_origin;
    return step;
}

DBU Database::getDBUSite(int step){
    return sites_origin + step*sites_step;
}

int Database::getRowDBU(DBU y){
     int step = (y-rows_origin)/rows_step;
    // eq2: lx = step*sites_origin;
    return step;
}

DBU Database::getDBURow(int step){
    return rows_origin + step*rows_step;
}

void Database::getIntersectedRows(utils::BoxT<DBU>& box,std::vector<int>& rows){
    int row_lx = (box.ly()-rows_origin)/rows_step;
    int row_hx = (box.hy()-rows_origin)/rows_step;
    for(auto i = row_lx; i < row_hx; i++)
        rows.push_back(i);
}//end getIntersectedSites

void Database::initRtrees(){
    auto rsynService = getRsynService();
    cell_rtree.clear();
    empty_rtree.clear();
    row_rtree.clear();

    for (auto& cell : cells){
        // log() << "cell: " << cell.getName() << ", "
        //       << "X: " << cell.rsynInstance.getX() << ", Y: " 
        //       << cell.rsynInstance.getY() << ", bds: "
        //       << cell.rsynInstance.getBounds() << std::endl;

        boostBox box(boostPoint(cell.rsynInstance.getBounds().getLower().x,
                                cell.rsynInstance.getBounds().getLower().y),
                     boostPoint(cell.rsynInstance.getBounds().getUpper().x,
                                cell.rsynInstance.getBounds().getUpper().y));
        cell_rtree.insert({box, cell.idx});
        // cellsToIdx[cell.getName()] =cell.idx; 
        // cell.initConnectedCells();
    }
    // log() << "phyRow..."<< std::endl;
    int emptyIdx = 0;
    int rowIdx = 0;
    for (Rsyn::PhysicalRow row : rsynService.physicalDesign.allPhysicalRows()){
        // log() << row.getName() << ", "<< row.getBounds() << std::endl;

        auto row_box = utils::BoxT<DBU>(row.getBounds().getLower().x,
                                        row.getBounds().getLower().y,
                                        row.getBounds().getUpper().x,
                                        row.getBounds().getUpper().y);
        boostBox row_box_bgt(boostPoint(row.getBounds().getLower().x,
                                        row.getBounds().getLower().y),
                            boostPoint( row.getBounds().getUpper().x,
                                        row.getBounds().getUpper().y));
    
        row_rtree.insert({row_box_bgt, rowIdx});
        rowsRsyn.push_back(row);
        rowIdx++;
        
        std::vector<db::Cell> cells;
        getCellsInBox(row_box,cells,10);

        std::sort(cells.begin(), cells.end(), [&](const db::Cell& lhs, const db::Cell& rhs) {
            return lhs.rsynInstance.getBounds().getLower().x <
                   rhs.rsynInstance.getBounds().getLower().x;
        });
        
        
        auto start = row.getBounds().getLower().x;
        auto stop = row.getBounds().getUpper().x;
        for (int i = 0; i < cells.size(); i++) {
            int emptySpace = cells[i].rsynInstance.getBounds().getLower().x - start;
            if (emptySpace > 0) {
                
                boostBox box(boostPoint(start,
                                cells[i].rsynInstance.getBounds().getLower().y),
                            boostPoint(cells[i].rsynInstance.getBounds().getLower().x,
                                cells[i].rsynInstance.getBounds().getUpper().y));
                empty_rtree.insert({box, emptyIdx}); 
                emptyIdx ++;
            }
            start = cells[i].rsynInstance.getBounds().getUpper().x;
        }//end for 
        int emptySpace = stop - start;
        if (emptySpace > 0) {
            
            boostBox box(boostPoint(start,
                            row.getBounds().getLower().y),
                        boostPoint(stop,
                            row.getBounds().getUpper().y));
            empty_rtree.insert({box, emptyIdx}); 
            emptyIdx ++;
        }//end if

        // std::vector<utils::BoxT<DBU>> emptys;
        // getEmptySpacesInBox(row_box,emptys,10);
    }//end row for loop
}//end initRtrees



void Database::getRowsInBox(utils::BoxT<DBU>& box,std::vector<utils::BoxT<DBU>>& rows, DBU eps){
    // log() << "empty space query..." << std::endl;
    boostBox rtreeQueryBox(
        boostPoint(box.lx()+eps,
                   box.ly()+eps),
        boostPoint(box.hx()-eps,
                   box.hy()-eps));
     vector<std::pair<boostBox, int>> queryResults;
        
    row_rtree.query(bgi::intersects(rtreeQueryBox), std::back_inserter(queryResults));
    
    for (const auto& queryResult : queryResults) {
        const auto& b = queryResult.first;
        
        
        auto box = utils::BoxT<DBU>(bg::get<bg::min_corner, 0>(b),
                                    bg::get<bg::min_corner, 1>(b),
                                    bg::get<bg::max_corner, 0>(b),
                                    bg::get<bg::max_corner, 1>(b));
        
        rows.push_back(box);
    }//end for
}//end getEmptySpacesInBox

void Database::getCellsInBox(utils::BoxT<DBU>& box,std::vector<db::Cell>& cells,DBU eps){
     boostBox rtreeQueryBox(
        boostPoint(box.lx()+eps,
                   box.ly()+eps),
        boostPoint(box.hx()-eps,
                   box.hy()-eps));
    vector<std::pair<boostBox, int>> queryResults;
    cell_rtree.query(bgi::intersects(rtreeQueryBox), std::back_inserter(queryResults));
    for (const auto& queryResult : queryResults) {
        auto cell_idx = queryResult.second;
        auto cell = database.cells[cell_idx];
        cells.push_back(cell);
    }//end for 

}//end getCellsInBox

void Database::markPinAndObsOccupancy() {
    bool debug = false;

    if(db::setting.debug){    
        if (db::setting.dbVerbose >= +db::VerboseLevelT::MIDDLE) {
            log() << "Mark pin & obs occupancy on RouteGrid ..." << std::endl;
        }
    }
    vector<std::pair<BoxOnLayer, int>> fixedMetalVec;

    // STEP 1: get fixed objects
    // Mark pins associated with nets
    for (const auto& net : nets) {
        for (const auto& accessBoxes : net.pinAccessBoxes) {
            for (const auto& box : accessBoxes) {
                fixedMetalVec.emplace_back(box, net.idx);
            }
        }
    }
    // Mark dangling pins
    // minor TODO: port?
    const Rsyn::Session session;
    const Rsyn::PhysicalDesign& physicalDesign =
        static_cast<Rsyn::PhysicalService*>(session.getService("rsyn.physical"))->getPhysicalDesign();
    const DBU libDBU = physicalDesign.getDatabaseUnits(Rsyn::LIBRARY_DBU);
    unsigned numUnusedPins = 0;
    unsigned numObs = 0;
    unsigned numSNetObs = 0;

    // write instances in csv file
    // std::stringstream ss;
    // ss << "inst_name,x,y,width,height" << std::endl;


    // // movement loop
    // if (debug) {
    //     for (Rsyn::Instance instance : rsynService.module.allInstances()) {
    //         if (instance.getType() != Rsyn::CELL) continue;
    //         // phCell
    //         Rsyn::Cell cell = instance.asCell();
    //         Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cell);
    //         Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cell);
    //         rsynService.physicalDesign.placeCell(phCell,10,10,cell.getOrientation());
    //         //     log() << "after move move inst name: " << instance.getName() 
    //         //         << " ( " << instance.getX() << ", " << instance.getY() << ") "<< std::endl; 
    //         // libPin
    //     }//end movement loop
    // }
    

    for (Rsyn::Instance instance : rsynService.module.allInstances()) {

        if (instance.getType() != Rsyn::CELL) continue;
        // phCell
        Rsyn::Cell cell = instance.asCell();
        Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cell);
        Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cell);
        const DBUxy origin(static_cast<DBU>(std::round(phLibCell.getMacro()->originX() * libDBU)),
                           static_cast<DBU>(std::round(phLibCell.getMacro()->originY() * libDBU)));


        // if (debug) {
        //     log() << "before move inst name: " << instance.getName() 
        //           << " ( " << instance.getX() << ", " << instance.getY() << ") "<< std::endl; 

        //     // add write instances in csv file
        //     log() << "Writing guides to file..." << std::endl;
        //     ss << instance.getName() << ", " << instance.getX() << ", "
        //        << instance.getY() << ", " << phLibCell.getWidth() << ", " 
        //        << phLibCell.getHeight() << std::endl;
        //     //end write instances in csv file
        
        // }//end debug

        
        for (Rsyn::Pin pin : instance.allPins(false)) {
            if (!pin.getNet()) {  // no associated net
                Rsyn::PhysicalLibraryPin phLibPin = rsynService.physicalDesign.getPhysicalLibraryPin(pin);
                vector<BoxOnLayer> accessBoxes;
                Net::getPinAccessBoxes(phLibPin, phCell, accessBoxes, origin);
                for (const auto& box : accessBoxes) {

                    fixedMetalVec.emplace_back(box, OBS_NET_IDX);
                }
                ++numUnusedPins;
            }
        }
        // libObs
        DBUxy displacement = phCell.getPosition() + origin;
        auto transform = phCell.getTransform();
        for (const Rsyn::PhysicalObstacle& phObs : phLibCell.allObstacles()) {
            if (phObs.getLayer().getType() != Rsyn::PhysicalLayerType::ROUTING) continue;
            const int layerIdx = phObs.getLayer().getRelativeIndex();
            for (auto bounds : phObs.allBounds()) {
                bounds.translate(displacement);
                bounds = transform.apply(bounds);
                const BoxOnLayer box(layerIdx, getBoxFromRsynBounds(bounds));
                fixedMetalVec.emplace_back(box, OBS_NET_IDX);
                ++numObs;
            }
        }
    }


    // write instances to csv file
    // std::string file_name_csv = "instances.csv";
    // std::ofstream fout(file_name_csv);
    // fout << ss.str();
    // fout.close();

    // Mark special nets
    // bool debug = true;
    for (Rsyn::PhysicalSpecialNet specialNet : rsynService.physicalDesign.allPhysicalSpecialNets()) {
        for (const DefWireDscp& wire : specialNet.getNet().clsWires) {
            for (const DefWireSegmentDscp& segment : wire.clsWireSegments) {
                int layerIdx =
                    rsynService.physicalDesign.getPhysicalLayerByName(segment.clsLayerName).getRelativeIndex();
                const DBU width = segment.clsRoutedWidth;
                DBUxy pos;
                DBU ext = 0;
                for (unsigned i = 0; i != segment.clsRoutingPoints.size(); ++i) {
                    const DefRoutingPointDscp& pt = segment.clsRoutingPoints[i];
                    const DBUxy& nextPos = pt.clsPos;
                    const DBU nextExt = pt.clsHasExtension ? pt.clsExtension : 0;
                    if (i >= 1) {
                        for (unsigned dim = 0; dim != 2; ++dim) {
                            if (pos[dim] == nextPos[dim]) continue;
                            const DBU l = pos[dim] < nextPos[dim] ? pos[dim] - ext : nextPos[dim] - nextExt;
                            const DBU h = pos[dim] < nextPos[dim] ? nextPos[dim] + nextExt : pos[dim] + ext;
                            BoxOnLayer box(layerIdx);
                            box[dim].Set(l, h);
                            box[1 - dim].Set(pos[1 - dim] - width / 2, pos[1 - dim] + width / 2);
                            if(debug){
                                // if(box.layerIdx==3){
                                    log() << "pt.ClsPos.x : " << pt.clsPos.x/2000.0 
                                          << ", y: " << pt.clsPos.y/2000.0 << std::endl;
                                    log() << "specialNet: " << box << std::endl;
                                // }
                            }
                            
                            fixedMetalVec.emplace_back(box, OBS_NET_IDX);
                            ++numSNetObs;
                            break;
                        }
                    }
                    pos = nextPos;
                    ext = nextExt;
                    if (!pt.clsHasVia) continue;
                    const Rsyn::PhysicalVia& via = rsynService.physicalDesign.getPhysicalViaByName(pt.clsViaName);
                    const int botLayerIdx = via.getBottomLayer().getRelativeIndex();
                    for (const Rsyn::PhysicalViaGeometry& geo : via.allBottomGeometries()) {
                        Bounds bounds = geo.getBounds();
                        bounds.translate(pos);
                        const BoxOnLayer box(botLayerIdx, getBoxFromRsynBounds(bounds));
                        if(debug){
                            if(box.layerIdx==3){
                                log() << "specialNet: " << box << std::endl;
                            }
                        }
                        fixedMetalVec.emplace_back(box, OBS_NET_IDX);
                        ++numSNetObs;
                    }
                    const int topLayerIdx = via.getTopLayer().getRelativeIndex();
                    for (const Rsyn::PhysicalViaGeometry& geo : via.allTopGeometries()) {
                        Bounds bounds = geo.getBounds();
                        bounds.translate(pos);
                        const BoxOnLayer box(topLayerIdx, getBoxFromRsynBounds(bounds));
                        if(debug){
                            if(box.layerIdx==3){
                                log() << "specialNet: " 
                                      << specialNet.getNet().clsName 
                                      << ", box: " << box << std::endl;
                            }
                        }
                        fixedMetalVec.emplace_back(box, OBS_NET_IDX);
                        ++numSNetObs;
                    }
                    if (via.hasViaRule()) {
                        const utils::PointT<int> numRowCol =
                            via.hasRowCol() ? utils::PointT<int>(via.getNumCols(), via.getNumRows())
                                            : utils::PointT<int>(1, 1);
                        BoxOnLayer botBox(botLayerIdx);
                        BoxOnLayer topBox(topLayerIdx);
                        for (unsigned dimIdx = 0; dimIdx != 2; ++dimIdx) {
                            const Dimension dim = static_cast<Dimension>(dimIdx);
                            const DBU origin = via.hasOrigin() ? pos[dim] + via.getOrigin(dim) : pos[dim];
                            const DBU botOff =
                                via.hasOffset() ? origin + via.getOffset(Rsyn::BOTTOM_VIA_LEVEL, dim) : origin;
                            const DBU topOff =
                                via.hasOffset() ? origin + via.getOffset(Rsyn::TOP_VIA_LEVEL, dim) : origin;
                            const DBU length =
                                (via.getCutSize(dim) * numRowCol[dim] + via.getSpacing(dim) * (numRowCol[dim] - 1)) / 2;
                            const DBU botEnc = length + via.getEnclosure(Rsyn::BOTTOM_VIA_LEVEL, dim);
                            const DBU topEnc = length + via.getEnclosure(Rsyn::TOP_VIA_LEVEL, dim);
                            botBox[dim].Set(botOff - botEnc, botOff + botEnc);
                            topBox[dim].Set(topOff - topEnc, topOff + topEnc);
                        }
                        if(debug){
                            if(botBox.layerIdx==3){
                                log() << "specialNet botBox: " << botBox << std::endl;
                            }
                        }
                        fixedMetalVec.emplace_back(botBox, OBS_NET_IDX);
                        if(debug){
                            if(topBox.layerIdx==3){
                                log() << "specialNet botBox: " << topBox << std::endl;
                            }
                        }
                        fixedMetalVec.emplace_back(topBox, OBS_NET_IDX);
                        numSNetObs += 2;
                    }
                    if (layerIdx == botLayerIdx)
                        layerIdx = topLayerIdx;
                    else if (layerIdx == topLayerIdx)
                        layerIdx = botLayerIdx;
                    else {
                        log() << "Error: Special net " << specialNet.getNet().clsName << " via " << pt.clsViaName
                              << " on wrong layer " << layerIdx << std::endl;
                        break;
                    }
                }
            }
        }
    }
    // Stat
    vector<int> layerNumFixedObjects(getLayerNum(), 0);
    for (const auto& fixedMetal : fixedMetalVec) {
        layerNumFixedObjects[fixedMetal.first.layerIdx]++;
    }
    // Print
    if(db::setting.debug){
        if (setting.dbVerbose >= +db::VerboseLevelT::MIDDLE) {
            log() << "The number of unused pins is " << numUnusedPins << std::endl;
            log() << "The number of OBS is " << numObs << std::endl;
            log() << "The number of special net OBS is " << numSNetObs << std::endl;
            log() << "The number of fixed objects on each layers:" << std::endl;
            for (unsigned i = 0; i < getLayerNum(); i++) {
                if (layerNumFixedObjects[i] > 0) log() << getLayer(i).name << ": " << layerNumFixedObjects[i] << std::endl;
            }
        }
        log() << std::endl;
    }

    // STEP 2: mark
    if(db::setting.debug){
        if (setting.dbVerbose >= +db::VerboseLevelT::MIDDLE) {
            printlog("mark fixed metal rtrees...");
        }
    }

    markFixedMetalBatch(fixedMetalVec, 0, fixedMetalVec.size());

    // logFixedMetals(fixedMetalVec);
}

void Database::initMTSafeMargin() {
    for (auto& layer : layers) {
        layer.mtSafeMargin = max({layer.minAreaMargin, layer.confLutMargin, layer.fixedMetalQueryMargin});
        if(db::setting.debug){
            if (db::setting.dbVerbose >= +db::VerboseLevelT::MIDDLE) {
            printlog(layer.name,
                     "mtSafeMargin = max {",
                     layer.minAreaMargin,
                     layer.confLutMargin,
                     layer.fixedMetalQueryMargin,
                     "} =",
                    layer.mtSafeMargin);
            }
        }
        
    }
}

void Database::getGridPinAccessBoxes(const Net& net, vector<vector<db::GridBoxOnLayer>>& gridPinAccessBoxes) const {
    gridPinAccessBoxes.resize(net.numOfPins());
    for (unsigned pinIdx = 0; pinIdx != net.numOfPins(); ++pinIdx) {
        vector<vector<db::GridBoxOnLayer>> pins(getLayerNum());
        for (const db::BoxOnLayer& pinAccessBox : net.pinAccessBoxes[pinIdx]) {
            int dir = getLayerDir(pinAccessBox.layerIdx);
            DBU pitch = getLayer(pinAccessBox.layerIdx).pitch;
            // pinForbidRegion
            auto pinForbidRegion = getMetalRectForbidRegion(pinAccessBox, AggrParaRunSpace::DEFAULT);
            const db::GridBoxOnLayer& gridPinForbidRegion = rangeSearch(pinForbidRegion);
            if (isValid(gridPinForbidRegion)) {
                pins[pinAccessBox.layerIdx].push_back(gridPinForbidRegion);
            }
            // One-pitch extension
            auto pinExtension = pinAccessBox;
            for (int d = 0; d < 2; ++d) {
                pinExtension[d].low -= pitch;
                pinExtension[d].high += pitch;
            }
            const db::GridBoxOnLayer& gridPinExtension = rangeSearch(pinExtension);
            for (int trackIdx = gridPinExtension.trackRange.low; trackIdx <= gridPinExtension.trackRange.high;
                 ++trackIdx) {
                for (int cpIdx = gridPinExtension.crossPointRange.low; cpIdx <= gridPinExtension.crossPointRange.high;
                     ++cpIdx) {
                    db::GridPoint pt(pinAccessBox.layerIdx, trackIdx, cpIdx);
                    if (!gridPinForbidRegion.includePoint(pt) && Dist(pinAccessBox, getLoc(pt)) <= pitch) {
                        pins[pinAccessBox.layerIdx].emplace_back(pinAccessBox.layerIdx,
                                                                 utils::IntervalT<int>{trackIdx, trackIdx},
                                                                 utils::IntervalT<int>{cpIdx, cpIdx});
                    }
                }
            }
        }

        // assign a relatively far grid access box if none (rarely happen)
        unsigned numBoxes = 0;
        for (const vector<db::GridBoxOnLayer>& pin : pins) {
            numBoxes += pin.size();
        }
        if (!numBoxes) {
            for (const db::BoxOnLayer& pinAccessBox : net.pinAccessBoxes[pinIdx]) {
                db::GridBoxOnLayer gridBox = rangeSearch(pinAccessBox);
                if (gridBox.trackRange.low > gridBox.trackRange.high) {
                    if (gridBox.trackRange.low == 0) {
                        gridBox.trackRange.high = 0;
                    } else {
                        gridBox.trackRange.low = gridBox.trackRange.high;
                    }
                }
                if (gridBox.crossPointRange.low > gridBox.crossPointRange.high) {
                    if (gridBox.crossPointRange.low == 0) {
                        gridBox.crossPointRange.high = 0;
                    } else {
                        gridBox.crossPointRange.low = gridBox.crossPointRange.high;
                    }
                }
                pins[pinAccessBox.layerIdx].push_back(gridBox);
            }
        }

        // slice
        gridPinAccessBoxes[pinIdx].clear();
        for (vector<db::GridBoxOnLayer>& pin : pins) {
            if (!pin.empty()) {
                db::GridBoxOnLayer::sliceGridPolygons(pin);
                for (const db::GridBoxOnLayer& box : pin) {
                    if (isValid(box)) {
                        gridPinAccessBoxes[pinIdx].push_back(box);
                    }
                }
            }
        }
        if (gridPinAccessBoxes[pinIdx].empty()) {
            log() << "Error: Net " << net.getName() << " Pin " << pinIdx << " has empty grid pin access boxes\n";
        }
    }
}

void Database::logCellLocations(int iter){
    bool debug = false;
    // 
    std::string file_name = db::setting.directory +  db::setting.benchmarkName+ ".cell."+std::to_string(iter) + ".csv";
    if(debug)log() << "logCellLocations file name: " << file_name << std::endl;
    std::ofstream file(file_name);
    std::stringstream stream;
    stream << "cell_name,x,y,w,h" << std::endl;

    for(auto cell : database.cells){
        auto box = cell.getCellBox();
        stream << cell.getName() 
               << "," << box.lx()
               << "," << box.ly()
               << "," << box.width()
               << "," << box.height()
                << std::endl;
    }
    file << stream.str();
    file.close();
}//end logCellLocations

// void Database::logFixedMetals(vector<std::pair<BoxOnLayer, int>>& fixedMetalVec){
void Database::logFixedMetals(int iter){
    std::string file_name = db::setting.directory + \
     db::setting.benchmarkName+ ".fixedMetals."+std::to_string(iter)+".csv";
    
    std::ofstream file(file_name);
    std::stringstream ss;
    ss << "l,xl,yl,xh,yh" << std::endl;



    for (int l_idx = 0; l_idx < database.getLayerNum(); l_idx++){
        auto rtree = database.getFixedMetals(l_idx);
        for(auto tmp : rtree){
            // ss << l_idx 
            //    << ", " << tmp.first 
            //    << ", " << tmp.second << std::endl;
            auto b = tmp.first;
            auto box = utils::BoxT<DBU>(bg::get<bg::min_corner, 0>(b),
                                    bg::get<bg::min_corner, 1>(b),
                                    bg::get<bg::max_corner, 0>(b),
                                    bg::get<bg::max_corner, 1>(b));
            // log() << tmp.second << std::endl;

            ss << std::to_string(l_idx) 
               << ", " << std::to_string(box.lx())
               << ", " << std::to_string(box.ly())
               << ", " << std::to_string(box.hx())
               << ", " << std::to_string(box.hy()) << std::endl;
               
        }
    }

    // for(auto fixedMetal : fixedMetalVec){
    //     auto box = fixedMetal.first;
    //     stream << box.layerIdx
    //            << "," << box.lx()
    //            << "," << box.ly()
    //            << "," << box.hx()
    //            << "," << box.hy()
    //            << std::endl;
    // }
    file << ss.str();
    file.close();
}//end logFixedMetals




void Database::logDie(){
    std::string file_name = db::setting.directory +  db::setting.benchmarkName+ ".die.csv";
    std::ofstream file(file_name);
    std::stringstream stream;
    stream << "bench_name,die_xl,die_yl,die_xh,die_yh" << std::endl;
    auto dieBound = rsynService.physicalDesign.getPhysicalDie().getBounds();
    
    stream << db::setting.benchmarkName
            << "," << dieBound.getLower().x
            << "," << dieBound.getLower().y
            << "," << dieBound.getUpper().x
            << "," << dieBound.getUpper().y
            << std::endl;
    
    file << stream.str();
    file.close();
}

void Database::logPatternRoute(int iter){
    std::string file_name = db::setting.directory +  db::setting.benchmarkName
        + ".patternroute."+std::to_string(iter)+".csv";
    std::ofstream file(file_name);
    std::string head;
    head = "net_name,l,xl,yl,xh,yh,type,cost\n";
    head = head + patternRoute_stream;
    file << head;
    file.close();
    // clean the pattern stream in each iteration
    patternRoute_stream = "";
}

void Database::logAstar(int iter){
    std::string file_name = db::setting.directory +  db::setting.benchmarkName
        + ".astar."+std::to_string(iter)+".csv";
    std::ofstream file(file_name);
    std::string head;
    head = "idx,net_name,l,xl,yl,xh,yh,eBackward,eForward,eDown,eUp,eBackwardCost,eForwardCost,eDownCost,eUpCost\n";
    head = head + astar_stream;
    file << head;
    file.close();
    // clean the pattern stream in each iteration
    astar_stream = "";
}

void Database::logAstarCoraseGrid(int iter){
    std::string file_name = db::setting.directory +  db::setting.benchmarkName
        + ".astar.coarsegrid."+std::to_string(iter)+".csv";
    std::ofstream file(file_name);
    std::string head;
    head = "idx,net_name,l,xl,yl,xh,yh,eBackward,eForward,eDown,eUp,eBackwardCost,eForwardCost,eDownCost,eUpCost\n";
    head = head + astar_stream_coarse;
    file << head;
    file.close();
    // clean the pattern stream in each iteration
    astar_stream_coarse = "";
}

void Database::logLayers(){
    std::string file_name = db::setting.directory +  db::setting.benchmarkName
        + ".layers.csv";
    std::ofstream file(file_name);       
    std::stringstream stream;
    stream << "name,idx,direction,width,minArea" << std::endl;

    for(auto layer : layers){
        stream << layer.name
               << "," << layer.idx
               << "," << layer.direction
               << "," << layer.width
               << "," << layer.minArea
               << std::endl;
    }
    file << stream.str();
    file.close();

}//end logLayers

void Database::writeDEF(const std::string& filename){
    DEFControlParser defParser;
    DefDscp def;
    def.clsDesignName = rsynService.design.getName();
    def.clsDatabaseUnits = rsynService.physicalDesign.getDatabaseUnits(Rsyn::DESIGN_DBU);
    def.clsHasDatabaseUnits = true;
    def.clsDieBounds = rsynService.physicalDesign.getPhysicalDie().getBounds();
    def.clsHasDieBounds = true;
    def.clsRows.reserve(rsynService.physicalDesign.getNumRows());
    for (Rsyn::PhysicalRow phRow : rsynService.physicalDesign.allPhysicalRows()) {
        def.clsRows.push_back(DefRowDscp());
        DefRowDscp& defRow = def.clsRows.back();
        defRow.clsName = phRow.getName();
        defRow.clsSite = phRow.getSiteName();
        defRow.clsOrigin = phRow.getOrigin();
        defRow.clsStepX = phRow.getStep(X);
        defRow.clsStepY = 0;  // phRow.getStep(Y);
        defRow.clsNumX = phRow.getNumSites(X);
        defRow.clsNumY = phRow.getNumSites(Y);
        defRow.clsOrientation = Rsyn::getPhysicalOrientation(phRow.getSiteOrientation());
    }  // end for

    // write def tracks
    def.clsTracks.reserve(rsynService.physicalDesign.getNumPhysicalTracks());
    for (const Rsyn::PhysicalTracks& phTrack : rsynService.physicalDesign.allPhysicalTracks()) {
        def.clsTracks.emplace_back();
        DefTrackDscp& defTrack = def.clsTracks.back();
        defTrack.clsDirection = Rsyn::getPhysicalTrackDirectionDEF(phTrack.getDirection());
        defTrack.clsLocation = phTrack.getLocation();
        defTrack.clsSpace = phTrack.getSpace();
        int numLayers = phTrack.getNumberOfLayers();
        defTrack.clsLayers.reserve(numLayers);
        for (Rsyn::PhysicalLayer phLayer : phTrack.allLayers()) defTrack.clsLayers.push_back(phLayer.getName());
        defTrack.clsNumTracks = phTrack.getNumberOfTracks();
    }  // end for

    // write def vias
    unsigned numVias = 0;
    for (const Rsyn::PhysicalVia& phVia : rsynService.physicalDesign.allPhysicalVias()) {
        if (phVia.isViaDesign()) ++numVias;
    }  // end for
    def.clsVias.reserve(numVias);
    for (const Rsyn::PhysicalVia& phVia : rsynService.physicalDesign.allPhysicalVias()) {
        if (!phVia.isViaDesign()) continue;
        def.clsVias.emplace_back();
        DefViaDscp& defVia = def.clsVias.back();
        defVia.clsName = phVia.getName();
        defVia.clsBottomLayer = phVia.getBottomLayer().getName();
        defVia.clsCutLayer = phVia.getCutLayer().getName();
        defVia.clsTopLayer = phVia.getTopLayer().getName();
        defVia.clsGeometries.emplace(defVia.clsBottomLayer, std::deque<DefViaGeometryDscp>());
        for (const Rsyn::PhysicalViaGeometry& phGeo : phVia.allBottomGeometries()) {
            DefViaGeometryDscp defGeo;
            defGeo.clsBounds = phGeo.getBounds();
            defVia.clsGeometries[defVia.clsBottomLayer].push_back(defGeo);
        }  // end for
        defVia.clsGeometries.emplace(defVia.clsCutLayer, std::deque<DefViaGeometryDscp>());
        for (const Rsyn::PhysicalViaGeometry& phGeo : phVia.allCutGeometries()) {
            DefViaGeometryDscp defGeo;
            defGeo.clsBounds = phGeo.getBounds();
            defVia.clsGeometries[defVia.clsCutLayer].push_back(defGeo);
        }  // end for
        defVia.clsGeometries.emplace(defVia.clsTopLayer, std::deque<DefViaGeometryDscp>());
        for (const Rsyn::PhysicalViaGeometry& phGeo : phVia.allTopGeometries()) {
            DefViaGeometryDscp defGeo;
            defGeo.clsBounds = phGeo.getBounds();
            defVia.clsGeometries[defVia.clsTopLayer].push_back(defGeo);
        }  // end for
        if (phVia.hasViaRule()) {
            defVia.clsHasViaRule = true;
            defVia.clsViaRuleName = phVia.getViaRule().getName();
            defVia.clsXCutSize = phVia.getCutSize(X);
            defVia.clsYCutSize = phVia.getCutSize(Y);
            defVia.clsXCutSpacing = phVia.getSpacing(X);
            defVia.clsYCutSpacing = phVia.getSpacing(Y);
            defVia.clsXBottomEnclosure = phVia.getEnclosure(Rsyn::BOTTOM_VIA_LEVEL, X);
            defVia.clsYBottomEnclosure = phVia.getEnclosure(Rsyn::BOTTOM_VIA_LEVEL, Y);
            defVia.clsXTopEnclosure = phVia.getEnclosure(Rsyn::TOP_VIA_LEVEL, X);
            defVia.clsYTopEnclosure = phVia.getEnclosure(Rsyn::TOP_VIA_LEVEL, Y);
        }
        if (phVia.hasRowCol()) {
            defVia.clsHasRowCol = true;
            defVia.clsNumCutRows = phVia.getNumRows();
            defVia.clsNumCutCols = phVia.getNumCols();
        }
        if (phVia.hasOrigin()) {
            defVia.clsHasOrigin = true;
            defVia.clsXOffsetOrigin = phVia.getOrigin(X);
            defVia.clsYOffsetOrigin = phVia.getOrigin(Y);
        }
        if (phVia.hasOffset()) {
            defVia.clsHasOffset = true;
            defVia.clsXBottomOffset = phVia.getOffset(Rsyn::BOTTOM_VIA_LEVEL, X);
            defVia.clsYBottomOffset = phVia.getOffset(Rsyn::BOTTOM_VIA_LEVEL, Y);
            defVia.clsXTopOffset = phVia.getOffset(Rsyn::TOP_VIA_LEVEL, X);
            defVia.clsYTopOffset = phVia.getOffset(Rsyn::TOP_VIA_LEVEL, Y);
        }
        if (phVia.hasPattern()) {
            defVia.clsHasPattern = true;
            defVia.clsPattern = phVia.getPattern();
        }
    }  // end for

    def.clsComps.reserve(rsynService.design.getNumInstances(Rsyn::CELL));
    for (Rsyn::Instance instance : rsynService.module.allInstances()) {
        if (instance.getType() != Rsyn::CELL) continue;

        Rsyn::Cell cell = instance.asCell();  // minor TODO: hack, assuming that the instance is a cell
        Rsyn::PhysicalCell ph = rsynService.physicalDesign.getPhysicalCell(cell);
        def.clsComps.push_back(DefComponentDscp());
        DefComponentDscp& defComp = def.clsComps.back();
        defComp.clsName = cell.getName();
        defComp.clsMacroName = cell.getLibraryCellName();
        defComp.clsPos = ph.getPosition();
        defComp.clsIsFixed = instance.isFixed();
        defComp.clsOrientation = Rsyn::getPhysicalOrientation(ph.getOrientation());
        defComp.clsIsPlaced = ph.isPlaced();
    }  // end for

    // write def special nets
    def.clsSpecialNets.reserve(rsynService.physicalDesign.getNumPhysicalSpecialNets());
    for (Rsyn::PhysicalSpecialNet phSpecialNet : rsynService.physicalDesign.allPhysicalSpecialNets()) {
        def.clsSpecialNets.push_back(phSpecialNet.getNet());
    }
    //----------------------------------------------------------------
    int numNets = rsynService.design.getNumNets();
    int i = 0;
    def.clsNets.reserve(numNets);
    for (Rsyn::Net net : rsynService.module.allNets()) {
        switch(net.getUse()) {
            case Rsyn::POWER:
                continue;
            case Rsyn::GROUND:
                continue;
            default:
                break;
        }
        def.clsNets.push_back(DefNetDscp());
        DefNetDscp& defNet = def.clsNets.back();
        defNet.clsName = net.getName();
        defNet.clsConnections.reserve(net.getNumPins());
        for (Rsyn::Pin pin : net.allPins()) {
            if (!pin.isPort()) continue;
            defNet.clsConnections.push_back(DefNetConnection());
            DefNetConnection& netConnection = defNet.clsConnections.back();
            netConnection.clsComponentName = "PIN";
            netConnection.clsPinName = pin.getInstanceName();
        }  // end for
        for (Rsyn::Pin pin : net.allPins()) {
            if (pin.isPort()) continue;
            defNet.clsConnections.push_back(DefNetConnection());
            DefNetConnection& netConnection = defNet.clsConnections.back();
            netConnection.clsComponentName = pin.getInstanceName();
            netConnection.clsPinName = pin.getName();
        }  // end for


        auto wire_descps = net.getWires();
        for(auto tmp : wire_descps )
            defNet.clsWires.push_back(tmp);
        // --- Net wire write ----------
        // for (const DefWireDscp& wire : wire_descps) {
        //     std::cout << "wire.clsWireSegments size: " << wire.clsWireSegments.size() << std::endl;
        //     for (const DefWireSegmentDscp& segment : wire.clsWireSegments) {
        //         int layerIdx =
        //             rsynService.physicalDesign.getPhysicalLayerByName(segment.clsLayerName).getRelativeIndex();
        //         const DBU width = segment.clsRoutedWidth;
        //         DBUxy pos;
        //         DBU ext = 0;
        //         std::cout << "segment.clsRoutingPoints: " << segment.clsRoutingPoints.size() << std::endl;
        //         for (unsigned i = 0; i != segment.clsRoutingPoints.size(); ++i) {
        //             const DefRoutingPointDscp& pt = segment.clsRoutingPoints[i];
        //             const DBUxy& nextPos = pt.clsPos;
        //             const DBU nextExt = pt.clsHasExtension ? pt.clsExtension : 0;
        //             if (i >= 1) {
        //                 for (unsigned dim = 0; dim != 2; ++dim) {
        //                     if (pos[dim] == nextPos[dim]) continue;
        //                     const DBU l = pos[dim] < nextPos[dim] ? pos[dim] - ext : nextPos[dim] - nextExt;
        //                     const DBU h = pos[dim] < nextPos[dim] ? nextPos[dim] + nextExt : pos[dim] + ext;
        //                     BoxOnLayer box(layerIdx);
        //                     box[dim].Set(l, h);
        //                     box[1 - dim].Set(pos[1 - dim] - width / 2, pos[1 - dim] + width / 2);
        //                     // fixedMetalVec.emplace_back(box, OBS_NET_IDX);
        //                     // initRouteDR.push_back(box);
        //                     log() << "wire_box: " << box << std::endl;
        //                     // ++numSNetObs;
        //                     break;
        //                 }
        //             }
        //             pos = nextPos;
        //             ext = nextExt;
        //             log() << "pt: " << nextPos << std::endl;
        //             log() << "has extension: " <<pt.clsHasExtension << std::endl;
        //             log() << "pt.clsHasVia: " << pt.clsHasVia << std::endl;
        //             log() << "pt.clsHasRectangle: " << pt.clsHasRectangle << std::endl;
        //             if(pt.clsHasRectangle){
        //                 log() << "rect: " << pt.clsRect << std::endl;
        //             }
        //             if (!pt.clsHasVia) continue;
        //             const Rsyn::PhysicalVia& via = rsynService.physicalDesign.getPhysicalViaByName(pt.clsViaName);
        //             const int botLayerIdx = via.getBottomLayer().getRelativeIndex();
        //             for (const Rsyn::PhysicalViaGeometry& geo : via.allBottomGeometries()) {
        //                 Bounds bounds = geo.getBounds();
        //                 bounds.translate(pos);
        //                 const BoxOnLayer box(botLayerIdx, getBoxFromRsynBounds(bounds));
        //                 log() << "via_box: " << box << std::endl;
        //                 // fixedMetalVec.emplace_back(box, OBS_NET_IDX);
        //                 // ++numSNetObs;
        //             }
        //             const int topLayerIdx = via.getTopLayer().getRelativeIndex();
        //             for (const Rsyn::PhysicalViaGeometry& geo : via.allTopGeometries()) {
        //                 Bounds bounds = geo.getBounds();
        //                 bounds.translate(pos);
        //                 const BoxOnLayer box(topLayerIdx, getBoxFromRsynBounds(bounds));
        //                 log() << "via_toplayer_box: " << box << std::endl;
        //                 // fixedMetalVec.emplace_back(box, OBS_NET_IDX);
        //                 // ++numSNetObs;
        //             }
        //             if (via.hasViaRule()) {
        //                 const utils::PointT<int> numRowCol =
        //                     via.hasRowCol() ? utils::PointT<int>(via.getNumCols(), via.getNumRows())
        //                                     : utils::PointT<int>(1, 1);
        //                 BoxOnLayer botBox(botLayerIdx);
        //                 BoxOnLayer topBox(topLayerIdx);
        //                 for (unsigned dimIdx = 0; dimIdx != 2; ++dimIdx) {
        //                     const Dimension dim = static_cast<Dimension>(dimIdx);
        //                     const DBU origin = via.hasOrigin() ? pos[dim] + via.getOrigin(dim) : pos[dim];
        //                     const DBU botOff =
        //                         via.hasOffset() ? origin + via.getOffset(Rsyn::BOTTOM_VIA_LEVEL, dim) : origin;
        //                     const DBU topOff =
        //                         via.hasOffset() ? origin + via.getOffset(Rsyn::TOP_VIA_LEVEL, dim) : origin;
        //                     const DBU length =
        //                         (via.getCutSize(dim) * numRowCol[dim] + via.getSpacing(dim) * (numRowCol[dim] - 1)) / 2;
        //                     const DBU botEnc = length + via.getEnclosure(Rsyn::BOTTOM_VIA_LEVEL, dim);
        //                     const DBU topEnc = length + via.getEnclosure(Rsyn::TOP_VIA_LEVEL, dim);
        //                     botBox[dim].Set(botOff - botEnc, botOff + botEnc);
        //                     topBox[dim].Set(topOff - topEnc, topOff + topEnc);
        //                 }
        //                 log() << "botBox: " << botBox << std::endl;
        //                 log() << "topBox: " << topBox << std::endl;
        //                 // fixedMetalVec.emplace_back(botBox, OBS_NET_IDX);
        //                 // fixedMetalVec.emplace_back(topBox, OBS_NET_IDX);
        //                 // numSNetObs += 2;
        //             }
        //             if (layerIdx == botLayerIdx)
        //                 layerIdx = topLayerIdx;
        //             else if (layerIdx == topLayerIdx)
        //                 layerIdx = botLayerIdx;
        //             else {
        //                 // log() << "Error: Special net " << specialNet.getNet().clsName << " via " << pt.clsViaName
        //                 //     << " on wrong layer " << layerIdx << std::endl;
        //                 break;
        //             }
        //         }
        //     }
        // }//end loop wiredescp

        // --- End Net Wire Write 

        // const db::Net& dbNet = nets[i++];
        // defNet.clsWires.clear();
        // if (!dbNet.defWireSegments.empty()) {
        //     defNet.clsWires.emplace_back();
        // defNet.clsWires.back().clsWireSegments = move(dbNet.defWireSegments);
        // vector<DefWireSegmentDscp>& defWireSegments = defNet.clsWires.back().clsWireSegments;
        // std::unordered_map<std::tuple<string, Dimension, DBU>, vector<std::pair<DBU, bool>>> tracks;
        //     for (const DefWireSegmentDscp& seg : dbNet.defWireSegments) {
        //         if (seg.clsRoutingPoints.size() == 1) {
        //             defWireSegments.push_back(seg);
        //             continue;
        //         }
        //         const string& layerName = seg.clsLayerName;
        //         const DBUxy& xy0 = seg.clsRoutingPoints[0].clsPos;
        //         const DBUxy& xy1 = seg.clsRoutingPoints[1].clsPos;
        //         for (unsigned dim = 0; dim != 2; ++dim) {
        //             if (xy0[dim] == xy1[dim]) {
        //                 const std::tuple<string, Dimension, DBU> key =
        //                     std::make_tuple(layerName, static_cast<Dimension>(dim), xy0[dim]);
        //                 std::unordered_map<std::tuple<string, Dimension, DBU>, vector<std::pair<DBU, bool>>>::iterator it =
        //                     tracks.find(key);
        //                 if (it == tracks.end()) {
        //                     it = tracks.emplace(key, vector<std::pair<DBU, bool>>()).first;
        //                 }
        //                 it->second.emplace_back(std::min(xy0[1 - dim], xy1[1 - dim]), true);
        //                 it->second.emplace_back(std::max(xy0[1 - dim], xy1[1 - dim]), false);
        //             }
        //         }
        //     }
        //     for (const std::pair<std::tuple<string, Dimension, DBU>, vector<std::pair<DBU, bool>>>& p : tracks) {
        //         vector<std::pair<DBU, bool>> pts = p.second;
        //         std::sort(pts.begin(), pts.end());
        //         unsigned isWire = 0;
        //         DBU start = std::numeric_limits<DBU>::has_infinity ? -std::numeric_limits<DBU>::infinity()
        //                                                         : std::numeric_limits<DBU>::lowest();
        //         const Dimension dim = std::get<1>(p.first);
        //         for (const std::pair<DBU, bool>& pt : pts) {
        //             if (isWire && pt.first != start) {
        //                 defWireSegments.emplace_back();
        //                 DefWireSegmentDscp& segment = defWireSegments.back();

        //                 segment.clsRoutingPoints.resize(2);
        //                 segment.clsRoutingPoints[0].clsPos[dim] = std::get<2>(p.first);
        //                 segment.clsRoutingPoints[1].clsPos[dim] = std::get<2>(p.first);
        //                 segment.clsRoutingPoints[0].clsPos[1 - dim] = start;
        //                 segment.clsRoutingPoints[1].clsPos[1 - dim] = pt.first;
        //                 segment.clsLayerName = std::get<0>(p.first);
        //             }
        //             if (pt.second) {
        //                 --isWire;
        //             } else {
        //                 ++isWire;
        //             }
        //             start = pt.first;
        //         }
        //     }
        // }
    }  // end for


    //----------------------------------------------------------------




    int numPorts = rsynService.module.getNumPorts(Rsyn::IN) + rsynService.module.getNumPorts(Rsyn::OUT);
    def.clsPorts.reserve(numPorts);
    for (Rsyn::Port port : rsynService.module.allPorts()) {
        Rsyn::PhysicalPort phPort = rsynService.physicalDesign.getPhysicalPort(port);
        def.clsPorts.push_back(DefPortDscp());
        DefPortDscp& defPort = def.clsPorts.back();
        defPort.clsName = port.getName();
        defPort.clsNetName = port.getInnerPin().getNetName();
        if (port.getDirection() == Rsyn::IN)
            defPort.clsDirection = "INPUT";
        else if (port.getDirection() == Rsyn::OUT)
            defPort.clsDirection = "OUTPUT";

        defPort.clsLocationType = "FIXED";
        defPort.clsOrientation = Rsyn::getPhysicalOrientation(phPort.getOrientation());
        defPort.clsLayerName = phPort.getLayer().getName();
        defPort.clsLayerBounds = phPort.getBounds();
        defPort.clsPos = phPort.getPosition();

    }  // end for

    log() << "write Def..." << filename << std::endl;
    defParser.writeFullDEF(filename.c_str(), def);
}//emd writeDef


}  // namespace db

MTStat runJobsMT(int numJobs, const std::function<void(int)>& handle) {
    int numThreads = min(numJobs, db::setting.numThreads);
    MTStat mtStat(max(1, db::setting.numThreads));
    if (numThreads <= 1) {
        utils::timer threadTimer;
        for (int i = 0; i < numJobs; ++i) {
            handle(i);
        }
        mtStat.durations[0] = threadTimer.elapsed();
    } else {
        int globalJobIdx = 0;
        std::mutex mtx;
        utils::timer threadTimer;
        auto thread_func = [&](int threadIdx) {
            int jobIdx;
            while (true) {
                mtx.lock();
                jobIdx = globalJobIdx++;
                mtx.unlock();
                if (jobIdx >= numJobs) {
                    mtStat.durations[threadIdx] = threadTimer.elapsed();
                    break;
                }
                handle(jobIdx);
            }
        };

        std::thread threads[numThreads];
        for (int i = 0; i < numThreads; i++) {
            threads[i] = std::thread(thread_func, i);
        }
        for (int i = 0; i < numThreads; i++) {
            threads[i].join();
        }
    }
    return mtStat;
}
