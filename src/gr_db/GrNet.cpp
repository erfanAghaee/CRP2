#include "GrNet.h"
#include "GrDatabase.h"
#include <random>
#include <iostream>

namespace gr {

void GrNet::init(const GCellGrid& gcellGrid) {

    
    max_edge_cost = 0;
    initPinAccessBoxes(gcellGrid);

    // get overlap points
    for (int p1 = 0; p1 < numOfPins(); p1++) {
        for (int p2 = p1 + 1; p2 < numOfPins(); p2++) {
            for (const auto& point1 : pinAccessBoxes[p1]) {
                for (const auto& point2 : pinAccessBoxes[p2]) {
                    if (point1 == point2) ovlpPoints.insert(point1);
                }
            }
        }
    }

    // judge if it is one pin net
    std::unordered_map<gr::GrPoint, int> pointSet;
    for (const auto& pin : pinAccessBoxes) {
        for (const auto& point : pin) pointSet[point]++;
    }

    for (const auto& pair : pointSet) {
        if (pair.second == numOfPins()) {
            isOnePin = true;
            break;
        }
    }

    for (const auto& pinBox : pinAccessBoxes) {
        for (const auto& grpnt : pinBox) {
            boundingBox[X].Update(grpnt[X]);
            boundingBox[Y].Update(grpnt[Y]);
        }
    }
}

void GrNet::update(const GCellGrid& gcellGrid) {
    isOnePin = false;
    // ovlpPoints.clear();
    updatePinAccessBoxes(gcellGrid);

    // get overlap points
    for (int p1 = 0; p1 < numOfPins(); p1++) {
        for (int p2 = p1 + 1; p2 < numOfPins(); p2++) {
            for (const auto& point1 : pinAccessBoxes[p1]) {
                for (const auto& point2 : pinAccessBoxes[p2]) {
                    if (point1 == point2) ovlpPoints.insert(point1);
                }
            }
        }
    }

    // judge if it is one pin net
    std::unordered_map<gr::GrPoint, int> pointSet;
    for (const auto& pin : pinAccessBoxes) {
        for (const auto& point : pin) pointSet[point]++;
    }

    // log() << "grNet update: " << getName()
    //       << ", numOfPins: " << numOfPins() << std::endl;

    for (const auto& pair : pointSet) {
        if (pair.second == numOfPins()) {
            isOnePin = true;
            break;
        }
    }

    for (const auto& pinBox : pinAccessBoxes) {
        for (const auto& grpnt : pinBox) {
            boundingBox[X].Update(grpnt[X]);
            boundingBox[Y].Update(grpnt[Y]);
        }
    }
}

// merged pins with same coor
vector<vector<PointOnLayer>> GrNet::getMergedPinAccessBoxes(
    const std::function<PointOnLayer(GrPoint)>& pointHash) const {
    vector<vector<PointOnLayer>> mergedPinAccessBoxes;

    vector<vector<PointOnLayer>> hashedPinAccessBoxes(numOfPins());
    for (int p = 0; p < numOfPins(); p++)
        for (const auto& point : pinAccessBoxes[p]) hashedPinAccessBoxes[p].push_back(pointHash(point));

    vector<vector<int>> pinConn(numOfPins());
    for (int p1 = 0; p1 < numOfPins(); p1++) {
        for (int p2 = p1 + 1; p2 < numOfPins(); p2++) {
            bool overlap = false;
            for (const auto& point1 : hashedPinAccessBoxes[p1]) {
                if (overlap) break;
                for (const auto& point2 : hashedPinAccessBoxes[p2]) {
                    if (overlap) break;
                    if (point1 == point2) {
                        pinConn[p1].push_back(p2);
                        pinConn[p2].push_back(p1);
                        overlap = true;
                    }
                }
            }
        }
    }
    vector<bool> visited(numOfPins(), false);
    for (int p = 0; p < numOfPins(); p++) {
        if (visited[p]) continue;
        std::queue<int> q;
        q.push(p);
        vector<PointOnLayer> points;
        while (!q.empty()) {
            int node = q.front();
            q.pop();
            copy(hashedPinAccessBoxes[node].begin(), hashedPinAccessBoxes[node].end(), std::back_inserter(points));
            visited[node] = true;
            for (auto child : pinConn[node])
                if (!visited[child]) q.push(child);
        }
        mergedPinAccessBoxes.resize(mergedPinAccessBoxes.size() + 1);
        mergedPinAccessBoxes.back() = move(points);
    }
    for (auto& points : mergedPinAccessBoxes) {
        std::sort(points.begin(), points.end(), [&](const PointOnLayer& lhs, const PointOnLayer& rhs) {
            if (lhs.layerIdx != rhs.layerIdx) {
                return lhs.layerIdx < rhs.layerIdx;
            } else {
                if (lhs[X] != rhs[X]) {
                    return lhs[X] < rhs[X];
                } else {
                    return lhs[Y] < rhs[Y];
                }
            }
        });
        points.erase(std::unique(points.begin(), points.end()), points.end());
    }

    return mergedPinAccessBoxes;
}

void GrNet::initPinAccessBoxes(const GCellGrid& gcellGrid) {
    bool debug = false;
    // if(getName() == "net1233") debug = false;
    // transform coor to grPoint, construct pinAccessBoxes
    if(debug) {
        log() << "initPinAccessBoxes...(grNet)" << getName() << std::endl;
    }
    // log() << "grNet initPinAccessBoxes: " << getName() << std::endl;

    pinAccessBoxes.resize(numOfPins());
    for (int i = 0; i < numOfPins(); i++) {
        const auto& boxes = dbNet.pinAccessBoxes[i];
        std::unordered_set<GrPoint> pointSet;

        

        DBU smallestVio = std::numeric_limits<DBU>::max();
        vector<const db::BoxOnLayer*> smallestBoxes;

        for (const auto& box : boxes) {
            if(debug)
                log() << "i: " << i    
                    << ", box: " << box << std::endl;
            int vio = database.getOvlpFixedMetalArea(box, dbNet.idx);
            if (vio <= smallestVio) {
                if (vio < smallestVio) smallestBoxes.clear();
                smallestVio = vio;
                smallestBoxes.push_back(&box);
            }

            if (vio == 0) {
                auto grBox = gcellGrid.rangeSearchGCell(box);
                // if(debug)
                //     log() << "i: " << i 
                //           << ", box: " << box 
                //           << ", grBox: " << grBox << std::endl;
                for (int x = grBox[X].low; x <= grBox[X].high; x++)
                    for (int y = grBox[Y].low; y <= grBox[Y].high; y++) {
                        pointSet.emplace(box.layerIdx, x, y);
                        if(debug)
                            log() << "box: " << box 
                                << ", box.layerIdx: " << box.layerIdx
                                << ", gcell_x: " << x
                                << ", gcell_y: " << y << std::endl;
                              
                    }
            }
        }
        // all have vio, add those with smallest vio
        if (pointSet.empty()) {
            for (auto box : smallestBoxes) {
                auto grBox = gcellGrid.rangeSearchGCell(*box);
                for (int x = grBox[X].low; x <= grBox[X].high; x++)
                    for (int y = grBox[Y].low; y <= grBox[Y].high; y++){
                        pointSet.emplace(box->layerIdx, x, y);
                    } 
            }
        }

        for (auto& point : pointSet){
            pinAccessBoxes[i].push_back(point);
            if(debug)
                log() << "pt.x: " << point.x
                      << ", pt.y: " << point.y << std::endl;
        }
    }
    debug = false;
}//end initPinAccessBoxes

void GrNet::updatePinAccessBoxes(const GCellGrid& gcellGrid) {
    bool debug = false;
    // if(getName() == "net1233") debug = false;
    if(debug) {
        log() << "updatePinAccessBoxes...(grNet)" << getName() << std::endl;
    }
    // transform coor to grPoint, construct pinAccessBoxes
    pinAccessBoxes.clear();
    pinAccessBoxes.resize(numOfPins());
    for (int i = 0; i < numOfPins(); i++) {
        const auto& boxes = dbNet.pinAccessBoxes[i];
        std::unordered_set<GrPoint> pointSet;

        

        DBU smallestVio = std::numeric_limits<DBU>::max();
        vector<const db::BoxOnLayer*> smallestBoxes;

        for (const auto& box : boxes) {
            if(debug)
                log() << "i: " << i    
                    << ", box: " << box << std::endl;
            int vio = database.getOvlpFixedMetalArea(box, dbNet.idx);
            if (vio <= smallestVio) {
                if (vio < smallestVio) smallestBoxes.clear();
                smallestVio = vio;
                smallestBoxes.push_back(&box);
            }

            if (vio == 0) {
                auto grBox = gcellGrid.rangeSearchGCell(box);
                if(debug)
                    log() << "i: " << i    
                          << ", box: " << box                      
                          << ", grBox: " << grBox << std::endl;
                for (int x = grBox[X].low; x <= grBox[X].high; x++)
                    for (int y = grBox[Y].low; y <= grBox[Y].high; y++) pointSet.emplace(box.layerIdx, x, y);
            }
        }
        // all have vio, add those with smallest vio
        if (pointSet.empty()) {
            for (auto box : smallestBoxes) {
                auto grBox = gcellGrid.rangeSearchGCell(*box);
                for (int x = grBox[X].low; x <= grBox[X].high; x++)
                    for (int y = grBox[Y].low; y <= grBox[Y].high; y++) pointSet.emplace(box->layerIdx, x, y);
            }
        }

        for (auto& point : pointSet) {
            pinAccessBoxes[i].push_back(point);
            if(debug)
                log() << "pt.x: " << point.x
                      << ", pt.y: " << point.y << std::endl;
        }
    }
    debug = false;
}//update pinAccesboxes

void GrNetlist::init(const GCellGrid& gcellGrid) {
    nets.clear();
    for (int i = 0, sz = database.nets.size(); i < sz; i++) nets.emplace_back(i);
    log() << "database.nets.size()" << database.nets.size() << std::endl;
    runJobsMT(database.nets.size(), [&](int i) { nets[i].init(gcellGrid); });
    
}

void GrNetlist::update(const GCellGrid& gcellGrid) {
    // for (int i = 0, sz = database.nets.size(); i < sz; i++) nets.emplace_back(i);
    runJobsMT(database.nets.size(), [&](int i) { nets[i].update(gcellGrid); });
}

void GrNet::postOrderVisitGridTopo(const std::function<void(std::shared_ptr<GrSteiner>)>& visit) const {
    for (const std::shared_ptr<GrSteiner>& tree : gridTopo) {
        GrSteiner::postOrder(tree, visit);
    }
}

void GrNet::preOrderVisitGridTopo(const std::function<void(std::shared_ptr<GrSteiner>)>& visit) const {
    for (const std::shared_ptr<GrSteiner>& tree : gridTopo) {
        GrSteiner::preOrder(tree, visit);
    }
}

DBU GrNet::getWirelength() const {
    DBU wirelength = 0;
    postOrderVisitGridTopo([&](std::shared_ptr<gr::GrSteiner> node) {
        auto parent = node;
        for (auto child : parent->children) {
            if (parent->layerIdx == child->layerIdx) wirelength += grDatabase.getDist(*parent, *child);
        }
    });
    return wirelength;
}

int GrNet::getCriticalCell(){
    DBU path_cost = 0;
    DBU wl_cost = 0;
    DBU wl_absolute_cost = 0;
    DBU via_cost = 0;
    bool debug = false;


    std::vector<int> cells_idx;
    for (auto pin : dbNet.rsynPins){
        std::string inst_name = pin.getInstanceName();
        int cell_idx = database.cellsToIdx[inst_name];
        // auto cell = database.cells[cell_idx];
        cells_idx.push_back(cell_idx);
    }   


    // for(auto idx : cells_idx){
    //     log() << "connected cells: " << database.cells[idx].getName() << std::endl;;
    // }


    // for(auto pinAccessBox: pinAccessBoxes){
    //     log () << "pin boxs ..." << std::endl;
    //     for(auto pin : pinAccessBox ){
    //         log() << "pin: [x: " << pin.x << ", y: " << pin.y 
    //               << ", l: " << pin.layerIdx << "]" << std::endl;
    //     }
    // }


    if(debug)
            log() << "getCriticalCell ... " << std::endl;
    auto pathCostFun = [&](std::shared_ptr<gr::GrSteiner> node) {
        auto parent = node;
        for (auto child : parent->children) {
            if (parent->layerIdx == child->layerIdx) {
                std::shared_ptr<gr::GrSteiner> lower, upper;
                if ((*parent)[X] < (*child)[X] || (*parent)[Y] < (*child)[Y]) {
                    lower = parent;
                    upper = child;
                } else {
                    lower = child;
                    upper = parent;
                }

                gr::GrEdge edge(gr::GrPoint(lower->layerIdx, lower->x , lower->y ),
                                gr::GrPoint(upper->layerIdx, upper->x, upper->y) );

                DBU edge_cost = grDatabase.getWireCost(edge);
                path_cost += edge_cost;
                wl_cost += edge_cost;
                wl_absolute_cost += grDatabase.getDist(*parent, *child);

                if(debug)
                    log() << "lower: [x: " << lower->x << ", y: " << lower->y 
                        << ", l: " << lower->layerIdx
                        << "], upper: [x: " << upper->x << ", y: " << upper->y 
                        << ", l: " << upper->layerIdx << "]"
                        << ", edge_cost: " << edge_cost << std::endl;

                // gr::GrEdge edge_test(gr::GrPoint(2, 52, 71),
                //                      gr::GrPoint(2, 55, 71) );

                // DBU edge_cost_test = grDatabase.getWireCost(edge_test);

                // if(debug)
                //     log() << "edge_cost_test: " << edge_cost_test << std::endl;


                // wl_cost += grDatabase.getDist(*parent, *child);

                // wireRouteGuides.emplace_back(lower->layerIdx,
                //                             utils::IntervalT<int>((*lower)[X], (*upper)[X]),
                //                             utils::IntervalT<int>((*lower)[Y], (*upper)[Y]));
            } else {
               
                auto lower_layer = std::min(parent->layerIdx,child->layerIdx);
                auto upper_layer = std::max(parent->layerIdx,child->layerIdx);

                // viaRouteGuides.emplace_back(parent->layerIdx,
                //                                 utils::IntervalT<int>((*parent)[X], (*parent)[X]),
                //                                 utils::IntervalT<int>((*parent)[Y], (*parent)[Y]));
                // viaRouteGuides.emplace_back(child->layerIdx,
                //                                 utils::IntervalT<int>((*child)[X], (*child)[X]),
                //                                 utils::IntervalT<int>((*child)[Y], (*child)[Y]));
                DBU via_pt = grDatabase.getStackViaCost(gr::GrPoint(lower_layer, parent->x , parent->y),
                                                std::abs(upper_layer-lower_layer));
                path_cost +=via_pt;
                    
                via_cost +=via_pt;
                if(debug)
                    log() << "parent: [x: " << parent->x << ", y: " << parent->y 
                        << ", l: " << parent->layerIdx
                        << "], child: [x: " << child->x << ", y: " << child->y 
                        << ", l: " << child->layerIdx << "]" 
                        <<", via_cost: " << via_pt << std::endl;
            }
        }
    };

    vector<DBU> path_costs;
    int i = 0;
    for(auto gridTopoTmp : gridTopoHist){
        path_cost = 0;
        wl_cost = 0;
        via_cost = 0;
        wl_absolute_cost = 0;
        if(debug)
            log() << "gridtopo#" << i << std::endl;
        for (const std::shared_ptr<gr::GrSteiner>& tree : gridTopoTmp) {
            gr::GrSteiner::postOrder(tree,pathCostFun );
        }
        path_costs.push_back(path_cost);
        if(debug){
            log() << "gridTopo path cost: " << path_cost << std::endl;
            log() << "gridTopo wl cost: " << wl_cost << std::endl;
            log() << "gridTopo via cost: " << via_cost << std::endl;   
            log() << "gridTopo wl_absolute_cost: " << wl_absolute_cost  << std::endl;   
            
        }
        
        i++;
    }


    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> rand_choice(0,cells_idx.size()-1); // distribution in range [1, 6]


    // return cells_idx[rand_choice(rng)];
    return cells_idx[0];


    // if(path_costs.size()>0){
    //     if(path_costs[0] - path_costs[path_costs.size() - 1] < 0){
    //         if(debug)
    //             log() << "grNet: " << getName() << ", cost increased... " << std::endl;
    //         critical_net = true;
    //     }
    // }
    // log() << "path_cost: " << path_cost << std::endl;
    // return path_cost;
}// end getCriticalCell


void GrNet::updateMaxEdgeCost(){
    bool debug = false;
    double max_edge_cost = 0;

    // if(getName() == "net610") debug = false;
    if(debug) log() << "updateMaxEdgeCost net: "  << getName() << std::endl;

    auto pathCostFun = [&](std::shared_ptr<gr::GrSteiner> node) {
        auto parent = node;
        for (auto child : parent->children) {
            if (parent->layerIdx == child->layerIdx) {
                std::shared_ptr<gr::GrSteiner> lower, upper;
                if ((*parent)[X] < (*child)[X] || (*parent)[Y] < (*child)[Y]) {
                    lower = parent;
                    upper = child;
                } else {
                    lower = child;
                    upper = parent;
                }
                if(debug){
                    log() << std::endl << std::endl;
                }

                gr::GrEdge edge(gr::GrPoint(lower->layerIdx, lower->x , lower->y ),
                                gr::GrPoint(upper->layerIdx, upper->x, upper->y) );

                DBU edge_cost = grDatabase.getWireCost(edge,debug);
                if(edge_cost > max_edge_cost)
                    max_edge_cost = edge_cost;
                
                if(debug){
                    log() << "lower: [x: " << lower->x << ", y: " << lower->y 
                        << ", l: " << lower->layerIdx
                        << "], upper: [x: " << upper->x << ", y: " << upper->y 
                        << ", l: " << upper->layerIdx << "]"
                        << ", edge_cost: " << edge_cost << std::endl;
                }
                    


            } else {
               if(debug){
                    log() << std::endl << std::endl;
                }
                auto lower_layer = std::min(parent->layerIdx,child->layerIdx);
                auto upper_layer = std::max(parent->layerIdx,child->layerIdx);

                // viaRouteGuides.emplace_back(parent->layerIdx,
                //                                 utils::IntervalT<int>((*parent)[X], (*parent)[X]),
                //                                 utils::IntervalT<int>((*parent)[Y], (*parent)[Y]));
                // viaRouteGuides.emplace_back(child->layerIdx,
                //                                 utils::IntervalT<int>((*child)[X], (*child)[X]),
                //                                 utils::IntervalT<int>((*child)[Y], (*child)[Y]));
                DBU via_pt = grDatabase.getStackViaCost(gr::GrPoint(lower_layer, parent->x , parent->y),
                                                std::abs(upper_layer-lower_layer),debug);

                if(debug)
                    log() << "parent: [x: " << parent->x << ", y: " << parent->y 
                        << ", l: " << parent->layerIdx
                        << "], child: [x: " << child->x << ", y: " << child->y 
                        << ", l: " << child->layerIdx << "]" 
                        <<", via_cost: " << via_pt << std::endl;
            }
        }
    };

    if(! (gridTopoHist.size() > 0) ){
        max_edge_cost = 0;
    } else{
        auto gridTopoTmp = gridTopoHist[gridTopoHist.size()-1];
        for (const std::shared_ptr<gr::GrSteiner>& tree : gridTopoTmp) {
            gr::GrSteiner::postOrder(tree,pathCostFun );
        }
    }
}


void GrNet::checkRouteHist(){
    DBU path_cost = 0;
    DBU wl_cost = 0;
    DBU via_cost = 0;
    DBU wl_absolute_cost = 0;
    bool debug = false;


    std::vector<int> cells_idx;
    for (auto pin : dbNet.rsynPins){
        std::string inst_name = pin.getInstanceName();
        int cell_idx = database.cellsToIdx[inst_name];
        // auto cell = database.cells[cell_idx];
        cells_idx.push_back(cell_idx);
    }   

    if(debug)
        for(auto idx : cells_idx){
            log() << "connected cells: " << database.cells[idx].getName() << std::endl;;
        }

    if(debug)
        for(auto pinAccessBox: pinAccessBoxes){
            log () << "pin boxs ..." << std::endl;
            for(auto pin : pinAccessBox ){
                log() << "pin: [x: " << pin.x << ", y: " << pin.y 
                    << ", l: " << pin.layerIdx << "]" << std::endl;
            }
        }


    if(debug)
            log() << "getCriticalCell ... " << std::endl;
    auto pathCostFun = [&](std::shared_ptr<gr::GrSteiner> node) {
        auto parent = node;
        for (auto child : parent->children) {
            if (parent->layerIdx == child->layerIdx) {
                std::shared_ptr<gr::GrSteiner> lower, upper;
                if ((*parent)[X] < (*child)[X] || (*parent)[Y] < (*child)[Y]) {
                    lower = parent;
                    upper = child;
                } else {
                    lower = child;
                    upper = parent;
                }

                if(debug)
                    log() << "lower: [x: " << lower->x << ", y: " << lower->y 
                        << ", l: " << lower->layerIdx
                        << "], upper: [x: " << upper->x << ", y: " << upper->y 
                        << ", l: " << upper->layerIdx << "]"<< std::endl;


                gr::GrEdge edge(gr::GrPoint(lower->layerIdx, lower->x , lower->y ),
                                gr::GrPoint(upper->layerIdx, upper->x, upper->y) );

                DBU edge_cost = grDatabase.getWireCost(edge);
                path_cost += edge_cost;
                wl_cost += edge_cost;
                wl_absolute_cost += grDatabase.getDist(*parent, *child);

                // wireRouteGuides.emplace_back(lower->layerIdx,
                //                             utils::IntervalT<int>((*lower)[X], (*upper)[X]),
                //                             utils::IntervalT<int>((*lower)[Y], (*upper)[Y]));
            } else {

                if(debug)
                    log() << "parent: [x: " << parent->x << ", y: " << parent->y 
                        << ", l: " << parent->layerIdx
                        << "], child: [x: " << child->x << ", y: " << child->y 
                        << ", l: " << child->layerIdx << "]" <<std::endl;
                
                auto lower_layer = std::min(parent->layerIdx,child->layerIdx);
                auto upper_layer = std::max(parent->layerIdx,child->layerIdx);

                // viaRouteGuides.emplace_back(parent->layerIdx,
                //                                 utils::IntervalT<int>((*parent)[X], (*parent)[X]),
                //                                 utils::IntervalT<int>((*parent)[Y], (*parent)[Y]));
                // viaRouteGuides.emplace_back(child->layerIdx,
                //                                 utils::IntervalT<int>((*child)[X], (*child)[X]),
                //                                 utils::IntervalT<int>((*child)[Y], (*child)[Y]));

                DBU via_pt = grDatabase.getStackViaCost(gr::GrPoint(lower_layer, parent->x , parent->y),
                                                std::abs(upper_layer-lower_layer));
                path_cost +=via_pt;
                    
                via_cost +=via_pt;
            }
        }
    };

    vector<DBU> path_costs;
    vector<DBU> wl_costs;
    vector<DBU> via_costs;
    vector<DBU> wl_absolute_costs;

    int i = 0;
    for(auto gridTopoTmp : gridTopoHist){
        path_cost = 0;
        wl_cost = 0;
        via_cost = 0;
        wl_absolute_cost = 0;
        if(debug)
            log() << "gridtopo#" << i << std::endl;
        for (const std::shared_ptr<gr::GrSteiner>& tree : gridTopoTmp) {
            gr::GrSteiner::postOrder(tree,pathCostFun );
        }
        path_costs.push_back(path_cost);
        wl_costs.push_back(wl_cost);
        via_costs.push_back(via_cost);
        wl_absolute_costs.push_back(wl_absolute_cost);

        if(debug){
            log() << "gridTopo path cost: " << path_cost << std::endl;
            log() << "gridTopo wl cost: " << wl_cost << std::endl;
            log() << "gridTopo wl_absolute_cost: " << wl_absolute_cost << std::endl;
            log() << "gridTopo via cost: " << via_cost << std::endl;    
        }
        
        i++;
    }

    if(path_costs.size()>0){
        // if(path_costs[0] - path_costs[path_costs.size() - 1] < 0){
        //     if(debug)
        //         log() << "grNet: " << getName() << ",path cost increased... " << std::endl;
        //     is_critical_net = true;
        // }

        if(wl_absolute_costs[0] - wl_absolute_costs[wl_absolute_costs.size() - 1] < 0){
            if(debug)
                log() << "grNet: " << getName() << ",wl_absolute_costs increased... " << std::endl;
            is_critical_net = true;
        }

    }
    // log() << "path_cost: " << path_cost << std::endl;
    // return path_cost;
}
void GrNet::getRouteSites(){
    bool debug = false;

    if(debug)
        log() << "get route sites: " << getName() << std::endl;

    std::vector<std::pair<std::tuple<int,int,int>,double>> wireMap_tpls;
    std::vector<std::pair<std::tuple<int,int,int>,double>> viaMap_tpls;
    for (const auto& guide : wireRouteGuides) grDatabase.removeWireRelax(guide,wireMap_tpls);
    const auto& viaGuides = viaRouteGuides;
    for (int g1 = 0; g1 < viaGuides.size(); g1++) {
        for (int g2 = g1 + 1; g2 < viaGuides.size(); g2++) {
            if (abs(viaGuides[g1].layerIdx - viaGuides[g2].layerIdx) != 1) continue;

            auto xIntvl = viaGuides[g1][X].IntersectWith(viaGuides[g2][X]);
            auto yIntvl = viaGuides[g1][Y].IntersectWith(viaGuides[g2][Y]);

            if (xIntvl.IsValid() && yIntvl.IsValid()){
                grDatabase.removeViaRelax({min(viaGuides[g1].layerIdx, viaGuides[g2].layerIdx), xIntvl, yIntvl},viaMap_tpls);
            }
                // removeVia({min(viaGuides[g1].layerIdx, viaGuides[g2].layerIdx), xIntvl, yIntvl});
                
                
        }
    }


    std::vector<utils::BoxT<DBU>> wires;
    std::vector<std::vector<utils::BoxT<DBU>>> rows_matrix;
    std::vector<std::vector<int>> sites_matrix;


    auto pathCostFun = [&](std::shared_ptr<gr::GrSteiner> node) {
        auto parent = node;
        for (auto child : parent->children) {
            if (parent->layerIdx == child->layerIdx) {
                std::shared_ptr<gr::GrSteiner> lower, upper;
                if ((*parent)[X] < (*child)[X] || (*parent)[Y] < (*child)[Y]) {
                    lower = parent;
                    upper = child;
                } else {
                    lower = child;
                    upper = parent;
                }

                gr::GrEdge edge(gr::GrPoint(lower->layerIdx, lower->x , lower->y ),
                                gr::GrPoint(upper->layerIdx, upper->x, upper->y) );

                         

                DBU edge_cost = grDatabase.getWireCostRelax(edge,wireMap_tpls,viaMap_tpls,false);

                auto lx =  grDatabase.getCoor((*lower)[X], X);
                auto ly =  grDatabase.getCoor((*lower)[Y], Y);
                auto hx =  grDatabase.getCoor((*upper)[X] + 1, X);
                auto hy =  grDatabase.getCoor((*upper)[Y]+ 1, Y);
                auto box_edge = utils::BoxT<DBU>(lx,ly,hx,hy);

                std::vector<utils::BoxT<DBU>> rows;
                std::vector<int> sites;
                database.getRowsInBox(box_edge, rows, 1000);
                database.getIntersectedSites(box_edge,sites);

                rows_matrix.push_back(rows);
                sites_matrix.push_back(sites);

                // ss << database.getLayer(guide.layerIdx).name << std::endl;

                if(debug){
                    log() << "(*lower)[X]: " << (*lower)[X] 
                        << ", (*upper)[X]: "  << (*upper)[X]
                        << ", (*lower)[Y]: "  << (*lower)[Y]
                        << ", (*upper)[Y]: "  << (*upper)[Y] << std::endl;
                    log() << "box_edge: " << box_edge << ",layerIdx: " << edge.getLayerIdx() << std::endl;
                }
                    
                    


                // wireRouteGuides.emplace_back(lower->layerIdx,
                //                             utils::IntervalT<int>((*lower)[X], (*upper)[X]),
                //                             utils::IntervalT<int>((*lower)[Y], (*upper)[Y]));
                if(debug)
                    log() << "lower: [x: " << lower->x << ", y: " << lower->y 
                        << ", l: " << lower->layerIdx
                        << "], upper: [x: " << upper->x << ", y: " << upper->y 
                        << "l: " << upper->layerIdx << "]"
                        << ", edge_cost: " << edge_cost
                        << std::endl;
            } else {
               
                auto lower_layer = std::min(parent->layerIdx,child->layerIdx);
                auto upper_layer = std::max(parent->layerIdx,child->layerIdx);

                auto lx =  grDatabase.getCoor((*parent)[X], X);
                auto ly =  grDatabase.getCoor((*parent)[Y], Y);
                auto hx =  grDatabase.getCoor((*parent)[X] + 1, X);
                auto hy =  grDatabase.getCoor((*parent)[Y]+ 1, Y);
                auto box_via = utils::BoxT<DBU>(lx,ly,hx,hy);

                // viaRouteGuides.emplace_back(parent->layerIdx,
                //                                 utils::IntervalT<int>((*parent)[X], (*parent)[X]),
                //                                 utils::IntervalT<int>((*parent)[Y], (*parent)[Y]));
                // viaRouteGuides.emplace_back(child->layerIdx,
                //                                 utils::IntervalT<int>((*child)[X], (*child)[X]),
                //                                 utils::IntervalT<int>((*child)[Y], (*child)[Y]));
                DBU via_pt =  grDatabase.getStackViaCostRelax(gr::GrPoint(lower_layer, parent->x , parent->y),
                                                std::abs(upper_layer-lower_layer),wireMap_tpls,viaMap_tpls,false);
                // DBU via_pt =  grDatabase.getStackViaCost(gr::GrPoint(lower_layer, parent->x , parent->y),
                //                                 std::abs(upper_layer-lower_layer),true);
                   
                if(debug){
                    log() << "parent: [x: " << parent->x << ", y: " << parent->y 
                        << ", l: " << parent->layerIdx
                        << ", child: [x: " << child->x << ", y: " << child->y 
                        << "l: " << child->layerIdx << "]" 
                        << ", via_pt: " << via_pt
                        << std::endl;
                    log() << "box_via: " << box_via 
                          << ", lower_layer: " << lower_layer
                          << ", upper_layer: " << upper_layer << std::endl;
                }
                    
            }
        }
    };


    for (const std::shared_ptr<gr::GrSteiner>& tree : gridTopo) {
        gr::GrSteiner::postOrder(tree,pathCostFun );
    }

    log() << "row and site matrix:" << std::endl;
    for(auto rows: rows_matrix ){
        log() << "new edge row" << std::endl;
        for(auto row: rows){
            log() << "row: " << row << std::endl;
        }
    }
    for(auto sites: sites_matrix ){
        log() << "new edge site" << std::endl;
        for(auto site: sites){
            log() << "row: " << site << std::endl;
        }
    }
}

DBU GrNet::getPathCost() const {
    DBU path_cost = 0;
    bool debug = false;
    auto pathCostFun = [&](std::shared_ptr<gr::GrSteiner> node) {
        if(debug)
            log() << "genTopoGuide ... " << std::endl;
        auto parent = node;
        for (auto child : parent->children) {
            if (parent->layerIdx == child->layerIdx) {
                std::shared_ptr<gr::GrSteiner> lower, upper;
                if ((*parent)[X] < (*child)[X] || (*parent)[Y] < (*child)[Y]) {
                    lower = parent;
                    upper = child;
                } else {
                    lower = child;
                    upper = parent;
                }

                gr::GrEdge edge(gr::GrPoint(lower->layerIdx, lower->x , lower->y ),
                                gr::GrPoint(upper->layerIdx, upper->x, upper->y) );

                DBU edge_cost = grDatabase.getWireCost(edge);
                path_cost += edge_cost;

                // wireRouteGuides.emplace_back(lower->layerIdx,
                //                             utils::IntervalT<int>((*lower)[X], (*upper)[X]),
                //                             utils::IntervalT<int>((*lower)[Y], (*upper)[Y]));
                if(debug)
                    log() << "lower: [x: " << lower->x << ", y: " << lower->y 
                        << ", l: " << lower->layerIdx
                        << "], upper: [x: " << upper->x << ", y: " << upper->y 
                        << "l: " << upper->layerIdx << "]"
                        << ", edge_cost: " << edge_cost
                        << std::endl;
            } else {
               
                auto lower_layer = std::min(parent->layerIdx,child->layerIdx);
                auto upper_layer = std::max(parent->layerIdx,child->layerIdx);

                // viaRouteGuides.emplace_back(parent->layerIdx,
                //                                 utils::IntervalT<int>((*parent)[X], (*parent)[X]),
                //                                 utils::IntervalT<int>((*parent)[Y], (*parent)[Y]));
                // viaRouteGuides.emplace_back(child->layerIdx,
                //                                 utils::IntervalT<int>((*child)[X], (*child)[X]),
                //                                 utils::IntervalT<int>((*child)[Y], (*child)[Y]));
                DBU via_pt =  grDatabase.getStackViaCost(gr::GrPoint(lower_layer, parent->x , parent->y),
                                                std::abs(upper_layer-lower_layer));
                path_cost += via_pt;
                   
                if(debug)
                    log() << "parent: [x: " << parent->x << ", y: " << parent->y 
                        << ", l: " << parent->layerIdx
                        << ", child: [x: " << child->x << ", y: " << child->y 
                        << "l: " << child->layerIdx << "]" 
                        << ", via_pt: " << via_pt
                        << std::endl;
            }
        }
    };


    for (const std::shared_ptr<gr::GrSteiner>& tree : gridTopo) {
        gr::GrSteiner::postOrder(tree,pathCostFun );
    }
    // log() << "path_cost: " << path_cost << std::endl;
    return path_cost;
}//end getPathCost

DBU GrNet::getLongestWire() const {
    DBU longest_wire = 0;
    
    auto pathCostFun = [&](std::shared_ptr<gr::GrSteiner> node) {
        auto parent = node;
        for (auto child : parent->children) {
            if (parent->layerIdx == child->layerIdx) {
                std::shared_ptr<gr::GrSteiner> lower, upper;
                if ((*parent)[X] < (*child)[X] || (*parent)[Y] < (*child)[Y]) {
                    lower = parent;
                    upper = child;
                } else {
                    lower = child;
                    upper = parent;
                }

                gr::GrEdge edge(gr::GrPoint(lower->layerIdx, lower->x , lower->y ),
                                gr::GrPoint(upper->layerIdx, upper->x, upper->y) );

                // DBU edge_cost = grDatabase.getWireCost(edge);
                // path_cost += edge_cost;
                auto wire_len = grDatabase.getDist(*parent, *child);
                if(wire_len > longest_wire){
                    longest_wire = wire_len;
                }
                // wireRouteGuides.emplace_back(lower->layerIdx,
                //                             utils::IntervalT<int>((*lower)[X], (*upper)[X]),
                //                             utils::IntervalT<int>((*lower)[Y], (*upper)[Y]));
                
            } else {
               
                auto lower_layer = std::min(parent->layerIdx,child->layerIdx);
                auto upper_layer = std::max(parent->layerIdx,child->layerIdx);

                
                // DBU via_pt =  grDatabase.getStackViaCost(gr::GrPoint(lower_layer, parent->x , parent->y),
                //                                 std::abs(upper_layer-lower_layer));
                // path_cost += via_pt;
            }
        }
    };


    for (const std::shared_ptr<gr::GrSteiner>& tree : gridTopo) {
        gr::GrSteiner::postOrder(tree,pathCostFun );
    }
    // log() << "path_cost: " << path_cost << std::endl;
    return longest_wire;
}//end getLongestWire


double GrNet::getEdgeHighestCost() const{
    double edge_highest_cost = 0;
    std::vector<std::pair<std::tuple<int,int,int>,double>> wireMap_tpls;
    std::vector<std::pair<std::tuple<int,int,int>,double>> viaMap_tpls;
    for (const auto& guide : wireRouteGuides) grDatabase.removeWireRelax(guide,wireMap_tpls);
    const auto& viaGuides = viaRouteGuides;
    for (int g1 = 0; g1 < viaGuides.size(); g1++) {
        for (int g2 = g1 + 1; g2 < viaGuides.size(); g2++) {
            if (abs(viaGuides[g1].layerIdx - viaGuides[g2].layerIdx) != 1) continue;

            auto xIntvl = viaGuides[g1][X].IntersectWith(viaGuides[g2][X]);
            auto yIntvl = viaGuides[g1][Y].IntersectWith(viaGuides[g2][Y]);

            if (xIntvl.IsValid() && yIntvl.IsValid()){
                grDatabase.removeViaRelax({min(viaGuides[g1].layerIdx, viaGuides[g2].layerIdx), xIntvl, yIntvl},viaMap_tpls);
            }
                // removeVia({min(viaGuides[g1].layerIdx, viaGuides[g2].layerIdx), xIntvl, yIntvl});
                
                
        }
    }


    std::vector<utils::BoxT<DBU>> wires;


    auto pathCostFun = [&](std::shared_ptr<gr::GrSteiner> node) {
        auto parent = node;
        for (auto child : parent->children) {
            if (parent->layerIdx == child->layerIdx) {
                std::shared_ptr<gr::GrSteiner> lower, upper;
                if ((*parent)[X] < (*child)[X] || (*parent)[Y] < (*child)[Y]) {
                    lower = parent;
                    upper = child;
                } else {
                    lower = child;
                    upper = parent;
                }

                gr::GrEdge edge(gr::GrPoint(lower->layerIdx, lower->x , lower->y ),
                                gr::GrPoint(upper->layerIdx, upper->x, upper->y) );

                         

                DBU edge_cost = grDatabase.getWireCostRelax(edge,wireMap_tpls,viaMap_tpls,false);
                if(edge_cost > edge_highest_cost){
                    edge_highest_cost = edge_cost;
                }
                
            } else {
               
                auto lower_layer = std::min(parent->layerIdx,child->layerIdx);
                auto upper_layer = std::max(parent->layerIdx,child->layerIdx);

                
                    
            }
        }
    };


    for (const std::shared_ptr<gr::GrSteiner>& tree : gridTopo) {
        gr::GrSteiner::postOrder(tree,pathCostFun );
    }

    return edge_highest_cost;


}//end getEdgeHighestCost


void GrNet::checkObidenceGrVsDr(){
    log() << "checkObidenceGrVsDr net: " << getName() << std::endl;
    // if(getName() != "net610") return;
    log() << "dbNet.initRouteDR.size(): " << dbNet.initRouteDR.size() << std::endl;
    log() << "dbNet.initRouteGuides.size(): " << dbNet.initRouteGuides.size() << std::endl;
    // wrong way wiring
    // 0: vertical
    // 1: horizontal
    // check wrong way wiring
    // for(auto box : dbNet.initRouteDR){
    //     // log() << "drBox: " << box << std::endl;
        
    //     int layer_dir = database.getLayerDir(box.layerIdx);
    //     int wire_dir = box.width() < box.height() ? 0 : 1;

    //     if(layer_dir != wire_dir) {
    //         log() << "net: " << getName() << " has wrong way wiring " << std::endl;
    //     }
    // }


    // initRouteGuides
    // for(int i = 0; i < database.nets.size();i++){
    std::vector<RTree> local_guides_rtree; 
    local_guides_rtree.resize(database.getLayerNum());
    int j = 0;
    for(auto guide : dbNet.initRouteGuides){

        boostBox box(boostPoint(guide.lx(),
                                guide.ly()),
                    boostPoint( guide.hx(),
                                guide.hy()));
        int layerIdx = guide.layerIdx;
        local_guides_rtree[layerIdx].insert({box, j});
        j++;
    }//end guide loop 
    
    for(auto drBox : dbNet.initRouteDR){
        vector<std::pair<boostBox, int>> queryResults;
        // log() << "drBox: " << drBox << std::endl;
        boostBox rtreeQueryBox(boostPoint(drBox.lx(),
                                drBox.ly()),
                    boostPoint( drBox.hx(),
                                drBox.hy()));
        auto layer_idx = drBox.layerIdx;
        local_guides_rtree[layer_idx].query(bgi::intersects(rtreeQueryBox), 
                                    std::back_inserter(queryResults));  
        
        bool isInside = false;

        log() << "layer_idx: " << layer_idx 
              << ", queryResults size: " << queryResults.size() << std::endl;

        for (const auto& queryResult : queryResults) {
            auto b = queryResult.first;

            auto guideBox = utils::BoxT<DBU>(bg::get<bg::min_corner, 0>(b),
                                bg::get<bg::min_corner, 1>(b),
                                bg::get<bg::max_corner, 0>(b),
                                bg::get<bg::max_corner, 1>(b));
            
            log() << "guideBox: " << guideBox << std::endl;
            log() << "drBox: " << drBox << std::endl;

            if((drBox.lx() >= guideBox.lx()) &&
                (drBox.ly() >= guideBox.ly()) && 
                (drBox.hx() <= guideBox.hx()) && 
                (drBox.hy() <= guideBox.hy())){
                    isInside = true;
                    break;
                }

            
                

        }//end database.nets.size()  

        if(!isInside){
            log() << "out of guide net: " << dbNet.getName() << std::endl;
            log() << "drBox: " << drBox << std::endl;
            // log() << "guideBox: " << guideBox << std::endl;
            return;
        }
            

    }//end drBox loop
    // }
        

        

        
        

}//end checkObidenceGrVsDr


void GrNet::updateNetFeatures(){
    // bool debug = false;

    // if(debug){
    //     log() << "update Net features" << std::endl;
    //     log() << "net: " << dbNet.getName() << std::endl;
    calcHPWL();
    // }
    // updateHPWLNet();
}//end updateNetFeatures

void GrNet::calcHPWL(){
    bool debug = false;
    std::vector<db::Cell> cells;
    std::vector<int> cells_idx;
    std::vector<utils::BoxT<DBU>> cell_boxs;
    getNetCells(cells_idx);
    for(int cell_idx : cells_idx){
        cells.push_back(database.cells[cell_idx]);
    }
    std::vector<DBU> x_positions;
    std::vector<DBU> y_positions;
    for (auto cell_tmp : cells){
        if(debug) log() << "cell: " << cell_tmp.getName() << std::endl;
        // if(rsynInstance.getName() == cell_tmp.getName()){
        // cell_boxs.push_back(box);
        // }else{
            // auto box_tmp = utils::BoxT<DBU>(
            //     cell_tmp.rsynInstance.getBounds().getLower().x +
            //              rsynInstance.getBounds().getWidth()/2.0,
            //     cell_tmp.rsynInstance.getBounds().getLower().y + 
            //     cell_tmp.rsynInstance.getBounds().getHeight()/2.0,
            //     cell_tmp.rsynInstance.getBounds().getUpper().x,
            //     cell_tmp.rsynInstance.getBounds().getUpper().y
            // );
            // cell_boxs.push_back(box_tmp);
        // }  
                 
        // x_positions.push_back(cell_tmp.rsynInstance.getBounds().getLower().x +
        //                       cell_tmp.rsynInstance.getBounds().getWidth()/2.0);
        // y_positions.push_back(cell_tmp.rsynInstance.getBounds().getLower().y +
        //                       cell_tmp.rsynInstance.getBounds().getHeight()/2.0);

        x_positions.push_back(cell_tmp.rsynInstance.getBounds().getLower().x);
        x_positions.push_back(cell_tmp.rsynInstance.getBounds().getUpper().x);
        y_positions.push_back(cell_tmp.rsynInstance.getBounds().getLower().y);
        y_positions.push_back(cell_tmp.rsynInstance.getBounds().getUpper().y);
                            
        
    }
    // std::vector<DBU> x_positions;
    // std::vector<DBU> y_positions;
    // for(auto cell_box : cell_boxs){
    //     x_positions.push_back(cell_box.lx());
    //     x_positions.push_back(cell_box.hx());
    //     y_positions.push_back(cell_box.ly());
    //     y_positions.push_back(cell_box.hy());
    // }

    DBU hpwl = getHPWLPositions(x_positions,y_positions);

    if(debug){
        log() << "net: " << getName() 
              << ", hpwl: " << hpwl << std::endl;
    }

     if(debug){
        for(auto x : x_positions)
            log() << x << ",";
        log() << std::endl;
        for(auto y : y_positions)
            log() << y << ",";
        log() << std::endl;

    }

    feature_.hpwl_ = hpwl;

    // log() << "hpwl: " << hpwl << std::endl;


}


void GrNet::getNetCells(std::vector<int>& cells_idx){
    for (auto pin : dbNet.rsynPins){
        std::string inst_name = pin.getInstanceName();
        int cell_idx = database.cellsToIdx[inst_name];
        cells_idx.push_back(cell_idx);
    }   
}

}  // namespace gr