#include "Cell.h"

#include <fstream>
#include <random>

#include "Setting.h"
#include "Database.h"
#include "gr_db/GrDatabase.h"
#include "single_net/SingleNetRouter.h"
#include "single_net/InitRoute.h"
#include "single_net/InitRouteCell.h"
#include "flute/flute.h"
#include "place/Legalizer.h"

std::mutex cellFluteMutex;


using namespace std;

extern "C" {
Tree flute(int d, DTYPE x[], DTYPE y[], int acc);
}

namespace db {

struct hash_tuple_cell {  // hash binary tuple
    template <class T>
    size_t operator()(const tuple<T, T>& tup) const {
        auto hash1 = hash<T>{}(get<0>(tup));
        auto hash2 = hash<T>{}(get<1>(tup));
        return hash1 ^ hash2;
    }
};

CellBase::~CellBase() {}



void CellBase::print(ostream& os) const {
    log() << "print cellbase ..." << std::endl;
}

Cell::Cell(int i, Rsyn::Instance instance, RsynService& rsynService,bool debug_t) {
    idx = i;
    rsynInstance = instance;
    cost = 0;
    isFixed = false;
    is_critical_cell = false;
    debug = debug_t;
}

void Cell::initConnectedCells(){
    std::vector<db::Net> nets;
    std::vector<int> nets_idx;
    
    getCellNetsIdx(nets_idx);

    for(auto idx : nets_idx){
        nets.push_back(database.nets[idx]);
    }

    for(auto net : nets){
        std::vector<int> cells_idx;
        getNetCells(net.idx,cells_idx);
        for(int cell_idx : cells_idx){
            if(cell_idx != idx)
                connectd_cells.insert(cell_idx);
        }
    }

}

void Cell::initCandidates(){
    // add current position as placement candidates
    auto box_init_position = 
            utils::BoxT<DBU>(rsynInstance.getBounds().getLower().x,
                                rsynInstance.getBounds().getLower().y,
                                rsynInstance.getBounds().getUpper().x,
                                rsynInstance.getBounds().getUpper().y);

    
    if(!isRedundantCandidate(box_init_position))
        placement_candidates.emplace_back(box_init_position,0);

}//end initCandidates

// check if the placement candidate already is in the 
// placement candidates or not? If yes, this candidate is redundant and 
// will not added to the placement_candidates vector.
bool Cell::isRedundantCandidate(utils::BoxT<DBU> box){
    for(auto cand_pair : placement_candidates){
        auto box_cand = cand_pair.first;
        if(box.isSame(box_cand))
            return true;
    }

    if(!isValidPostion(box)){
        return true;
    }
    // is valid candidate check
    // auto site_min  = database.getRowDBU(box.lx());
    // auto row_min   = database.getSiteDBU(box.ly());
    // auto cell_width = getCellSiteWidth();

    // if(site_min + cell_width > database.sites_num){
    //     return true;
    // }

    // if(site_min < 0 || row_min < 0 ){
    //     return true;
    // }
    // if(row_min > database.rows_num){
    //     return true;
    // }

    return false;
}

// This will be added in the case this cell is 
// overlapping with another cell. 
void Cell::addOvrlpCandidates(utils::BoxT<DBU> box){
    if(isFixed) return;
    initCandidates();
    if(!isRedundantCandidate(box))
        placement_candidates.emplace_back(box,0);
}

// calculate the cost of all new placement candidates.
void Cell::calcCostPlacementCandidates(){
    if(placement_candidates.size() <= 1) return;

    // Temp
    

    for(auto& pl_cad : placement_candidates){
        // log() << "cell: " << getName() 
        //       << ", pl_ca: " << pl_cad.first << std::endl;
        calcCostPlacementCandidate(pl_cad.first);
    }
    // for(int i = 0; i< placement_candidates.size() ; i++){
    //     // if(database.suspected_cells_dict.find(idx) != database.suspected_cells_dict.end())
    //     //     placement_candidates[i].second = -1*cost_eqs[i].getCost();     
    //     // else
    //     // if(i==0)
    //     //     placement_candidates[i].second = (feature.degree_connected_cells *cost_eqs[i].getCost())-1;     
    //     // else
    //         placement_candidates[i].second = feature.degree_connected_cells *cost_eqs[i].getCost();     
    //     // log() << "cell: " << getName() 
    //     //       << ", box: " << placement_candidates[i].first 
    //     //       << ", cost: " << placement_candidates[i].second << std::endl;
    // }
    // // normalizeCostEq();  
    std::vector<std::pair<utils::BoxT<DBU>,double>> valid_placement_candidates;
    for(int i = 0; i< placement_candidates.size() ; i++){
        auto box_init_position = 
        utils::BoxT<DBU>(placement_candidates[i].first.lx(),
                         placement_candidates[i].first.ly(),
                         placement_candidates[i].first.hx(),
                         placement_candidates[i].first.hy());

    
        double cost = feature.degree_connected_cells *cost_eqs[i].getCost();
        // if(!isRedundantCandidate(box_init_position))
        // if(!cost_eqs[i].has_vio)
            valid_placement_candidates.emplace_back(box_init_position,cost);
    }

    placement_candidates.clear();

    for(auto pl : valid_placement_candidates){
        placement_candidates.push_back(pl);
    }
    // normalizeCostEq();  

}

void Cell::normalizeCostEq(){
    double alpha = 1,beta = 1,gamma = 1;
    double init_wl = 0,init_via = 0,init_cong = 0;

    auto costEq_max_wl = *std::max_element(cost_eqs.begin(), cost_eqs.end(),
            [] (CostEq& lhs, CostEq& rhs) {
            return lhs.wl_cost < rhs.wl_cost;
    });
    auto costEq_min_wl = *std::max_element(cost_eqs.begin(), cost_eqs.end(),
            [] (CostEq& lhs, CostEq& rhs) {
            return lhs.wl_cost > rhs.wl_cost;
    });
    auto costEq_max_via = *std::max_element(cost_eqs.begin(), cost_eqs.end(),
            [] (CostEq& lhs, CostEq& rhs) {
            return lhs.via_cost < rhs.via_cost;
    });
    auto costEq_min_via = *std::max_element(cost_eqs.begin(), cost_eqs.end(),
            [] (CostEq& lhs, CostEq& rhs) {
            return lhs.via_cost > rhs.via_cost;
    });
     auto costEq_max_cong = *std::max_element(cost_eqs.begin(), cost_eqs.end(),
            [] (CostEq& lhs, CostEq& rhs) {
            return lhs.cong_cost < rhs.cong_cost;
    });
    auto costEq_min_cong = *std::max_element(cost_eqs.begin(), cost_eqs.end(),
            [] (CostEq& lhs, CostEq& rhs) {
            return lhs.cong_cost > rhs.cong_cost;
    });

    for(int i = 0; i < cost_eqs.size() ; i++){
        auto& cost_eq = cost_eqs[i];
        double wl_tmp = 0, via_tmp = 0, cong_tmp = 0; 
        if(costEq_max_wl.wl_cost != costEq_min_wl.wl_cost)
            wl_tmp = (cost_eq.wl_cost - costEq_min_wl.wl_cost)/ 
                (costEq_max_wl.wl_cost - costEq_min_wl.wl_cost);
        else
            wl_tmp = 1;

        if(costEq_max_via.via_cost != costEq_min_via.via_cost)
            via_tmp = (cost_eq.via_cost - costEq_min_via.via_cost)/ 
                (costEq_max_via.via_cost - costEq_min_via.via_cost);
        else
            via_tmp = 1;

        if(costEq_max_cong.cong_cost != costEq_min_cong.cong_cost)
            cong_tmp = (cost_eq.cong_cost - costEq_min_cong.cong_cost)/ 
                (costEq_max_cong.cong_cost - costEq_min_cong.cong_cost);
        else
            cong_tmp = 1;
        
        placement_candidates[i].second = cost_eq.wl_abs_cost;//wl_tmp;// + via_tmp + cong_tmp;
        // placement_candidates[i].second = cong_tmp;
        //     + cost_eq.via_cost;
        // log() << "cost_wl " << cost_eq.wl_cost << std::endl;
        // if(i == 0){
        //     init_wl = cost_eq.wl_cost;
        //     init_via = cost_eq.via_cost;
        //     init_cong = cost_eq.cong_cost;
        //     placement_candidates[i].second = alpha*init_wl;// + beta + gamma;
        // }else{
        //     placement_candidates[i].second = 
        //         (alpha*(cost_eq.wl_cost-init_wl)/init_wl);
        //         //  +
        //         // (beta*(cost_eq.via_cost-init_via)/init_via);
        //         //  +
        //         // (alpha*(cost_eq.cong_cost-init_cong)/init_cong);
        // }
    }//end for 
}


void Cell::getCellNetsIdx(std::vector<int>& nets_idx){
    auto rsynService = database.getRsynService();
    for (Rsyn::Pin pin : rsynInstance.allPins(false)) {
        auto netRsyn = pin.getNet();
        if(netRsyn){
            int idx = database.netsToIdx[netRsyn.getName()];
            nets_idx.push_back(idx);
        }
    }//end for allpins
}//end getCellNetsIdx

void Cell::getNetCells(int net_idx,std::vector<int>& cells_idx){
    auto net = database.nets[net_idx];
    for (auto pin : net.rsynPins){
        std::string inst_name = pin.getInstanceName();
        int cell_idx = database.cellsToIdx[inst_name];
        cells_idx.push_back(cell_idx);
    }   
}//end getNetCells

void Cell::getPinAccessBoxes(Rsyn::PhysicalPort phPort, vector<BoxOnLayer>& accessBoxes) {
    bool debug = false;
    auto displacement = phPort.getPosition();
    auto bounds = phPort.getBounds();
    Bounds dummyCellBounds(displacement, displacement);
    Rsyn::PhysicalTransform transform(dummyCellBounds, phPort.getOrientation());
    bounds.translate(displacement);
    if(debug){
        log() << phPort.getName() << std::endl;
        // for(auto bd : bounds){
            log() << bounds.getLower().x 
                  << ", " << bounds.getLower().y
                  << ", " << bounds.getUpper().x
                  << ", " << bounds.getUpper().y
                  << std::endl;
        // }
            
    }

    bounds = transform.apply(bounds);



    accessBoxes.emplace_back(phPort.getLayer().getRelativeIndex(), getBoxFromRsynBounds(bounds));
}// getPinAccessBoxes

void Cell::getPinAccessBoxes(Rsyn::PhysicalLibraryPin phLibPin,
                            Rsyn::PhysicalCell phCell,
                            vector<BoxOnLayer>& accessBoxes,
                            utils::BoxT<DBU>& new_position,
                            const DBUxy& origin) {

    bool debug = false;
    if (!phLibPin.hasPinGeometries()) {
        log() << "Warning: pin of " << phCell.getName() << " has no pinGeometries" << std::endl;
        return;
    }
    auto rsynService = database.getRsynService();
    
    auto orientation_mov = getMoveOrientation(new_position);
    const DBUxy displacement = phCell.getPosition() + origin;
    const DBUxy displacement_mov = DBUxy({new_position.lx(),new_position.ly()}) + origin; 

    auto transform = phCell.getTransform();
    Bounds rsynBds(displacement_mov.x,displacement_mov.y,
                  displacement_mov.x + new_position.width(),
                  displacement_mov.y + new_position.height());
    Rsyn::PhysicalTransform transform_mov(rsynBds, orientation_mov);
    

    // for (Rsyn::PhysicalPinGeometry phPinGeo : phLibPin.allPinGeometries()) {
    // TODO: check why multiple PinGeometry on 8t4 inst60849
    auto phPinGeo = phLibPin.allPinGeometries()[0];
    // if(debug){
    //     log() << phLibPin.getName() << std::endl;
    // }
    for (Rsyn::PhysicalPinLayer phPinLayer : phPinGeo.allPinLayers()) {
        if (!phPinLayer.hasRectangleBounds()) {
            log() << "Warning: pin has no RectangleBounds" << std::endl;
            continue;
        }
        int layerIdx = phPinLayer.getLayer().getRelativeIndex();
        for (auto bounds : phPinLayer.allBounds()) {
            // log() << "bds: " << bounds << std::endl;
            if(getName() == phCell.getInstance().getName()){
                // log() << "displacement_mov: " << displacement_mov << std::endl;
                bounds.translate(displacement_mov);
                bounds = transform_mov.apply(bounds);
                // log() << "final_bounds_mov: " << bounds << std::endl;
            }else{
                // log() << "displacement: " << displacement << std::endl;
                bounds.translate(displacement);
                bounds = transform.apply(bounds);
                // log() << "final_bounds: " << bounds << std::endl;
            }
            if(debug){
                
                // for(auto bd : bounds){
                    log() << bounds.getLower().x 
                        << ", " << bounds.getLower().y
                        << ", " << bounds.getUpper().x
                        << ", " << bounds.getUpper().y
                        << std::endl;
                // }
                    
            }
            accessBoxes.emplace_back(layerIdx, getBoxFromRsynBounds(bounds));
        }
    }
    // for(auto ac : accessBoxes){
    //     log() << "ac: " << ac << std::endl;
    // }
}

void Cell::initPinAccessBoxes(Rsyn::Pin rsynPin,
                             RsynService& rsynService,
                             vector<BoxOnLayer>& accessBoxes,
                             utils::BoxT<DBU>& new_position,
                             const DBU libDBU){
    bool debug = false;
    if(debug)
        log() << "rsynPin: " << rsynPin.getName() << std::endl;
    if (rsynPin.isPort()) {
        if(debug){
            log() << "is port..." << std::endl;
        }
        Rsyn::PhysicalPort phPort = rsynService.physicalDesign.getPhysicalPort(rsynPin.getPort());
        getPinAccessBoxes(phPort, accessBoxes);
        return;
    }

    // PhysicalLibraryPin
    Rsyn::PhysicalLibraryPin phLibPin = rsynService.physicalDesign.getPhysicalLibraryPin(rsynPin);

    // PhysicalCell
    Rsyn::Instance instance = rsynPin.getInstance();
    if (instance.getType() != Rsyn::CELL) {
        log() << "Warning: pin is not on a cell " << rsynPin.getNetName() << " " << rsynPin.getInstanceName()
              << std::endl;
        return;
    }
    Rsyn::Cell cell = instance.asCell();
    Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cell);
    Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cell);
    const DBUxy origin(static_cast<DBU>(std::round(phLibCell.getMacro()->originX() * libDBU)),
                       static_cast<DBU>(std::round(phLibCell.getMacro()->originY() * libDBU)));
    // fill accessBoxes
    getPinAccessBoxes(phLibPin, phCell, accessBoxes, new_position, origin);
}//end getPinAccessBoxes

void Cell::calcCostPlacementCandidate(utils::BoxT<DBU>& box){
    bool debug = false;
    auto rsynService = database.getRsynService();
    const Rsyn::Session session;
    const Rsyn::PhysicalDesign& physicalDesign =
        static_cast<Rsyn::PhysicalService*>(session.getService("rsyn.physical"))->getPhysicalDesign();
    const DBU libDBU = physicalDesign.getDatabaseUnits(Rsyn::LIBRARY_DBU);
    
    std::vector<db::Net> nets;
    std::vector<int> nets_idx;
    getCellNetsIdx(nets_idx);
    for(auto idx : nets_idx)
        nets.push_back(database.nets[idx]);

    CostEq cost_eq;

    auto old_box = utils::BoxT<DBU>(
                rsynInstance.getBounds().getLower().x,
                rsynInstance.getBounds().getLower().y,
                rsynInstance.getBounds().getUpper().x,
                rsynInstance.getBounds().getUpper().y
            );

    bool isSameBox = old_box.isSame(box);
    // debug=true;
    // if(debug)
    //     log() << "cell: " << getName() << ", new_box: " << box << std::endl;
    // debug=false;

    for(auto net : nets){
        if(net.getName() == "net26604"){
            debug = true;
        }else{
            debug = false;
        }


        if(debug){
            log() << "cell: " << getName() << std::endl;
            log() << "net: " << net.getName() << std::endl;
        }
        //----------------PinAccessBoxes--------------
        vector<vector<BoxOnLayer>> pinAccessBoxes; 
        if(database.policy_set.find("updatePinAccessBoxes") == database.policy_set.end()){
            int tgPin_idx = updatePinAccessBoxes(net.idx,box,pinAccessBoxes);
        }
        std::vector<utils::BoxT<DBU>> cell_boxs;
        if(database.policy_set.find("updateCellBoxes") == database.policy_set.end()){
            updateCellBoxes(net.idx,box,cell_boxs);
        }
        //--------------------------------------------
        // calcAlignmentCost(pinAccessBoxes,tgPin_idx,cost_eq);
        // calcAlignmentSiteBasedCost(pinAccessBoxes,tgPin_idx,cost_eq);
        //--------------------------------------------
        // if(!isSameBox){

        // if(database.policy_set.find("hybridCostEstimation") == database.policy_set.end()){
        //     if(net.numOfPins() > 5){
                calcFluteRouting3DCost(net.idx,pinAccessBoxes,cell_boxs,cost_eq);
        //     }else{
        //         calcHPWLCost(cell_boxs,cost_eq);
        //     }
        // }else{
        //     if(database.policy_set.find("calcFluteRouting3DCost") == database.policy_set.end()){
        //         calcFluteRouting3DCost(net.idx,pinAccessBoxes,cell_boxs,cost_eq);
        //     }
        //     if(database.policy_set.find("calcHPWLCost") == database.policy_set.end()){
        //         calcHPWLCost(cell_boxs,cost_eq);
        //     }

        // }

        if(debug){
            log() << "cost_eq: " << cost_eq.getCost() << std::endl;
        }

        
        
        // }else{
        //     cost_eq.path_cost = grDatabase.nets[net.idx].getPathCost();
        //     cost_eq.wl_abs_cost = grDatabase.nets[net.idx].getWirelength();
        // }
            
        // //--------------------------------------------
        // calcHPWLCost(cell_boxs,cost_eq);
    }//end for 
    
    cost_eqs.push_back(cost_eq);
    
}//end calcCostPlacementCandidate

void Cell::calcHPWLCost(
    std::vector<utils::BoxT<DBU>>& cell_boxs
    ,CostEq& cost_eq){
    std::vector<DBU> x_positions;
    std::vector<DBU> y_positions;

    bool debug = false;

    // if(debug) {
    //     log() << "calcHPWLCost ..." << std::endl;
    //     log() << "cell name: " << getName() << std::endl;
    //     for(auto tmp : cell_boxs){
    //         log() << tmp << std::endl;
    //     }
    // }

    for(auto cell_box : cell_boxs){
        x_positions.push_back(cell_box.lx());
        x_positions.push_back(cell_box.hx());
        y_positions.push_back(cell_box.ly());
        y_positions.push_back(cell_box.hy());
    }

    if(debug){
        for(auto x : x_positions)
            log() << x << ",";
        log() << std::endl;
        for(auto y : y_positions)
            log() << y << ",";
        log() << std::endl;

    }

    DBU hpwl = getHPWLPositions(x_positions,y_positions);
    if(debug)
        log() << "hpwl: " << hpwl << std::endl;
    cost_eq.updateHPWLCost(hpwl);
    cost_eq.updateHybridCost(hpwl);
}//end calcHPWLCost

void Cell::calcFluteRouting3DCost(int net_idx
    ,vector<vector<BoxOnLayer>>& pinAccessBoxes
    ,std::vector<utils::BoxT<DBU>>& cell_boxs
    ,CostEq& cost_eq){
    bool debug_tmp = false;
    // if(grDatabase.nets[net_idx].getName() == "net3144") debug_tmp = true;
    // if(grDatabase.nets[net_idx].getName() == "net2564") debug_tmp = true;

    bool relax_mode = true;
    //temp
    // grDatabase.removeNet(grDatabase.nets[net_idx]);
    // grDatabase.nets[net_idx].gridTopo.clear();
    // allNetStatus[net_idx] = db::RouteStatus::FAIL_UNPROCESSED;
    //end temp
    routeCell::InitRouteCell router_estimator(grDatabase.nets[net_idx],cell_boxs,pinAccessBoxes,relax_mode,debug_tmp);
    router_estimator.plan_fluteOnly();
    router_estimator.getRoutingOrder();
    if(db::setting.patternRouteMemorization){
        router_estimator.patternRouteMemo();
    }else{
        router_estimator.patternRoute();
    }
    bool has_vio = router_estimator.buildTopo();
    debug_tmp = false;
    if(debug_tmp){
        log() << "grDatabase.nets[net_idx]: " << grDatabase.nets[net_idx].getName() << std::endl;
        // log() << "wl_cost: " << router_estimator.wl_cost << std::endl;
        // log() << "wl_abs_cost: " << router_estimator.wl_abs_cost  << std::endl;
        // log() << "path_cost: " << router_estimator.path_cost << std::endl;
        // log() << "path_cost_grNet: " << grDatabase.nets[net_idx].getPathCost()<< std::endl;
        log() << "has_vio: " << has_vio << std::endl;
    }
    
    int penalty = 1;
    // if(has_vio)
    //     penalty = 1000000;
    
    cost_eq.updateCongestionCost(penalty*router_estimator.cong_cost);
    cost_eq.updateWireLengthAbsCost(penalty*router_estimator.wl_abs_cost);
    cost_eq.updateWireLengthCost(penalty*router_estimator.wl_cost);
    cost_eq.updateViaAbsCost(penalty*router_estimator.via_abs_cost);
    cost_eq.updateViaCost(penalty*router_estimator.via_cost);
    cost_eq.updatePathCost(penalty*router_estimator.path_cost);
    // hybrid cost
    cost_eq.updateHybridCost(penalty*router_estimator.wl_abs_cost);
    cost_eq.updateHasVio(has_vio);
}//end 

void Cell::calcAlignmentCost(vector<vector<BoxOnLayer>>& pinAccessBoxes,
            int tgPin_idx,CostEq& cost_eq){
    DBU overlap_pinArea = getPinIoU(pinAccessBoxes,tgPin_idx);
    DBU total_pinArea = getPinTotalArea(pinAccessBoxes,tgPin_idx);
    auto alignment_cost = total_pinArea-overlap_pinArea;
    cost_eq.updateAlignmentCost(alignment_cost);
}//end calcAlignmentCost

void Cell::calcAlignmentSiteBasedCost(vector<vector<BoxOnLayer>>& pinAccessBoxes,int tgPin_idx,CostEq& cost_eq){
    // How many pins the target cells has overlap with
    // log() << "calcAlignmentSiteBasedCost cell: " << getName() << std::endl;
    // log() << "pinAccessBoxes size: " << pinAccessBoxes << std::endl;
    int num_ovrlp_pins = 0;
    float alignment_cost = 0;
    std::set<int> pin_tg_sites;
    int j =0;
    for (auto boxs : pinAccessBoxes) {
        if (j!= tgPin_idx) {
            j++;
            continue;
        }
        for(auto box : boxs){
            auto tmpbox = utils::BoxT<DBU>(box.lx(), box.ly(),box.hx(), box.hy());  
            std::vector<int> tmp_sites;
            // log() << "tmpbox: " << tmpbox << std::endl;
            getIntesectedSites(tmpbox,tmp_sites);
            for(auto site : tmp_sites){
                // log() << "site: " << site << std::endl;
                pin_tg_sites.insert(site);
            }
                
        }
        break;
    }

    // for(auto tmp : pin_tg_sites){
    //     log() << "pin_tg_sites: " << tmp << std::endl;
    // }

    int i = 0;
    for(auto boxs : pinAccessBoxes){
        if(i == tgPin_idx) {
            i++;
            continue;
        }
        std::set<int> pin_sites;
        for(auto box : boxs){
            auto tmpbox = utils::BoxT<DBU>(box.lx(), box.ly(),box.hx(), box.hy());  
            std::vector<int> tmp_sites;
            getIntesectedSites(tmpbox,tmp_sites);
            for(auto site : tmp_sites)
                pin_sites.insert(site);
        }//end for boxs
        // for(auto tmp : pin_sites){
        //     log() << "pin_sites: " << tmp << std::endl;
        // }

        bool is_algined = false;
        int num_aligned = 0;
        for(auto site_tg : pin_tg_sites){
            if(pin_sites.find(site_tg) != pin_sites.end()){
                is_algined = true;
                num_aligned++;
            }
        }
        int num_total_sites_pin = pin_sites.size();
        float alignment_percentage = 0;
        if(num_total_sites_pin != 0)
            alignment_percentage = float(num_aligned)/float(num_total_sites_pin);
        // log() << "num_aligned: " << num_aligned << std::endl;
        // log() << "num_total_sites_pin: " << num_total_sites_pin << std::endl;
        // log() << "alignment_percentage: " << alignment_percentage << std::endl;

        alignment_cost += alignment_percentage;

    }//end for pinAccessBoxs

    // log() << "alignment_cost: " << alignment_cost << std::endl;

    cost_eq.updateAlignmentCost(alignment_cost*100);

    
}//end calcAlignmentSiteBasedCost

int Cell::updatePinAccessBoxes(int net_idx,utils::BoxT<DBU>& box,vector<vector<BoxOnLayer>>& pinAccessBoxes){
    auto rsynService = database.getRsynService();
    const Rsyn::Session session;
    const Rsyn::PhysicalDesign& physicalDesign =
        static_cast<Rsyn::PhysicalService*>(session.getService("rsyn.physical"))->getPhysicalDesign();
    const DBU libDBU = physicalDesign.getDatabaseUnits(Rsyn::LIBRARY_DBU);
    int tgPin_idx = 0;
    int tmpIdx = 0;
    auto net = database.nets[net_idx];

    bool debug = false;

    // if(net.getName() == "net121")
    //     debug = true;

    for (auto rsynPin : net.rsynNet.allPins()) {
        //rsynPin, rsynService
        //  rsynPins.push_back(RsynPin);
        pinAccessBoxes.emplace_back();
        // if(RsynPin.getInstance().getName() != getName())
        // if(debug)
        // log() << "inst pins: " << rsynPin.getInstance().getName() << std::endl;
        // log() << "initPinAccessBoxes..." << std::endl;
        initPinAccessBoxes(rsynPin, rsynService, pinAccessBoxes.back(), box,libDBU);

        

        if( rsynPin.getInstance().getName() == getName()){
            tgPin_idx = tmpIdx;
        }
        tmpIdx++;
    }

    if(debug){
        for(auto vec_pins : pinAccessBoxes){
            log() << "new pin"<< std::endl;
            for(auto pin_box : vec_pins){
                log() << "pin_box: " << pin_box << std::endl;
            }
        }

        log() << "tgPin_idx: " << tgPin_idx << std::endl;
    }

    return tgPin_idx;
}//end updatePinAccessBoxes

void Cell::updateCellBoxes(int net_idx, utils::BoxT<DBU>& box, 
          std::vector<utils::BoxT<DBU>>& cell_boxs){

    auto net = database.nets[net_idx];
    std::vector<int> cells_idx;
    std::vector<db::Cell> cells;
    getNetCells(net.idx,cells_idx);
    for(int cell_idx : cells_idx){
        cells.push_back(database.cells[cell_idx]);
    }
    for (auto cell_tmp : cells){
        if(rsynInstance.getName() == cell_tmp.getName()){
            cell_boxs.push_back(box);
        }else{
            auto box_tmp = utils::BoxT<DBU>(
                cell_tmp.rsynInstance.getBounds().getLower().x,
                cell_tmp.rsynInstance.getBounds().getLower().y,
                cell_tmp.rsynInstance.getBounds().getUpper().x,
                cell_tmp.rsynInstance.getBounds().getUpper().y
            );
            cell_boxs.push_back(box_tmp);
        }           
    }
    
}//end updateCellBoxes

int Cell::getHPWLPositions(std::vector<DBU> xs,std::vector<DBU> ys){
    int max_x = *std::max_element(xs.begin(), xs.end());
    int min_x = *std::min_element(xs.begin(), xs.end());
    int max_y = *std::max_element(ys.begin(), ys.end());
    int min_y = *std::min_element(ys.begin(), ys.end());

    return 2.5*std::abs(max_x-min_x) + std::abs(max_y-min_y);
}
void Cell::getConnectedCells(std::vector<int>& cells_idx){
    std::vector<db::Net> nets;
    std::vector<int> nets_idx;

    getCellNetsIdx(nets_idx);


    for(auto idx : nets_idx)
        nets.push_back(database.nets[idx]);

    for(auto net : nets){
        std::vector<int> cells_idx_tmp;        
        
        getNetCells(net.idx,cells_idx_tmp);
        for(auto tmp : cells_idx_tmp){
            cells_idx.push_back(tmp);
        }
    }

}

void Cell::updatePinBoundingBox(){
    auto rsynService = database.getRsynService();
    const Rsyn::Session session;
    const Rsyn::PhysicalDesign& physicalDesign =
    static_cast<Rsyn::PhysicalService*>(session.getService("rsyn.physical"))->getPhysicalDesign();
    const DBU libDBU = physicalDesign.getDatabaseUnits(Rsyn::LIBRARY_DBU);
    auto dieBound = rsynService.physicalDesign.getPhysicalDie().getBounds();

    std::vector<db::Net> nets;
    std::vector<int> nets_idx;

    getCellNetsIdx(nets_idx);


    for(auto idx : nets_idx)
        nets.push_back(database.nets[idx]);

    for(auto net : nets){
        std::vector<int> cells_idx;        
        
        getNetCells(net.idx,cells_idx);

        std::vector<Rsyn::Pin> rsynPins;
        for (auto rsynPin : net.rsynNet.allPins()) {
            Rsyn::Instance instance = rsynPin.getInstance();

            rsynPins.push_back(rsynPin);
            
            Rsyn::PhysicalLibraryPin phLibPin = rsynService.physicalDesign.getPhysicalLibraryPin(rsynPin);

            if (instance.getType() != Rsyn::CELL) {
                log() << "Warning: pin is not on a cell " << rsynPin.getNetName() << " " << rsynPin.getInstanceName()
                    << std::endl;
                return;
            }
            Rsyn::Cell cell = instance.asCell();
            Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cell);
            Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cell);
            const DBUxy origin(static_cast<DBU>(std::round(phLibCell.getMacro()->originX() * libDBU)),
                               static_cast<DBU>(std::round(phLibCell.getMacro()->originY() * libDBU)));
        
            if (!phLibPin.hasPinGeometries()) {
                log() << "Warning: pin of " << phCell.getName() << " has no pinGeometries" << std::endl;
                return;
            }
            const DBUxy displacement = phCell.getPosition() + origin;

            auto transform = phCell.getTransform();
    
            // for (Rsyn::PhysicalPinGeometry phPinGeo : phLibPin.allPinGeometries()) {
            // TODO: check why multiple PinGeometry on 8t4 inst60849
            auto phPinGeo = phLibPin.allPinGeometries()[0];
            int i = 0;
            for (Rsyn::PhysicalPinLayer phPinLayer : phPinGeo.allPinLayers()) {
                if (!phPinLayer.hasRectangleBounds()) {
                    log() << "Warning: pin has no RectangleBounds" << std::endl;
                    continue;
                }
                int layerIdx = phPinLayer.getLayer().getRelativeIndex();
                // log() << "i: " << i << std::endl;
                std::vector<DBU> xs;
                std::vector<DBU> ys;
                for (auto bounds : phPinLayer.allBounds()) {
                    // log() << "before translate bounds: " << bounds << std::endl;
                    bounds.translate(displacement);
                    bounds = transform.apply(bounds);                    
                    xs.push_back(bounds.getLower().x);
                    xs.push_back(bounds.getUpper().x);
                    ys.push_back(bounds.getLower().y);
                    ys.push_back(bounds.getUpper().y);
                }
                int max_x = *std::max_element(xs.begin(), xs.end());
                int min_x = *std::min_element(xs.begin(), xs.end());
                int max_y = *std::max_element(ys.begin(), ys.end());
                int min_y = *std::min_element(ys.begin(), ys.end());

                utils::BoxT<DBU> bd_pin(min_x,min_y,max_x,max_y);
                std::string pin_name = instance.getName();// + "_" + rsynPin.getName(); 
                pin_bd_dict[pin_name] = bd_pin;
            }
        }

    }//end nets

}//end updatePinBoundingBox

void Cell::addAlignmentCandidates(std::vector<utils::BoxT<DBU>>& placement_candidates_alignment){
    auto rsynService = database.getRsynService();
    auto dieBound = rsynService.physicalDesign.getPhysicalDie().getBounds();
    Rsyn::Cell cellRsyn = rsynInstance.asCell();
    Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cellRsyn);
    Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cellRsyn);
    int cell_width_movable = phLibCell.getWidth();

    auto cell_lx_movable = rsynInstance.getBounds().getLower().x;
    auto cell_ly_movable = rsynInstance.getBounds().getLower().y;
    auto cell_hx_movable = rsynInstance.getBounds().getUpper().x;
    auto cell_hy_movable = rsynInstance.getBounds().getUpper().y;

    utils::BoxT<DBU> cell_box_movable(cell_lx_movable,cell_ly_movable,
                                      cell_hx_movable,cell_hy_movable);
    int cell_min_site_movable = getSiteDBU(cell_lx_movable);
    int cell_max_site_movable = getSiteDBU(cell_hx_movable);
    int cell_width_site_movable = cell_max_site_movable - cell_min_site_movable;
    

    utils::BoxT<DBU> pin_box_movable;
    for(auto pair : pin_bd_dict){
        if(pair.first == getName()){
            auto pin_bd = pair.second;
            pin_box_movable = pin_bd;
        }
    }



    // DBU cur_pin_width = pin_box_movable.width();

    // int cell_height = phLibCell.getHeight();
    DBU  eps = 10;
    // log() << "all pin_bd_box cell: " << getName() << std::endl;

    // log() << "cell_min_site_movable: " << cell_min_site_movable << std::endl;
    // log() << "cell_max_site_movable: " << cell_max_site_movable << std::endl;
    // log() << "database.sites_num: " << database.sites_num << std::endl;

    std::vector<int> pin_sites_movable;
    getIntesectedSites(pin_box_movable,pin_sites_movable);
    int pin_min_site_movable = *std::min_element(pin_sites_movable.begin(),
                                             pin_sites_movable.end());
    int cell_pin_dist_site_movable = pin_min_site_movable - cell_min_site_movable;
    
    // for(auto tmp : pin_sites_movable){
    //     log() << "pin_sites_movable: " << tmp << std::endl;
    // }

    for(auto pair : pin_bd_dict){
        if(pair.first == getName()) continue;

        std::vector<int> pin_sites;
        auto pin_bd = pair.second;
        getIntesectedSites(pin_bd,pin_sites);
        // for(auto tmp : pin_sites){
        //     log() << "pin_sites: " << tmp << std::endl;
        // }
        int pin_max_site = *std::max_element(pin_sites.begin(), pin_sites.end());
        int pin_min_site = *std::min_element(pin_sites.begin(), pin_sites.end());

        for(int i = pin_min_site - (pin_sites_movable.size()-1); i <= pin_max_site;i++){
            int new_cell_site_lx = i-cell_pin_dist_site_movable;
            int new_cell_site_hx = new_cell_site_lx + cell_width_site_movable - 1; 
            if(new_cell_site_lx < 0 || new_cell_site_hx > database.sites_num-1) continue;
            // log() << "new_cell_site_lx: " << new_cell_site_lx << std::endl;
            auto new_cell_xmin = getDBUSite(new_cell_site_lx);
            // log() << "new_cell_xmin: " << new_cell_xmin/database.libDBU << std::endl;
            // log() << "new_cell_xmax: " << (new_cell_xmin+cell_width_movable)/database.libDBU << std::endl;
            
            placement_candidates_alignment.emplace_back(new_cell_xmin, cell_ly_movable , 
                                                        new_cell_xmin+cell_width_movable, cell_hy_movable);
        }//end for 

        // DBU pin_bd_width = pin_bd.width();
        // DBU max_step = 0;
        // DBU step_size = 0;
        // log() << "pin_bd_width: " << pin_bd_width 
        //       << ", cur_pin_width: " << cur_pin_width << std::endl;
        // if(pin_bd_width > cur_pin_width ){
        //     max_step = pin_bd_width/cur_pin_width;
        //     step_size = cur_pin_width/2;
        // }else{
        //     max_step = cur_pin_width/pin_bd_width;
        //     step_size = pin_bd_width/2;
        // }
        // log() << "step_size: " << step_size << std::endl;
        // log() << "max_step: " << max_step << std::endl;
        // log() << "pin_bd: " << pin_bd << std::endl;
        // log() << "cur_pin_box: " << cur_pin_box << std::endl;
        // log() << "cur_lx_cell: " << cur_lx_cell/database.libDBU << std::endl;

        // auto lx_pin = pin_bd.lx();
        // auto ly_pin = pin_bd.ly();
        // auto hx_pin = pin_bd.hx();
        // auto hy_pin = pin_bd.hy();
        // auto rows_fixCell = getRowsNear(lx_pin,ly_pin);    
        // // log() << "rows_fixCell size: " << rows_fixCell.size() << std::endl;
        // for(auto rowBox_fixCell : rows_fixCell){
        //     // log() << "rowBox_fixCell: " << rowBox_fixCell << std::endl;
        //     float row_height = std::abs(rowBox_fixCell.hy()-rowBox_fixCell.ly());
        //     for(int i = 1; i <= number_of_alignment_rows; i++ ){
        //         auto  rows_fixCell_neigh = getRowsNear(lx_pin,ly_pin-i*row_height);    
        //         // log() << "rows_fixCell_neigh size: " << rows_fixCell_neigh.size() << std::endl;
        //         for(auto rowBox_fixCell_neigh : rows_fixCell_neigh ){
        //             // log() << "rowBox_fixCell_neigh: " << rowBox_fixCell_neigh << std::endl;
        //             DBU step = lx_pin;
        //             for(int j = 0 ; j <= max_step; j++){
        //                 log() << "lx_pin: " << lx_pin/database.libDBU << std::endl;
        //                 step = lx_pin + j*step_size;
        //                 log() << "step: " << step/database.libDBU << std::endl;
        //                 DBU disp = 0;
        //                 DBU xmin = 0;
        //                 if(lx_pin <= cur_pin_box.lx()) {
        //                     disp = cur_pin_box.lx() - step;
        //                     xmin = cur_lx_cell - disp;
        //                     log() << "cur_lx_cell: " << cur_lx_cell
        //                           << "disp: " << disp
        //                           << "xmin: " << xmin << std::endl;
        //                     if(xmin < dieBound.getLower().x){
        //                         xmin = dieBound.getLower().x;
        //                     }
        //                 }else{
        //                     disp = cur_pin_box.lx() - step;
        //                     xmin = cur_lx_cell - disp;
        //                     log() << "cur_lx_cell: " << cur_lx_cell
        //                           << "disp: " << disp
        //                           << "xmin: " << xmin << std::endl;
        //                     if(xmin > dieBound.getUpper().x){
        //                         xmin = dieBound.getUpper().x-cell_width;
        //                     }
        //                 }
        //                 log() << "disp: " << disp/database.libDBU << std::endl;
        //                 log() << "xmin: " << xmin/database.libDBU << std::endl;
                        
        //                 placement_candidates_alignment.emplace_back(xmin, rowBox_fixCell_neigh.ly() , 
        //                                                     xmin+cell_width, rowBox_fixCell_neigh.hy());
        //             }
                    
        //             // placement_candidates_alignment.emplace_back(lx+(cell_width/2), rowBox_fixCell_neigh.ly() , 
        //             //                                         lx+(cell_width/2)+cell_width, rowBox_fixCell_neigh.hy());
        //         }//end rows box neighbour fixed row box
        //     }
        // }//end rows loop
    }// end pin bounding box loop
}//end add alignment candidates

void Cell::addRandomNetBoxCellCandidates(){
    int cell_width = rsynInstance.getBounds().getWidth();
    std::vector<utils::BoxT<DBU>> placement_candidates_rnd;
    //---------------find the net box---------------------
    utils::BoxT<DBU> net_box;
    getNetsBox(net_box);
    std::vector<db::Cell> cells_in_netBox;
    database.getCellsInBox(net_box,cells_in_netBox,10);
    //----------------------------------------------------


    // add random positions
    for(int i = 0; i < 5; i++){
        auto rand = getRandom(net_box);
        auto rand_rows = getRowsNear(rand.first,rand.second);
        for (auto row_box : rand_rows){
            placement_candidates_rnd.emplace_back(rand.first,
                row_box.ly() , rand.first+cell_width, row_box.hy());
        }
    }

    for(auto tmp : placement_candidates_rnd){
        collectCandidatePositions(tmp);
    }

}//end addRandomCellCandidates

void Cell::addEqualWidthNetBoxCellCandidates(){   
    int cell_width = rsynInstance.getBounds().getWidth();
    std::vector<utils::BoxT<DBU>> placement_candidates_rnd;
    //---------------find the net box---------------------
    utils::BoxT<DBU> net_box;
    getNetsBox(net_box);
    std::vector<db::Cell> cells_in_netBox;
    database.getCellsInBox(net_box,cells_in_netBox,10);
    //----------------------------------------------------


    // add random positions
    for(auto cell_tmp : cells_in_netBox){
        if(cell_tmp.rsynInstance.getBounds().getWidth() == cell_width){
            auto rand_rows = getRowsNear(cell_tmp.rsynInstance.getBounds().getLower().x,
                                         cell_tmp.rsynInstance.getBounds().getLower().y);
            for (auto row_box : rand_rows){
                placement_candidates_rnd.emplace_back(cell_tmp.rsynInstance.getBounds().getLower().x,
                    row_box.ly() , cell_tmp.rsynInstance.getBounds().getLower().x+cell_width,
                    row_box.hy());
            }//end for 
        }//end if 
    }//end for 

    for(auto tmp : placement_candidates_rnd){
        collectCandidatePositions(tmp);
    }


}//end addEqualWidthNetBoxCellCandidates
void Cell::addWireDependentCellCandidates(){
    std::vector<int> nets_idx;
    getCellNetsIdx(nets_idx);

    for(auto net_idx : nets_idx){
        if(grDatabase.nets[net_idx].getName() != "net1237") continue;
        log() << "grDatabase.nets[net_idx]: " << grDatabase.nets[net_idx].getName() << std::endl;
        // grDatabase.removeNetRelax(grDatabase.nets[net_idx]);
        // grDatabase.removeNet(grDatabase.nets[net_idx]);
        grDatabase.nets[net_idx].getRouteSites();
    }

}//end addWireDependentCellCandidates

void Cell::addPinAlignmentCellCandidates(){  
    if(isFixed) return;
    if(!is_critical_cell) return; 

    updatePinBoundingBox();

    std::vector<utils::BoxT<DBU>> placement_candidates_alignment;
    addAlignmentCandidates(placement_candidates_alignment);

    for(auto tmp : placement_candidates_alignment){
        collectCandidatePositions(tmp);
    }

    printPlacementCandidates();
}// addPinAlignmentCellCandidates


void Cell::addManualCellCandidates(){
    if(isFixed) return;
    if(!is_critical_cell) return; 
    auto rsynService = database.getRsynService();
    Rsyn::Cell cellRsyn = rsynInstance.asCell();
    Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cellRsyn);
    Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cellRsyn);
    int cell_width = phLibCell.getWidth();
    bool manual = true;
    DBU new_x = 41.8*database.libDBU;
    DBU new_y = 41.04*database.libDBU;

    if(manual){
        if(getName() == "inst4382"){
            std::vector<utils::BoxT<DBU>> placement_candidates_manual;
            auto bd_ly = rsynInstance.getBounds().getLower().y;
            auto bd_hy = rsynInstance.getBounds().getUpper().y;
            auto height = bd_hy-bd_ly;
            placement_candidates_manual.emplace_back(new_x, new_y
                                        , new_x+cell_width, new_y+height);
            for(auto tmp : placement_candidates_manual){
                collectCandidatePositions(tmp);
            }
        }
        

    }else{
        std::vector<int> cells_idx;
        getConnectedCells(cells_idx);
        for(auto cell_idx : cells_idx){
            auto cell = database.cells[cell_idx];
            if(cell.getName() == getName() ) continue;

            
            auto bds = cell.rsynInstance.getBounds();

            auto rows = getRowsNear(bds.getLower().x,bds.getLower().y);

            std::vector<utils::BoxT<DBU>> placement_candidates_manual;

            for (auto row_box : rows){
                placement_candidates_manual.emplace_back(bds.getUpper().x, row_box.ly()
                                        , bds.getUpper().x+cell_width, row_box.hy());
            }

            for(auto tmp : placement_candidates_manual){
                collectCandidatePositions(tmp);
            }
        }

    }



}//end addManualCellCandidates

void Cell::addMedianCellCandidates(){
    if(isFixed) return;
    if(!is_critical_cell) return;

    // log() << "add Median cell: " << getName() << std::endl;

    auto rsynService = database.getRsynService();

    auto median = getMedianCell();
    auto median_x = median.first;
    auto median_y = median.second;
    
    Rsyn::Cell cellRsyn = rsynInstance.asCell();
    Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cellRsyn);
    Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cellRsyn);
    int cell_width = phLibCell.getWidth();
    int cell_height = phLibCell.getHeight();
    

    // log() << "median_x: " << double(median_x)/database.libDBU << ", median_y: " << double(median_y)/database.libDBU << std::endl;
    auto rows = getRowsNear(median_x,median_y);
    

    std::vector<utils::BoxT<DBU>> placement_candidates_median;
    // log() << "row and sites for cell: " << getName() << std::endl;
    for (auto row_box : rows){
        // utils::BoxT<DBU> median_box(median_x, row_box.ly() , median_x+cell_width, row_box.hy());
        // std::vector<int> sites;
        // database.getIntesectedSites(median_box,sites);
            
        // log() << "row: " << row_box << std::endl;
        // for(auto site : sites){
        //     log() << "site: " << site 
        //           << ", dbu: " <<double(getDBUSite(site))/database.libDBU << std::endl;
        // }
        
        placement_candidates_median.emplace_back(median_x, row_box.ly() , median_x+cell_width, row_box.hy());
        placement_candidates_median.emplace_back(
            rsynInstance.getBounds().getLower().x, row_box.ly(),
            rsynInstance.getBounds().getUpper().x, row_box.hy());
        placement_candidates_median.emplace_back(
            median_x, rsynInstance.getBounds().getLower().y,
            median_x+cell_width, rsynInstance.getBounds().getUpper().y);
    }

    // add random positions
    // for(int i = 0; i < 100; i++){
    //     auto rand = getRandom();
    //     auto rand_rows = getRowsNear(rand.first,rand.second);
    //     for (auto row_box : rand_rows){
    //         // log() << "ovrlp rows: " << row_box << std::endl;
    //         placement_candidates_median.emplace_back(rand.first,
    //             row_box.ly() , rand.first+cell_width, row_box.hy());
    //     }
    // }

    for(auto tmp : placement_candidates_median){
        collectCandidatePositions(tmp);
    }

}//end addMedianCellCandidates

void Cell::addGCellBasedCandidates(){
    if(isFixed) return;
    if(!is_critical_cell) return;

    log() << "addGCellBasedCandidates cell: " << getName() << std::endl;

    auto rsynService = database.getRsynService();
    Rsyn::Cell cellRsyn = rsynInstance.asCell();
    Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cellRsyn);
    Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cellRsyn);
    int cell_width = phLibCell.getWidth();
    int cell_height = phLibCell.getHeight();

    // database.getCellsInBox(new_box,cells,10);
    auto cell_init_position = utils::BoxT<DBU>(rsynInstance.getBounds().getLower().x,
                            rsynInstance.getBounds().getLower().y,
                            rsynInstance.getBounds().getUpper().x,
                            rsynInstance.getBounds().getUpper().y);
    auto cell_init_position_boxOnLayer = db::BoxOnLayer(0, cell_init_position);
    // cellMutex.lock();
    auto grBox = grDatabase.rangeSearchGCell(cell_init_position_boxOnLayer);
 
    std::vector<int> grBox_xs;
    std::vector<int> grBox_ys;
     
    for (int x = grBox[X].low; x <= grBox[X].high; x++)
        for (int y = grBox[Y].low; y <= grBox[Y].high; y++) {
            grBox_xs.push_back(x);
            grBox_ys.push_back(y);
                log() << ", gcell_x: " << x
                      << ", gcell_y: " << y << std::endl;
        }

    int grBox_lx = *std::min_element(grBox_xs.begin(),grBox_xs.end());
    int grBox_ly = *std::min_element(grBox_ys.begin(),grBox_ys.end());
    int grBox_hx = *std::max_element(grBox_xs.begin(),grBox_xs.end());
    int grBox_hy = *std::max_element(grBox_ys.begin(),grBox_ys.end());

    log() << "grBox_lx: "    << grBox_lx
          << ", grBox_ly: "  << grBox_ly
          << ", grBox_hx: "  << grBox_hx
          << ", grBox_hy: "  << grBox_hy << std::endl;

    int ext_grBox_lx = 0;//
    // ((grBox_lx-db::setting.legalizationWindowSize) <= 0) ? 
    //                     0 : (grBox_lx-db::setting.legalizationWindowSize);
    int ext_grBox_ly = 0;
    // ((grBox_ly-db::setting.legalizationWindowSize) <= 0) ? 
    //                     0 : (grBox_ly-db::setting.legalizationWindowSize);
    int ext_grBox_hx = 0;
    // ((grBox_hx+db::setting.legalizationWindowSize) >= grDatabase.getNumGrPoint(X)) ? (grDatabase.getNumGrPoint(X)
    //                     ) : (grBox_hx+db::setting.legalizationWindowSize);
    int ext_grBox_hy = 0;
    // ((grBox_hy+db::setting.legalizationWindowSize) >= grDatabase.getNumGrPoint(Y)) ? (grDatabase.getNumGrPoint(Y)
    //                     ) : (grBox_hy+db::setting.legalizationWindowSize);

    
    utils::BoxT<DBU> legalize_box(grDatabase.getCoor(ext_grBox_lx, X),grDatabase.getCoor(ext_grBox_ly, Y),
                                 grDatabase.getCoor(ext_grBox_hx, X),grDatabase.getCoor(ext_grBox_hy, Y));
    log() << "getCoor(grBox_lx, X): "<< double(grDatabase.getCoor(ext_grBox_lx, X))/database.libDBU << std::endl;
    log() << "getCoor(grBox_ly, Y): "<< double(grDatabase.getCoor(ext_grBox_ly, Y))/database.libDBU << std::endl;
    log() << "getCoor(grBox_hx, X): "<< double(grDatabase.getCoor(ext_grBox_hx, X))/database.libDBU << std::endl;
    log() << "getCoor(grBox_hy, Y): "<< double(grDatabase.getCoor(ext_grBox_hy, Y))/database.libDBU << std::endl;  
    // legalize(cell_init_position,new_box);

    

}//end addGCellBasedCandidates
void Cell::addILPLegalizedCellCandidates(){
    if(isFixed) return;
    if(!is_critical_cell) return;
    bool debug = false;
    // if(getName() == "inst4948"){
    //     log() << "legalizer cell: " << getName() << std::endl;
        Legalizer legalizer(idx);
        if(getName() == "g55131_u0"){
            legalizer.debug_global = true;
            debug = true;
        }
        legalizer.debug_global_all = true;
            
        if(database.policy_set.find("legalizer") == database.policy_set.end())
            legalizer.run();
        
        for(auto sol_tuple : legalizer.legalizer_sols){
            int cell_idx = std::get<0>(sol_tuple);
            int new_row  = std::get<1>(sol_tuple);
            int new_site = std::get<2>(sol_tuple);

            
            auto cell = database.cells[cell_idx];

            auto cell_box = cell.getCellBox();
            auto row_dbu  = database.getDBURow(new_row);
            auto site_dbu = database.getDBUSite(new_site);

            if(debug){
                log() << "legalizer_sols: "
                      << cell.getName()
                      << ", " << site_dbu
                      << ", " << row_dbu << std::endl;
            }
            auto old_site = database.getSiteDBU(cell_box.lx());
            auto old_row  = database.getRowDBU(cell_box.ly());

            if(cell_idx == idx){    
                utils::BoxT<DBU> new_box(site_dbu,row_dbu,
                                        site_dbu+cell_box.width(),row_dbu+cell_box.height());
                if(debug)
                    log() << "cell: " << cell.getName() << ", new_box: " << new_box << std::endl;
                if(!isRedundantCandidate(new_box))
                    placement_candidates.push_back(std::make_pair(new_box,0));
            }else{
                // continue;
                utils::BoxT<DBU> new_box(site_dbu,row_dbu,
                                        site_dbu+cell_box.width(),row_dbu+cell_box.height());
                if(debug)
                    log() << "ovrlp_cell: " << cell.getName() << ", new_box: " << new_box << std::endl;
                ovrlp_cells_.push_back(std::make_pair(cell.idx,new_box));
            }
        }
    // }
    
}

bool Cell::isMultiRowCell(){
    auto box = getCellBox();
    auto height = std::abs(box.hy()-box.ly());
    if(database.rows_step < height) return true;
    return false;
}


void Cell::addMedianCellCandidatesSiteBased(){
    if(isFixed) return;
    if(!is_critical_cell) return;

    log() << "add Median cell: " << getName() << std::endl;

    auto rsynService = database.getRsynService();

    auto median = feature.median_x_y_cell;//getMedianCell();
    auto median_x = median.first;
    auto median_y = median.second;
    
    Rsyn::Cell cellRsyn = rsynInstance.asCell();
    Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cellRsyn);
    Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cellRsyn);
    int cell_width = phLibCell.getWidth();
    int cell_height = phLibCell.getHeight();
    
    auto median_x_site_step = getSiteDBU(median_x);
    auto median_y_site_step = getSiteDBU(median_y);

    median_x = getDBUSite(median_x_site_step);
    median_y = getDBUSite(median_y_site_step);

    // log() << "median_x: " << double(median_x)/database.libDBU << ", median_y: " << double(median_y)/database.libDBU << std::endl;
    auto rows = getRowsNear(median_x,median_y);

    std::vector<utils::BoxT<DBU>> placement_candidates_median;
    // log() << "row and sites for cell: " << getName() << std::endl;
    for (auto row_box : rows){
        placement_candidates_median.emplace_back(median_x, row_box.ly() , median_x+cell_width, row_box.hy());
        placement_candidates_median.emplace_back(
            rsynInstance.getBounds().getLower().x, row_box.ly(),
            rsynInstance.getBounds().getUpper().x, row_box.hy());
        placement_candidates_median.emplace_back(
            median_x, rsynInstance.getBounds().getLower().y,
            median_x+cell_width, rsynInstance.getBounds().getUpper().y);
    }

    log() << "placement_candidates_median..." << std::endl;
    for(auto pl : placement_candidates_median){
        log() << pl << std::endl;
    }


    // add random positions
    // for(int i = 0; i < 100; i++){
    //     auto rand = getRandom();
    //     auto rand_rows = getRowsNear(rand.first,rand.second);
    //     for (auto row_box : rand_rows){
    //         // log() << "ovrlp rows: " << row_box << std::endl;
    //         placement_candidates_median.emplace_back(rand.first,
    //             row_box.ly() , rand.first+cell_width, row_box.hy());
    //     }
    // }

    for(auto tmp : placement_candidates_median){
        collectCandidatePositions(tmp);
    }

}//end addMedianCellCandidates




std::vector<utils::BoxT<DBU>> Cell::getRowsNear(double median_x, double median_y){
    auto rsynService = database.getRsynService();
    Rsyn::Cell cellRsyn = rsynInstance.asCell();
    Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cellRsyn);
    Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cellRsyn);
    int cell_width = phLibCell.getWidth();
    int cell_height = phLibCell.getHeight();

    auto tmpbox = utils::BoxT<DBU>(median_x, median_y,median_x+cell_width, median_y+cell_height/2);
    
    auto grBox = db::BoxOnLayer(0, tmpbox);
    // cellMutex.lock();
    auto gcellBox = grDatabase.rangeSearchGCell(grBox);
    // cellMutex.unlock();

    std::vector<utils::BoxT<DBU>> rows;
    database.getRowsInBox(tmpbox, rows, 1000);
    return rows;
}



void Cell::collectCandidatePositions(utils::BoxT<DBU>& placement_candidate){
    // debug = false;
    std::vector<std::vector<utils::BoxT<DBU>>> solutions_list;
    auto rsynService = database.getRsynService();
    auto dieBound = rsynService.physicalDesign.getPhysicalDie().getBounds();
    // dieRegion = getBoxFromRsynBounds(dieBound);
    std::vector<utils::BoxT<DBU>> sol;

    palceInSite(placement_candidate);

    std::vector<db::Cell> ovrlp_cells;

    // log() << "placement_candidate: " << placement_candidate << std::endl;

    database.getCellsInBox(placement_candidate,ovrlp_cells,10);


    // for(auto tmp : ovrlp_cells){
    //     log() << "ovrlp cells: " << tmp.getName() << std::endl;
    // }

    bool isEmpty = (ovrlp_cells.size() > 0) ? false : true;

    // log() << "isEmpty: " << isEmpty << std::endl;

    if (isEmpty){
        sol.emplace_back(placement_candidate);
        solutions_list.push_back(sol);
        sol.clear();
    }else{
        bool isSame = isSameCellVsOvrlpCells(ovrlp_cells);

        // log() << "isSame: " << isSame << std::endl;

        if (isSame) {
            sol.emplace_back(placement_candidate);
            solutions_list.push_back(sol);
            sol.clear();
        }else{
            utils::BoxT<DBU> new_box;
            // eps 10 here has the risk of error please review it later.
            // std::vector<utils::BoxT<DBU>> sol;
            ovrlp_cells.clear();
            sol.clear();
            if(getCellsInBoxNeigh(placement_candidate,ovrlp_cells,new_box,10,"M"))
                sol.emplace_back(new_box);
            ovrlp_cells.clear();
            // if(getCellsInBoxNeigh(placement_candidate,ovrlp_cells,new_box,10,"N"))
            //     sol.emplace_back(new_box);
            // ovrlp_cells.clear();
            // if(getCellsInBoxNeigh(placement_candidate,ovrlp_cells,new_box,10,"S"))
            //     sol.emplace_back(new_box);
            // ovrlp_cells.clear();
            // if(getCellsInBoxNeigh(placement_candidate,ovrlp_cells,new_box,10,"W"))
            //     sol.emplace_back(new_box);
            // ovrlp_cells.clear();
            // if(getCellsInBoxNeigh(placement_candidate,ovrlp_cells,new_box,10,"E"))
            //     sol.emplace_back(new_box);
            // ovrlp_cells.clear();
            // if(getCellsInBoxNeigh(placement_candidate,ovrlp_cells,new_box,10,"SE"))
            //     sol.emplace_back(new_box);
            // ovrlp_cells.clear();
            // if(getCellsInBoxNeigh(placement_candidate,ovrlp_cells,new_box,10,"SW"))
            //     sol.emplace_back(new_box);
            // ovrlp_cells.clear();
            // if(getCellsInBoxNeigh(placement_candidate,ovrlp_cells,new_box,10,"NW"))
            //     sol.emplace_back(new_box);
            // ovrlp_cells.clear();
            // if(getCellsInBoxNeigh(placement_candidate,ovrlp_cells,new_box,10,"NE"))
            //     sol.emplace_back(new_box);


            solutions_list.push_back(sol);
            sol.clear();
        }
    }
    for(auto cells_row_cand: solutions_list){
        for(auto cell_cand : cells_row_cand){
            if(!isRedundantCandidate(cell_cand))
                placement_candidates.emplace_back(cell_cand,0);            
        }
    }

}//end collectCandidatePositions

bool Cell::getCellsInBoxNeigh(utils::BoxT<DBU>& cur_box,std::vector<db::Cell>& cells,
                            utils::BoxT<DBU>& new_box,DBU eps,std::string mode){

    if(mode == "M") {
        new_box = utils::BoxT<DBU>(cur_box.lx(),cur_box.ly(),
                                   cur_box.hx(),cur_box.hy());                                
    }else if(mode == "N") {
        new_box = utils::BoxT<DBU>(cur_box.lx(),cur_box.ly()+ cur_box.height(),
                                       cur_box.hx(),cur_box.hy()+ cur_box.height());
    }else if(mode == "S") {
        new_box = utils::BoxT<DBU>(cur_box.lx(),cur_box.ly()- cur_box.height(),
                                       cur_box.hx(),cur_box.hy()- cur_box.height());
    
    }else if(mode == "E") {
        new_box = utils::BoxT<DBU>(cur_box.lx()+cur_box.width(),cur_box.ly(),
                                       cur_box.hx()+cur_box.width(),cur_box.hy());
    }else if(mode == "W") {
        new_box = utils::BoxT<DBU>(cur_box.lx()-cur_box.width(),cur_box.ly(),
                                       cur_box.hx()-cur_box.width(),cur_box.hy());
    }else if(mode == "NE") {
        new_box = utils::BoxT<DBU>(cur_box.lx()+cur_box.width(),cur_box.ly()+ cur_box.height(),
                                   cur_box.hx()+cur_box.width(),cur_box.hy()+ cur_box.height());
    }else if(mode == "NW") {
        new_box = utils::BoxT<DBU>(cur_box.lx()-cur_box.width(),cur_box.ly()+ cur_box.height(),
                                       cur_box.hx()-cur_box.width(),cur_box.hy()+ cur_box.height());

    }else if(mode == "SW") {
        new_box = utils::BoxT<DBU>(cur_box.lx()-cur_box.width(),cur_box.ly()-cur_box.height(),
                                       cur_box.hx()-cur_box.width(),cur_box.hy()-cur_box.height());

    }else if(mode == "SE") {
        new_box = utils::BoxT<DBU>(cur_box.lx()+cur_box.width(),cur_box.ly()-cur_box.height(),
                                    cur_box.hx()+cur_box.width(),cur_box.hy()-cur_box.height());  
    }
            
    // log() << "new box " << new_box << "mode: " << mode << std::endl;
    bool isValid = isValidPostion(new_box);
    // log() << "isValid " << isValid << std::endl;
    if (!isValid){
        return false;
    }
    // eps 10 here has the risk of error please review it later.
    database.getCellsInBox(new_box,cells,10);
     auto box_init_position = 
            utils::BoxT<DBU>(rsynInstance.getBounds().getLower().x,
                                rsynInstance.getBounds().getLower().y,
                                rsynInstance.getBounds().getUpper().x,
                                rsynInstance.getBounds().getUpper().y);

    // if(cells.size() > 1) return false;       

    // for(auto cell : cells){
    //     log() << "move_cell: " << getName() 
    //           << ", overlap_cell: " << cell.getName()
    //           << ", isFixed: " << cell.isFixed << std::endl;
    // }
    
    for(auto cell : cells){
        // log() << "move_cell: " << getName() 
        //       << ", cell ovrlp: " << cell.getName() << ", isFixed: "
        //       << cell.isFixed << std::endl;
        if ( (connectd_cells.find(cell.idx) == connectd_cells.end()) && 
            (cell.rsynInstance.getWidth() == rsynInstance.getWidth()) && 
            (!cell.isFixed)){
            ovrlp_cells_.push_back(std::make_pair(cell.idx,box_init_position));
        }else{
            return false;
        }
        
    }

    // log() << "good to add solution list" << std::endl;
    
    return true;

}

int Cell::getCellSiteWidth(){
    return std::abs(database.getSiteDBU(rsynInstance.getBounds().getUpper().x)-
                    database.getSiteDBU(rsynInstance.getBounds().getLower().x));

}

bool Cell::isSameCellVsOvrlpCells(std::vector<db::Cell>& cells){
    for(auto& cell_tmp : cells){
        if(cell_tmp.getName() != rsynInstance.getName()){
            return false;
        }
    }
    return true;
}//end isSameCellVsOvrlpCells

bool Cell::isValidPostion(utils::BoxT<DBU>& box){
    auto rsynService = database.getRsynService();
    auto dieBound = rsynService.physicalDesign.getPhysicalDie().getBounds();
    bool isValid = false;

    isValid = ((dieBound.getLower().x <= box.lx()) &&
               (dieBound.getLower().y <= box.ly()) &&
               (dieBound.getUpper().x >= box.hx()) &&
               (dieBound.getUpper().y >= box.hy())) ? true : false;
    

    return isValid;
}


void Cell::getDistOptPoint(){
    auto rsynService = database.getRsynService();
    std::vector<db::Net> nets;
    std::vector<int> nets_idx;
  
    getCellNetsIdx(nets_idx);

    for(auto idx : nets_idx){
        // log() << "nets_idx: " << idx << std::endl;
        nets.push_back(database.nets[idx]);
    }
        
    std::vector<int> cells_idx;
    std::vector<db::Cell> cells;
    
    // log() << "connected nets: " << std::endl;
    for(auto net : nets){
        getNetCells(net.idx,cells_idx);
        // log() << "net.idx: " << net.idx<< std::endl;
    }
    // log() << "cells_idx..." << std::endl;
    
    for(int cell_idx : cells_idx){
        // log() << "cell_idx: " << cell_idx<< std::endl;
        cells.push_back(database.cells[cell_idx]);
    }
    
    std::vector<double> x_positions;
    std::vector<double> y_positions;
    
    for (auto cell_tmp : cells){
        x_positions.push_back(cell_tmp.rsynInstance.getBounds().getLower().x);
        y_positions.push_back(cell_tmp.rsynInstance.getBounds().getLower().y);
    }
 
    
    if(x_positions.empty() || y_positions.empty()){
        // log() << "no median " << std::endl;
        return;
    }

    std::nth_element(x_positions.begin(), x_positions.begin() + x_positions.size()/2, x_positions.end());
    auto median_x = x_positions[x_positions.size()/2];
    std::nth_element(y_positions.begin(), y_positions.begin() + y_positions.size()/2, y_positions.end());
    auto median_y = y_positions[y_positions.size()/2];
    // log() << "median_x: " << median_x/database.libDBU 
    //       << ", median_y: " << median_y/database.libDBU  << std::endl;

    cost =  std::abs(median_x - rsynInstance.getBounds().getLower().x)+
            std::abs(median_y - rsynInstance.getBounds().getLower().y);

    // log() << "inner cell: " << getName() << "cost: " << cost << std::endl;
}//end getDistOptPoint

void Cell::getNetsBox(utils::BoxT<DBU>& box){
    auto rsynService = database.getRsynService();
    std::vector<db::Net> nets;
    std::vector<int> nets_idx;
  
    getCellNetsIdx(nets_idx);

    for(auto idx : nets_idx){
        // log() << "nets_idx: " << idx << std::endl;
        nets.push_back(database.nets[idx]);
    }
        
    std::vector<int> cells_idx;
    std::vector<db::Cell> cells;
    
    // log() << "connected nets: " << std::endl;
    for(auto net : nets){
        getNetCells(net.idx,cells_idx);
        // log() << "net.idx: " << net.idx<< std::endl;
    }
    // log() << "cells_idx..." << std::endl;
    
    for(int cell_idx : cells_idx){
        // log() << "cell_idx: " << cell_idx<< std::endl;
        cells.push_back(database.cells[cell_idx]);
    }
    
    std::vector<DBU> xs;
    std::vector<DBU> ys;
    
    for (auto cell_tmp : cells){
        xs.push_back(cell_tmp.rsynInstance.getBounds().getLower().x);
        ys.push_back(cell_tmp.rsynInstance.getBounds().getLower().y);
    }
    DBU max_x = *std::max_element(xs.begin(), xs.end());
    DBU min_x = *std::min_element(xs.begin(), xs.end());
    DBU max_y = *std::max_element(ys.begin(), ys.end());
    DBU min_y = *std::min_element(ys.begin(), ys.end());
    box.Set(min_x,min_y,max_x,max_y);
}//end getNetsBox


void Cell::getDistOptPointV2(){
    auto rsynService = database.getRsynService();
    std::vector<db::Net> nets;
    std::vector<int> nets_idx;

    // log() << "getDistOptPointV2 cell: " << rsynInstance.getName() << std::endl;
    
    getCellNetsIdx(nets_idx);

    for(auto idx : nets_idx){
        // log() << "nets_idx: " << idx << std::endl;
        nets.push_back(database.nets[idx]);
    }
        
    
    
    // log() << "connected nets: " << std::endl;
    std::vector<double> med_xs;
    std::vector<double> med_ys;
    for(auto net : nets){
        std::vector<double> x_positions;
        std::vector<double> y_positions;
        std::vector<int> cells_idx;
        std::vector<db::Cell> cells;
        getNetCells(net.idx,cells_idx);
        for(int cell_idx : cells_idx){
            if(cell_idx == idx) continue;
            cells.push_back(database.cells[cell_idx]);
        }
        for (auto cell_tmp : cells){
            x_positions.push_back(cell_tmp.rsynInstance.getBounds().getLower().x);
            y_positions.push_back(cell_tmp.rsynInstance.getBounds().getLower().y);
        }
        if(x_positions.empty() || y_positions.empty()){
            // log() << "no median " << std::endl;
            return;
        }

        // log() << "x_positions ..." << std::endl;
        // for(auto x : x_positions)
        //     log() << "x: " << x << std::endl;

        // log() << "y_positions ..." << std::endl;
        // for(auto y : y_positions)
        //     log() << "y: " << y << std::endl;

        std::nth_element(x_positions.begin(), x_positions.begin() + x_positions.size()/2, x_positions.end());
        auto median_x = x_positions[x_positions.size()/2];
        std::nth_element(y_positions.begin(), y_positions.begin() + y_positions.size()/2, y_positions.end());
        auto median_y = y_positions[y_positions.size()/2];
        // log() << "median_x: " << median_x/database.libDBU 
        //       << ", median_y: " << median_y/database.libDBU  << std::endl;

        med_xs.push_back(median_x);
        med_ys.push_back(median_y);
    }

    double len_x = med_xs.size();
    double len_y = med_ys.size();
    if(len_x == 0 || len_y == 0){
        cost = 0;
        return;
    }

    double avg_x = accumulate( med_xs.begin(), med_xs.end(), 0.0) / len_x; 
    double avg_y = accumulate( med_ys.begin(), med_ys.end(), 0.0) / len_y; 


    cost =  std::abs(avg_x - rsynInstance.getBounds().getLower().x)+
            std::abs(avg_y - rsynInstance.getBounds().getLower().y);

}//end getDistOptPoint

void Cell::updateMedianCellV2(){

}//end updateMedianCellV2

void Cell::updateMedianCell(){
    auto rsynService = database.getRsynService();
    std::vector<db::Net> nets;
    std::vector<int> nets_idx;

    // log() << "getDistOptPointV2 cell: " << rsynInstance.getName() << std::endl;
    
    getCellNetsIdx(nets_idx);

    for(auto idx : nets_idx){
        // log() << "nets_idx: " << idx << std::endl;
        nets.push_back(database.nets[idx]);
    }
        
    // log() << "updateMedianCell: " << getName() << std::endl;
    
    // log() << "connected nets: " << std::endl;
    std::vector<double> med_xs;
    std::vector<double> med_ys;
    for(auto net : nets){
        // log() << "net: " << net.getName() << std::endl;
        std::vector<double> x_positions;
        std::vector<double> y_positions;
        std::vector<int> cells_idx;
        std::vector<db::Cell> cells;
        getNetCells(net.idx,cells_idx);
        for(int cell_idx : cells_idx){
            if(cell_idx == idx) continue;
            cells.push_back(database.cells[cell_idx]);
        }
        for (auto cell_tmp : cells){
            // log() << "cell_tmp: " << cell_tmp.getName() << std::endl;
            x_positions.push_back(cell_tmp.rsynInstance.getBounds().getLower().x);
            y_positions.push_back(cell_tmp.rsynInstance.getBounds().getLower().y);
        }
        if(x_positions.empty() || y_positions.empty()){
            // log() << "no median " << std::endl;
            median_x_y_cell = std::make_pair(0,0);
            return;
        }

        // log() << "x_positions ..." << std::endl;
        // for(auto x : x_positions)
        //     log() << "x: " << x << std::endl;

        // log() << "y_positions ..." << std::endl;
        // for(auto y : y_positions)
        //     log() << "y: " << y << std::endl;

        // for(int i = 0; i < x_positions.size() ; i++){
        //     log() << "x: " << double(x_positions[i])/database.libDBU
        //           << ", y: " << double(y_positions[i])/database.libDBU << std::endl;
        // }


        std::nth_element(x_positions.begin(), x_positions.begin() + x_positions.size()/2, x_positions.end());
        auto median_x = x_positions[x_positions.size()/2];
        std::nth_element(y_positions.begin(), y_positions.begin() + y_positions.size()/2, y_positions.end());
        auto median_y = y_positions[y_positions.size()/2];
        // log() << "median_x: " << double(median_x)/database.libDBU
        //       << ", median_y: " << double(median_y)/database.libDBU  << std::endl;

        med_xs.push_back(median_x);
        med_ys.push_back(median_y);
    }

    double len_x = med_xs.size();
    double len_y = med_ys.size();
    if(len_x == 0 || len_y == 0){
        median_x_y_cell = std::make_pair(0,0);
        return;
    }

    double avg_x = accumulate( med_xs.begin(), med_xs.end(), 0.0) / len_x; 
    double avg_y = accumulate( med_ys.begin(), med_ys.end(), 0.0) / len_y; 

    

    median_x_y_cell = std::make_pair(avg_x,avg_y);
    
}//end updateMedianCell


std::pair<double,double> Cell::getMedianBoxs(std::vector<utils::BoxT<DBU>>& boxs){
    std::vector<double> med_xs;
    std::vector<double> med_ys;
    std::vector<double> x_positions;
    std::vector<double> y_positions;

    for(auto box : boxs){
        x_positions.push_back(box.lx());
        y_positions.push_back(box.ly());
    }

    if(x_positions.empty() || y_positions.empty()){
        // log() << "no median " << std::endl;
        return std::make_pair(0,0);
        
    }

    std::nth_element(x_positions.begin(), x_positions.begin() + x_positions.size()/2, x_positions.end());
    auto median_x = x_positions[x_positions.size()/2];
    std::nth_element(y_positions.begin(), y_positions.begin() + y_positions.size()/2, y_positions.end());
    auto median_y = y_positions[y_positions.size()/2];
    
    return std::make_pair(median_x,median_y);

}//end getMedianBoxs

void Cell::updateMedianPin(){
    auto rsynService = database.getRsynService();        
    const Rsyn::Session session;
    const Rsyn::PhysicalDesign& physicalDesign =
        static_cast<Rsyn::PhysicalService*>(session.getService("rsyn.physical"))->getPhysicalDesign();
    const DBU libDBU = physicalDesign.getDatabaseUnits(Rsyn::LIBRARY_DBU);
    std::vector<db::Net> nets;
    std::vector<int> nets_idx;

    // log() << "getDistOptPointV2 cell: " << rsynInstance.getName() << std::endl;
    
    getCellNetsIdx(nets_idx);

    for(auto idx : nets_idx){
        // log() << "nets_idx: " << idx << std::endl;
        nets.push_back(database.nets[idx]);
    }
        
    
    // log() << "updateMedianPin cell: " << getName() << std::endl;
    // log() << "connected nets: " << std::endl;
    std::vector<double> med_xs;
    std::vector<double> med_ys;      
    for(auto net : nets){
        // log() << "net: " << net.getName() << std::endl;
        std::vector<double> x_positions;
        std::vector<double> y_positions;
        std::vector<int> cells_idx;
        std::vector<db::Cell> cells;
        
        auto box = getCellBox();

        for (auto rsynPin : net.rsynNet.allPins()) {
            vector<BoxOnLayer> pinAccessBoxes; 
            std::vector<double> x_acs;
            std::vector<double> y_acs;
            if(rsynPin.getInstance().getName() == getName()) continue;
            //rsynPin, rsynService
            //  rsynPins.push_back(RsynPin);
            // pinAccessBoxes.emplace_back();
            // log() << "instance: " << rsynPin.getInstance().getName() 
            //       << ", pin: " << rsynPin.getName() << std::endl;
            // if(RsynPin.getInstance().getName() != getName())
            // log() << "initPinAccessBoxes..." << std::endl;
            initPinAccessBoxes(rsynPin, rsynService, pinAccessBoxes, box,libDBU);
            // log() << "pinAccessboxes: "  << std::endl;
            // for(auto tmp : pinAccessBoxes){
            //     log() << "pinAccess: " << tmp << std::endl;
            // }
            std::vector<utils::BoxT<DBU>> acBoxs;
            for(auto tmpBoxOnLayer : pinAccessBoxes){
                x_acs.push_back(tmpBoxOnLayer.lx());
                x_acs.push_back(tmpBoxOnLayer.hx());
                y_acs.push_back(tmpBoxOnLayer.ly());
                y_acs.push_back(tmpBoxOnLayer.hy());
            }

            auto min_x_acs = *std::min_element(x_acs.begin(),x_acs.end());
            auto min_y_acs = *std::min_element(y_acs.begin(),y_acs.end());

            // log() << "median_acBoxs: [x: " << double(median_acBoxs.first)/database.libDBU
            //       << ", y: " << double(median_acBoxs.second)/database.libDBU << " ]" << std::endl;

            // log() << "min_x_acs: " << double(min_x_acs)/database.libDBU
            //       << ", min_y_acs: " << double(min_y_acs)/database.libDBU << std::endl;

            // for(auto acBox : pinAccessBoxes){
            x_positions.push_back(min_x_acs);
            y_positions.push_back(min_y_acs);
            // }

        }

        if(x_positions.empty() || y_positions.empty()){
            // log() << "no median " << std::endl;
            median_x_y_pin = std::make_pair(0,0);
            return;
        }

        // for(int i = 0; i < x_positions.size() ; i++){
        //     log() << "x: " << double(x_positions[i])/database.libDBU
        //           << ", y: " << double(y_positions[i])/database.libDBU << std::endl;
        // }



        std::nth_element(x_positions.begin(), x_positions.begin() + x_positions.size()/2, x_positions.end());
        auto median_x = x_positions[x_positions.size()/2];
        std::nth_element(y_positions.begin(), y_positions.begin() + y_positions.size()/2, y_positions.end());
        auto median_y = y_positions[y_positions.size()/2];
        // log() << "median_x: " << double(median_x)/database.libDBU
        //       << ", median_y: " << double(median_y)/database.libDBU  << std::endl;

        med_xs.push_back(median_x);
        med_ys.push_back(median_y);
    }

    double len_x = med_xs.size();
    double len_y = med_ys.size();
    if(len_x == 0 || len_y == 0){
        median_x_y_pin = std::make_pair(0,0);
        return;
    }

    

    double avg_x = accumulate( med_xs.begin(), med_xs.end(), 0.0) / len_x; 
    double avg_y = accumulate( med_ys.begin(), med_ys.end(), 0.0) / len_y; 
    median_x_y_pin = std::make_pair(avg_x,avg_y);

    // log() << "median_x_y_pin: [x: " << double(median_x_y_pin.first)/database.libDBU
    //     << ", y: " << double(median_x_y_pin.second)/database.libDBU << " ]" << std::endl;

    
}//end updateMedianCell


std::pair<double,double> Cell::getRandom(){
    auto rsynService = database.getRsynService();
    auto dieBound = rsynService.physicalDesign.getPhysicalDie().getBounds();

    auto die_lx = dieBound.getLower().x;
    auto die_ly = dieBound.getLower().y;
    auto die_hx = dieBound.getUpper().x;
    auto die_hy = dieBound.getUpper().y;


    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist_x(die_lx,die_hx); // distribution in range [1, 6]
    std::uniform_int_distribution<std::mt19937::result_type> dist_y(die_ly,die_hy); // distribution in range [1, 6]

    auto x = dist_x(rng);
    auto y = dist_y(rng);

    // log() << "random: [x: " << x << ", y: " << y << "]"<<std::endl;
               
    
    return std::make_pair(x,y);

}//end getRandom

std::pair<double,double> Cell::getRandom(utils::BoxT<DBU>& box){
    auto rsynService = database.getRsynService();
    auto dieBound = rsynService.physicalDesign.getPhysicalDie().getBounds();

    auto die_lx = dieBound.getLower().x;
    auto die_ly = dieBound.getLower().y;
    auto die_hx = dieBound.getUpper().x;
    auto die_hy = dieBound.getUpper().y;

    auto box_lx = box.lx();
    auto box_ly = box.ly();
    auto box_hx = box.hx();
    auto box_hy = box.hy();

    if(box_ly < die_ly) box_ly = die_ly;
    if(box_lx < die_lx) box_lx = die_lx;
    if(box_hy > die_hy) box_hy = die_hy;
    if(box_hx > die_hx) box_hx = die_hx;




    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist_x(box_lx,box_hx); // distribution in range [1, 6]
    std::uniform_int_distribution<std::mt19937::result_type> dist_y(box_ly,box_hy); // distribution in range [1, 6]

    auto x = dist_x(rng);
    auto y = dist_y(rng);

    // log() << "random: [x: " << x << ", y: " << y << "]"<<std::endl;
               
    
    return std::make_pair(x,y);

}//end getRandom

Rsyn::PhysicalOrientation Cell::getMoveOrientation(utils::BoxT<DBU>& box){
     int eps = 10;
    // get row 
    boostBox rtreeQueryBox(boostPoint(box.lx()+eps,
                                        box.ly()+eps),
                            boostPoint(box.hx()-eps,
                                        box.hy()-eps));
    vector<std::pair<boostBox, int>> queryResults;
    
    database.row_rtree.query(bgi::intersects(rtreeQueryBox), std::back_inserter(queryResults));
    for (const auto& queryResult : queryResults) {
        auto row_idx = queryResult.second;        
        auto row = database.rowsRsyn[row_idx];
        return row.getSiteOrientation();
        // break; 
    }//end for
    return Rsyn::ORIENTATION_INVALID;
}//end getMoveOrientation

DBU Cell::getWirelength(){
    std::vector<gr::GrNet> nets;
    std::vector<int> nets_idx;
    
    getCellNetsIdx(nets_idx);

    for(auto idx : nets_idx){
        nets.push_back(grDatabase.nets[idx]);
    }
    DBU wl = 0;
    for(auto net : nets){
        wl += net.getWirelength();
    }
    return wl;

}//end getWirelength

void Cell::updateCellFeatures(){
    if(database.policy_set.find("updateMedianCell") == database.policy_set.end())
        updateMedianCell();
    if(database.policy_set.find("updateMedianPin") == database.policy_set.end())
        updateMedianPin();
    int wl = getWirelengthLast();    
    std::vector<int> nets_idx;
    getCellNetsIdx(nets_idx);
    connectd_cells.size();


    feature.median_x_y_cell = median_x_y_cell;
    feature.median_x_y_pin = median_x_y_pin;
    feature.wl = wl;
    feature.degree_connected_nets = nets_idx.size();
    feature.degree_connected_cells = connectd_cells.size();
    feature.congestion = 0;

    if(database.suspected_cells_dict.find(idx) != database.suspected_cells_dict.end()){
        feature.is_suspected_cell = true;
    }else{
        feature.is_suspected_cell = false;
    }

    if(database.policy_set.find("getMaxEdgeCost") == database.policy_set.end())
        for(auto net_idx : nets_idx){
            feature.max_edge_cost += grDatabase.nets[net_idx].getMaxEdgeCost();
        }

    if(database.policy_set.find("connected_nets_degrees") == database.policy_set.end())
        for(auto idx : nets_idx){
            // log() << "net: " << database.nets[idx].getName() 
            //       << "num pins: " << database.nets[idx].rsynPins.size() << std::endl;
            feature.connected_nets_degrees.push_back(std::make_pair(
                idx,database.nets[idx].rsynPins.size()
            ));
            
        }

    
    
}

DBU Cell::getWirelengthLast(){
    std::vector<gr::GrNet> nets;
    std::vector<int> nets_idx;
    
    getCellNetsIdx(nets_idx);

    for(auto idx : nets_idx){
        nets.push_back(grDatabase.nets[idx]);
    }
    DBU wl = 0;
    for(auto net : nets){
        wl += net.getWirelengthLast();
    }
    return wl;

}//end getWirelength


DBU Cell::getPinIoU(vector<vector<BoxOnLayer>>& pinAccessBoxes,int tgPin_idx){
    DBU area = 0;
    // box of unions
    vector<utils::BoxT<DBU>> box_us;
    for (auto boxs : pinAccessBoxes) {
        utils::BoxT<DBU> box_u;
        for(auto box : boxs){
            auto tmpbox = utils::BoxT<DBU>(box.lx(), box.ly(),box.hx(), box.hy());            
            box_u = box_u.UnionWith(tmpbox);
        }
        box_us.push_back(box_u);
    }
    auto box_u_tg = box_us[tgPin_idx];

    Bounds box_u_tg_x = Bounds(box_u_tg.lx(),0,box_u_tg.hx(),1);
    Bounds box_u_tg_y = Bounds(0,box_u_tg.ly(),1,box_u_tg.hy());

    
    // vertical alignment 
    for(int i = 0; i < box_us.size(); i++){
        if(i == tgPin_idx) continue;
        // utils::BoxT<DBU> box_u_x = utils::BoxT<DBU>(box_us[i].lx(),0,box_us[i].hx(),1);
        // utils::BoxT<DBU> box_u_y = utils::BoxT<DBU>(0,box_us[i].ly(),1,box_us[i].hy());
        auto box_u = box_us[i];
        Bounds box_u_bds = Bounds(box_u.lx(),box_u.ly(),box_u.hx(),box_u.hy());
        Bounds box_u_x = Bounds(box_u_bds.getLower().x,0,box_u_bds.getUpper().x,1);
        Bounds box_u_y = Bounds(0,box_u_bds.getLower().y,1,box_u_bds.getUpper().y);
        // auto box_x = box_u_tg_x.IntersectWith(box_u_x);
        // auto box_y = box_u_tg_y.IntersectWith(box_u_y);
        auto box_x_bds = box_u_tg_x.overlapRectangle(box_u_x);
        auto box_y_bds = box_u_tg_y.overlapRectangle(box_u_y);

        utils::BoxT<DBU> box_x = utils::BoxT<DBU>(box_x_bds.getLower().x,box_x_bds.getLower().y
                                                 ,box_x_bds.getUpper().x,box_x_bds.getUpper().y);
        utils::BoxT<DBU> box_y = utils::BoxT<DBU>(box_y_bds.getLower().x,box_y_bds.getLower().y
                                                 ,box_y_bds.getUpper().x,box_y_bds.getUpper().y);

        area = area + box_x.area() + box_y.area();
    }
    return area;
}//end getIoUX 

DBU Cell::getPinTotalArea(vector<vector<BoxOnLayer>>& pinAccessBoxes,int tgPin_idx){
    DBU area = 0;
    // box of unions
    vector<utils::BoxT<DBU>> box_us;
    for (auto boxs : pinAccessBoxes) {
        utils::BoxT<DBU> box_u;
        for(auto box : boxs){
            auto tmpbox = utils::BoxT<DBU>(box.lx(), box.ly(),box.hx(), box.hy());
            box_u = box_u.UnionWith(tmpbox);
        }
        box_us.push_back(box_u);
    }

    // vertical alignment 
    for(int i = 0; i < box_us.size(); i++){
        if(i == tgPin_idx) continue;

        auto box_u = box_us[i];
        area = area + box_u.area();
    }
    return area;
}//end getPinTotalArea

void Cell::palceInSite(utils::BoxT<DBU>& box){
    // eq1: step = (cur_loc*sites_origin)/sites_step;
    int step = (box.lx()-database.sites_origin)/database.sites_step;
    // eq2: lx = step*sites_origin;
    DBU new_lx = database.sites_origin + step*database.sites_step;
    DBU box_width = box.width();
    box.Set(new_lx,box.ly(),new_lx + box_width, box.hy());
}//end palceInSite

int Cell::getSiteDBU(DBU x){
     int step = (x-database.sites_origin)/database.sites_step;
    // eq2: lx = step*sites_origin;
    return step;
}

DBU Cell::getDBUSite(int step){
    return database.sites_origin + step*database.sites_step;
}

void Cell::getIntesectedSites(utils::BoxT<DBU>& box,std::vector<int>& sites){
    int site_lx = (box.lx()-database.sites_origin)/database.sites_step;
    int site_hx = (box.hx()-database.sites_origin)/database.sites_step;
    for(auto i = site_lx; i < site_hx; i++)
        sites.push_back(i);
}//end getIntersectedSites

int Cell::getRowDBU(DBU y){
     int step = (y-database.rows_origin)/database.rows_step;
    // eq2: lx = step*sites_origin;
    return step;
}

DBU Cell::getDBURow(int step){
    return database.rows_origin + step*database.rows_step;
}

void Cell::getIntesectedRows(utils::BoxT<DBU>& box,std::vector<int>& rows){
    int row_lx = (box.ly()-database.rows_origin)/database.rows_step;
    int row_hx = (box.hy()-database.rows_origin)/database.rows_step;
    for(auto i = row_lx; i < row_hx; i++)
        rows.push_back(i);
}//end getIntersectedSites


void Cell::printPlacementCandidates(){
    log() << "placement_candidates cell: " << getName() << std::endl;
    for(auto tmp : placement_candidates){
        log() << "pl_cand: " << tmp.first << ", cost: " << tmp.second << std::endl;
    }
}//end printPlacementCandidates

utils::BoxT<DBU> Cell::getCellBox(){
    auto cell_box = utils::BoxT<DBU>( rsynInstance.getBounds().getLower().x,
                                rsynInstance.getBounds().getLower().y,
                                rsynInstance.getBounds().getUpper().x,
                                rsynInstance.getBounds().getUpper().y);
    return cell_box;

}//end getCellBox

int Cell::getCellRow(){
    int step = 0;
    auto cell_box = getCellBox();
    step = getRowDBU(cell_box.ly());

    return step;

}//end getCellRow
int Cell::getCellSite(){
    int step = 0;
    auto cell_box = getCellBox();
    step = getSiteDBU(cell_box.lx());

    return step;

}//end getCellSite


void CellList::init(RsynService& rsynService) {
    cells.clear();
    cells.reserve(rsynService.design.getNumInstances());

    for (Rsyn::Instance instance : rsynService.module.allInstances()) {
        if (instance.getType() != Rsyn::CELL) continue;
        cells.emplace_back(cells.size(), instance, rsynService);
    }//end movement loop
}



}  // namespace db
