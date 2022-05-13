#include "Placer.h"
#include <random>
#include <iostream>
#include <execution>
#include "solver.h"

void Placer::runMTISPD(std::vector<int>& netsToRoute,int cellWidth, int cellHeight){
    congMap.init(cellWidth, cellHeight);

    log() << "init cells and update median of each..." << std::endl;
    initCellsUpdateFeatures();
    log() << "end init cells and update median of each!" << std::endl;

    log() << "selectCriticalCells..." << std::endl;
    selectCriticalCells();
    log() << "end selectCriticalCells..." << std::endl;

   

}//end runMTISPD

void Placer::runMT(std::vector<int>& netsToRoute,int cellWidth, int cellHeight){
    bool log_debug = false;
    // return;
    congMap.init(cellWidth, cellHeight);
    // profile_time_str << "congMap_" << std::to_string(iter) << "," << std::to_string(profile_time.getTimer()) << std::endl;
    profile_time_str << "congMap" << "," << std::to_string(profile_time.getTimer()) << std::endl;
    // log() << "setting.name: " << db::setting.name << std::endl;

    // log() << "after routing path cost " << std::endl;
    // for(auto grNet : grDatabase.nets){
    //     log() << "grNet: " << grNet.getName() << ", path_cost: " << grNet.getPathCost() << std::endl;
    // }
    if(db::setting.debug){
        log() << "init cells and update median of each..." << std::endl;
    }
    initCellsUpdateFeatures();
    if(db::setting.debug){
        log() << "end init cells and update median of each!" << std::endl;
    }
    if(log_debug)
        logCellsFeature();


    // profile_time_str << "initCellsUpdateFeatures_" << std::to_string(iter) << "," << std::to_string(profile_time.getTimer()) << std::endl;
    profile_time_str << "initCellsUpdateFeatures" << "," << std::to_string(profile_time.getTimer()) << std::endl;

    if(db::setting.debug){
        log() << "selectCriticalCells..." << std::endl;
    }        
    selectCriticalCells();
    if(db::setting.debug){
        log() << "end selectCriticalCells..." << std::endl;
    }
    if(log_debug)
        logCellsCritical();
    // profile_time_str << "selectCriticalCells_" << std::to_string(iter) << "," << std::to_string(profile_time.getTimer()) << std::endl;
    profile_time_str << "selectCriticalCells"<< "," << std::to_string(profile_time.getTimer()) << std::endl;
    
    // return;
    if(db::setting.debug){
        log() << "initCellsCandidate..." << std::endl;
    }
    initCellsCandidate();
    if(db::setting.debug){
        log() << "end initCellsCandidate..." << std::endl;
    }
    // profile_time_str << "initCellsCandidate_" << std::to_string(iter) << "," << std::to_string(profile_time.getTimer()) << std::endl;
    profile_time_str << "initCellsCandidate"  << "," << std::to_string(profile_time.getTimer()) << std::endl;
    
    if(database.policy_set.find("calcCellsCostMT") == database.policy_set.end()){
        if(db::setting.debug){
            log() << "calcCellsCostMT..." << std::endl;
        }
        calcCellsCostMT();
        if(db::setting.debug){
            log() << "end calcCellsCostMT..." << std::endl;
        }
    }

    // if(log_debug)
    logCellsCandidates();
    // profile_time_str << "calcCellsCostMT_" << std::to_string(iter) << "," << std::to_string(profile_time.getTimer()) << std::endl;
    profile_time_str << "calcCellsCostMT"<< "," << std::to_string(profile_time.getTimer()) << std::endl;
    if(db::setting.debug){
        log() << "findBestEntryMT..." << std::endl;
    }
    findBestEntryMTV2();
    if(db::setting.debug){
        log() << "end findBestEntryMT..." << std::endl;
    }
    // profile_time_str << "findBestEntryMT_" << std::to_string(iter) << "," << std::to_string(profile_time.getTimer()) << std::endl;
    profile_time_str << "findBestEntryMT"<< "," << std::to_string(profile_time.getTimer()) << std::endl;

    // log() << "logRefinePlacement..." << std::endl;
    // logRefinePlacement();
    // log() << "end logRefinePlacement..." << std::endl;
    if(log_debug)
        logCellsMoved();
    if(db::setting.debug){
        log() << "refinePlacement..." << std::endl;
    }
    refinePlacement(netsToRoute);   
    if(db::setting.debug){
        log() << "end refinePlacement..." << std::endl;
    }
    // profile_time_str << "refinePlacement_" << std::to_string(iter) << "," << std::to_string(profile_time.getTimer()) << std::endl;
    profile_time_str << "refinePlacement"<< "," << std::to_string(profile_time.getTimer()) << std::endl;

    if(db::setting.debug){
        log() << "updateRoutingDataBase..." << std::endl;
    }
    updateRoutingDataBase(netsToRoute);
    if(db::setting.debug){
        log() << "end updateRoutingDataBase..." << std::endl;
    }
    // profile_time_str << "updateRoutingDataBase_" << std::to_string(iter) << "," << std::to_string(profile_time.getTimer()) << std::endl;
    profile_time_str << "updateRoutingDataBase"  << "," << std::to_string(profile_time.getTimer()) << std::endl;



    if(db::setting.debug){
        std::string file_name_csv = db::setting.outputFile + "_legalizer_critical_"+std::to_string(iter)+".csv";
        std::ofstream fout(file_name_csv);
        fout << database.ss_legalizer.str();
        fout.close(); 
    }

}//end runMT


void Placer::initCellsCandidate(){
    if(db::setting.debug){
        log() << "cells size: " << cells.size() << std::endl;
    }
    runJobsMT(cells.size(), [&](int i) { cells[i].initCandidates(); });
    // if(database.policy_set.find("addMedianCellCandidates") == database.policy_set.end()){
    //     if(db::setting.debug){
    //         log() << "addMedianCellCandidates..." << std::endl;
    //     }
    //     runJobsMT(cells.size(), [&](int i) { cells[i].addMedianCellCandidates(); });
    //     if(db::setting.debug){
    //         log() << "end addMedianCellCandidates..." << std::endl;
    //     }
    // }
        
    // runJobsMT(cells.size(), [&](int i) { cells[i].addRandomNetBoxCellCandidates(); });
    // // runJobsMT(cells.size(), [&](int i) { cells[i].addMedianCellCandidatesSiteBased(); });
    // // runJobsMT(cells.size(), [&](int i) { cells[i].addGCellBasedCandidates(); });
    // // runJobsMT(cells.size(), [&](int i) { cells[i].addWireDependentCellCandidates(); });
    // // log() << "addPinAlignmentCellCandidates..." << std::endl;
    // // runJobsMT(cells.size(), [&](int i) { cells[i].addPinAlignmentCellCandidates(); });
    // // log() << "end addPinAlignmentCellCandidates..." << std::endl;
    // // log() << "addManualCellCandidates..." << std::endl;
    // // runJobsMT(cells.size(), [&](int i) { cells[i].addManualCellCandidates(); });
    // // log() << "end addManualCellCandidates..." << std::endl;
    // // log() << "addLegalizedCellCandidates..." << std::endl;
    // // runJobsMT(cells.size(), [&](int i) { cells[i].addLegalizedCellCandidates(); });
    // // log() << "end addLegalizedCellCandidates..." << std::endl;
    // // log() << "addILPLegalizedCellCandidates..." << std::endl;
    // // runJobsMT(cells.size(), [&](int i) { cells[i].addILPLegalizedCellCandidates(); });
    // // log() << "end addILPLegalizedCellCandidates..." << std::endl;
    // // log() << "addBacktrackingLegalizedCellCandidates..." << std::endl;
    // // runJobsMT(cells.size(), [&](int i) { cells[i].addBacktrackingLegalizedCellCandidates(); });
    // // log() << "end addBacktrackingLegalizedCellCandidates..." << std::endl;
    
    
    if(database.policy_set.find("addILPLegalizedCellCandidates") == database.policy_set.end()){
        if(db::setting.debug){
            log() << "addILPLegalizedCellCandidates..." << std::endl;
        }
        runJobsMT(cells.size(), [&](int i) { cells[i].addILPLegalizedCellCandidates(); });

        // #pragma omp parallel for
        // for(int i = 0; i < cells.size(); i++){
        //     cells[i].addILPLegalizedCellCandidates(); 
        // }
        if(db::setting.debug){
            log() << "end addILPLegalizedCellCandidates..." << std::endl;
        }
    }
        


    // for (auto cell : cells){
    //     log() << "cell: " << cell.getName() << std::endl;
    //     auto ovrlp_pair_cells = cell.ovrlp_cells_;
    //     for(auto pair : ovrlp_pair_cells ){
            
    //         log() << "cellOvrlp: " << database.cells[pair.first].getName() 
    //               << ", tg box: " << pair.second << std::endl;
            
    //     }
    // }
    // return;
    bool debug = false;
    for (auto cell : cells){
        if(cell.getName() == "g61949_u3"){
            log() << "cell: " << cell.getName() << std::endl;
            debug = true;
        }else{
            debug = false;
        }
            
        std::set<int> seen_cells;
        auto cell_name = cell.getName();
        auto ovrlp_pair_cells = cell.ovrlp_cells_;
        for(auto pair : ovrlp_pair_cells ){
            if(seen_cells.find(pair.first) == seen_cells.end() && 
                (database.cells[pair.first].getName() != cell_name)){
                // log() << "cellOvrlp: " << database.cells[pair.first].getName() 
                //   << ", tg box: " << pair.second << std::endl;
                seen_cells.insert(pair.first);
                // auto cell_ovrlp = database.cells[pair.first];
                auto& cell_ovrlp = cells[pair.first];
                // cell_ovrlp.initCandidates();
                if(debug || cell_ovrlp.getName() == "g61949_u3"){

                    log() << "overlap_cell: " << cell_ovrlp.getName() 
                          << pair.second 
                          << ", in cell: " << cell.getName() << std::endl;

                }
                    
                cell_ovrlp.addOvrlpCandidates(pair.second);
                // cells.push_back(cell_ovrlp);

            }
            
        }
    }//end for 

    

}//end initCellsCandidate



Placer::Placer(CongestionMap& congMap,vector<vector<int>>& routeTableTmp,int iter_t,utils::timer& profile_time_,
        std::stringstream& profile_time_str_) 
    : congMap(congMap), routeTable(routeTableTmp), iter(iter_t)
    , profile_time(profile_time_),profile_time_str(profile_time_str_) {
    // readLUT();  // read flute LUT    
    // for (auto& net : grDatabase.nets){
    //     netsToIdx[net.getName()] =net.dbNet.idx;
    // } 
    // initRtree();
    debug_global = false;
}

void Placer::checkGlobalRoutingResults(){
    // std::set<int> nets_survey;
    // for(int i = 0; i < routeTable.size(); i++){
    //     auto net_idx = database.nets[i].idx;
    //     // ss << net.getName() << ",";
    //     for(int j = 0 ; j < routeTable[i].size(); j++){
    //         auto size = routeTable[i].size()-1;
    //         if(routeTable[i][size] - routeTable[i][0] > 0){
    //             nets_survey.insert(net_idx);
    //         }
    //     }
    // }//end for 

    // log() << "critical nets which I found..." << std::endl;
    // for(auto net_idx : nets_survey){        
    //     log() << "net name: " << database.nets[net_idx].getName() << std::endl;
    //     grDatabase.nets[net_idx].getCriticalCell();
    // }

    std::vector<gr::GrNet> gr_nets;
    for(auto gr_net : grDatabase.nets){        
        gr_nets.push_back(gr_net);
        // log() << "net name: " << gr_net.getName() << std::endl;
        // gr_net.checkRouteHist();
    }
    bool debug = false;
    // if(iter == 0)
    runJobsMT(gr_nets.size(), [&](int i) { gr_nets[i].checkRouteHist(); });
    runJobsMT(gr_nets.size(), [&](int i) { gr_nets[i].checkObidenceGrVsDr(); });
    // runJobsMT(10, [&](int i) { gr_nets[i].checkObidenceGrVsDr(); });
    // else{
    //     std::vector<int> nets_tmp;
    //     for (auto& net : grDatabase.nets) nets_tmp.push_back(net.dbNet.idx);
    //     sort(nets_tmp.begin(), nets_tmp.end(), [&](int id1, int id2) {
    //         return grDatabase.nets[id1].getWirelength() < grDatabase.nets[id2].getWirelength();
    //     });
    //     double percentage = 0.1;
    //     for(int i =0; i < percentage*nets_tmp.size();i++){
    //         if(debug)
    //             log() << "nets_tmp" << gr_nets[nets_tmp[i]].getName() << std::endl;
    //         gr_nets[nets_tmp[i]].is_critical_net = true;
    //     }
    // }
    //  if(debug){
    //     log() << "new critical nets ... " << std::endl;
    //     for(int i = 0; i < gr_nets.size(); i++){        
    //         log() << "net: " << gr_nets[i].getName() 
    //                 << ", wl: " << gr_nets[i].getWirelength()
    //                 << ", is_critical: " << gr_nets[i].is_critical_net << std::endl;
    //     }        
    // }

    for(int i = 0; i < gr_nets.size(); i++){        
        grDatabase.nets[i].is_critical_net = gr_nets[i].is_critical_net; 
    }
    gr_nets.clear();

    // for(auto& gr_net : grDatabase.nets){        
    //     // log() << "net name: " << gr_net.getName() << std::endl;
    //     gr_net.checkRouteHist();
    // }

    // log() << "critical nets which I found..." << std::endl;
    // std::set<int> critical_cells_set;
    for(auto& gr_net : grDatabase.nets){   
        if(gr_net.is_critical_net){
            // log() << "net name: " << gr_net.getName() << std::endl;
            auto cell_idx = gr_net.getCriticalCell();
            database.cells[cell_idx].is_critical_cell = true;
            // grDatabase.critical_cells.insert(cell_idx);
            // log() << "critical cell: " << database.cells[cell_idx].getName() << std::endl;
            // break;
        }
        // gr_net.getCriticalCell();
    }

    // for(auto idx : critical_cells)
    //     grDatabase.critical_cells.insert(idx);



    // log() << "checkGlobalRoutingResults ... " << std::endl;
    // for(auto net : grDatabase.nets){
    //     log() << "net: " << net.getName() << std::endl;
    //     log() << "topo size: " << net.gridTopo.size() << std::endl;
    //     log() << "path cost: " << net.getPathCost() << std::endl;
    // }
    
}//end checkGlobalRoutingResults

void Placer::initCellsUpdateFeatures(){
    bool debug = false;

    for(int i = 0; i < database.cells.size() ; i++){
        cells.push_back(database.cells[i]);        
    }//end loop 
    if(db::setting.debug){
        log() << "update updateWirelengthHist" << std::endl;
    }
    for(auto& grNet : grDatabase.nets){
        if(database.policy_set.find("updateWirelengthHist") == database.policy_set.end())
            grNet.updateWirelengthHist();
        if(database.policy_set.find("updateMaxEdgeCost") == database.policy_set.end())
            grNet.updateMaxEdgeCost();
        if(database.policy_set.find("updateNetFeatures") == database.policy_set.end())
            grNet.updateNetFeatures();
    }
    if(database.policy_set.find("calcNetsHPWL") == database.policy_set.end())
        calcNetsHPWL();
    
    
    if(db::setting.debug){
        log() << "end updateWirelengthHist" << std::endl;
    }

    // runJobsMT(cells.size(), [&](int i) { cells[i].updateMedianCell(); });
    if(database.policy_set.find("updateCellFeatures") == database.policy_set.end())
        runJobsMT(cells.size(), [&](int i) { cells[i].updateCellFeatures(); });
   

    for(auto cell : cells){
        database.cells[cell.idx].median_x_y_cell = cell.median_x_y_cell;
        database.cells[cell.idx].median_x_y_pin = cell.median_x_y_pin;
        database.cells[cell.idx].feature = cell.feature;
    }//end for


    // bool log_features = false;

    // if(log_features)
    //     logCellsFeature();


}//end initCells


void Placer::selectCriticalCellsNetDegree(){

}//end selectCriticalCellsNetDegree

void Placer::selectCriticalCellsISPD(){
    std::set<std::string> critical_cells_set;

}//end selectCriticalCellsISPD

void Placer::selectCriticalCells(){    
    // std::set<std::string> critical_cells_set;//= {"inst4132","inst5195",
                                             // "inst5333","inst2591","inst6050", "inst6458"};
    // std::set<std::string> critical_cells_set;//= {
        // // "inst8861", 
        // "inst6050"}; 
        // "inst7112", 
        // "inst6299", 
        // "inst5320", 
        // "inst5092", 
        // "inst150 ", 
        // "inst87  " };
                                                
    vector<std ::string> filter_cells_name;
    boost::split(filter_cells_name, db::setting.critical_cells_set, boost::is_any_of(","));



    for(auto cell_name : filter_cells_name){
       critical_cells_set.insert(cell_name);
    }

    // for(auto cell_tmp : critical_cells_set){
    //     std::cout << "cell_tmp: " << cell_tmp << std::endl;
    // }

    // std::set<std::string> critical_cells_set = {"inst4132","inst5195","inst2591","inst6050"};
    // std::set<std::string> critical_cells_set = {"inst4132"};
    // checkGlobalRoutingResults();
    // // for(int i = 0; i < database.cells.size() ; i++){
    // //     critical_cells_set.insert(database.cells[i].getName());
    // // }
    if(db::setting.critical_cells_set == "-1") 
        schellingSegerationModel(critical_cells_set);
    // MWIS(critical_cells_set);

    // for(auto pair : database.suspected_cells_dict){
    //     auto cell_idx = pair.first;
    //     critical_cells_set.insert(database.cells[cell_idx].getName());
    // }

    // log() << "size of critical_cells_set: " << critical_cells_set.size() << std::endl;
    // for(auto cell_name : critical_cells_set){
    //     log() << "cell: " << cell_name << std::endl;
    // }
    // return;


    //------------ Temporary Critical Cells ----------------
    // 1- Scenario#1: constant set of critical cells for debugging
    // std::set<std::string> critical_cells_set = {"inst8059"};
    // 2- Scenario#2: randomly select critical_Cells
    // std::random_device dev;
    // std::mt19937 rng(dev());
    // // // std::set<std::string> critical_cells_set;
    // for(int i = 0; i < database.cells.size()/10; i++){
    //     std::uniform_int_distribution<std::mt19937::result_type> rnd_fun(0,database.cells.size()-1); // distribution in range [0, cell set]
    //     // log() << "random number: " << rnd_fun(rng); << std::endl;
    //     // auto cell_name = database.cells[rnd_val].getName();
    //     critical_cells_set.insert(database.cells[rnd_fun(rng)].getName());
    // }
    // 3- Scenario#3: All cells one-by-one
        // std::set<std::string> critical_cells_set;
    // critical_cells_set.insert("inst4132");
 
    // log() << "critical_cells_set: " << std::endl;
    // for(auto cell_name : critical_cells_set){
    //     log() << "critical_cell: " << cell_name << std::endl;
    // }

    

    for(int i = 0; i < database.cells.size() ; i++){
        auto cell_name = database.cells[i].getName();
        if(critical_cells_set.find(cell_name) == critical_cells_set.end()){
            database.cells[i].is_critical_cell = false;
            cells[i].is_critical_cell = false;
        }else{
            database.cells[i].is_critical_cell = true;
            cells[i].is_critical_cell = true;
        }
    }



    // for(int i : critical_cells){
    //     cells.push_back(database.cells[i]);
    // }


    // for(int i : grDatabase.critical_cells){
    //     if(!cells[i].isFixed){
    //         auto connected_cells_idx = cells[i].connectd_cells;
    //         for(auto cell_idx : connected_cells_idx){
    //             // log() << "connected: " << cells[cell_idx].getName() << std::endl;
    //             cells[cell_idx].isFixed = true;
    //             database.cells[cell_idx].isFixed = true;
    //         }

    //     }
    // }

    bool debug = false;
    // for(auto cell : cells){
    //     log() << "cell: " << cell.getName() 
    //           << ", critical: " << cell.is_critical_cell << std::endl;
    // }


    if(debug)
        log() << "critical cells iter: " << iter << std::endl;
    for(int i = 0; i < database.cells.size() ; i++){
        if(debug)
        {
            
            if(database.cells[i].is_critical_cell)
                log() << "cell: " << cells[i].getName() << std::endl;
        }
            
        if(!cells[i].isFixed && cells[i].is_critical_cell){
            auto connected_cells_idx = cells[i].connectd_cells;
            for(auto cell_idx : connected_cells_idx){
                // log() << "connected: " << cells[cell_idx].getName() << std::endl;
                if(cells[i].getName() != cells[cell_idx].getName()){
                    cells[cell_idx].isFixed = true;
                    database.cells[cell_idx].isFixed = true;
                }
            }

        }
    }

    if(debug){
        
        for(auto cell : cells){
            log() << "cell: " << cell.getName()
                  << "box: " << cell.getCellBox()
                  << ", isFixed: " << cell.isFixed
                  << ", isCritical: " << cell.is_critical_cell << std::endl;
        }

    }


    if(db::setting.debug){
        int num_critical_cells = 0;
        int num_fixed_cells = 0;
        int num_movable_cells = 0;
        std::stringstream ss;
        ss << "inst,lx,ly,hx,hy,status" << std::endl;

        for(auto cell : cells){
            int status = 0;
            if(cell.isFixed){
                status=1;
                num_fixed_cells++;
            } 
            if(cell.is_critical_cell) {
                num_critical_cells++;
                status=2;
            }

            auto box = cell.getCellBox();
            ss << cell.getName() 
               << "," << box.lx()
               << "," << box.ly()
               << "," << box.hx()
               << "," << box.hy()
               << "," << status
                << std::endl;
        }

        std::string file_name_csv = db::setting.outputFile + "_placer_criticalCells_"+std::to_string(iter)+".csv";
        std::ofstream fout(file_name_csv);
        fout << ss.str();
        fout.close();

        log() << "num_critical_cells: " << num_critical_cells
              << ", num_fixed_cells: "  << num_fixed_cells
              << ", num_movable_cells: "  << cells.size() - num_fixed_cells << std::endl;
    }


    bool log_critical=false;

    if(log_critical)
        logCellsCritical();
}//selectCriticalCells

void Placer::updateRoutingDataBase(std::vector<int>& netsToRoute){
    if(db::setting.debug){
        log() << "################################################################" << std::endl;
        log() << "update routing database ..." << std::endl;
    }
    auto rsynService = database.getRsynService();

    // database.clear();
    // update pin access box of the net
    for (auto idx : netsToRoute){
        auto& net = database.nets[idx];
        net.updatePinAccessBoxes(rsynService);
    }



    // ----------------
    database.initRtrees(); 
    database.markPinAndObsOccupancy();
    database.initMTSafeMargin();
    grDatabase.update();

    // log() << "grDatabase this: " << &grDatabase << std::endl;

    for(auto& gr_net : grDatabase.nets){
        gr_net.is_critical_net = false;
    }
    for(auto& cell : database.cells){
        cell.is_critical_cell = false;
        cell.isFixed = false;
    }


    
}//end 

void Placer::writeToCSV(){
    logInstances();
    // logCongestion();
}//end writeToCSV

void Placer::logCongestion(){
    congMap.logCSV("log_cong.csv");

    // get edgecost and via cost
    auto numPointsX = ceil(grDatabase.getNumGrPoint(X) / (double)congMap.getCellWidth());
    auto numPointsY = ceil(grDatabase.getNumGrPoint(Y) / (double)congMap.getCellHeight());
    int numPoints = numPointsX * numPointsY * database.getLayerNum();

       
    for (int l = 0; l < database.getLayerNum() - 1; l++) {
        for (int x = 0; x < numPointsX; x++) {
            for (int y = 0; y < numPointsY; y++) {
                log() << "l: " << l << ",x: " 
                      << x << ",y: " << y << ", resourceUsage: "
                      << congMap.getRsrcUsage(l, x, y) << std::endl;
            }
        }
    }


    // 3. Add inter-layer connection
    for (int l = 0; l < database.getLayerNum() - 1; l++) {
        for (int x = 0; x < numPointsX; x++) {
            for (int y = 0; y < numPointsY; y++) {
                log() << "l: " << l << ",x: " 
                      << x << ",y: " << y << ", viaCost: "
                      << congMap.getCrsnViaCost({l, x, y}) << std::endl;
            }
        }
    }

    // 4. add intra-layer connection
    for (int l = 0; l < database.getLayerNum(); l++) {
        if (database.getLayerDir(l) == X) {
            for (int x = 0; x < numPointsX; x++) {
                for (int y = 0; y < numPointsY - 1; y++) {
                    log() << "l: " << l << ",x: " 
                      << x << ",y: " << y << "," << y+1 <<", edgeCost: "
                      << congMap.getCrsnEdgeCost({l, x, y}, {l, x, y + 1}) << std::endl;
                }
            }
        } else {
            for (int x = 0; x < numPointsX - 1; x++) {
                for (int y = 0; y < numPointsY; y++) {
                    log() << "l: " << l << ",x: " 
                      << x << ", " << x+1 <<",y: " << y << ", edgeCost: "
                      << congMap.getCrsnEdgeCost({l, x, y}, {l, x + 1, y}) << std::endl;
                }
            }
        }
    }



}

void Placer::logInstances(){
    auto rsynService = database.getRsynService();
    // write instances in csv file
    std::stringstream ss;
    ss << "inst_name,x,y,width,height" << std::endl;
    for (Rsyn::Instance instance : rsynService.module.allInstances()) {
        // if (instance.getName() != "inst5638") continue;
        if (instance.getType() != Rsyn::CELL) continue;
        // phCell
        Rsyn::Cell cell = instance.asCell();
        Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cell);
        Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cell);
        ss << instance.getName() << ", " << instance.getX() << ", "
            << instance.getY() << ", " << phLibCell.getWidth() << ", " 
            << phLibCell.getHeight() << std::endl;
    }
    std::string file_name_csv = "instances.csv";
    std::ofstream fout(file_name_csv);
    fout << ss.str();
    fout.close();
}

void Placer::refinePlacement(std::vector<int>& netsToRoute){
    // logCellLocations();
    std::stringstream ss;
    bool log_moved=false;

    if(log_moved){
        ss << "cell_name,init_x,init_y,new_x,new_y" << std::endl;
    }
    
    std::string file_name = db::setting.outputFile+ ".cell_"+std::to_string(iter) + ".moved.csv";
    auto rsynService = database.getRsynService();
    std::set<int> netsToRoute_set;
    int i = 0;
    bool debug = false;
    // log() << "cellsToMove: " << cellsToMove.size() << std::endl;
    database.total_cells_moved += cellsToMove.size();

    int max_disp_x = 0;
    int total_disp_x = 0;
    int max_disp_y = 0;
    int total_disp_y = 0;

    
    for(auto cell_box_pair: cellsToMove){
        auto cell = database.cells[cell_box_pair.first];
        auto box = cell_box_pair.second;
        database.critical_cells_hist_moved.insert(cell.idx);

        if(db::setting.debug){
            // if(!(cell.getName() == "inst2977" || cell.getName() == "inst2916")) continue;
            log() << "cell: " << cell.getName() 
            
                << ", is moved from " << cell.rsynInstance.getBounds()
                << ", to: " << box << std::endl;
        }
        if(log_moved){
            ss  << cell.getName() 
                << "," << cell.rsynInstance.getBounds().getLower().x/database.libDBU
                << "," << cell.rsynInstance.getBounds().getLower().y/database.libDBU
                << "," << box.lx()/database.libDBU
                << "," << box.ly()/database.libDBU << std::endl;
        }


        int disp_x =  std::abs(cell.rsynInstance.getBounds().getLower().x - box.lx());
        int disp_y =  std::abs(cell.rsynInstance.getBounds().getLower().y - box.ly());
        total_disp_x = disp_x + total_disp_x;
        total_disp_y = disp_y + total_disp_y;


        if(max_disp_x < disp_x){
            max_disp_x = disp_x;
        }
        if(max_disp_y < disp_y){
            max_disp_y = disp_y;
        }

        int eps = 10;

        Rsyn::Cell cellRsyn = cell.rsynInstance.asCell();
        Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cellRsyn);
        Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cellRsyn);

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

            // cell orientation must change to row orientation
            // the current version is not right
            // rsynService.physicalDesign.placeCell(phCell,
            //             cell_box_pair.second.lx(),
            //             cell_box_pair.second.ly(),
            //             cellRsyn.getOrientation());
            rsynService.physicalDesign.placeCell(phCell,
                        cell_box_pair.second.lx(),
                        cell_box_pair.second.ly(),
                        row.getSiteOrientation());


            // break; 
        }//end for

        
        for (Rsyn::Pin pin : cell.rsynInstance.allPins(false)) {
            auto net = pin.getNet();
            if(net){
                int idx = database.netsToIdx[net.getName()];
                netsToRoute_set.insert(idx);
                // log() << "nets to reroute after movement " << net.getName() << ", "
                //     << idx << std::endl;
            }
            
        }//end for allpins
    }//end for 

    for(auto tmp : netsToRoute_set){
        netsToRoute.push_back(tmp);
    }

    int num_pin = 0;
    for(auto net_idx : netsToRoute){
        num_pin = num_pin + database.nets[net_idx].numOfPins();
    }
    int avg_connected_cells= 0;
    for(auto cell_box_pair: cellsToMove){
        auto cell = database.cells[cell_box_pair.first];
        // int connected_cells = cell.feature.degree_connected_cells;
        avg_connected_cells = avg_connected_cells 
                            + cell.feature.degree_connected_cells;

    }

    avg_connected_cells = avg_connected_cells/cellsToMove.size();
    float avg_disp_x = float(total_disp_x)/float(cellsToMove.size());
    float avg_disp_y = float(total_disp_y)/float(cellsToMove.size());


    auto dieBound = rsynService.physicalDesign.getPhysicalDie().getBounds();
    DBU die_lx = dieBound.getLower().x;
    DBU die_ly = dieBound.getLower().y;
    DBU die_hx = dieBound.getUpper().x;
    DBU die_hy = dieBound.getUpper().y;

    if(log_moved){
        std::ofstream fout(file_name);
        fout << ss.str();
        fout.close();
    }//end log moved

    try{
        log() << "nets characterestics that are moved"<<std::endl;
        log() << "num_nets: " << netsToRoute.size() 
            << " out of " << database.nets.size() 
            << " ( " << float(netsToRoute.size()*100)/float(database.nets.size())
            << " %) " << std::endl;
    }catch (const std::exception& e){
        log() << "divide by zero!" << std::endl;
    }

    try{
        if(netsToRoute.size() != 0){
            log() << "average num pins: " << num_pin/netsToRoute.size() << std::endl;
            log() << "max_disp_x: " << max_disp_x 
                << " die width " << std::abs(die_lx-die_hx) 
                << " ( "  << float(max_disp_x*100)/float(std::abs(die_lx-die_hx))
                << " %) " << std::endl;
        }
        
    }catch (const std::exception& e){
        log() << "divide by zero!" << std::endl;
    }
    
    try{
        log() << "max_disp_y: " << max_disp_y 
            << " die height " << std::abs(die_ly-die_hy) 
            << " ( "  << float(max_disp_y*100)/float(std::abs(die_ly-die_hy))
            << " %) " << std::endl;
    }catch (const std::exception& e){
        log() << "divide by zero!" << std::endl;
    }

    try{
        log() << "avg_disp_x: " << avg_disp_x
            << " ( " << (avg_disp_x*100)/float(std::abs(die_lx-die_hx))
            << " %) "<< std::endl;
    }catch (const std::exception& e){
        log() << "divide by zero!" << std::endl;
    }
        
    try{
        log() << "avg_disp_y: " << avg_disp_y
            << " ( " << (avg_disp_y*100)/float(std::abs(die_ly-die_hy))
            << " %) "<< std::endl;
    }catch (const std::exception& e){
        log() << "divide by zero!" << std::endl;
    }

    try{
        log() << "avg cells degree: " << avg_connected_cells << std::endl;
        log() << "TotalCellsToMove: " << cellsToMove.size()
            << " out of " << database.cells.size() 
            << " ( "  << float(cellsToMove.size()*100)/float(database.cells.size())
            << " %) " << std::endl;
    }catch (const std::exception& e){
        log() << "divide by zero!" << std::endl;
    }

    

}//end refineplacement


void Placer::calcCellsCostMT(){
    runJobsMT(cells.size(), [&](int i) { cells[i].calcCostPlacementCandidates(); });

    bool log_candidates = false;

    if(log_candidates)
        logCellsCandidates();
}


void Placer::findBestEntryMT(){
    bool debug = false;
    std::vector<double> weights;
    // idx of weight that have conflict with each other.
    std::vector<std::vector<int>> conflict_matrix;
    std::vector<std::vector<int>> conflict_matrix_connected;
    std::vector<std::vector<int>> legalize_matrix;

    //idx weight to idx cells and new box position
    std::unordered_map<int,std::pair<int,utils::BoxT<DBU>>> weightIdxToCell;
    // cell index to weight indexes
    std::unordered_map<int,std::vector<int>> cellIdxToWeightIdxs;
    int i = 0;
    // log() << "construct ILP model ..." << std::endl;
    // construct ILP model and find add constraints of having one solution
    // for each cell. 
    for(auto& cell : cells){
        std::vector<int> each_cell_has_one_location_constraint;
        // if(cell.placement_candidates.size() <= 1 ) continue;

        for(auto candidate_pair : cell.placement_candidates){
            auto box = candidate_pair.first;
            auto cost = candidate_pair.second;
            weights.push_back(cost);
            auto idx = cell.idx;
            weightIdxToCell[i] = std::make_pair(idx,box);
            if(cellIdxToWeightIdxs.find(idx) != cellIdxToWeightIdxs.end()){
                auto& vec_weight_idxs = cellIdxToWeightIdxs[idx];
                vec_weight_idxs.push_back(i);
            }else{
                std::vector<int> vec_weight_idxs;
                vec_weight_idxs.push_back(i);
                cellIdxToWeightIdxs[idx] = vec_weight_idxs;
            }
            each_cell_has_one_location_constraint.push_back(i);
            if(debug)
                log() << "ILPidx: " << i <<", cell: " << cell.getName() << ", cell.candidate: " << box
                    << ", cost: " << cost << std::endl;
            i++;
            
        }
        conflict_matrix.push_back(each_cell_has_one_location_constraint);
    }//end for 

    // add connected cells constraint
    // for(auto cell : database.cells){
    //     std::vector<int> one_cell_each_net_constraint;
    //     auto connectd_cells = cell.connectd_cells;
    //     for(auto connected_cell_idx : connectd_cells){
    //         auto vec_weight_idxs = cellIdxToWeightIdxs[connected_cell_idx];
    //         for(auto weight_idx: vec_weight_idxs){
    //             one_cell_each_net_constraint.push_back(weight_idx);
    //         }
    //     }

    //     auto vec_weight_idxs = cellIdxToWeightIdxs[cell.idx];
    //     for(auto weight_idx: vec_weight_idxs){
    //         one_cell_each_net_constraint.push_back(weight_idx);
    //     }
    //     conflict_matrix_connected.push_back(one_cell_each_net_constraint);
    // }

    // for(auto net : database.nets){
    //     std::vector<int> one_cell_each_net_constraint;
    //     std::vector<int> cells_idx;
    //     database.getNetCells(net.idx,cells_idx);
    //     for(auto cell_idx : cells_idx){
    //         auto vec_weight_idxs = cellIdxToWeightIdxs[cell_idx];
    //         for(int i = 1; i < vec_weight_idxs.size(); i++){
    //             one_cell_each_net_constraint.push_back(vec_weight_idxs[i]);
    //         }
    //     }
    //     if(one_cell_each_net_constraint.size() > 1)
    //         conflict_matrix_connected.push_back(one_cell_each_net_constraint);
    // }//end for 


    //---------------------------------------------
    // legalization constraints
    RTree local_candidates_rtree; 
    for(auto pair_idx_cell : weightIdxToCell){
        auto vars_idx = pair_idx_cell.first;
        auto cellIdxToBox_pair = pair_idx_cell.second;
        boostBox box(boostPoint(cellIdxToBox_pair.second.lx(),
                                cellIdxToBox_pair.second.ly()),
                     boostPoint(cellIdxToBox_pair.second.hx(),
                                cellIdxToBox_pair.second.hy()));
        local_candidates_rtree.insert({box, vars_idx});
    }//end for 

  
    for(auto pair_idx_cell : weightIdxToCell){
        auto vars_idx = pair_idx_cell.first;
        auto cellIdxToBox_pair = pair_idx_cell.second;
        auto cell_idx = cellIdxToBox_pair.first;
        auto box = cellIdxToBox_pair.second;

        // std::vector<int> cells_ovrlp_constraint;
        // cells_ovrlp_constraint.push_back(vars_idx);
        int eps = 10;
        boostBox rtreeQueryBox(boostPoint(box.lx()+eps,
                                          box.ly()+eps),
                               boostPoint(box.hx()-eps,
                                          box.hy()-eps));
        vector<std::pair<boostBox, int>> queryResults;
        local_candidates_rtree.query(bgi::intersects(rtreeQueryBox), std::back_inserter(queryResults));
        for (const auto& queryResult : queryResults) {
            std::vector<int> cells_ovrlp_constraint;
            cells_ovrlp_constraint.push_back(vars_idx);

            auto conf_vars_idx = queryResult.second;
            auto conf_cell_idx = weightIdxToCell[conf_vars_idx].first;
            // log() << "cell_idx: " << cell_idx
            //       << ", conf_vars_idx: " << conf_vars_idx 
            //       << ", conf_cell_idx: " << conf_cell_idx
            //       << std::endl ;

            if( cell_idx != conf_cell_idx){
                cells_ovrlp_constraint.push_back(conf_vars_idx);
                legalize_matrix.push_back(cells_ovrlp_constraint);
            }
        }//end for 

        // log() << "cells_ovrlp_constraint: " << cells_ovrlp_constraint.size() << std::endl;

        // if(cells_ovrlp_constraint.size() > 1)
        //     legalize_matrix.push_back(cells_ovrlp_constraint);
    }

    // 
    // log() << "legalize matrix: " << std::endl;
    // for(auto tmp_row : legalize_matrix){
    //     log() << "new tmp: " << std::endl;
    //     for(auto tmp : tmp_row) {
    //         log() << "tmp: " << tmp << std::endl;
    //     }
        
    // }

    //---------------------------------------------


    std::vector<int> sols; 
    std::string model_name = "model_" + std::to_string(iter) + ".lp";
    ILPSolver(weights,weightIdxToCell,conflict_matrix,conflict_matrix_connected,legalize_matrix,sols,false,model_name);

    // log() << "ILPSolver solution..." <<std::endl;
    for(auto idx : sols){
        auto cell_box_pair = weightIdxToCell[idx];

        auto old_box = utils::BoxT<DBU>(
            database.cells[cell_box_pair.first].rsynInstance.getBounds().getLower().x,
            database.cells[cell_box_pair.first].rsynInstance.getBounds().getLower().y,
            database.cells[cell_box_pair.first].rsynInstance.getBounds().getUpper().x,
            database.cells[cell_box_pair.first].rsynInstance.getBounds().getUpper().y
        );

        if(!old_box.isSame(cell_box_pair.second)){
            if(debug){
                log() << "sol: " << idx << ", " 
                      << database.cells[cell_box_pair.first].getName() 
                      << ", " << cell_box_pair.second << std::endl;
            }
            cellsToMove.push_back(std::make_pair(cell_box_pair.first,cell_box_pair.second));
        }//end if 

    }//end for 



}//end findBestEntry


void Placer::getOvrlpConflicts(
        std::vector<CandWeight>& weights
    ,   std::vector<std::vector<int>>& ovrlps
){

    double eps = 0.01;

    bool debug = false || debug_global;
    
    std::set<std::pair<int,int>> rs_loc;
    RTree cand_rtree;

    for(int i = 0 ; i < weights.size() ; i++){
        auto box_cell = weights[i].box;
        boostBox box(boostPoint(box_cell.lx()+eps,box_cell.ly()+eps),\
                     boostPoint(box_cell.hx()-eps,box_cell.hy()-eps));
        cand_rtree.insert({box, i});
    }//end for

    if(debug){
        std::stringstream ss;

        ss << "w_idx,inst,r,s"<<std::endl;

        for(auto cand : cand_rtree){
            const auto& b = cand.first;
            auto idx = cand.second;
            auto name = database.cells[weights[idx].cell_idx].getName();
                
            auto box = utils::BoxT<DBU>(bg::get<bg::min_corner, 0>(b),
                                        bg::get<bg::min_corner, 1>(b),
                                        bg::get<bg::max_corner, 0>(b),
                                        bg::get<bg::max_corner, 1>(b));


            ss << std::to_string(idx)
            << "," << name 
            << "," << std::to_string(box.lx())
            << "," << std::to_string(box.ly())
            << std::endl;
        
        }


        std::string file_name_csv = db::setting.outputFile + "_legRtree_placer.csv";
        std::ofstream fout(file_name_csv);
        fout << ss.str();
        fout.close();
    }
    

  


    std::set<std::string> encode_idx;
    std::set<std::pair<int,int>> conflict_set;

    for(int i = 0; i < weights.size() ; i++)
    {
        
        // conflict_set.insert(i);
        auto box_cell = weights[i].box;
        int cell_idx = weights[i].cell_idx;
        std::string cell_name = database.cells[cell_idx].getName();


        double xl = box_cell.lx()+eps;
        double yl = box_cell.ly()+eps;
        double xh = box_cell.hx()-eps;
        double yh = box_cell.hy()-eps;

        boostBox rtreeQueryBox(
        boostPoint(xl,
                   yl),
        boostPoint(xh,
                   yh));
        vector<std::pair<boostBox, int>> queryResults;


        auto box_query = utils::BoxT<DBU>(bg::get<bg::min_corner, 0>(rtreeQueryBox),
                                    bg::get<bg::min_corner, 1>(rtreeQueryBox),
                                    bg::get<bg::max_corner, 0>(rtreeQueryBox),
                                    bg::get<bg::max_corner, 1>(rtreeQueryBox));
        if(debug){
            // log() << "i: " << i
            //   << ", cell_name: " << cell_name
            //   << ", xl: "    <<xl
            //   << ", yl: " << yl
            //   << ", xh: " << xh
            //   << ", yh: " << yh
            //   << ", eps: " << eps << std::endl;
        }
        

        // log() << "box: " << box_query << std::endl;
        
        cand_rtree.query(bgi::intersects(rtreeQueryBox), std::back_inserter(queryResults));

        
        // std::vector<int> ovrlp_set;
        std::set<int> ovrlp_set;
        


        // if(queryResults.size() <= 1){
        //     ovrlp_vec.push_back({});
        //     continue;
        // }

        for (const auto& queryResult : queryResults) {
            const auto& b = queryResult.first;
            auto j = queryResult.second;
            
            
            auto box = utils::BoxT<DBU>(bg::get<bg::min_corner, 0>(b),
                                        bg::get<bg::min_corner, 1>(b),
                                        bg::get<bg::max_corner, 0>(b),
                                        bg::get<bg::max_corner, 1>(b));
            // CellOvrlp cell_ovrlp;
            // // ovrlp_set.insert(idx);
            // if(ovrlp_set.find(weights[idx].idx) == ovrlp_set.end()){
            //     ovrlp_set.insert(weights[idx].idx);
            
            int cell_idx1 = weights[j].cell_idx;
            auto box_cell1 = weights[j].box;
            std::string cell_name1 = database.cells[cell_idx1].getName();

            if(debug){
                // log() 
                //     << "cell_name: " << cell_name1
                //     << ", r: " <<r1
                //     << ", s: " << s1
                //     << ", w: " << w1
                //     << ", j: " << j
                //     << std::endl;
            }

            int max = std::max(i,j);
            int min = std::min(i,j);
            
            if(conflict_set.find(std::make_pair(min,max)) == conflict_set.end()){
                std::vector<int > ovrlp_vec;
                if(cell_idx1 != cell_idx){
                    conflict_set.insert(std::make_pair(min,max));
                    ovrlp_vec.push_back(min);
                    ovrlp_vec.push_back(max);
                    ovrlps.push_back(ovrlp_vec);
                }//end if 
            }//end if 



            
        }//end for
        
       
    }

    

}//end getOvrlpConflicts

void Placer::ILPSolverV2(
                  std::vector<CandWeight>& weights
                , std::vector<std::vector<int>>& ovrlps
                , std::vector<std::vector<int>>& constraints
                , std::vector<int>& sol){
    bool debug = false || debug_global;
    
    vector<int> selected;
    // std::vector<int> sol;
    selected.resize(weights.size());
    std::string log_name = "placer_ilp_"+std::to_string(iter)+".lp";
    bool maxmimize = false;
    
    IloEnv env;

    try {
        IloModel model(env);
        IloIntArray weightsSet(env);
        IloIntVarArray vars(env);
        IloRangeArray con(env);


        for (auto tmp : weights) {
            weightsSet.add(tmp.cost);
            // tmp_idx++;
        }

        int counter = 0;
        

        for (int i = 0 ; i < weights.size(); i++) {
            vars.add(IloIntVar(env, 0, 1));

            if(debug){
                int cell_idx = weights[i].cell_idx; 
                auto box  = weights[i].box;
            

                std::string name = database.cells[cell_idx].getName() 
                                + "_" + std::to_string(box.lx())
                                + "_" + std::to_string(box.ly())
                                + "_" + std::to_string(box.hx())
                                + "_" + std::to_string(box.hy())
                                + "_" + std::to_string(i);
                vars[i].setName(name.c_str());
            }
            
        }

        // naming
        // // only one location for each cell constraint
        int j = 0;
        int counter_cons = 0;

        for (auto constraint_row : constraints) {
            IloNumExpr expr_constraint(env);
            for(auto var_idx : constraint_row){
                expr_constraint += vars[var_idx];    
            }  
            con.add(expr_constraint == 1);  
            if(debug){
                std::string const_name = "c_" + std::to_string(counter_cons);
                con[counter_cons].setName(const_name.c_str());
                counter_cons++;  
            }
            expr_constraint.end();    
            // if(j==100)break;
            // j++;
            // break;
        }


        for (auto conflict_row : ovrlps) {
            // if(j=!72){
            //     j++;
            //     continue;
            // }
            IloNumExpr expr_constraint(env);
            for(auto var_idx : conflict_row){
                expr_constraint += vars[var_idx];    
            }  
            con.add(expr_constraint <= 1);  
            if(debug){
                std::string const_name = "ovrlp_" + std::to_string(counter_cons);
                con[counter_cons].setName(const_name.c_str());
                counter_cons++;  
            }
            expr_constraint.end();    
            // if(j==100)break;
            // j++;
            // break;
        }  

        model.add(con);        

        IloNumExpr expression(env);
        for (int i = 0; i < weights.size(); i++) {
            expression += vars[i] * weightsSet[i];
        }
        if(!maxmimize){
            model.add(IloMinimize(env, expression));
        }else{
            model.add(IloMaximize(env, expression));
        }

        IloCplex cplex(model);
        if(!debug)
            cplex.setOut(env.getNullStream());


        cplex.setParam(IloCplex::Param::TimeLimit, db::setting.legalizerOptimzerTimeLimit);
        cplex.setParam(IloCplex::Param::Threads, db::setting.numThreads);
        
        // cplex.setOut(env.getNullStream());
        
        // if(db::setting.debugLegalizerExportILPModel)
        if(debug)
            cplex.exportModel(log_name.c_str());
        // return;
        
        
     
        if (!cplex.solve()) {
            log() << "Failed to optimize.\n";
            env.end();
            // std::exit(1);
            return;
        }
     


        IloNumArray vals(env);
    
        cplex.getValues(vals, vars);
        
        for (int i = 0; i < selected.size(); i++) {
            try {
                if(debug)
                    log() << "vals[" << i << "]: " << vals[i] << std::endl;
                if (vals[i] == 1) {
                    selected[i] = 1;
                    sol.push_back(i);
                }
            } catch (...) {
                log() << "Exception on index " << i << "\n";
            }
        }

        
    } catch (IloException& e) {
        log() << "Concert exception caught: " << e << endl;
    } catch (...) {
        log() << "Unknown exception caught" << endl;
    }
    env.end();                   

}//end ILPSolverV2


void Placer::findBestEntryMTV2(){
    crp::Solver solver;

    bool debug = false;
    std::vector<CandWeight> weights;
    std::vector<std::vector<int>> constraints;
    std::vector<std::vector<int>> ovrlps;
    std::vector<int> sols;

   
    int num_movable_cells = 0;
   
    // construct weights
    for(auto& cell : cells){
        if(cell.placement_candidates.size() > 1){
            num_movable_cells++;
        }//end if num_movable_cells
        int start_idx = weights.size();
        for(auto candidate_pair : cell.placement_candidates){
            auto box = candidate_pair.first;
            auto cost = candidate_pair.second;
            CandWeight cand_weight;
            cand_weight.cell_idx = cell.idx;

            cand_weight.box = box;// cell box
            cand_weight.cost = cost;
            weights.push_back(cand_weight);
        }        
        int stop_idx = weights.size();
        std::vector<int> single_position_cons;
        for(int i = start_idx; i < stop_idx;i++){
            single_position_cons.push_back(i);
        }
        constraints.push_back(single_position_cons);
        
    }//end for 

    logPlacementWeights(weights);

    getOvrlpConflicts(weights,ovrlps);

    log() << "num_movable_cells: " << num_movable_cells << std::endl; 


    ILPSolverV2(
                  weights
                , ovrlps
                , constraints
                , sols);


    for(auto idx : sols){
        auto cell_idx = weights[idx].cell_idx;
        auto new_box = weights[idx].box;

        auto old_box = utils::BoxT<DBU>(
            database.cells[cell_idx].rsynInstance.getBounds().getLower().x,
            database.cells[cell_idx].rsynInstance.getBounds().getLower().y,
            database.cells[cell_idx].rsynInstance.getBounds().getUpper().x,
            database.cells[cell_idx].rsynInstance.getBounds().getUpper().y
        );

        if(!old_box.isSame(new_box)){
            if(debug){
                log() << "sol: " << idx << ", " 
                      << database.cells[cell_idx].getName() 
                      << ", " << new_box << std::endl;
            }
            cellsToMove.push_back(std::make_pair(cell_idx,new_box));
        }//end if 

    }//end for

    // add connected cells constraint
    // for(auto cell : database.cells){
    //     std::vector<int> one_cell_each_net_constraint;
    //     auto connectd_cells = cell.connectd_cells;
    //     for(auto connected_cell_idx : connectd_cells){
    //         auto vec_weight_idxs = cellIdxToWeightIdxs[connected_cell_idx];
    //         for(auto weight_idx: vec_weight_idxs){
    //             one_cell_each_net_constraint.push_back(weight_idx);
    //         }
    //     }

    //     auto vec_weight_idxs = cellIdxToWeightIdxs[cell.idx];
    //     for(auto weight_idx: vec_weight_idxs){
    //         one_cell_each_net_constraint.push_back(weight_idx);
    //     }
    //     conflict_matrix_connected.push_back(one_cell_each_net_constraint);
    // }

    // for(auto net : database.nets){
    //     std::vector<int> one_cell_each_net_constraint;
    //     std::vector<int> cells_idx;
    //     database.getNetCells(net.idx,cells_idx);
    //     for(auto cell_idx : cells_idx){
    //         auto vec_weight_idxs = cellIdxToWeightIdxs[cell_idx];
    //         for(int i = 1; i < vec_weight_idxs.size(); i++){
    //             one_cell_each_net_constraint.push_back(vec_weight_idxs[i]);
    //         }
    //     }
    //     if(one_cell_each_net_constraint.size() > 1)
    //         conflict_matrix_connected.push_back(one_cell_each_net_constraint);
    // }//end for 


    // //---------------------------------------------
    // // legalization constraints
    // RTree local_candidates_rtree; 
    // for(auto pair_idx_cell : weightIdxToCell){
    //     auto vars_idx = pair_idx_cell.first;
    //     auto cellIdxToBox_pair = pair_idx_cell.second;
    //     boostBox box(boostPoint(cellIdxToBox_pair.second.lx(),
    //                             cellIdxToBox_pair.second.ly()),
    //                  boostPoint(cellIdxToBox_pair.second.hx(),
    //                             cellIdxToBox_pair.second.hy()));
    //     local_candidates_rtree.insert({box, vars_idx});
    // }//end for 

  
    // for(auto pair_idx_cell : weightIdxToCell){
    //     auto vars_idx = pair_idx_cell.first;
    //     auto cellIdxToBox_pair = pair_idx_cell.second;
    //     auto cell_idx = cellIdxToBox_pair.first;
    //     auto box = cellIdxToBox_pair.second;

    //     // std::vector<int> cells_ovrlp_constraint;
    //     // cells_ovrlp_constraint.push_back(vars_idx);
    //     int eps = 10;
    //     boostBox rtreeQueryBox(boostPoint(box.lx()+eps,
    //                                       box.ly()+eps),
    //                            boostPoint(box.hx()-eps,
    //                                       box.hy()-eps));
    //     vector<std::pair<boostBox, int>> queryResults;
    //     local_candidates_rtree.query(bgi::intersects(rtreeQueryBox), std::back_inserter(queryResults));
    //     for (const auto& queryResult : queryResults) {
    //         std::vector<int> cells_ovrlp_constraint;
    //         cells_ovrlp_constraint.push_back(vars_idx);

    //         auto conf_vars_idx = queryResult.second;
    //         auto conf_cell_idx = weightIdxToCell[conf_vars_idx].first;
    //         // log() << "cell_idx: " << cell_idx
    //         //       << ", conf_vars_idx: " << conf_vars_idx 
    //         //       << ", conf_cell_idx: " << conf_cell_idx
    //         //       << std::endl ;

    //         if( cell_idx != conf_cell_idx){
    //             cells_ovrlp_constraint.push_back(conf_vars_idx);
    //             legalize_matrix.push_back(cells_ovrlp_constraint);
    //         }
    //     }//end for 

        // log() << "cells_ovrlp_constraint: " << cells_ovrlp_constraint.size() << std::endl;

        // if(cells_ovrlp_constraint.size() > 1)
        //     legalize_matrix.push_back(cells_ovrlp_constraint);
    // }

    // 
    // log() << "legalize matrix: " << std::endl;
    // for(auto tmp_row : legalize_matrix){
    //     log() << "new tmp: " << std::endl;
    //     for(auto tmp : tmp_row) {
    //         log() << "tmp: " << tmp << std::endl;
    //     }
        
    // }

    //---------------------------------------------


    // std::vector<int> sols; 
    // std::string model_name = "model_" + std::to_string(iter) + ".lp";
    // ILPSolver(weights,weightIdxToCell,conflict_matrix,conflict_matrix_connected,legalize_matrix,sols,false,model_name);

    // // log() << "ILPSolver solution..." <<std::endl;
    // for(auto idx : sols){
    //     auto cell_box_pair = weightIdxToCell[idx];

    //     auto old_box = utils::BoxT<DBU>(
    //         database.cells[cell_box_pair.first].rsynInstance.getBounds().getLower().x,
    //         database.cells[cell_box_pair.first].rsynInstance.getBounds().getLower().y,
    //         database.cells[cell_box_pair.first].rsynInstance.getBounds().getUpper().x,
    //         database.cells[cell_box_pair.first].rsynInstance.getBounds().getUpper().y
    //     );

    //     if(!old_box.isSame(cell_box_pair.second)){
    //         if(debug){
    //             log() << "sol: " << idx << ", " 
    //                   << database.cells[cell_box_pair.first].getName() 
    //                   << ", " << cell_box_pair.second << std::endl;
    //         }
    //         cellsToMove.push_back(std::make_pair(cell_box_pair.first,cell_box_pair.second));
    //     }//end if 

    // }//end for 



}//end findBestEntryV2




void Placer::ILPSolver(std::vector<double>& weights
        ,std::unordered_map<int,std::pair<int,utils::BoxT<DBU>>>& weightIdxToCell
        ,std::vector<std::vector<int>>& conflict_matrix
        ,std::vector<std::vector<int>>& conflict_matrix_connected
        ,std::vector<std::vector<int>>& legalize_matrix
        ,std::vector<int>& sol
        ,bool maxmimize
        ,std::string log_name 
         ){
    // log() << "cplex ..." << std::endl;
    /*
        var_type should be:
            IloNumVar   : for LP
            IloBoolVar  : for ILP
    */
    // using var_type = IloBoolVar;
    // IloEnv env;
    bool debug = false;
    vector<int> selected;
    // vector<int> selectedIndexes;
    selected.resize(weights.size());
    // std::vector<int> weights = {1,2,3,5};
    IloEnv env;
    try {
        IloModel model(env);
        IloIntArray weightsSet(env);
        // IloNumArray weightsSet(env);
        IloIntVarArray vars(env);
        IloRangeArray con(env);

        for (auto i : weights) {
            // weightsSet.add(std::max(i, 1));
            weightsSet.add(i);
        }

        int counter = 0;
        // for (int i : selected) {
        //     vars.add(IloIntVar(env, 0, 1));
        //     std::string name = "x_" + std::to_string(counter);
        //     vars[i].setName(name.c_str());
        //     counter++;
        // }
        if(debug)
            log() << "weights idx ..." << std::endl;


        for (int i = 0 ; i < weights.size(); i++) {
            vars.add(IloIntVar(env, 0, 1));
            // std::string cell_name = cellsToIdx[]
            auto cell_box_pair = weightIdxToCell[i];
            auto cell_name = database.cells[cell_box_pair.first].getName();

            

            auto new_pos = cell_box_pair.second;
            double x_u = new_pos.lx()/database.libDBU;
            double y_u = new_pos.ly()/database.libDBU;
            std::stringstream stream_x;
            std::stringstream stream_y;
            std::stringstream stream;
            stream_x << std::fixed << std::setprecision(3) << x_u;
            stream_y << std::fixed << std::setprecision(3) << y_u;
            std::string s_x = stream_x.str();
            std::string s_y = stream_y.str();

            std::string name = cell_name + "_" 
                + s_x + "_" + s_y + "_"
                + std::to_string(i);
            if(debug){
              log() << "i: " << i << ", name: " << name << std::endl;
            }
            
            vars[i].setName(name.c_str());
        }

        // naming
        // only one location for each cell constraint
        int j = 0;
        int counter_cons = 0;
        for (auto conflict_row : conflict_matrix) {
            // if(j=!72){
            //     j++;
            //     continue;
            // }
            IloNumExpr expr_constraint(env);
            for(auto var_idx : conflict_row){
                expr_constraint += vars[var_idx];    
            }  
            con.add(expr_constraint == 1);  
            std::string const_name = "c_" + std::to_string(counter_cons);
            con[counter_cons].setName(const_name.c_str());
            counter_cons++;  
            expr_constraint.end();    
            // if(j==100)break;
            // j++;
            // break;
        }

        for (auto conflict_row : conflict_matrix_connected) {
            // if(j=!72){
            //     j++;
            //     continue;
            // }
            IloNumExpr expr_constraint(env);
            for(auto var_idx : conflict_row){
                expr_constraint += vars[var_idx];    
            }  
            con.add(expr_constraint <= 1);  
            std::string const_name = "cn_" + std::to_string(counter_cons);
            con[counter_cons].setName(const_name.c_str());
            counter_cons++;  
            expr_constraint.end();    
            // if(j==100)break;
            // j++;
            // break;
        }
        // model.add(con);

        // legalization constraints
        // for (auto conflict_row : legalize_matrix) {
        for (int i = 0; i < legalize_matrix.size(); i++) {
            auto conflict_row = legalize_matrix[i];
            IloNumExpr expr_constraint(env);
            for(auto var_idx : conflict_row){
                expr_constraint += vars[var_idx];    
            }  
            con.add(expr_constraint <= 1);  
            std::string const_name = "leg_" + std::to_string(counter_cons);
            con[counter_cons].setName(const_name.c_str());
            counter_cons++;  
            expr_constraint.end();    
            
        }
        
        model.add(con);
        

        IloNumExpr expression(env);
        for (int i = 0; i < weights.size(); i++) {
            expression += vars[i] * weightsSet[i];
        }
        if(!maxmimize){
            model.add(IloMinimize(env, expression));
        }else{
            model.add(IloMaximize(env, expression));
        }

        IloCplex cplex(model);
        // cplex.setOut(env.getNullStream());

        cplex.setParam(IloCplex::Param::Threads, db::setting.numThreads);
        cplex.setOut(env.getNullStream());
        if(db::setting.debug)
            cplex.exportModel(log_name.c_str());
        // log() << "before solve" << std::endl;
        if (!cplex.solve()) {
            env.error() << "iter: " << std::to_string(iter) << std::endl;
            env.error() << "Failed to optimize.\n";
            throw (-1);
        }
        
        IloNumArray vals(env);
        
        cplex.getValues(vals, vars);
        
        for (int i = 0; i < selected.size(); i++) {
            try {
                if(debug)
                    log() << "vals[" << i << "]: " << vals[i] << std::endl;
                if (vals[i] == 1) {
                    selected[i] = 1;
                    sol.push_back(i);
                }
            } catch (...) {
                cerr << "Exception on index " << i << "\n";
            }
        }
    } catch (IloException& e) {
        cerr << "Concert exception caught: " << e << endl;
    } catch (...) {
        cerr << "Unknown exception caught" << endl;
    }
    env.end();
    // return selectedIndexes;
}//end ILPAlg


void Placer::MWIS(std::set<std::string>& critical_cells_set){
    std::vector<double> weights_tmp;
    std::unordered_map<int,std::pair<int,utils::BoxT<DBU>>> weightIdxToCell_tmp;
    std::vector<std::vector<int>> conflict_matrix_tmp;
    std::vector<std::vector<int>> conflict_matrix_connected_tmp;
    std::vector<std::vector<int>> legalize_matrix_tmp;
    std::vector<int> sols_tmp; 

    // int i = 0;
    // log() << "construct ILP model ..." << std::endl;
    // construct ILP model and find add constraints of having one solution
    // for each cell. 
    // weightIdxToCell_tmp.resize(cells.size());
    int i = 0;
    for(auto& cell : cells){
        std::vector<int> connected_cells;
        // log() << "cell: " << cell.getName() << ", idx: " << cell.idx << std::endl;
        auto box = utils::BoxT<DBU>(cell.rsynInstance.getBounds().getLower().x,
                                    cell.rsynInstance.getBounds().getLower().y
                                   ,cell.rsynInstance.getBounds().getUpper().x,
                                    cell.rsynInstance.getBounds().getUpper().y);
        weightIdxToCell_tmp[cell.idx] = std::make_pair(cell.idx,box);
        auto cost = cell.getWirelength();
        if(cost == 0){
            weights_tmp.push_back(1);
        }else{
            weights_tmp.push_back(cost);
        }
        
        connected_cells.push_back(cell.idx);

        for(auto connected_cell_idx : cell.connectd_cells){
            connected_cells.push_back(connected_cell_idx);
                // log() << "connected_cell_idx: " << connected_cell_idx << std::endl;
            // log() << "ILPidx: " << i <<", cell: " << cell.getName() << ", cell.candidate: " << box
            //       << ", cost: " << cost << std::endl;
        }
        // i++;
        // if(i==5) break;
        conflict_matrix_tmp.push_back(connected_cells);
        
    }//end for 

    // for(auto tmps: conflict_matrix_tmp){
    //     log() << "new_conflict" << std::endl;
    //     for(auto tmp : tmps){
    //         log() << "tmp" << tmp << std::endl;
    //     }
    // }
    for(auto tmp : weights_tmp){
        log() << "w: " << tmp << std::endl;
    }

    ILPSolver(weights_tmp,weightIdxToCell_tmp,conflict_matrix_tmp,conflict_matrix_connected_tmp,legalize_matrix_tmp,sols_tmp,true,"modelMWIS.lp");


    for(auto sol : sols_tmp){
        // log() << "critical cell: " << cells[sol].getName() << std::endl;
        critical_cells_set.insert(cells[sol].getName());
    }
}//end MWIS


void Placer::schellingSegerationModel(std::set<std::string>& critical_cells_set){
    bool debug = false;
    std::vector<int> nets_tmp;
    std::vector<int> cells_criticality_order;
    // The percentage of critical cells to select and move
    // 0.01 of 8734 is 87 crticial cell.
    float critical_cells_percentage =  db::setting.percCriticalCells;


    // for (auto& net : grDatabase.nets) nets_tmp.push_back(net.dbNet.idx);
    for (auto& cell : database.cells) cells_criticality_order.push_back(cell.idx);
    

    //  runJobsMT(cells.size(), [&](int i) { cells[i].calcCostPlacementCandidates(); });
    // sort(nets_tmp.begin(), nets_tmp.end(), [&](int id1, int id2) {
    //     return grDatabase.nets[id1].getWirelength() > grDatabase.nets[id2].getWirelength();
    // });

    if(db::setting.debug){
        log() << "sort cells_criticality_order" << std::endl;
    }
    // sort(cells_criticality_order.begin(), cells_criticality_order.end(), [&](int id1, int id2) {
    //     return database.cells[id1].getWirelengthLast() > database.cells[id2].getWirelengthLast();
    // });
    sort(cells_criticality_order.begin(), cells_criticality_order.end(), [&](int id1, int id2) {
        return database.cells[id1].feature.getCost() > database.cells[id2].feature.getCost();
    });
    if(db::setting.debug){
        log() << "end sort cells_criticality_order" << std::endl;
    }
  
    std::random_device dev;
    std::mt19937 rng(dev());
    // std::uniform_int_distribution<std::mt19937::result_type> rnd_fun(0,1);
    std::uniform_real_distribution<> rnd_fun(0,1);
    
    if(db::setting.debug){
        log() << "add cells_criticality_order" << std::endl;
    }
    for(int i = 0; i < cells_criticality_order.size(); i ++){
        std::set<int> connected_cells_set;
        std::vector<int> nets_idx;
        auto cell = database.cells[cells_criticality_order[i]];
        cell.getCellNetsIdx(nets_idx);
        if(cell.isMultiRowCell()) continue;
        if(nets_idx.size() == 0) continue;
        

        for(auto connected_cell_idx : cell.connectd_cells){
            connected_cells_set.insert(connected_cell_idx);
        }//end for 
        bool unique_cell_in_net = false;
        for(auto tmp : connected_cells_set){
            if(critical_cells_set.find(database.cells[tmp].getName()) != critical_cells_set.end()){
                unique_cell_in_net = true;
                break;
            }
        }
        if(unique_cell_in_net) continue;

        auto cell_idx = cell.idx;
        // if we could not find the critical cell in current critical cell set. 
        if(critical_cells_set.find(cell.getName()) == critical_cells_set.end()){
            int was_in_critical_cells = 0;
            int was_in_critical_cells_moved = 0;
            if(database.critical_cells_hist.find(cell.idx) != database.critical_cells_hist.end()){
                was_in_critical_cells = 1;
                if(database.critical_cells_hist_moved.find(cell.idx) != database.critical_cells_hist_moved.end()){
                    was_in_critical_cells_moved = 1;
                }
            }
            
            // DBU prev_cost = database.cost_hist[database.cost_hist.size()-2];
            // DBU cur_cost = database.cost_hist[database.cost_hist.size()-1];
            float temperature = 1;
            double simulated_annealing_exp = std::exp((-1*(was_in_critical_cells+was_in_critical_cells_moved))/temperature);
            double rnd_val = rnd_fun(rng);
            // log() << "prev_cost: " << prev_cost << ", cur_cost: " << cur_cost << std::endl;
            // log() << "simulated_annealing_exp: " << simulated_annealing_exp
            //       << ", rnd_val: " << rnd_val << std::endl;
            bool debug= true;
            if(simulated_annealing_exp > rnd_val){
                if(debug){
                    if(cell.getName() == "inst31497")
                        continue;
                    else
                        critical_cells_set.insert(cell.getName());
                }else{
                    critical_cells_set.insert(cell.getName());
                }
                
            }

            // if(exp^)
        }

        if(critical_cells_set.size() > critical_cells_percentage*database.cells.size()){
            break;
        }

    }//end for 
    if(db::setting.debug){
        log() << "end cells_criticality_order" << std::endl;
    }


    for(auto cell_name : critical_cells_set){
        database.critical_cells_hist.insert(database.cellsToIdx[cell_name]);
    }

    // sort(nets_tmp.begin(), nets_tmp.end(), [&](int id1, int id2) {
    //     return grDatabase.nets[id1].getWirelength() > grDatabase.nets[id2].getWirelength();
    // });

    // // log() << "list of angry nets..." << std::endl;
    // // for(auto idx : nets_tmp){
    // //     if(debug){
    // //         log() << "net: "  << grDatabase.nets[idx].getName() 
    // //               << ", wl: " << grDatabase.nets[idx].getWirelength()  << std::endl;
    // //     }
    // // }

    // // log() << "list of angry cells..." << std::endl;
    // // for(auto idx : cells_tmp){
    // //     if(debug){
    // //         log() << "cell: " << database.cells[idx].getName() 
    // //               << ", wl: " << database.cells[idx].getWirelength() << std::endl;
    // //     }
    // //     // if(database.cells[idx].getName() == "inst4132"){
    // //     //     critical_cells_set.insert("inst4132");
    // //     // }   
    // // }

    // // log() << "list of angry wires" << std::endl;
    // // std::vector<int,DBU> nets_longestWire;
    // // for(auto net_idx : nets_tmp){
    // //     auto long_wire = grDatabase.nets[net_idx].getLongestWire();
    // //     nets_longestWire.push_back(std::make_pair(net_idx,long_wire));
    // // }

    // sort(nets_tmp.begin(), nets_tmp.end(), [&](int id1, int id2) {
    //         return grDatabase.nets[id1].getLongestWire() > grDatabase.nets[id2].getLongestWire();
    //     });

    // // for(auto net_idx : nets_tmp){
    // //     log() << "grNet: " << grDatabase.nets[net_idx].getName()
    // //           << ", wireLen: " << grDatabase.nets[net_idx].getLongestWire() << std::endl;
    // // }

    // // log() << "list of highest edge cost" << std::endl;

    // sort(nets_tmp.begin(), nets_tmp.end(), [&](int id1, int id2) {
    //         return grDatabase.nets[id1].getEdgeHighestCost() > grDatabase.nets[id2].getEdgeHighestCost();
    //     });

    // for(auto net_idx : nets_tmp){
    //     log() << "grNet: " << grDatabase.nets[net_idx].getName()
    //           << ", wireLen: " << grDatabase.nets[net_idx].getEdgeHighestCost() << std::endl;
    // }

}//end schellingSegerationModel


void Placer::logRefinePlacement(){
    auto rsynService = database.getRsynService();
    const Rsyn::Session session;
    const Rsyn::PhysicalDesign& physicalDesign =
    static_cast<Rsyn::PhysicalService*>(session.getService("rsyn.physical"))->getPhysicalDesign();
    const DBU libDBU = physicalDesign.getDatabaseUnits(Rsyn::LIBRARY_DBU);
    auto dieBound = rsynService.physicalDesign.getPhysicalDie().getBounds();
    DBU die_lx = dieBound.getLower().x;
    DBU die_ly = dieBound.getLower().y;
    DBU die_hx = dieBound.getUpper().x;
    DBU die_hy = dieBound.getUpper().y;
    std::unordered_map<std::string, std::vector<std::pair<double,double>> > lib_abundance_dict;


    for(auto cell : database.cells){
        DBU x_center = (cell.rsynInstance.getBounds().getLower().x + cell.rsynInstance.getBounds().getUpper().x)/2;
        DBU y_center = (cell.rsynInstance.getBounds().getLower().y + cell.rsynInstance.getBounds().getUpper().y)/2;
        double x_norm = double(x_center)/double(die_hx-die_lx);
        double y_norm = double(y_center)/double(die_hy-die_ly);
        Rsyn::Cell cellRsyn = cell.rsynInstance.asCell();
        Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cellRsyn);
        Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cellRsyn);
        log() << "cell: " << cell.getName() 
              << ", lib: " << cellRsyn.getLibraryCellName() 
              << ", x_norm: " << x_norm << ", y_norm: " << y_norm << std::endl;

        if(lib_abundance_dict.find(cellRsyn.getLibraryCellName()) == lib_abundance_dict.end()){
            std::vector<std::pair<double,double>> tmp;
            tmp.push_back(std::make_pair(x_norm,y_norm));
            lib_abundance_dict[cellRsyn.getLibraryCellName()] = tmp;
        }else{
            auto& tmp = lib_abundance_dict[cellRsyn.getLibraryCellName()];
            tmp.push_back(std::make_pair(x_norm,y_norm));
        }
    }

    log() << "abundance: " << std::endl;

    for(auto lib_abundance_pair : lib_abundance_dict){
        auto lib_name = lib_abundance_pair.first;
        auto vec = lib_abundance_pair.second;
        log() << "lib: " << lib_name << std::endl;
        for(auto tmp : vec){
            log() << tmp.first << ", " << tmp.second << std::endl;
        }
    }

    log() << "number of logics: " << lib_abundance_dict.size() << std::endl;
    for(auto lib_abundance_pair : lib_abundance_dict){
        auto lib_name = lib_abundance_pair.first;
        auto vec = lib_abundance_pair.second;
        log() << "logic: " << lib_name << ", num: " << vec.size() << std::endl;
    }

    // std::vector<std::vector<int>> matrix_abnd;
    std::vector <std::vector<int> > vec2D(10, std::vector<int>(10, 0));
    // MX2XL
    for(auto lib_abundance_pair : lib_abundance_dict){
        auto lib_name = lib_abundance_pair.first;
        auto vec = lib_abundance_pair.second;
        if(lib_name == "MX2XL"){
            for(auto tmp : vec){
                int i = int(10*tmp.first);
                int j = int(10*tmp.second);
                vec2D[i][j] +=1;
            }   
        }
        // log() << "lib: " << lib_name << std::endl;
        
    }
    log() << "vec2D..." << std::endl;
    for(int i = 0; i < vec2D.size() ; i++){
        log() << vec2D[i][0] << ", " << vec2D[i][1] << ", " 
              << vec2D[i][2] << ", " << vec2D[i][3] << ", " << vec2D[i][4] << ", "
              << vec2D[i][5] << ", " << vec2D[i][6] << ", " 
              << vec2D[i][7] << ", " << vec2D[i][8] << ", " << vec2D[i][9] << std::endl;
    }

    

    // for(auto )

    // std::stringstream ss;
    // // for (Rsyn::Instance instance : rsynService.module.allInstances()) {
    // //     // if (instance.getName() != "inst5638") continue;
    // //     if (instance.getType() != Rsyn::CELL) continue;
    // //     // phCell
    // //     Rsyn::Cell cell = instance.asCell();
    // //     Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cell);
    // //     Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cell);
    // //     ss << instance.getName() << ", " << instance.getX() << ", "
    // //         << instance.getY() << ", " << phLibCell.getWidth() << ", " 
    // //         << phLibCell.getHeight() << std::endl;
    // // }
    // std::string file_name_csv = "instances.csv";
    // std::ofstream fout(file_name_csv);
    // fout << ss.str();
    // fout.close();

    std::vector<DBU> diff_x_vec;
    std::vector<DBU> diff_y_vec;

    for(auto cell_box_pair: cellsToMove){
        auto cell = database.cells[cell_box_pair.first];
        auto box = cell_box_pair.second;
        Rsyn::Cell cellRsyn = cell.rsynInstance.asCell();
        Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cellRsyn);
        Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cellRsyn);
        auto lx_old = cell.rsynInstance.getBounds().getLower().x;
        auto ly_old = cell.rsynInstance.getBounds().getLower().y;
        auto hx_old = cell.rsynInstance.getBounds().getUpper().x;
        auto hy_old = cell.rsynInstance.getBounds().getUpper().y;
        
        auto lx_new = box.lx();
        auto ly_new = box.ly();
        auto hx_new = box.hx();
        auto hy_new = box.hy();
        diff_x_vec.push_back(std::abs(lx_new-lx_old));
        diff_y_vec.push_back(std::abs(ly_new-ly_old));
        
        log() << cell.getName() 
              << ", " << cellRsyn.getLibraryCellName() 
              << ", " << std::to_string(iter) 
              << ",x_diff: " << std::abs(lx_new-lx_old)
              << ",y_diff: " << std::abs(ly_new-ly_old) 
              << std::endl;
    }

    
    DBU max_diff_x_vec = *std::max_element(diff_x_vec.begin(),diff_x_vec.end());
    DBU max_diff_y_vec = *std::max_element(diff_y_vec.begin(),diff_y_vec.end());

    log() << "max_x: " << max_diff_x_vec << ", max_y: " << max_diff_y_vec << std::endl;

}


void Placer::calcNetsHPWL(){
    DBU hpwl_total = 0;
    bool debug = false;
    std::string file_name = "net_hpwl_"+std::to_string(iter) + ".csv";
    
    std::ofstream file(file_name);
    std::stringstream stream;
    stream << "net,hpwl" << std::endl;

 
    for(auto& grNet : grDatabase.nets){
        auto hpwl_tmp = grNet.getHPWL();
        hpwl_total += hpwl_tmp;
        stream << grNet.getName() << "," << std::to_string(hpwl_tmp)<<std::endl;
    }
    log() << "total HPWL: " << hpwl_total << std::endl;
    database.db_hpwls.push_back(hpwl_total);
    file << stream.str();
    file.close();
}

void Placer::logCellLocations(){
    bool debug = false;
    std::string file_name = db::setting.benchmarkName+ ".cell_"+std::to_string(iter) + ".csv";
    if(debug)log() << "logCellLocations file name: " << file_name << std::endl;
    std::ofstream file(file_name);
    std::stringstream stream;
    stream << "cell_name,x,y,w,h" << std::endl;

    for(auto cell : cells){
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


void Placer::logCellsFeature(){
    bool debug = false;
    std::string file_name = db::setting.directory +  
        db::setting.benchmarkName+ ".cells.features."+std::to_string(iter) + ".csv";
    
    if(debug)log() << "logCellLocations file name: " << file_name << std::endl;
    std::ofstream file(file_name);
    std::stringstream stream;
    // medCx -> median cell x
    // medPx -> median pin x
    // cong -> congestion
    stream << "cell_name,xl,yl,xh,yh,medCx,medCy,medPx,medPy,wl,cong,deg_nets,deg_cells,suspected,maxEdgeCost"<< std::endl;

    for(auto cell : cells){
        auto box = cell.getCellBox();
        stream << cell.getName() 
            << "," << box.lx()
            << "," << box.ly()
            << "," << box.hx()
            << "," << box.hy()
            << "," << cell.feature.median_x_y_cell.first
            << "," << cell.feature.median_x_y_cell.second
            << "," << cell.feature.median_x_y_pin.first
            << "," << cell.feature.median_x_y_pin.second
            << "," << cell.feature.wl
            << "," << cell.feature.congestion
            << "," << cell.feature.degree_connected_nets
            << "," << cell.feature.degree_connected_cells
            << "," << cell.feature.is_suspected_cell
            << "," << cell.feature.max_edge_cost << std::endl;
    }
    file << stream.str();
    file.close();

}//end logCellLocations


void Placer::logCellsCritical(){
    bool debug = false;
    std::string file_name = db::setting.directory +  
        db::setting.benchmarkName+ ".cells.critical."+std::to_string(iter) + ".csv";
    
    if(debug)log() << "logCellLocations file name: " << file_name << std::endl;
    std::ofstream file(file_name);
    std::stringstream stream;
    // medCx -> median cell x
    // medPx -> median pin x
    // cong -> congestion
    stream << "cell_name,xl,yl,xh,yh,medCx,medCy,medPx,medPy,wl,cong,deg_nets,deg_cells,suspected,maxEdgeCost"<< std::endl;

    for(auto cell_idx : critical_cells){
        auto cell = database.cells[cell_idx];
        auto box = cell.getCellBox();
        stream << cell.getName() 
            << "," << box.lx()
            << "," << box.ly()
            << "," << box.hx()
            << "," << box.hy()
            << "," << cell.feature.median_x_y_cell.first
            << "," << cell.feature.median_x_y_cell.second
            << "," << cell.feature.median_x_y_pin.first
            << "," << cell.feature.median_x_y_pin.second
            << "," << cell.feature.wl
            << "," << cell.feature.congestion
            << "," << cell.feature.degree_connected_nets
            << "," << cell.feature.degree_connected_cells
            << "," << cell.feature.is_suspected_cell
            << "," << cell.feature.max_edge_cost << std::endl;
    }
    file << stream.str();
    file.close();
}//end logCellsCritical

void Placer::logCellsCandidates(){
    bool debug = false;
    std::string file_name = db::setting.directory +  
        db::setting.benchmarkName+ ".cells.candidates."+std::to_string(iter) + ".csv";
    
    if(debug)log() << "logCellLocations file name: " << file_name << std::endl;
    std::ofstream file(file_name);
    std::stringstream stream;
    // medCx -> median cell x
    // medPx -> median pin x
    // cong -> congestion
    stream << "cell_name,xl,yl,xh,yh,type,cost"<< std::endl;

    for(auto cell : cells){
        for(auto pl_cad : cell.placement_candidates){
            auto box = pl_cad.first;
            stream << cell.getName() 
                << "," << box.lx()
                << "," << box.ly()
                << "," << box.hx()
                << "," << box.hy()
                << "," << pl_cad.second << std::endl;
        }//end for 
    }//end for    

    file << stream.str();
    file.close();
}//end logCellsCandidates

void Placer::logCellsMoved(){
    bool debug = false;
    std::string file_name = db::setting.directory +  
        db::setting.benchmarkName+ ".cells.moved."+std::to_string(iter) + ".csv";
    
    if(debug)log() << "logCellLocations file name: " << file_name << std::endl;
    std::ofstream file(file_name);
    std::stringstream stream;
    // medCx -> median cell x
    // medPx -> median pin x
    // cong -> congestion
    stream << "cell_name,xl,yl,xh,yh,newXl,newYl,newXh,newYh"<< std::endl;

    for(auto cell_box_pair : cellsToMove){
        auto cell = database.cells[cell_box_pair.first];
        auto oldBox = cell.getCellBox();
        auto newBox = cell_box_pair.second;
       
        stream << cell.getName() 
            << "," << oldBox.lx()
            << "," << oldBox.ly()
            << "," << oldBox.hx()
            << "," << oldBox.hy()
            << "," << newBox.lx()
            << "," << newBox.ly()
            << "," << newBox.hx()
            << "," << newBox.hy()
            << std::endl;
       
    }//end for    

    file << stream.str();
    file.close();

}//end logCellsMoved


void Placer::logPlacementWeights(std::vector<CandWeight>& weights){
    std::stringstream ss;
    auto rsynService = database.getRsynService();
       
    ss << "weight_idx,cell_name,xl,yl,xh,yh,cost" << std::endl;    

    for(int i = 0; i < weights.size(); i++){
        auto cell_idx = weights[i].cell_idx;
        auto box = weights[i].box;
        auto cost = weights[i].cost;
    

        ss << i
           << "," << database.cells[cell_idx].getName()
           << "," << box.lx()
           << "," << box.ly()
           << "," << box.hx()
           << "," << box.hy()
           << "," << cost << std::endl;


    }//end for loop


    
    std::string file_name_csv_b = db::setting.directory +  db::setting.benchmarkName 
                                + ".placement.weights.csv";
    std::ofstream fout_blockage(file_name_csv_b);
    fout_blockage << ss.str();
    fout_blockage.close();
}//end logLegalizerBoard



