#include "Router.h"
#include "flute/flute.h"
#include "Scheduler.h"
#include "single_net/InitRoute.h"
#include "place/Placer.h"

extern "C" {
void readLUT();
}

const MTStat& MTStat::operator+=(const MTStat& rhs) {
    auto dur = rhs.durations;
    std::sort(dur.begin(), dur.end());
    if (durations.size() < dur.size()) {
        durations.resize(dur.size(), 0.0);
    }
    for (int i = 0; i < dur.size(); ++i) {
        durations[i] += dur[i];
    }
    return *this;
}

ostream& operator<<(ostream& os, const MTStat mtStat) {
    double minDur = std::numeric_limits<double>::max(), maxDur = 0.0, avgDur = 0.0;
    for (double dur : mtStat.durations) {
        minDur = min(minDur, dur);
        maxDur = max(maxDur, dur);
        avgDur += dur;
    }
    avgDur /= mtStat.durations.size();
    os << "#threads=" << mtStat.durations.size() << " (dur: min=" << minDur << ", max=" << maxDur << ", avg=" << avgDur
       << ")";
    return os;
}
Router::Router() {
    routeTable.resize(database.nets.size() + 1);
    readLUT();  // read flute LUT
   
    boost::split(filter_routers, db::setting.rrrRouters, boost::is_any_of(","));
    boost::split(filter_routers_applys, db::setting.rrrRoutersApply, boost::is_any_of(","));
    

    log() << "router in each iteration: "  << std::endl;
    for(int i = 0; i < filter_routers.size() ;i++){
        if(filter_routers_applys[i] == "1")
            log() << "router: " << filter_routers[i] << std::endl;
    }

    if(filter_routers.size() != db::setting.rrrIterLimit){
        log() << "Error: number of iteration of routers (rrrIters) is not equal to rrrRoutes!" << std::endl;
        // std::exit(1);
    }
    
}

// void Router::runISPD() {
//     bool debug = false;
//     log() << "run ISPD engine ..." << std::endl;


//     coef_stream << "iter,unitWireCostRaw,unitViaCostRaw,unitShortVioCostRaw,unitShortVioCost,\
//                     unitShortVioCostDiscounted,unitViaCost,step,db::setting.rrrInitVioCostDiscount,\
//                     db::setting.rrrIterLimit,db::setting.initLogisticSlope,db::setting.rrrFadeCoeff,\
//                     unitViaMultiplier,logisticSlope,wireCapDiscount" << std::endl;
//     std::string coef_stream_file_name = db::setting.outputFile + ".gr.coef.csv";
//     std::ofstream coef_stream_file(coef_stream_file_name);

//     log() << "coef: " << coef_stream_file_name << std::endl;

//     logCoef();

//     filterNets();


//     allNetStatus.resize(database.nets.size(), db::RouteStatus::FAIL_UNPROCESSED);
//     vector<int> netsToRoute;
//     int num_iter=0;

//     if(debug){
//         log() << "before looop ... " << std::endl;
//         grDatabase.statHistCost();
//     }
        
//     for(int i = 0; i <= num_iter; i++){  
//         // //-----
//         // iter=0;
//         // netsToRoute = getNetsToRoute();
//         // ripup(netsToRoute);
//         // Placer place(congMap,routeTable,i);
//         // place.updateRoutingDataBase(netsToRoute);
//         // updateCostInit();
//         // netsToRoute.clear();
//         // if(debug){
//         //     log() << "before applyRoute#" << i << std::endl;
//         //     grDatabase.statHistCost();
//         // }
//         // // Route
//         log() << "start update" << std::endl;
//         database.init();
//         grDatabase.init();
//         log() << "stop update" << std::endl;

//         applyOnlyRoute(netsToRoute);  

//         if(num_iter-1 == i) break;
//         utils::timer profile_time;
//         std::stringstream profile_time_str;
//         Placer place(congMap,routeTable,i,profile_time,profile_time_str);
//         place.runMTISPD(netsToRoute,cellWidth, cellHeight);
//         printStat();      
//         if(debug){
//             log() << "after applyRoute#" << i << std::endl;
//             grDatabase.statHistCost();
//         }
//         // if(num_iter==i)
//         //     break;
//         // applyPlacement(netsToRoute,i);
//     }
    
    

    


//     // postprocessing
//     if(db::setting.postProcessing){
//         log() << "postProcessing..." << std::endl;
//         for (auto& net : grDatabase.nets) {
//             if(filter_nets.size() != 0){
//                 if(filter_nets.find(net.dbNet.idx) != filter_nets.end()){
//                     GuideGenerator guideGen(net);
//                     guideGen.genPatchGuides();
//                 }
//             }else{
//                 GuideGenerator guideGen(net);
//                 guideGen.genPatchGuides();
//             }
            
//         }
//         log() << "end postProcessing..." << std::endl;
//     }

//     if(db::setting.debug){    
//         if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) guideGenStat.print();

//         log() << std::endl;
//         log() << "################################################################" << std::endl;
//         database.setUnitVioCost(1);  // set the cost back to without discount
//         log() << "Finish all RRR iterations and PostRoute" << std::endl;
//         log() << "MEM: cur=" << utils::mem_use::get_current() << "MB, peak=" << utils::mem_use::get_peak() << "MB"
//             << std::endl;

//         printStat();
//     }

//     // close log files
//     coef_stream_file << coef_stream.str();
//     coef_stream_file.close();
    

// }//end runISPD

// void Router::run() {
//     allNetStatus.resize(database.nets.size(), db::RouteStatus::FAIL_UNPROCESSED);
//     for (iter = 0; iter < db::setting.rrrIterLimit; iter++) {
//         log() << std::endl;
//         log() << "################################################################" << std::endl;
//         log() << "Start RRR iteration " << iter << std::endl;
//         log() << std::endl;
//         db::routeStat.clear();
//         guideGenStat.reset();

//         vector<int> netsToRoute;
//         getNetsToRoute(netsToRoute,iter);
//         if (netsToRoute.empty()) {
//             if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) {
//                 log() << "No net is identified for this iteration of RRR." << std::endl;
//                 log() << std::endl;
//             }
//             break;
//         }
//         sortNets(netsToRoute);  // Note: only effective when doing mazeroute sequentially

//         updateCost(iter);
//         grDatabase.statHistCost();

//         if (iter > 0) {
//             ripup(netsToRoute);
//             congMap.init(cellWidth, cellHeight);
//         }

//         if(iter == 0)
//             routeApprx(netsToRoute,PATTERNROUTE);
//         else 
//             routeApprx(netsToRoute,ASTAR);

//         log() << std::endl;
//         log() << "Finish RRR iteration " << iter << std::endl;
//         log() << "MEM: cur=" << utils::mem_use::get_current() << "MB, peak=" << utils::mem_use::get_peak() << "MB"
//               << std::endl;
//         if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) db::routeStat.print();
//     }
//     // postprocessing
//     for (auto& net : grDatabase.nets) {
//         GuideGenerator guideGen(net);
//         guideGen.genPatchGuides();
//     }
//     if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) guideGenStat.print();

//     log() << std::endl;
//     log() << "################################################################" << std::endl;
//     database.setUnitVioCost(1);  // set the cost back to without discount
//     log() << "Finish all RRR iterations and PostRoute" << std::endl;
//     log() << "MEM: cur=" << utils::mem_use::get_current() << "MB, peak=" << utils::mem_use::get_peak() << "MB"
//           << std::endl;

//     printStat();
// }


void Router::run() {
    allNetStatus.resize(database.nets.size(), db::RouteStatus::FAIL_UNPROCESSED);
    utils::timer profile_time;
    std::stringstream profile_time_str;
    // logAll();
    // return;
    
    
    
    // it is used to route specific nets
    filterNets();

    vector<int> netsToRoute;
    getNetsToRoute(netsToRoute,0);
    netFeatureExtraction2D(netsToRoute);
    

    // return;

    if(!db::setting.refinePlacement){
        
    // ripupRerouteCUGR(netsToRoute); 
        ripupRerouteCRP(netsToRoute); 
    
    }else{
        vector<int> netsToRoute;
        vector<int> netsToRouteByPatternRoute;
        vector<int> netsToRouteByAStar;
        // ripupReroute(netsToRoute);
        ripupRerouteCRP(netsToRoute); 
        updateCostPlacement(0);
        netsToRoute.clear();
        for(int i = 0; i < db::setting.numRefinePlacement; i++){
            printRouteStart("placement",i);
            applyPlacement(netsToRoute,i, profile_time, profile_time_str);
            
            // route(netsToRoute,PATTERNROUTE);
            // route(netsToRoute,ASTAR);

            // need to select right router 
            classifyNets(netsToRoute, netsToRouteByPatternRoute, netsToRouteByAStar );
        
            // if(filter_routers[iter] == "patternroute" && filter_routers_applys[iter] =="1"){
            // route by patternRoute
            ripup(netsToRoute);
            congMap.init(cellWidth, cellHeight);
            routeApprx(netsToRoute, PATTERNROUTE);
        
            updateRouteTable(); 
            printStat();
            logAll();        

        // }
            // if(filter_routers[iter] == "astar" && filter_routers_applys[iter] =="1" ) {
            // route by Astar
            // ripup(netsToRouteByAStar);
            
            congMap.init(cellWidth, cellHeight);
            routeApprx(netsToRouteByAStar, ASTAR);
            updateRouteTable(); 
            printStat();
            logAll();

            netsToRoute.clear();
            logAll();
        }//end refinePlacement loop

    }//end refinePlacement


    // postprocessing
    if(db::setting.postProcessing){
        if(db::setting.debug){
            log() << "postProcessing..." << std::endl;
        }
        for (auto& net : grDatabase.nets) {
            if(filter_nets.size() != 0){
                if(filter_nets.find(net.dbNet.idx) != filter_nets.end()){
                    GuideGenerator guideGen(net);
                    guideGen.genPatchGuides();
                }
            }else{
                GuideGenerator guideGen(net);
                guideGen.genPatchGuides();
            }
            
        }
        if(db::setting.debug){
            log() << "end postProcessing..." << std::endl;
        }
    }

    logAll();

    
 

    // if(!db::setting.refinePlacement){
        
    //     applyOnlyRoute(netsToRoute);        
        

    //     // reroute nets that still have violation
    //     for (auto& net : grDatabase.nets)
    //         if (grDatabase.hasVio(net)) netsToRoute.push_back(net.dbNet.idx);


    //     for(int idx :netsToRoute ){
    //         std::vector<int> netsToRouteTmp2;
    //         netsToRouteTmp2.push_back(idx);
    //         routeAStarSeq(netsToRouteTmp2);
    //     }
        

    //     db::routeStat.clear();
    //     guideGenStat.reset();
    //     congMap.init(cellWidth, cellHeight);   

    //     netsToRoute.clear();
    //     for (auto& net : grDatabase.nets){
    //         if (grDatabase.hasVio(net)) {
    //             netsToRoute.push_back(net.dbNet.idx);
    //             // log() << "violation " << net.getName() << std::endl;
    //         }
    //     }
            
        
            
       



            
    // }else{
    //     for(int j = 0; j < db::setting.numGlobalRouting; j++){
    //         log() << "global routing + refine Placement..." << std::endl;
    //         // 1- apply cugr flute route & maze route
    //         applyOnlyRoute(netsToRoute); 
    //         printStat();
    //         debug=true;
    //         if(debug){
    //             int layerIdx_tmp = 3;
    //             int gridline_tmp = 366;
    //             int cp_tmp = 277;

    //             log() << "l: " << layerIdx_tmp
    //                 << ", gridline: " << gridline_tmp
    //                 << ", cp: " << cp_tmp
    //                 << ", wire usage: " << grDatabase.getWireUsage(layerIdx_tmp, gridline_tmp, cp_tmp)
    //                 << ", fixed usage: " << grDatabase.getFixedUsage(layerIdx_tmp, gridline_tmp, cp_tmp)
    //                 << ", Tracks: " << grDatabase.getNumTracks(layerIdx_tmp, gridline_tmp) << std::endl;
    //         }   
    //         debug=false;
    //         // profile_time_str << "applyOnlyRoute_" << std::to_string(j) << "," << std::to_string(profile_time.getTimer()) << std::endl;
    //         profile_time_str << "applyOnlyRoute"  << "," << std::to_string(profile_time.getTimer()) << std::endl;

    //         if(db::setting.debug){
    //             logReport(report_tracker,j,-1);
    //         }
    //         // profile_time_str << "logReport_" << std::to_string(j) << "," << std::to_string(profile_time.getTimer()) << std::endl;

    //         DBU total_wl = 0;
    //         for(auto net : grDatabase.nets){
    //             total_wl += net.getWirelength();
    //         }
    //         database.cost_hist[database.cost_hist.size()-1] = total_wl;
    //         // profile_time_str << "WLHistUpdate_" << std::to_string(j)  << ","<< std::to_string(profile_time.getTimer()) << std::endl;
            
    //         if(db::setting.debug){
    //             log() << "#########################" << std::endl;
    //             log() << "State After Only Route..." << std::endl;
    //             printStat();
    //             log() << "#########################" << std::endl;
    //         }
    //         // profile_time_str << "printStat_" << std::to_string(j) << "," << std::to_string(profile_time.getTimer()) << std::endl;

    //         // log() << "before placement and route wirelength..." << std::endl;
    //         // for(auto net : grDatabase.nets){
    //         //     log() << "net: " << net.getName() << ", wl: " << net.getWirelength() << std::endl;
    //         // }

    //         // logRouteTable("routeTable.csv");
    //         // printCongMap("after_routing.csv");
            
            

    //         // placement and route
    //         for(int i = 0; i < db::setting.numRefinePlacement; i++){
                
    //             log() << "################################################################" << std::endl;
    //             log() << "Start Placer " << i << std::endl;
                
                
    //             applyPlacement(netsToRoute,i, profile_time, profile_time_str);
    //             debug=true;
    //             if(debug){
    //                 int layerIdx_tmp = 3;
    //                 int gridline_tmp = 366;
    //                 int cp_tmp = 277;

    //                 log() << "l: " << layerIdx_tmp
    //                     << ", gridline: " << gridline_tmp
    //                     << ", cp: " << cp_tmp
    //                     << ", wire usage: " << grDatabase.getWireUsage(layerIdx_tmp, gridline_tmp, cp_tmp)
    //                     << ", fixed usage: " << grDatabase.getFixedUsage(layerIdx_tmp, gridline_tmp, cp_tmp)
    //                     << ", Tracks: " << grDatabase.getNumTracks(layerIdx_tmp, gridline_tmp) << std::endl;
    //             }   
    //             debug=false;
    //             // profile_time_str << "applyPlacement_" << std::to_string(j) << "_"<< \
    //             //     std::to_string(i)<< "," << std::to_string(profile_time.getTimer()) << std::endl;
    //             profile_time_str << "applyPlacement" << "," << std::to_string(profile_time.getTimer()) << std::endl;
    //             if(db::setting.debug){
    //                 log() << "after placement and route wirelength..." << std::endl;
    //             }

    //             // calc nets characterestics that 
    //             // are targeted to be rerouted


    //             if(db::setting.debug){
    //                 logNetsToRoute(netsToRoute);
    //             }
    //             // profile_time_str << "logNetsToRoute_" << std::to_string(j) << "_"<< \
    //             //     std::to_string(i)<< "," << std::to_string(profile_time.getTimer()) << std::endl;
    //             profile_time_str << "logNetsToRoute" << "," << std::to_string(profile_time.getTimer()) << std::endl;
    //             // for(auto net : grDatabase.nets){
    //             //     log() << "net: " << net.getName() << ", wl: " << net.getWirelength() << std::endl;
    //             // }
    //             // printCongMap("after_placement.csv");
    //             applyFluteRoute3D(netsToRoute);
    //             debug=true;
    //             if(debug){
    //                 int layerIdx_tmp = 3;
    //                 int gridline_tmp = 366;
    //                 int cp_tmp = 277;

    //                 log() << "l: " << layerIdx_tmp
    //                     << ", gridline: " << gridline_tmp
    //                     << ", cp: " << cp_tmp
    //                     << ", wire usage: " << grDatabase.getWireUsage(layerIdx_tmp, gridline_tmp, cp_tmp)
    //                     << ", fixed usage: " << grDatabase.getFixedUsage(layerIdx_tmp, gridline_tmp, cp_tmp)
    //                     << ", Tracks: " << grDatabase.getNumTracks(layerIdx_tmp, gridline_tmp) << std::endl;
    //             }   
    //             debug=false;
    //             // profile_time_str << "applyFluteRoute3D_" << std::to_string(j) << "_"<< \
    //             //     std::to_string(i)<< "," << std::to_string(profile_time.getTimer()) << std::endl;
    //             profile_time_str << "applyFluteRoute3D" << "," << std::to_string(profile_time.getTimer()) << std::endl;
    //             // applyOnlyRoutePlacement(netsToRoute);
    //             // applyOnlyRoute(netsToRoute);
    //             // netsToRoute.clear();
    //             netsToRoute.clear();
    //             // database.suspected_cells_dict.clear();
    //             // applyFluteRoute3D(netsToRoute);
    //             // for(auto net_idx : netsToRoute){
    //             //     log() << "need to reroute due to violation: "
    //             //           << database.nets[net_idx].getName() << std::endl;
    //             // }
                
    //             // applyOnlyRoutePlacement(netsToRoute);
    //             // netsToRoute.clear();
    //             // if(db::setting.debug){
    //                 printStat();
    //             // }
    //             // profile_time_str << "printStat_" << std::to_string(j) << "_" << \
    //             //     std::to_string(i)<< "," << std::to_string(profile_time.getTimer()) << std::endl;
    //             profile_time_str << "printStat"<< "," << std::to_string(profile_time.getTimer()) << std::endl;
                
    //             // log() << "MEM: cur=" << utils::mem_use::get_current() << "MB, peak=" << utils::mem_use::get_peak() << "MB"
    //             //       << std::endl;
    //             total_wl = 0;
    //             for(auto net : grDatabase.nets){
    //                 total_wl += net.getWirelength();
    //             }
    //             // profile_time_str << "calcTotalWL_" << std::to_string(j) << "_"<< \
    //             //     std::to_string(i)<< "," << std::to_string(profile_time.getTimer()) << std::endl;
    //             profile_time_str << "calcTotalWL" << "," << std::to_string(profile_time.getTimer()) << std::endl;
    //             database.cost_hist.push_back(total_wl);
    //             if(db::setting.debug){
    //                 logReport(report_tracker,j,i);
    //             }

    //             database.logCellLocations(i+1);

    //             investigateRoute(i,netsToRoute);

    //         }//end for
    //         // profile_time_str << "CRP_" << std::to_string(j) << "," << std::to_string(profile_time.getTimer()) << std::endl;
    //         profile_time_str << "CRP" << "," << std::to_string(profile_time.getTimer()) << std::endl;
    //         // 1- apply cugr flute route & maze route
            
    //         if(db::setting.debug){
    //             log() << "last netsToRoute size: "  << netsToRoute.size() << std::endl;
    //             log() << "end global routing + refine Placement..." << std::endl;
    //         }
    //     }//end for iterate the flow
        
    // }//end if-else refineplacement

    // // log() << "nets wirelength..." << std::endl;
    // // for(auto net : grDatabase.nets){
    // //     log() << "net: " << net.getName() << ", wl: " << net.getWirelength() << std::endl;
    // // }
    
    // // postprocessing
    // if(db::setting.postProcessing){
    //     if(db::setting.debug){
    //         log() << "postProcessing..." << std::endl;
    //     }
    //     for (auto& net : grDatabase.nets) {
    //         if(filter_nets.size() != 0){
    //             if(filter_nets.find(net.dbNet.idx) != filter_nets.end()){
    //                 GuideGenerator guideGen(net);
    //                 guideGen.genPatchGuides();
    //             }
    //         }else{
    //             GuideGenerator guideGen(net);
    //             guideGen.genPatchGuides();
    //         }
            
    //     }
    //     if(db::setting.debug){
    //         log() << "end postProcessing..." << std::endl;
    //     }
    // }

    // grDatabase.logNets(iter);
    // grDatabase.logVio(iter);
    // grDatabase.logCongestion(iter);

    // profile_time_str << "postProcessing" << "," << std::to_string(profile_time.getTimer()) << std::endl;

    
    // if(db::setting.debug){    
    //     if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) guideGenStat.print();
    // }

    // if(db::setting.debug){
    //     log() << std::endl;
    //     log() << "################################################################" << std::endl;
    //     database.setUnitVioCost(1);  // set the cost back to without discount
    //     log() << "Finish all RRR iterations and PostRoute" << std::endl;
    //     log() << "MEM: cur=" << utils::mem_use::get_current() << "MB, peak=" << utils::mem_use::get_peak() << "MB"
    //         << std::endl;
    // }

    // printStat();
    // // logReport(report_tracker);
    // if(db::setting.debug){
    //     logReport(report_tracker,-1,-1);
    // }

    // if(db::setting.debug){
    //     std::ofstream fout(db::setting.benchmarkName + "_gr");
    //     fout << report_tracker.str();
    //     fout.close();

    //     std::ofstream fout_refinePlacement(db::setting.benchmarkName + "_refinePlacment");
    //     fout_refinePlacement << database.refinePlacement_report_tracker.str();
    //     fout_refinePlacement.close();

    // }

    // std::ofstream profileTime_file(db::setting.outputFile + ".profileTime.csv");
    // profileTime_file << profile_time_str.str();
    // profileTime_file.close();

    // if(database.db_hpwls.size() > 0)
    //     log() << "improve hpwl: " << double(database.db_hpwls[database.db_hpwls.size() -1 ]) / double(database.db_hpwls[0])  << std::endl;


    // // //-----Test net congestion ----
    // // for(auto net : grDatabase.nets){
    // //     grDatabase.getNetCongestion(net);
    // // }//end for net
    // //-----end of testing net congestion ----
}

// void Router::investigateRoute(int iter_i,std::vector<int>& netsToRoute){
//     std::vector<int> nets_violated; 
//     std::vector<gr::GrEdge> edges_viol; 
//     std::vector<gr::GrPoint> vias_viol; 
//     netsToRoute.clear();
//     for (auto& net : grDatabase.nets){
//          if (grDatabase.getVioReport(net,edges_viol,vias_viol,viol_dict)) {
//         // if (grDatabase.getNumVio(net,false)) {
//             // if(edges_viol.size()>0){
//             //     log() << "violated edge size: " << edges_viol.size() << std::endl;
//             //     for(auto edge : edges_viol){
//             //         log() << "n1: " << edge.lowerGrPoint() 
//             //                 << ", n2: " << edge.upperGrPoint() << std::endl;
//             //         auto lx = grDatabase.getCoor(edge.lowerGrPoint().x, X);
//             //         auto ly = grDatabase.getCoor(edge.lowerGrPoint().y, Y);
//             //         auto hx = grDatabase.getCoor(edge.upperGrPoint().x+1, X);
//             //         auto hy = grDatabase.getCoor(edge.upperGrPoint().y+1, Y);

//             //         DBUxy violation_pt((lx+hx)/2.0,(ly+hy)/2.0);

//             //         log() << "lx: "   << double(lx)/database.libDBU
//             //                 << ", ly: " << double(ly)/database.libDBU
//             //                 << ", hx: " << double(hx)/database.libDBU
//             //                 << ", hy: " << double(hy)/database.libDBU << std::endl;           

//             //         // log() << "edge: " << edge <<", dir: " << dir << std::endl;
//             //     }
//             // }

//             // if(vias_viol.size()>0){
//             //     log() << "violated vias size: " << vias_viol.size() << std::endl;
//             //     for(auto viaTmp : vias_viol){
//             //         log() << "viaTmp: " << viaTmp << std::endl;
//             //     }
//             // }

//             // getSuspectedCellsToViolation(net,edges_viol,vias_viol);
//             netsToRoute.push_back(net.dbNet.idx);
//             log() << "net: " << net.getName() 
//                   << ", has violation in iteration: " << iter_i << std::endl;

            
//         }
//     } 


//     auto net_idx_set = viol_dict[std::make_tuple(3,369,257)];


//     std::set<std::string> nets_viols;

//     for(auto tuple_pair : viol_dict){
//         auto tuple = tuple_pair.first;
//         auto net_idx_set = tuple_pair.second;

//         int layerIdx =  std::get<0>(tuple);
//         int gridline =  std::get<1>(tuple);
//         int cp =  std::get<2>(tuple);
//         log() << "layerIdx: " << layerIdx
//               << ", gridline: " << gridline 
//               << ", cp: " << cp << std::endl;
//         for(auto net_idx: net_idx_set){
//             log() << database.nets[net_idx].getName()  << " ";
//             nets_viols.insert(database.nets[net_idx].getName() );
//             // log() << "ok net: " <<net_idx  << std::endl;
//         }
//         log() << std::endl;
//     }


//     log() << "size of viol_dicts: " << viol_dict.size() << std::endl;
//     log() << "nets need to be routed: " <<  nets_viols.size() << std::endl;
//     for(auto net_str : nets_viols){
//         log() << net_str << std::endl;
//     }
    
    

// }//end investigateRoute

// void Router::getSuspectedCellsToViolation(const gr::GrNet& net
//                               , std::vector<gr::GrEdge>& edges_viol
//                               , std::vector<gr::GrPoint>& vias_viol){
//     std::vector<int> cells_idx;
//     database.getNetCells(net.dbNet.idx,cells_idx);
//     log() << "net with violation: " << net.dbNet.getName() << std::endl;
//     for(auto edge : edges_viol){
//         sort(cells_idx.begin(), cells_idx.end(), [&](int id1, int id2) {
//             auto cell_center_id1 = database.cells[id1].rsynInstance.getBounds().computeCenter();
//             auto cell_center_id2 = database.cells[id2].rsynInstance.getBounds().computeCenter();
//             auto lx = grDatabase.getCoor(edge.lowerGrPoint().x, X);
//             auto ly = grDatabase.getCoor(edge.lowerGrPoint().y, Y);
//             auto hx = grDatabase.getCoor(edge.upperGrPoint().x+1, X);
//             auto hy = grDatabase.getCoor(edge.upperGrPoint().y+1, Y);
//             auto middle_x = (lx+hx)/2.0;
//             auto middle_y = (ly+hy)/2.0;
//             auto md_id1 = std::abs(cell_center_id1.x - middle_x) + std::abs(cell_center_id1.y - middle_y);
//             auto md_id2 = std::abs(cell_center_id2.x - middle_x) + std::abs(cell_center_id2.y - middle_y);
//             return  md_id1 >= md_id2 ;
//         });
//         for(auto cell_idx : cells_idx) {
//             auto cell = database.cells[cell_idx];
//             if(cell.isMultiRowCell()) {
//                 log() << "multiRow cell: " << cell.getName() << std::endl;
//                 continue;
//             }
//             // if(cell.getCellHeight())
//             log() << "suspected cell: " << database.cells[cell_idx].getName() << std::endl;
//             // auto cell_center = database.cells[cell_idx].rsynInstance.getBounds().computeCenter();
//             // min(viaGuides[g1].layerIdx, viaGuides[g2].layerIdx)
//             auto cell_box =  cell.getCellBox();
//             auto cell_boxLayer = db::BoxOnLayer(0,cell_box);
//             // cellMutex.lock();
//             auto grBox = grDatabase.rangeSearchGCell(cell_boxLayer);
        
//             std::vector<int> grBox_xs;
//             std::vector<int> grBox_ys;
            
//             for (int x = grBox[X].low; x <= grBox[X].high; x++)
//                 for (int y = grBox[Y].low; y <= grBox[Y].high; y++) {
//                     grBox_xs.push_back(x);
//                     grBox_ys.push_back(y);
//                         log() << ", gcell_x: " << x
//                             << ", gcell_y: " << y << std::endl;
//                 }

//             int grBox_lx = *std::min_element(grBox_xs.begin(),grBox_xs.end());
//             int grBox_ly = *std::min_element(grBox_ys.begin(),grBox_ys.end());
//             int grBox_hx = *std::max_element(grBox_xs.begin(),grBox_xs.end());
//             int grBox_hy = *std::max_element(grBox_ys.begin(),grBox_ys.end());

//             log() << "grBox_lx: "  << grBox_lx
//                 << ", grBox_ly: "  << grBox_ly
//                 << ", grBox_hx: "  << grBox_hx
//                 << ", grBox_hy: "  << grBox_hy << std::endl;

//             log() << "getCoor(grBox_lx, X): "<< double(grDatabase.getCoor(grBox_lx, X))/database.libDBU << std::endl;
//             log() << "getCoor(grBox_ly, Y): "<< double(grDatabase.getCoor(grBox_ly, Y))/database.libDBU << std::endl;
//             log() << "getCoor(grBox_hx, X): "<< double(grDatabase.getCoor(grBox_hx+1, X))/database.libDBU << std::endl;
//             log() << "getCoor(grBox_hy, Y): "<< double(grDatabase.getCoor(grBox_hy+1, Y))/database.libDBU << std::endl;  

//             utils::BoxT<DBU> penalty_box(grDatabase.getCoor(grBox_lx, X)
//                                 , grDatabase.getCoor(grBox_ly, Y)
//                                 , grDatabase.getCoor(grBox_hx+1, X)
//                                 , grDatabase.getCoor(grBox_hy+1, Y));

//             if(database.suspected_cells_dict.find(cell_idx) == database.suspected_cells_dict.end())
//                 database.suspected_cells_dict[cell_idx] = penalty_box;
//             else{
//                 log() << "alreay have penalty_box cell: " << database.cells[cell_idx].getName() << std::endl;
//             }
//         }
//     }//end for edge_viol
    


// }//end getSuspectedCells

// void Router::applyFluteRoute3D(vector<int>& netsToRoute){
//     if(db::setting.debug){
//         log() << std::endl;
//         log() << "----------------------------------------------------------------" << std::endl;
//         log() << "Start RRR iteration " << iter << std::endl;
//         log() << std::endl;
//     }
//     db::routeStat.clear();
//     guideGenStat.reset();
//     bool debug = false;
//     // netsToRoute.clear();
//     if(debug){
//         // log() << "netsToRoute Size: " <<netsToRoute.size() << std::endl;
//         if(db::setting.debug){
//             log() << "nets to route after movement ..." << std::endl;
//         }
//         // for(auto net : database.nets){
//         //     log() << "net idx: " << net.idx << ", net_name " << net.getName() << std::endl;
//         //     // netsToRoute.push_back(net.idx);
//         // }
//         if(db::setting.debug){
//             for(auto net_idx : netsToRoute){
//                 log() << "net: " << database.nets[net_idx].getName() 
//                     << ",idx: " << net_idx << std::endl;
//             }
//         }

//         // log() << "before ripup wirelength..." << std::endl;
//         // for(auto net : grDatabase.nets){
//         //     log() << "net: " << net.getName() << ", wl: " << net.getWirelength() << std::endl;
//         // }
//     }
    
    

//     if (netsToRoute.empty()) {
//         if(db::setting.debug){
//             if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) {
//                 log() << "No net is identified for this iteration of RRR." << std::endl;
//                 log() << std::endl;
//             }
//         }
//         return;
//     }
//     sortNets(netsToRoute);  // Note: only effective when doing mazeroute sequentially
//     // printCongMap("after_placement_before_ripup.csv");
//     // grDatabase.logRoutedViaMap("routeViaMap_before_ripup.csv");
//     ripup(netsToRoute);        
//     congMap.init(cellWidth, cellHeight);
//     // grDatabase.logRoutedViaMap("routedViaMap_after_gr.csv");
//     // // printCongMap("after_placement_after_ripup.csv");
//     // grDatabase.logRoutedViaMap("routeViaMap_after_ripup.csv");

//     // log() << "before routeApprxCell wirelength..." << std::endl;
//     // for(auto net : grDatabase.nets){
//     //     log() << "net: " << net.getName() << ", wl: " << net.getWirelength() << std::endl;
//     // }

//     routeApprxCell(netsToRoute);   

//     // log() << "after routeApprxCell wirelength..." << std::endl;
//     // for(auto net : grDatabase.nets){
//     //     log() << "net: " << net.getName() << ", wl: " << net.getWirelength() << std::endl;
//     // }

//     if(db::setting.debug){
//         log() << std::endl;
//         log() << "Finish RRR iteration " << iter << std::endl;
//         log() << "MEM: cur=" << utils::mem_use::get_current() << "MB, peak=" << utils::mem_use::get_peak() << "MB"
//             << std::endl;
//         if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) db::routeStat.print();
//     }   


//     updateRouteTable();
// }

void Router::applyPlacement(vector<int>& netsToRoute,int iter_t,utils::timer& profile_time,std::stringstream& profile_time_str){
    netsToRoute.clear();
    Placer place(congMap,routeTable,iter_t, profile_time,profile_time_str);
    place.runMT(netsToRoute,cellWidth, cellHeight);
}

// void Router::applyOnlyRoutePlacement(vector<int>& netsToRoute){
//     bool debug = false;
//     for (iter = 0; iter < db::setting.rrrIterLimit; iter++) {
//         if(db::setting.debug){
//             log() << std::endl;
//             log() << "################################################################" << std::endl;
//             log() << "Start RRR iteration " << iter << std::endl;
//             log() << std::endl;
//         }
//         db::routeStat.clear();
//         guideGenStat.reset();
//         if(debug)
//             log() << "netsToRoute Size: " <<netsToRoute.size() << std::endl;
        
//         // netsToRoute = getNetsToRoute();

//         if(debug){
//             log() << "iteration: " << iter << std::endl;
//             for(auto net_idx : netsToRoute){
//                 log() << "net_name_to_reroute: " << database.nets[net_idx].getName() << std::endl;
//             }
//         }
 
//         if (netsToRoute.empty()) {
//             if(db::setting.debug){
//                 if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) {
//                     log() << "No net is identified for this iteration of RRR." << std::endl;
//                     log() << std::endl;
//                 }
//             }
//             break;
//         }
//         sortNets(netsToRoute);  // Note: only effective when doing mazeroute sequentially

//         updateCost();
//         grDatabase.statHistCost();
//         // grDatabase.logRoutedViaMap("routedViaMap_after_gr.csv");

//         if (iter > 0 ) {
//             ripup(netsToRoute);
//             congMap.init(cellWidth, cellHeight);
//         }
        
        

//         routeApprx(netsToRoute);

//         if(db::setting.debug){
//             log() << std::endl;
//             log() << "Finish RRR iteration " << iter << std::endl;
//             log() << "MEM: cur=" << utils::mem_use::get_current() << "MB, peak=" << utils::mem_use::get_peak() << "MB"
//                 << std::endl;
//             if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) db::routeStat.print();
//         }
//         updateRouteTable();      

//         break;  
//     }
// }// end applyOnlyRoute

void Router::ripupRerouteCUGR(vector<int>& netsToRoute){
    bool debug = false;
    utils::timer profile_time;
    std::stringstream profile_time_str;
    

    // injectCongestion();
    for (int iter = 0; iter < db::setting.rrrIterLimit; iter++) {        
    // for (int iter = 0; iter < 1; iter++) {        
        std::vector<int> netsToRoute;
        printRouteStart("routing",iter);
        routeStateClear();
        getNetsToRoute(netsToRoute,iter);
        sortNets(netsToRoute);  // Note: only effective when doing mazeroute sequentially
        updateCost(iter);
        
        if(filter_routers[iter] == "patternroute" && filter_routers_applys[iter] =="1"){
            ripup(netsToRoute);
            congMap.init(cellWidth, cellHeight);
            routeApprx(netsToRoute, PATTERNROUTE);
            
            updateRouteTable(); 
            printStat();
            logAll();

            // if(db::setting.refinePlacement){
            //     // get all nets with violation after pattern route
            //     netsToRoute.clear();
            //     for (auto& net : grDatabase.nets)
            //         if (grDatabase.hasVio(net)) netsToRoute.push_back(net.dbNet.idx);

            //     if(debug){
            //         for(auto ii : netsToRoute){
            //             log() << "net: " << database.nets[ii].getName() << std::endl;
            //         }
            //     }
                

            //     for(int i = 0; i < db::setting.numRefinePlacement; i++){
            //         printRouteStart("placement",i);
            //         applyPlacement(netsToRoute,i, profile_time, profile_time_str);
            //         route(netsToRoute,PATTERNROUTE);
            //         // route(netsToRoute,ASTAR);
            //         // netsToRoute.clear();
                    
            //         updateRouteTable(); 
            //         printStat();
            //         logAll();
            //     }//end loop
            // }//end refinePlacement

        

        }else if(filter_routers[iter] == "astar" && filter_routers_applys[iter] =="1" ) {
            ripup(netsToRoute);
            congMap.init(cellWidth, cellHeight);
            routeApprx(netsToRoute, ASTAR);
            updateRouteTable(); 
            printStat();
            logAll();
        }

        if(db::setting.refinePlacement && iter == 3){
            // get all nets with violation after pattern route
            netsToRoute.clear();
            for (auto& net : grDatabase.nets)
                if (grDatabase.hasVio(net)) netsToRoute.push_back(net.dbNet.idx);

            if(debug){
                for(auto ii : netsToRoute){
                    log() << "net: " << database.nets[ii].getName() << std::endl;
                }
            }
            

            for(int i = 0; i < db::setting.numRefinePlacement; i++){
                printRouteStart("placement",i);
                applyPlacement(netsToRoute,i, profile_time, profile_time_str);
                route(netsToRoute,PATTERNROUTE);
                // route(netsToRoute,ASTAR);
                // netsToRoute.clear();
                
                updateRouteTable(); 
                printStat();
                logAll();
            }//end loop
        }//end refinePlacement
            

        printRouteEnd("routing",iter);

        

        // break;       
    }
}// end applyOnlyRoute




void Router::ripupRerouteCRP(vector<int>& netsToRoute){
    bool debug = false;
    utils::timer profile_time;
    std::stringstream profile_time_str;
    

    // injectCongestion();
    for (int iter = 0; iter < db::setting.rrrIterLimit; iter++) {        
    // for (int iter = 0; iter < 1; iter++) {        
        std::vector<int> netsToRoute;
        std::vector<int> netsToRouteByPatternRoute;
        std::vector<int> netsToRouteByAStar;

        printRouteStart("routing",iter);
        routeStateClear();
        getNetsToRoute(netsToRoute,iter);
        // sortNets(netsToRoute);  // Note: only effective when doing mazeroute sequentially
        updateCost(iter);

        classifyNets(netsToRoute, netsToRouteByPatternRoute, netsToRouteByAStar );
        
        // if(filter_routers[iter] == "patternroute" && filter_routers_applys[iter] =="1"){
        // route by patternRoute
        ripup(netsToRoute);
        congMap.init(cellWidth, cellHeight);
        routeApprx(netsToRouteByPatternRoute, PATTERNROUTE);
        
        updateRouteTable(); 
        printStat();
        logAll();        

        // }
        // if(filter_routers[iter] == "astar" && filter_routers_applys[iter] =="1" ) {
        // route by Astar
        ripup(netsToRouteByAStar);
        congMap.init(cellWidth, cellHeight);
        routeApprx(netsToRouteByAStar, ASTAR);
        updateRouteTable(); 
        printStat();
        logAll();
        // }

       
            

        printRouteEnd("routing",iter);

        

        break;       
    }
    // return;



    for (int iter = 1; iter < db::setting.rrrIterLimit; iter++) {        
    // for (int iter = 0; iter < 1; iter++) {        
        std::vector<int> netsToRoute;

        printRouteStart("routing",iter);
        routeStateClear();
        getNetsToRoute(netsToRoute,iter);
        sortNets(netsToRoute);  // Note: only effective when doing mazeroute sequentially
        updateCost(iter);

        // classifyNets(netsToRoute, netsToRouteByPatternRoute, netsToRouteByAStar );
        
        // if(filter_routers[iter] == "patternroute" && filter_routers_applys[iter] =="1"){
        // route by patternRoute
        // ripup(netsToRoute);
        // congMap.init(cellWidth, cellHeight);
        // routeApprx(netsToRouteByPatternRoute, PATTERNROUTE);
        
        // updateRouteTable(); 
        // printStat();
        // logAll();        

        // }
        // if(filter_routers[iter] == "astar" && filter_routers_applys[iter] =="1" ) {
        // route by Astar
        ripup(netsToRoute);
        congMap.init(cellWidth, cellHeight);
        routeApprx(netsToRoute, ASTAR);
        updateRouteTable(); 
        printStat();
        logAll();
    }

       
            

    //     printRouteEnd("routing",iter);

        

         
    // }

}// end applyOnlyRoute





void Router::route(const vector<int>& netsToRoute, RouterName routerName){
    log() << "route ..." << std::endl;
    routeStateClear();
    // getNetsToRoute(netsToRoute,iter);
    // sortNets(netsToRoute);  // Note: only effective when doing mazeroute sequentially
    // updateCost(iter);
    
    
    ripup(netsToRoute);
    congMap.init(cellWidth, cellHeight);
    
            
    if(routerName == PATTERNROUTE){
        routeApprx(netsToRoute, PATTERNROUTE);
    }else if (routerName == ASTAR){
        routeApprx(netsToRoute, ASTAR);
    }
    

    // printRouteEnd("routing",iter);

    updateRouteTable(); 
    printStat();
    logAll();
    log() << "end route!" << std::endl;
}//end route



void Router::ripup(const vector<int>& netsToRoute) {
    for (auto id : netsToRoute) {
        grDatabase.removeNet(grDatabase.nets[id]);
        grDatabase.nets[id].gridTopo.clear();
        allNetStatus[id] = db::RouteStatus::FAIL_UNPROCESSED;
    }
}
// void Router::updateCostInit(){
//     grDatabase.update();
//     grDatabase.setLogisticSlope(1);
// }//end updateCostInit


void Router::updateCost(int iter) {
    if (iter > 0) {
        // apply in rrr stages
        grDatabase.addHistCost();
        grDatabase.fadeHistCost();

        grDatabase.setUnitViaMultiplier(
            max(100 / pow(5, iter - 1), 4.0));  // note: enlarge unit via cost to avoid extra use of vias
        grDatabase.setLogisticSlope(db::setting.initLogisticSlope * pow(2, iter));
    }

    if (db::setting.rrrIterLimit > 1) {
        double step = (1.0 - db::setting.rrrInitVioCostDiscount) / (db::setting.rrrIterLimit - 1);
        // if(debug)
        //     log() << "step: " << step << ", iter: " << iter 
        //           << "db::setting.rrrInitVioCostDiscount: " << db::setting.rrrInitVioCostDiscount << std::endl;
        database.setUnitVioCost(db::setting.rrrInitVioCostDiscount + step * iter);
    }

    grDatabase.statHistCost();
}//end updateCost

void Router::updateCostPlacement(int iter) {
    
    grDatabase.setUnitViaMultiplier(4);        
    grDatabase.setLogisticSlope(db::setting.initLogisticSlope);
    database.setUnitVioCost(db::setting.rrrInitVioCostDiscount);
    
}//end updateCost



vector<vector<int>> Router::getBatches(vector<SingleNetRouter>& routers, const vector<int>& netsToRoute) {
    vector<int> batch(netsToRoute.size());
    for (int i = 0; i < netsToRoute.size(); i++) batch[i] = i;

    runJobsMT(batch.size(), [&](int jobIdx) {
        auto& router = routers[batch[jobIdx]];

        const auto mergedPinAccessBoxes = grDatabase.nets[netsToRoute[jobIdx]].getMergedPinAccessBoxes(
            [](const gr::GrPoint& point) { return gr::PointOnLayer(point.layerIdx, point[X], point[Y]); });
        utils::IntervalT<int> xIntvl, yIntvl;
        for (auto& points : mergedPinAccessBoxes) {
            for (auto& point : points) {
                xIntvl.Update(point[X]);
                yIntvl.Update(point[Y]);
            }
        }
        router.guides.emplace_back(0, xIntvl, yIntvl);
    });

    Scheduler scheduler(routers);
    const vector<vector<int>>& batches = scheduler.schedule();

    if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) {
        log() << "Finish multi-thread scheduling" << ((db::setting.numThreads == 0) ? " using simple mode" : "")
              << ". There will be " << batches.size() << " batches for " << netsToRoute.size() << " nets." << std::endl;
        log() << std::endl;
    }

    return batches;
}

void Router::routeApprx(const vector<int>& netsToRoute, RouterName routerName) {
    bool debug = true;


    if (routerName == PATTERNROUTE) {
        // log() << "netsToRoute before flute: " << netsToRoute.size() << std::endl;
        // log() << "before flute netsToRoute size: " << netsToRoute.size() << std::endl;
        fluteAllAndRoute(netsToRoute);
        // log() << "after flute netsToRoute size: " << netsToRoute.size() << std::endl;
        // log() << "netsToRoute after flute: " << netsToRoute.size() << std::endl;
        // grDatabase.logNets(iter);
        // grDatabase.logVio(iter);
        // grDatabase.logCongestion(iter);
    } else if (routerName == ASTAR){
        
        vector<SingleNetRouter> routers;
        routers.reserve(netsToRoute.size());
        for (auto id : netsToRoute) routers.emplace_back(grDatabase.nets[id]);

        
        vector<vector<int>> batches = getBatches(routers, netsToRoute);
        // log() << "netsToRoute before maze: " << netsToRoute.size() << std::endl;
        for (const vector<int>& batch : batches) {
            runJobsMT(batch.size(), [&](int jobIdx) {
                auto& router = routers[batch[jobIdx]];
                router.planMazeRoute(congMap);
            });

            for (auto jobIdx : batch) {
                auto& router = routers[jobIdx];
                router.mazeRoute();
                router.finish();

                int netIdx = netsToRoute[jobIdx];
                congMap.update(grDatabase.nets[netIdx]);
                allNetStatus[netIdx] = router.status;
                // database.astar_stream = database.astar_stream + router.stream_str;
                // database.astar_stream_coarse = database.astar_stream_coarse + router.stream_coarseGrid_str;
            }
        }
        
        // grDatabase.logNets(iter);
        // grDatabase.logVio(iter);
        // grDatabase.logCongestion(iter);
        // log() << "netsToRoute after maze: " << netsToRoute.size() << std::endl;
    }
}

// void Router::routeApprxCell(const vector<int>& netsToRoute) {
    
//     fluteAllAndRouteCell(netsToRoute);
        
// }

void Router::fluteAllAndRoute(const vector<int>& netsToRoute) {
    vector<SingleNetRouter> routers;
    vector<InitRoute> initRouters;
    routers.reserve(netsToRoute.size());
    initRouters.reserve(netsToRoute.size());
    for (auto id : netsToRoute) routers.emplace_back(grDatabase.nets[id]);
    for (auto id : netsToRoute) initRouters.emplace_back(grDatabase.nets[id]);

    grDatabase.init2DMaps(grDatabase);

    for (auto& router : initRouters){
        if (router.grNet.needToRoute()) router.plan_fluteOnly();
    }
        
    if(db::setting.debug){
        printlog("finish planning");
    }

    for (int i = 0; i < db::setting.edgeShiftingIter; i++) {
        grDatabase.edge_shifted = 0;
        for (auto& router : initRouters){
            if(db::setting.escapeRouteBlockage){
                if (router.grNet.needToRoute()) router.edge_shift2d_blockage(router.getRouteNodes());
            }else{
                if (router.grNet.needToRoute()) router.edge_shift2d(router.getRouteNodes());
            }
        }
            
        if(db::setting.debug){
            log() << "Total Number of edges in Iter " << i << " : " << grDatabase.tot_edge
                << ". Edge Shifted: " << grDatabase.edge_shifted << "("
                << 100.0 * grDatabase.edge_shifted / grDatabase.tot_edge << "%)" << std::endl;
        }
    }
    for (auto& router : initRouters)
        if (router.grNet.needToRoute()) router.getRoutingOrder();


    for (int i = 0; i < routers.size(); i++) {
        auto& router = routers[i];
        router.initRoutePattern(initRouters[i]);
        router.finish();
        allNetStatus[netsToRoute[i]] = router.status;

        // pattern route stream for debugging purposes. 
        database.patternRoute_stream = database.patternRoute_stream + initRouters[i].stream.str();
    }
    if(db::setting.debug){
        printlog("finish pattern route");
    }
}

void Router::fluteAllAndRouteCell(const vector<int>& netsToRoute) {
    bool debug= true;

    for(int id : netsToRoute)
        log() << "net to route: " << grDatabase.nets[id].getName() << std::endl;

    vector<SingleNetRouter> routers;
    vector<InitRoute> initRouters;
    routers.reserve(netsToRoute.size());
    initRouters.reserve(netsToRoute.size());
    for (auto id : netsToRoute) routers.emplace_back(grDatabase.nets[id]);
    for (auto id : netsToRoute) initRouters.emplace_back(grDatabase.nets[id],true);

    grDatabase.init2DMaps(grDatabase);

    for (auto& router : initRouters)
        if (router.grNet.needToRoute()) router.plan_fluteOnly();
    if(db::setting.debug){
        printlog("finish planning");
    }

    for (int i = 0; i < db::setting.edgeShiftingIter; i++) {
        grDatabase.edge_shifted = 0;
        for (auto& router : initRouters)
            if (router.grNet.needToRoute()) router.edge_shift2d(router.getRouteNodes());
        if(db::setting.debug){
            log() << "Total Number of edges in Iter " << i << " : " << grDatabase.tot_edge
                << ". Edge Shifted: " << grDatabase.edge_shifted << "("
                << 100.0 * grDatabase.edge_shifted / grDatabase.tot_edge << "%)" << std::endl;
        }
    }
    for (auto& router : initRouters)
        if (router.grNet.needToRoute()) router.getRoutingOrder();
    if(db::setting.debug){
        printlog("finish edge shifting");
    }

    for (int i = 0; i < routers.size(); i++) {
        auto& router = routers[i];
        router.initRoutePattern(initRouters[i]);
        router.finish();
        allNetStatus[netsToRoute[i]] = router.status;

        // getVioReport
        std::vector<gr::GrEdge> edges_viol; 
        std::vector<gr::GrPoint> vias_viol; 
        if(grDatabase.getVioReport(router.grNet,edges_viol,vias_viol,viol_dict)){
            log() << "net: " << router.grNet.getName()
                  << ", has violation!" << std::endl;
        }
        // end getVioReport

    }
    if(db::setting.debug){
        printlog("finish pattern route");
    }
}

void Router::sortNets(vector<int>& netsToRoute) {
    sort(netsToRoute.begin(), netsToRoute.end(), [&](int id1, int id2) {
        return grDatabase.nets[id1].boundingBox.hp() < grDatabase.nets[id2].boundingBox.hp();
    });
}


void Router::classifyNets(vector<int>& netsToRoute, vector<int>& netsToRouteByPatternRoute, vector<int>& netsToRouteByAStar ) {
    for(auto idx : netsToRoute){
        if(grDatabase.nets[idx].routeType == SIMPLE){
            netsToRouteByPatternRoute.push_back(idx);
        }else if(grDatabase.nets[idx].routeType == COMPLEX){
            netsToRouteByAStar.push_back(idx);
        }
    }

    sortNets(netsToRouteByPatternRoute);
    sortNets(netsToRouteByAStar);

}//end classifyNets

void Router::getNetsToRoute(vector<int>& netsToRoute,int iter) {
    if (iter == 0) {
        for (auto& net : grDatabase.nets) {
            if(filter_nets.size() != 0){
                if(filter_nets.find(net.dbNet.idx) != filter_nets.end())
                    netsToRoute.push_back(net.dbNet.idx);   
            }
            else{
                netsToRoute.push_back(net.dbNet.idx);   
            }
        } 
    } else {
        for (auto& net : grDatabase.nets)
            if (grDatabase.hasVio(net)) netsToRoute.push_back(net.dbNet.idx);
    }
}//end getNetsToRoute

void Router::printStat() {
    log() << std::endl;
    log() << "----------------------------------------------------------------" << std::endl;
    db::routeStat.print();
    grDatabase.printAllUsageAndVio();
    double hpwl_total = 0;
    for(auto& grNet : grDatabase.nets){
        auto hpwl_tmp = grNet.getHPWL();
        hpwl_total += hpwl_tmp;
    }
    log() << "total HPWL: " << hpwl_total << std::endl;
    log() << "----------------------------------------------------------------" << std::endl;
    log() << std::endl;
}

void Router::printCongMap(std::string log_name){
    congMap.logCSV(log_name);
}


void Router::logDatabase(std::string name){
    log() << "Writing database to csv file..." << std::endl;

    log() << "database.unitWireCostRaw: " <<  database.unitWireCostRaw << std::endl;
    log() << "database.unitViaCostRaw: " <<  database.unitViaCostRaw << std::endl;
    log() << "database.unitShortVioCostRaw: " <<  database.unitShortVioCostRaw << std::endl;
    log() << "database.unitViaCost: " <<  database.unitViaCost << std::endl;
    
    for (auto tmp : database.unitShortVioCost){
        log() << "unitShortVioCost: " << tmp << std::endl;
    }
    for (auto tmp : database.unitShortVioCostDiscounted){
        log() << "unitShortVioCostDiscounted: " << tmp << std::endl;
    }
   

    std::stringstream ss;
    ss << "layer,idx,lx,ly,hx,hy" << std::endl;

    // std::cout << "get RsrcUsage: " << getRsrcUsage(1,1,1) << std::endl;
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
               << ", " << std::to_string(tmp.second)
               << ", " << std::to_string(box.lx())
               << ", " << std::to_string(box.ly())
               << ", " << std::to_string(box.hx())
               << ", " << std::to_string(box.hy()) << std::endl;
               
        }
    }

    std::ofstream fout(name);
    fout << ss.str();
    fout.close();

}//end logDatabase


void Router::logReport(std::stringstream& ss,int iter_router,int iter_refine_placement){    
    // std::stringstream ss;
    // ss << "name,wl,via,shorts" << std::endl;

    // std::cout << "get RsrcUsage: " << getRsrcUsage(1,1,1) << std::endl;
    double viaNum = grDatabase.getTotViaNum();
    vector<double> buckets = {
        -1, 0, 0.3, 0.6, 0.8, 0.9, 1, 1.1, 1.3, 1.5, 2, 3};  // the i-th bucket: buckets[i] <= x < buckets[i+1]
    // Wire
    vector<int> routedWireUsageGrid;
    vector<DBU> routedWireUsageLength;
    double wireLength = grDatabase.getAllWireUsage(buckets, routedWireUsageGrid, routedWireUsageLength);


    wireLength /= double(database.getLayer(1).pitch);

    double numShort = grDatabase.getAllVio();


    vector<std::string> items = {"wirelength", "# vias", "short"};
    vector<double> metrics = {wireLength, viaNum, numShort};
    vector<double> weights = {db::setting.weightWirelength, db::setting.weightViaNum, db::setting.weightShortArea};
    double totalScore = 0;
    for (int i = 0; i < items.size(); ++i) {
        totalScore += metrics[i] * weights[i];
    }

    std::string name = db::setting.benchmarkName;
    ss << name 
        << ", " << std::to_string(wireLength)
        << ", " << std::to_string(viaNum)
        << ", " << std::to_string(numShort)
        << ", " << std::to_string(iter_router)
        << ", " << std::to_string(iter_refine_placement)
        << ", " << std::to_string(totalScore) << std::endl;               
    

}


void Router::updateRouteTable(){
    for(auto& grnet : grDatabase.nets){
        routeTable[grnet.dbNet.idx].push_back(grnet.getPathCost());
        // routeTable[grnet.dbNet.idx].push_back(grnet.getWirelength());
        grnet.gridTopoHist.push_back(grnet.gridTopo);
    }

}

void Router::logRouteTable(std::string name){
    std::stringstream ss;
    int num_iter = 0;
    for(auto tmp : routeTable){
        num_iter = tmp.size();
        break;
    }
    std::string title_string = "nets,";
    ss << title_string;
    for(int i = 0; i < num_iter; i++){
        ss << "itr" << std::to_string(i);
        if(i != num_iter-1)
            ss << ",";
    }
    ss << std::endl;

    for(int i = 0; i < routeTable.size(); i++){
        auto net = database.nets[i];
        ss << net.getName() << ",";
        for(int j = 0 ; j < routeTable[i].size(); j++){
            ss << std::to_string(routeTable[i][j]);
            if(j != routeTable[i].size()-1)
                ss << ",";
        }
        ss << std::endl;
        
    }//end for 
    

    std::ofstream fout(name);
    fout << ss.str();
    fout.close();

}//end logDatabase


void Router::logNetsToRoute(vector<int>& netsToRoute){
    if(db::setting.debug){
        log() << "Log netsToBeRoute..." << std::endl;
    }

    std::stringstream ss;
    ss << "net" << std::endl;    
    for(auto net_idx : netsToRoute){
        ss << database.nets[net_idx].getName() << std::endl;
            
    }

    std::string file_name_csv = db::setting.outputFile + "_netToRoute.csv";
    std::ofstream fout(file_name_csv);
    fout << ss.str();
    fout.close();
}//end logNetsToRoute

void Router::filterNets(){
    bool debug = false;
    filter_nets.clear();
    // temp I will fix this part later 
    vector<std ::string> filter_nets_name;
    vector<std ::string> filter_nets_NumPin;
    set<std ::string> filter_nets_name_set;
    set<int> filter_nets_NumPin_set;
    boost::split(filter_nets_name, db::setting.filter_nets_name, boost::is_any_of(","));
    boost::split(filter_nets_NumPin, db::setting.filter_nets_NumPin, boost::is_any_of(","));


    for(auto net_name : filter_nets_name){
       filter_nets_name_set.insert(net_name);
    }

    for(auto num_pin : filter_nets_NumPin){
       filter_nets_NumPin_set.insert(std::stoi( num_pin ));
    }

    bool isEmptyFilterNetsName = filter_nets_name_set.size() == 0 ? true : false;
    bool isEmptyFilterNetsNumPin = filter_nets_NumPin_set.size() == 0 ? true : false;

    // filter nets
    for(auto net : database.nets){
        if(!isEmptyFilterNetsName || !isEmptyFilterNetsNumPin){
            if(filter_nets_name_set.find(net.getName()) != filter_nets_name_set.end()){
                filter_nets.insert(net.idx);
            }
                
            if(filter_nets_NumPin_set.find(net.numOfPins()) != filter_nets_NumPin_set.end()){
                filter_nets.insert(net.idx);
            }
                
        }//end empty check

    }//end for net 
}//end filterNets

// void Router::logCoef(int iter){
//     // coef_stream << std::to_string(iter)
//     //             << "," << std::to_string(database.unitWireCostRaw)
//     //             << "," << std::to_string(database.unitViaCostRaw)
//     //             << "," << std::to_string(database.unitShortVioCostRaw)
//     //             << "," << std::to_string(0) //database.unitShortVioCost)
//     //             << "," << std::to_string(0) //database.unitShortVioCostDiscounted)
//     //             << "," << std::to_string(database.unitViaCost)
//     //             << "," << std::to_string(0) // step)
//     //             << "," << std::to_string(db::setting.rrrInitVioCostDiscount )
//     //             << "," << std::to_string(db::setting.rrrIterLimit )
//     //             << "," << std::to_string(db::setting.initLogisticSlope )
//     //             << "," << std::to_string(db::setting.rrrFadeCoeff )
//     //             << "," << std::to_string(grDatabase.unitViaMultiplier)
//     //             << "," << std::to_string(grDatabase.logisticSlope) 
//     //             << "," << std::to_string(grDatabase.wireCapDiscount) 
//     //             << std::endl;
    
// }//end logCoef


void Router::gridMapReport(){
    std::stringstream ss;

    ss << "l,gridline,cp,wireusage,fixedusage,tracks,elx,ely,ehx,ehy" << std::endl;
    // \
    // ,lx,ly,hx,hy" << std::endl;


    // write instances to csv file
    for (int layerIdx = 0; layerIdx < database.getLayerNum(); ++layerIdx) {
        Dimension dir = database.getLayerDir(layerIdx);

        for (int gridline = 0; gridline < grDatabase.getNumGrPoint(dir); gridline++) {
            for (int cp = 0; cp < grDatabase.getNumGrEdge(layerIdx); cp++) {
                // double numWire = getWireUsage(layerIdx, gridline, cp);
                // double usage = (numWire + getFixedUsage({layerIdx, gridline, cp})) / getNumTracks(layerIdx, gridline);
                gr::GrEdge tempEdge(layerIdx,gridline,cp);  

                // log() << "l: " <<layerIdx
                //       << ", gridline: " << gridline
                //       << ", cp: " << cp
                //       << ", e.lGrP: " << tempEdge.lowerGrPoint()
                //       << ", e.uGrP: " << tempEdge.upperGrPoint()
                //       << ", e.lx: " << grDatabase.getCoor(tempEdge.lowerGrPoint().x, X)/database.libDBU
                //       << ", e.ly: " << grDatabase.getCoor(tempEdge.lowerGrPoint().y, Y)/database.libDBU
                //       << ", e.hx: " << grDatabase.getCoor(tempEdge.upperGrPoint().x+1, X)/database.libDBU
                //       << ", e.hy: " << grDatabase.getCoor(tempEdge.upperGrPoint().y+1, Y)/database.libDBU
                //       << std::endl;

                auto elx = grDatabase.getCoor(tempEdge.lowerGrPoint().x, X)/database.libDBU;
                auto ely = grDatabase.getCoor(tempEdge.lowerGrPoint().y, Y)/database.libDBU;
                auto ehx = grDatabase.getCoor(tempEdge.upperGrPoint().x+1, X)/database.libDBU;
                auto ehy = grDatabase.getCoor(tempEdge.upperGrPoint().y+1, Y)/database.libDBU;


                ss << std::to_string(layerIdx)
                   << "," << std::to_string(gridline)
                   << "," << std::to_string(cp)
                   << "," << std::to_string(grDatabase.getWireUsage(layerIdx, gridline, cp))
                   << "," << std::to_string(grDatabase.getFixedUsage({layerIdx, gridline, cp}))
                   << "," << std::to_string(grDatabase.getNumTracks(layerIdx, gridline))
                   << "," << std::to_string(elx)
                   << "," << std::to_string(ely)
                   << "," << std::to_string(ehx)
                   << "," << std::to_string(ehy) << std::endl;
                //    << std::to_string(lx)
                //    << std::to_string(ly)
                //    << std::to_string(hx)
                //    << std::to_string(hy) << std::endl;
            }
        }
    }


    std::ofstream file(db::setting.outputFile + ".grid.csv");
    file << ss.str();
    file.close();
}

void Router::visualiseCircuit(){
    // grDatabase.printGrid();
    // gridMapReport();
}//end visualiseCircuit


// void Router::routeAStarSeq(vector<int>& netsToRoute){
//     // iter++;
//     if(!grDatabase.hasVio(grDatabase.nets[netsToRoute[0]])) return;

//     sortNets(netsToRoute);  // Note: only effective when doing mazeroute sequentially

//     updateCost();
//     grDatabase.statHistCost();

//     if(db::setting.debug)
//         logCoef();
//     // grDatabase.logRoutedViaMap("routedViaMap_after_gr.csv");

    
//     ripup(netsToRoute);
//     congMap.init(cellWidth, cellHeight);     
    

//     vector<SingleNetRouter> routers;
//     routers.reserve(netsToRoute.size());
//     for (auto id : netsToRoute) routers.emplace_back(grDatabase.nets[id]);


//     for(auto router : routers){
//         router.planMazeRoute(congMap);
//         router.mazeRoute();
//         router.finish();
//         int netIdx = router.grNet.dbNet.idx;
//         congMap.update(netIdx);
//         allNetStatus[netIdx] = router.status;
//     }

     
                
            




// }//end routeAStarSeq

void Router::logAll(){
    if(db::setting.logAll){
        log() << "log nets/cells/vio/congestion ..." << std::endl;
        if(log_iter == 0){
            database.logDie();
            database.logLayers();
            grDatabase.logGCellGrid();
        }
        
        database.logCellLocations(log_iter);
        database.logFixedMetals(log_iter);
        grDatabase.logNets(log_iter);
        grDatabase.logVio(log_iter);
        // grDatabase.logCongestion(log_iter);
        grDatabase.logCoef(log_iter);
        database.logPatternRoute(log_iter);
        database.logAstar(log_iter);
        database.logAstarCoraseGrid(log_iter);
        log_iter++;
        log() << "end log!" << std::endl;
    }
}


void Router::printRouteStart(std::string txt,int iter){
    log() << std::endl;
    log() << "################################################################" << std::endl;
    log() << txt << " , iter: " << iter << std::endl;
    log() << std::endl;
}//end printRouteStart
void Router::printRouteEnd(std::string txt,int iter){
    log() << std::endl;
    log() << "Finish RRR iteration " << iter << std::endl;
    log() << "MEM: cur=" << utils::mem_use::get_current() << "MB, peak=" << utils::mem_use::get_peak() << "MB"
        << std::endl;
    if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) db::routeStat.print();
}


void Router::routeStateClear(){
    db::routeStat.clear();
    guideGenStat.reset();
}//end routeStateClear


void Router::injectCongestion(){

    struct fakeCong{
        utils::IntervalT<int> xx;
        utils::IntervalT<int> yy;
        int l; //layer
        int usage;
    };

    std::vector<fakeCong> congestion = {
        //  xx,yy,l,usage
        {utils::IntervalT<int>(5, 7),utils::IntervalT<int>(9,9),2,6},
        {utils::IntervalT<int>(5, 7),utils::IntervalT<int>(8,8),2,6},
        {utils::IntervalT<int>(5, 7),utils::IntervalT<int>(7,7),2,6},
        {utils::IntervalT<int>(5, 7),utils::IntervalT<int>(6,6),2,6},
        {utils::IntervalT<int>(5, 7),utils::IntervalT<int>(5,5),2,6},
        {utils::IntervalT<int>(5, 6),utils::IntervalT<int>(4,4),2,6},
        {utils::IntervalT<int>(5, 6),utils::IntervalT<int>(3,3),2,6},
        {utils::IntervalT<int>(5, 6),utils::IntervalT<int>(2,2),2,6},
        {utils::IntervalT<int>(5, 6),utils::IntervalT<int>(1,1),2,6},
        {utils::IntervalT<int>(5, 6),utils::IntervalT<int>(0,0),2,6}
    };

    
    for (auto cong : congestion){      
        gr::GrBoxOnLayer grBox(cong.l,cong.xx,cong.yy);
        for(int i = 0; i < cong.usage; i++){
            grDatabase.useWire(grBox);        
        }
    }
    

}//end fakeViolationGenerator


void Router::netFeatureExtraction2D(vector<int>& netsToRoute){
    vector<SingleNetRouter> routers;
    vector<InitRoute> initRouters;
    routers.reserve(netsToRoute.size());
    initRouters.reserve(netsToRoute.size());

    grDatabase.init2DMaps(grDatabase);

    for (auto id : netsToRoute) routers.emplace_back(grDatabase.nets[id]);
    for (auto id : netsToRoute) initRouters.emplace_back(grDatabase.nets[id]);
    

    runJobsMT(initRouters.size(), [&](int i) { initRouters[i].plan_fluteOnly(); });
    runJobsMT(initRouters.size(), [&](int i) { initRouters[i].extractFeature(initRouters[i].getRouteNodes()); });
           
    
    
    int i = 0; 
    int j = 0;
    for(auto& router: initRouters){
        if(router.routeType==SIMPLE){
            i++;
            router.grNet.routeType = SIMPLE;
        }else if(router.routeType==COMPLEX){
            j++;
            router.grNet.routeType = COMPLEX;
        }
    }
    log() << "Total: " << netsToRoute.size() 
          << ", SIMPLE: " << i 
          << ", COMPLEX: " << j 
          << std::endl;


    
    
            

}//end netFeatureExtraction2D