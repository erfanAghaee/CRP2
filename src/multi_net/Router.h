#pragma once

#include "db/Database.h"
#include "gr_db/GrDatabase.h"
#include "single_net/SingleNetRouter.h"
#include "CongestionMap.h"

class Router {
public:
    Router();
    void run();
    // void runISPD();

    void printCongMap(std::string log_name);
private:
    int iter;
    
    vector<db::RouteStatus> allNetStatus;
    vector<vector<int>> routeTable;
    CongestionMap congMap;
    int cellWidth = 5, cellHeight = 5;
    // int cellWidth = 1, cellHeight = 1;

    void getNetsToRoute(vector<int>& netsToRoute,int iter);
    void sortNets(vector<int>& netsToRoute);
    void classifyNets(vector<int>& netsToRoute, vector<int>& netsToRouteByPatternRoute, vector<int>& netsToRouteByAStar );
    vector<vector<int>> getBatches(vector<SingleNetRouter>& routers, const vector<int>& netsToRoute);

    void routeApprx(const vector<int>& netsToRoute, RouterName routeName);
    // void routeApprxCell(const vector<int>& netsToRoute);
    // void routeAStarSeq( vector<int>& netsToRoute);
    void fluteAllAndRoute(const vector<int>& netsToRoute);
    void fluteAllAndRouteCell(const vector<int>& netsToRoute);

    void ripup(const vector<int>& netsToRoute);
    void updateCost(int iter);
    void updateCostPlacement(int iter);

    // void updateCostInit();
    void updateRouteTable();

    void ripupRerouteCUGR(vector<int>& netsToRoute);
    void ripupRerouteCRP(vector<int>& netsToRoute);
    void route(const vector<int>& netsToRoute, RouterName routeName);
    // void applyOnlyRoutePlacement(vector<int>& netsToRoute);
    void applyPlacement(vector<int>& netsToRoute,int iter_t,utils::timer& profile_time,std::stringstream& profile_time_str);
    void applyFluteRoute3D(vector<int>& netsToRoute);
    // void investigateRoute(int iter_i,std::vector<int>& netsToRoute);

    void printStat();

    // void getSuspectedCellsToViolation(const gr::GrNet& net
    //                           , std::vector<gr::GrEdge>& edges_viol
    //                           , std::vector<gr::GrPoint>& vias_viol);


    void logDatabase(std::string name);
    void logRouteTable(std::string name);
    void logReport(std::stringstream& ss,int iter_router,int iter_refine_placement);
    void logNetsToRoute(vector<int>& netsToRoute);
    void netFeatureExtraction2D(vector<int>& netsToRoute);

    void gridMapReport();


    void logAll();


    void visualiseCircuit();

    void filterNets();
    

    std::set<int> filter_nets;
    std::vector<std::string> filter_routers;
    std::vector<std::string> filter_routers_applys;


    // void logCoef();
    std::stringstream coef_stream;
    
    // l, gridline, cp -> vector of net_idx
    std::unordered_map<std::tuple<int,int,int>,std::set<int>> viol_dict;


    void printRouteStart(std::string txt, int iter);
    void printRouteEnd(std::string txt, int iter);
    void routeStateClear();

    void injectCongestion();
  

    int log_iter=0;
    

    
    
};
