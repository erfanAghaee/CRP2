#pragma once

#include "global.h"
#include "flute/flute.h"
#include "gr_db/GrDatabase.h"
#include "GenGuide.h"
#include <map>
namespace routeCell{

extern "C" {
Tree flute(int d, DTYPE x[], DTYPE y[], int acc);
}

class RouteNode {
public:
    int x;
    int y;
    int idx;
    vector<int> pinLayers;    // pins' layers within the node
    vector<int> pinIdxs;      // pins' indexes within the node
    std::set<int> toConnect;  // nodes connecting to this one
    vector<int> childIdxs;

    // PatternRoute
    vector<db::CostT> exitCosts;  // cost exiting the node from different layers
    // cost of entering the node on different layers from different children
    std::unordered_map<int, vector<db::CostT>> enterCosts;  // (childIdx->layerIdx->cost)
    // topo of entering the node on different layers from different children
    std::unordered_map<int, vector<vector<gr::GrPoint>>> enterEdges;  // (childIdx->layerIdx->edge)
    // which layer should a child edge enter if exiting from a specific layer
    std::unordered_map<int, vector<int>> exitEnterLayers;  // (childIdx->exit layer-> enter layer)

    int degree() { return toConnect.size(); }

    std::string str() {  // return string description
        return std::to_string(idx) + "(" + std::to_string(x) + ", " + std::to_string(y) + ")" +
               (pinIdxs.size() == 0 ? "Steiner" : " ");
    }

    friend inline std::ostream& operator<<(std::ostream& os, const RouteNode& rhs) {
        os << "[x: " << rhs.x << ", y: " << rhs.y << ", idx: " << rhs.idx << "]";
        return os;
    }

    RouteNode() : x(-1), y(-1) {}
    RouteNode(int nx, int ny) : x(nx), y(ny) {}
    RouteNode(std::tuple<int, int> loc) : x(std::get<0>(loc)), y(std::get<1>(loc)) {}
    
    void printExitCost(){
        log() << "printExitCost..."<< std::endl;
        for(int i =0; i < exitCosts.size() ; i++){
            log()   << "idx: " << idx 
                    << ", exitCosts[" << i << "]: " << exitCosts[i] << std::endl;
        }
    }
    void printEnterCosts(){
        log() << "printEnterCosts..."<< std::endl;
        for(auto enterCostTmp : enterCosts){
            auto idx = enterCostTmp.first;
            auto cost_vec = enterCostTmp.second;
            for(int j = 0; j < cost_vec.size() ; j++){
                log() << "idx: " << idx << ", enterCosts[" << j
                      << "]: " << cost_vec[j] << std::endl;
            }
        }
    }
    void printExitEnterLayers(){
        log() << "printExitEnterLayers..."<< std::endl;
        for(auto exitEnterLayersTmp : exitEnterLayers){
            auto idx = exitEnterLayersTmp.first;
            auto cost_vec = exitEnterLayersTmp.second;
            for(int j = 0; j < cost_vec.size() ; j++){
                log() << "idx: " << idx << ", exitEnterLayers[" << j
                      << "]: " << cost_vec[j] << std::endl;
            }
        }
    }
    void printEnterEdges(){
        log() << "printEnterEdges..."<< std::endl;
        for(auto enterEdgeTmp : enterEdges){
            auto idx = enterEdgeTmp.first;
            auto mat_gr_pts = enterEdgeTmp.second;
            for(auto vec_gr_pts : mat_gr_pts){
                for(auto gr_pts : vec_gr_pts){
                    log() << "pts: " << gr_pts << std::endl;
                }
            }
        }
    }


    void print(){
        log() << "exitCosts.size(): " << exitCosts.size() << std::endl;
        log() << "enterCosts: " << enterCosts.size() << std::endl;
        log() << "enterEdges: " << enterEdges.size() << std::endl;
        log() << "exitEnterLayers: " << exitEnterLayers.size() << std::endl;
        log() << str() << std::endl;
        printExitCost();
        printEnterCosts();
        printExitEnterLayers();
        printEnterEdges();

        
    }//end print
};

struct RouteEdge {
    int from;
    int to;
};

class InitRouteCell {
public:
    InitRouteCell(gr::GrNet &grNetData,std::vector<utils::BoxT<DBU>>& cellBoxsTmp, vector<vector<db::BoxOnLayer>>& pinAccessBoxesTmp ,
    bool relax_m = false,bool debug = false) :
     grNet(grNetData),cellBoxs(cellBoxsTmp),pinAccessBoxes(pinAccessBoxesTmp)
     ,debug(debug) {
        // grNet.gridTopo.clear();
        cong_cost = 0;
        wl_cost = 0;
        via_cost = 0;
        via_abs_cost = 0;
        wl_abs_cost = 0;
        path_cost = 0;
        relax = relax_m;
    }

    void patternRoute();  // pattern routing
    void patternRouteMT();  // pattern routing
    void patternRouteMemo();  // pattern routing
    
    bool buildTopo();

    void plan_fluteOnly();
    // void costEstimator(std::map<int, RouteNode> &routeNodes);
    void getRoutingOrder();
    static void sliceGuides(vector<gr::GrBoxOnLayer>& guides, bool mergeAdj = false);
    // void addUsage2D(RouteNode &u, RouteNode &v, double usage = 1);
//     void removeUsage2D(RouteNode &u, RouteNode &v, double usage = 1);

//     std::map<int, RouteNode> &getRouteNodes();

//     db::RouteStatus status;
    gr::GrNet &grNet;
    std::vector<utils::BoxT<DBU>>& cellBoxs;

// private:
//     GuideGenerator guideGen;
    std::map<int, RouteNode> routeNodes;
    vector<RouteEdge> routeEdges;

    db::CostT getBufferedWireCost(gr::GrEdge edge);
    db::CostT getBufferedWireCostRelax(gr::GrEdge edge, std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
                                    , std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls);
    std::unordered_map<gr::GrEdge, db::CostT> wireCostBuffer;

    void LShape(const RouteEdge &edge, std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
                                    , std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls);

    void runFlute();
    void runFluteV2();

    void getGrPinAccessBoxes(const gr::GCellGrid& gcellGrid);
    DBU wl_cost;
    DBU cong_cost;
    DBU via_cost;
    DBU via_abs_cost;
    DBU wl_abs_cost;
    DBU path_cost;
    bool debug;
    vector<vector<db::BoxOnLayer>>& pinAccessBoxes; 
    vector<vector<gr::GrPoint>> grPinAccessBoxes;
    bool relax;
    
};


};

