#pragma once

#include "global.h"
#include "flute/flute.h"
#include "gr_db/GrDatabase.h"
#include "GenGuide.h"
#include <map>




struct hash_tuple {  // hash binary tuple
    template <class T>
    size_t operator()(const tuple<T, T>& tup) const {
        auto hash1 = hash<T>{}(get<0>(tup));
        auto hash2 = hash<T>{}(get<1>(tup));
        return hash1 ^ hash2;
    }
};

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

class InitRoute {
public:
    InitRoute(gr::GrNet &grNetData, bool relax_m = false) : grNet(grNetData), guideGen(grNet), status(db::RouteStatus::SUCC_NORMAL) {
        grNet.gridTopo.clear();
        debug = false;
        relax = relax_m;
        // if(grNet.getName() == "net1233") debug = false;
    }

    void patternRoute();  // pattern routing
    void patternRouteMT();  // pattern routing
    void buildTopo();

    void plan_fluteOnly();
    void edge_shift2d(std::map<int, RouteNode> &routeNodes);
    // try to avoid blockage
    void edge_shift2d_blockage(std::map<int, RouteNode> &routeNodes);
    void getRoutingOrder();
    void addUsage2D(RouteNode &u, RouteNode &v, double usage = 1);
    void removeUsage2D(RouteNode &u, RouteNode &v, double usage = 1);

    std::map<int, RouteNode> &getRouteNodes();

    db::RouteStatus status;
    gr::GrNet &grNet;

    // stream for debugging
    std::stringstream stream;
    std::stringstream stream_time;
    unordered_map<tuple<int, int>, vector<int>, hash_tuple> loc2Pins;
    float net_ctrx;
    float net_ctry;
    

private:
    GuideGenerator guideGen;
    std::map<int, RouteNode> routeNodes;
    vector<RouteEdge> routeEdges;

    db::CostT getBufferedWireCost(gr::GrEdge edge,bool debug=false);
    std::unordered_map<gr::GrEdge, db::CostT> wireCostBuffer;

    void LShape(const RouteEdge &edge);

    void runFlute();
    std::pair<float,float> getNetCenter();
    void getPinCenter(vector<tuple<int, int>>& pinCenters, float net_ctrx,float net_ctry);
    void constructRouteNodes(Tree& flutetree
                            , int degree
                            , unordered_map<tuple<int, int>, vector<int>, hash_tuple>& loc2Pins
                            , int node_cnt);
    void routeNodePostProcessing(Tree& flutetree
                            , int degree
                            , unordered_map<tuple<int, int>, vector<int>, hash_tuple>& loc2Pins
                            , int node_cnt);

    
    void getLoc2Pins( unordered_map<tuple<int, int>, vector<int>, hash_tuple>& loc2Pins
                    , vector<tuple<int, int>>& pinCenters
                    );
    
    void logRouteNodes();
    
    void logRoute();

    bool debug;
    bool relax;

    
};
