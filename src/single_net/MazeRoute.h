#pragma once

#include "GridGraph.h"
#include "multi_net/CongestionMap.h"

enum MAZE_TYPE { COARSEGRID = 0, FINEGRID = 1};

class Solution {
public:
    db::CostT cost;
    int vertex;
    std::shared_ptr<Solution> prev;

    Solution(db::CostT c, int v, const std::shared_ptr<Solution> &p) : cost(c), vertex(v), prev(p) {}

    friend ostream &operator<<(ostream &os, const Solution &sol);
};

class MazeRoute {
public:
    MazeRoute(gr::GrNet &grNetData, MAZE_TYPE type_in) 
        : grNet(grNetData) 
        , type(type_in)
    {
        iter = 0;
        int cellWidth = 1;
        int cellHeight = 1;
    }

    void constructGridGraph(const vector<gr::GrBoxOnLayer> &guides);
    void constructGridGraph(const CongestionMap& congMap);
    db::RouteStatus run();

    int iter;

    std::string getGridGraphStreamString(){return graph.stream.str();}
    std::string getCoarseGridGraphStreamString(){return graph.stream_coarse.str();}

private:
    gr::GrNet &grNet;
    GridGraph graph;
    // this can be finegrid or corasegrid
    MAZE_TYPE type;

    vector<db::CostT> vertexCosts;              // min cost upper bound for each vertex
    vector<std::shared_ptr<Solution>> pinSols;  // best solution for each pin
    vector<vector<gr::PointOnLayer>> mergedPinAccessBoxes;

    db::RouteStatus route(int startPin);
    void getResult();
    int cellWidth; 
    int cellHeight; 
};
