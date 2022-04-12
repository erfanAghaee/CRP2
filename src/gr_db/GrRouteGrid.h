#pragma once

#include "GCell.h"
#include "GrGeoPrimitive.h"
#include "GrNet.h"

namespace gr {
class GrRouteGrid : public GCellGrid {
public:
    using UsageT = double;
    using UsageMapT = vector<vector<vector<UsageT>>>;
    int edge_shifted =0;
    int tot_edge = 0;

    void init();
    void update();
    void clear();
    // void reset();

    void useNet(const GrNet& net);
    void removeNet(GrNet& net);
    void removeNetRelax(GrNet& net,std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
                        ,std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls);

    db::CostT getViaCost(const GrPoint& via,bool debug = false) const;
    db::CostT getViaCostRelax(const GrPoint& via
                              ,std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
                              ,std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls
                              ,bool debug = false) const;
    db::CostT getStackViaCost(const GrPoint& via, int height,bool debug = false) const;
    db::CostT getStackViaCostRelax(const GrPoint& via, int height
                                ,std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
                                ,std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls
                                    ,bool debug = false) const;
    db::CostT getWireCost(const GrEdge& edge,bool debug = false) const;
    db::CostT getWireCostRelax(const GrEdge& edge
                                ,std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
                                ,std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls
                                ,bool debug = false) const;
    db::CostT getHistCost(int layerIdx, int gridline, int cp) const;

    double getCellResource(const GrPoint& point) const;  // Note: simplified version

    bool hasVio(const GrNet& net, bool hasCommit = true) const;
    bool hasVio(const GrEdge& edge, bool hasCommit = true) const;
    bool hasVioRelax(const GrEdge& edge,std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
                    ,std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls, bool hasCommit = true) const;
    bool hasVio(const GrPoint& via, bool hasCommit = true) const;

    void print() const;
    void printAllUsageAndVio() const;

    double getWirelength() const;

    void setViaCapDiscount(double discount = 1);
    void setWireCapDiscount(double discount = 1);
    void setUnitViaMultiplier(double multiplier = 1);
    void setLogisticSlope(double slope = 1);

    double getUnitViaMultiplier(){return unitViaMultiplier; }

    double getLogisticSlope() const { return logisticSlope; }

    // for ripup and reroute
    void addHistCost();
    void fadeHistCost();
    void statHistCost() const;

    friend class GrRouteGrid2D;

public:
    // only pref-dir is modeled
    UsageMapT routedWireMap;  // model cross-cell routing
    UsageMapT routedViaMap;
    vector<vector<vector<std::pair<int, DBU>>>> fixedMetalMap;  // layer, x, y, (# blocked tracks, avg blocked length)
    vector<vector<vector<double>>> histWireUsageMap;


    db::CostT getWireDistCost(const GrEdge& edge) const;
    db::CostT getWireShortCost(const GrEdge& edge, bool debug=false) const;
    db::CostT getWireShortCostRelax(const GrEdge& edge
                                    , std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
                                    , std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls
                                    , bool debug=false) const;
    db::CostT getViaShortCost(const GrPoint& via, bool debug=false) const;
    db::CostT getViaShortCostRelax(const GrPoint& via
                                    , std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
                                    , std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls
                                    , bool debug=false) const;

    double viaCapDiscount = 1;
    double wireCapDiscount = 1;
    double unitViaMultiplier = 1;
    double logisticSlope = 1;
public:
    std::pair<double, double> printAllUsage() const;
    double printAllVio() const;
    double getAllVio() const;

    double getAllWireUsage(const vector<double>& buckets,
                           vector<int>& wireUsageGrid,
                           vector<DBU>& wireUsageLength) const;
    void getAllInCellUsage(const vector<double>& buckets, vector<int>& viaUsage) const;
    double getTotViaNum() const;
public:
    void markFixedMetals();
public:
    void markFixed(int layerIdx, int gridline, int cp, int num_track, DBU avg_length);

    bool getVioReport(const GrNet& net,std::vector<GrEdge>& edges_viol
                              , std::vector<GrPoint>& vias_viol
                              , std::unordered_map<std::tuple<int,int,int>,std::set<int>>& viol_dict
                              , bool hasCommit=true) const;
    // bool getVioReportRelax(const GrNet& net,std::vector<GrEdge>& edges_viol
    //                           , std::vector<GrPoint>& vias_viol, bool hasCommit=true) const;
    bool getVioReportRelax(const GrNet& net
                           , std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
                           , std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls 
                           , vector<gr::GrBoxOnLayer>& wireRouteGuides
                           , vector<gr::GrBoxOnLayer>& viaRouteGuides) const;                              
                              

    double getNumVio(const GrNet& net, bool hasCommit) const;
    double getNetCongestion(const GrNet& net) const;
    double getEdgeCongestion(const GrEdge& edge) const;
    double getGrPointCongestion(const GrPoint& via) const;
    double getNumVio(const GrEdge& edge, double selfUsage) const;
    double getNumVioRelax(const GrEdge& edge
                          ,std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
                          ,std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls, double selfUsage) const;
    double getNumVio(const GrPoint& via, double selfUsage) const;
    // double getNumVioRelax(const GrPoint& via,std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls) const;

    double getFixedUsage(int layerIdx, int gridline, int cp) const;
    double getFixedUsage(const GrEdge& edge) const;
    double getWireUsage(int layerIdx, int gridline, int cp) const;
    double getWireUsage(const GrEdge& edge) const;
    double getWireUsageRelax(int layerIdx, int gridline, int cp
            , std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls) const;
    double getWireUsageRelax(const GrEdge& edge, std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls) const;

    double getViaUsage(int layerIdx, int x, int y) const;
    double getViaUsage(const GrPoint& via) const;
    
    double getViaUsageRelax(int layerIdx, int x, int y,std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls) const;
    double getViaUsageRelax(const GrPoint& via,std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls) const;

    double getCellResource(int layerIdx, int x, int y) const;

    double getWireCapacity(const GrEdge& edge) const;
    double getInCellArea(const GrPoint& point) const;

    double getInCellUsedArea(const GrPoint& point) const;
    double getInCellUsedAreaRelax(const GrPoint& point
                                    , std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
                                    , std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls) const;
    double getFixedUsedArea(const GrEdge& edge) const;  // fixed used area of the area of an edge
    DBU getFixedLength(const GrEdge& edge) const;
    double getInCellViaNum(const GrPoint& point) const;
    double getInCellViaNumRelax(const GrPoint& point, std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls) const;
    // double getUnitViaArea(const GrPoint& point, int side) const;

    void removeWire(const GrBoxOnLayer& box);
    void removeWireRelax(const GrBoxOnLayer& box, std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls);
    void removeVia(const GrBoxOnLayer& box);
    void removeViaRelax(const GrBoxOnLayer& box, std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls);

    void useWire(const GrBoxOnLayer& box);
    void useWirePessimistic(const GrBoxOnLayer& box);
    void useWire(int layerIdx, int gridline, int cp, double usage = 1);
    void useVia(const GrBoxOnLayer& box);
    void useVia(int layerIdx, int x, int y, double usage = 1);
    void useHistWire(int layerIdx, int gridline, int cp, double usage);
};

class GrRouteGrid2D {
public:
    using UsageT = double;
    using UsageMapT = vector<vector<vector<UsageT>>>;

    void init2DMaps(const GrRouteGrid& routeGrid);
    void useWire2D(int dir, int gridline, int cp, double usage = 1);
    void removeUsage2D(int dir, int gridline, int cp, double usage = 1);
    double getCost2D(int dir, int gridline, int cp) const;
    

private:
    UsageMapT wireUsageMap2D;   // 2d, first dim show X,Y
    UsageMapT fixedMetalMap2D;  // 2d, first dim show X,Y
    UsageMapT capacityMap2D;    // 2d, first dim show X,Y
};

}  // namespace gr
