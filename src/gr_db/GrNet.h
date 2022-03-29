#pragma once

#include "db/Database.h"
#include "GrGeoPrimitive.h"
#include "GCell.h"
#include "GridTopo.h"

namespace gr {

class NetFeatures{
    public:
        NetFeatures(){
            rsmt_ = -1;
            hpwl_ = -1;
            // median_x_y_cell = std::make_pair(0,0);
            // median_x_y_pin = std::make_pair(0,0);
            // wl = 0;
            // congestion = 0;
            // degree_connected_nets = 0;
            // degree_connected_cells = 0;
            // is_suspected_cell = false;
            // max_edge_cost = 0;
        }
        // NetFeatures(const NetFeatures& cf){
        //     median_x_y_cell.first = cf.median_x_y_cell.first;
        //     median_x_y_cell.second = cf.median_x_y_cell.second;
        //     median_x_y_pin.first = cf.median_x_y_pin.first;
        //     median_x_y_pin.second = cf.median_x_y_pin.second;
        //     wl = cf.wl;
        //     congestion = cf.congestion;
        //     degree_connected_nets = cf.degree_connected_nets;
        //     degree_connected_cells = cf.degree_connected_cells;
        // }

        // ~NetFeatures(){}

        // double getCost(){
        //     double violation_penalty = 1000000;
        //     if(!is_suspected_cell)
        //         return wl*(1+congestion);
        //     return violation_penalty*(wl*(1+congestion));
        // }

        // std::pair<double,double> median_x_y_cell;
        // std::pair<double,double> median_x_y_pin;
        // DBU wl;
        // double congestion;
        // int degree_connected_nets;
        // int degree_connected_cells;
        // bool is_suspected_cell;
        // double max_edge_cost;
        // //net idx to net degree
        // std::vector<std::pair<int,int>> connected_nets_degrees;
        DBU rsmt_;
        DBU hpwl_;
};//end class NetFeatures

class GrNet {
public:
    db::Net& dbNet;

    bool is_critical_net;

    vector<vector<GrPoint>> pinAccessBoxes;
    std::unordered_set<GrPoint> ovlpPoints;
    GrBoxOnLayer boundingBox;

    vector<std::shared_ptr<GrSteiner>> gridTopo;
    vector<vector<std::shared_ptr<GrSteiner>>> gridTopoHist;

    vector<GrBoxOnLayer> wireRouteGuides;
    vector<GrBoxOnLayer> viaRouteGuides;
    vector<GrBoxOnLayer> patchRouteGuides;

    GrNet(int i) : dbNet(database.nets[i]) { is_critical_net = false;}
    void init(const GCellGrid& gcellGrid);
    void update(const GCellGrid& gcellGrid);

    unsigned numOfPins() const { return dbNet.numOfPins(); }
    const std::string& getName() const { return dbNet.getName(); }

    void postOrderVisitGridTopo(const std::function<void(std::shared_ptr<GrSteiner>)>& visit) const;
    void preOrderVisitGridTopo(const std::function<void(std::shared_ptr<GrSteiner>)>& visit) const;
    DBU getWirelength() const;

    DBU getPathCost() const;
    DBU getLongestWire() const;
    double getEdgeHighestCost() const;

    int getCriticalCell();
    void checkRouteHist();
    void getRouteSites();

    vector<vector<PointOnLayer>> getMergedPinAccessBoxes(const std::function<PointOnLayer(GrPoint)>& pointHash) const;

    bool needToRoute() { return !isOnePin; }
    DBU getWirelengthLast(){
        if(wl_hist.size() == 0) return 0;
        return wl_hist[wl_hist.size()-1];
    }
    void updateWirelengthHist(){
        auto wl = getWirelength();
        wl_hist.push_back(wl);
    }

    void updateMaxEdgeCost();

    void checkObidenceGrVsDr();

    double getMaxEdgeCost() {return max_edge_cost;}

    void updateNetFeatures();
    NetFeatures feature_;


    void getNetCells(std::vector<int>& cells_idx);
    void calcHPWL();


    int getHPWLPositions(std::vector<DBU> xs,std::vector<DBU> ys){
        int max_x = *std::max_element(xs.begin(), xs.end());
        int min_x = *std::min_element(xs.begin(), xs.end());
        int max_y = *std::max_element(ys.begin(), ys.end());
        int min_y = *std::min_element(ys.begin(), ys.end());

        return 2.5*std::abs(max_x-min_x) + std::abs(max_y-min_y);
    }

    DBU getHPWL(){return feature_.hpwl_;}

private:
    void initPinAccessBoxes(const GCellGrid& gcellGrid);
    void updatePinAccessBoxes(const GCellGrid& gcellGrid);

    bool isOnePin = false;
    std::vector<DBU> wl_hist;
    double max_edge_cost;
};

class GrNetlist {
public:
    vector<GrNet> nets;

    void init(const GCellGrid& gcellGrid);
    void update(const GCellGrid& gcellGrid);
};
}  // namespace gr
