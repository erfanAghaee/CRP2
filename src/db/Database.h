#pragma once

#include "RsynService.h"
#include "RouteGrid.h"
#include "Net.h"
#include "Cell.h"
#include "Setting.h"
#include "Stat.h"
#include "LookupTbs.h"

class MTStat {
public:
    vector<double> durations;
    MTStat(int numOfThreads = 0) : durations(numOfThreads, 0.0) {}
    const MTStat& operator+=(const MTStat& rhs);
    friend ostream& operator<<(ostream& os, const MTStat mtStat);
};

namespace db {

class Database : public RouteGrid, public NetList, public CellList {
public:
    std::string debug_net = "net941";
    std::vector<int> db_hpwls;

    // streams for debugging
    std::string patternRoute_stream="";
    std::string astar_stream="";
    std::string astar_stream_coarse="";

    LookupTbs lookup_tb;
    utils::BoxT<DBU> dieRegion;
    std::unordered_map<std::string,int> netsToIdx;
    std::unordered_map<std::string,int> cellsToIdx;

    void init();
    void initBlockagesRTree();
    void initRsynService();
    void clear() { RouteGrid::clear(); }

    // get girdPinAccessBoxes
    // TODO: better way to differetiate same-layer and diff-layer girdPinAccessBoxes
    void getGridPinAccessBoxes(const Net& net, vector<vector<db::GridBoxOnLayer>>& gridPinAccessBoxes) const;
    RsynService& getRsynService() {return rsynService;}
    void writeDEF(const std::string& filename);

    //using for placement
    RTree cell_rtree;
    RTree blockage_rtree;
    RTree empty_rtree;
    RTree row_rtree;
    std::vector<Rsyn::PhysicalRow> rowsRsyn;
    DBU sites_origin = 0;
    DBU sites_step = 0;
    int sites_num = 0;
    DBU rows_origin = 0;
    DBU rows_step = 0;
    int rows_num = 0;

    std::set<int> critical_cells_hist;
    std::set<int> critical_cells_hist_moved;
    std::vector<DBU> cost_hist;

    std::vector<std::vector<std::string>> critical_cellsList;

    void getRowsInBox(utils::BoxT<DBU>& box,std::vector<utils::BoxT<DBU>>& rows, DBU eps);
    void getCellsInBox(utils::BoxT<DBU>& box,std::vector<db::Cell>& cells,DBU eps);
    void initRtrees();
    void initSites();
    void initIllegalPlacementBoxs();
    void getNetCells(int net_idx,std::vector<int>& cells_idx);
    std::vector<utils::BoxT<DBU>> queryBlockageRTree(utils::BoxT<DBU>& box,DBU eps);

    void getIntersectedSites(utils::BoxT<DBU>& box,std::vector<int>& sites);
    void getIntersectedRows(utils::BoxT<DBU>& box,std::vector<int>& rows);
    int getSiteDBU(DBU x);
    DBU getDBUSite(int step);
    int getRowDBU(DBU y);
    DBU getDBURow(int step);

    DBU getRowStep() {return rows_step;}
    DBU getSiteStep() {return sites_step;}


    void logCellLocations(int iter);
    void logDie();
    // void logFixedMetals(vector<std::pair<BoxOnLayer, int>>& fixedMetalVec);
    void logFixedMetals(int iter);
    void logPatternRoute(int iter);
    void logAstar(int iter);
    void logAstarCoraseGrid(int iter);
    void logLayers();
    void addToIllegalPlacementBox(utils::BoxT<DBU>& box);

    
    

    // function for filtering the functions in the felow. 
    void initPolicy();
    
    // report strings
    std::stringstream refinePlacement_report_tracker;    
    std::stringstream cost_report_tracker;
    std::unordered_map<int,utils::BoxT<DBU>> suspected_cells_dict;
    int total_cells_moved;
    // row -> set of sites
    std::unordered_map<int,std::set<int>> illegal_placement_map;
    
    
    
    std::set<std::string> policy_set;

    // for debugging legalizer critical area 
    std::stringstream ss_legalizer;
    std::vector<utils::BoxT<DBU>> illegal_placement_boxs;
    int origin_offset_die;

    DBU libDBU;
    


private:
    RsynService rsynService;
    
    

//temporary put markPinAndObs and initMTsafe to public
public:

    // mark pin and obstacle occupancy on RouteGrid
    void markPinAndObsOccupancy();

    // init safe margin for multi-thread
    void initMTSafeMargin();
};

}  //   namespace db

extern db::Database database;

namespace std {

//  hash function for Dimension
template <>
struct hash<Dimension> {
    std::size_t operator()(const Dimension d) const { return (hash<unsigned>()(d)); }
};

//  hash function for std::tuple<typename t0, typename t1, typename t2>
template <typename t0, typename t1, typename t2>
struct hash<std::tuple<t0, t1, t2>> {
    std::size_t operator()(const std::tuple<t0, t1, t2>& t) const {
        return (hash<t0>()(std::get<0>(t)) ^ hash<t1>()(std::get<1>(t)) ^ hash<t2>()(std::get<2>(t)));
    }
};

}  // namespace std

MTStat runJobsMT(int numJobs, const std::function<void(int)>& handle);
