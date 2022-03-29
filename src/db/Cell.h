#pragma once

#include "RsynService.h"
#include "GeoPrimitive.h"
#include "Setting.h"

#ifndef IL_STD
#define IL_STD
#endif
#include <cstring>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN


// 

namespace db {

// extern "C" {
// Tree flute(int d, DTYPE x[], DTYPE y[], int acc);
// }

class CellFeatures{
    public:
        CellFeatures(){
            median_x_y_cell = std::make_pair(0,0);
            median_x_y_pin = std::make_pair(0,0);
            wl = 0;
            congestion = 0;
            degree_connected_nets = 0;
            degree_connected_cells = 0;
            is_suspected_cell = false;
            max_edge_cost = 0;
        }
        CellFeatures(const CellFeatures& cf){
            median_x_y_cell.first = cf.median_x_y_cell.first;
            median_x_y_cell.second = cf.median_x_y_cell.second;
            median_x_y_pin.first = cf.median_x_y_pin.first;
            median_x_y_pin.second = cf.median_x_y_pin.second;
            wl = cf.wl;
            congestion = cf.congestion;
            degree_connected_nets = cf.degree_connected_nets;
            degree_connected_cells = cf.degree_connected_cells;
        }

        ~CellFeatures(){}

        double getCost(){
            double violation_penalty = 1000000;
            if(!is_suspected_cell)
                return wl*(1+congestion);
            return violation_penalty*(wl*(1+congestion));
        }

        std::pair<double,double> median_x_y_cell;
        std::pair<double,double> median_x_y_pin;
        DBU wl;
        double congestion;
        int degree_connected_nets;
        int degree_connected_cells;
        bool is_suspected_cell;
        double max_edge_cost;
        //net idx to net degree
        std::vector<std::pair<int,int>> connected_nets_degrees;
};//end class CellFeatures



class CostEq{
public:
    CostEq()
    {
        wl_cost = 0;
        cong_cost = 0;
        via_cost = 0;
        via_abs_cost = 0;
        wl_abs_cost = 0;
        path_cost = 0;
        alignment_cost = 0;
        hpwl_cost = 0;
        hybrid_cost = 0;
        has_vio = false;
    };
    void updateAlignmentCost(DBU x) {alignment_cost += x;}
    void updateWireLengthCost(DBU x) {wl_cost += x;}
    void updateViaCost(DBU x) {via_cost += x;}
    void updateCongestionCost(DBU x) {cong_cost += x;}
    void updateViaAbsCost(DBU x) {via_abs_cost += x;}
    void updateWireLengthAbsCost(DBU x) {wl_abs_cost += x;}
    void updatePathCost(DBU x) {path_cost += x;}
    void updateHPWLCost(DBU x) {hpwl_cost += x;}
    void updateHybridCost(DBU x) {hybrid_cost += x;}
    void updateHasVio(bool vio) {
        if(vio==true){
            has_vio = true;
        }
    }//end updateHasVio

    DBU getCost(){
        // cost equation can be defined here
        // return alignment_cost;
        if(db::setting.costFunction == "alignment_cost"){
            return alignment_cost;
        }else if (db::setting.costFunction == "wl_cost"){
            return wl_cost;
        }else if (db::setting.costFunction == "via_cost"){
            return via_cost;
        }else if (db::setting.costFunction == "cong_cost"){
            return cong_cost;
        }else if (db::setting.costFunction == "via_abs_cost"){
            return via_abs_cost;
        }else if (db::setting.costFunction == "wl_abs_cost"){
            return wl_abs_cost;
        }else if (db::setting.costFunction == "path_cost"){
            return path_cost;
        }else if (db::setting.costFunction == "hpwl_cost"){
            return hpwl_cost;
        }else if (db::setting.costFunction == "path_hpwl_cost"){
            return hpwl_cost+path_cost;
        }else if (db::setting.costFunction == "hybrid_cost"){
            return hybrid_cost;
        }
        return path_cost;
    }

    DBU wl_cost;
    DBU cong_cost;
    DBU via_cost;
    DBU via_abs_cost;
    DBU wl_abs_cost;
    DBU path_cost;
    DBU alignment_cost;
    DBU hpwl_cost;
    DBU hybrid_cost;
    bool has_vio;
};


class CellBase {
public:
    ~CellBase();

    int idx;
    Rsyn::Instance rsynInstance;
    const std::string& getName() const { return rsynInstance.getName(); }

    // pins
    // vector<Rsyn::Pin> rsynPins;
    // vector<vector<BoxOnLayer>> pinAccessBoxes;  // (pinIdx, accessBoxIdx) -> BoxOnLayer
    // unsigned numOfPins() const noexcept { return pinAccessBoxes.size(); }
    // BoxOnLayer getMaxAccessBox(int pinIdx) const;

    // // route guides
    // vector<BoxOnLayer> routeGuides;
    // vector<GridBoxOnLayer> gridRouteGuides;

    void print(ostream& os = std::cout) const;
};


class Cell : public CellBase {
public:
    Cell(int i, Rsyn::Instance cell, RsynService& rsynService,bool debug_t = false);
    bool debug;
    // more cell guide information
    // the cell cost
    double cost;
    bool isFixed;
    bool is_critical_cell;
    std::pair<double,double> median_x_y_cell;
    std::pair<double,double> median_x_y_pin;
    

    // placement candidate is the pair of box and cost
    std::vector<std::pair<utils::BoxT<DBU>,double>> placement_candidates;

    void calcCostPlacementCandidates();
    void initCandidates();
    void initConnectedCells();
    void addOvrlpCandidates(utils::BoxT<DBU> box);
    // median of all cells
    void getDistOptPoint();
    // median of net boxs
    void getDistOptPointV2();

    DBU getWirelength();
    DBU getWirelengthLast();

    void updateCellFeatures();
    CellFeatures feature;

    // cell index to box position
    vector<std::pair<int,utils::BoxT<DBU>>> ovrlp_cells_;
    void addMedianCellCandidates();
    void addLegalizedCellCandidates();
    void addILPLegalizedCellCandidates();
    void addBacktrackingLegalizedCellCandidates();
    void addMedianCellCandidatesSiteBased();
    void addPinAlignmentCellCandidates();
    void addManualCellCandidates();
    void addRandomNetBoxCellCandidates();
    void addEqualWidthNetBoxCellCandidates();
    void addWireDependentCellCandidates();
    void addGCellBasedCandidates();
    void legalize(utils::BoxT<DBU>& init_box,utils::BoxT<DBU>& new_box);
    void abacuslegalizer(utils::BoxT<DBU>& init_box);
    void backtrackingLegalizer( std::vector<int>& legalize_rows
                            , std::vector<int>& legalize_sites
                            , std::vector<std::pair<int,int>>& movable_cells_width
                            , std::vector<std::vector<int>>& blockage_matrix);
    void getNetsBox(utils::BoxT<DBU>& box);
    void updateMedianCell();
    void updateMedianCellV2();
    void updateMedianPin();

    std::set<int> connectd_cells;
    std::unordered_map<std::string,utils::BoxT<DBU>> pin_bd_dict;
// protected:

    bool isRedundantCandidate(utils::BoxT<DBU> box);

    bool isSameCellVsOvrlpCells(std::vector<db::Cell>& cells);
    bool isValidPostion(utils::BoxT<DBU>& box);
    bool getCellsInBoxNeigh(utils::BoxT<DBU>& cur_box,std::vector<db::Cell>& cells,
                            utils::BoxT<DBU>& new_box,DBU eps,std::string mode);
    int getCellSiteWidth();                            
    void collectCandidatePositions(utils::BoxT<DBU>& placement_candidate);
    void collectLegalizedCandidatePositions(utils::BoxT<DBU>& placement_candidate);
    
    void normalizeCostEq();

    std::pair<double,double> getMedianCell() {return median_x_y_cell;}
    std::pair<double,double> getMedianPin() {return median_x_y_pin;}
    std::pair<double,double> getMedianBoxs(std::vector<utils::BoxT<DBU>>& boxs);
    
    std::vector<utils::BoxT<DBU>> getRowsNear(double x, double y);
    std::pair<double,double> getRandom();
    std::pair<double,double> getRandom(utils::BoxT<DBU>& box);
    void getCellNetsIdx(std::vector<int>& nets_idx);
    void calcCostPlacementCandidate(utils::BoxT<DBU>& box);
    int getHPWLPositions(std::vector<DBU> xs,std::vector<DBU> ys);
    void getNetCells(int net_idx,std::vector<int>& cells_idx);

    void initPinAccessBoxes(Rsyn::Pin rsynPin,
                             RsynService& rsynService,
                             vector<BoxOnLayer>& accessBoxes,
                             utils::BoxT<DBU>& new_position,
                             const DBU libDBU);

    void getPinAccessBoxes(Rsyn::PhysicalPort phPort, vector<BoxOnLayer>& accessBoxes);
    
    void getPinAccessBoxes(Rsyn::PhysicalLibraryPin phLibPin,
                            Rsyn::PhysicalCell phCell,
                            vector<BoxOnLayer>& accessBoxes,
                            utils::BoxT<DBU>& new_position,
                            const DBUxy& origin);
    // Get intersect over Union X
    DBU getPinIoU(vector<vector<BoxOnLayer>>& pinAccessBoxes,int tgPin_idx);
    DBU getPinTotalArea(vector<vector<BoxOnLayer>>& pinAccessBoxes,int tgPin_idx);
    
    Rsyn::PhysicalOrientation getMoveOrientation(utils::BoxT<DBU>& box);
    void getConnectedCells(std::vector<int>& cells_idx);
    void updatePinBoundingBox();
    void addAlignmentCandidates(std::vector<utils::BoxT<DBU>>& placement_candidates_alignment);
    // any given box will be validated according to site
    void palceInSite(utils::BoxT<DBU>& box);
    void getIntesectedSites(utils::BoxT<DBU>& box,std::vector<int>& sites);
    int getSiteDBU(DBU x);
    DBU getDBUSite(int step);
    void getIntesectedRows(utils::BoxT<DBU>& box,std::vector<int>& rows);
    int getRowDBU(DBU y);
    DBU getDBURow(int step);

    int updatePinAccessBoxes(int net_idx,utils::BoxT<DBU>& box,vector<vector<BoxOnLayer>>& pinAccessBoxes);
    void updateCellBoxes(int net_idx, utils::BoxT<DBU>& box, 
          std::vector<utils::BoxT<DBU>>& cell_boxs);

    void calcAlignmentCost(vector<vector<BoxOnLayer>>& pinAccessBoxes,int tgPin_idx,CostEq& cost_eq);
    void calcAlignmentSiteBasedCost(vector<vector<BoxOnLayer>>& pinAccessBoxes,int tgPin_idx,CostEq& cost_eq);
    void calcFluteRouting3DCost(int net_idx
                                    ,vector<vector<BoxOnLayer>>& pinAccessBoxes
                                    ,std::vector<utils::BoxT<DBU>>& cell_boxs
                                    ,CostEq& cost_eq);
    void calcHPWLCost(std::vector<utils::BoxT<DBU>>& cell_boxs
                            ,CostEq& cost_eq);
    std::vector<CostEq> cost_eqs;

    void printPlacementCandidates();
    void ILPSolverLegalizer(std::vector<double>& weights
        ,std::unordered_map<int,std::pair<int,utils::BoxT<DBU>>>& weightIdxToCell
        ,std::vector<std::vector<int>>& conflict_matrix
        ,std::vector<std::vector<int>>& conflict_matrix_connected
        ,std::vector<std::vector<int>>& legalize_matrix
        ,std::vector<int>& sol
        ,bool maxmimize
        ,std::string log_name 
         );

    utils::BoxT<DBU> getCellBox();
    int getCellRow();
    int getCellSite();

    // Legalizer Functions
    bool getLegalizeBox(utils::BoxT<DBU>& new_box, utils::BoxT<DBU>& legalize_box);
    void getCellsInsideLegalizeBox(std::vector<db::Cell>& cells
                                  ,utils::BoxT<DBU>& legalize_box
                                  ,std::set<int>& cells_inside);
    int getLegalizerBlockageMatrix(std::vector<db::Cell>& cells
                                  ,std::set<int>& cells_inside_set  
                                  ,std::vector<int>& legalize_rows
                                  ,std::vector<int>& legalize_sites
                                  ,std::vector<std::vector<int>>& blockage_matrix);
    void logLegalizerBlockageMatrix(std::vector<std::vector<int>>& blockage_matrix);
    int getLegalizerInitCellPositions(utils::BoxT<DBU>& new_box
                            ,std::set<int>& cells_inside_set
                            ,std::vector<int>& movable_cells
                            ,std::unordered_map<int,std::pair<int,int>>& cells_init_position);
    void getLegalizerCostCube(std::vector<double>& weights
                            , std::vector<int>& legalize_rows
                            , std::vector<int>& legalize_sites
                            , std::vector<int>& movable_cells
                            , std::unordered_map<int,std::pair<int,int>>& cells_init_position
                            , std::vector<std::vector<int>>& blockage_matrix
                            , std::vector<int>& cells_width_site
                            , std::unordered_map<int,int>& weightToCellIdxDict
                            , std::unordered_map<int,std::tuple<int,int,int>>& weigthToCellRowSite
                            , std::vector<std::vector<std::vector<int>>>& cube_cost
                            , std::vector<std::vector<std::vector<int>>>& cube_weightIdx);

    void logCubeWeight( std::vector<int>& movable_cells,std::vector<std::vector<std::vector<int>>>& cube_weightIdx);
    void logCubeCost( std::vector<int>& movable_cells
                    , std::vector<std::vector<std::vector<int>>>& cube_cost
                    , std::vector<std::vector<std::vector<int>>>& cube_weightIdx);
    
    void getLegalizerOvrlpConflictMatrix(
                    std::vector<std::vector<std::vector<int>>>& cube_cost
                    , std::vector<std::vector<std::vector<int>>>& cube_weightIdx
                    , std::vector<int>& cells_width_site
                    , std::vector<int>& legalize_rows
                    , std::vector<int>& legalize_sites
                    , std::vector<std::vector<int>>& conflict_matrix);
    
    void getLegalizerSinglePositionConflictMatrix(std::vector<std::vector<std::vector<int>>>& cube_cost
                                                , std::vector<std::vector<std::vector<int>>>& cube_weightIdx
                                                , std::vector<int>& movable_cells
                                                , std::vector<std::vector<int>>& conflict_matrix);


    void logLegalizerConflictMatrix(std::vector<std::vector<int>>& conflict_matrix);


    void getLegalizerSolutionILPSolver(std::vector<double>& weights
                            , std::vector<std::vector<int>>& conflict_matrix_ovrlp
                            , std::vector<std::vector<int>>& conflict_matrix_single_position
                            , std::unordered_map<int,int>& weightToCellIdxDict
                            , std::vector<int>& sol);

    bool isMultiRowCell();

};

class CellList {
public:
    vector<Cell> cells;
    

    void init(RsynService& rsynService);
};

}  // namespace db