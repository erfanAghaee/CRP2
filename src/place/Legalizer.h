#pragma once

#include "db/RsynService.h"
#include "db/GeoPrimitive.h"
#include "db/Setting.h"
#include "db/Database.h"
#include "gr_db/GrDatabase.h"
// #include "single_net/SingleNetRouter.h"
// #include "db/RsynService.h"
// #include "db/Cell.h"
// #include "db/GeoPrimitive.h"
// #include "flute/flute.h"
// Magic tricks to have CPLEX behave well:
// source https://github.com/alberto-santini/cplex-example
#ifndef IL_STD
#define IL_STD
#endif
#include <cstring>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
// End magic tricks

namespace db {

struct cellWrap{
        int idx,r,s;
        double cost;
    };
struct CellOvrlp{
    int i , j;
};

struct Segment{
    int r; // row_idx
    int start,stop;
    std::vector<std::vector<int>> orders;
};//end Segment struct


class Legalizer {
public:
    Legalizer(int cell_idx) 
        : cell_idx_ (cell_idx)
    {
        int min_site_    = -1;
        int max_site_    = -1;
        int site_offset_ = -1;
        int row_offset_  = -1;
        int min_row_    = -1;
        int max_row_    = -1;
        debug_global = false;
        debug_global_all = false;
    }
    ~Legalizer(){}
    void run();
    
    std::vector<std::tuple<int,int,int>> legalizer_sols;

    bool debug_global;
    bool debug_global_all;

private:
    bool getLegalizeBox(utils::BoxT<DBU>& box, utils::BoxT<DBU>& legalize_box);
    void getMovableCellsInsideLegalizeBox(std::vector<db::Cell>& cells
                                ,utils::BoxT<DBU>& legalize_box
                                ,std::set<int>& cells_inside);
    bool isLegalizedRowsValid(std::vector<int>& legalize_rows);
    bool isLegalizedSitesValid(std::vector<int>& legalize_sites);
    int getLegalizerBlockageMatrix(std::vector<db::Cell>& cells
                                  ,std::set<int>& cells_inside_set  
                                  ,std::vector<int>& legalize_rows
                                  ,std::vector<int>& legalize_sites
                                  ,std::vector<std::vector<int>>& blockage_matrix);
    
    void updateMinMaxSites(std::vector<int>& legalize_sites);
    void updateMinMaxRows(std::vector<int>& legalize_rows);
    
    // legalizers 
    void legalizerBacktracking(  std::vector<int>& movable_cells
                               , std::vector<int>& legalize_rows
                               , std::vector<int>& legalize_sites
                               , std::vector<std::vector<int>>& blockage_matrix);

    
    void getAvailableSpaceInEachRow(std::vector<std::vector<int>>& blockage_matrix
                                    , std::vector<std::pair<int,int>>& rowIdxToAvailableSpace);


    bool isEnoughSpaceForCellInRow(int x
                                  , int r
                                  , std::vector<int>& movable_cells
                                  , std::vector<int>& rowIdxToSpaceUsedInRow
                                  , std::vector<std::pair<int,int>>& rowIdxToAvailableSpace);

    void updateSpaceUsedInRow(int x
                            , int r
                            , std::vector<int>& movable_cells
                            , std::vector<std::vector<int>>& cellRow_table
                            , std::vector<int>& rowIdxToSpaceUsedInRow
                            , bool increase);


    void collectSolutionBT(std::vector<int>& movable_cells
                                ,std::vector<int>& legalize_rows
                                ,std::vector<std::vector<int>>& cellRow_table
                                ,std::vector<std::pair<int,int>>& cellRow_sols);
    
    
    bool isSolutionBT(int placed_cells,std::vector<int>& movable_cells);

    void getAllPossibleArrangementsOfCellsInDifferentRows(std::vector<int>& movable_cells
                                                        , std::vector<int>& legalize_rows
                                                        , std::vector<std::pair<int,int>>& queue_bt
                                                        , std::vector<int>& rowIdxToSpaceUsedInRow
                                                        , std::vector<std::pair<int,int>>& rowIdxToAvailableSpace
                                                        , std::vector<std::vector<int>>& cellRow_table
                                                        , std::vector<std::vector<std::pair<int,int>>>& cellRow_sols_table);

    void legalizerILP(std::vector<int>& movable_cells
                    , std::vector<int>& legalize_rows
                    , std::vector<int>& legalize_sites
                    , std::vector<std::vector<int>>& blockage_matrix);


    void legalizerILPV2(std::vector<int>& movable_cells
                    , std::vector<int>& legalize_rows
                    , std::vector<int>& legalize_sites
                    , std::vector<std::vector<int>>& blockage_matrix);                
    
    void getLegalizerILPCostCube(std::vector<double>& weights
                        , std::vector<int>& legalize_rows
                        , std::vector<int>& legalize_sites
                        , std::vector<int>& movable_cells
                        , std::vector<std::vector<int>>& blockage_matrix
                        , std::unordered_map<int,int>& weightToCellIdxDict
                        , std::unordered_map<int,std::tuple<int,int,int>>& weigthToCellRowSite
                        , std::vector<std::vector<std::vector<int>>>& cube_cost
                        , std::vector<std::vector<std::vector<int>>>& cube_weightIdx);
    void getLegalizerILPCostCubeV2(std::vector<double>& weights
                        , std::vector<int>& legalize_rows
                        , std::vector<int>& legalize_sites
                        , std::vector<int>& movable_cells
                        , std::vector<std::vector<int>>& blockage_matrix
                        , std::vector<cellWrap>& cell_wraps
                        // , std::unordered_map<int,int>& weightToCellIdxDict
                        // , std::unordered_map<int,std::tuple<int,int,int>>& weigthToCellRowSite
                        , std::vector<std::vector<std::vector<int>>>& cube_cost
                        , std::vector<std::vector<std::vector<int>>>& cube_weightIdx);

    void getLegalizerILPCostCubeV3(std::vector<cellWrap>& weights
                    , std::vector<int>& legalize_rows
                    , std::vector<int>& legalize_sites
                    , std::vector<int>& movable_cells
                    , std::vector<std::vector<int>>& blockage_matrix
                    // , std::unordered_map<int,int>& weightToCellIdxDict
                    // , std::unordered_map<int,std::tuple<int,int,int>>& weigthToCellRowSite
                    // , std::vector<std::vector<std::vector<int>>>& cube_cost
                    // , std::vector<std::vector<std::vector<int>>>& cube_weightIdx
                    );
    void findAllPermutations(std::vector<cellWrap>& weights
                , std::vector<int>& legalize_rows
                , std::vector<int>& legalize_sites
                , std::vector<int>& movable_cells
                , std::vector<std::vector<int>>& mat_spr
                , std::vector<std::vector<int>>& blockage_matrix
                // , std::unordered_map<int,int>& weightToCellIdxDict
                // , std::unordered_map<int,std::tuple<int,int,int>>& weigthToCellRowSite
                // , std::vector<std::vector<std::vector<int>>>& cube_cost
                // , std::vector<std::vector<std::vector<int>>>& cube_weightIdx
                );

    void findAllPermutationsV2(std::vector<cellWrap>& weights
                , std::vector<int>& legalize_rows
                , std::vector<int>& legalize_sites
                , std::vector<int>& movable_cells
                , std::vector<std::vector<int>>& mat_spr
                , std::vector<std::vector<int>>& blockage_matrix
                // , std::unordered_map<int,int>& weightToCellIdxDict
                // , std::unordered_map<int,std::tuple<int,int,int>>& weigthToCellRowSite
                // , std::vector<std::vector<std::vector<int>>>& cube_cost
                // , std::vector<std::vector<std::vector<int>>>& cube_weightIdx
                );
    

    void getLegalizerILPOvrlpConflictMatrix(
                std::vector<std::vector<std::vector<int>>>& cube_cost
                , std::vector<std::vector<std::vector<int>>>& cube_weightIdx
                , std::vector<int>& movable_cells
                , std::vector<int>& legalize_rows
                , std::vector<int>& legalize_sites
                , std::vector<std::vector<int>>& conflict_matrix);
        
    void getLegalizerILPOvrlpConflictMatrixV2(
                std::vector<cellWrap>& weights
                , std::vector<std::vector<int>>& ovrlps);

    void getLegalizerILPSinglePositionConflictMatrix(
                std::vector<std::vector<std::vector<int>>>& cube_cost
                , std::vector<std::vector<std::vector<int>>>& cube_weightIdx
                , std::vector<int>& movable_cells
                , std::vector<std::vector<int>>& conflict_matrix);
    
    void getLegalizerILPSinglePositionConflictMatrixV3(
                std::vector<cellWrap>& weights
                , std::vector<int>& movable_cells
                , std::vector<std::vector<int>>& constraints);
    
    void getLegalizerSolutionILPSolver(
                  std::vector<double>& weights
                , std::vector<std::vector<int>>& conflict_matrix_ovrlp
                , std::vector<std::vector<int>>& conflict_matrix_single_position
                , std::unordered_map<int,int>& weightToCellIdxDict
                , std::unordered_map<int,std::tuple<int,int,int>>& weigthToCellRowSite
                , std::vector<int>& sol);

    void getLegalizerSolutionILPSolverV2(
                  std::vector<cellWrap>& weights
                , std::vector<int>& legalize_rows
                , std::vector<int>& legalize_sites
                , std::vector<std::vector<int>>& ovrlps
                , std::vector<std::vector<int>>& constraints
                , std::vector<int>& sol);

    void getLegalizerSolutionILPSolverSolutionPool(
                  std::vector<double>& weights
                , std::vector<std::vector<int>>& conflict_matrix_ovrlp
                , std::vector<std::vector<int>>& conflict_matrix_single_position
                , std::unordered_map<int,int>& weightToCellIdxDict
                , std::vector<std::vector<int>>& sols);

    // prints 
    void printMovableCells(std::vector<int>& movable_cells);
    void printLegalizeRows(std::vector<int>& legalize_rows);
    void printLegalizeSites(std::vector<int>& legalize_sites);
    void logLegalizerBlockageMatrix(std::vector<std::vector<int>>& blockage_matrix);
    void logCubeCost( std::vector<int>& movable_cells
                    , std::vector<std::vector<std::vector<int>>>& cube_cost
                    , std::vector<std::vector<std::vector<int>>>& cube_weightIdx);
    void logCubeWeight( std::vector<int>& movable_cells,std::vector<std::vector<std::vector<int>>>& cube_weightIdx);


    double getCandidateCost(db::Cell& cell, int row_idx, int site_idx);


    void logLegalizerBoard(std::vector<std::vector<int>>& blockage_matrix
                                , std::vector<int>& legalize_rows
                                , std::vector<int>& legalize_sites);
    void logLegalizerWeights(std::vector<cellWrap>& weights
                            , std::vector<int>& legalize_rows
                            , std::vector<int>& legalize_sites);
   void logLegalizerOvrlps(std::vector<cellWrap>& weights
                        , std::vector<int>& legalize_rows
                        , std::vector<int>& legalize_sites
                        ,std::vector<std::vector<int>>& ovrlps);
    void logLegalizerSols(std::vector<cellWrap>& weights
                                , std::vector<int>& legalize_rows
                                , std::vector<int>& legalize_sites
                                , std::vector<int>& sols);
    //variables 
    int cell_idx_;
    int min_site_;
    int max_site_;
    int min_row_;
    int max_row_;
    int site_offset_;
    int row_offset_;
    // cell_idx->row->site
    
};//end legalizer class

}//end namespace db


