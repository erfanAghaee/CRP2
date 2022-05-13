#pragma once

#include "db/Database.h"
#include "gr_db/GrDatabase.h"
#include "single_net/SingleNetRouter.h"
#include "db/RsynService.h"
#include "db/Cell.h"
#include "db/GeoPrimitive.h"
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
struct CandWeight{
    int cell_idx;
    utils::BoxT<DBU> box;// cell box
    double cost;

};//

class Placer {
public:
    Placer(CongestionMap& congMap,vector<vector<int>>& routeTableTmp,int iter_t,utils::timer& profile_time,
        std::stringstream& profile_time_str);
    void runMT(std::vector<int>& netsToRoute,int cellWidth, int cellHeight);
    void runMTISPD(std::vector<int>& netsToRoute,int cellWidth, int cellHeight);

public:
    // MT
    void initCellsCandidate();
    void findBestEntryMT();
    void findBestEntryMTV2();
    void calcCellsCostMT();
    void initCellsUpdateFeatures();
    void selectCriticalCells();
    void selectCriticalCellsNetDegree();
    void selectCriticalCellsISPD();
    void updateRoutingDataBase(std::vector<int>& netsToRoute);

    void refinePlacement(std::vector<int>& netsToRoute);
    // double getCellPositionCost(db::Cell& cell, utils::BoxT<DBU> box);

    void writeToCSV();
    void logInstances();
    void logCongestion();
    void logPlacementWeights(std::vector<CandWeight>& weights);
    // void reportCells(std::string filename);
    

    
    void ILPSolver(std::vector<double>& weights
                  ,std::unordered_map<int,std::pair<int,utils::BoxT<DBU>>>& weightIdxToCell
                  ,std::vector<std::vector<int>>& conflict_matrix
                  ,std::vector<std::vector<int>>& conflict_matrix_connected
                  ,std::vector<std::vector<int>>& legalize_matrix
                  ,std::vector<int>& sol
                  ,bool maxmimize
                  ,std::string log_name = "model.lp"
                   );
    
    // solving Maximum Weight indepndent set (MWIS) to select cells to move.
    void MWIS(std::set<std::string>& critical_cells_set);
    


    void checkGlobalRoutingResults();
    void schellingSegerationModel(std::set<std::string>& critical_cells_set);// happiness and unhappiness of each cell


    void logRefinePlacement();

    // functions for extract net features
    void calcNetsHPWL();

    void logCellLocations();


    void getOvrlpConflicts(
            std::vector<CandWeight>& weights
        ,   std::vector<std::vector<int>>& ovrlps
    );
    void ILPSolverV2(
                  std::vector<CandWeight>& weights
                , std::vector<std::vector<int>>& ovrlps
                , std::vector<std::vector<int>>& constraints
                , std::vector<int>& sol);

    void logCellsFeature();
    void logCellsCritical();
    void logCellsCandidates();
    void logCellsMoved();
    

    // // cells index in database to new box for movement
    std::vector<std::pair<int,utils::BoxT<DBU>>> cellsToMove;
    CongestionMap& congMap;
    std::vector<db::Cell> cells;
    vector<vector<int>>& routeTable;
    std::set<int> critical_cells;
    int iter;
    utils::timer& profile_time;
    std::stringstream& profile_time_str;
    bool debug_global;
    std::set<std::string> critical_cells_set;


    
    
    
};


