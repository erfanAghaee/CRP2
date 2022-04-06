#pragma once

#include "global.h"

namespace db {

BETTER_ENUM(VerboseLevelT, int, LOW = 0, MIDDLE = 1, HIGH = 2);

// global setting
class Setting {
public:
    // basic
    std::string directory="./";
    std::string lookuptbs_dir="./..";
    std::string outputFile;
    std::string outputFileDef;
    bool refinePlacement;
    bool postProcessing;
    std::string name;
    std::string benchmarkName;
    std::string costFunction;
    int numThreads = 1;  // 0 for simple scheduling
    int numRefinePlacement = 1;
    int numGlobalRouting = 1;
    int legalizationSolutionPoolIntensity = 0;
    int legalizationWindowSizeRow = 1;
    int legalizationWindowSizeSite = 1;
    int legalizationMaxNumCellsToLegalize = 1;
    int legalizationMaxNumSolutionPool = 1;
    bool legalizationActivateSolutionPool = false;
    float percCriticalCells = 0.001;
    double legalizerOptimzerTimeLimit = 1;
    double patchPinRegions= false;
    double patchLongSegments= false;
    double patchVioCells= false;
    int gcellSize=3000;
    int tat = std::numeric_limits<int>::max();

    bool hasInputGuide = false;
    std::string inputGuide = "";

    // multi_net
    VerboseLevelT multiNetVerbose = VerboseLevelT::MIDDLE;
    bool multiNetScheduleSortAll = true;
    bool multiNetScheduleSort = true;
    bool multiNetScheduleReverse = true;
    int multiNetSelectViaTypesIter = 3;
    int rrrIterLimit = 4;
    bool rrrWriteEachIter = false;
    double rrrInitVioCostDiscount = 0.1;
    double rrrFadeCoeff = 0.01;  // should be <= 0.5 to make sure fade/(1-fade) <= 1

    int edgeShiftingIter = 2;

    // single_net
    VerboseLevelT singleNetVerbose = VerboseLevelT::MIDDLE;
    bool fixOpenBySST = true;
    double unitSqrtViaUsage = 1.5;   // # of tracks that unit sqrt(via) uses
    double initLogisticSlope = 1.0;  // the slope of the logistic func in cost func

    // db
    VerboseLevelT dbVerbose = VerboseLevelT::MIDDLE;
    int maxNumWarnForEachRouteStatus = 5;
    bool dbWriteDebugFile = false;
    double dbInitHistUsageForPinAccess = 0.1;

    //  Metric weights of ISPD 2018 Contest
    //  Wirelength unit is M2 pitch
    //  1. basic objective
    static constexpr double weightWirelength = 0.5;
    static constexpr int weightViaNum = 4;
    //  2. routing preference
    static constexpr int weightOutOfGuideWirelength = 1;
    static constexpr int weightOutOfGuideViaNum = 1;
    static constexpr double weightOffTrackWirelength = 0.5;
    static constexpr int weightOffTrackViaNum = 1;
    static constexpr int weightWrongWayWirelength = 1;
    static constexpr int wirelenCostWeight = 4;

    //  3. violation
    //  normalized by square of M2 pitch
    static constexpr int weightShortArea = 500;
    //  including wire spacing, eol spacing, cut spacing
    static constexpr int weightSpaceVioNum = 500;
    static constexpr int weightMinAreaVioNum = 500;

    // 4. debug 
    bool debugLegalizerExportILPModel = false;
    bool debug=false;
    std::string runtimeImprvCRPPolicies="";


    bool patchPinRegionsMode=false;

    std::string filter_nets_name="";
    std::string filter_nets_NumPin="";
    std::string critical_cells_set="";


    bool patternRouteMemorization=false;
    bool featureExtractionParallel=false;
    bool findAllPermutations=false;
    bool RLOutput=false;
    bool moveToRemoveViol=false;


    void makeItSilent();
    void adapt();
};

extern Setting setting;

}  //   namespace db
