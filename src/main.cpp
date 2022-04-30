#include "db/Database.h"
#include "global.h"
#include "multi_net/Router.h"
#include "gr_db/GrDatabase.h"

#include <cmath>



double LARGE_NUM = 100000000;

void generateLookUpTable(){
    int len = 3;
    // permutation string
    std::stringstream perm_ss;

    // for(int len = 0; len < total_len; len++){
        std::vector<int> orders;
        std::set<std::string> orders_set;
        std::vector<std::vector<int>> orders_2d;

        
        int expected_results = 0;
        for(int i = 0; i <= len ; i++){
            expected_results += std::tgamma(len+1)/std::tgamma(len-i+1);

        }

        for(int i = 0; i < len; i++){
            orders.push_back(i);
        }

        for(int n =0; n <= orders.size(); n++){
            do{

                //Display the current permutation
                std::vector<int> order_tmp;
                std::string txt2 = "";
                for(int i=0;i<n;i++) {
                    txt2 = txt2 + std::to_string(orders[i]) + " ";
                }
            
                orders_set.insert(txt2);
                //     cout << orders[i] << " ";
                // cout << endl;
            }while(std::next_permutation(orders.begin(), orders.end()));
        }//end for 


        
            
        for(auto txt : orders_set){
            // log() << txt << std::endl;
            perm_ss << txt << std::endl;
        }

        // log() << "len: " << len << std::endl;
        // log() << "size: " << orders_set.size() << std::endl;
        // log() << "expected_results: " << expected_results << std::endl;

        
    // }
    
    std::ofstream perm_file("./../lookupTable/perm3.txt");
    perm_file << perm_ss.str();
    perm_file.close();
    
        
    
}//end generateLookUpTable


// -----------------------------------------------------------------------------

void signalHandler(int signum) {
    std::cout << "Signal (" << signum << ") received. Exiting...\n";

    // cleanup and close up stuff here

    std::exit(signum);
}

// -----------------------------------------------------------------------------

void runISPD18Flow(const boost::program_options::variables_map& vm) {
    // db::setting.makeItSilent();

    Rsyn::Session session;

    // Parse options
    // required
    std::string lefFile = vm.at("lef").as<std::string>();
    std::string defFile = vm.at("def").as<std::string>();
    db::setting.numThreads = vm.at("threads").as<int>();
    db::setting.legalizationWindowSizeRow = vm.at("legalizationWindowSizeRow").as<int>();
    db::setting.legalizationWindowSizeSite = vm.at("legalizationWindowSizeSite").as<int>();
    db::setting.legalizationMaxNumCellsToLegalize = vm.at("legalizationMaxNumCellsToLegalize").as<int>();
    db::setting.legalizationMaxNumSolutionPool = vm.at("legalizationMaxNumSolutionPool").as<int>();
    db::setting.legalizationActivateSolutionPool = vm.at("legalizationActivateSolutionPool").as<bool>();
    db::setting.legalizationSolutionPoolIntensity = vm.at("legalizationSolutionPoolIntensity").as<int>();
    db::setting.outputFile = vm.at("output").as<std::string>();
    db::setting.outputFileDef = vm.at("outputDef").as<std::string>();
    db::setting.refinePlacement = vm.at("refinePlacement").as<bool>();
    db::setting.costFunction = vm.at("costFunction").as<std::string>();
    db::setting.numGlobalRouting = vm.at("numGlobalRouting").as<int>();
    db::setting.numRefinePlacement = vm.at("numRefinePlacement").as<int>();
    db::setting.percCriticalCells = vm.at("percCriticalCells").as<float>();
    db::setting.postProcessing = vm.at("postProcessing").as<bool>();
    db::setting.benchmarkName = vm.at("benchmarkName").as<std::string>();
    db::setting.inputGuide = vm.at("inputGuide").as<std::string>();
    db::setting.hasInputGuide = vm.at("hasInputGuide").as<bool>();
    db::setting.patchPinRegions = vm.at("patchPinRegions").as<bool>();
    db::setting.patchLongSegments = vm.at("patchLongSegments").as<bool>();
    db::setting.patchVioCells = vm.at("patchVioCells").as<bool>();
    db::setting.legalizerOptimzerTimeLimit = vm.at("legalizerOptimzerTimeLimit").as<double>();
    db::setting.gcellSize = vm.at("gcellSize").as<int>();
    db::setting.debug = vm.at("debug").as<bool>();
    
    
    
    
    // optional
    if (vm.count("tat")) {
        db::setting.tat = vm.at("tat").as<int>();
    }
    // multi_net
    if (vm.count("multiNetVerbose")) {
        db::setting.multiNetVerbose =
            db::VerboseLevelT::_from_string(vm.at("multiNetVerbose").as<std::string>().c_str());
    }
    if (vm.count("multiNetScheduleSortAll")) {
        db::setting.multiNetScheduleSortAll = vm.at("multiNetScheduleSortAll").as<bool>();
    }
    if (vm.count("multiNetScheduleReverse")) {
        db::setting.multiNetScheduleReverse = vm.at("multiNetScheduleReverse").as<bool>();
    }
    if (vm.count("multiNetScheduleSort")) {
        db::setting.multiNetScheduleSort = vm.at("multiNetScheduleSort").as<bool>();
    }
    if (vm.count("rrrIters")) {
        db::setting.rrrIterLimit = vm.at("rrrIters").as<int>();
    }
    if (vm.count("rrrWriteEachIter")) {
        db::setting.rrrWriteEachIter = vm.at("rrrWriteEachIter").as<bool>();
    }
    if (vm.count("rrrInitVioCostDiscount")) {
        db::setting.rrrInitVioCostDiscount = vm.at("rrrInitVioCostDiscount").as<double>();
    }
    if (vm.count("rrrFadeCoeff")) {
        db::setting.rrrFadeCoeff = vm.at("rrrFadeCoeff").as<double>();
    }
    if (vm.count("edgeShiftingIter")) {
        db::setting.edgeShiftingIter = vm.at("edgeShiftingIter").as<int>();
    }
    // single_net
    if (vm.count("fixOpenBySST")) {
        db::setting.fixOpenBySST = vm.at("fixOpenBySST").as<bool>();
    }
    // db
    if (vm.count("dbVerbose")) {
        db::setting.dbVerbose = db::VerboseLevelT::_from_string(vm.at("dbVerbose").as<std::string>().c_str());
    }
    if (vm.count("dbInitHistUsageForPinAccess")) {
        db::setting.dbInitHistUsageForPinAccess = vm.at("dbInitHistUsageForPinAccess").as<double>();
    }
    if (vm.count("debugLegalizerExportILPModel")) {
        db::setting.debugLegalizerExportILPModel = vm.at("debugLegalizerExportILPModel").as<bool>();
    }
    if (vm.count("patchPinRegionsMode")) {
        db::setting.patchPinRegionsMode = vm.at("patchPinRegionsMode").as<bool>();
    }
    if (vm.count("filter_nets_name")) {
        db::setting.filter_nets_name = vm.at("filter_nets_name").as<std::string>();
    }
    if (vm.count("filter_nets_NumPin")) {
        db::setting.filter_nets_NumPin = vm.at("filter_nets_NumPin").as<std::string>();
    }
    if (vm.count("runtimeImprvCRPPolicies")) {
        db::setting.runtimeImprvCRPPolicies = vm.at("runtimeImprvCRPPolicies").as<std::string>();
    }
    if (vm.count("patternRouteMemorization")) {
        db::setting.patternRouteMemorization = vm.at("patternRouteMemorization").as<bool>();
    }
    if (vm.count("featureExtractionParallel")) {
        db::setting.featureExtractionParallel = vm.at("featureExtractionParallel").as<bool>();
    }
    if (vm.count("critical_cells_set")) {
        db::setting.critical_cells_set = vm.at("critical_cells_set").as<std::string>();
    }
    if (vm.count("findAllPermutations")) {
        db::setting.findAllPermutations = vm.at("findAllPermutations").as<bool>();
    }
    if (vm.count("RLOutput")) {
        db::setting.RLOutput = vm.at("RLOutput").as<bool>();
    }

    if (vm.count("moveToRemoveViol")) {
        db::setting.moveToRemoveViol = vm.at("moveToRemoveViol").as<bool>();
    }
    if (vm.count("directory")) {
        db::setting.directory = vm.at("directory").as<std::string>();
    }
    if (vm.count("lookuptbs_dir")) {
        db::setting.lookuptbs_dir = vm.at("lookuptbs_dir").as<std::string>();
    }
    if (vm.count("rrrRouters")) {
        db::setting.rrrRouters = vm.at("rrrRouters").as<std::string>();
    }
    if (vm.count("rrrRoutersApply")) {
        db::setting.rrrRoutersApply = vm.at("rrrRoutersApply").as<std::string>();
    }
    if (vm.count("logAll")) {
        db::setting.logAll = vm.at("logAll").as<bool>();
    }
    



    
    
    
    


    utils::timer total_runtime;

    vector<std ::string> strs;
    boost::split(strs, db::setting.outputFile, boost::is_any_of("."));
    db::setting.name = strs[0];

    // Read benchmarks
    Rsyn::ISPD2018Reader reader;
    const Rsyn::Json params = {
        {"lefFile", lefFile},
        {"defFile", defFile},
        {"inputGuide", db::setting.inputGuide}
    };
    log() << std::endl;
    if (db::setting.dbVerbose >= +db::VerboseLevelT::HIGH) {
        log() << "################################################################" << std::endl;
        log() << "Start reading benchmarks" << std::endl;
    }
    reader.load(params);
    if (db::setting.dbVerbose >= +db::VerboseLevelT::HIGH) {
        log() << "Finish reading benchmarks" << std::endl;
        log() << "MEM: cur=" << utils::mem_use::get_current() << "MB, peak=" << utils::mem_use::get_peak() << "MB"
              << std::endl;
        log() << std::endl;
    }

    // generateLookUpTable();
    // return;

    // Route
    database.initRsynService();
    database.init();
    db::setting.adapt();
    grDatabase.init();

    log() << "finish reading benchmark" << std::endl;

    Router router;
    router.run();
    // router.runISPD();
    // router.printCongMap();

    // Router router1;
    // router1.run();
    

    // Router router;
    // router.run();

    grDatabase.writeGuides(db::setting.outputFile);
    // grDatabase.writeGuidesCSV(db::setting.outputFile);
    // grDatabase.reportGR(db::setting.outputFile);
    // grDatabase.reportCells(db::setting.outputFile);

    if(db::setting.RLOutput)
        grDatabase.reportRL(db::setting.outputFile);

    database.writeDEF(db::setting.outputFileDef);

    database.clear();
    grDatabase.clear();

    

    log() << "Total Time: "<< total_runtime.getTimer()
          <<", MEM: cur=" << utils::mem_use::get_current() << "MB, peak=" << utils::mem_use::get_peak() << "MB"
          << std::endl;
    log() << std::endl;
}

// -----------------------------------------------------------------------------

int main(int argc, char* argv[]) {

    signal(SIGINT, signalHandler);
    signal(SIGTERM, signalHandler);

    Rsyn::Session::init();

    printlog("------------------------------------------------------------------------------");
    printlog("                     Global Routing & RefinePlacement                      ");
    printlog("        Developer: Erfan Aghaeekiasaraee                                     ");
    printlog("        Supervisors: Dr. Laleh Behjat, Dr. David Westwick                                         ");
    printlog("        Affiliation: University of Calgary                     ");
    printlog("        Note: The source code is developed on CUGR                     ");
    printlog("------------------------------------------------------------------------------");

    std::cout << std::boolalpha;  // set std::boolalpha to std::cout

    try {
        using namespace boost::program_options;
        options_description desc{"Options"};
        // clang-format off
        desc.add_options()
                ("help", "Help message.")
                ("lef", value<std::string>()->required(), "Input .lef file")
                ("def", value<std::string>()->required(), "Input .def file.")
                ("threads", value<int>()->required(), "# of threads")
                ("output", value<std::string>()->required(), "Output file name")
                ("outputDef", value<std::string>()->required(), "Output def file name")
                ("refinePlacement", value<bool>()->required(), "refinePlacement is required")
                ("costFunction", value<std::string>()->required(), "costFunction is required")
                ("numRefinePlacement", value<int>()->required(), "NumRefinePlacement iterations is required")
                ("legalizationWindowSizeRow", value<int>()->required(), "legalizationWindowSizeRow is required")
                ("legalizationWindowSizeSite", value<int>()->required(), "legalizationWindowSizeSite is required")
                ("legalizationMaxNumCellsToLegalize", value<int>()->required(), "legalizationMaxNumCellsToLegalize is required")
                ("legalizationMaxNumSolutionPool", value<int>()->required(), "legalizationMaxNumSolutionPool is required")
                ("legalizationActivateSolutionPool", value<bool>()->required(), "legalizationActivateSolutionPool is required")
                ("legalizationSolutionPoolIntensity", value<int>()->required(), "legalizationSolutionPoolIntensity is required")
                ("numGlobalRouting", value<int>()->required(), "NumGlobalRouting iterations is required")
                ("percCriticalCells", value<float>()->required(), "percCriticalCells (Percentage of critical cells in each iteration) is required")
                ("postProcessing", value<bool>()->required(), "postProcessing is required")
                ("benchmarkName", value<std::string>()->required(), "benchmarkName is required")
                ("inputGuide", value<std::string>()->required(),"inputGuide is required")
                ("hasInputGuide", value<bool>()->required(),"hasInputGuide is required")
                ("patchPinRegions", value<bool>()->required(),"patchPinRegions is required")
                ("patchLongSegments", value<bool>()->required(),"patchLongSegments is required")
                ("patchVioCells", value<bool>()->required(),"patchVioCells is required")
                ("legalizerOptimzerTimeLimit", value<double>()->required(),"legalizerOptimzerTimeLimit is required")
                ("gcellSize", value<int>()->required(),"gcellSize is required")
                ("debug", value<bool>()->required(),"debug is required")
                // optional
                ("tat", value<int>(), "Runtime limit (sec)")
                ("multiNetVerbose", value<std::string>())
                ("multiNetScheduleSortAll", value<bool>())
                ("multiNetScheduleReverse", value<bool>())
                ("multiNetScheduleSort", value<bool>())
                ("rrrIters", value<int>())
                ("rrrWriteEachIter", value<bool>())
                ("rrrInitVioCostDiscount", value<double>())
                ("rrrFadeCoeff", value<double>())
                ("edgeShiftingIter", value<int>())
                ("fixOpenBySST", value<bool>())
                ("dbVerbose", value<std::string>())
                ("dbInitHistUsageForPinAccess", value<double>())
                ("debugLegalizerExportILPModel", value<bool>())
                ("patchPinRegionsMode", value<bool>())
                ("filter_nets_name", value<std::string>())
                ("filter_nets_NumPin", value<std::string>())
                ("runtimeImprvCRPPolicies", value<std::string>())
                ("patternRouteMemorization", value<bool>())
                ("featureExtractionParallel", value<bool>())
                ("critical_cells_set", value<std::string>())
                ("findAllPermutations", value<bool>())     
                ("RLOutput", value<bool>())  
                ("moveToRemoveViol", value<bool>())         
                ("directory", value<std::string>())    
                ("lookuptbs_dir", value<std::string>())    
                ("rrrRouters", value<std::string>())    
                ("rrrRoutersApply", value<std::string>())    
                ("logAll", value<bool>())
                ;
        // clang-format on
        variables_map vm;
        store(command_line_parser(argc, argv)
                  .options(desc)
                  .style(command_line_style::unix_style | command_line_style::allow_long_disguise)
                  .run(),
              vm);
        notify(vm);
        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return 0;
        }
        for (const auto& option : desc.options()) {
            if (vm.count(option->long_name())) {
                std::string name = option->description().empty() ? option->long_name() : option->description();
                log() << std::left << std::setw(18) << name << ": ";
                const auto& value = vm.at(option->long_name()).value();
                if (auto v = boost::any_cast<double>(&value)) {
                    std::cout << *v;
                } else if (auto v = boost::any_cast<int>(&value)) {
                    std::cout << *v;
                } else if (auto v = boost::any_cast<std::string>(&value)) {
                    std::cout << *v;
                } else if (auto v = boost::any_cast<bool>(&value)) {
                    std::cout << *v;
                } else {
                    std::cout << "unresolved type";
                }
                std::cout << std::endl;
            }
        }
        runISPD18Flow(vm);
    } catch (const boost::program_options::error& e) {
        printlog(e.what());
    }

    printlog("---------------------------------------------------------------------------");
    printlog("                               Terminated...                               ");
    printlog("---------------------------------------------------------------------------");

    return 0;
}


