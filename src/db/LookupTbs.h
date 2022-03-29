#pragma once
#include "global.h"

namespace db{
class LookupTbs{
public:
    LookupTbs()
    {}
    ~LookupTbs(){}

    void run(std::string& bench_name){
        log() << "init lookup tables..."<< std::endl;

        readPermutations();
        readIellgalPlacementBoxes(bench_name);
        readOriginOffsetDie(bench_name);
        
        

    }//end runfunction

    void readPermutations(){
        for(int ii = 1; ii <= 3; ii++ ){
            std::string file_name = "./../lookupTable/perm"+\
                                     std::to_string(ii)+".txt";
            std::ifstream file(file_name);
            std::vector<std::vector<int>> mat;
            std::string str; 
            if(file){
                while (std::getline(file, str))
                {
                    // Process str
                    
                    // log() << str << std::endl;
                    vector<std ::string> strs;
                    boost::split(strs, str, boost::is_any_of(" "));
                    std::vector<int> tmps;
                    for(int i=0; i < strs.size()-1; i++ ){
                        auto str_tmp = strs[i];
                        // log() << "str_tmp: " << str_tmp << std::endl;
                        tmps.push_back(std::stoi(str_tmp));
                    }
                    if(tmps.size()>0)
                        mat.push_back(tmps);
                }
            }else{
                log() << "could not find the file_name" << std::endl;
                std::exit(1);
            }
            perms[ii] = mat;

        }//end loop
    }//end 

    void readIellgalPlacementBoxes(std::string& bench_name){
        
        std::string file_name = "./../lookupTable/illegal_placement_boxs/"+\
                                    bench_name+".illegal"+".txt";
        std::ifstream file(file_name);
        std::string str; 
        if(file){
            while (std::getline(file, str))
            {
                illegal_placement_boxs_strings.push_back(str);
            }
        }else{
            log() << "could not find the file_name" << std::endl;
            std::exit(1);
        }
        

        
    }//end 

    void readOriginOffsetDie(std::string& bench_name){
        
        std::string file_name = "./../lookupTable/origin_offset_die/"+\
                                    bench_name+".offset"+".txt";
        std::ifstream file(file_name);
        std::string str; 
        if(file){
            while (std::getline(file, str))
            {
                origin_offset_die = std::stoi(str);
                break;
            }
        }else{
            log() << "could not find the file_name" << std::endl;
            std::exit(1);
        }
        

        
    }//end 

    // keep all permutations
    std::unordered_map<int,std::vector<std::vector<int>>> perms;
    std::vector<std::string> illegal_placement_boxs_strings;
    int origin_offset_die;


};//end class 
};//end namespace 


