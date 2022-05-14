#include "Legalizer.h"
#include <random>



namespace db {

// struct hash_tuple_leg {  // hash binary tuple
//     template <class T>
//     size_t operator()(const tuple<T, T>& tup) const {
//         auto hash1 = hash<T>{}(get<0>(tup));
//         auto hash2 = hash<T>{}(get<1>(tup));
//         return hash1 ^ hash2;
//     }
// };




void Legalizer::run(){
    bool debug = false || debug_global;

    if(debug)
        log() << "legalize cell: " << database.cells[cell_idx_].getName() << std::endl;

    utils::BoxT<DBU> legalize_box;
    std::vector<db::Cell> total_cells;
    std::set<int> movable_cells_set;
    std::vector<int> movable_cells;
    std::vector<int> legalize_rows;
    std::vector<int> legalize_sites;
    std::vector<std::vector<int>> blockage_matrix;

    movable_cells_set.insert(cell_idx_);

    utils::BoxT<DBU> cell_box = database.cells[cell_idx_].getCellBox();

    // 1- init legalize box
    bool isLegalized_boxValid = getLegalizeBox(cell_box,legalize_box);
    if(!isLegalized_boxValid) return;


    // 2- init legalize rows and sites
    database.getIntersectedRows(legalize_box,legalize_rows);
    database.getIntersectedSites(legalize_box,legalize_sites); 
    if(!isLegalizedRowsValid(legalize_rows) || 
       !isLegalizedSitesValid(legalize_sites) ) 
        return;
    
    // 3- update minMaxRowsSites
    updateMinMaxRows(legalize_rows);
    updateMinMaxSites(legalize_sites);

    // print legalizeRowsSites
    if(debug){
        log() << "legalize box: " << legalize_box << std::endl;
        printLegalizeRows(legalize_rows);
        printLegalizeSites(legalize_sites);
        auto cell_width = database.cells[cell_idx_].getCellSiteWidth();
    }

    
    // 4- get initCellsInBox
    database.getCellsInBox(legalize_box,total_cells,10);

    // 5- getMovableCellsInsideLegalizeBox
    getMovableCellsInsideLegalizeBox(total_cells,legalize_box,movable_cells_set);


    // print movable cells and involved cells
    if(debug){
        log() << "involved cells" << std::endl;
        for(auto cell : total_cells){
            log() << "cell: " << cell.getName()
                << ", isFixed: " << cell.isFixed 
                << ", s: "  << cell.getCellSite()
                << ", r: "  << cell.getCellRow()
                << ", W: "  << cell.getCellSiteWidth()
                << ", box: "<< cell.getCellBox() << std::endl;
        }
    }
    

    

    int total_empty_space = 
        getLegalizerBlockageMatrix(   total_cells
                                    , movable_cells_set
                                    , legalize_rows
                                    , legalize_sites
                                    , blockage_matrix);

    
    // print movable_cells_str
    if(db::setting.debug){
        auto min_row = *std::min_element(legalize_rows.begin(),legalize_rows.end());
        auto max_row = *std::max_element(legalize_rows.begin(),legalize_rows.end());
        auto min_site = *std::min_element(legalize_sites.begin(),legalize_sites.end());
        auto max_site = *std::max_element(legalize_sites.begin(),legalize_sites.end());

        auto min_row_dbu = database.getDBURow(min_row);
        auto max_row_dbu = database.getDBURow(max_row+1);
        auto min_site_dbu = database.getDBUSite(min_site);
        auto max_site_dbu = database.getDBUSite(max_site+1);

        std::string movable_cells_str="";
        std::string total_cells_str="";

        for(auto cell_tmp : total_cells){
            total_cells_str=total_cells_str+"|"+cell_tmp.getName();
        }

        for(auto idx_tmp : movable_cells_set){
            movable_cells_str=movable_cells_str+"|"+database.cells[idx_tmp].getName();
        }


        database.ss_legalizer << database.cells[cell_idx_].getName()
                      << "," <<min_site_dbu
                      << "," <<min_row_dbu
                      << "," <<max_site_dbu
                      << "," <<max_row_dbu
                      << "," <<total_cells_str
                      << "," <<movable_cells_str << std::endl;                
    }

    
    
    int total_cells_width = 0;
    for(auto cell_idx : movable_cells_set){
        auto cell = database.cells[cell_idx];
        if(!cell.isFixed) {
            movable_cells.push_back(cell_idx);
            total_cells_width += cell.getCellSiteWidth();
        }
    }//end for 

    if(debug){
        log() << "total_empty_space: " << total_empty_space << std::endl;
        log() << "total_required_cell_space: " << total_cells_width << std::endl;
        printMovableCells(movable_cells);
    }
    if(total_empty_space < total_cells_width ) return;

  
    if(database.policy_set.find("legalizerILP") == database.policy_set.end())
        legalizerILPV2( movable_cells
                    , legalize_rows
                    , legalize_sites
                    , blockage_matrix);
}//end run function


bool Legalizer::getLegalizeBox(utils::BoxT<DBU>& new_box, utils::BoxT<DBU>& legalize_box){
    bool debug = false || debug_global;
    // log() << "cell_box: " << new_box << std::endl;
    auto cell_site = database.cells[cell_idx_].getCellSite();
    auto cell_row = database.cells[cell_idx_].getCellRow();
    auto cell_width = database.cells[cell_idx_].getCellSiteWidth();
    auto site_step = cell_width * db::setting.legalizationWindowSizeSite;
    // auto site_step = db::setting.legalizationWindowSizeSite;
    if(debug){
        log() << "cell_site: " << cell_site
          << ", cell_row: " << cell_row 
          << ", site_step: " << site_step << std::endl;
        log() << "cell_site dbu: " << database.getDBUSite(cell_site)
          << ", cell_row dbu: " << database.getDBURow(cell_row)
          << std::endl;
    }
    


    if(cell_site < 0 || cell_row < 0) return false;

    int min_site = 0;
    int max_site = 0;
    int min_row = 0;
    int max_row = 0;

    if(cell_site - site_step <= 0){
        min_site = 0;
        max_site = cell_site + 2*site_step + cell_width;
    }else if (cell_site + cell_width + site_step >= database.sites_num){
        max_site = database.sites_num;
        min_site = max_site - 2*site_step - cell_width;
    }else{
        min_site = cell_site - site_step;
        max_site = cell_site + cell_width + site_step;
    }

    if(max_site >= database.sites_num) max_site = database.sites_num;
    if(min_site <= 0) min_site = 0;
   

    // log() << "site_step: " << site_step
    //       << "min_site: "  << min_site 
    //       << ", max_site: "  << max_site << std::endl;

    // if(cell_site - db::setting.legalizationWindowSizeSite/2 <= 0){
    //     min_site = 0;
    //     max_site = min_site + db::setting.legalizationWindowSizeSite;
    // }else if (cell_site + db::setting.legalizationWindowSizeSite/2 >= database.sites_num){
    //     max_site = database.sites_num;
    //     min_site = max_site - db::setting.legalizationWindowSizeSite;
    // }else{
    //     min_site = cell_site - db::setting.legalizationWindowSizeSite/2;
    //     max_site = cell_site + db::setting.legalizationWindowSizeSite/2;
    // }

    if(cell_row - db::setting.legalizationWindowSizeRow/2 <= 0){
        min_row = 0;
        max_row = min_row + db::setting.legalizationWindowSizeRow;
    }else if (cell_row + db::setting.legalizationWindowSizeRow/2 >= database.rows_num){
        max_row = database.rows_num;
        min_row = max_row - db::setting.legalizationWindowSizeRow;
    }else{
        min_row = cell_row - db::setting.legalizationWindowSizeRow/2;
        max_row = cell_row + db::setting.legalizationWindowSizeRow/2 + 1;
    }
    

    auto min_row_dbu = database.getDBURow(min_row);
    auto max_row_dbu = database.getDBURow(max_row);
    auto min_site_dbu = database.getDBUSite(min_site);
    auto max_site_dbu = database.getDBUSite(max_site);

    legalize_box.Set(min_site_dbu,min_row_dbu,max_site_dbu,max_row_dbu);

    // auto boundary_box_onLayer = db::BoxOnLayer(0, new_box);
    // auto grBox = grDatabase.rangeSearchGCell(boundary_box_onLayer);
    
    // auto dbu_site_min = database.getDBUSite(0);
    // auto dbu_site_max = database.getDBUSite(database.sites_num);
    // auto dbu_row_min = database.getDBURow(0);
    // auto dbu_row_max = database.getDBURow(database.rows_num);

    // log() << "database.sites_num: " << database.sites_num
    //       << ", dbu_site_max: " << double(dbu_site_max)/database.libDBU
    //       << ", dbu_site_min: " << double(dbu_site_min)/database.libDBU
    //       << ", dbu_row_min: " << double(dbu_row_min)/database.libDBU
    //       << ", dbu_row_max: " << double(dbu_row_max)/database.libDBU
    //       << ", database.rows_num: " << database.rows_num;
     

    // int ext_grBox_lx = 0 ;//grBox.lx()-db::setting.legalizationWindowSize;      
    // int ext_grBox_ly = 0 ;//grBox.ly()-db::setting.legalizationWindowSize;
    // int ext_grBox_hx = 0 ;//grBox.hx()+db::setting.legalizationWindowSize;
    // int ext_grBox_hy = 0 ;//grBox.hy()+db::setting.legalizationWindowSize+1;
    
   

    // auto lx = grDatabase.getCoor(ext_grBox_lx, X);
    // auto ly = grDatabase.getCoor(ext_grBox_ly, Y);
    // auto hx = grDatabase.getCoor(ext_grBox_hx, X);
    // auto hy = grDatabase.getCoor(ext_grBox_hy, Y);

    // if(hx < dbu_site_min) return false;
    // if(hy < dbu_row_min) return false;
    // if(lx > dbu_site_max) return false;
    // if(ly > dbu_row_max) return false;

    // if(lx < dbu_site_min) lx = dbu_site_min;
    // if(ly < dbu_row_min) ly = dbu_row_min;
    // if(hx > dbu_site_max) hx = dbu_site_max;
    // if(hy > dbu_row_max) hy = dbu_row_max;

    // auto row = database.cells[cell_idx_].getCellRow();

    // int min_site = 0;
    // int max_site = 0;
    // int min_row = 0;
    // int max_row = 0;

    // int min_site_dbu = 0;
    // int max_site_dbu = 0;
    // int min_row_dbu  = 0;
    // int max_row_dbu  = 0;

    // auto min_row = 0;
    // auto max_row = database.rows_num;

    // if(row == 0){
    //     min_row = 0;
    //     max_row = min_row + 2;
    // }else if (row == database.rows_num-1){
    //     max_row = row + 1;
    //     min_row = row - 2;
    // }else{
    //     min_row = row-1;
    //     max_row = row+2;
    // }

    // // log() << "min_row: " << min_row << ", max_row: " << max_row << std::endl;

    // auto min_row_dbu = database.getDBURow(min_row);
    // auto max_row_dbu = database.getDBURow(max_row);

    return true;
}//end getLegalizerBox

void Legalizer::updateMinMaxSites(std::vector<int>& legalize_sites){
    min_site_ = *std::min_element(legalize_sites.begin(),legalize_sites.end());
    max_site_ = *std::max_element(legalize_sites.begin(),legalize_sites.end());    
    site_offset_ = min_site_;
}//end updateMinMaxSites
void Legalizer::updateMinMaxRows(std::vector<int>& legalize_rows){
    min_row_ = *std::min_element(legalize_rows.begin(),legalize_rows.end());
    max_row_ = *std::max_element(legalize_rows.begin(),legalize_rows.end());    
    row_offset_ = min_row_;
}//end updateMinMaxRows


void Legalizer::getMovableCellsInsideLegalizeBox(std::vector<db::Cell>& cells
            ,utils::BoxT<DBU>& legalize_box
            ,std::set<int>& cells_inside){

    

    // std::set<int> cells_inside_tmp;
    std::vector<int> cells_inside_tmp;
    std::set<int> connected_cells_freeze;
    for(auto cell : cells){
        auto cell_box = cell.getCellBox();
        // log() << "inside cell loop: " << cell.getName() 
        //       << "leg_box: " << legalize_box
        //       << "cell_box: " << cell_box << std::endl;
        if((cell_box.lx() >= legalize_box.lx()) &&
           (cell_box.ly() >= legalize_box.ly()) &&
           (cell_box.hx() <= legalize_box.hx()) &&
           (cell_box.hy() <= legalize_box.hy()) ){
            //    log() << "added " << std::endl;
            //    cells_inside_tmp.insert(cell.idx);
                cells_inside_tmp.push_back(cell.idx);
           }
    }

    int maximum_number_of_cells_to_legalize = db::setting.legalizationMaxNumCellsToLegalize-1;
    int i = 0;

    std::sort(cells_inside_tmp.begin(),cells_inside_tmp.end(),[](int idx1,int idx2){
        return database.cells[idx1].getCellSiteWidth() < database.cells[idx1].getCellSiteWidth();
    });

     

    // old
    if(maximum_number_of_cells_to_legalize >= 1 ||
       db::setting.legalizationMaxNumCellsToLegalize == -1 )
        for(auto cell_idx : cells_inside_tmp){
            // log() << "inside box: " << database.cells[cell_idx].getName() 
            //       << ", isFixed: " << database.cells[cell_idx].isFixed
            //       << std::endl;
            if(connected_cells_freeze.find(cell_idx) != connected_cells_freeze.end())
                continue;

            if((!database.cells[cell_idx].isFixed) && (cell_idx != cell_idx_)){
                cells_inside.insert(cell_idx);
                
                auto connected_cells = database.cells[cell_idx].connectd_cells;
                // log() << "cell: " << database.cells[cell_idx].getName() << std::endl;
                // log() << "connected_cells size " << connected_cells.size() << std::endl;
                for(auto c : connected_cells){
                    connected_cells_freeze.insert(c);
                }

                i++;
                if(db::setting.legalizationMaxNumCellsToLegalize != -1)
                    if(i == maximum_number_of_cells_to_legalize){
                        break;
                    }
            }
        }


        // for(auto cell_idx : cells_inside_tmp){
        //     // log() << "inside box: " << database.cells[cell_idx].getName() 
        //     //       << ", isFixed: " << database.cells[cell_idx].isFixed
        //     //       << std::endl;
        //     if(connected_cells_freeze.find(cell_idx) != connected_cells_freeze.end())
        //         continue;

        //     if((!database.cells[cell_idx].isFixed) && (cell_idx != cell_idx_)){
        //         cells_inside.insert(cell_idx);
                
        //         auto connected_cells = database.cells[cell_idx].connectd_cells;
        //         // log() << "cell: " << database.cells[cell_idx].getName() << std::endl;
        //         // log() << "connected_cells size " << connected_cells.size() << std::endl;
        //         for(auto c : connected_cells){
        //             connected_cells_freeze.insert(c);
        //         }

                
        //     }
        // }
}//end getCellsInsideLegalizeBox

bool Legalizer::isLegalizedRowsValid(std::vector<int>& legalize_rows){
    if(legalize_rows.size() <= 0) return false;
    for(auto tmp : legalize_rows){
        if(tmp < 0){
            return false;
        } 
    }

    return true;
}//end isLegalizedRowsValid
bool Legalizer::isLegalizedSitesValid(std::vector<int>& legalize_sites){
    
    if(legalize_sites.size() <= 0) return false;

    for(auto tmp : legalize_sites){
        if(tmp < 0) {
            return false;
        }
    }
    return true;
}//end isLegalizedSitesValid

int Legalizer::getLegalizerBlockageMatrix(std::vector<db::Cell>& cells
                                  ,std::set<int>& cells_inside_set  
                                  ,std::vector<int>& legalize_rows
                                  ,std::vector<int>& legalize_sites
                                  ,std::vector<std::vector<int>>& blockage_matrix){
      
    blockage_matrix.resize(legalize_rows.size(), std::vector<int>(legalize_sites.size()));

    // Hard coded and temporary by erfan
    // illegal area to be place cell.
    // RTree illegal_placement_area_rtree;
    // int eps = 10;
    // boostBox box1(boostPoint(193600+eps,256200+eps),
    //              boostPoint(313600-eps,292200-eps));
    

    // utils::BoxT<DBU> illegal_box1(193600,256200,313600,292200);
    // utils::BoxT<DBU> illegal_box2(96200,196200,179600,234600);
    // utils::BoxT<DBU> illegal_box3(194800,133800,315600,174600);
    // utils::BoxT<DBU> illegal_box4(80000,76200,136800,109800);

    // std::vector<utils::BoxT<DBU>> illegal_boxs;
    // illegal_boxs.push_back(illegal_box1);
    // illegal_boxs.push_back(illegal_box2);
    // illegal_boxs.push_back(illegal_box3);
    // illegal_boxs.push_back(illegal_box4);
 

    // remove rows and sites that have overlap with 
    // fixed metal layers in upper layers.

    
    


    // illegal_placement_area_rtree.insert({box1, 0});

    // boostBox box1(boostPoint(193600+eps,256200+eps),
    //              boostPoint(313600-eps,292200-eps));
    // illegal_placement_area_rtree.insert({box1, 0});


    // end hard coded
    bool debug = false || debug_global;
    if(debug){
        log() << "legalize_rows: " << legalize_rows.size() 
              << "legalize_sites: " << legalize_sites.size() << std::endl;
    }

    int min_site = *std::min_element(legalize_sites.begin(),legalize_sites.end());
    int max_site = *std::max_element(legalize_sites.begin(),legalize_sites.end());
    int min_row = *std::min_element(legalize_rows.begin(),legalize_rows.end());
    int max_row = *std::max_element(legalize_rows.begin(),legalize_rows.end());
    int site_offset = min_site;
    int row_offset = min_row;
    int total_empty_width = 0;
    
    for(auto cell : cells){
        auto cell_box = cell.getCellBox();
        auto cell_row = cell.getCellRow();


        if(cells_inside_set.find(cell.idx) != cells_inside_set.end()){
            if(!cell.isFixed) {
                continue;
            }
        }
        

        auto lx_site = database.getSiteDBU(cell_box.lx());
        auto hx_site = database.getSiteDBU(cell_box.hx());

        for(int i = lx_site ; i < hx_site;i++){
            // if((min_site <= i) && (i <= max_site) &&
            //     (min_row <= cell_row) && (cell_row <= max_row)
            // ){
                // log() << "cell_row-row_offset: " << cell_row-row_offset
                //       << ", i-site_offset: " << i-site_offset << std::endl;
            if((cell_row-row_offset < legalize_rows.size()) &&
               (cell_row-row_offset >= 0) &&
               (i-site_offset < legalize_sites.size()) &&
               (i-site_offset >= 0)){
                    blockage_matrix[cell_row-row_offset][i-site_offset] = -1;
               }
                   
            // }
        }
    }

    for(int i = 0; i < blockage_matrix.size(); i++){
        if(database.illegal_placement_map.find(i+row_offset) != database.illegal_placement_map.end()){
            for(int j = 0; j < blockage_matrix[i].size(); j++){
                if(database.illegal_placement_map[i+row_offset].find(j+site_offset) 
                   != database.illegal_placement_map[i+row_offset].end()){
                    blockage_matrix[i][j] = -1;
                }
            }        
        }

        // add to avoid placement in special nets
        
    }
    
    for(auto tmps : blockage_matrix){
        for(auto tmp : tmps){
            if(tmp == 0)
                total_empty_width +=1;
        }
    }
    return total_empty_width;
}//end getLegalizerBlockageMatrix

void Legalizer::getAvailableSpaceInEachRow(std::vector<std::vector<int>>& blockage_matrix
                                , std::vector<std::pair<int,int>>& rowIdxToAvailableSpace){

    for(int r = 0; r < blockage_matrix.size();r++){
        auto sites = blockage_matrix[r];
        int row_empty_space = 0;
        for(auto site : sites){
            if(site == 0)
                row_empty_space +=1;
        }
        rowIdxToAvailableSpace.push_back(std::make_pair(r,row_empty_space));
    }//end for 
                                    

}//end getAvailableSpaceInEachRow

bool Legalizer::isEnoughSpaceForCellInRow(int x
                            , int r
                            , std::vector<int>& movable_cells
                            , std::vector<int>& rowIdxToSpaceUsedInRow
                            , std::vector<std::pair<int,int>>& rowIdxToAvailableSpace){
    bool debug = false;
    int cell_idx = movable_cells[x];

                     
    if(rowIdxToSpaceUsedInRow[r] +
        database.cells[cell_idx].getCellSiteWidth() > 
        rowIdxToAvailableSpace[r].second)
            return false;
    return true;
}//end isEnoughSpaceForCellInRow

void  Legalizer::updateSpaceUsedInRow(int x
                            , int r
                            , std::vector<int>& movable_cells
                            , std::vector<std::vector<int>>& cellRow_table
                            , std::vector<int>& rowIdxToSpaceUsedInRow
                            , bool increase){
    int cell_idx = movable_cells[x];
    if(increase){
        cellRow_table[x][r] = 1;
        rowIdxToSpaceUsedInRow[r] =
            rowIdxToSpaceUsedInRow[r] +
            database.cells[cell_idx].getCellSiteWidth();
    }else{
        cellRow_table[x][r] = 0;
        rowIdxToSpaceUsedInRow[r] =
            rowIdxToSpaceUsedInRow[r] -
            database.cells[cell_idx].getCellSiteWidth();
    }
        
}//end 

bool Legalizer::isSolutionBT(int placed_cells,std::vector<int>& movable_cells){
    if(placed_cells == movable_cells.size())
        return true;
    return false;
}//end isSolution

void Legalizer::collectSolutionBT(std::vector<int>& movable_cells
                                ,std::vector<int>& legalize_rows
                                ,std::vector<std::vector<int>>& cellRow_table
                                ,std::vector<std::pair<int,int>>& cellRow_sols){
    
    for(int i = 0; i < cellRow_table.size(); i++){
        for(int j = 0; j < cellRow_table[i].size();j++){
            if(cellRow_table[i][j]==1){
                auto cell_idx = movable_cells[i];
                int row_idx = legalize_rows[j];
                cellRow_sols.push_back(std::make_pair(cell_idx,row_idx));
            }
        }
    }                                
}//end collectioSolutionBT


void Legalizer::getAllPossibleArrangementsOfCellsInDifferentRows(
      std::vector<int>& movable_cells
    , std::vector<int>& legalize_rows
    , std::vector<std::pair<int,int>>& queue_bt
    , std::vector<int>& rowIdxToSpaceUsedInRow
    , std::vector<std::pair<int,int>>& rowIdxToAvailableSpace
    , std::vector<std::vector<int>>& cellRow_table
    , std::vector<std::vector<std::pair<int,int>>>& cellRow_sols_table
){
    bool debug = false;
    int x_bt = 0;
    int r_bt = 0;
    int placed_cells = 0;
    bool found_solution = false;
    bool isFound = false;
    int sol_i = 0;
     while(true){
        for(int r = r_bt; r < legalize_rows.size(); r++){
            found_solution = false;            
            if(isEnoughSpaceForCellInRow(x_bt, r,movable_cells
                                        , rowIdxToSpaceUsedInRow
                                        , rowIdxToAvailableSpace)){
                updateSpaceUsedInRow(x_bt,r,movable_cells
                                    , cellRow_table
                                    , rowIdxToSpaceUsedInRow
                                    , true);
                placed_cells++;
                isFound = true;
                queue_bt.push_back(std::make_pair(x_bt,r));
                // if found a solution
                if(placed_cells == movable_cells.size()){
                    std::vector<std::pair<int,int>> cellRow_sols;
                    collectSolutionBT( movable_cells
                                , legalize_rows
                                , cellRow_table
                                , cellRow_sols);
                    cellRow_sols_table.push_back(cellRow_sols);

                    sol_i++;

                    found_solution = true;
                }
                break;
            }//end if 
        }//end loop-r

        if(found_solution){
            if(queue_bt.empty()){
                break;
            }
            auto x_r_pair = queue_bt.back();
            queue_bt.pop_back();
            x_bt = x_r_pair.first;
            r_bt = x_r_pair.second+1;
            int r_bt_old = x_r_pair.second;
            updateSpaceUsedInRow(x_bt,r_bt_old,movable_cells
                    , cellRow_table
                    , rowIdxToSpaceUsedInRow
                    , false);
            isFound = false;
            placed_cells--;
        }else{
            if(isFound){
                x_bt++;
                r_bt = 0;
                isFound=false;
            }else{
                if(queue_bt.empty()){
                    break;
                }
                auto x_r_pair = queue_bt.back();
                queue_bt.pop_back();
                x_bt = x_r_pair.first;
                r_bt = x_r_pair.second+1;
                int r_bt_old = x_r_pair.second;
                
                updateSpaceUsedInRow(x_bt,r_bt_old,movable_cells
                    , cellRow_table
                    , rowIdxToSpaceUsedInRow
                    , false);
                isFound = false;
                placed_cells--;
                
            }// end if isFound
        }//end if found_solution
    }//end while
}//end getAllPossibleArrangementsofCellsInDifferentRows


void Legalizer::legalizerBacktracking(  std::vector<int>& movable_cells
                            , std::vector<int>& legalize_rows
                            , std::vector<int>& legalize_sites
                            , std::vector<std::vector<int>>& blockage_matrix){
    std::vector<std::pair<int,int>> rowIdxToAvailableSpace;
    std::vector<int> rowIdxToSpaceUsedInRow;
    std::vector<std::vector<std::pair<int,int>>> cellRow_sols_table;
    // cell_idx to row
    std::vector<std::pair<int,int>> queue_bt;
    // cell-row table
    std::vector<std::vector<int>> cellRow_table;


    rowIdxToSpaceUsedInRow.resize(legalize_rows.size());
    cellRow_table.resize(movable_cells.size(), std::vector<int>(legalize_rows.size()));
    
    getAvailableSpaceInEachRow(blockage_matrix,rowIdxToAvailableSpace);


    getAllPossibleArrangementsOfCellsInDifferentRows(movable_cells
                                                    ,  legalize_rows
                                                    ,  queue_bt
                                                    ,  rowIdxToSpaceUsedInRow
                                                    ,  rowIdxToAvailableSpace
                                                    ,  cellRow_table
                                                    ,  cellRow_sols_table);
                                                    
    log() << "cellRow_sols_table size: " << cellRow_sols_table.size() << std::endl;
    for(int i =0; i < cellRow_sols_table.size() ; i++){
        log() << "sol_" << i << std::endl;
        for(int j = 0; j < cellRow_sols_table[i].size(); j++ ){
            log() << "sol: cell: " << database.cells[cellRow_sols_table[i][j].first].getName()
                  << ", row: " << cellRow_sols_table[i][j].second << std::endl;

        }//end for j
    }//end for i

}//end LegalizerBacktracking algorithm


void Legalizer::printMovableCells(std::vector<int>& movable_cells){
    for(auto cell_idx : movable_cells){
        log() << "mov_cell: " << database.cells[cell_idx].getName() << std::endl;
    }
}
void Legalizer::printLegalizeRows(std::vector<int>& legalize_rows){
    log() << "legalize_rows: " << std::endl;
    for(auto row : legalize_rows){
        log() << "row: " << row << "-> " << database.getDBURow(row)<< std::endl;
    }
}
void Legalizer::printLegalizeSites(std::vector<int>& legalize_sites){
        log() << "legalize_sites: " << std::endl;
        for(auto site : legalize_sites){
            log() << "site: " << site << "-> " << database.getDBUSite(site)<< std::endl;
        }
}

void Legalizer::legalizerILPV2(std::vector<int>& movable_cells
                    , std::vector<int>& legalize_rows
                    , std::vector<int>& legalize_sites
                    , std::vector<std::vector<int>>& blockage_matrix){
    std::vector<cellWrap> weights;
    std::vector<std::vector<int>> ovrlps;
    std::vector<std::vector<int>> constraints;
    std::vector<int> sols;


    if(db::setting.debug)
        std::cout << "cell: " << database.cells[cell_idx_].getName() << std::endl;
    
    bool debug = false || debug_global;
    if(db::setting.debug)
        std::cout << "getLegalizerILPCostCubeV3" << std::endl;
    getLegalizerILPCostCubeV3( weights
                    , legalize_rows
                    , legalize_sites
                    , movable_cells
                    , blockage_matrix);

    if(debug){
        logLegalizerWeights(weights,legalize_rows,legalize_sites);
    }
    
        
    // return;  

    // if(debug)
    //     for(auto weight : weights){
    //         log() << "cell: " << database.cells[weight.idx].getName()
    //             << ", r: " << weight.r
    //             << ", s: " << weight.s
    //             << ", cost: " << weight.cost << std::endl;
    //     }

    
    if(db::setting.debug)
        std::cout << "getLegalizerILPOvrlpConflictMatrixV2" << std::endl;
    getLegalizerILPOvrlpConflictMatrixV2(
                    weights
                    , ovrlps
                    );

    if(debug){
        // std::vector<cellWrap> weights;
        // std::vector<std::vector<int>> ovrlps;
        logLegalizerOvrlps(weights,legalize_rows,legalize_sites,ovrlps);
    }

    if(debug){
        log() << "ovrlps: " << ovrlps.size() << std::endl;
        for(auto ovrlp_vec : ovrlps){
            log() << "ovrlp: " << std::endl;
            for(auto ovrlp_idx : ovrlp_vec){
                int s = weights[ovrlp_idx].s;
                int r = weights[ovrlp_idx].r;
                int w = database.cells[weights[ovrlp_idx].idx].getCellSiteWidth();
                auto name = database.cells[weights[ovrlp_idx].idx].getName();
                
                log() << "name: " << name
                      << ", s: " << s
                      << ", r: " << r
                      << ", w: "  << w
                      << std::endl;
            }//end loop
            
        }//end loop
    }


    if(db::setting.debug)
        std::cout << "getLegalizerILPSinglePositionConflictMatrixV3" << std::endl;
    getLegalizerILPSinglePositionConflictMatrixV3(
              weights
            , movable_cells
            , constraints);
        
    if(db::setting.debug)
        std::cout << "getLegalizerSolutionILPSolverV2" << std::endl;
    getLegalizerSolutionILPSolverV2(
                weights
                , legalize_rows
                , legalize_sites
                , ovrlps
                , constraints
                , sols);

    if(debug){
        logLegalizerSols( weights
                        , legalize_rows
                        , legalize_sites
                        , sols);
    }
    

    if(db::setting.debug)
        log() << "sols size out: " << sols.size() << std::endl;



    for(auto sol : sols){
        int cell_idx = weights[sol].idx;
        int new_row  = legalize_rows[weights[sol].r];
        int new_site = legalize_sites[weights[sol].s];;
        legalizer_sols.push_back(std::make_tuple(cell_idx,new_row,new_site));
    }//end sols 
    
    if(db::setting.debug)
        log() << "legalizer_sols: " << legalizer_sols.size() << std::endl;
}//end legalizerILPV2


void Legalizer::legalizerILP(std::vector<int>& movable_cells
                , std::vector<int>& legalize_rows
                , std::vector<int>& legalize_sites
                , std::vector<std::vector<int>>& blockage_matrix){

    std::vector<double> weights;
    std::unordered_map<int,int> weightToCellIdxDict;
    std::unordered_map<int,std::tuple<int,int,int>> weigthToCellRowSite;
    std::vector<std::vector<std::vector<int>>> cube_cost;
    std::vector<std::vector<std::vector<int>>> cube_weightIdx;
    std::vector<std::vector<int>> conflict_matrix_ovrlp;
    std::vector<std::vector<int>> conflict_matrix_single_position;
    std::vector<int> sols;
    std::vector<std::vector<int>> sols_tb;
    std::vector<cellWrap> cell_wraps;
    
    bool debug = true;

    
    if(debug) log() << "getLegalizerILPCostCube" << std::endl;
    if(database.policy_set.find("getLegalizerILPCostCube") == database.policy_set.end())
        // getLegalizerILPCostCubeV2( weights
        //                     , legalize_rows
        //                     , legalize_sites
        //                     , movable_cells
        //                     , blockage_matrix
        //                     , cell_wraps
        //                     // , weightToCellIdxDict
        //                     // , weigthToCellRowSite
        //                     , cube_cost
        //                     , cube_weightIdx);

        getLegalizerILPCostCube( weights
                            , legalize_rows
                            , legalize_sites
                            , movable_cells
                            , blockage_matrix
                            , weightToCellIdxDict
                            , weigthToCellRowSite
                            , cube_cost
                            , cube_weightIdx);
        
    if(debug) log() << "end getLegalizerILPCostCube" << std::endl;

    if(debug) log() << "getLegalizerILPOvrlpConflictMatrix" << std::endl;

    if(database.policy_set.find("getLegalizerILPOvrlpConflictMatrix") == database.policy_set.end())
        getLegalizerILPOvrlpConflictMatrix(
                    cube_cost
                    , cube_weightIdx
                    , movable_cells
                    , legalize_rows
                    , legalize_sites
                    , conflict_matrix_ovrlp);
    if(debug) log() << "end getLegalizerILPOvrlpConflictMatrix" << std::endl;


    if(debug) log() << "getLegalizerILPSinglePositionConflictMatrix" << std::endl;
    if(database.policy_set.find("getLegalizerILPSinglePositionConflictMatrix") == database.policy_set.end())
        getLegalizerILPSinglePositionConflictMatrix(
                    cube_cost
                    , cube_weightIdx
                    , movable_cells
                    , conflict_matrix_single_position);
    if(debug) log() << "end getLegalizerILPSinglePositionConflictMatrix" << std::endl;
    // // logs
    // logLegalizerBlockageMatrix(blockage_matrix);
    // logCubeCost(  movable_cells
    //            , cube_cost
    //            , cube_weightIdx);
    // logCubeWeight(  movable_cells
    //             , cube_weightIdx);
    // return;
    if(!db::setting.legalizationActivateSolutionPool){
        if(debug) log() << "getLegalizerSolutionILPSolver" << std::endl;
        getLegalizerSolutionILPSolver(
                    weights
                    , conflict_matrix_ovrlp
                    , conflict_matrix_single_position
                    , weightToCellIdxDict
                    , weigthToCellRowSite
                    , sols);
        if(debug) log() << "end getLegalizerSolutionILPSolver" << std::endl;
        for(auto sol : sols){
            int cell_idx = std::get<0>(weigthToCellRowSite[sol]);
            int new_row  = std::get<1>(weigthToCellRowSite[sol]) + row_offset_;
            int new_site = std::get<2>(weigthToCellRowSite[sol]) + site_offset_;
            legalizer_sols.push_back(std::make_tuple(cell_idx,new_row,new_site));
        }//end sols 
    }else{
        if(debug) log() << "getLegalizerSolutionILPSolverSolutionPool" << std::endl;
        getLegalizerSolutionILPSolverSolutionPool(
                        weights
                        , conflict_matrix_ovrlp
                        , conflict_matrix_single_position
                        , weightToCellIdxDict
                        , sols_tb);
            

        for(auto sols_row : sols_tb){
            for(auto sol : sols_row){
                int cell_idx = std::get<0>(weigthToCellRowSite[sol]);
                int new_row  = std::get<1>(weigthToCellRowSite[sol]) + row_offset_;
                int new_site = std::get<2>(weigthToCellRowSite[sol]) + site_offset_;
                legalizer_sols.push_back(std::make_tuple(cell_idx,new_row,new_site));
            }
        }//end sols 
        if(debug) log() << "end getLegalizerSolutionILPSolverSolutionPool" << std::endl;
    }//end if-else
    


    


}//end legalizerILP

void Legalizer::getLegalizerILPCostCube(std::vector<double>& weights
                        , std::vector<int>& legalize_rows
                        , std::vector<int>& legalize_sites
                        , std::vector<int>& movable_cells
                        , std::vector<std::vector<int>>& blockage_matrix
                        , std::unordered_map<int,int>& weightToCellIdxDict
                        , std::unordered_map<int,std::tuple<int,int,int>>& weigthToCellRowSite
                        , std::vector<std::vector<std::vector<int>>>& cube_cost
                        , std::vector<std::vector<std::vector<int>>>& cube_weightIdx){
    // weight idx to cell idx

    cube_weightIdx.resize(movable_cells.size());
    for (int r = 0;r< movable_cells.size(); r++) {        
        cube_weightIdx[r].resize(legalize_rows.size(), std::vector<int>(legalize_sites.size()));
    }

    int weight_idx=0;
    for(int movable_cells_idx = 0;movable_cells_idx < movable_cells.size(); movable_cells_idx++){
        auto cell_idx = movable_cells[movable_cells_idx];
        auto cell = database.cells[cell_idx];
        
        auto medianTmp = cell.feature.median_x_y_pin;//getMedianPin();
        int median_siteTmp = database.getSiteDBU(medianTmp.first);
        int median_rowTmp  = database.getRowDBU(medianTmp.second);
        // log() << "cell name: " << cell.getName()  
        //       << ", median: " << medianTmp 
        //       << ", median_siteTmp: " << median_siteTmp 
        //       << ", median_rowTmp: " << median_rowTmp << std::endl;
        auto cell_width = cell.getCellSiteWidth();
        std::vector<std::vector<int>> rows_vec;
        for(int r = 0; r < legalize_rows.size();r++){
            std::vector<int> sites_vec;
            for(int s = 0; s < legalize_sites.size(); s++){
                int cost_median = 0;
                int cost = 0;
                int row_idx  = legalize_rows[r];
                int site_idx = legalize_sites[s];
                int row_dbu  = database.getDBURow(row_idx);
                int site_dbu = database.getDBUSite(site_idx);
                utils::BoxT<DBU> new_cell_box(site_dbu,row_dbu
                                        ,site_dbu+(database.sites_step*cell_width)
                                        ,row_dbu+database.rows_step);
                // int site_idx_offset = site_idx - site_offset_;
                // int row_idx_offset = row_idx - row_offset_;
                
                // auto median = cell.feature.median_x_y_pin;//getMedianPin();
                auto median = cell.feature.median_x_y_pin;//getMedianPin();
                int median_site = database.getSiteDBU(median.first);
                int median_row  = database.getRowDBU(median.second);
                int diff_x = database.sites_step*std::abs(median_site-site_idx);
                int diff_y = database.rows_step*std::abs(median_row-row_idx);
                cost_median = diff_x + diff_y;
                // if(median_row == row_idx)
                //     cost = 2.5*cost_median;
                // else
                if(database.policy_set.find("suspected_cells_dict") == database.policy_set.end())
                    if(database.suspected_cells_dict.find(cell.idx) != database.suspected_cells_dict.end() ){
                        auto penalty_box = database.suspected_cells_dict[cell.idx];
                        if(penalty_box.HasIntersectWith(new_cell_box))
                            cost = 500*cost_median;
                    }else{
                        cost = cost_median;
                    }
                        
                    
                    
                // log() << "cost: " << cost << ", max_site: " << max_site_ << std::endl;
                if(database.policy_set.find("blockage_matrix[r][ss]") == database.policy_set.end())
                    for(int ss = s; ss < s+cell_width;ss++ ){
                        // log() << "r: " << r
                        //       << ", ss: " << ss 
                        //       << ", blockage_matrix[r][ss]: " << blockage_matrix[r][ss]  << std::endl;
                        if(blockage_matrix[r][ss] == -1){
                            cost = -1;
                            break;
                        }
                    }
                
                if(site_idx+cell_width > max_site_+1){
                    cost = -1;
                }
                // log() << "cost_after site_idx+cell_width_sites > max_site: "  << cost << std::endl;


                sites_vec.push_back(cost);
                // log() << "r: " << r
                //       << ", s: " << s 
                //       << ", weight_idx: " << weight_idx
                //       << ", cost: " << cost << std::endl;
                weights.push_back(cost);

                if(database.policy_set.find("cube_weightIdx[movable_cells_idx][r][s]") == database.policy_set.end())
                    cube_weightIdx[movable_cells_idx][r][s] = weight_idx;
                weightToCellIdxDict[weight_idx] = cell.idx;
                // weigthToCellRowSite[weight_idx] = std::make_tuple(cell.idx,row_idx,site_idx);
                if(database.policy_set.find("weigthToCellRowSite[weight_idx]") == database.policy_set.end())
                    weigthToCellRowSite[weight_idx] = std::make_tuple(cell.idx,r,s);
                weight_idx++;
            }//end site
            rows_vec.push_back(sites_vec);
        }//end for 
        cube_cost.push_back(rows_vec);
    }//end for 
}//end getLegalizerILPCostCube

// void Legalizer::findAllPermutations(std::vector<cellWrap>& weights
//                 , std::vector<int>& legalize_rows
//                 , std::vector<int>& legalize_sites
//                 , std::vector<int>& movable_cells
//                 , std::vector<std::vector<int>>& mat_spr
//                 , std::vector<std::vector<int>>& blockage_matrix
//                 // , std::unordered_map<int,int>& weightToCellIdxDict
//                 // , std::unordered_map<int,std::tuple<int,int,int>>& weigthToCellRowSite
//                 // , std::vector<std::vector<std::vector<int>>>& cube_cost
//                 // , std::vector<std::vector<std::vector<int>>>& cube_weightIdx
//                 ){
//     auto rsynService = database.getRsynService();
//     bool debug = false || debug_global;
//     if(debug)
//         log() << "findAllPermutations..." << std::endl;
//     std::stringstream ss;
//     ss << "inst,r,s,w,h,cost" << std::endl;  
//     // find segments

//     std::vector<std::vector<std::pair<int,int>>> segs;

//     segs.resize(legalize_rows.size());

    
//     for(int r = 0; r < mat_spr.size() ; r++){
//         std::vector<int> seqs;
//         for(int s = 0; s < mat_spr[r].size(); s++){
//             seqs.push_back(mat_spr[r][s]);
//             if(s+1 < mat_spr[r].size()){
//                 if(std::abs(mat_spr[r][s]-mat_spr[r][s+1]) != 1){
//                     int min = *std::min_element(seqs.begin(),seqs.end());
//                     int max = *std::max_element(seqs.begin(),seqs.end());
//                     segs[r].push_back(std::make_pair(min,max));
//                     seqs.clear();
//                 }//end if 
//             }else{
//                 int min = *std::min_element(seqs.begin(),seqs.end());
//                 int max = *std::max_element(seqs.begin(),seqs.end());
//                 segs[r].push_back(std::make_pair(min,max));
//                 seqs.clear();

//             }//end if 
//         }//end for 
//     }//end for 


//     if(debug){
//         log() << "seqs ..." << std::endl;
//         for(int r = 0; r < segs.size() ; r++){
//             log() << "r: " << r 
//                   << "->" << database.getDBURow(r) << std::endl;
//             auto seqs = segs[r];
//             std::string txt = "";
//             for(auto pair : seqs){
//                 txt = txt 
//                     + "(" 
//                     + std::to_string(pair.first) 
//                     + ","
//                     + std::to_string(pair.second)  
//                     + ") ->"
//                     + "(" 
//                     + std::to_string(database.getDBUSite(pair.first)) 
//                     + ","
//                     + std::to_string(database.getDBUSite(pair.second))  
//                     + ")";
//             }
//             log() << txt << std::endl;
//         }
//     }//end if 

    
//     std::vector<int> orders;
//     std::vector<std::vector<int>> orders_2d; 
//     std::vector<int> cell_w;


//     for(int i = 0; i < movable_cells.size(); i++){
        
//         orders.push_back(i);
//         auto cell = database.cells[movable_cells[i]];
//         cell_w.push_back(cell.getCellSiteWidth());
//         if(debug){
//             log() << "cell: " << cell.getName()
//               << ", i: " << i
//               << ", w: " << cell_w[cell_w.size() -1 ] << std::endl;
//         }
        
//     }

//     for(int i = 0; i < orders.size(); i++){
//         orders_2d.push_back(std::vector<int>{i});
//     }

//     for(int n =2; n <= orders.size(); n++){
//         do{
//             std::string txt2 = "";
//             for(int ii = 0; ii < orders.size(); ii++){
//                 txt2 = txt2 + std::to_string(orders[ii]) + " ";
//             }
//             log() << "current order: " << txt2 << std::endl;
            
                

//             //Display the current permutation
//             std::vector<int> order_tmp;
//             for(int i=0;i<n;i++) 
//                 order_tmp.push_back(orders[i]);
//             orders_2d.push_back(order_tmp);
//             //     cout << orders[i] << " ";
//             // cout << endl;
//         }while(std::next_permutation(orders.begin(), orders.end()));
//     }//end for 

//     if(debug){
//         for(auto ords : orders_2d){
//             std::string txt = "";
//             for(auto tmp : ords){
//                 txt = txt + std::to_string(tmp) + " ";
//             }
//             log() << txt << std::endl;
//         }
//     }


//     std::vector<Segment> segments;


//     for(int r = 0; r < segs.size(); r++){
//         for(int i = 0; i < segs[r].size(); i++){
//             // get orders 
//             int start = segs[r][i].first;
//             int stop = segs[r][i].second;
//             Segment seg_obj;
//             seg_obj.r = r;
//             seg_obj.start = start;
//             seg_obj.stop = stop;
            
//             for(int ii = 0; ii < orders_2d.size(); ii++){
//                 int sum = 0;
//                 std::string txt = "";
//                 for(int jj = 0; jj < orders_2d[ii].size(); jj++){
//                     txt = txt + std::to_string(orders_2d[ii][jj]) 
//                         + " " + std::to_string(cell_w[orders_2d[ii][jj]])
//                         + " " ;
//                     sum += cell_w[orders_2d[ii][jj]];
//                 }//end for
//                 if(debug){
//                     log() << "order: " << txt << std::endl; 
//                     log()  << "r: " << r
//                         << ", start: " << seg_obj.start
//                         << ", stop: " << seg_obj.stop
//                         << ", sum: " << sum << std::endl;
//                 }
                
//                 if(sum <= std::abs(start-stop)+1){
//                     seg_obj.orders.push_back(orders_2d[ii]);
//                 }
//             }//end for

//             segments.push_back(seg_obj);
//         }//end for 
//     }//end for

    
//     if(debug){
//         for(auto seg_obj : segments){
            
//             log() << "r: " << seg_obj.r
//                   << ", start: " << seg_obj.start
//                   << ", stop: " << seg_obj.stop
//                   << ", possible orders: " << std::endl;


//             for(auto ords : seg_obj.orders){
//                 std::string txt = "";
//                 for(auto elem : ords ){
//                     txt = txt + std::to_string(elem) + " ";
//                 }
//                 log() << txt << std::endl;
//             }
            
//         }
//     }//end if debug 

//     if(debug){
//         log() << "cell indexing ... " << std::endl;
//     }//end if 

//     for(auto seg : segments){
//         int r = seg.r;
//         auto orders_2d = seg.orders;
//         int seg_w = std::abs(seg.start - seg.stop)+1;
//         if(debug){
//             log() << "seg.r: " << r
//                   << ", seg.w: " << seg_w << std::endl;
//             for(auto ords : orders_2d){
//                 std::string txt = "";
//                 for(auto elem : ords ){
//                     txt = txt + std::to_string(elem) + " ";
//                 }
//                 log() << txt << std::endl;
//             }//end for 
//         }//end if 
        

//         for(auto order : orders_2d){
//             int cells_w_sum = 0;
//             std::string txt = "";


//             for(int idx : order){
//                 txt = txt + std::to_string(idx) + " ";
//                 cells_w_sum += cell_w[idx];
//             }

//             if(debug){
//                 log() << "order: " << txt << std::endl;
//                 log() << "cell_w_sum: " << cells_w_sum << std::endl;
//             }
                

//             int empty_w = seg_w-cells_w_sum;

//             if(debug)
//                 log() << "empty_w: " << empty_w << std::endl;


//             int num_line_seg = order.size();

//             int num_empty_spot = num_line_seg + 1;

//             int empty_offset = empty_w / num_empty_spot;
//             int reminder = empty_w % num_empty_spot;

            
//             int n = num_line_seg + num_empty_spot;
//             int idx = seg.start;

//             bool emptySw = true;
//             int i_order = 0;

//             if(debug){
//                 log() << "empty_offset: " << empty_offset
//                       << ", reminder: " << reminder
//                       << ", n: " << n << std::endl;
//             }

//             for(int i = 0; i < n; i++){
//                 if(idx >= seg.stop)
//                     continue;

//                 if(debug){
//                     log() << "idx: " << idx << std::endl;
//                 }
//                 if(emptySw){
//                     if(reminder > 0){
//                         std::random_device dev;
//                         std::mt19937 rng(dev());
//                         std::uniform_int_distribution<std::mt19937::result_type> rand_dist(0,reminder); // distribution in range [1, 6]
//                         int rnd_por = rand_dist(rng);
                        
//                         if(debug){
//                             log() << "rnd_por: " << rnd_por
//                                   << ", reminder: " << reminder 
//                                   << ", emptySw: " << emptySw << std::endl;
//                         }
                        
                        
//                         reminder = reminder - rnd_por;
//                     }

//                     idx = idx + empty_offset + 0;//reminder;
//                 }else{
//                     cellWrap cell_wrap;
                    
//                     auto cell = database.cells[movable_cells[order[i_order]]];
//                     cell_wrap.idx = cell.idx;
//                     cell_wrap.r = seg.r;
//                     cell_wrap.s = idx;

//                     // find cost 
//                     int row_idx  = legalize_rows[cell_wrap.r];
//                     int site_idx = legalize_sites[cell_wrap.s];
                                        
//                     // auto median = cell.feature.median_x_y_pin;//getMedianPin();
                    
                    
//                     // double cost_median = diff_x + diff_y;
//                     double cost_median = getCandidateCost(cell,row_idx,site_idx);
//                     if(cost_median > 0)
//                         cell_wrap.cost = cost_median;
//                     else
//                         cell_wrap.cost = 1;


//                     if(debug){
                            
//                             Rsyn::Cell cellRsyn = cell.rsynInstance.asCell();
//                             Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cellRsyn);
//                             Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cellRsyn);
//                             int cell_width_tmp = phLibCell.getWidth();
//                             int cell_height_tmp = phLibCell.getHeight();

//                             log() << cell.getName() 
//                                 << ",r: " << std::to_string(cell_wrap.r)
//                                 << ",s: " << std::to_string(cell_wrap.s)
//                                 << ",w: " << std::to_string(cell_width_tmp)
//                                 << ",h: " << std::to_string(cell_height_tmp)
//                                 << ",cost: " << std::to_string(cost_median) << std::endl;

//                             ss << cell.getName() 
//                                 << "," << std::to_string(cell_wrap.r)
//                                 << "," << std::to_string(cell_wrap.s)
//                                 << "," << std::to_string(cell_width_tmp)
//                                 << "," << std::to_string(cell_height_tmp)
//                                 << "," << std::to_string(cost_median) << std::endl;
                            
//                         }

//                     weights.push_back(cell_wrap);
//                     idx = idx + cell_w[order[i_order]];
//                     i_order++;
//                 }

//                 emptySw = emptySw xor true;

//                 if(debug){
//                     log() << std::endl;
//                 }
                
//             }//end for 



//         }//end order loop

//     }//end segments loop
    
//     if(debug){
//         std::string file_name_csv =db::setting.directory +  db::setting.benchmarkName 
//                                   + ".legalizer.csv";
//         // log() << "file_name_csv " << file_name_csv << std::endl;
//         std::ofstream fout(file_name_csv);
//         fout << ss.str();
//         fout.close();
//     }


// }//end findAllpermutation


void Legalizer::findAllPermutationsV2(std::vector<cellWrap>& weights
                , std::vector<int>& legalize_rows
                , std::vector<int>& legalize_sites
                , std::vector<int>& movable_cells
                , std::vector<std::vector<int>>& mat_spr
                , std::vector<std::vector<int>>& blockage_matrix
                // , std::unordered_map<int,int>& weightToCellIdxDict
                // , std::unordered_map<int,std::tuple<int,int,int>>& weigthToCellRowSite
                // , std::vector<std::vector<std::vector<int>>>& cube_cost
                // , std::vector<std::vector<std::vector<int>>>& cube_weightIdx
                ){
    auto rsynService = database.getRsynService();
    bool debug = false || debug_global;
    if(debug)
        log() << "findAllPermutations..." << std::endl;
    std::stringstream ss;
    // ss << "inst,r,s,w,h,cost" << std::endl;  
    ss << "inst_name,xl,yl,xh,yh,cost" << std::endl;  
    // find segments

    std::vector<std::vector<std::pair<int,int>>> segs;

    segs.resize(legalize_rows.size());

    // construct segs
    constructSegments(mat_spr,segs);

    


    if(debug){
        log() << "seqs ..." << std::endl;
        for(int r = 0; r < segs.size() ; r++){
            log() << "r: " << r 
                  << "->" << database.getDBURow(legalize_rows[r]) << std::endl;
            auto seqs = segs[r];
            std::string txt = "";
            for(auto pair : seqs){
                txt = txt 
                    + "(" 
                    + std::to_string(pair.first) 
                    + ","
                    + std::to_string(pair.second)  
                    + ") ->"
                    + "(" 
                    + std::to_string(database.getDBUSite(legalize_sites[pair.first])) 
                    + ","
                    + std::to_string(database.getDBUSite(legalize_sites[pair.second]))  
                    + ")";
            }
            log() << txt << std::endl;
        }
    }//end if 

    
    std::vector<int> orders;
    std::vector<std::vector<int>> orders_2d(database.lookup_tb.perms[movable_cells.size()]); 
    // sum of cells_w in each permutation
    std::vector<int> orders_2d_sum;
    std::vector<int> cell_w;

    assignMovableCellsToOrders( movable_cells, orders, cell_w);

    if(debug){
        log() << "database.lookup_tb.perms size: " 
              << database.lookup_tb.perms.size() << std::endl;
        for(auto ords : database.lookup_tb.perms[movable_cells.size()]){
            std::string txt = "";
            for(auto tmp : ords){
                txt = txt + std::to_string(tmp) + " ";
            }
            log() << txt << std::endl;
        }
    }


    std::vector<Segment> segments;

    initSegments( segs
                , orders_2d
                , cell_w
                , segments);

    

    

    
    if(debug){
        for(auto seg_obj : segments){
            
            log() << "r: " << seg_obj.r
                  << ", start: " << seg_obj.start
                  << ", stop: " << seg_obj.stop
                  << ", possible orders: " << std::endl;


            for(auto ords : seg_obj.orders){
                std::string txt = "";
                for(auto elem : ords ){
                    txt = txt + std::to_string(elem) + " ";
                }
                log() << txt << std::endl;
            }
            
        }
    }//end if debug 

    if(debug){
        log() << "cell indexing ... " << std::endl;
    }//end if 

   
    assignCellsToSegments(segments
                        , movable_cells
                        , legalize_rows
                        , legalize_sites
                        , cell_w
                        , weights
                        , ss);
    
    
    if(debug){
        std::string file_name_csv = db::setting.directory +  db::setting.benchmarkName 
             + ".legalizer.csv";
        // log() << "file_name_csv " << file_name_csv << std::endl;
        std::ofstream fout(file_name_csv);
        fout << ss.str();
        fout.close();
    }


}//end findAllpermutation



double Legalizer::getCandidateCost(db::Cell& cell, int row_idx, int site_idx){
    auto median = cell.feature.median_x_y_pin;//getMedianPin();
    // if(debug_global)
        // log() << "cell: "  << cell.getName() << ", "
        //     << "median.first: " << median.first << ", median.second" << median.second << std::endl;
    int median_site = database.getSiteDBU(median.first);
    int median_row  = database.getRowDBU(median.second);
    int diff_x = database.sites_step*std::abs(median_site-site_idx);
    int diff_y = database.rows_step*std::abs(median_row-row_idx);
    double cost = 2.5*diff_x + diff_y;


    // std::random_device dev;
    // std::mt19937 rng(dev());
    // std::uniform_int_distribution<std::mt19937::result_type> rand_dist(1,1000); // distribution in range [1, 6]
    // cost =  1;//rand_dist(rng);

    // if(cost <= 0)
    //     cost = 1;

    return cost;
}//end getCandidateCost


void Legalizer::getLegalizerILPCostCubeV3(std::vector<cellWrap>& weights
                        , std::vector<int>& legalize_rows
                        , std::vector<int>& legalize_sites
                        , std::vector<int>& movable_cells
                        , std::vector<std::vector<int>>& blockage_matrix
                        // , std::unordered_map<int,int>& weightToCellIdxDict
                        // , std::unordered_map<int,std::tuple<int,int,int>>& weigthToCellRowSite
                        // , std::vector<std::vector<std::vector<int>>>& cube_cost
                        // , std::vector<std::vector<std::vector<int>>>& cube_weightIdx
                        ){



    // weight idx to cell idx
    bool debug = false || debug_global || debug_global_all;

    auto rsynService = database.getRsynService();

    // matrix sparse
    std::vector<std::vector<int>> mat_spr;
    mat_spr.resize(blockage_matrix.size(),std::vector<int>{});

    if(debug){
        log() << "legalize_rows: " << legalize_rows.size() << std::endl;
        log() << "legalize_sites: " << legalize_sites.size() << std::endl;
    }

    if(debug){
        logLegalizerBoard(blockage_matrix,legalize_rows,legalize_sites);
    }
    

    for(int i = 0; i < blockage_matrix.size(); i++)
        for(int j = 0; j < blockage_matrix[i].size(); j++)
            if(blockage_matrix[i][j] != -1){
                mat_spr[i].push_back(j);
                // log() << "i: " << i << ", j: " << j << ", b: " << blockage_matrix[i][j] << std::endl;
            }

    // return;
    if(db::setting.findAllPermutations){
        findAllPermutationsV2( weights
                    , legalize_rows
                    , legalize_sites
                    , movable_cells
                    , mat_spr
                    , blockage_matrix);
    }else{
        std::stringstream ss;
        ss << "inst,r,s,w,h,cost" << std::endl;   

        int weight_idx=0;
        for(int movable_cells_idx = 0;movable_cells_idx < movable_cells.size(); movable_cells_idx++){
            auto cell_idx = movable_cells[movable_cells_idx];
            auto cell = database.cells[cell_idx];
            
            auto medianTmp = cell.feature.median_x_y_pin;//getMedianPin();
            int median_siteTmp = database.getSiteDBU(medianTmp.first);
            int median_rowTmp  = database.getRowDBU(medianTmp.second);
            auto cell_width = cell.getCellSiteWidth();
            int cell_width_tmp = 0;
            int cell_height_tmp = 0;
            // temporary need to be removed
            // width, height
            if(debug){
                Rsyn::Cell cellRsyn = cell.rsynInstance.asCell();
                Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cellRsyn);
                Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cellRsyn);
                cell_width_tmp = phLibCell.getWidth();
                cell_height_tmp = phLibCell.getHeight();
            }


            std::vector<std::vector<int>> rows_vec;
            for(int r = 0; r < mat_spr.size();r++){
                std::vector<int> sites_vec;
                // for(int s = 0; s < mat_spr[r].size(); s++){
                for(auto s : mat_spr[r]){
                    int cost_median = 0;
                    int cost = 0;
                    int row_idx  = legalize_rows[r];
                    int site_idx = legalize_sites[s];
                    
                    if(site_idx+cell_width > max_site_+1){
                        cost = -1;
                    }

                    if(cost != -1){
                        if(database.policy_set.find("blockage_matrix[r][ss]") == database.policy_set.end())
                        for(int ss = s; ss < s+cell_width;ss++ ){
                            // log() << "r: " << r
                            //       << ", ss: " << ss 
                            //       << ", blockage_matrix[r][ss]: " << blockage_matrix[r][ss]  << std::endl;
                            if(blockage_matrix[r][ss] == -1){
                                cost = -1;
                                break;
                            }
                        }
                    }
                    
                    // log() << "cost: " << cost << ", max_site: " << max_site_ << std::endl;
                    
                    
                    
                    // log() << "cost_after site_idx+cell_width_sites > max_site: "  << cost << std::endl;

                    if(cost != -1){

                        int row_dbu  = database.getDBURow(row_idx);
                        int site_dbu = database.getDBUSite(site_idx);
                        utils::BoxT<DBU> new_cell_box(site_dbu,row_dbu
                                                ,site_dbu+(database.sites_step*cell_width)
                                                ,row_dbu+database.rows_step);
                        // int site_idx_offset = site_idx - site_offset_;
                        // int row_idx_offset = row_idx - row_offset_;
                        
                        // auto median = cell.feature.median_x_y_pin;//getMedianPin();
                        auto median = cell.feature.median_x_y_pin;//getMedianPin();
                        int median_site = database.getSiteDBU(median.first);
                        int median_row  = database.getRowDBU(median.second);
                        int diff_x = database.sites_step*std::abs(median_site-site_idx);
                        int diff_y = database.rows_step*std::abs(median_row-row_idx);
                        cost_median = diff_x + diff_y;
                        // if(median_row == row_idx)
                        //     cost = 2.5*cost_median;
                        // else
                        if(database.policy_set.find("suspected_cells_dict") == database.policy_set.end())
                            if(database.suspected_cells_dict.find(cell.idx) != database.suspected_cells_dict.end() ){
                                auto penalty_box = database.suspected_cells_dict[cell.idx];
                                if(penalty_box.HasIntersectWith(new_cell_box))
                                    cost = 500*cost_median;
                            }else{
                                cost = cost_median;
                            }
                    }


                    sites_vec.push_back(cost);
                    // log() << "r: " << r
                    //       << ", s: " << s 
                    //       << ", weight_idx: " << weight_idx
                    //       << ", cost: " << cost << std::endl;

                    if(cost != -1){
                            if(debug){
                                ss << cell.getName() 
                                    << "," << std::to_string(r)
                                    << "," << std::to_string(s)
                                    << "," << std::to_string(cell_width_tmp)
                                    << "," << std::to_string(cell_height_tmp)
                                    << "," << std::to_string(cost) << std::endl;
                            }
                            if(cost != 0)
                                weights.push_back({cell.idx,r,s,double(cost)});
                            else
                                weights.push_back({cell.idx,r,s,double(1)});
                    }
                        

                    

                    // if(database.policy_set.find("cube_weightIdx[movable_cells_idx][r][s]") == database.policy_set.end())
                    //     cube_weightIdx[movable_cells_idx][r][s] = weight_idx;
                    // weightToCellIdxDict[weight_idx] = cell.idx;
                    // // weigthToCellRowSite[weight_idx] = std::make_tuple(cell.idx,row_idx,site_idx);
                    // if(database.policy_set.find("weigthToCellRowSite[weight_idx]") == database.policy_set.end())
                    //     weigthToCellRowSite[weight_idx] = std::make_tuple(cell.idx,r,s);
                    // weight_idx++;
                }//end site
                // rows_vec.push_back(sites_vec);
            }//end for 
            // cube_cost.push_back(rows_vec);
        }//end for 

        if(debug){
            std::string file_name_csv = db::setting.outputFile + "_legalizer.csv";
            // log() << "file_name_csv " << file_name_csv << std::endl;
            std::ofstream fout(file_name_csv);
            fout << ss.str();
            fout.close();
        }

    }// end if(db::setting.findAllPermutations)

}//end getLegalizerILPCostCube

void Legalizer::getLegalizerILPCostCubeV2(std::vector<double>& weights
                        , std::vector<int>& legalize_rows
                        , std::vector<int>& legalize_sites
                        , std::vector<int>& movable_cells
                        , std::vector<std::vector<int>>& blockage_matrix
                        , std::vector<cellWrap>& cell_wraps
                        // , std::unordered_map<int,int>& weightToCellIdxDict
                        // , std::unordered_map<int,std::tuple<int,int,int>>& weigthToCellRowSite
                        , std::vector<std::vector<std::vector<int>>>& cube_cost
                        , std::vector<std::vector<std::vector<int>>>& cube_weightIdx){
    // weight idx to cell idx
    utils::timer cube_cost_timer;
    bool debug = true;

    // matrix sparse
    std::vector<std::vector<int>> mat_spr;
    mat_spr.resize(blockage_matrix.size());

    for(int i = 0; i < blockage_matrix.size(); i++)
        for(int j = 0; j < blockage_matrix[i].size(); j++)
            if(blockage_matrix[i][j] != -1)
                mat_spr[i].push_back(j);


    // init space 
    std::stringstream ss_blockage;
    ss_blockage << "r,s,val,w,h" << std::endl;    
    for(int r = 0; r < legalize_rows.size() ; r++)
        for(int s = 0; s < legalize_sites.size() ; s++){
            int row_idx  = legalize_rows[r];
            int site_idx = legalize_sites[s];
            int row_dbu  = database.getDBURow(row_idx);
            int site_dbu = database.getDBUSite(site_idx);

            log() << "row_um: " <<  row_dbu/database.libDBU
                  << "site_um: " <<  site_dbu/database.libDBU << std::endl;

            ss_blockage << std::to_string(r)
                        << "," << std::to_string(s)
                        << "," << std::to_string(blockage_matrix[r][s])
                        << "," << std::to_string(database.getSiteStep())
                        << "," << std::to_string(database.getRowStep()) << std::endl;
        }

    
    std::string file_name_csv_b = db::setting.outputFile + "_initLegalizer.csv";
    log() << "file_name_csv " << file_name_csv_b << std::endl;
    std::ofstream fout_blockage(file_name_csv_b);
    fout_blockage << ss_blockage.str();
    fout_blockage.close();
    


    std::stringstream ss;
    ss << "inst,r,s,w,h,cost" << std::endl;    
    

    cube_weightIdx.resize(movable_cells.size());
    for (int r = 0;r< movable_cells.size(); r++) {        
        cube_weightIdx[r].resize(legalize_rows.size(), std::vector<int>(legalize_sites.size()));
    }

    

    // if(debug){
    //     log() << "movable_cells.size(): " << movable_cells.size() << std::endl;
    //     log() << "legalize_rows.size(): " << legalize_rows.size() << std::endl;
    //     log() << "legalize_sites.size(): " << legalize_sites.size() << std::endl;

    // }

    // int x = 1;
    // for(int i = 0 ; i < 3; i++)
    //     for(int j = 0 ; j < 5; j++)
    //         for(int k = 0 ; k < 40; k++)
    //             x = x + 1;

    // int x = 1;
    // for(int i = 0 ; i < 3; i++)
    //     for(int j = 0 ; j < 200; j++)
    //         // for(int k = 0 ; k < 40; k++)
    //             x = x + 1;

    
    // std::vector<cellWrap> cell_wraps;

    auto rsynService = database.getRsynService();
    

    int weight_idx=0;
    for(int movable_cells_idx = 0;movable_cells_idx < movable_cells.size(); movable_cells_idx++){
        auto cell_idx = movable_cells[movable_cells_idx];
        auto cell = database.cells[cell_idx];
        auto medianTmp = cell.feature.median_x_y_pin;//getMedianPin();
        int median_siteTmp = database.getSiteDBU(medianTmp.first);
        int median_rowTmp  = database.getRowDBU(medianTmp.second);

        // width, height
        Rsyn::Cell cellRsyn = cell.rsynInstance.asCell();
        Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cellRsyn);
        Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cellRsyn);
        int cell_width_tmp = phLibCell.getWidth();
        int cell_height_tmp = phLibCell.getHeight();

        
        log() << "cell name: " << cell.getName()  
              << ", median: " << medianTmp 
              << ", median_siteTmp: " << median_siteTmp 
              << ", median_rowTmp: " << median_rowTmp << std::endl;
        auto cell_width = cell.getCellSiteWidth();
        std::vector<std::vector<int>> rows_vec;
        for(int r = 0; r < mat_spr.size();r++){
            std::vector<int> sites_vec;
            for(int s = 0; s < mat_spr[r].size(); s++){
                int cost_median = 0;
                int cost = 0;
                int row_idx  = legalize_rows[r];
                int site_idx = legalize_sites[s];
                int row_dbu  = database.getDBURow(row_idx);
                int site_dbu = database.getDBUSite(site_idx);
                utils::BoxT<DBU> new_cell_box(site_dbu,row_dbu
                                        ,site_dbu+(database.sites_step*cell_width)
                                        ,row_dbu+database.rows_step);
                // int site_idx_offset = site_idx - site_offset_;
                // int row_idx_offset = row_idx - row_offset_;
                auto median = cell.feature.median_x_y_pin;//getMedianPin();
                int median_site = database.getSiteDBU(median.first);
                int median_row  = database.getRowDBU(median.second);

                // //temporary and need to be removed
                // // width, height
                // Rsyn::Cell cellRsyn = cell.rsynInstance.asCell();
                // Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cellRsyn);
                // Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cellRsyn);
                // int cell_width_tmp = phLibCell.getWidth();
                // int cell_height_tmp = phLibCell.getHeight();
                // // end tmp


                int diff_x = database.sites_step*std::abs(median_site-site_idx);
                int diff_y = database.rows_step*std::abs(median_row-row_idx);
                cost_median = diff_x + diff_y;
                if(database.suspected_cells_dict.find(cell.idx) != database.suspected_cells_dict.end() ){
                    auto penalty_box = database.suspected_cells_dict[cell.idx];
                    if(penalty_box.HasIntersectWith(new_cell_box))
                        cost = 500*cost_median;
                }else{
                    cost = cost_median;
                }//end if suspected_cells

                for(int ss = s; ss < s+cell_width;ss++ ){
                    // log() << "r: " << r
                    //       << ", ss: " << ss 
                    //       << ", blockage_matrix[r][ss]: " << blockage_matrix[r][ss]  << std::endl;
                    if(blockage_matrix[r][ss] == -1){
                        cost = -1;
                        break;
                    }
                }//end check overlap with blockages
                if(site_idx+cell_width > max_site_+1){
                    cost = -1;
                }
                // log() << "cost_after site_idx+cell_width_sites > max_site: "  << cost << std::endl;


                sites_vec.push_back(cost);
                // log() << "r: " << r
                //       << ", s: " << s 
                //       << ", weight_idx: " << weight_idx
                //       << ", cost: " << cost << std::endl;
                weights.push_back(cost);

                // ss << inst,r,s,w,h,cost
                ss << cell.getName() 
                   << "," << std::to_string(r)
                   << "," << std::to_string(s)
                   << "," << std::to_string(cell_width_tmp)
                   << "," << std::to_string(cell_height_tmp)
                   << "," << std::to_string(cost) << std::endl;

                cell_wraps.push_back({cell.idx,r,s});
                // cube_weightIdx[movable_cells_idx][r][s] = weight_idx;
                // weightToCellIdxDict[weight_idx] = cell.idx;
                // // weigthToCellRowSite[weight_idx] = std::make_tuple(cell.idx,row_idx,site_idx);
                
                // weigthToCellRowSite[weight_idx] = std::make_tuple(cell.idx,r,s);
                weight_idx++;

            }//end legalize_sites loop
            rows_vec.push_back(sites_vec);
        }//end legalize_rows loop
        cube_cost.push_back(rows_vec);
    }

    
    

    std::string file_name_csv = db::setting.outputFile + "_legalizer.csv";
    log() << "file_name_csv " << file_name_csv << std::endl;
    std::ofstream fout(file_name_csv);
    fout << ss.str();
    fout.close();
}//end getLegalizerILPCostCube


void Legalizer::getLegalizerILPOvrlpConflictMatrix(
      std::vector<std::vector<std::vector<int>>>& cube_cost
    , std::vector<std::vector<std::vector<int>>>& cube_weightIdx
    , std::vector<int>& movable_cells
    , std::vector<int>& legalize_rows
    , std::vector<int>& legalize_sites
    , std::vector<std::vector<int>>& conflict_matrix){
    // overlap conflicts
    
    for(int c = 0; c < cube_weightIdx.size() ; c++){
        for(int r = 0; r < cube_weightIdx[c].size(); r++){
            for(int s = 0; s < cube_weightIdx[c][r].size(); s++){
                if(cube_cost[c][r][s] == -1) continue;
                
                int weight_idx = cube_weightIdx[c][r][s];
                int cell_idx = movable_cells[c];
                int site_width = database.cells[cell_idx].getCellSiteWidth();
                // log() << "c: " << c << ", r: " << r << ", s: " << s << std::endl;
                // log() << "cell: " << database.cells[cell_idx].getName()
                //       << "c_w: " << site_width << std::endl;
                std::vector<int> conflict_matrix_row;
                conflict_matrix_row.push_back(weight_idx);
                for(int cc = 0; cc < cube_weightIdx.size() ; cc++){
                    if(c == cc ) continue;
                    for(int rr = 0; rr < cube_weightIdx[cc].size(); rr++){
                        if(r != rr) continue;
                        for(int ss = s; ss < s+site_width; ss++){
                            // log() << "cc: " << cc 
                            //           << ", rr: " << rr 
                            //           << ", ss: " << ss << std::endl;
                            // log() << "cube_cost[cc][rr][ss]: " << cube_cost[cc][rr][ss] << std::endl;
                            if(ss < legalize_sites.size()){
                                if(cube_cost[cc][rr][ss] == -1) continue;
                                int weight_idx_neigh = cube_weightIdx[cc][rr][ss];
                                // log() << "added" << std::endl;
                                conflict_matrix_row.push_back(weight_idx_neigh);
                                conflict_matrix.push_back(conflict_matrix_row);
                                conflict_matrix_row.pop_back();
                            }//end if
                        }//end loop ss
                    }//end loop rr
                }//end loop cc
                // if(conflict_matrix_row.size() > 1 )
                //     conflict_matrix.push_back(conflict_matrix_row);
            }//end loop s
        }//end loop r
    }// end loop c
}//end getLegalizerOvrlpConflictMatrix

void Legalizer::getLegalizerILPOvrlpConflictMatrixV2(
                std::vector<cellWrap>& weights
                , std::vector<std::vector<int>>& ovrlps){
    // overlap conflicts
    // int a = 0;
    double eps = 0.01;
    // bool debug = true;
    // log() << "weights.size(): "<<weights.size() << std::endl;
    // for(int i = 0; i < weights.size(); i++){
    //     for(int j = 0; j < weights.size(); j++){
    //         if((weights[i].idx != weights[j].idx)){
    //             if(weights[i].r == weights[j].r){
    //                 int s1 = weights[i].s;
    //                 int s2 = weights[j].s;
                    
                    
    //                 int w1 = database.cells[weights[i].idx].getCellSiteWidth();
    //                 int w2 = database.cells[weights[j].idx].getCellSiteWidth();

    //                 bool ovrlp = !(((s1+w1) < s2+eps) ||
    //                              ((s2+w2) < s1+eps));

    //                 if(ovrlp){
    //                     ovrlps.push_back({i,j});
    //                 }
 
    //             }
    //         }
    //     }//end loop 
    // }//end loop
    bool debug = false || debug_global;
    
    std::set<std::pair<int,int>> rs_loc;
    RTree cand_rtree;

    for(int i = 0 ; i < weights.size() ; i++){
        int s = weights[i].s;
        int r = weights[i].r;
        rs_loc.insert(std::make_pair(r,s));
        int w = database.cells[weights[i].idx].getCellSiteWidth();
        int h = 1;

        boostBox box(boostPoint(s+eps,r+eps),boostPoint(s+w-eps,r+h-eps));
        cand_rtree.insert({box, i});
    }//end for

    if(debug){
        std::stringstream ss;

        ss << "w_idx,inst,r,s"<<std::endl;

        for(auto cand : cand_rtree){
            const auto& b = cand.first;
            auto idx = cand.second;
            auto name = database.cells[weights[idx].idx].getName();
                
            auto box = utils::BoxT<DBU>(bg::get<bg::min_corner, 0>(b),
                                        bg::get<bg::min_corner, 1>(b),
                                        bg::get<bg::max_corner, 0>(b),
                                        bg::get<bg::max_corner, 1>(b));


            ss << std::to_string(idx)
            << "," << name 
            << "," << std::to_string(box.lx())
            << "," << std::to_string(box.ly())
            << std::endl;
        
        }


        std::string file_name_csv = db::setting.outputFile + "_legRtree.csv";
        std::ofstream fout(file_name_csv);
        fout << ss.str();
        fout.close();
    }
    

  


    std::set<std::string> encode_idx;
    std::set<std::pair<int,int>> conflict_set;

    for(int i = 0; i < weights.size() ; i++)
    {
        
        // conflict_set.insert(i);

        int r = weights[i].r;
        int s = weights[i].s;
        int cell_idx = weights[i].idx;
        int w = database.cells[cell_idx].getCellSiteWidth();
        std::string cell_name = database.cells[cell_idx].getName();


        double xl = s+eps;
        double yl = r+eps;
        double xh = s+w-eps;
        double yh = r+1-eps;

        boostBox rtreeQueryBox(
        boostPoint(xl,
                   yl),
        boostPoint(xh,
                   yh));
        vector<std::pair<boostBox, int>> queryResults;


        auto box_query = utils::BoxT<DBU>(bg::get<bg::min_corner, 0>(rtreeQueryBox),
                                    bg::get<bg::min_corner, 1>(rtreeQueryBox),
                                    bg::get<bg::max_corner, 0>(rtreeQueryBox),
                                    bg::get<bg::max_corner, 1>(rtreeQueryBox));
        if(debug){
            // log() << "i: " << i
            //   << ", cell_name: " << cell_name
            //   << ", xl: "    <<xl
            //   << ", yl: " << yl
            //   << ", xh: " << xh
            //   << ", yh: " << yh
            //   << ", eps: " << eps << std::endl;
        }
        

        // log() << "box: " << box_query << std::endl;
        
        cand_rtree.query(bgi::intersects(rtreeQueryBox), std::back_inserter(queryResults));

        
        // std::vector<int> ovrlp_set;
        std::set<int> ovrlp_set;
        


        // if(queryResults.size() <= 1){
        //     ovrlp_vec.push_back({});
        //     continue;
        // }

        for (const auto& queryResult : queryResults) {
            const auto& b = queryResult.first;
            auto j = queryResult.second;
            
            
            auto box = utils::BoxT<DBU>(bg::get<bg::min_corner, 0>(b),
                                        bg::get<bg::min_corner, 1>(b),
                                        bg::get<bg::max_corner, 0>(b),
                                        bg::get<bg::max_corner, 1>(b));
            // CellOvrlp cell_ovrlp;
            // // ovrlp_set.insert(idx);
            // if(ovrlp_set.find(weights[idx].idx) == ovrlp_set.end()){
            //     ovrlp_set.insert(weights[idx].idx);
            int r1 = weights[j].r;
            int s1 = weights[j].s;
            int cell_idx1 = weights[j].idx;
            int w1 = database.cells[cell_idx1].getCellSiteWidth();
            std::string cell_name1 = database.cells[cell_idx1].getName();

            if(debug){
                // log() 
                //     << "cell_name: " << cell_name1
                //     << ", r: " <<r1
                //     << ", s: " << s1
                //     << ", w: " << w1
                //     << ", j: " << j
                //     << std::endl;
            }

            int max = std::max(i,j);
            int min = std::min(i,j);
            
            if(conflict_set.find(std::make_pair(min,max)) == conflict_set.end()){
                std::vector<int > ovrlp_vec;
                if(cell_idx1 != cell_idx){
                    conflict_set.insert(std::make_pair(min,max));
                    ovrlp_vec.push_back(min);
                    ovrlp_vec.push_back(max);
                    ovrlps.push_back(ovrlp_vec);
                }//end if 
            }//end if 
                

            
            // }


            
        }//end for
        

        

        // if(debug){
        //     log() << "ovrlp idxs: "<< std::endl;
        //     std::string txt="";
        //     for(auto tmp : conflict_set){
        //         txt = txt + std::to_string(tmp) + "_";
        //     }
            

            
        //     log() <<txt <<  std::endl;
        // }

        // std::string encoded_txt = "";
        // for(auto tmp : conflict_set){
        //     encoded_txt = encoded_txt + std::to_string(tmp) + "_";
        // }
        // encoded_txt.pop_back();

        // if(encode_idx.find(encoded_txt) == encode_idx.end())
        // {
            
        //     encode_idx.insert(encoded_txt);
        // }
        


        // for(auto tmp : ovrlp_set)
        //     ovrlps.push_back(cell_ovrlp);
    }

    
    // for(auto code : encode_idx){
    //     vector<std::string> str_idxs;
    //     std::vector<int> ovrlp_vec;
    
    //     boost::split(str_idxs, code, boost::is_any_of("_"));

    //     if(debug)
    //         log() << "code: " << code << std::endl;

    //     for(auto str_idx : str_idxs){
    //         int idx = std::atoi(str_idx.c_str());
    //         // if(debug){
    //         int r = weights[idx].r;
    //         int s = weights[idx].s;
    //         int cell_idx = weights[idx].idx;
    //         int w = database.cells[cell_idx].getCellSiteWidth();
    //         std::string cell_name = database.cells[cell_idx].getName();
    //         ovrlp_vec.push_back(idx);

    //         if(debug){
    //             log() << ", cell_name: " << cell_name
    //             << ", r: "    <<r
    //             << ", s: " << s
    //             << ", w: " << w
    //             << "i: " << idx
    //             << std::endl;
    //         }
            
    //         // }
    //     }
    //     ovrlps.push_back(ovrlp_vec);
    // }
        

    // if(debug){
    //     log() << "ovrlps size: " << ovrlps.size() << std::endl;
    //     log() << "weights size: " << weights.size() << std::endl;
    // }
    
    
    
    

    // construct rtree
//  
// cell_rtree.insert({box, cell.idx});
    
}//end getLegalizerOvrlpConflictMatrix

void Legalizer::getLegalizerILPSinglePositionConflictMatrix(
      std::vector<std::vector<std::vector<int>>>& cube_cost
    , std::vector<std::vector<std::vector<int>>>& cube_weightIdx
    , std::vector<int>& movable_cells
    , std::vector<std::vector<int>>& conflict_matrix){
    // one site for each cell conflict
    
    for(int c = 0; c < cube_weightIdx.size() ; c++){
        auto cell_idx = movable_cells[c];
        auto cell = database.cells[cell_idx];
        std::vector<int> conflict_matrix_row;
        for(int r = 0; r < cube_weightIdx[c].size(); r++){
            for(int s = 0; s < cube_weightIdx[c][r].size(); s++){
                if(cube_cost[c][r][s] == -1 ) continue;
                conflict_matrix_row.push_back(cube_weightIdx[c][r][s]);
            }
        }
        conflict_matrix.push_back(conflict_matrix_row);
    }
}// end getLegalizerSinglePositionConflictMatrix


void Legalizer::getLegalizerILPSinglePositionConflictMatrixV3(
            std::vector<cellWrap>& weights
            , std::vector<int>& movable_cells
            , std::vector<std::vector<int>>& constraints){
    int num_cells = movable_cells.size();
    std::map<int,int> cellIdxToConstraintMap;

    for(int i = 0; i< movable_cells.size(); i++){
        cellIdxToConstraintMap[movable_cells[i]] = i;
    }



    constraints.resize(num_cells);

    for(int i = 0; i < weights.size(); i++){
        int cell_idx = weights[i].idx;
        constraints[cellIdxToConstraintMap[cell_idx]].push_back(i);
    }

}//end getLegalizerILPSinglePositionConflictMatrixV2

// void Legalizer::getLegalizerILPSinglePositionConflictMatrix(
//       std::vector<std::vector<std::vector<int>>>& cube_cost
//     , std::vector<std::vector<std::vector<int>>>& cube_weightIdx
//     , std::vector<int>& movable_cells
//     , std::vector<std::vector<int>>& conflict_matrix){
//     // one site for each cell conflict
    
//     for(int c = 0; c < cube_weightIdx.size() ; c++){
//         auto cell_idx = movable_cells[c];
//         auto cell = database.cells[cell_idx];
//         std::vector<int> conflict_matrix_row;
//         for(int r = 0; r < cube_weightIdx[c].size(); r++){
//             for(int s = 0; s < cube_weightIdx[c][r].size(); s++){
//                 if(cube_cost[c][r][s] == -1 ) continue;
//                 conflict_matrix_row.push_back(cube_weightIdx[c][r][s]);
//             }
//         }
//         conflict_matrix.push_back(conflict_matrix_row);
//     }
// }// end getLegalizerSinglePositionConflictMatrix


void Legalizer::getLegalizerSolutionILPSolverV2(
                  std::vector<cellWrap>& weights
                , std::vector<int>& legalize_rows
                , std::vector<int>& legalize_sites
                , std::vector<std::vector<int>>& ovrlps
                , std::vector<std::vector<int>>& constraints
                , std::vector<int>& sol){
    
    bool debug = false || debug_global;
    
    vector<int> selected;
    // std::vector<int> sol;
    selected.resize(weights.size());
    std::string log_name = "legalizer_ilp_"+database.cells[cell_idx_].getName()+".lp";
    bool maxmimize = false;
    
    IloEnv env;

    try {
        IloModel model(env);
        IloIntArray weightsSet(env);
        IloIntVarArray vars(env);
        IloRangeArray con(env);


        for (auto tmp : weights) {
            weightsSet.add(tmp.cost);
            // tmp_idx++;
        }

        int counter = 0;
        

        for (int i = 0 ; i < weights.size(); i++) {
            vars.add(IloIntVar(env, 0, 1));

            if(debug){
                int cell_idx = weights[i].idx; 
                int new_row  = legalize_rows[weights[i].r];
                int new_site = legalize_sites[weights[i].s];
                auto row_dbu = database.getDBURow(new_row)/database.libDBU;
                auto site_dbu = database.getDBUSite(new_site)/database.libDBU;



                std::string name = database.cells[cell_idx].getName() 
                                // + "_" + std::to_string(site_dbu)
                                // + "_" + std::to_string(row_dbu)
                                + "_" + std::to_string(weights[i].r)
                                + "_" + std::to_string(weights[i].s)
                                + "_" + std::to_string(i);
                vars[i].setName(name.c_str());
            }
            
        }

        // naming
        // // only one location for each cell constraint
        int j = 0;
        int counter_cons = 0;

        for (auto constraint_row : constraints) {
            IloNumExpr expr_constraint(env);
            for(auto var_idx : constraint_row){
                expr_constraint += vars[var_idx];    
            }  
            con.add(expr_constraint == 1);  
            if(debug){
                std::string const_name = "c_" + std::to_string(counter_cons);
                con[counter_cons].setName(const_name.c_str());
                counter_cons++;  
            }
            expr_constraint.end();    
            // if(j==100)break;
            // j++;
            // break;
        }


        for (auto conflict_row : ovrlps) {
            // if(j=!72){
            //     j++;
            //     continue;
            // }
            IloNumExpr expr_constraint(env);
            for(auto var_idx : conflict_row){
                expr_constraint += vars[var_idx];    
            }  
            con.add(expr_constraint <= 1);  
            if(debug){
                std::string const_name = "ovrlp_" + std::to_string(counter_cons);
                con[counter_cons].setName(const_name.c_str());
                counter_cons++;  
            }
            expr_constraint.end();    
            // if(j==100)break;
            // j++;
            // break;
        }


        // for (int i = 0 ; i < weights.size(); i++) {
        //     if(weights[i] != -1) continue;
        //     IloNumExpr expr_constraint(env);
            
        //     expr_constraint += vars[i];    
            
        //     con.add(expr_constraint == 0); 
        //     std::string const_name = "cons_" + std::to_string(counter_cons);
        //     con[counter_cons].setName(const_name.c_str());
        //     counter_cons++;  
        //     expr_constraint.end();    
        //     // if(j==100)break;
        //     // j++;
        //     // break;
        // }

        model.add(con);        

        IloNumExpr expression(env);
        for (int i = 0; i < weights.size(); i++) {
            expression += vars[i] * weightsSet[i];
        }
        if(!maxmimize){
            model.add(IloMinimize(env, expression));
        }else{
            model.add(IloMaximize(env, expression));
        }

        IloCplex cplex(model);
        if(!debug)
            cplex.setOut(env.getNullStream());


        cplex.setParam(IloCplex::Param::TimeLimit, db::setting.legalizerOptimzerTimeLimit);
        cplex.setParam(IloCplex::Param::Threads, db::setting.numThreads);
        if(db::setting.legalizationActivateSolutionPool){
            cplex.setParam(IloCplex::Param::MIP::Pool::Intensity, db::setting.legalizationSolutionPoolIntensity);
            cplex.setParam(IloCplex::Param::MIP::Limits::Populate, db::setting.legalizationMaxNumSolutionPool);
        }
        
        // cplex.setOut(env.getNullStream());
        
        // if(db::setting.debugLegalizerExportILPModel)
        if(debug)
            cplex.exportModel(log_name.c_str());
        // return;
        
        
        if(db::setting.legalizationActivateSolutionPool){
            if (!cplex.populate()) {
                log() << "Failed to optimize.\n";
                env.end();
                // std::exit(1);
                return;
            }
        }else{
            if (!cplex.solve()) {
                log() << "Failed to optimize.\n";
                env.end();
                // std::exit(1);
                return;
            }
        }

        if(db::setting.legalizationActivateSolutionPool){
            int numsol = cplex.getSolnPoolNsolns();
            log() << "numSol: " << numsol << std::endl;
        
        
            for (int j = 0; j < numsol; j++) {
                IloNumArray vals(env);
                cplex.getValues(vals, vars,j);
                // std::vector<int> sol;
                for (int i = 0; i < vals.getSize(); i++) {
                    try {
                        // log() << "vals[" << i << "]: " << vals[i] << std::endl;
                        if (vals[i] == 1) {
                            selected[i] = 1;
                            sol.push_back(i);
                        }
                    } catch (...) {
                        log() << "Exception on index " << i << "\n";
                    }
                }
                // log() << "sol size: " << sol.size() << std::endl;
                // if(db::setting.legalizationMaxNumSolutionPool != -1)
                //     if(j == db::setting.legalizationMaxNumSolutionPool)
                //         break;
            }//end numsol loop

        }else{
            IloNumArray vals(env);
        
            cplex.getValues(vals, vars);
            
            for (int i = 0; i < selected.size(); i++) {
                try {
                    if(debug)
                        log() << "vals[" << i << "]: " << vals[i] << std::endl;
                    if (vals[i] == 1) {
                        selected[i] = 1;
                        sol.push_back(i);
                    }
                } catch (...) {
                    log() << "Exception on index " << i << "\n";
                }
            }
        }//end if(db::setting.legalizationActivateSolutionPool)
        
    } catch (IloException& e) {
        log() << "Concert exception caught: " << e << endl;
    } catch (...) {
        log() << "Unknown exception caught" << endl;
    }
    env.end();
    // return selectedIndexes;

}//end getLegalizerSolutionILPSolverV2

void Legalizer::getLegalizerSolutionILPSolver(
      std::vector<double>& weights
    , std::vector<std::vector<int>>& conflict_matrix_ovrlp
    , std::vector<std::vector<int>>& conflict_matrix_single_position
    , std::unordered_map<int,int>& weightToCellIdxDict
    , std::unordered_map<int,std::tuple<int,int,int>>& weigthToCellRowSite
    , std::vector<int>& sol){
    //####################ILP Model####################################
    vector<int> selected;
    // std::vector<int> sol;
    selected.resize(weights.size());
    std::string log_name = "legalizer_ilp_"+database.cells[cell_idx_].getName()+".lp";
    bool maxmimize = false;
    
    IloEnv env;

    try {
        IloModel model(env);
        IloIntArray weightsSet(env);
        // IloNumArray weightsSet(env);
        IloIntVarArray vars(env);
        IloRangeArray con(env);
        int tmp_idx = 0;

        for (auto i : weights) {
            // log() << "idx: " << tmp_idx 
            //       << ", weight: " << i << std::endl;
            // weightsSet.add(std::max(i, 1));
            // if(i == -1) continue;
            weightsSet.add(i);
            tmp_idx++;
        }

        int counter = 0;
        // for (int i : selected) {
        //     vars.add(IloIntVar(env, 0, 1));
        //     std::string name = "x_" + std::to_string(counter);
        //     vars[i].setName(name.c_str());
        //     counter++;
        // }

        for (int i = 0 ; i < weights.size(); i++) {
            // if(weights[i]==-1) continue;
            vars.add(IloIntVar(env, 0, 1));
            int cell_idx = std::get<0>(weigthToCellRowSite[i]);
            int new_row  = std::get<1>(weigthToCellRowSite[i]) + row_offset_;
            int new_site = std::get<2>(weigthToCellRowSite[i]) + site_offset_;
            auto row_dbu = database.getDBURow(new_row)/database.libDBU;
            auto site_dbu = database.getDBUSite(new_site)/database.libDBU;

            // log() << "cell_idx: " << cell_idx
            //       << "new_row: " << new_row
            //       << "new_site: " << new_site
            //       << "row_dbu: " << row_dbu
            //       << "site_dbu: " << site_dbu << std::endl;

            std::string name = database.cells[weightToCellIdxDict[i]].getName() 
                               + "_" + std::to_string(site_dbu)
                               + "_" + std::to_string(row_dbu)
                               + "_" + std::to_string(i);
            vars[i].setName(name.c_str());
        }

        // naming
        // only one location for each cell constraint
        int j = 0;
        int counter_cons = 0;

        for (auto conflict_row : conflict_matrix_single_position) {
            IloNumExpr expr_constraint(env);
            for(auto var_idx : conflict_row){
                expr_constraint += vars[var_idx];    
            }  
            con.add(expr_constraint == 1);  
            std::string const_name = "c_" + std::to_string(counter_cons);
            con[counter_cons].setName(const_name.c_str());
            counter_cons++;  
            expr_constraint.end();    
            // if(j==100)break;
            // j++;
            // break;
        }


        for (auto conflict_row : conflict_matrix_ovrlp) {
            // if(j=!72){
            //     j++;
            //     continue;
            // }
            IloNumExpr expr_constraint(env);
            for(auto var_idx : conflict_row){
                expr_constraint += vars[var_idx];    
            }  
            con.add(expr_constraint <= 1);  
            std::string const_name = "ovrlp_" + std::to_string(counter_cons);
            con[counter_cons].setName(const_name.c_str());
            counter_cons++;  
            expr_constraint.end();    
            // if(j==100)break;
            // j++;
            // break;
        }


        for (int i = 0 ; i < weights.size(); i++) {
            if(weights[i] != -1) continue;
            IloNumExpr expr_constraint(env);
            
            expr_constraint += vars[i];    
            
            con.add(expr_constraint == 0); 
            std::string const_name = "cons_" + std::to_string(counter_cons);
            con[counter_cons].setName(const_name.c_str());
            counter_cons++;  
            expr_constraint.end();    
            // if(j==100)break;
            // j++;
            // break;
        }

        model.add(con);        

        IloNumExpr expression(env);
        for (int i = 0; i < weights.size(); i++) {
            if(weightsSet[i] == -1)
                expression += 0 * weightsSet[i];
            else
                expression += vars[i] * weightsSet[i];
        }
        if(!maxmimize){
            model.add(IloMinimize(env, expression));
        }else{
            model.add(IloMaximize(env, expression));
        }

        IloCplex cplex(model);
        cplex.setOut(env.getNullStream());
        cplex.setParam(IloCplex::Param::TimeLimit, db::setting.legalizerOptimzerTimeLimit);
        cplex.setParam(IloCplex::Param::Threads, db::setting.numThreads);
        
        // cplex.setOut(env.getNullStream());
        
        if(db::setting.debugLegalizerExportILPModel)
            cplex.exportModel(log_name.c_str());
        
        
        if (!cplex.solve()) {
            log() << "Failed to optimize.\n";
            env.end();
            // std::exit(1);
            return;
        }

        
        IloNumArray vals(env);
        
        cplex.getValues(vals, vars);
        
        for (int i = 0; i < selected.size(); i++) {
            try {
                // log() << "vals[" << i << "]: " << vals[i] << std::endl;
                if (vals[i] == 1) {
                    selected[i] = 1;
                    sol.push_back(i);
                }
            } catch (...) {
                log() << "Exception on index " << i << "\n";
            }
        }
    } catch (IloException& e) {
        log() << "Concert exception caught: " << e << endl;
    } catch (...) {
        log() << "Unknown exception caught" << endl;
    }
    env.end();
    // return selectedIndexes;


}//end getLegalizerSolutionILPSolver


void Legalizer::getLegalizerSolutionILPSolverSolutionPool(
      std::vector<double>& weights
    , std::vector<std::vector<int>>& conflict_matrix_ovrlp
    , std::vector<std::vector<int>>& conflict_matrix_single_position
    , std::unordered_map<int,int>& weightToCellIdxDict
    , std::vector<std::vector<int>>& sols){
    //####################ILP Model####################################
    vector<int> selected;
    // std::vector<int> sol;
    selected.resize(weights.size());
    std::string log_name = "legalizer_ilp_"+database.cells[cell_idx_].getName()+".lp";
    bool maxmimize = false;
    
    IloEnv env;
    try {
        IloModel model(env);
        IloIntArray weightsSet(env);
        // IloNumArray weightsSet(env);
        IloIntVarArray vars(env);
        IloRangeArray con(env);
        int tmp_idx = 0;
        for (auto i : weights) {
            // log() << "idx: " << tmp_idx 
            //       << ", weight: " << i << std::endl;
            // weightsSet.add(std::max(i, 1));
            // if(i == -1) continue;
            weightsSet.add(i);
            tmp_idx++;
        }

        int counter = 0;
        // for (int i : selected) {
        //     vars.add(IloIntVar(env, 0, 1));
        //     std::string name = "x_" + std::to_string(counter);
        //     vars[i].setName(name.c_str());
        //     counter++;
        // }
        for (int i = 0 ; i < weights.size(); i++) {
            // if(weights[i]==-1) continue;
            vars.add(IloIntVar(env, 0, 1));
            std::string name = database.cells[weightToCellIdxDict[i]].getName() + "_" + std::to_string(i);
            vars[i].setName(name.c_str());
        }

        // naming
        // only one location for each cell constraint
        int j = 0;
        int counter_cons = 0;
        for (auto conflict_row : conflict_matrix_single_position) {
            IloNumExpr expr_constraint(env);
            for(auto var_idx : conflict_row){
                expr_constraint += vars[var_idx];    
            }  
            con.add(expr_constraint == 1);  
            std::string const_name = "c_" + std::to_string(counter_cons);
            con[counter_cons].setName(const_name.c_str());
            counter_cons++;  
            expr_constraint.end();    
            // if(j==100)break;
            // j++;
            // break;
        }

        for (auto conflict_row : conflict_matrix_ovrlp) {
            // if(j=!72){
            //     j++;
            //     continue;
            // }
            IloNumExpr expr_constraint(env);
            for(auto var_idx : conflict_row){
                expr_constraint += vars[var_idx];    
            }  
            con.add(expr_constraint <= 1);  
            std::string const_name = "ovrlp_" + std::to_string(counter_cons);
            con[counter_cons].setName(const_name.c_str());
            counter_cons++;  
            expr_constraint.end();    
            // if(j==100)break;
            // j++;
            // break;
        }

        for (int i = 0 ; i < weights.size(); i++) {
            if(weights[i] != -1) continue;
            IloNumExpr expr_constraint(env);
            
            expr_constraint += vars[i];    
            
            con.add(expr_constraint == 0); 
            std::string const_name = "cons_" + std::to_string(counter_cons);
            con[counter_cons].setName(const_name.c_str());
            counter_cons++;  
            expr_constraint.end();    
            // if(j==100)break;
            // j++;
            // break;
        }

        model.add(con);        

        IloNumExpr expression(env);
        for (int i = 0; i < weights.size(); i++) {
            if(weightsSet[i] == -1)
                expression += 0 * weightsSet[i];
            else
                expression += vars[i] * weightsSet[i];
        }
        if(!maxmimize){
            model.add(IloMinimize(env, expression));
        }else{
            model.add(IloMaximize(env, expression));
        }

        IloCplex cplex(model);
        cplex.setOut(env.getNullStream());

        cplex.setParam(IloCplex::Param::Threads, db::setting.numThreads);
        cplex.setParam(IloCplex::Param::MIP::Pool::Intensity, db::setting.legalizationSolutionPoolIntensity);
        cplex.setParam(IloCplex::Param::MIP::Limits::Populate, db::setting.legalizationMaxNumSolutionPool);
        cplex.setParam(IloCplex::Param::TimeLimit, db::setting.legalizerOptimzerTimeLimit);
        
        if (!cplex.populate()) {
            log() << "Failed to optimize.\n";
            env.end();
            // std::exit(1);
            return;
        }

        int numsol = cplex.getSolnPoolNsolns();
        // log() << "numSol: " << numsol << std::endl;
        
        
        for (int j = 0; j < numsol; j++) {
            IloNumArray vals(env);
            cplex.getValues(vals, vars,j);
            std::vector<int> sol;
            for (int i = 0; i < vals.getSize(); i++) {
                try {
                    // log() << "vals[" << i << "]: " << vals[i] << std::endl;
                    if (vals[i] == 1) {
                        selected[i] = 1;
                        sol.push_back(i);
                    }
                } catch (...) {
                    log() << "Exception on index " << i << "\n";
                }
            }
            sols.push_back(sol);
            // if(db::setting.legalizationMaxNumSolutionPool != -1)
            //     if(j == db::setting.legalizationMaxNumSolutionPool)
            //         break;
        }//end numsol loop
    } catch (IloException& e) {
        log() << "Concert exception caught: " << e << endl;
    } catch (...) {
        log() << "Unknown exception caught" << endl;
    }
    env.end();
    // return selectedIndexes;


}//end getLegalizerSolutionILPSolver

void Legalizer::logLegalizerBlockageMatrix(std::vector<std::vector<int>>& blockage_matrix){
    // ##### Write blockage Matrix ######
    auto rsynService = database.getRsynService();
    // write instances in csv file
    std::stringstream ss;

    for(int r = 0; r < blockage_matrix.size();r++){
        for(int s = 0; s < blockage_matrix[r].size(); s++){
            ss << " ( r: "      << std::to_string(r)
               << ", s: "       << std::to_string(s)
               << ", v[r][s]: " << std::to_string(blockage_matrix[r][s]) << ") ," << std::endl;
            // log() << " ( r: "      << std::to_string(r)
            //       << ", s: "       << std::to_string(s)
            //       << ", v[r][s]: " << std::to_string(blockage_matrix[r][s]) << ") ," << std::endl;
            
        }
        ss << std::endl;
    }


    std::string file_name_csv = "blockage_matrix_"+database.cells[cell_idx_].getName()+".csv";
    std::ofstream fout(file_name_csv);
    fout << ss.str();
    fout.close();
}//end logLegalizerBlockageMatrix


void Legalizer::logCubeCost( std::vector<int>& movable_cells
                    , std::vector<std::vector<std::vector<int>>>& cube_cost
                    , std::vector<std::vector<std::vector<int>>>& cube_weightIdx){
    log() << "cube_cost: " << std::endl;
    for(int c = 0; c < cube_cost.size() ; c++){
        std::stringstream ss1;
        auto cell_idx = movable_cells[c];
        auto cell = database.cells[cell_idx];
        for(int r = 0; r < cube_cost[c].size(); r++){
            for(int s = 0; s < cube_cost[c][r].size(); s++){
                ss1 << std::to_string(cube_weightIdx[c][r][s])<< ":" << std::to_string(cube_cost[c][r][s])  <<", ";
            }
            ss1 << std::endl;
        }
        std::string file_name_csv1 = "cost_matrix_"+cell.getName()+".csv";
        std::ofstream fout1(file_name_csv1);
        fout1 << ss1.str();
        fout1.close();
    }   
}//end logCubeCost

void Legalizer::logCubeWeight( std::vector<int>& movable_cells,std::vector<std::vector<std::vector<int>>>& cube_weightIdx){
    log() << "cube_weightIdx: " << std::endl;
    for(int c = 0; c < cube_weightIdx.size() ; c++){
        std::stringstream ss1;
        auto cell_idx = movable_cells[c];
        auto cell = database.cells[cell_idx];
        for(int r = 0; r < cube_weightIdx[c].size(); r++){
            for(int s = 0; s < cube_weightIdx[c][r].size(); s++){
                ss1 << std::to_string(cube_weightIdx[c][r][s]) << ", ";
            }
            ss1 << std::endl;
        }
        std::string file_name_csv1 = "weightIdx_"+cell.getName()+".csv";
        std::ofstream fout1(file_name_csv1);
        fout1 << ss1.str();
        fout1.close();
    }
}//end logCubeWeight


void Legalizer::logLegalizerBoard(std::vector<std::vector<int>>& blockage_matrix
                                , std::vector<int>& legalize_rows
                                , std::vector<int>& legalize_sites){
    std::stringstream ss_blockage;
    // ss_blockage << "r,s,val,w,h" << std::endl;    
    ss_blockage << "xl,yl,xh,yh,cost" << std::endl;    
    for(int r = 0; r < legalize_rows.size() ; r++)
        for(int s = 0; s < legalize_sites.size() ; s++){
            int row_idx  = legalize_rows[r];
            int site_idx = legalize_sites[s];
            int row_dbu  = database.getDBURow(row_idx);
            int site_dbu = database.getDBUSite(site_idx);

            // if(debug){
            //     // log() << "r: " << r << ", s: " << s << std::endl;
            //     // log() << "row_idx: " << row_idx << ", site_idx: " << site_idx << std::endl;
            //     log() << "row_um: " <<  double(row_dbu)/double(database.libDBU)
            //         << ", site_um: " <<  double(site_dbu)/double(database.libDBU) << std::endl;
            // }
                

            auto xl = database.getDBUSite(legalize_sites[s]);
            auto yl = database.getDBURow(legalize_rows[r]);
            auto xh = xl + database.getSiteStep();
            auto yh = yl + database.getRowStep();
            

            ss_blockage << xl
                        << "," << yl
                        << "," << xh
                        << "," << yh
                        << "," << std::to_string(blockage_matrix[r][s])
                        << std::endl;
        }

    
    std::string file_name_csv_b = db::setting.directory +  db::setting.benchmarkName 
                                + ".legalizer.board.csv";
    if(debug_global_all){
        std::string cell_name = database.cells[cell_idx_].getName();
        file_name_csv_b = db::setting.directory +  db::setting.benchmarkName 
                                + ".legalizer.board."+cell_name+".csv";
    }
    
    std::ofstream fout_blockage(file_name_csv_b);
    fout_blockage << ss_blockage.str();
    fout_blockage.close();
}//end logLegalizerBoard

void Legalizer::logLegalizerWeights(std::vector<cellWrap>& weights
                                , std::vector<int>& legalize_rows
                                , std::vector<int>& legalize_sites){
    std::stringstream ss;
    auto rsynService = database.getRsynService();
       
    ss << "weight_idx,cell_name,xl,yl,xh,yh,cost" << std::endl;    

    for(int i = 0; i < weights.size(); i++){
        auto cell_wrap = weights[i];
        auto cell = database.cells[cell_wrap.idx];
        Rsyn::Cell cellRsyn = cell.rsynInstance.asCell();
        Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cellRsyn);
        Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cellRsyn);
        int cell_width_tmp = phLibCell.getWidth();
        int cell_height_tmp = phLibCell.getHeight();

        auto xl = database.getDBUSite(legalize_sites[cell_wrap.s]);
        auto yl = database.getDBURow(legalize_rows[cell_wrap.r]);
        auto xh = xl + cell_width_tmp;
        auto yh = yl + cell_height_tmp;

        ss << i
           << "," << cell.getName()
           << "," << xl
           << "," << yl
           << "," << xh
           << "," << yh
           << "," << cell_wrap.cost << std::endl;


    }//end for loop


    
    std::string file_name_csv_b = db::setting.directory +  db::setting.benchmarkName 
                                + ".legalizer.weights.csv";
    if(debug_global_all){
        std::string cell_name = database.cells[cell_idx_].getName();
        file_name_csv_b = db::setting.directory +  db::setting.benchmarkName 
                                + ".legalizer.weights."+cell_name+".csv";
    }

    std::ofstream fout_blockage(file_name_csv_b);
    fout_blockage << ss.str();
    fout_blockage.close();
}//end logLegalizerBoard


void Legalizer::logLegalizerOvrlps(std::vector<cellWrap>& weights
                                  , std::vector<int>& legalize_rows
                                  , std::vector<int>& legalize_sites
                                  ,std::vector<std::vector<int>>& ovrlps){
    std::stringstream ss;
    auto rsynService = database.getRsynService();
       
    ss << "weight_idx1,cell_name1,xl1,yl1,xh1,yh1,weight_idx2,cell_name2,xl2,yl2,xh2,yh2" << std::endl;    

    auto getWeightString = [&](int i){
        std::string str="";
        auto cell_wrap = weights[i];
        auto cell = database.cells[cell_wrap.idx];
        Rsyn::Cell cellRsyn = cell.rsynInstance.asCell();
        Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cellRsyn);
        Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cellRsyn);
        int cell_width_tmp = phLibCell.getWidth();
        int cell_height_tmp = phLibCell.getHeight();

        auto xl = database.getDBUSite(legalize_sites[cell_wrap.s]);
        auto yl = database.getDBURow(legalize_rows[cell_wrap.r]);
        auto xh = xl + cell_width_tmp;
        auto yh = yl + cell_height_tmp;

        str = std::to_string(i)
           + "," + cell.getName()
           + "," + std::to_string(xl)
           + "," + std::to_string(yl)
           + "," + std::to_string(xh)
           + "," + std::to_string(yh);
        return str;
    };

    for(auto ovrlp_row : ovrlps){
        if(ovrlp_row.size() == 2){
            std::string txt = "";
            bool first = true;
            for(auto ovrlp : ovrlp_row){
                if(first){
                    txt = getWeightString(ovrlp);
                    first = false;
                }else{
                    txt = txt + "," + getWeightString(ovrlp);
                }            
            }//end loop
            ss << txt << std::endl;
        }else{
            log() << "ovrlp size must be 2, wrong logging" << std::endl;
            std::exit(1);
        }//end if
        

    }//end loop 

    
    std::string file_name_csv_b = db::setting.directory +  db::setting.benchmarkName 
                                + ".legalizer.overlaps.csv";
    if(debug_global_all){
        std::string cell_name = database.cells[cell_idx_].getName();
        file_name_csv_b = db::setting.directory +  db::setting.benchmarkName 
                                + ".legalizer.overlaps."+cell_name+".csv";
    }                            
    std::ofstream fout_blockage(file_name_csv_b);
    fout_blockage << ss.str();
    fout_blockage.close();

}//end logLegalizerOvrlps


void Legalizer::logLegalizerSols(std::vector<cellWrap>& weights
                                , std::vector<int>& legalize_rows
                                , std::vector<int>& legalize_sites
                                , std::vector<int>& sols
                                ){
    std::stringstream ss;
    auto rsynService = database.getRsynService();
       
    ss << "weight_idx,cell_name,xl,yl,xh,yh" << std::endl;    

    auto getWeightString = [&](int i){
        std::string str="";
        auto cell_wrap = weights[i];
        auto cell = database.cells[cell_wrap.idx];
        Rsyn::Cell cellRsyn = cell.rsynInstance.asCell();
        Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cellRsyn);
        Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cellRsyn);
        int cell_width_tmp = phLibCell.getWidth();
        int cell_height_tmp = phLibCell.getHeight();

        auto xl = database.getDBUSite(legalize_sites[cell_wrap.s]);
        auto yl = database.getDBURow(legalize_rows[cell_wrap.r]);
        auto xh = xl + cell_width_tmp;
        auto yh = yl + cell_height_tmp;

        str = std::to_string(i)
           + "," + cell.getName()
           + "," + std::to_string(xl)
           + "," + std::to_string(yl)
           + "," + std::to_string(xh)
           + "," + std::to_string(yh);
        return str;
    };

    for(auto sol : sols){
        ss << getWeightString(sol) << std::endl;
    }

    
    std::string file_name_csv_b = db::setting.directory +  db::setting.benchmarkName 
                                + ".legalizer.solution.csv";

    if(debug_global_all){
        std::string cell_name = database.cells[cell_idx_].getName();
        file_name_csv_b = db::setting.directory +  db::setting.benchmarkName 
                                + ".legalizer.solution."+cell_name+".csv";
    }   

    std::ofstream fout_blockage(file_name_csv_b);
    fout_blockage << ss.str();
    fout_blockage.close();

}//end logLegalizerSols


void Legalizer::constructSegments(std::vector<std::vector<int>>& mat_spr
                          , std::vector<std::vector<std::pair<int,int>>>& segs){

    for(int r = 0; r < mat_spr.size() ; r++){
        std::vector<int> seqs;
        for(int s = 0; s < mat_spr[r].size(); s++){
            seqs.push_back(mat_spr[r][s]);
            if(s+1 < mat_spr[r].size()){
                if(std::abs(mat_spr[r][s]-mat_spr[r][s+1]) != 1){
                    int min = *std::min_element(seqs.begin(),seqs.end());
                    int max = *std::max_element(seqs.begin(),seqs.end());
                    segs[r].push_back(std::make_pair(min,max));
                    seqs.clear();
                }//end if 
            }else{
                int min = *std::min_element(seqs.begin(),seqs.end());
                int max = *std::max_element(seqs.begin(),seqs.end());
                segs[r].push_back(std::make_pair(min,max));
                seqs.clear();

            }//end if 
        }//end for 
    }//end for 

}// constructSegments

void Legalizer::assignMovableCellsToOrders( std::vector<int>& movable_cells
                                   , std::vector<int>& orders
                                   , std::vector<int>& cell_w){
    for(int i = 0; i < movable_cells.size(); i++){
        
        orders.push_back(i);
        auto cell = database.cells[movable_cells[i]];
        cell_w.push_back(cell.getCellSiteWidth());
        // if(debug){
        //     log() << "cell: " << cell.getName()
        //       << ", i: " << i
        //       << ", w: " << cell_w[cell_w.size() -1 ] << std::endl;
        // }
        
    }
}//assignMovableCellsToOrders

void Legalizer::initSegments(std::vector<std::vector<std::pair<int,int>>>& segs
                , std::vector<std::vector<int>>& orders_2d
                , std::vector<int>& cell_w
                , std::vector<Segment>& segments){
    for(int r = 0; r < segs.size(); r++){
        for(int i = 0; i < segs[r].size(); i++){
            // get orders 
            int start = segs[r][i].first;
            int stop = segs[r][i].second;
            Segment seg_obj;
            seg_obj.r = r;
            seg_obj.start = start;
            seg_obj.stop = stop;
            
            for(int ii = 0; ii < orders_2d.size(); ii++){
                int sum = 0;
                // std::string txt = "";
                for(int jj = 0; jj < orders_2d[ii].size(); jj++){
                    // if(debug){
                    //     txt = txt + " (ord: "+std::to_string(orders_2d[ii][jj]) 
                    //         + ", w: " + std::to_string(cell_w[orders_2d[ii][jj]])
                    //         + ") " ;
                    // }
                    sum += cell_w[orders_2d[ii][jj]];
                }//end for
                // if(debug){
                //     log() << "order: " << txt << std::endl; 
                //     log()  << "r: " << r
                //         << ", start: " << seg_obj.start
                //         << ", stop: " << seg_obj.stop
                //         << ", sum: " << sum << std::endl;
                // }
                
                if(sum <= std::abs(start-stop)+1){
                    seg_obj.orders.push_back(orders_2d[ii]);
                }
            }//end for

            segments.push_back(seg_obj);
        }//end for 
    }//end for

}//initSegments

void Legalizer::assignCellsToSegments(std::vector<Segment>& segments
                        , std::vector<int>& movable_cells
                        , std::vector<int>& legalize_rows
                        , std::vector<int>& legalize_sites
                        , std::vector<int>& cell_w
                        , std::vector<cellWrap>& weights
                        , std::stringstream& ss){
    auto printSegs = [](int r, int seg_w, std::vector<std::vector<int>>& seg_orders_2d){
        log() << "seg.r: " << r
                << ", seg.w: " << seg_w << std::endl;
        for(auto ords : seg_orders_2d){
            std::string txt = "";
            for(auto elem : ords ){
                txt = txt + std::to_string(elem) + " ";
            }
            log() << txt << std::endl;
        }//end for 
    };


    bool debug = true || debug_global;    
    auto rsynService = database.getRsynService();
    
    for(auto seg : segments){
        int r = seg.r;
        auto seg_orders_2d = seg.orders;
        int seg_w = std::abs(seg.start - seg.stop)+1;
        if(debug){
            printSegs(r,seg_w,seg_orders_2d);
        }//end if 

        
        

        for(auto order : seg_orders_2d){
            // order idx -> r, s
            std::unordered_map<int,std::pair<int,int>> weights_dict;

            initWeightIdx(weights_dict,order);
            initAssignment(weights_dict,cell_w,order,seg);

            int cells_w_sum = 0;
            std::string txt = "";

            // get the total sum of cells' width
            for(int idx : order){
                txt = txt + std::to_string(idx) + " ";
                cells_w_sum += cell_w[idx];
            }

            if(debug){
                log() << "order: " << txt << std::endl;
                log() << "cell_w_sum: " << cells_w_sum << std::endl;
            }
                
            int empty_space = seg_w-cells_w_sum;

            if(debug)
                log() << "empty_space: " << empty_space << std::endl;

            insertGapBetweenCells(weights_dict
                                  , empty_space
                                  , cell_w
                                  , order
                                  , seg);

            initCost( weights_dict
                    , empty_space
                    , movable_cells
                    , cell_w
                    , order
                    , legalize_rows
                    , legalize_sites
                    , seg
                    , weights);

        }//end order loop

    }//end segments loop                            

}//assignCellsToSegments

void Legalizer::initWeightIdx(std::unordered_map<int,std::pair<int,int>>& weight_dict, std::vector<int>& order){
    for(int i : order){
        weight_dict[i] = std::make_pair(-1,-1);
    }//end for
}//end initWeightIdx
void Legalizer::initAssignment(std::unordered_map<int,std::pair<int,int>>& weight_dict
                        ,std::vector<int>& cell_w
                        ,std::vector<int>& order
                        ,Segment& seg){
    int idx = seg.start;
    for(int i : order){
        weight_dict[i]=std::make_pair(idx,idx+cell_w[i]-1);
        idx = idx + cell_w[i];
    }//end loop                        
}//end initAssignment
int Legalizer::getRandomInt(int rng_val){
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> rand_dist(0,rng_val); // distribution in range [1, 6]
    return rand_dist(rng);
}//end getRandomInt
void Legalizer::insertGapBetweenCells(std::unordered_map<int,std::pair<int,int>>& weight_dict
                                , int empty_space
                                , std::vector<int>& cell_w
                                , std::vector<int>& order
                                , Segment& seg){
    if(empty_space <= 0)
        return;
    
    int reminder = empty_space;
    int idx = seg.start;
    for(int i : order){
        auto& weight_tmp = weight_dict[i];
        int rnd = getRandomInt(reminder);
        weight_tmp = std::make_pair(idx+rnd,idx+cell_w[i]+rnd-1);
        idx = weight_tmp.second+1;
        reminder = reminder - rnd;
    }//end loop
}//end insertGapBetweenCells
void Legalizer::initCost(std::unordered_map<int,std::pair<int,int>>& weight_dict
                , int empty_space
                , std::vector<int>& movable_cells
                , std::vector<int>& cell_w
                , std::vector<int>& order
                , std::vector<int>& legalize_rows
                , std::vector<int>& legalize_sites
                , Segment& seg
                , std::vector<cellWrap>& weights){
    for(int i : order){
        cellWrap cell_wrap;
        auto& weight_tmp = weight_dict[i];
        auto cell = database.cells[movable_cells[i]];
        cell_wrap.idx = cell.idx;
        cell_wrap.r = seg.r;
        cell_wrap.s = weight_tmp.first;

        // find cost 
        int row_idx  = legalize_rows[cell_wrap.r];
        int site_idx = legalize_sites[cell_wrap.s];
                            
        // auto median = cell.feature.median_x_y_pin;//getMedianPin();
        
        
        // double cost_median = diff_x + diff_y;
        double cost_median = getCandidateCost(cell,row_idx,site_idx);
        if(cost_median > 0)
            cell_wrap.cost = cost_median;
        else
            cell_wrap.cost = 1;

            weights.push_back(cell_wrap);

    }//end loop
}//end initCost

}//end namespace db



// auto printSegs = [](int r, int seg_w, std::vector<std::vector<int>>& seg_orders_2d){
//         log() << "seg.r: " << r
//                 << ", seg.w: " << seg_w << std::endl;
//         for(auto ords : seg_orders_2d){
//             std::string txt = "";
//             for(auto elem : ords ){
//                 txt = txt + std::to_string(elem) + " ";
//             }
//             log() << txt << std::endl;
//         }//end for 
//     };


//     auto initAssignment = [](){

//     };

//     bool debug = true || debug_global;    
//     auto rsynService = database.getRsynService();
    
//     for(auto seg : segments){
//         int r = seg.r;
//         auto seg_orders_2d = seg.orders;
//         int seg_w = std::abs(seg.start - seg.stop)+1;
//         if(debug){
//             printSegs(r,seg_w,seg_orders_2d);
//         }//end if 

        
        

//         for(auto order : seg_orders_2d){
//             // int cell idx -> r, s
//             std::unordered_map<int,std::pair<int,int>> weights_dict;
            
//             int cells_w_sum = 0;
//             std::string txt = "";

//             // get the total sum of cells' width
//             for(int idx : order){
//                 txt = txt + std::to_string(idx) + " ";
//                 cells_w_sum += cell_w[idx];
//             }

//             if(debug){
//                 log() << "order: " << txt << std::endl;
//                 log() << "cell_w_sum: " << cells_w_sum << std::endl;
//             }
                
//             int empty_w = seg_w-cells_w_sum;

//             if(debug)
//                 log() << "empty_w: " << empty_w << std::endl;


//             int num_line_seg = order.size();
//             int num_empty_spot = num_line_seg + 1;
//             int empty_offset = empty_w / num_empty_spot;
//             int reminder = empty_w % num_empty_spot;
//             int n = num_line_seg + num_empty_spot;
//             int idx = seg.start;
//             bool emptySw = true;
//             int i_order = 0;

//             if(debug){
//                 log() << "empty_offset: " << empty_offset
//                       << ", reminder: " << reminder
//                       << ", n: " << n << std::endl;
//             }

            

//             for(int i = 0; i < n; i++){
//                 if(idx > seg.stop)
//                     continue;

//                 if(debug){
//                     log() << "idx: " << idx << std::endl;
//                 }
//                 if(emptySw){
//                 // if(false){
//                     if(reminder > 0){
//                         std::random_device dev;
//                         std::mt19937 rng(dev());
//                         std::uniform_int_distribution<std::mt19937::result_type> rand_dist(0,reminder); // distribution in range [1, 6]
//                         int rnd_por = rand_dist(rng);
                        
//                         if(debug){
//                             log() << "rnd_por: " << rnd_por
//                                   << ", reminder: " << reminder 
//                                   << ", emptySw: " << emptySw << std::endl;
//                         }
                        
                        
//                         reminder = reminder - rnd_por;
//                     }

//                     idx = idx + empty_offset + reminder;
//                 }else{

                    
//                     cellWrap cell_wrap;

//                     // log() << "i_order: " << i_order << ", orer: " << order[i_order]
//                     //       << ", movable_cells[order[i_order]]: " << movable_cells[order[i_order]]
//                     //       << std::endl;
//                     // log() << "movable_cells size: " << movable_cells.size() << std::endl;

//                     // auto cell = database.cells[movable_cells[order[0]]];
//                     auto cell = database.cells[movable_cells[order[i_order]]];
//                     // continue;
//                     cell_wrap.idx = cell.idx;
//                     cell_wrap.r = seg.r;
//                     cell_wrap.s = idx;

//                     // find cost 
//                     int row_idx  = legalize_rows[cell_wrap.r];
//                     int site_idx = legalize_sites[cell_wrap.s];
                                        
//                     // auto median = cell.feature.median_x_y_pin;//getMedianPin();
                    
                    
//                     // double cost_median = diff_x + diff_y;
//                     double cost_median = getCandidateCost(cell,row_idx,site_idx);
//                     if(cost_median > 0)
//                         cell_wrap.cost = cost_median;
//                     else
//                         cell_wrap.cost = 1;


                    

                    
//                     if(debug){
                            
//                             Rsyn::Cell cellRsyn = cell.rsynInstance.asCell();
//                             Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cellRsyn);
//                             Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cellRsyn);
//                             int cell_width_tmp = phLibCell.getWidth();
//                             int cell_height_tmp = phLibCell.getHeight();

//                             log() << cell.getName() 
//                                 << ",r: " << std::to_string(cell_wrap.r)
//                                 << ",s: " << std::to_string(cell_wrap.s)
//                                 << ",w: " << std::to_string(cell_width_tmp)
//                                 << ",h: " << std::to_string(cell_height_tmp)
//                                 << ",cost: " << std::to_string(cost_median) << std::endl;

//                             // ss << cell.getName() 
//                             //     << "," << std::to_string(cell_wrap.r)
//                             //     << "," << std::to_string(cell_wrap.s)
//                             //     << "," << std::to_string(cell_width_tmp)
//                             //     << "," << std::to_string(cell_height_tmp)
//                             //     << "," << std::to_string(cost_median) << std::endl;
//                             auto xl = database.getDBUSite(legalize_sites[cell_wrap.s]);
//                             auto yl = database.getDBURow(legalize_rows[cell_wrap.r]);
//                             auto xh = xl + cell_width_tmp;
//                             auto yh = yl + cell_height_tmp;
//                             ss  << cell.getName() 
//                                 << "," << xl
//                                 << "," << yl
//                                 << "," << xh
//                                 << "," << yh
//                                 << "," << std::to_string(cost_median)
//                                 << std::endl;
//                         }
                    
//                     weights.push_back(cell_wrap);
//                     idx = idx + cell_w[order[i_order]];
//                     i_order++;
//                 }

//                 emptySw = emptySw xor true;

//                 if(debug){
//                     log() << std::endl;
//                 }
                
//             }//end for 



//         }//end order loop

//     }//end segments loop                    