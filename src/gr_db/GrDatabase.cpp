#include "GrDatabase.h"
#include <fstream>

gr::GrDatabase grDatabase;

namespace gr {
void GrDatabase::init() {
    if(db::setting.debug){
        log() << "this: " << this << std::endl;
    }
    GrRouteGrid::init();
    
    GrNetlist::init(*this);
    if(db::setting.debug){
        log() << "after GrRouteGrid: " << this << std::endl;
        GrRouteGrid::print();
    }
}

void GrDatabase::update() {
    // GrRouteGrid::clear();
    GrRouteGrid::update();
    // log() << "this: " << this << std::endl;
    GrNetlist::update(*this);

    // log() << "*this pointer " << this << std::endl;
    if(db::setting.debug){
        GrRouteGrid::print();
    }

    // setUnitViaMultiplier(1);
            
    // setLogisticSlope(1);
}


void GrDatabase::writeGuides(std::string filename) {
    log() << "Writing guides to file..." << std::endl;

    std::stringstream ss;

    auto printGrGuides = [&](const vector<GrBoxOnLayer>& guides) {
        for (const auto& guide : guides) {
            // log() << "guide: " << guide << std::endl;
            // ss << getCoor(guide[X].low, X) << " ";
            // ss << getCoor(guide[Y].low, Y) << " ";
            // ss << getCoor(guide[X].high + 1, X) << " ";
            // ss << getCoor(guide[Y].high + 1, Y) << " ";
            // ss << database.getLayer(guide.layerIdx).name << std::endl;
            ss << getCoor(guide[X].low, X) << " ";
            ss << getCoor(guide[Y].low, Y) << " ";
            ss << getCoor(guide[X].high + 1, X) << " ";
            ss << getCoor(guide[Y].high + 1, Y) << " ";
            ss << database.getLayer(guide.layerIdx).name << std::endl;
        }
    };

    for (const auto& net : grDatabase.nets) {
        // log() << "net: " << net.getName() << std::endl;
        ss << net.getName() << std::endl;
        ss << "(" << std::endl;
        // log() << "routeGuide " << std::endl;
        printGrGuides(net.wireRouteGuides);
        // log() << "viaGuide " << std::endl;
        printGrGuides(net.viaRouteGuides);
        // log() << "patchGuide     " << std::endl;
        printGrGuides(net.patchRouteGuides);
        ss << ")" << std::endl;
    }

    std::ofstream fout(filename);
    fout << ss.str();
    fout.close();
}

void GrDatabase::writeGuidesCSV(std::string filename) {
    log() << "Writing guides to file..." << std::endl;

    std::stringstream ss;
    ss << "net,box" << std::endl;

    auto printGrGuides = [&](const vector<GrBoxOnLayer>& guides) {
        for (const auto& guide : guides) {
            ss << getCoor(guide[X].low, X) << "_";
            ss << getCoor(guide[Y].low, Y) << "_";
            ss << getCoor(guide[X].high + 1, X) << "_";
            ss << getCoor(guide[Y].high + 1, Y) << "_";
            ss << database.getLayer(guide.layerIdx).name << "|";
        }
    };

    for (const auto& net : grDatabase.nets) {
        ss << net.getName();
        ss << ", \"";
        printGrGuides(net.wireRouteGuides);
        printGrGuides(net.viaRouteGuides);
        printGrGuides(net.patchRouteGuides);
        ss << "\""<<std::endl;
    }
    std::string file_name_csv = filename + "_guides.csv";
    std::ofstream fout(file_name_csv);
    fout << ss.str();
    fout.close();
}

void GrDatabase::logGCellGrid(){
    
    
    std::string file_name = db::setting.directory +  db::setting.benchmarkName+ ".gcell.csv";
    
    std::ofstream file(file_name);
    std::stringstream stream;
    stream << "l,x,y,w,h,dir" << std::endl;

    for (int layerIdx = 0; (layerIdx + 1) < database.getLayerNum(); ++layerIdx) {
       Dimension dir = database.getLayerDir(layerIdx);
        for (int x = 0; x < getNumGrPoint(X); x++) {
            for (int y = 0; y < getNumGrPoint(Y); y++) {
                auto grPoint = GrPoint({layerIdx,x,y});
                stream << std::to_string(layerIdx)
                   << "," <<  getCoorIntvl(grPoint,X).low
                   << "," <<  getCoorIntvl(grPoint,Y).low
                   << "," <<  getCoorIntvl(grPoint,X).high
                   << "," <<  getCoorIntvl(grPoint,Y).high
                   << "," <<  dir
                   << std::endl;
            }
        }
    }
    file << stream.str();
    file.close();
}//end logCellLocations

void GrDatabase::reportGR(std::string filename) {
    log() << "Writing GR report to file..." << std::endl;

    std::stringstream ss;
    ss << "net,hpwl,instTerms,wl,path_cost,time_gr" << std::endl;

    // auto printGrGuides = [&](const vector<GrBoxOnLayer>& guides) {
    //     for (const auto& guide : guides) {
    //         ss << getCoor(guide[X].low, X) << "_";
    //         ss << getCoor(guide[Y].low, Y) << "_";
    //         ss << getCoor(guide[X].high + 1, X) << "_";
    //         ss << getCoor(guide[Y].high + 1, Y) << "_";
    //         ss << database.getLayer(guide.layerIdx).name << "|";
    //     }
    // };

    for (auto& net : grDatabase.nets) {
        net.calcHPWL();
        ss << net.getName();
        ss << "," << net.getHPWL();
        ss << "," << net.dbNet.numOfPins();
        ss << "," << net.getWirelength();
        ss << "," << net.getPathCost();
        ss << "," << net.dbNet.getTime();
        ss << std::endl;
    }
    std::string file_name_csv = filename + ".gr.report.csv";
    std::ofstream fout(file_name_csv);
    fout << ss.str();
    fout.close();
}

void GrDatabase::reportCells(std::string filename){
    log() << "Writing GR report to file..." << std::endl;
    

    std::stringstream ss;
    ss << "name,x,y,w,h,degree,logic,connected_cells,connected_nets" << std::endl;

  

    auto rsynService = database.getRsynService();

    for (auto& cellWrapper : database.cells) {
        auto& instance = cellWrapper.rsynInstance;
        auto cellWrapper_idx = cellWrapper.idx;
        Rsyn::Cell cell = instance.asCell();
        Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cell);
        Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cell);


        ss << instance.getName();
        ss << "," << instance.getX();
        ss << "," << instance.getY();
        ss << "," << phLibCell.getWidth();
        ss << "," << phLibCell.getHeight();
        ss << "," << cell.getLibraryCellName();


        std::vector<db::Net> nets;
        std::vector<int> nets_idx;
        std::set<int> connectd_cells;
        
        cellWrapper.getCellNetsIdx(nets_idx);

        for(auto idx : nets_idx){
            nets.push_back(database.nets[idx]);
        }

        for(auto net : nets){
            std::vector<int> cells_idx;
            cellWrapper.getNetCells(net.idx,cells_idx);
            for(int cell_idx : cells_idx){
                if(cell_idx != cellWrapper_idx)
                    connectd_cells.insert(cell_idx);
            }
        }

        ss << "," ;     
        ss << connectd_cells.size();

        ss << "," ;     
        for(auto cell_idx : connectd_cells){
            ss << database.cells[cell_idx].getName() << "_";
        }


        ss << "," ;     
        for(auto idx : nets_idx){
            ss << database.nets[idx].getName() << "_";
        }
        ss << std::endl;     




    }
    std::string file_name_csv = filename + ".cells.report.csv";
    std::ofstream fout(file_name_csv);
    fout << ss.str();
    fout.close();


}//end report Cells


void GrDatabase::reportRL(std::string filename){
    
    std::stringstream ss;
    ss << "GGridBoundaryIdx " << "0 0 " \
       << getNumGrPoint(X) << " "
       << getNumGrPoint(Y) << " " << std::endl;

    ss << "NumLayer " << database.getLayerNum() << std::endl;

    for (int layerIdx = 0; layerIdx < database.getLayerNum(); ++layerIdx) {
        double avg_num_tracks = 0;
        int i = 0;
        for (int x = 0; x < getNumGrPoint(X); x++) {
            for (int y = 0; y < getNumGrPoint(Y); y++) {
                avg_num_tracks += getInCellArea({layerIdx, x, y});
                i++;
            }
        }

        avg_num_tracks = round(avg_num_tracks/i);
        Dimension dir = database.getLayerDir(layerIdx);

        if(dir == 1){
            ss << "Lay M" << std::to_string(layerIdx) 
               << " H " << avg_num_tracks << std::endl;
        }else{
            ss << "Lay M" << std::to_string(layerIdx) 
               << " V " << avg_num_tracks << std::endl;
        }
    }

    int i = 0;
    for(auto net : grDatabase.nets){
        // if(net.numOfPins() == 2){
            i++;
        // }
    }

    ss << "NumNets " << i << std::endl;

    for(auto net : grDatabase.nets){
        // if(net.numOfPins() == 2){
            ss << net.getName() 
               << " " << net.numOfPins() << std::endl;

            for(auto pins_pt : net.pinAccessBoxes){
                for(auto pt : pins_pt){
                    ss << pt.getPrefIdx() << " ";
                }
                ss << std::endl;
            }
        // }
    }

    ss << "routing " << std::endl;

    for(auto net : grDatabase.nets){
        ss << net.getName() << std::endl;
        // Read wire guides and via guides
        const auto& guides = net.wireRouteGuides;
        for (const auto& guide : guides) {
            ss << guide.lx() 
            << ", " << guide.ly() 
            << ", " << guide.hx()+1 
            << ", " << guide.hy()+1
            << ", " << guide.layerIdx << std::endl; 
        }



        

        const auto& viaGuides = net.viaRouteGuides;
        for (const auto& guide : viaGuides) {
            ss << guide.lx() 
            << ", " << guide.ly() 
            << ", " << guide.hx()+1 
            << ", " << guide.hy()+1
            << ", " << guide.layerIdx << std::endl; 
        }


    }



    std::string file_name_csv = filename + ".RLoutput.txt";
    std::ofstream fout(file_name_csv);
    fout << ss.str();
    fout.close();


}//end reportRL


}  // namespace gr
