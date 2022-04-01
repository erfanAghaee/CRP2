#pragma once

#include "db/Database.h"
#include "GrRouteGrid.h"
#include "global.h"
#include "GrNet.h"

namespace gr {
class GrDatabase : public GrRouteGrid, public GrRouteGrid2D, public GrNetlist {
public:
    void init();
    void update();
    void writeGuides(std::string filename);
    void writeGuidesCSV(std::string filename);
    void reportGR(std::string filename);
    void reportCells(std::string filename);
    void reportRL(std::string filename);

    
    void logGCellGrid();
    void logNets(int iter);
    void logVio(int iter);

    // std::set<int> critical_cells;
private:
};

}  // namespace gr

extern gr::GrDatabase grDatabase;