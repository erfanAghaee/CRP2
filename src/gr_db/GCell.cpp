#include "GCell.h"

namespace gr {
void GCellGrid::clear(){
    numTracks.clear();
    grid.clear();
}//end clear

void GCellGrid::init() {
    db::RsynService rsynService;
    clear();
    rsynService.init();
    const Rsyn::Session session;
    Rsyn::PhysicalDesign physicalDesign =
        static_cast<Rsyn::PhysicalService *>(session.getService("rsyn.physical"))->getPhysicalDesign();
    const DBU libDBU = physicalDesign.getDatabaseUnits(Rsyn::LIBRARY_DBU);

    grid.resize(2);
    grid[0].push_back(database.dieRegion[X].low);
    grid[1].push_back(database.dieRegion[Y].low);


    if (physicalDesign.allPhysicalGCell().empty()) {
        for (unsigned direction = 0; direction != 2; ++direction) {
             const int lo = database.dieRegion[direction].low;
             const int hi = database.dieRegion[direction].high;
             for (int i = lo + db::setting.gcellSize; i + db::setting.gcellSize < hi; i += db::setting.gcellSize) {
                 grid[direction].push_back(i);
             }
         }
    } else {
        for (const Rsyn::PhysicalGCell &rsynGCell : physicalDesign.allPhysicalGCell()) {
            const int location = rsynGCell.getLocation();
            const int step = rsynGCell.getStep();
            const int numStep = rsynGCell.getNumTracks();
            const Dimension direction = rsynGCell.getDirection() == Rsyn::PhysicalGCellDirection::VERTICAL ? X : Y;

            for (int i = 1; i < numStep; ++i) {
                grid[direction].push_back(location + step * i);
            }
        }
    }

    for (unsigned direction = 0; direction != 2; ++direction) {
        sort(grid[direction].begin(), grid[direction].end());
        if (grid[direction].back() < database.dieRegion[direction].high) {
            grid[direction].push_back(database.dieRegion[direction].high);
        } else if (grid[direction].back() > database.dieRegion[direction].high) {
            log() << "Warning: grid line " << grid[direction].back() << " exceeds die region " << database.dieRegion[direction].high << std::endl;
        }
    }

    numTracks.resize(database.getLayerNum());
    for (int i = 0; i < database.getLayerNum(); i++) {
        Dimension dir = database.getLayerDir(i);
        numTracks[i].resize(grid[dir].size() - 1);
        for (unsigned g = 0; g < grid[dir].size() - 1; g++) {
            utils::IntervalT<DBU> coorIntvl(grid[dir][g], grid[dir][g + 1]);
            if (grid[dir][g] >= database.getLayer(i).lastTrackLoc()) {
                numTracks[i][g] = numTracks[i][g - 1];
            } else {
                auto trackIntvl = database.rangeSearchTrack(i, coorIntvl);
                bool includeBnd = database.getLayer(i).tracks[trackIntvl.high].location == coorIntvl.high;
                numTracks[i][g] = (g == 0 ? 0 : numTracks[i][g - 1]) + trackIntvl.range() + 1 - includeBnd;
            }
        }
    }

    // print();
}


void GCellGrid::initV2() {
    db::RsynService rsynService;
    clear();
    rsynService.init();
    const Rsyn::Session session;
    Rsyn::PhysicalDesign physicalDesign =
        static_cast<Rsyn::PhysicalService *>(session.getService("rsyn.physical"))->getPhysicalDesign();
    const DBU libDBU = physicalDesign.getDatabaseUnits(Rsyn::LIBRARY_DBU);

    grid.resize(2);
    grid[0].push_back(database.dieRegion[X].low);
    grid[1].push_back(database.dieRegion[Y].low);

    // if (physicalDesign.allPhysicalGCell().empty()) {
        for (unsigned direction = 0; direction != 2; ++direction) {
             const int lo = database.dieRegion[direction].low-database.origin_offset_die;//+2120
             const int hi = database.dieRegion[direction].high+db::setting.gcellSize;
             for (int i = lo + db::setting.gcellSize; i + db::setting.gcellSize < hi; i += db::setting.gcellSize) {
                // log() << "i: "  << i << std::endl;
                // if(i < database.dieRegion[direction].low ){
                //     grid[direction].push_back(database.dieRegion[direction].low);
                //     continue;
                // }
                // if(i > database.dieRegion[direction].high ){
                //     grid[direction].push_back(database.dieRegion[direction].high);
                //     continue;
                // }
                grid[direction].push_back(i);
                 
                 
                 
             }
         }
    // } else {
    //     for (const Rsyn::PhysicalGCell &rsynGCell : physicalDesign.allPhysicalGCell()) {
    //         const int location = rsynGCell.getLocation();
    //         const int step = rsynGCell.getStep();
    //         const int numStep = rsynGCell.getNumTracks();
    //         const Dimension direction = rsynGCell.getDirection() == Rsyn::PhysicalGCellDirection::VERTICAL ? X : Y;

    //         for (int i = 1; i < numStep; ++i) {
    //             grid[direction].push_back(location + step * i);
    //         }
    //     }
    // }

    for (unsigned direction = 0; direction != 2; ++direction) {
        sort(grid[direction].begin(), grid[direction].end());
        if (grid[direction].back() < database.dieRegion[direction].high) {
            grid[direction].push_back(database.dieRegion[direction].high);
        } else if (grid[direction].back() > database.dieRegion[direction].high) {
            log() << "Warning: grid line " << grid[direction].back() << " exceeds die region " << database.dieRegion[direction].high << std::endl;
        }
    }

    numTracks.resize(database.getLayerNum());
    for (int i = 0; i < database.getLayerNum(); i++) {
        Dimension dir = database.getLayerDir(i);
        numTracks[i].resize(grid[dir].size() - 1);
        for (unsigned g = 0; g < grid[dir].size() - 1; g++) {
            utils::IntervalT<DBU> coorIntvl(grid[dir][g], grid[dir][g + 1]);
            if (grid[dir][g] >= database.getLayer(i).lastTrackLoc()) {
                numTracks[i][g] = numTracks[i][g - 1];
            } else {
                auto trackIntvl = database.rangeSearchTrack(i, coorIntvl);
                bool includeBnd = database.getLayer(i).tracks[trackIntvl.high].location == coorIntvl.high;
                numTracks[i][g] = (g == 0 ? 0 : numTracks[i][g - 1]) + trackIntvl.range() + 1 - includeBnd;
            }
        }
    }
}//initV2

GrBoxOnLayer GCellGrid::rangeSearchGCell(const db::BoxOnLayer &box) const {
    return GrBoxOnLayer(box.layerIdx, rangeSearchGCell(box[X], X), rangeSearchGCell(box[Y], Y));
}

utils::IntervalT<int> GCellGrid::rangeSearchGCell(const utils::IntervalT<DBU> &intvl, Dimension dir) const {
    int lo_idx = lower_bound(grid[dir].begin(), grid[dir].end(), intvl.low) - grid[dir].begin();
    if (grid[dir][lo_idx] > intvl.low) lo_idx--;
    int hi_idx = lower_bound(grid[dir].begin(), grid[dir].end(), intvl.high) - grid[dir].begin();
    hi_idx--;

    return utils::IntervalT<int>(lo_idx, hi_idx);
}

vector<vector<GrEdge>> GCellGrid::rangeSearchGCellEdge(const db::BoxOnLayer &box) const {
    Dimension dir = database.getLayerDir(box.layerIdx);

    auto grBox = rangeSearchGCell(box);

    vector<vector<GrEdge>> edges;
    edges.resize(grBox[dir].range() + 1);

    int jMin = max(grBox[1 - dir].low, 0);
    int jMax = min(grBox[1 - dir].high, getNumGrPoint(1 - dir) - 2);
    for (int i = grBox[dir].low; i < grBox[dir].high; i++) {
        for (int j = jMin; j <= jMax; j++) {
            int lx, ly, hx, hy;
            if (dir == X) {
                lx = hx = i;
                ly = j;
                hy = j + 1;
            } else {
                ly = hy = i;
                lx = j;
                hx = j + 1;
            }
            GrPoint p1(box.layerIdx, lx, ly), p2(box.layerIdx, hx, hy);
            edges[i].emplace_back(p1, p2);
        }
    }

    return edges;
}

utils::IntervalT<int> GCellGrid::getTrackIntvl(const GrPoint &point) const {
    return getTrackIntvl(point.layerIdx, point[database.getLayerDir(point.layerIdx)]);
}

utils::IntervalT<int> GCellGrid::getTrackIntvl(int layerIdx, int idx) const {
    return {idx == 0 ? 0 : numTracks[layerIdx][idx - 1], numTracks[layerIdx][idx] - 1};
}

void GCellGrid::print() const {
    log() << "========= gcell grid info ============" << std::endl;
    printflog("#gcell=%d*%d=%d, #edge=%d/%d\n",
              getNumGrLine(X) - 1,
              getNumGrLine(Y) - 1,
              (getNumGrLine(X) - 1) * (getNumGrLine(Y) - 1),
              (getNumGrLine(X) - 1) * (getNumGrLine(Y) - 2),
              (getNumGrLine(X) - 2) * (getNumGrLine(Y) - 1));
    for (int l = 0; l < database.getLayerNum(); l++) {
        int totTrackNum = 0;
        for (unsigned i = 0; i < numTracks[l].size(); i++) totTrackNum += getNumTracks(l, i);

        printlog(database.getLayer(l).name, ": avg #track=", totTrackNum / numTracks[l].size());
    }
}

void GCellGrid::printGrid() const {
    log() << "gridX" << std::endl;
    for(int i = 0; i < grid[X].size(); i++){
        log() << "i: " << i 
              << ", gridx: "  << grid[X][i] 
              << ", coord: " << getX(i)/database.libDBU << std::endl;
    }

    log() << "gridY" << std::endl;
    for(int i = 0; i < grid[Y].size(); i++){
        log() << "j: " << i 
              << ", gridy: "  << grid[Y][i] 
              << ", coord: " << getY(i)/database.libDBU << std::endl;
    }
}

}  // namespace gr
