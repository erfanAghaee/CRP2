#include "GrRouteGrid.h"
#include "db/Database.h"
#include "GrDatabase.h"
#include "db/Setting.h"

namespace gr {
void GrRouteGrid::init() {
    clear();
    int numLayers = database.getLayerNum();
    routedWireMap.resize(numLayers);
    fixedMetalMap.resize(numLayers);
    histWireUsageMap.resize(numLayers);
    routedViaMap.resize(numLayers);

    GCellGrid::initV2();
    // GCellGrid::init();

    for (int l = 0; l < numLayers; l++) {
        auto dir = database.getLayerDir(l);
        routedWireMap[l].resize(getNumGrPoint(dir), vector<UsageT>(getNumGrEdge(l)));
        fixedMetalMap[l].resize(getNumGrPoint(dir), vector<std::pair<int, DBU>>(getNumGrEdge(l)));
        histWireUsageMap[l].resize(getNumGrPoint(dir), vector<double>(getNumGrEdge(l), 0));
        routedViaMap[l].resize(getNumGrPoint(X), vector<UsageT>(getNumGrPoint(Y)));
    }

    markFixedMetals();
}
void GrRouteGrid::update() {
    int numLayers = database.getLayerNum();
    routedWireMap.resize(numLayers);
    fixedMetalMap.resize(numLayers);
    histWireUsageMap.resize(numLayers);
    routedViaMap.resize(numLayers);

    // // GCellGrid::init();

    for (int l = 0; l < numLayers; l++) {
        auto dir = database.getLayerDir(l);
        routedWireMap[l].resize(getNumGrPoint(dir), vector<UsageT>(getNumGrEdge(l)));
        fixedMetalMap[l].resize(getNumGrPoint(dir), vector<std::pair<int, DBU>>(getNumGrEdge(l)));
        histWireUsageMap[l].resize(getNumGrPoint(dir), vector<double>(getNumGrEdge(l), 0));
        routedViaMap[l].resize(getNumGrPoint(X), vector<UsageT>(getNumGrPoint(Y)));
    }

    markFixedMetals();
}

void GrRouteGrid::clear() {
    routedWireMap.clear();
    fixedMetalMap.clear();
    histWireUsageMap.clear();
    routedViaMap.clear();
    viaCapDiscount = 1;
    wireCapDiscount = 1;
    unitViaMultiplier = 1;
    logisticSlope = 1;
}

// void GrRouteGrid::reset() {
// }//end reset

void GrRouteGrid::markFixed(int layerIdx, int gridline, int cp, int num_track, DBU avg_length) {
    fixedMetalMap[layerIdx][gridline][cp] = std::make_pair(num_track, avg_length);
}

void GrRouteGrid::useWire(int layerIdx, int gridline, int cp, double usage) {
    // log() << "useWire Map: " 
    //       << ", layerIdx: " << layerIdx
    //       << ", gridline: " << gridline
    //       << ", cp: " << cp 
    //       << ", usage: " << usage 
    //       << ", cur_routedWireMap: " << routedWireMap[layerIdx][gridline][cp] << std::endl;
    
    routedWireMap[layerIdx][gridline][cp] += usage;
}

void GrRouteGrid::useVia(int layerIdx, int x, int y, double usage) { 
        // log() << "useVia Map: " 
        //   << ", layerIdx: " << layerIdx
        //   << ", x: " << x
        //   << ", y: " << y 
        //   << ", usage: " << usage 
        //   << ", cur_routedViaMap: " << routedViaMap[layerIdx][x][y] << std::endl;
    routedViaMap[layerIdx][x][y] += usage; }

void GrRouteGrid::useWire(const GrBoxOnLayer& box) {
    // log() << "box: " << box << std::endl;
    int layerIdx = box.layerIdx;
    auto dir = database.getLayerDir(layerIdx);
    double usage = 1.0 / (box[dir].range() + 1);
    for (int gridline = box[dir].low; gridline <= box[dir].high; gridline++)
        for (int cp = box[1 - dir].low; cp < box[1 - dir].high; cp++) useWire(layerIdx, gridline, cp, usage);
}

void GrRouteGrid::useWirePessimistic(const GrBoxOnLayer& box) {
    // log() << "box: " << box << std::endl;
    // log() << getCoor(box[X].low, X) << ", "
    //       << getCoor(box[Y].low, Y) << ", "
    //       << getCoor(box[X].high + 1, X) << ", "
    //       << getCoor(box[Y].high + 1, Y) << std::endl;
    int layerIdx = box.layerIdx;
    auto dir = database.getLayerDir(layerIdx);
    double usage = 1.0 / (box[dir].range() + 1);
    for (int gridline = box[dir].low; gridline <= box[dir].high; gridline++)
        for (int cp = box[1 - dir].low; cp < box[1 - dir].high; cp++) useWire(layerIdx, gridline, cp, usage);
}

void GrRouteGrid::useVia(const GrBoxOnLayer& box) {
    double usage = 1.0 / ((box[X].range() + 1) * (box[Y].range() + 1));
    for (int x = box[X].low; x <= box[X].high; x++)
        for (int y = box[Y].low; y <= box[Y].high; y++) useVia(box.layerIdx, x, y, usage);
}

void GrRouteGrid::useNet(const GrNet& net) {
    // log() << "useNet: " << net.getName() << std::endl;
    for (const auto& guide : net.wireRouteGuides) useWirePessimistic(guide);

    const auto& viaGuides = net.viaRouteGuides;
    for (int g1 = 0; g1 < viaGuides.size(); g1++) {
        for (int g2 = g1 + 1; g2 < viaGuides.size(); g2++) {
            if (abs(viaGuides[g1].layerIdx - viaGuides[g2].layerIdx) != 1) continue;

            auto xIntvl = viaGuides[g1][X].IntersectWith(viaGuides[g2][X]);
            auto yIntvl = viaGuides[g1][Y].IntersectWith(viaGuides[g2][Y]);

            if (xIntvl.IsValid() && yIntvl.IsValid())
                useVia({min(viaGuides[g1].layerIdx, viaGuides[g2].layerIdx), xIntvl, yIntvl});
        }
    }
}

// Note: may accumulate the size of each recorder

void GrRouteGrid::removeWire(const GrBoxOnLayer& box) {
    int layerIdx = box.layerIdx;
    auto dir = database.getLayerDir(layerIdx);
    double usage = 1.0 / (box[dir].range() + 1);
    for (int gridline = box[dir].low; gridline <= box[dir].high; gridline++)
        for (int cp = box[1 - dir].low; cp < box[1 - dir].high; cp++) {
            // log() << "layerIdx: " << layerIdx 
            //       << ", gridline: " << gridline
            //       << ", cp: " << cp << std::endl;
            routedWireMap[layerIdx][gridline][cp] -= usage;
        }
}

void GrRouteGrid::removeWireRelax(const GrBoxOnLayer& box, std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls) {
    // std::vector<std::pair<std::tuple<int,int,int>,double>> wireMap_tpls;
    int layerIdx = box.layerIdx;
    auto dir = database.getLayerDir(layerIdx);
    double usage = 1.0 / (box[dir].range() + 1);
    for (int gridline = box[dir].low; gridline <= box[dir].high; gridline++)
        for (int cp = box[1 - dir].low; cp < box[1 - dir].high; cp++) {

            std::tuple<int,int,int> tpl(layerIdx,gridline,cp);
            wireMap_tpls.push_back(std::make_pair(tpl,usage));
            // routedWireMap[layerIdx][gridline][cp] -= usage;
        }
}

void GrRouteGrid::removeViaRelax(const GrBoxOnLayer& box, std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls) {
    double usage = 1.0 / ((box[X].range() + 1) * (box[Y].range() + 1));
    for (int x = box[X].low; x <= box[X].high; x++)
        for (int y = box[Y].low; y <= box[Y].high; y++) {
            // log() << "box.layerIdx: " << box.layerIdx 
            //       << ", x: " << x
            //       << ", y: " << y << std::endl;
            std::tuple<int,int,int> tpl(box.layerIdx,x,y);
            viaMap_tpls.push_back(std::make_pair(tpl,usage));
            // routedViaMap[box.layerIdx][x][y] -= usage;
        }
    // for (auto& pair : routedViaMap[box.layerIdx][x][y])
    //     if (pair.first == netIdx) pair.second = 0;
}

void GrRouteGrid::removeVia(const GrBoxOnLayer& box) {
    double usage = 1.0 / ((box[X].range() + 1) * (box[Y].range() + 1));
    for (int x = box[X].low; x <= box[X].high; x++)
        for (int y = box[Y].low; y <= box[Y].high; y++) {
            // log() << "box.layerIdx: " << box.layerIdx 
            //       << ", x: " << x
            //       << ", y: " << y << std::endl;
            routedViaMap[box.layerIdx][x][y] -= usage;
        }
    // for (auto& pair : routedViaMap[box.layerIdx][x][y])
    //     if (pair.first == netIdx) pair.second = 0;
}

void GrRouteGrid::removeNet(GrNet& net) {
    // log() << "net: " << net.getName() << std::endl;
    // log() << "remove wire: " <<  std::endl;
    for (const auto& guide : net.wireRouteGuides) removeWire(guide);

    // log() << "remove via: " <<  std::endl;
    const auto& viaGuides = net.viaRouteGuides;
    for (int g1 = 0; g1 < viaGuides.size(); g1++) {
        for (int g2 = g1 + 1; g2 < viaGuides.size(); g2++) {
            if (abs(viaGuides[g1].layerIdx - viaGuides[g2].layerIdx) != 1) continue;

            auto xIntvl = viaGuides[g1][X].IntersectWith(viaGuides[g2][X]);
            auto yIntvl = viaGuides[g1][Y].IntersectWith(viaGuides[g2][Y]);

            if (xIntvl.IsValid() && yIntvl.IsValid())
                removeVia({min(viaGuides[g1].layerIdx, viaGuides[g2].layerIdx), xIntvl, yIntvl});
        }
    }

    net.wireRouteGuides.clear();
    net.viaRouteGuides.clear();
}

void GrRouteGrid::removeNetRelax(GrNet& net,std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
    ,std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls){

    
    for (const auto& guide : net.wireRouteGuides){
        // log() << "guide: " << guide << std::endl;
        removeWireRelax(guide,wireMap_tpls);
    } 

    // log() << "remove via: " <<  std::endl;
    const auto& viaGuides = net.viaRouteGuides;
    for (int g1 = 0; g1 < viaGuides.size(); g1++) {
        for (int g2 = g1 + 1; g2 < viaGuides.size(); g2++) {
            if (abs(viaGuides[g1].layerIdx - viaGuides[g2].layerIdx) != 1) continue;

            auto xIntvl = viaGuides[g1][X].IntersectWith(viaGuides[g2][X]);
            auto yIntvl = viaGuides[g1][Y].IntersectWith(viaGuides[g2][Y]);

            if (xIntvl.IsValid() && yIntvl.IsValid()){
                removeViaRelax({min(viaGuides[g1].layerIdx, viaGuides[g2].layerIdx), xIntvl, yIntvl},viaMap_tpls);
            }
                // removeVia({min(viaGuides[g1].layerIdx, viaGuides[g2].layerIdx), xIntvl, yIntvl});
                
                
        }
    }

    // net.wireRouteGuides.clear();
    // net.viaRouteGuides.clear();

    bool debug = false;
    if(debug){
        log() << "wire usage" << std::endl;
        for(auto wire_map_pair : wireMap_tpls){
        auto tpls = wire_map_pair.first;
        auto usage = wire_map_pair.second;
        auto layerIdx = std::get<0>(tpls);
        auto gridline = std::get<1>(tpls);
        auto cp = std::get<2>(tpls);
        log() << "layerIdx: " << layerIdx 
                << ", gridline: " << gridline
                << ", cp: " << cp 
                << ", usage: " << usage 
                << std::endl;
        }
        log() << "via usage" << std::endl;
        for(auto via_map_pair : viaMap_tpls){
            auto tpls = via_map_pair.first;
            auto usage = via_map_pair.second;
            auto layerIdx = std::get<0>(tpls);
            auto x = std::get<1>(tpls);
            auto y = std::get<2>(tpls);
            log() << "layerIdx: " << layerIdx 
                    << ", x: " << x
                    << ", y: " << y 
                    << ", usage: " << usage 
                    << std::endl;
        }
    }
    
}//end removeNetRelax

double GrRouteGrid::getWireCapacity(const GrEdge& edge) const {
    if (abs(edge.u[X] - edge.v[X]) > 1 || abs(edge.u[Y] - edge.v[Y]) > 1)
        printlog("ERROR: in GrRouteGrid::getWireCapacity, edge len > 1");
    int layerIdx = edge.getLayerIdx();
    return getNumTracks(layerIdx, edge.u[database.getLayerDir(layerIdx)]) * wireCapDiscount;
}

double GrRouteGrid::getInCellArea(const GrPoint& point) const {
    int layerIdx = point.layerIdx;
    auto dir = database.getLayerDir(layerIdx);
    return getNumTracks(layerIdx, point[dir]);
}

double GrRouteGrid::getFixedUsage(const GrEdge& edge) const {  // get num of tracks blocked by fixed metal
    if (edge.getGrLen() != 1) printlog("Error");
    auto dir = database.getLayerDir(edge.getLayerIdx());
    return getFixedUsage(edge.getLayerIdx(), edge.u[dir], edge.u[1 - dir]);
}

DBU GrRouteGrid::getFixedLength(const GrEdge& edge) const {  // get avg length of the tracks blocked by fixed metal
    if (edge.getGrLen() != 1) printlog("Error");
    auto dir = database.getLayerDir(edge.getLayerIdx());
    return fixedMetalMap[edge.getLayerIdx()][edge.u[dir]][edge.u[1 - dir]].second;
}

double GrRouteGrid::getWireUsage(const GrEdge& edge) const {
    if (edge.getGrLen() != 1) printlog("Error");
    auto dir = database.getLayerDir(edge.getLayerIdx());
    return getWireUsage(edge.getLayerIdx(), edge.u[dir], edge.u[1 - dir]);
}

double GrRouteGrid::getWireUsageRelax(const GrEdge& edge,std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls) const {
    if (edge.getGrLen() != 1) printlog("Error");
    auto dir = database.getLayerDir(edge.getLayerIdx());
    return getWireUsageRelax(edge.getLayerIdx(), edge.u[dir], edge.u[1 - dir],wireMap_tpls);
}

double GrRouteGrid::getViaUsage(const GrPoint& via) const { return getViaUsage(via.layerIdx, via[X], via[Y]); }
double GrRouteGrid::getViaUsageRelax(const GrPoint& via
                                    ,std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls) const
                                     { return getViaUsageRelax(via.layerIdx, via[X], via[Y],viaMap_tpls); }

double GrRouteGrid::getCellResource(const GrPoint& point) const { return getCellResource(point.layerIdx, point.x, point.y); }

double GrRouteGrid::getInCellUsedArea(const GrPoint& point) const {
    // note: the area defined here = # of used tracks * avg_used_length (total used length of tracks in the gcell)
    // todo: consider boundary effect
    bool debug = false;
    double used_area = 0;
    int layerIdx = point.layerIdx;
    auto dir = database.getLayerDir(layerIdx);
    int low_x = point.x - (dir == Y);
    int low_y = point.y - (dir == X);
    auto low_edge = GrEdge({layerIdx, low_x, low_y}, point);
    if (low_x >= 0 && low_y >= 0) {
        if(debug)
        log() << "getInCellUsedArea: "
              << ", low_x: " << low_x
              << ", low_y: " << low_y 
              << ", getFixedUsage(low_edge): " << getFixedUsage(low_edge) 
              << ", getWireUsage(low_edge): " << getWireUsage(low_edge) << std::endl; 
        used_area += getFixedUsage(low_edge) + getWireUsage(low_edge);
    }
    int high_x = point.x + (dir == Y);
    int high_y = point.y + (dir == X);
    auto high_edge = GrEdge(point, {layerIdx, high_x, high_y});
    if (high_x < getNumGrPoint(X) && high_y < getNumGrPoint(Y)) {
        used_area += getFixedUsage(high_edge) + getWireUsage(high_edge);
        if(debug)
        log() << "getInCellUsedArea: "
              << ", high_x: " << high_x
              << ", high_y: " << high_y 
              << ", getFixedUsage(high_edge): " << getFixedUsage(high_edge) 
              << ", getWireUsage(high_edge): " << getWireUsage(high_edge) << std::endl;
    }



    return used_area / 2;
}

double GrRouteGrid::getInCellUsedAreaRelax(const GrPoint& point
                                    , std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
                                    , std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls) const {
    // note: the area defined here = # of used tracks * avg_used_length (total used length of tracks in the gcell)
    // todo: consider boundary effect
    double used_area = 0;
    int layerIdx = point.layerIdx;
    auto dir = database.getLayerDir(layerIdx);
    int low_x = point.x - (dir == Y);
    int low_y = point.y - (dir == X);
    auto low_edge = GrEdge({layerIdx, low_x, low_y}, point);
    if (low_x >= 0 && low_y >= 0) {
        used_area += getFixedUsage(low_edge) + getWireUsageRelax(low_edge,wireMap_tpls);
    }
    int high_x = point.x + (dir == Y);
    int high_y = point.y + (dir == X);
    auto high_edge = GrEdge(point, {layerIdx, high_x, high_y});
    if (high_x < getNumGrPoint(X) && high_y < getNumGrPoint(Y)) {
        used_area += getFixedUsage(high_edge) + getWireUsageRelax(high_edge,wireMap_tpls);
    }

    return used_area / 2;
}

double GrRouteGrid::getFixedUsedArea(const GrEdge& edge) const {
    // a little different to the area of a point, this is the area of an edge, which is half of 2 gcell's area combined
    return getFixedUsage(edge) * getFixedLength(edge);
}

double GrRouteGrid::getInCellViaNum(const GrPoint& point) const {  // get the num of vias passing through the gcell
    double num = 0;
    for (int side = -1; side <= 0;
         side++) {  // a cell is used by both the via from lower layer and that from higher layer
        // side = -1, from lower layer to current layer; side = 0 from current layer to upper layer
        auto via_point = GrPoint({point.layerIdx + side, point.x, point.y});
        if (via_point.layerIdx < 0 || via_point.layerIdx >= database.getLayerNum() - 1) continue;
        // log() << "via_point: " << via_point << std::endl;
        num += getViaUsage(via_point);
    }
    return num;
}

double GrRouteGrid::getInCellViaNumRelax(const GrPoint& point, std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls) const {  // get the num of vias passing through the gcell
    double num = 0;
    for (int side = -1; side <= 0;
         side++) {  // a cell is used by both the via from lower layer and that from higher layer
        // side = -1, from lower layer to current layer; side = 0 from current layer to upper layer
        auto via_point = GrPoint({point.layerIdx + side, point.x, point.y});
        if (via_point.layerIdx < 0 || via_point.layerIdx >= database.getLayerNum() - 1) continue;
        // log() << "via_point: " << via_point << std::endl;
        num += getViaUsageRelax(via_point,viaMap_tpls);
    }
    return num;
}

// double GrRouteGrid::getUnitViaArea(const GrPoint& point, int side) const {
//     // side = -1, from lower layer to current layer; side = 0 from current layer to upper layer
//     auto& viaType = database.getCutLayer(point.layerIdx).defaultViaType();
//     auto viaBox = (side == 1 ? viaType.top : viaType.bot);
//
//     auto dir = database.getLayerDir(point.layerIdx);
//     int nBlockedTrack = viaBox[dir].range() * grDatabase.getNumTracks(point.layerIdx, point[dir]) /
//     grDatabase.getCoorIntvl(point, dir).range() + 1; return viaBox[1 - dir].range() * nBlockedTrack;
// }

double GrRouteGrid::getFixedUsage(int layerIdx,
                                  int gridline,
                                  int cp) const {  // get num of tracks blocked by fixed metal
    return fixedMetalMap[layerIdx][gridline][cp].first;
}

double GrRouteGrid::getWireUsage(int layerIdx, int gridline, int cp) const {
    return routedWireMap[layerIdx][gridline][cp];
}
double GrRouteGrid::getWireUsageRelax(int layerIdx, int gridline, int cp,std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls) const {
    double usage_tpl = 0;
    for(auto wireMap_pair : wireMap_tpls){
        auto tpls = wireMap_pair.first;
        auto layerIdx_tmp = std::get<0>(tpls);
        auto gridline_tmp = std::get<1>(tpls);
        auto cp_tmp = std::get<2>(tpls);
        if( (layerIdx == layerIdx_tmp) && 
            (gridline == gridline_tmp) && 
            (cp == cp_tmp)){
                usage_tpl = wireMap_pair.second;
        }
    }
    auto res = ((routedWireMap[layerIdx][gridline][cp]-usage_tpl) >= 0 ) ? routedWireMap[layerIdx][gridline][cp]-usage_tpl : 0;
    return res;
}

double GrRouteGrid::getViaUsageRelax(int layerIdx, int x, int y,std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls) const { 
    double usage = 0;
    for(auto via_map_pair : viaMap_tpls){
        auto tpls = via_map_pair.first;
        auto layerIdx_tmp = std::get<0>(tpls);
        auto x_tmp = std::get<1>(tpls);
        auto y_tmp = std::get<2>(tpls);
        if( (x_tmp == x) &&
            (y_tmp == y) &&
            (layerIdx_tmp == layerIdx)
        ){
            usage  = via_map_pair.second;
            break;
        }
        
    }
    auto res = ((routedViaMap[layerIdx][x][y]-usage) >= 0) ? routedViaMap[layerIdx][x][y]-usage : 0;
    return res;
    }

double GrRouteGrid::getViaUsage(int layerIdx, int x, int y) const { 
    // log() << "layerIdx: " << x 
    //       << ", x: " << x 
    //       << ", y: " << y 
    //       << ", routedViaMap: " << routedViaMap[layerIdx][x][y] << std::endl;
    return routedViaMap[layerIdx][x][y];
    }

double GrRouteGrid::getCellResource(int layerIdx, int x, int y) const {
    auto layerDir = database.getLayerDir(layerIdx);
    double totalRsrc = grDatabase.getNumTracks(layerIdx, layerDir == X ? x : y);
    double cellUsage = 0;
    int low_x = x - (layerDir == Y);
    int low_y = y - (layerDir == X);
    if (low_x >= 0 && low_y >= 0) {
        auto low_edge = gr::GrEdge({layerIdx, low_x, low_y}, {layerIdx, x, y});
        cellUsage += grDatabase.getFixedUsage(low_edge) + grDatabase.getWireUsage(low_edge);
    }
    int high_x = x + (layerDir == Y);
    int high_y = y + (layerDir == X);
    if (high_x < grDatabase.getNumGrPoint(X) && high_y < grDatabase.getNumGrPoint(Y)) {
        auto high_edge = gr::GrEdge({layerIdx, x, y}, {layerIdx, high_x, high_y});
        cellUsage += grDatabase.getFixedUsage(high_edge) + grDatabase.getWireUsage(high_edge);
    }
    cellUsage /= 2;
    cellUsage += sqrt(grDatabase.getInCellViaNum({layerIdx, x, y})) * db::setting.unitSqrtViaUsage;
    return totalRsrc - cellUsage;
}

void GrRouteGrid::print() const { GCellGrid::print(); }

void GrRouteGrid::logWireUsage(std::string log_name){
    // UsageMapT routedWireMap;  // model cross-cell routing
    bool debug = false;

    UsageMapT routedViaMap;

    int numLayers = database.getLayerNum();
    if(debug){
        log() << "num Layers: " << numLayers << std::endl;
        for (int l = 0; l < numLayers; l++) {
            auto dir = database.getLayerDir(l);
            log() << "l: " << l 
                << ", dir: " << dir
                << ", getNumGrPoint(dir): " << getNumGrPoint(dir)
                << ", getNumGrEdge(l): " << getNumGrEdge(l) << std::endl;
        }

        log() << "Writing routedViaMap to csv file..." << std::endl;
    }
    
    

    std::stringstream ss;
    ss << "layer,gridLine,cp,numWire,fixedUsage,numTracks,usage,dist,cost,\
    wire2,fix2,via2,via_u,via_v,usage2,cap2" << std::endl;

    for (int layerIdx = 0; layerIdx < database.getLayerNum(); ++layerIdx) {
        Dimension dir = database.getLayerDir(layerIdx);

        for (int gridline = 0; gridline < getNumGrPoint(dir); gridline++) {
            for (int cp = 0; cp < getNumGrEdge(layerIdx); cp++) {
                GrEdge tempEdge(layerIdx, gridline, cp);

                double numWire = getWireUsage(layerIdx, gridline, cp);
                double fixedUsage = getFixedUsage({layerIdx, gridline, cp});
                double numTracks = getNumTracks(layerIdx, gridline);
                double usage = (numWire + getFixedUsage({layerIdx, gridline, cp})) / getNumTracks(layerIdx, gridline);
                DBU dist = (getCoor(cp + 2, 1 - dir) - getCoor(cp, 1 - dir)) / 2;
                double wireCost = getWireCost(tempEdge);

                if(debug)
                    log() << "edge: " << tempEdge
                        << ", wireCost: " << wireCost << std::endl;

                // log() << "LayerIdx: " << layerIdx
                //       << ", dir: " << dir
                //       << ", getNumGrPoint(dir): "  << getNumGrPoint(dir) 
                //       << ", getNumGrEdge(layerIdx): " << getNumGrEdge(layerIdx)
                //       << std::endl;

                double numWire2 = getWireUsage(tempEdge);
                double fixedUsage2 = getFixedUsage(tempEdge);
                double viaUsage2 = sqrt((getInCellViaNum(tempEdge.u) + getInCellViaNum(tempEdge.v)) / 2) * db::setting.unitSqrtViaUsage;
                double via_u = getInCellViaNum(tempEdge.u);
                double via_v = getInCellViaNum(tempEdge.v);
                double usage2 = numWire2+fixedUsage2+viaUsage2;
                double wireCap2 = getWireCapacity(tempEdge);

                ss << std::to_string(layerIdx)
                   << "," << std::to_string(gridline)
                   << "," << std::to_string(cp)
                   << "," << std::to_string(numWire)
                   << "," << std::to_string(fixedUsage)
                   << "," << std::to_string(numTracks)
                   << "," << std::to_string(usage)
                   << "," << std::to_string(dist)
                   << "," << std::to_string(wireCost)
                //     2
                   << "," << std::to_string(numWire2)
                   << "," << std::to_string(fixedUsage2)
                   << "," << std::to_string(viaUsage2)
                   << "," << std::to_string(via_u)
                   << "," << std::to_string(via_v)
                   << "," << std::to_string(usage2)
                   << "," << std::to_string(wireCap2)

                   << std::endl;
                // int bucketIdx = buckets.size() - 1;
                // while (buckets[bucketIdx] >= usage) --bucketIdx;
                // bucketIdx = max(bucketIdx, 0);
                // wireUsageGrid[bucketIdx]++;
                // wireUsageLength[bucketIdx] += dist;
                // wirelength += dist * numWire;
            }
        }
    }

    

    // std::cout << "get RsrcUsage: " << getRsrcUsage(1,1,1) << std::endl;
    // for (int l_idx = 0; l_idx < database.getLayerNum(); l_idx++){
    //     auto dir = database.getLayerDir(l_idx);
    //     for (int g = 0; g < getNumGrPoint(dir); ++g)
    //         for (int cp = 0; cp < getNumGrEdge(l_idx); ++cp)
    //         {
    //             log() << std::to_string(l_idx) << ", " << std::to_string(g)
    //                << ", " << std::to_string(cp) << ", ";
    //             log() << routedViaMap[l_idx][g][cp]<< std::endl;

    //             // log() << std::to_string(l_idx) << ", " << std::to_string(x_idx)
    //             //       << ", " << std::to_string(y_idx) << ", "
    //             //       << std::to_string(grDatabase.getCellResource({l_idx, x_idx, y_idx})) << ", "
    //             //       << rsrcMap[l_idx][x_idx][y_idx]<< std::endl;
    //         }

    // }

    // log() <<  "routedViaMap.size():  "<< routedViaMap.size()  << std::endl;
        
    
    // std::string file_name_csv = "congestionMap.csv";
    std::ofstream fout(log_name);
    fout << ss.str();
    fout.close();
}






void GrRouteGrid::logRoutedViaMap(std::string log_name){
    // UsageMapT routedWireMap;  // model cross-cell routing
    UsageMapT routedViaMap;

    log() << "Writing routedViaMap to csv file..." << std::endl;

    std::stringstream ss;
    ss << "layer,x,y,via_usage,via_cost" << std::endl;

    for (int layerIdx = 0; (layerIdx + 1) < database.getLayerNum(); ++layerIdx) {
        for (int x = 0; x < getNumGrPoint(X); x++) {
            for (int y = 0; y < getNumGrPoint(Y); y++) {
                ss << std::to_string(layerIdx)
                   << "," <<  std::to_string(x)
                   << "," <<  std::to_string(y)
                   << "," <<  std::to_string(getViaUsage(layerIdx, x, y))
                   << "," <<  std::to_string(getViaCost(gr::GrPoint(layerIdx,  x, y)))
                   << std::endl;
            }
        }
    }

    // std::cout << "get RsrcUsage: " << getRsrcUsage(1,1,1) << std::endl;
    // for (int l_idx = 0; l_idx < database.getLayerNum(); l_idx++){
    //     auto dir = database.getLayerDir(l_idx);
    //     for (int g = 0; g < getNumGrPoint(dir); ++g)
    //         for (int cp = 0; cp < getNumGrEdge(l_idx); ++cp)
    //         {
    //             log() << std::to_string(l_idx) << ", " << std::to_string(g)
    //                << ", " << std::to_string(cp) << ", ";
    //             log() << routedViaMap[l_idx][g][cp]<< std::endl;

    //             // log() << std::to_string(l_idx) << ", " << std::to_string(x_idx)
    //             //       << ", " << std::to_string(y_idx) << ", "
    //             //       << std::to_string(grDatabase.getCellResource({l_idx, x_idx, y_idx})) << ", "
    //             //       << rsrcMap[l_idx][x_idx][y_idx]<< std::endl;
    //         }

    // }

    // log() <<  "routedViaMap.size():  "<< routedViaMap.size()  << std::endl;
        
    
    // std::string file_name_csv = "congestionMap.csv";
    std::ofstream fout(log_name);
    fout << ss.str();
    fout.close();
}

void GrRouteGrid::printAllUsageAndVio() const {
    const int width = 10;
    auto wlVia = printAllUsage();
    double numShort = printAllVio();
    log() << "--- Estimated Scores ---" << std::endl;
    vector<std::string> items = {"wirelength", "# vias", "short"};
    vector<double> metrics = {wlVia.first, wlVia.second, numShort};
    vector<double> weights = {db::setting.weightWirelength, db::setting.weightViaNum, db::setting.weightShortArea};
    double totalScore = 0;
    for (int i = 0; i < items.size(); ++i) {
        totalScore += metrics[i] * weights[i];
    }
    log() << std::setw(width) << "item"
          << " | " << std::setw(width + 2) << "metric"
          << " | " << std::setw(width) << "weight"
          << " | " << std::setw(width + 2) << "score"
          << " | " << std::setw(width) << "\%" << std::endl;
    for (int i = 0; i < items.size(); ++i) {
        double score = metrics[i] * weights[i];
        log() << std::setw(width) << items[i] << " | " << std::setw(width + 2) << metrics[i] << " | "
              << std::setw(width) << weights[i] << " | " << std::setw(width + 2) << score << " | " << std::setw(width)
              << score / totalScore << std::endl;
    }
    log() << "total score = " << totalScore << std::endl;
}

double GrRouteGrid::getAllWireUsage(const vector<double>& buckets,
                                    vector<int>& wireUsageGrid,
                                    vector<DBU>& wireUsageLength) const {
    double wirelength = 0;
    wireUsageGrid.assign(buckets.size(), 0);
    wireUsageLength.assign(buckets.size(), 0);
    for (int layerIdx = 0; layerIdx < database.getLayerNum(); ++layerIdx) {
        Dimension dir = database.getLayerDir(layerIdx);

        for (int gridline = 0; gridline < getNumGrPoint(dir); gridline++) {
            for (int cp = 0; cp < getNumGrEdge(layerIdx); cp++) {
                double numWire = getWireUsage(layerIdx, gridline, cp);
                double usage = (numWire + getFixedUsage({layerIdx, gridline, cp})) / getNumTracks(layerIdx, gridline);
                DBU dist = (getCoor(cp + 2, 1 - dir) - getCoor(cp, 1 - dir)) / 2;
                int bucketIdx = buckets.size() - 1;
                while (buckets[bucketIdx] >= usage) --bucketIdx;
                bucketIdx = max(bucketIdx, 0);
                wireUsageGrid[bucketIdx]++;
                wireUsageLength[bucketIdx] += dist;
                wirelength += dist * numWire;
            }
        }
    }

    return wirelength;
}

double GrRouteGrid::getWirelength() const {
    double wirelength = 0;
    for (auto& net : grDatabase.nets) wirelength += net.getWirelength();
    return wirelength;
}

void GrRouteGrid::getAllInCellUsage(const vector<double>& buckets, vector<int>& inCellUsage) const {
    inCellUsage.assign(buckets.size(), 0);
    for (int layerIdx = 0; layerIdx < database.getLayerNum(); ++layerIdx) {
        for (int x = 0; x < getNumGrPoint(X); x++) {
            for (int y = 0; y < getNumGrPoint(Y); y++) {
                double usage = getInCellUsedArea({layerIdx, x, y}) / getInCellArea({layerIdx, x, y});
                int bucketIdx = buckets.size() - 1;
                while (buckets[bucketIdx] >= usage) --bucketIdx;
                bucketIdx = max(bucketIdx, 0);
                inCellUsage[bucketIdx]++;
            }
        }
    }
}

double GrRouteGrid::getTotViaNum() const {
    double viaNum = 0;
    for (int layerIdx = 0; (layerIdx + 1) < database.getLayerNum(); ++layerIdx) {
        for (int x = 0; x < getNumGrPoint(X); x++) {
            for (int y = 0; y < getNumGrPoint(Y); y++) {
                viaNum += getViaUsage(layerIdx, x, y);
            }
        }
    }
    return viaNum;
}

std::pair<double, double> GrRouteGrid::printAllUsage() const {
    const int width = 10;
    vector<double> buckets = {
        -1, 0, 0.3, 0.6, 0.8, 0.9, 1, 1.1, 1.3, 1.5, 2, 3};  // the i-th bucket: buckets[i] <= x < buckets[i+1]

    auto getRangeStr = [](const vector<double>& buckets, int i) {
        std::string range;
        if (i == 0) {
            range = "      " + std::to_string_with_precision(0.0, 2) + " ";
        } else if ((i + 1) < buckets.size()) {
            range = "(" + std::to_string_with_precision(buckets[i], 2) + "~" +
                    std::to_string_with_precision(buckets[i + 1], 2) + "]";
        } else {
            range = "(" + std::to_string_with_precision(buckets[i], 2) + "~inf" + ")";
        }
        return range;
    };

    // Wire
    vector<int> routedWireUsageGrid;
    vector<DBU> routedWireUsageLength;
    double wireLength = getAllWireUsage(buckets, routedWireUsageGrid, routedWireUsageLength);
    log() << "total wireLength: " << wireLength << std::endl;
    log() << "--- Wire Usage ---" << std::endl;
    log() << std::setw(width) << "usage"
          << " | " << std::setw(width) << "   grid   "
          << " | " << std::setw(width) << "  length  " << std::endl;
    for (int i = 0; i < buckets.size(); ++i) {
        if (routedWireUsageGrid[i] == 0 && routedWireUsageLength[i] == 0) continue;

        log() << std::setw(width) << getRangeStr(buckets, i) << " | " << std::setw(width) << routedWireUsageGrid[i]
              << " | " << std::setw(width) << routedWireUsageLength[i] / double(database.getLayer(1).pitch)
              << std::endl;
    }
    wireLength /= double(database.getLayer(1).pitch);

    log() << "database.getLayer(1).pitch: " << database.getLayer(1).pitch << std::endl;
    log() << "wireLength after pitch div: " << wireLength << std::endl;
    // in-Cell
    vector<int> routedViaUsage;
    getAllInCellUsage(buckets, routedViaUsage);
    log() << "--- in-Cell Usage ---" << std::endl;
    log() << std::setw(width) << "usage"
          << " | " << std::setw(width) << "routed" << std::endl;
    for (int i = 0; i < buckets.size(); ++i) {
        if (routedViaUsage[i] == 0) continue;

        log() << std::setw(width) << getRangeStr(buckets, i) << " | " << std::setw(width) << routedViaUsage[i]
              << std::endl;
    }

    double viaNum = getTotViaNum();
    return {wireLength, viaNum};
}

double GrRouteGrid::printAllVio() const {
    // log() << "printAllVio" << std::endl;
    const int width = 10;
    auto sumVec = [](const vector<int>& vec) {
        int sum = 0;
        for (int val : vec) {
            sum += val;
        }
        return sum;
    };

    // Wire violations
    vector<double> shortLen(database.getLayerNum(), 0.0);
    for (int layerIdx = 0; layerIdx < database.getLayerNum(); ++layerIdx) {
        Dimension dir = database.getLayerDir(layerIdx);

        for (int gridline = 0; gridline < getNumGrPoint(dir); gridline++) {
            for (int cp = 0; cp < getNumGrEdge(layerIdx); cp++) {
                double numWire = getWireUsage(layerIdx, gridline, cp) + getFixedUsage(layerIdx, gridline, cp);
                DBU dist = (getCoor(cp + 2, 1 - dir) - getCoor(cp, 1 - dir)) / 2;

                double overflow = max(0.0, numWire - getNumTracks(layerIdx, gridline));

                if(overflow>0){
                    GrEdge tempEdge(layerIdx,gridline,cp);         
                    auto lx = getCoor(tempEdge.lowerGrPoint().x, X)/2000.0;
                    auto ly = getCoor(tempEdge.lowerGrPoint().y, Y)/2000.0;
                    auto hx = getCoor(tempEdge.upperGrPoint().x+1, X)/2000.0;
                    auto hy = getCoor(tempEdge.upperGrPoint().y+1, Y)/2000.0;

                    log() << "vio layer: " << layerIdx
                          << ", dir: " << dir
                          << ", gridline: " << gridline
                          << ", cp: " << cp
                          << ", wireusage: " << getWireUsage(layerIdx, gridline, cp)
                          << ", fixedusage: " << getFixedUsage(layerIdx, gridline, cp)
                          << ", tracks: " << getNumTracks(layerIdx, gridline) 
                          << ", lx: " << lx 
                          << ", ly: " << ly 
                          << ", hx: " << hx 
                          << ", hy: " << hy 
                          << std::endl;


                }

                shortLen[layerIdx] += overflow * dist;
            }
        }
    }

    log() << "--- Wire-Wire Short Vios ---" << std::endl;
    log() << std::setw(width) << "usage"
          << " | " << std::setw(width * 2 + 3) << "      short area     " << std::endl;
    log() << std::setw(width) << "layer"
          << " | " << std::setw(width) << "wire-wire" << std::endl;
    double routedShortArea = 0;
    for (int i = 0; i < database.getLayerNum(); ++i) {
        if (shortLen[i] == 0) continue;
        const auto& layer = database.getLayer(i);
        double routedArea = double(shortLen[i]) * layer.width / database.getLayer(1).pitch / database.getLayer(1).pitch;
        log() << std::setw(width) << database.getLayer(i).name << " | " << std::setw(width) << routedArea << std::endl;
        routedShortArea += routedArea;
    }
    log() << std::setw(width) << "SumW"
          << " | " << std::setw(width) << routedShortArea << std::endl;

    // Via violations
    vector<double> shortNum(database.getLayerNum() - 1, 0.0);
    for (int layerIdx = 0; layerIdx < database.getLayerNum() - 1; ++layerIdx) {
        for (int x = 0; x < getNumGrPoint(X); x++) {
            for (int y = 0; y < getNumGrPoint(Y); y++) {
                double overflow = getNumVio(GrPoint({layerIdx, x, y}), 0);
                shortNum[layerIdx] += overflow;
            }
        }
    }

    log() << "--- Via-Via Short Vios ---" << std::endl;
    log() << std::setw(width) << "usage"
          << " | " << std::setw(width * 2 + 3) << "      #short     " << std::endl;
    log() << std::setw(width) << "layer"
          << " | " << std::setw(width) << "via-via" << std::endl;
    double routedShortViaNum = 0;
    for (int i = 0; i < database.getLayerNum() - 1; ++i) {
        if (shortNum[i] == 0) continue;
        log() << std::setw(width) << database.getCutLayer(i).name << " | " << std::setw(width) << shortNum[i]
              << std::endl;
        routedShortViaNum += shortNum[i];
    }
    log() << std::setw(width) << "SumW"
          << " | " << std::setw(width) << routedShortViaNum << std::endl;
    routedShortArea += routedShortViaNum;
    return routedShortArea;
}//end printAllVio

double GrRouteGrid::getAllVio() const {
    const int width = 10;
    auto sumVec = [](const vector<int>& vec) {
        int sum = 0;
        for (int val : vec) {
            sum += val;
        }
        return sum;
    };

    // Wire violations
    vector<double> shortLen(database.getLayerNum(), 0.0);
    for (int layerIdx = 0; layerIdx < database.getLayerNum(); ++layerIdx) {
        Dimension dir = database.getLayerDir(layerIdx);

        for (int gridline = 0; gridline < getNumGrPoint(dir); gridline++) {
            for (int cp = 0; cp < getNumGrEdge(layerIdx); cp++) {
                double numWire = getWireUsage(layerIdx, gridline, cp) + getFixedUsage(layerIdx, gridline, cp);
                DBU dist = (getCoor(cp + 2, 1 - dir) - getCoor(cp, 1 - dir)) / 2;

                double overflow = max(0.0, numWire - getNumTracks(layerIdx, gridline));
                shortLen[layerIdx] += overflow * dist;
            }
        }
    }

  
    double routedShortArea = 0;
    for (int i = 0; i < database.getLayerNum(); ++i) {
        if (shortLen[i] == 0) continue;
        const auto& layer = database.getLayer(i);
        double routedArea = double(shortLen[i]) * layer.width / database.getLayer(1).pitch / database.getLayer(1).pitch;
        routedShortArea += routedArea;
    }
  
    // Via violations
    vector<double> shortNum(database.getLayerNum() - 1, 0.0);
    for (int layerIdx = 0; layerIdx < database.getLayerNum() - 1; ++layerIdx) {
        for (int x = 0; x < getNumGrPoint(X); x++) {
            for (int y = 0; y < getNumGrPoint(Y); y++) {
                double overflow = getNumVio(GrPoint({layerIdx, x, y}), 0);
                shortNum[layerIdx] += overflow;
            }
        }
    }

   
    double routedShortViaNum = 0;
    for (int i = 0; i < database.getLayerNum() - 1; ++i) {
        if (shortNum[i] == 0) continue;
        routedShortViaNum += shortNum[i];
    }
    routedShortArea += routedShortViaNum;
    return routedShortArea;
}

void GrRouteGrid::markFixedMetals() {
    log() << "markFixedMetals... " << std::endl;
    bool debug = true;

    if(debug){
        log() << "markFixedMetals" << std::endl;
        db::BoxOnLayer box(1,
                        189.7,//23.35,
                        289.7,//4.3,
                        190.4,//23.36,
                        291.0);//4.6);
        auto tmps = database.getOvlpFixedMetals(box, -2);

        log() << "overlap size: " << tmps.size() <<  std::endl;
        for(auto tmp : tmps){
            log() << "overlaps: " << tmp.first << std::endl;
        }
    }


    for (int l = 0; l < database.getLayerNum(); l++) {
        std::unordered_map<std::pair<int, int>,
                           vector<std::pair<utils::IntervalT<int>, DBU>>,
                           boost::hash<std::pair<int, int>>>
            markingBuffer;  // (gridline, cp) -> (interval,netIdx)

        Dimension dir = database.getLayerDir(l);
        const RTree& tree = database.getFixedMetals(l);
        auto bit = tree.qbegin(bgi::satisfies([](auto const&) { return true; })), eit = tree.qend();
        for (auto iter = bit; iter != eit; iter++) {
            const auto& pair = *iter;
            int netIdx = pair.second;

            db::BoxOnLayer box(l,
                               bg::get<bg::min_corner, 0>(pair.first),
                               bg::get<bg::min_corner, 1>(pair.first),
                               bg::get<bg::max_corner, 0>(pair.first),
                               bg::get<bg::max_corner, 1>(pair.first));

            debug = false;
            if(debug){
                log() << "net: " << database.nets[netIdx].getName() 
                     << ", iter box: " << box << std::endl;
            }
            debug = false;
            // compute the forbid region of a fixed metal
            db::AggrParaRunSpace aggr = db::AggrParaRunSpace::DEFAULT;
            if (database.getLayer(0).parallelLength.size() <= 1) {
                // hack for ISPD'18 test cases
                aggr = db::AggrParaRunSpace::LARGER_WIDTH;
                if (min(box.width(), box.height()) == database.getLayer(box.layerIdx).width &&
                    database.getOvlpFixedMetals(box, -2).size() == 1) {
                    aggr = db::AggrParaRunSpace::DEFAULT;
                }
            } else {
                // hack for ISPD'19 test cases
                aggr = db::AggrParaRunSpace::LARGER_LENGTH;
            }
            auto forbidRegion = database.getMetalRectForbidRegion(box, aggr);
            auto gridBox = database.rangeSearch(forbidRegion, aggr == db::AggrParaRunSpace::LARGER_WIDTH);  // TODO: change to false
            if (!database.isValid(gridBox)) continue;
            box = database.getLoc(gridBox);

            // log() << "grroutegrid markFixedMetals rangeSearchGCell before (box): " << box << std::endl;

            auto grBox = rangeSearchGCell(box);

            // log() << "grroutegrid markFixedMetals rangeSearchGCell after (grBox): " << grBox << std::endl;
            auto trackIntvl = database.rangeSearchTrack(l, box[dir]);
            if (!trackIntvl.IsValid()) {
                continue;
            }

            // mark wire usage
            const int iMin = max(grBox[dir].low, 0);
            const int iMax = min(grBox[dir].high, getNumGrPoint(dir) - 1);
            const int jMin = max(grBox[1 - dir].low - 1, 0);
            const int jMax = min(grBox[1 - dir].high, getNumGrPoint(1 - dir) - 2);

            

            for (int i = iMin; i <= iMax; ++i) {
                for (int j = jMin; j <= jMax; ++j) {
                    utils::IntervalT<DBU> gcellIntvl1 = {getCoor(j, 1 - dir), getCoor(j + 1, 1 - dir)};
                    utils::IntervalT<DBU> gcellIntvl2 = {getCoor(j + 1, 1 - dir), getCoor(j + 2, 1 - dir)};
                    utils::IntervalT<DBU> edgeIntvl = {gcellIntvl1.center(), gcellIntvl2.center()};

                    auto blocked_length = box[1 - dir].IntersectWith(edgeIntvl).range();
                    if (blocked_length > 0) {
                        auto gcellTrackIntvl = getTrackIntvl(l, i);
                        auto blockedIntvl = gcellTrackIntvl.IntersectWith(trackIntvl);
                        if (blockedIntvl.IsValid()){
                            markingBuffer[std::make_pair(i, j)].emplace_back(
                                std::make_pair(blockedIntvl, blocked_length));

                            bool debug=false;
                            if(debug){
                                if(l==3 && 
                                    i == 366
                                    & j == 277){
                                    log() << "l" << l << ", box: " << box << std::endl;
                                    log() << "iMin: " << iMin 
                                        << ", iMax: " << iMax 
                                        << ", jMin: " << jMin 
                                        << ", jMax: " << jMax << std::endl;
                                }
                            }
                        }
                            
                    }
                }
            }
        }

        for (auto& buf : markingBuffer) {
            bool debug = false;
            // if(l==3 && buf.first.first == 366 & buf.first.second == 277){
            //     debug=true;
            // }
            auto gcellTrackIntvl = getTrackIntvl(l, buf.first.first);
            vector<DBU> trackBlocked(gcellTrackIntvl.range() + 1, 0);  // blocked track length
            for (auto& pair : buf.second) {
                for (int t = pair.first.low; t <= pair.first.high; t++){
                    if(debug){
                        log() << "t: " << t 
                              << ", pair.second: " << pair.second << std::endl;
                    }
                    trackBlocked[t - gcellTrackIntvl.low] += pair.second;
                }
                    
            }
            int num_blocked = 0;
            DBU avg_blocked_len = 0;

            for (auto& len : trackBlocked) {
                if (len > 0) {
                    num_blocked++;
                    avg_blocked_len += len;
                }
            }
            avg_blocked_len /= num_blocked;


            debug = false;
            if(debug){
                if(l==3 && buf.first.first == 366 & buf.first.second == 277){
                    log() << "num_blocked: " << num_blocked
                          << ", trackInterval: " << gcellTrackIntvl << std::endl;
                }
            }
            

            markFixed(l, buf.first.first, buf.first.second, num_blocked, avg_blocked_len);
        }
    }
}

db::CostT GrRouteGrid::getViaCost(const GrPoint& via,bool debug) const {
    return database.getUnitViaCost() * (unitViaMultiplier + getViaShortCost(via,debug));
}

db::CostT GrRouteGrid::getViaCostRelax(const GrPoint& via
                                    ,std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
                                    ,std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls
                                    ,bool debug) const {
    return database.getUnitViaCost() * (unitViaMultiplier + getViaShortCostRelax(via,wireMap_tpls,viaMap_tpls,debug));
}

db::CostT GrRouteGrid::getStackViaCost(const GrPoint& via, int height,bool debug) const {
    db::CostT cost = 0;
    for (int i = 0; i < height; i++) {
        cost += getViaCost(gr::GrPoint(via.layerIdx + i, via.x, via.y),debug);
    }
    return cost;
};

db::CostT GrRouteGrid::getStackViaCostRelax(const GrPoint& via, int height
    ,std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
    ,std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls
    ,bool debug) const {
    db::CostT cost = 0;
    for (int i = 0; i < height; i++) {
        cost += getViaCostRelax(gr::GrPoint(via.layerIdx + i, via.x, via.y),wireMap_tpls,viaMap_tpls,debug);
    }
    return cost;
};

db::CostT GrRouteGrid::getWireDistCost(const GrEdge& edge) const {
    Dimension dir = database.getLayerDir(edge.getLayerIdx());
    return getDist(edge.lowerGrPoint(), edge.upperGrPoint(), 1 - dir);
}

// logistic cost, ref: nctugr (8)
// " With a growing number of rip-up and rerouting iterations,
// the value of slope also grows to free the edge usage constraint. "

db::CostT GrRouteGrid::getWireShortCost(const GrEdge& edge, bool debug) const {
    debug = false;
    // Note: didn't consider occurrence penalty
    int layerIdx = edge.getLayerIdx();
    auto dir = database.getLayerDir(layerIdx);
    db::CostT cost = 0;


    int gridline = edge.u[dir];
    for (int i = edge.u[1 - dir]; i < edge.v[1 - dir]; i++) {
        GrEdge tempEdge(layerIdx, gridline, i);

        auto fixed_used = getFixedUsage(tempEdge);
        auto wire_used = getWireUsage(tempEdge);

        if(debug){
            log() << "edge: " << tempEdge << std::endl;
            log() << "tmp_edge: [l: " << layerIdx << ", gridline: " << gridline 
                  << "i: " << i
                  << "], fixed_used: " << fixed_used 
                  << ", wire_used: " << wire_used << std::endl;
        }

        DBU expected_of_len = fixed_used * getFixedLength(tempEdge) + wire_used * getWireDistCost(tempEdge);

        if(debug){
            log() << "expected_of_len before division: " << expected_of_len << std::endl;
        }

        expected_of_len /= max(grDatabase.getNumTracks(layerIdx, edge.u[dir]), 1);

        if(debug){
            log() << "expected_of_len after division: " << expected_of_len << std::endl;
        }



        auto demand = fixed_used + wire_used + 1;

        if(debug){
            log() << "demand before sqrt: " << demand << std::endl;
        }

        demand += sqrt((getInCellViaNum(tempEdge.u) + getInCellViaNum(tempEdge.v)) / 2) * db::setting.unitSqrtViaUsage;

        if(debug){
            log() << "demand after sqrt: " << demand << std::endl;
        }

        auto capacity = getWireCapacity(tempEdge);

        if(debug){
            log() << "capacity: " << capacity << std::endl;
        }
        

        cost += expected_of_len / (1.0 + exp(-logisticSlope * (demand - capacity)));

        if(debug){
            log() << "expected_of_len / (1.0 + exp(-logisticSlope * (demand - capacity))): "
                  << expected_of_len / (1.0 + exp(-logisticSlope * (demand - capacity))) << std::endl;
        }
    }
    if(debug)
        log() << "total_cost: " << cost 
            << ", cost * database.getUnitShortCost(layerIdx): " << std::endl;

    return cost * database.getUnitShortCost(layerIdx);
}

db::CostT GrRouteGrid::getWireShortCostRelax(const GrEdge& edge
                                            , std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
                                            , std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls
                                            , bool debug) const {
    debug = false;
    // Note: didn't consider occurrence penalty
    int layerIdx = edge.getLayerIdx();
    auto dir = database.getLayerDir(layerIdx);

    db::CostT cost = 0;

    int gridline = edge.u[dir];
    for (int i = edge.u[1 - dir]; i < edge.v[1 - dir]; i++) {
        GrEdge tempEdge(layerIdx, gridline, i);

        // double usage_tpl = 0;
        // for(auto wireMap_pair : wireMap_tpls){
        //     auto tpls = wireMap_pair.first;
        //     auto layerIdx_tmp = std::get<0>(tpls);
        //     auto gridline_tmp = std::get<1>(tpls);
        //     auto cp_tmp = std::get<2>(tpls);
        //     if( (layerIdx == layerIdx_tmp) && 
        //         (gridline == gridline_tmp) && 
        //         (i == cp_tmp)){
        //             usage_tpl = wireMap_pair.second;
        //     }
        // }

        auto fixed_used = getFixedUsage(tempEdge);
        auto wire_used = getWireUsageRelax(tempEdge,wireMap_tpls);
        // wire_used = wire_used - usage_tpl;

        if(debug){
            log() << "edge: " << tempEdge << std::endl;
            log() << "tmp_edge: [l: " << layerIdx << ", gridline: " << gridline 
                  << "i: " << i
                  << "], fixed_used: " << fixed_used 
                  << ", wire_used: " << wire_used << std::endl;
                  
        }

        DBU expected_of_len = fixed_used * getFixedLength(tempEdge) + wire_used * getWireDistCost(tempEdge);

        if(debug){
            log() << "expected_of_len before division: " << expected_of_len << std::endl;
        }

        expected_of_len /= max(grDatabase.getNumTracks(layerIdx, edge.u[dir]), 1);

        if(debug){
            log() << "expected_of_len after division: " << expected_of_len << std::endl;
        }



        auto demand = fixed_used + wire_used + 1;

        if(debug){
            log() << "demand before sqrt: " << demand << std::endl;
            log() << "tempEdge.u: " << tempEdge.u << std::endl;
            log() << "tempEdge.v: " << tempEdge.v << std::endl;
            log() << "getInCellViaNumRelax(tempEdge.u,viaMap_tpls): " << getInCellViaNumRelax(tempEdge.u,viaMap_tpls) << std::endl;
            log() << "getInCellViaNumRelax(tempEdge.v,viaMap_tpls): " << getInCellViaNumRelax(tempEdge.v,viaMap_tpls) << std::endl;
        }

        demand += sqrt((getInCellViaNumRelax(tempEdge.u,viaMap_tpls) + getInCellViaNumRelax(tempEdge.v,viaMap_tpls)) / 2) * db::setting.unitSqrtViaUsage;

        if(debug){
            log() << "demand after sqrt: " << demand << std::endl;
        }

        auto capacity = getWireCapacity(tempEdge);

        if(debug){
            log() << "capacity: " << capacity << std::endl;
        }
        

        cost += expected_of_len / (1.0 + exp(-logisticSlope * (demand - capacity)));

        if(debug){
            log() << "expected_of_len / (1.0 + exp(-logisticSlope * (demand - capacity))): "
                  << expected_of_len / (1.0 + exp(-logisticSlope * (demand - capacity))) << std::endl;
        }
    }
    if(debug)
        log() << "total_cost: " << cost 
            << ", cost * database.getUnitShortCost(layerIdx): " << std::endl;

    return cost * database.getUnitShortCost(layerIdx);
}

db::CostT GrRouteGrid::getViaShortCost(const GrPoint& via, bool debug) const {
    double cost = 0;
    for (int side = 0; side <= 1; side++) {
        auto side_point = GrPoint({via.layerIdx + side, via.x, via.y});
        double incell_used = getInCellUsedArea(side_point);
        double incell_area = getInCellArea(side_point);
        if(debug){
            log() << "viaCost: [l: " << via.layerIdx + side 
                  << ", x: " << via.x 
                  << ", y: " << via.y  
                  << "], incellUsedArea: " << incell_used
                  << ", incell_area: " << incell_area
                  << ", logisticSlope: " << logisticSlope << std::endl;
        }

        cost += 1.0 / (1.0 + exp(-logisticSlope * (incell_used - incell_area)));

        if(debug){
            log() << "1.0 / (1.0 + exp(-logisticSlope * (incell_used - incell_area))): "
                  << 1.0 / (1.0 + exp(-logisticSlope * (incell_used - incell_area)))
                  << ", cost: " << cost << std::endl;
        }
    }
    if(debug){
            log() << "total cost: " << cost << std::endl;
        }
    return cost;
    // return cost;
}

db::CostT GrRouteGrid::getViaShortCostRelax(const GrPoint& via
                                    , std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
                                    , std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls
                                    , bool debug) const {
    debug = false;
    double cost = 0;
    for (int side = 0; side <= 1; side++) {
        auto side_point = GrPoint({via.layerIdx + side, via.x, via.y});
        double incell_used = getInCellUsedAreaRelax(side_point,wireMap_tpls,viaMap_tpls);
        double incell_area = getInCellArea(side_point);
        if(debug){
            log() << "viaCost: [l: " << via.layerIdx + side 
                  << ", x: " << via.x 
                  << ", y: " << via.y  
                  << "], incellUsedArea: " << incell_used
                  << ", incell_area: " << incell_area
                  << ", logisticSlope: " << logisticSlope << std::endl;
        }

        cost += 1.0 / (1.0 + exp(-logisticSlope * (incell_used - incell_area)));

        if(debug){
            log() << "1.0 / (1.0 + exp(-logisticSlope * (incell_used - incell_area))): "
                  << 1.0 / (1.0 + exp(-logisticSlope * (incell_used - incell_area)))
                  << ", cost: " << cost << std::endl;
        }
    }
    if(debug){
            log() << "total cost: " << cost << std::endl;
        }
    return cost;
    // return cost;
}

db::CostT GrRouteGrid::getWireCost(const GrEdge& edge,bool debug) const {
    // debug = false;  
    if (edge.getLayerIdx() == 0) return LARGE_NUM;
    if(debug){
        log() << "getWireDistCost(edge): " << getWireDistCost(edge)
              << "getWireShortCost(edge): " << getWireShortCost(edge) << std::endl;
    }
    return getWireDistCost(edge) + getWireShortCost(edge,debug);
    // return db::setting.wirelenCostWeight*getWireDistCost(edge) + getWireShortCost(edge);
}

db::CostT GrRouteGrid::getWireCostRelax(const GrEdge& edge
                                        ,std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
                                        ,std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls
                                        ,bool debug) const {
                                                      
    if (edge.getLayerIdx() == 0) return LARGE_NUM;
    if(debug){
        log() << "getWireDistCost(edge): " << getWireDistCost(edge)
              << "getWireShortCost(edge): " << getWireShortCostRelax(edge,wireMap_tpls,viaMap_tpls) << std::endl;
    }
    return getWireDistCost(edge) + getWireShortCostRelax(edge,wireMap_tpls,viaMap_tpls,debug);
    // return db::setting.wirelenCostWeight*getWireDistCost(edge) + getWireShortCost(edge);
}

bool GrRouteGrid::hasVio(const GrNet& net, bool hasCommit) const { return getNumVio(net, hasCommit) > 0; }

bool GrRouteGrid::hasVio(const GrEdge& edge, bool hasCommit) const { return getNumVio(edge, !hasCommit) > 0; }
bool GrRouteGrid::hasVioRelax(const GrEdge& edge,std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
    ,std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls, bool hasCommit) const { 
        return getNumVioRelax(edge,wireMap_tpls,viaMap_tpls, !hasCommit) > 0; 
        }

bool GrRouteGrid::hasVio(const GrPoint& via, bool hasCommit) const { return getNumVio(via, !hasCommit) > 0; }

double GrRouteGrid::getNumVio(const GrNet& net, bool hasCommit) const {
    double numVio = 0;

    // log() << "getNumVio net: " << net.getName()
    //       << ", wl: " << net.getWirelength() <<  std::endl;

    const auto& guides = net.wireRouteGuides;
    for (const auto& guide : guides) {
        // log() << "getNumVio guide: " << guide << std::endl;
        auto dir = database.getLayerDir(guide.layerIdx);
        double usage = 1.0 / (guide[dir].range() + 1);
        for (int gridline = guide[dir].low; gridline <= guide[dir].high; gridline++) {
            GrEdge tempEdge = (dir == X) ? GrEdge(guide.layerIdx, gridline, guide[Y].low, guide[Y].high)
                                         : GrEdge(guide.layerIdx, gridline, guide[X].low, guide[X].high);
            auto numVioTmp = getNumVio(tempEdge, hasCommit ? 0 : usage);
            numVio = numVio + numVioTmp;

            // log() << "numVioTmpWireGuide: "  << numVioTmp << std::endl;


        }
    }

    const auto& viaGuides = net.viaRouteGuides;
    for (int g1 = 0; g1 < viaGuides.size(); g1++) {
        for (int g2 = g1 + 1; g2 < viaGuides.size(); g2++) {
            if (abs(viaGuides[g1].layerIdx - viaGuides[g2].layerIdx) != 1) continue;

            auto xIntvl = viaGuides[g1][X].IntersectWith(viaGuides[g2][X]);
            auto yIntvl = viaGuides[g1][Y].IntersectWith(viaGuides[g2][Y]);

            if (xIntvl.IsValid() && yIntvl.IsValid()) {
                double usage = 1.0 / ((xIntvl.range() + 1) * (yIntvl.range() + 1));

                for (int x = xIntvl.low; x <= xIntvl.high; x++)
                    for (int y = yIntvl.low; y <= yIntvl.high; y++){
                        int numVioTmp = getNumVio(GrPoint(min(viaGuides[g1].layerIdx, viaGuides[g2].layerIdx), x, y),
                                            hasCommit ? 0 : usage);
                        
                        numVio = numVio + numVioTmp;
                        // log() << "numVioTmpViaGuide: "  << numVioTmp << std::endl;
                    }
                        
            }
        }
    }

    // log() << "total_vio: " << numVio << std::endl;

    return numVio;
}

bool GrRouteGrid::getVioReport(const GrNet& net
                              , std::vector<GrEdge>& edges_viol
                              , std::vector<GrPoint>& vias_viol
                              , std::unordered_map<std::tuple<int,int,int>,std::set<int>>& viol_dict
                              , bool hasCommit) const
{

    

    double numVio = 0;
    // log() << "getVioReport" << std::endl;


    // log() << "getNumVio net: " << net.getName()
    //       << ", wl: " << net.getWirelength() <<  std::endl;

    const auto& guides = net.wireRouteGuides;

    for(auto box : guides){
        int layerIdx = box.layerIdx;
        auto dir = database.getLayerDir(layerIdx);
        double usage = 1.0 / (box[dir].range() + 1);
        for (int gridline = box[dir].low; gridline <= box[dir].high; gridline++)
            for (int cp = box[1 - dir].low; cp < box[1 - dir].high; cp++){
                double numWire = getWireUsage(layerIdx, gridline, cp) + getFixedUsage(layerIdx, gridline, cp);
                double overflow = max(0.0, numWire - getNumTracks(layerIdx, gridline));
                if(overflow>0)
                    numVio = numVio + 1;
                    if(overflow>0){
                    log() << "vio layer: " << layerIdx
                          << ", dir: " << dir
                          << ", gridline: " << gridline
                          << ", cp: " << cp
                          << ", numWire: " << numWire
                          << ", getNumTracks(layerIdx, gridline): " << getNumTracks(layerIdx, gridline) << std::endl;

                    if(viol_dict.find(std::make_tuple(layerIdx,gridline,cp)) != viol_dict.end()){
                        auto& nets_set = viol_dict[std::make_tuple(layerIdx,gridline,cp)];
                        nets_set.insert(net.dbNet.idx);
                    }else{
                        std::set<int> nets_set;
                        nets_set.insert(net.dbNet.idx);
                        viol_dict[std::make_tuple(layerIdx,gridline,cp)] = nets_set;
                    }
                }
            } 
    }//end for 

    // for (const auto& guide : guides) {
    //     // log() << "getNumVio guide: " << guide << std::endl;
    //     auto dir = database.getLayerDir(guide.layerIdx);
    //     double usage = 1.0 / (guide[dir].range() + 1);
    //     for (int gridline = guide[dir].low; gridline <= guide[dir].high; gridline++) {
    //         GrEdge tempEdge = (dir == X) ? GrEdge(guide.layerIdx, gridline, guide[Y].low, guide[Y].high)
    //                                      : GrEdge(guide.layerIdx, gridline, guide[X].low, guide[X].high);
    //         auto numVioTmp = getNumVio(tempEdge, hasCommit ? 0 : usage);
    //         if(numVioTmp > 0){
    //             // log() << "violation guide: " << guide << std::endl;
    //             edges_viol.push_back(tempEdge);
    //         }
    //         numVio = numVio + numVioTmp;

    //         // log() << "numVioTmpWireGuide: "  << numVioTmp 
    //         //       << "dir: " << dir << std::endl;


    //     }
    // }



    

    const auto& viaGuides = net.viaRouteGuides;
    for (int g1 = 0; g1 < viaGuides.size(); g1++) {
        for (int g2 = g1 + 1; g2 < viaGuides.size(); g2++) {
            if (abs(viaGuides[g1].layerIdx - viaGuides[g2].layerIdx) != 1) continue;

            auto xIntvl = viaGuides[g1][X].IntersectWith(viaGuides[g2][X]);
            auto yIntvl = viaGuides[g1][Y].IntersectWith(viaGuides[g2][Y]);

            if (xIntvl.IsValid() && yIntvl.IsValid()) {
                double usage = 1.0 / ((xIntvl.range() + 1) * (yIntvl.range() + 1));

                for (int x = xIntvl.low; x <= xIntvl.high; x++)
                    for (int y = yIntvl.low; y <= yIntvl.high; y++){
                        auto grPointTmp = GrPoint(min(viaGuides[g1].layerIdx, viaGuides[g2].layerIdx), x, y);
                        // int numVioTmp = getNumVio(GrPoint(min(viaGuides[g1].layerIdx, viaGuides[g2].layerIdx), x, y),
                                            // hasCommit ? 0 : usage);
                        int numVioTmp = getNumVio(grPointTmp,0);// hasCommit ? 0 : usage);
                        if(numVioTmp > 0){
                            
                            vias_viol.push_back(grPointTmp);
                        }
                        numVio = numVio + numVioTmp;
                        // log() << "numVioTmpViaGuide: "  << numVioTmp << std::endl;
                    }
                        
            }
        }
    }





    return (numVio>0)? true : false;

}//eng getVioReport

bool GrRouteGrid::getVioReportRelax(const GrNet& net
                           , std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
                           , std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls 
                           , vector<gr::GrBoxOnLayer>& wireRouteGuides
                           , vector<gr::GrBoxOnLayer>& viaRouteGuides) const
{

    double numVio = 0;


    // log() << "getNumVio net: " << net.getName()
    //       << ", wl: " << net.getWirelength() <<  std::endl;

    const auto& guides = wireRouteGuides;
    for (const auto& guide : guides) {
        // log() << "getNumVio guide: " << guide << std::endl;
        auto dir = database.getLayerDir(guide.layerIdx);
        double usage = 1.0 / (guide[dir].range() + 1);
        for (int gridline = guide[dir].low; gridline <= guide[dir].high; gridline++) {
            GrEdge tempEdge = (dir == X) ? GrEdge(guide.layerIdx, gridline, guide[Y].low, guide[Y].high)
                                         : GrEdge(guide.layerIdx, gridline, guide[X].low, guide[X].high);
            auto numVioTmp = getNumVioRelax(tempEdge,wireMap_tpls,viaMap_tpls, usage);
            // if(numVioTmp > 0){
            //     // log() << "violation guide: " << guide << std::endl;
            //     edges_viol.push_back(tempEdge);
            // }
            numVio = numVio + numVioTmp;

            // log() << "numVioTmpWireGuide: "  << numVioTmp 
            //       << "dir: " << dir << std::endl;


        }
    }



    

    // const auto& viaGuides = net.viaRouteGuides;
    // for (int g1 = 0; g1 < viaGuides.size(); g1++) {
    //     for (int g2 = g1 + 1; g2 < viaGuides.size(); g2++) {
    //         if (abs(viaGuides[g1].layerIdx - viaGuides[g2].layerIdx) != 1) continue;

    //         auto xIntvl = viaGuides[g1][X].IntersectWith(viaGuides[g2][X]);
    //         auto yIntvl = viaGuides[g1][Y].IntersectWith(viaGuides[g2][Y]);

    //         if (xIntvl.IsValid() && yIntvl.IsValid()) {
    //             double usage = 1.0 / ((xIntvl.range() + 1) * (yIntvl.range() + 1));

    //             for (int x = xIntvl.low; x <= xIntvl.high; x++)
    //                 for (int y = yIntvl.low; y <= yIntvl.high; y++){
    //                     auto grPointTmp = GrPoint(min(viaGuides[g1].layerIdx, viaGuides[g2].layerIdx), x, y);
    //                     // int numVioTmp = getNumVio(GrPoint(min(viaGuides[g1].layerIdx, viaGuides[g2].layerIdx), x, y),
    //                                         // hasCommit ? 0 : usage);
    //                     int numVioTmp = getNumVio(grPointTmp, hasCommit ? 0 : usage);
    //                     if(numVioTmp > 0){
                            
    //                         vias_viol.push_back(grPointTmp);
    //                     }
    //                     numVio = numVio + numVioTmp;
    //                     // log() << "numVioTmpViaGuide: "  << numVioTmp << std::endl;
    //                 }
                        
    //         }
    //     }
    // }





    return (numVio>0)? true : false;

}//end getVioReportRelax



double GrRouteGrid::getNumVio(const GrEdge& edge, double selfUsage) const {
    int layerIdx = edge.getLayerIdx();
    auto dir = database.getLayerDir(layerIdx);
    int gridline = edge.u[dir];
    // log() << "getNumVio Edge: " << std::endl;
    // log() << "edge_layerIdx: " << layerIdx 
    //       << ", dir: " << dir 
    //       << ", gridline: " << gridline << std::endl; 

    // log() << "edge.u: " << edge.u << std::endl;
    // log() << "edge.v: " << edge.v << std::endl;


    double numVio = 0;

    for (int i = edge.u[1 - dir]; i < edge.v[1 - dir]; i++) {
        GrEdge tempEdge(layerIdx, gridline, i);
        numVio += max(
            0.0,
            getWireUsage(tempEdge) + selfUsage + getFixedUsage(tempEdge) +
                sqrt((getInCellViaNum(tempEdge.u) + getInCellViaNum(tempEdge.v)) / 2) * db::setting.unitSqrtViaUsage -
                getWireCapacity(tempEdge));

        // log() << " tempEdge_layerIdx: " << layerIdx 
        //         << ", dir: " << i 
        //         << ", gridline: " << gridline << std::endl; 
        // log() << " getWireUsage(tempEdge): " << getWireUsage(tempEdge)
        //       << ", selfUsage: " << selfUsage
        //       << ", getFixedUsage(tempEdge): " << getFixedUsage(tempEdge)
        //       << ", sqrt(tempEdge.u) +(tempEdge.v)) / 2) * db::setting.unitSqrtViaUsage: " << sqrt((getInCellViaNum(tempEdge.u) + getInCellViaNum(tempEdge.v)) / 2) * db::setting.unitSqrtViaUsage
        //       << ", getWireCapacity(tempEdge): " << getWireCapacity(tempEdge) << std::endl;
        // log() << "db::setting.unitSqrtViaUsage: " << db::setting.unitSqrtViaUsage << std::endl;              
        // log() << "getInCellViaNum(tempEdge.u): " << getInCellViaNum(tempEdge.u) << std::endl;
        // log() << "getInCellViaNum(tempEdge.v): " << getInCellViaNum(tempEdge.v) << std::endl;
        // log() << "numVio: " << numVio << std::endl;
    }
    
    return numVio;
}

double GrRouteGrid::getNumVioRelax(const GrEdge& edge,std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
    ,std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls, double selfUsage) const {
    int layerIdx = edge.getLayerIdx();
    auto dir = database.getLayerDir(layerIdx);
    int gridline = edge.u[dir];
    // log() << "getNumVio Edge: " << std::endl;
    // log() << "edge_layerIdx: " << layerIdx 
    //       << ", dir: " << dir 
    //       << ", gridline: " << gridline << std::endl; 

    // log() << "edge.u: " << edge.u << std::endl;
    // log() << "edge.v: " << edge.v << std::endl;


    double numVio = 0;

    for (int i = edge.u[1 - dir]; i < edge.v[1 - dir]; i++) {
        GrEdge tempEdge(layerIdx, gridline, i);
        numVio += max(
            0.0,
            getWireUsageRelax(tempEdge,wireMap_tpls) + selfUsage + getFixedUsage(tempEdge) +
                sqrt((getInCellViaNumRelax(tempEdge.u,viaMap_tpls) + getInCellViaNumRelax(tempEdge.v,viaMap_tpls)) / 2) * db::setting.unitSqrtViaUsage -
                getWireCapacity(tempEdge));

        // numVio += max(
        //     0.0,
        //     getWireUsage(tempEdge) + selfUsage + getFixedUsage(tempEdge) - getWireCapacity(tempEdge));
        // numVio *= 0;
        // log() << " tempEdge_layerIdx: " << layerIdx 
        //         << ", dir: " << i 
        //         << ", gridline: " << gridline << std::endl; 
        // log() << " getWireUsage(tempEdge): " << getWireUsage(tempEdge)
        //       << ", selfUsage: " << selfUsage
        //       << ", getFixedUsage(tempEdge): " << getFixedUsage(tempEdge)
        //       << ", sqrt(tempEdge.u) +(tempEdge.v)) / 2) * db::setting.unitSqrtViaUsage: " << sqrt((getInCellViaNum(tempEdge.u) + getInCellViaNum(tempEdge.v)) / 2) * db::setting.unitSqrtViaUsage
        //       << ", getWireCapacity(tempEdge): " << getWireCapacity(tempEdge) << std::endl;
        // log() << "db::setting.unitSqrtViaUsage: " << db::setting.unitSqrtViaUsage << std::endl;              
        // log() << "getInCellViaNum(tempEdge.u): " << getInCellViaNum(tempEdge.u) << std::endl;
        // log() << "getInCellViaNum(tempEdge.v): " << getInCellViaNum(tempEdge.v) << std::endl;
        // log() << "numVio: " << numVio << std::endl;
    }
    
    return numVio;
}

double GrRouteGrid::getNumVio(const GrPoint& via, double selfUsage) const {
    double via_usage = getViaUsage(via) + selfUsage;
    double numVio = 0;
    for (int side = 0; side <= 1; side++) {
        auto side_point = GrPoint({via.layerIdx + side, via.x, via.y});
        auto incell_area = getInCellArea(side_point);
        auto incell_used_area = getInCellUsedArea(side_point);
        if (incell_area >= incell_used_area) {  // have remaining area

        } else {  // have no remaining area
            numVio += via_usage;
        }
    }
    // if (numVio != 0)
    // log() << numVio << std::endl;
    return numVio;
}

void GrRouteGrid::setViaCapDiscount(double discount) {
    if(db::setting.debug){
        printflog("viaCapDiscount change: %.2f->%.2f\n", viaCapDiscount, discount);
    }
    viaCapDiscount = discount;
}

void GrRouteGrid::setWireCapDiscount(double discount) {
    if(db::setting.debug){
        printflog("wireCapDiscount change: %.2f->%.2f\n", wireCapDiscount, discount);
    }
    wireCapDiscount = discount;
}

void GrRouteGrid::setUnitViaMultiplier(double multiplier) {
    if(db::setting.debug){
        printflog("unitViaMultiplier change: %.2f->%.2f\n", unitViaMultiplier, multiplier);
    }
    unitViaMultiplier = multiplier;
}

void GrRouteGrid::setLogisticSlope(double slope) {
    if(db::setting.debug){
        printflog("logisticSlope change: %.2f->%.2f\n", logisticSlope, slope);
    }
    logisticSlope = slope;
}

void GrRouteGrid::addHistCost() {
    if(db::setting.debug){
        if (db::setting.dbVerbose >= +db::VerboseLevelT::MIDDLE) {
            printlog("Add hist cost");
        }
        log() << "add hist cost " << std::endl;
    }
    for (int layerIdx = 0; layerIdx < database.getLayerNum(); ++layerIdx) {
        if(db::setting.debug){
            log() << "layerIdx: " << layerIdx << std::endl;
        }
        auto dir = database.getLayerDir(layerIdx);
        for (int g = 0; g < getNumGrPoint(dir); ++g) {
            for (int cp = 0; cp < getNumGrEdge(layerIdx); ++cp) {
                GrEdge tempEdge(layerIdx, g, cp);
                if (hasVio(tempEdge)) useHistWire(layerIdx, g, cp, 1);
            }
        }
    }
}

void GrRouteGrid::useHistWire(int layerIdx, int gridline, int cp, double usage) {
    histWireUsageMap[layerIdx][gridline][cp] += usage;
}

void GrRouteGrid::fadeHistCost() {
    if (db::setting.dbVerbose >= +db::VerboseLevelT::MIDDLE) {
        printlog("Fade hist cost by", db::setting.rrrFadeCoeff, "...");
    }
    for (int layerIdx = 0; layerIdx < database.getLayerNum(); ++layerIdx) {
        auto dir = database.getLayerDir(layerIdx);
        // wire
        for (int g = 0; g < getNumGrPoint(dir); ++g)
            for (int cp = 0; cp < getNumGrEdge(layerIdx); ++cp)
                histWireUsageMap[layerIdx][g][cp] *= db::setting.rrrFadeCoeff;
    }
}

void GrRouteGrid::statHistCost() const {
    if (db::setting.dbVerbose >= +db::VerboseLevelT::MIDDLE) {
        std::map<db::CostT, int> histWireUsage, histViaUsage;
        for (int layerIdx = 0; layerIdx < database.getLayerNum(); ++layerIdx) {
            auto dir = database.getLayerDir(layerIdx);
            // wire
            for (int g = 0; g < getNumGrPoint(dir); ++g)
                for (int cp = 0; cp < getNumGrEdge(layerIdx); ++cp) ++histWireUsage[histWireUsageMap[layerIdx][g][cp]];
            // via
        }
        printlog("Hist wire usage is", histWireUsage);
        printlog("Hist via usage is", histViaUsage);
    }
}

double GrRouteGrid::getNetCongestion(const GrNet& net) const{
    if(db::setting.debug){
        log() << "hey this is a net congestion...." << std::endl;
        log() << "net: "<< net.getName() << std::endl;
    }
    double netCongestion = 0;
    double viaCongestion = 0;
    double patchCongestion = 0;

    // log() << "getNumVio net: " << net.getName()
    //       << ", wl: " << net.getWirelength() <<  std::endl;

    const auto& guides = net.wireRouteGuides;
    std::vector<double> edgeCongestions;
    std::vector<double> viaCongestions;
    for (const auto& guide : guides) {
        if(db::setting.debug){
            log() << "getWire guide: " << guide << std::endl;
        }
        auto dir = database.getLayerDir(guide.layerIdx);
        double usage = 1.0 / (guide[dir].range() + 1);
        for (int gridline = guide[dir].low; gridline <= guide[dir].high; gridline++) {
            GrEdge tempEdge = (dir == X) ? GrEdge(guide.layerIdx, gridline, guide[Y].low, guide[Y].high)
                                         : GrEdge(guide.layerIdx, gridline, guide[X].low, guide[X].high);
            auto edgeCongestion = getEdgeCongestion(tempEdge);
            // netCongestion = netCongestion + edgeCongestion;
            edgeCongestions.push_back(edgeCongestion);
        }
    }


    const auto& viaGuides = net.viaRouteGuides;
    if(db::setting.debug){
        log() << "net.viaRouteGuides.size(): " << net.viaRouteGuides.size() << std::endl;
    }
    for (int g1 = 0; g1 < viaGuides.size(); g1++) {
        for (int g2 = g1 + 1; g2 < viaGuides.size(); g2++) {
            if (abs(viaGuides[g1].layerIdx - viaGuides[g2].layerIdx) != 1) continue;

            auto xIntvl = viaGuides[g1][X].IntersectWith(viaGuides[g2][X]);
            auto yIntvl = viaGuides[g1][Y].IntersectWith(viaGuides[g2][Y]);

            if (xIntvl.IsValid() && yIntvl.IsValid()) {
                double usage = 1.0 / ((xIntvl.range() + 1) * (yIntvl.range() + 1));

                for (int x = xIntvl.low; x <= xIntvl.high; x++)
                    for (int y = yIntvl.low; y <= yIntvl.high; y++){
                        int viaCongestionTmp = getGrPointCongestion(GrPoint(min(viaGuides[g1].layerIdx, viaGuides[g2].layerIdx)
                                                                    , x, y));
                        viaCongestions.push_back(viaCongestionTmp);
                        // numVio = numVio + numVioTmp;
                        // log() << "numVioTmpViaGuide: "  << numVioTmp << std::endl;
                    }
                        
            }
        }
    }

    // patchRouteGuides
    const auto& pathGuides = net.patchRouteGuides;
    std::vector<double> patchCongestions;
    for (const auto& guide : pathGuides) {
        if(db::setting.debug){
            log() << "getpatchGuide guide: " << guide << std::endl;
        }
        auto dir = database.getLayerDir(guide.layerIdx);
        double usage = 1.0 / (guide[dir].range() + 1);
        for (int gridline = guide[dir].low; gridline <= guide[dir].high; gridline++) {
            GrEdge tempEdge = (dir == X) ? GrEdge(guide.layerIdx, gridline, guide[Y].low, guide[Y].high)
                                         : GrEdge(guide.layerIdx, gridline, guide[X].low, guide[X].high);
            auto edgeCongestion = getEdgeCongestion(tempEdge);
            // netCongestion = netCongestion + edgeCongestion;
            patchCongestions.push_back(edgeCongestion);
        }
    }

    if(edgeCongestions.size() >= 1)
        netCongestion = *std::max_element(edgeCongestions.begin(), edgeCongestions.end(),
                [] (double lhs, double rhs) {
                return lhs < rhs;
        });

    if(patchCongestions.size() >= 1)
        patchCongestion = *std::max_element(patchCongestions.begin(), patchCongestions.end(),
                [] (double lhs, double rhs) {
                return lhs < rhs;
        });

    if(viaCongestions.size() >= 1)
        viaCongestion = *std::max_element(viaCongestions.begin(), viaCongestions.end(),
                [] (double lhs, double rhs) {
                return lhs < rhs;
        });

    if(db::setting.debug){
        log() << "net: "<< net.getName() 
          << ", netCongestion: " << netCongestion
          << ", viaCongestion: " << viaCongestion
          << ", patchCongestion: " << patchCongestion << std::endl;
    }

    return netCongestion;

}
double GrRouteGrid::getEdgeCongestion(const GrEdge& edge) const{
    int layerIdx = edge.getLayerIdx();
    auto dir = database.getLayerDir(layerIdx);
    int gridline = edge.u[dir];
    double selfUsage = 1;
    // log() << "getNumVio Edge: " << std::endl;
    // log() << "edge_layerIdx: " << layerIdx 
    //       << ", dir: " << dir 
    //       << ", gridline: " << gridline << std::endl; 

    // log() << "edge.u: " << edge.u << std::endl;
    // log() << "edge.v: " << edge.v << std::endl;


    double edgeCongestion = 0;
    std::vector<double> edgeCongestions;

    for (int i = edge.u[1 - dir]; i < edge.v[1 - dir]; i++) {
        GrEdge tempEdge(layerIdx, gridline, i);
        // auto edgeCongestionTmp = (getWireUsage(tempEdge) + selfUsage + getFixedUsage(tempEdge) +
        //         sqrt((getInCellViaNum(tempEdge.u) + getInCellViaNum(tempEdge.v)) / 2) * db::setting.unitSqrtViaUsage) /
        //         getWireCapacity(tempEdge);
         auto edgeCongestionTmp = (getWireUsage(tempEdge) + selfUsage + getFixedUsage(tempEdge)) /
                getWireCapacity(tempEdge);
        edgeCongestions.push_back(edgeCongestionTmp);        

        // log() << " tempEdge_layerIdx: " << layerIdx 
        //         << ", dir: " << i 
        //         << ", gridline: " << gridline << std::endl; 
        // log() << "edge_sum_usage: " << (getWireUsage(tempEdge) + selfUsage + getFixedUsage(tempEdge)) << std::endl;
        // log() << "getWireCapacity(tempEdge): " << getWireCapacity(tempEdge) << std::endl;
        // log() << " getWireUsage(tempEdge): " << getWireUsage(tempEdge)
        //       << ", selfUsage: " << selfUsage
        //       << ", getFixedUsage(tempEdge): " << getFixedUsage(tempEdge)
        //       << ", sqrt(tempEdge.u) +(tempEdge.v)) / 2) * db::setting.unitSqrtViaUsage: " << sqrt((getInCellViaNum(tempEdge.u) + getInCellViaNum(tempEdge.v)) / 2) * db::setting.unitSqrtViaUsage
        //       << ", getWireCapacity(tempEdge): " << getWireCapacity(tempEdge) << std::endl;
        // log() << "db::setting.unitSqrtViaUsage: " << db::setting.unitSqrtViaUsage << std::endl;              
        // log() << "getInCellViaNum(tempEdge.u): " << getInCellViaNum(tempEdge.u) << std::endl;
        // log() << "getInCellViaNum(tempEdge.v): " << getInCellViaNum(tempEdge.v) << std::endl;
        // log() << "numVio: " << numVio << std::endl;
    }
    if(edgeCongestions.size() >= 1)
        edgeCongestion = *std::max_element(edgeCongestions.begin(), edgeCongestions.end(),
                [] (double lhs, double rhs) {
                return lhs < rhs;
        });
    
    return edgeCongestion;

}//end getEdgeCongestion




double GrRouteGrid::getGrPointCongestion(const GrPoint& via) const {
    int selfUsage = 1;
    double via_usage = getViaUsage(via) + selfUsage;
    double grPointCongestion = 0;

    std::vector<double> grPointCongestions;
    for (int side = 0; side <= 1; side++) {
        auto side_point = GrPoint({via.layerIdx + side, via.x, via.y});
        auto incell_area = getInCellArea(side_point);
        auto incell_used_area = getInCellUsedArea(side_point);
        // if (incell_area >= incell_used_area) {  // have remaining area

        // } else {  // have no remaining area
        // //  numVio += via_usage;
        // //     grPointCongestionTmp = via_usage;
        grPointCongestions.push_back(via_usage);
        // }
        // log() << "via_usage: " << via_usage << std::endl;
    }

    if(grPointCongestions.size() >= 1)
        grPointCongestion = *std::max_element(grPointCongestions.begin(), grPointCongestions.end(),
                [] (double lhs, double rhs) {
                return lhs < rhs;
        });
    // if (numVio != 0)
    // log() << numVio << std::endl;
    // log() << "grPointCongestion: " << grPointCongestion << std::endl;
    return grPointCongestion;
}



void GrRouteGrid2D::init2DMaps(const GrRouteGrid& routeGrid) {
    int numLayers = database.getLayerNum();
    wireUsageMap2D.clear();
    fixedMetalMap2D.clear();
    capacityMap2D.clear();
    if(db::setting.debug){
        log() << "init2DMaps " << std::endl;
    }

    // log() << "before resize(2) wireUsageMap2D: "  << wireUsageMap2D.size() << std::endl;
    // log() << "before resize(2) fixedMetalMap2D: " << fixedMetalMap2D.size() << std::endl;
    // log() << "before resize(2) capacityMap2D: " << capacityMap2D.size() << std::endl;

    wireUsageMap2D.resize(2);
    fixedMetalMap2D.resize(2);
    capacityMap2D.resize(2);

    // log() << "after resize(2) wireUsageMap2D: "  << wireUsageMap2D.size() << std::endl;
    // log() << "after resize(2) fixedMetalMap2D: " << fixedMetalMap2D.size() << std::endl;
    // log() << "after resize(2) capacityMap2D: " << capacityMap2D.size() << std::endl;

    int xNumGrPoint = routeGrid.getNumGrPoint(X);
    int yNumGrPoint = routeGrid.getNumGrPoint(Y);
    int xNumGrEdge = yNumGrPoint - 1;
    int yNumGrEdge = xNumGrPoint - 1;
    wireUsageMap2D[X].resize(xNumGrPoint, vector<double>(xNumGrEdge, 0));
    fixedMetalMap2D[X].resize(xNumGrPoint, vector<double>(xNumGrEdge, 0));
    capacityMap2D[X].resize(xNumGrPoint, vector<double>(xNumGrEdge, 0));

    wireUsageMap2D[Y].resize(yNumGrPoint, vector<double>(yNumGrEdge, 0));
    fixedMetalMap2D[Y].resize(yNumGrPoint, vector<double>(yNumGrEdge, 0));
    capacityMap2D[Y].resize(yNumGrPoint, vector<double>(yNumGrEdge, 0));

    for (int l = 0; l < numLayers; l++) {
        auto dir = database.getLayerDir(l);
        if (dir == X) {
            for (int gridline = 0; gridline < xNumGrPoint; gridline++) {
                for (int cp = 0; cp < xNumGrEdge; cp++) {
                    GrEdge tempEdge(l, gridline, cp);
                    fixedMetalMap2D[X][gridline][cp] += routeGrid.getFixedUsage(tempEdge);
                    capacityMap2D[X][gridline][cp] += routeGrid.getWireCapacity(tempEdge);
                    
                    // log() << "xNumGrEdge tempEdge: " << ", l: " << l << ", gridline: "
                    //       << gridline << ", cp: "<< cp
                    //       << ", fixedUsage: " << routeGrid.getFixedUsage(tempEdge) 
                    //       << ", wireCapacity: " << routeGrid.getWireCapacity(tempEdge) << std::endl;
                }
            }
        } else {
            for (int gridline = 0; gridline < yNumGrPoint; gridline++) {
                for (int cp = 0; cp < yNumGrEdge; cp++) {
                    GrEdge tempEdge(l, gridline, cp);
                    fixedMetalMap2D[Y][gridline][cp] += routeGrid.getFixedUsage(tempEdge);
                    capacityMap2D[Y][gridline][cp] += routeGrid.getWireCapacity(tempEdge);

                    // log() << "yNumGrEdge tempEdge: " << ", l: " << l << ", gridline: "
                    //       << gridline << ", cp: "<< cp
                    //       << ", fixedUsage: " << routeGrid.getFixedUsage(tempEdge) 
                    //       << ", wireCapacity: " << routeGrid.getWireCapacity(tempEdge) << std::endl;
                }
            }
        }
    }
}

void GrRouteGrid2D::useWire2D(int dir, int gridline, int cp, double usage) {
    wireUsageMap2D[dir][gridline][cp] += usage;
}
void GrRouteGrid2D::removeUsage2D(int dir, int gridline, int cp, double usage) {
    wireUsageMap2D[dir][gridline][cp] -= usage;
}
double GrRouteGrid2D::getCost2D(int dir, int gridline, int cp) const {
    double cost = 0;

    
    // log() << "fixedMetalMap2D: " << fixedMetalMap2D.size() << std::endl;
    // log() << "wireUsageMap2D: " << wireUsageMap2D.size() << std::endl;
    // log() << "capacityMap2D: " << capacityMap2D.size() << std::endl;

    double fixed_usage = fixedMetalMap2D[dir][gridline][cp];
    double wire_usage = wireUsageMap2D[dir][gridline][cp];
    double cap = capacityMap2D[dir][gridline][cp];

    // log() << "fixedMetalMap2D: " << fixed_usage << std::endl;
    // log() << "wireUsageMap2D: " << wire_usage << std::endl;
    // log() << "capacityMap2D: " << cap << std::endl;

    double demand = fixed_usage + wire_usage;
    // log() << "demand: " << demand << std::endl;
    cost += 1 / (1.0 + exp(-1 * grDatabase.getLogisticSlope() * (demand - cap)));
    // log() << "grDatabase.getLogisticSlope(): " << grDatabase.getLogisticSlope() << std::endl;
    // log() << "cost: " << cost << std::endl;
    return cost;
}


void GrRouteGrid2D::logWireUsage2D(std::string log_name,const GrRouteGrid& routeGrid){
 

 

    std::stringstream ss;
    ss << "layer,gridLine,cp,wireUsageMap2D,fixedMetalMap2D,capacityMap2D,costMap2D" << std::endl;
    // ss << "layer,gridLine,cp" << std::endl;

    



    int numLayers = database.getLayerNum();
    int xNumGrPoint = routeGrid.getNumGrPoint(X);
    int yNumGrPoint = routeGrid.getNumGrPoint(Y);
    int xNumGrEdge = yNumGrPoint - 1;
    int yNumGrEdge = xNumGrPoint - 1;



    for (int l = 0; l < numLayers; l++) {
        auto dir = database.getLayerDir(l);
        if (l==2) 
            break;

        if (dir == X) {
            for (int gridline = 0; gridline < xNumGrPoint; gridline++) {
                for (int cp = 0; cp < xNumGrEdge; cp++) {                   
                    ss << std::to_string(1)
                       << "," << std::to_string(gridline)
                       << "," << std::to_string(cp)
                       << "," << std::to_string(wireUsageMap2D[X][gridline][cp])
                       << "," << std::to_string(fixedMetalMap2D[X][gridline][cp])
                       << "," << std::to_string(capacityMap2D[X][gridline][cp]) 
                       << "," << std::to_string(getCost2D(X,gridline,cp))
                       << std::endl;
                }
            }
        } 
        else {
            for (int gridline = 0; gridline < yNumGrPoint; gridline++) {
                for (int cp = 0; cp < yNumGrEdge; cp++) {
                    ss << std::to_string(0)
                       << "," << std::to_string(gridline)
                       << "," << std::to_string(cp)
                       << "," << std::to_string(wireUsageMap2D[Y][gridline][cp])
                       << "," << std::to_string(fixedMetalMap2D[Y][gridline][cp] )
                       << "," << std::to_string(capacityMap2D[Y][gridline][cp] ) 
                       << "," << std::to_string(getCost2D(Y,gridline,cp))
                       << std::endl;
                }
            }
        }
        
    }


    std::ofstream fout(log_name);
    fout << ss.str();
    fout.close();


    
}//end logWireUsage2D



}  // namespace gr
