#include "InitRouteCell.h"

using namespace std;

namespace routeCell{

mutex fluteMutex;

struct hash_tuple_cell {  // hash binary tuple
    template <class T>
    size_t operator()(const tuple<T, T>& tup) const {
        auto hash1 = hash<T>{}(get<0>(tup));
        auto hash2 = hash<T>{}(get<1>(tup));
        return hash1 ^ hash2;
    }
};

void InitRouteCell::runFluteV2(){
    // if(grNet.getName() == "net275" )debug = false;
    std::vector<std::pair<std::tuple<int,int,int>,double>> wireMap_tpls;
    std::vector<std::pair<std::tuple<int,int,int>,double>> viaMap_tpls;
    grDatabase.removeNetRelax(grNet,wireMap_tpls,viaMap_tpls);

    // if(grNet.getName()=="net121")
    //     debug = true;
    if(debug){
        log() << "grNet: " << grNet.getName() << ", initRoute ..." << std::endl;
    }

    getGrPinAccessBoxes(*&grDatabase);

    // if(grNet.getName() == "net275" )debug = false;
    if(debug){
        log() << "grNet: " << grNet.getName() << ", initRoute ..." << std::endl;
    }
    // get net center
    // float net_ctrz = 0;
    float net_ctrx = 0;
    float net_ctry = 0;


    

    // log() << "runFlute " << grNet.getName() << std::endl;
    int i = 0;
    for (auto& pinBoxes : grPinAccessBoxes) {
        // float pin_ctrz = 0;
        float pin_ctrx = 0;
        float pin_ctry = 0;
        for (auto& pinBox : pinBoxes) {
            pin_ctrx += pinBox.x;
            pin_ctry += pinBox.y;
            // log() << "pinBox.x: " << pinBox.x
            //       << ", pinBox.y: " << pinBox.y << std::endl;
        }
        pin_ctrx /= pinBoxes.size();
        pin_ctry /= pinBoxes.size();

        net_ctrx += pin_ctrx;
        net_ctry += pin_ctry;
    }
    net_ctrx /= grPinAccessBoxes.size();
    net_ctry /= grPinAccessBoxes.size();


    // log() << "net_ctrx: " << net_ctrx
    //       << ", net_ctry: " << net_ctry << std::endl;

    // get pin center
    vector<tuple<int, int>> pinCenters;

    for (auto& pinBoxes : grPinAccessBoxes) {
        float xCenter = 0;
        float yCenter = 0;
        for (auto& pinBox : pinBoxes) {
            xCenter += pinBox.x;
            yCenter += pinBox.y;
            if(debug)
                log() << "i: " << i 
                      << ", pinBox.x: " << pinBox.x << ", pinBox.y: " << pinBox.y << std::endl;
            // log() << "xCenter: " << xCenter << ", yCenter: " << yCenter << std::endl;
        }
        // log() << "pinBoxes.size(): " << pinBoxes.size() << std::endl;
        xCenter /= pinBoxes.size();
        yCenter /= pinBoxes.size();

        // for(auto tmp : pinCenters){
        // log() << "xCenter: " << xCenter << ", yCenter: " << yCenter << std::endl;
        // }
        // 3 kinds of accessibility (0) totally vio-free (1) one side vio-free (2) no side vio-free
        float min_dist[3];
        min_dist[0] = min_dist[1] = min_dist[2] = numeric_limits<float>::max();
        int best_box[3] = {-1, -1, -1};

        for (int pb = 0; pb < pinBoxes.size(); pb++) {
            // check pin box's accessibility
            auto& pb_x = pinBoxes[pb].x;
            auto& pb_y = pinBoxes[pb].y;
            int layer_idx = pinBoxes[pb].layerIdx;
            if(debug)
                log() << "pinsBox: [x: " <<  pb_x << ", y: " << pb_y << std::endl;
            int is_x = database.getLayerDir(layer_idx) == X ? 1 : 0;
            auto low_edge =
                gr::GrEdge({layer_idx, max(pb_x - (1 - is_x), 0), max(pb_y - is_x, 0)}, {layer_idx, pb_x, pb_y});
            auto high_edge = gr::GrEdge({layer_idx, pb_x, pb_y},
                                        {layer_idx,
                                         min(pb_x + (1 - is_x), grDatabase.getNumGrPoint(X) - 1),
                                         min(pb_y + is_x, grDatabase.getNumGrPoint(Y) - 1)});
            if(debug){
                log() << "low_edge: " << low_edge << std::endl;
                log() << "low_edge has viol: " << grDatabase.hasVio(low_edge, false)<< std::endl;
                log() << "high_edge: " << high_edge << std::endl;
                log() << "high_edge has viol: " << grDatabase.hasVio(high_edge, false)<< std::endl;
            }

            int pb_access = 0;
            if(relax){
                pb_access += int(grDatabase.hasVioRelax(low_edge,wireMap_tpls,viaMap_tpls, false));
                pb_access += int(grDatabase.hasVioRelax(high_edge,wireMap_tpls,viaMap_tpls, false));
            }else{
                pb_access += int(grDatabase.hasVio(low_edge, false));
                pb_access += int(grDatabase.hasVio(high_edge, false));
            }
             
            float dist = abs(pinBoxes[pb].x - net_ctrx) + abs(pinBoxes[pb].y - net_ctry);

            if(debug){
                log() << "pb: " << pb 
                      << ", pinBoxes[pb].x: " << pinBoxes[pb].x 
                      << ", pinBoxes[pb].y: " << pinBoxes[pb].y 
                      << ", net_ctrx: "  << net_ctrx
                      << ", net_ctry: "  << net_ctry 
                      << ", dist: " << dist << std::endl;

            }
            if (dist < min_dist[pb_access]) {
                min_dist[pb_access] = dist;
                best_box[pb_access] = pb;
            }
        }
        for (int ac = 0; ac < 3; ac++) {
            if (best_box[ac] != -1) {
                pinCenters.emplace_back(make_tuple(pinBoxes[best_box[ac]].x, pinBoxes[best_box[ac]].y));
                break;
            }
        }
        i++;
    }

    if(debug){
        log() << "pinCenters" << std::endl;
        for(auto tmp : pinCenters){
            log() << "x: " << std::get<0>(tmp) << ", y: " << std::get<1>(tmp) << std::endl;
        }
    }
    

    // location to pins
    unordered_map<tuple<int, int>, vector<int>, hash_tuple_cell> loc2Pins;
    for (int i = 0; i < pinCenters.size(); i++) {
        loc2Pins[pinCenters[i]].emplace_back(i);
    }

    // log() << "loc2Pins size: " << loc2Pins.size() << std::endl;

    // flute
    int degree = loc2Pins.size();
    int xs[100 * degree];
    int ys[100 * degree];
    int pt_cnt = 0;
    int node_cnt = 0;
    // if (loc2Pins.size() > MAX_DEGREE) {
    //     log() << "Warning: Net degree larger than MAXD(flute)" << endl;
    // }

    for (auto& loc2pin : loc2Pins) {
        const tuple<int, int>& loc = loc2pin.first;
        xs[pt_cnt] = get<0>(loc);
        ys[pt_cnt] = get<1>(loc);
        if(debug)
            log() << "x: " << get<0>(loc) << ", y: " << get<1>(loc) << std::endl;
        pt_cnt++;
    }
    if (degree >= 2) {
        fluteMutex.lock();
        Tree flutetree = flute(degree, xs, ys, ACCURACY);
        fluteMutex.unlock();
        // log() << "flutetree length: " << flutetree.length << std::endl;
        unordered_map<tuple<int, int>, int, hash_tuple_cell> loc2Node;  // location -> RouteNode index
        for (int i = 0; i < degree * 2 - 2; i++) {
            Branch& branch1 = flutetree.branch[i];
            Branch& branch2 = flutetree.branch[branch1.n];
            tuple<int, int> fluteEdge[2];
            fluteEdge[0] = make_tuple(branch1.x, branch1.y);
            fluteEdge[1] = make_tuple(branch2.x, branch2.y);
            // create route nodes
            for (int j = 0; j < 2; j++) {
                tuple<int, int>& nodeLoc = fluteEdge[j];
                if (loc2Node.find(nodeLoc) == loc2Node.end()) {
                    RouteNode node(fluteEdge[j]);
                    // pin info
                    if (loc2Pins.find(nodeLoc) != loc2Pins.end()) {
                        node.pinIdxs = loc2Pins[nodeLoc];
                        for (int pIdx : node.pinIdxs) {
                            int layerIdx = grNet.pinAccessBoxes[pIdx][0].layerIdx;
                            node.pinLayers.emplace_back(layerIdx);

                            // for (auto& box: grNet.pinAccessBoxes[pIdx]) {
                            //     if (box.layerIdx != layerIdx) {
                            //         log() << "Error: Pin across layers spotted" << endl;
                            //     }
                            // }
                        }
                    }
                    node.idx = node_cnt;
                    routeNodes[node_cnt++] = node;
                    loc2Node[nodeLoc] = routeNodes.size() - 1;
                }
            }
            // add connectivity info
            if (fluteEdge[0] != fluteEdge[1]) {
                int nodeIdx1 = loc2Node[fluteEdge[0]];
                int nodeIdx2 = loc2Node[fluteEdge[1]];
                routeNodes[nodeIdx1].toConnect.insert(nodeIdx2);
                routeNodes[nodeIdx2].toConnect.insert(nodeIdx1);
            }
        }
    } else if (degree == 1) {
        auto nodeLoc = make_tuple(xs[0], ys[0]);
        RouteNode node(nodeLoc);
        if (loc2Pins.find(nodeLoc) != loc2Pins.end()) {
            node.pinIdxs = loc2Pins[nodeLoc];
            for (int pIdx : node.pinIdxs) {
                int layerIdx = grNet.pinAccessBoxes[pIdx][0].layerIdx;
                node.pinLayers.emplace_back(layerIdx);
            }
        } else {
            log() << "Error: loc2Pins is empty." << std::endl;
        }
        node.idx = node_cnt;
        routeNodes[node_cnt++] = node;
        return;
    } else {
        log() << "Error: degree = 0." << std::endl;
    }
   

}//end runFluteV2   

void InitRouteCell::runFlute(){
    // get cell 
    // log() << "getFlute "<< std::endl;
    // std::map<int, RouteNode> routeNodes;
    // if(grNet.getName() != "net275") debug = false;
    // if(debug)
    //     log() << "grNet: " << grNet.getName() << ", initRouteCell..." <<std::endl;
    
    vector<vector<gr::GrPoint>> cellAccessBoxs;
    cellAccessBoxs.resize(cellBoxs.size());
    // convert box coordination from database DBU to GlobalGrid coordination
    for(int i = 0; i < cellBoxs.size() ; i++){
        std::unordered_set<gr::GrPoint> pointSet;
        const auto& box = cellBoxs[i];
        auto box_layer = db::BoxOnLayer(0, box);
        auto grBox = grDatabase.rangeSearchGCell(box_layer);
        // if(debug)
        //     log() << "i: " << i 
        //           << ", box: " << box 
        //           << ", grBox: " << grBox << std::endl;
        for (int x = grBox[X].low; x <= grBox[X].high; x++){
            for (int y = grBox[Y].low; y <= grBox[Y].high; y++) {
                pointSet.emplace(box_layer.layerIdx, x, y);
                    // log() << "box: " << box 
                    //       << ", box.layerIdx: " << box_layer.layerIdx
                    //       << ", gcell_x: " << x
                    //       << ", gcell_y: " << y << std::endl;
                            
                }
        }
        for (auto& point : pointSet) {
            // if(debug)
            //     log() << "box: " << box 
            //           << ", pt.x: " << point.x
            //           << ", pt.y: " << point.y << std::endl;
            cellAccessBoxs[i].push_back(point);
        }   
    }
    
    // get net center
    // float net_ctrz = 0;
    float net_ctrx = 0;
    float net_ctry = 0;

    // log() << "runFlute " << grNet.getName() << std::endl;

    for (auto& pinBoxes : cellAccessBoxs) {
        // float pin_ctrz = 0;
        float pin_ctrx = 0;
        float pin_ctry = 0;
        for (auto& pinBox : pinBoxes) {
            pin_ctrx += pinBox.x;
            pin_ctry += pinBox.y;
            // log() << "pinBox.x: " << pinBox.x
            //       << ", pinBox.y: " << pinBox.y << std::endl;
        }
        pin_ctrx /= pinBoxes.size();
        pin_ctry /= pinBoxes.size();

        net_ctrx += pin_ctrx;
        net_ctry += pin_ctry;
    }
    net_ctrx /= cellAccessBoxs.size();
    net_ctry /= cellAccessBoxs.size();


    // log() << "net_ctrx: " << net_ctrx
    //       << ", net_ctry: " << net_ctry << std::endl;

    // get pin center
    vector<tuple<int, int>> pinCenters;
    
    for (auto& pinBoxes : cellAccessBoxs) {
        float xCenter = 0;
        float yCenter = 0;
        for (auto& pinBox : pinBoxes) {
            xCenter += pinBox.x;
            yCenter += pinBox.y;
            // log() << "pinBox.x: " << pinBox.x << ", pinBox.y: " << pinBox.y << std::endl;
            // log() << "xCenter: " << xCenter << ", yCenter: " << yCenter << std::endl;
        }
        // log() << "pinBoxes.size(): " << pinBoxes.size() << std::endl;
        xCenter /= pinBoxes.size();
        yCenter /= pinBoxes.size();

        // for(auto tmp : pinCenters){
        //     log() << "xCenter: " << xCenter << ", yCenter: " << yCenter << std::endl;
        // }
        // 3 kinds of accessibility (0) totally vio-free (1) one side vio-free (2) no side vio-free
        float min_dist[3];
        min_dist[0] = min_dist[1] = min_dist[2] = numeric_limits<float>::max();
        int best_box[3] = {-1, -1, -1};

        for (int pb = 0; pb < pinBoxes.size(); pb++) {
            // check pin box's accessibility
            auto& pb_x = pinBoxes[pb].x;
            auto& pb_y = pinBoxes[pb].y;
            int layer_idx = pinBoxes[pb].layerIdx;
            int is_x = database.getLayerDir(layer_idx) == X ? 1 : 0;
            auto low_edge =
                gr::GrEdge({layer_idx, max(pb_x - (1 - is_x), 0), max(pb_y - is_x, 0)}, {layer_idx, pb_x, pb_y});
            auto high_edge = gr::GrEdge({layer_idx, pb_x, pb_y},
                                        {layer_idx,
                                         min(pb_x + (1 - is_x), grDatabase.getNumGrPoint(X) - 1),
                                         min(pb_y + is_x, grDatabase.getNumGrPoint(Y) - 1)});

            int pb_access = 0;
            pb_access += int(grDatabase.hasVio(low_edge, false));
            pb_access += int(grDatabase.hasVio(high_edge, false));

            float dist = abs(pinBoxes[pb].x - net_ctrx) + abs(pinBoxes[pb].y - net_ctrx);
            if (dist < min_dist[pb_access]) {
                min_dist[pb_access] = dist;
                best_box[pb_access] = pb;
            }
        }
        for (int ac = 0; ac < 3; ac++) {
            if (best_box[ac] != -1) {
                pinCenters.emplace_back(make_tuple(pinBoxes[best_box[ac]].x, pinBoxes[best_box[ac]].y));
                break;
            }
        }
    }

    // log() << "pinCenters" << std::endl;
    // for(auto tmp : pinCenters){
    //     log() << "x: " << std::get<0>(tmp) << ", y: " << std::get<1>(tmp) << std::endl;
    // }

    // location to pins
    unordered_map<tuple<int, int>, vector<int>, hash_tuple_cell> loc2Pins;
    for (int i = 0; i < pinCenters.size(); i++) {
        loc2Pins[pinCenters[i]].emplace_back(i);
    }

    // log() << "loc2Pins size: " << loc2Pins.size() << std::endl;

    // flute
    int degree = loc2Pins.size();
    int xs[100 * degree];
    int ys[100 * degree];
    int pt_cnt = 0;
    int node_cnt = 0;
    // if (loc2Pins.size() > MAX_DEGREE) {
    //     log() << "Warning: Net degree larger than MAXD(flute)" << endl;
    // }
    // log() << "loc2Pins " << std::endl;
    for (auto& loc2pin : loc2Pins) {
        const tuple<int, int>& loc = loc2pin.first;
        xs[pt_cnt] = get<0>(loc);
        ys[pt_cnt] = get<1>(loc);
        // if(debug)
        //     log() << "x: " << get<0>(loc) << ", y: " << get<1>(loc) << std::endl;
        pt_cnt++;
    }
    
    if (degree >= 2) {
        fluteMutex.lock();
        Tree flutetree = flute(degree, xs, ys, ACCURACY);
        // log() << "flutetree wl: " << flutetree.length << ",degree: " << degree << std::endl;
        fluteMutex.unlock();
        
        
        unordered_map<tuple<int, int>, int, hash_tuple_cell> loc2Node;  // location -> RouteNode index
        for (int i = 0; i < degree * 2 - 2; i++) {
            Branch& branch1 = flutetree.branch[i];
            Branch& branch2 = flutetree.branch[branch1.n];
            tuple<int, int> fluteEdge[2];
            fluteEdge[0] = make_tuple(branch1.x, branch1.y);
            fluteEdge[1] = make_tuple(branch2.x, branch2.y);
            // create route nodes
            for (int j = 0; j < 2; j++) {
                tuple<int, int>& nodeLoc = fluteEdge[j];
                if (loc2Node.find(nodeLoc) == loc2Node.end()) {
                    RouteNode node(fluteEdge[j]);
                    // pin info
                    if (loc2Pins.find(nodeLoc) != loc2Pins.end()) {
                        node.pinIdxs = loc2Pins[nodeLoc];
                        for (int pIdx : node.pinIdxs) {
                            int layerIdx = 0;// grNet.pinAccessBoxes[pIdx][0].layerIdx;
                            node.pinLayers.emplace_back(layerIdx);

                            // for (auto& box: grNet.pinAccessBoxes[pIdx]) {
                            //     if (box.layerIdx != layerIdx) {
                            //         log() << "Error: Pin across layers spotted" << endl;
                            //     }
                            // }
                        }
                    }
                    node.idx = node_cnt;
                    routeNodes[node_cnt++] = node;
                    loc2Node[nodeLoc] = routeNodes.size() - 1;
                }
            }
            // add connectivity info
            if (fluteEdge[0] != fluteEdge[1]) {
                int nodeIdx1 = loc2Node[fluteEdge[0]];
                int nodeIdx2 = loc2Node[fluteEdge[1]];
                routeNodes[nodeIdx1].toConnect.insert(nodeIdx2);
                routeNodes[nodeIdx2].toConnect.insert(nodeIdx1);
            }
        }
        // log() << "routeNodes: " << routeNodes.size() << std::endl;
        wl_cost = flutetree.length;
        // return flutetree.length;
    } else if (degree == 1) {
        // return 1;
        auto nodeLoc = make_tuple(xs[0], ys[0]);
        RouteNode node1(nodeLoc);
        if (loc2Pins.find(nodeLoc) != loc2Pins.end()) {
            node1.pinIdxs = loc2Pins[nodeLoc];
            for (int pIdx : node1.pinIdxs) {
                // int layerIdx = grNet.pinAccessBoxes[pIdx][0].layerIdx;
                int layerIdx = 0;
                node1.pinLayers.emplace_back(layerIdx);
            }
        } else {
            log() << "Error: loc2Pins is empty." << std::endl;
        }
        node1.idx = node_cnt;
        routeNodes[node_cnt++] = node1;
        // log() << "routeNodes: " << routeNodes.size() << std::endl;
        // return 1;
        wl_cost = 1;
    } else {
        log() << "Error: degree = 0." << std::endl;
    }
    // wl_cost = 1;
    // return 1;
}//end runFlute

void InitRouteCell::patternRoute() {
    bool debug = false;
    // if(grNet.getName() == database.debug_net) debug = false;
    // add relaxation to router estimator 
    std::vector<std::pair<std::tuple<int,int,int>,double>> wireMap_tpls;
    std::vector<std::pair<std::tuple<int,int,int>,double>> viaMap_tpls;
    for (const auto& guide : grNet.wireRouteGuides) grDatabase.removeWireRelax(guide,wireMap_tpls);
    const auto& viaGuides = grNet.viaRouteGuides;
    for (int g1 = 0; g1 < viaGuides.size(); g1++) {
        for (int g2 = g1 + 1; g2 < viaGuides.size(); g2++) {
            if (abs(viaGuides[g1].layerIdx - viaGuides[g2].layerIdx) != 1) continue;

            auto xIntvl = viaGuides[g1][X].IntersectWith(viaGuides[g2][X]);
            auto yIntvl = viaGuides[g1][Y].IntersectWith(viaGuides[g2][Y]);

            if (xIntvl.IsValid() && yIntvl.IsValid()){
                grDatabase.removeViaRelax({min(viaGuides[g1].layerIdx, viaGuides[g2].layerIdx), xIntvl, yIntvl},viaMap_tpls);
            }
        }
    }




    if(debug)
        log() << "patternRouteRefinePlacement grNet: " << grNet.getName() << std::endl;
    if(debug){
        log() << "grNet: " << grNet.getName() << std::endl;
        log() << "database.getUnitViaCost(): " << database.getUnitViaCost() << std::endl;
        log() << "unitSqrtViaUsage: " << db::setting.unitSqrtViaUsage << std::endl;
        log() << "grDatabase.getLogisticSlope(): " << grDatabase.getLogisticSlope() << std::endl;
        log() << "grDatabase.getUnitViaMultiplier(): " << grDatabase.getUnitViaMultiplier() << std::endl;
    }


    int layer_cnt = database.getLayerNum();
    for (auto& edge : routeEdges) {
        auto& fromNode = routeNodes[edge.from];
        auto& toNode = routeNodes[edge.to];

        if(debug){
            log() << "fromNode idx: " << fromNode.idx
                  << "[x: " << fromNode.x
                  << ",y: " << fromNode.y << "]" << " -> "
                  << "toNode idx: " << toNode.idx
                  << "[x: " << toNode.x
                  << ",y: " << toNode.y << "]" << std::endl;
        }
        // Initialize fromNode's exitCosts
        fromNode.exitCosts.resize(layer_cnt, numeric_limits<db::CostT>::max());
        int lowest_pin = layer_cnt;
        int highest_pin = -1;
        for (int pin_layer : fromNode.pinLayers) {
            lowest_pin = min(lowest_pin, pin_layer);
            highest_pin = max(highest_pin, pin_layer);
        }
        if (fromNode.childIdxs.size() == 0) {  // no children nodes, only via costs
            for (int layer_idx = 0; layer_idx < layer_cnt; layer_idx++) {
                int lowest_layer = min(lowest_pin, layer_idx);
                int highest_layer = max(highest_pin, layer_idx);

                db::CostT exit_cost = grDatabase.getStackViaCostRelax(gr::GrPoint(lowest_layer, fromNode.x, fromNode.y),
                                                                 highest_layer - lowest_layer,wireMap_tpls,viaMap_tpls);
                fromNode.exitCosts[layer_idx] = exit_cost;
            }
        } else {
            // initialize exitEnterLayers map
            for (int child_idx : fromNode.childIdxs) {
                fromNode.exitEnterLayers[child_idx] = vector<int>(layer_cnt, -1);
            }
            // enumerate all children layer combinations
            int num_children = fromNode.childIdxs.size();
            for (int enum_idx = 0; enum_idx < pow(layer_cnt, num_children); enum_idx++) {
                vector<int> enum_vec;
                int e_idx = enum_idx;
                for (int en = 0; en < num_children; en++) {
                    enum_vec.emplace_back(e_idx % layer_cnt);
                    e_idx /= layer_cnt;
                }
                db::CostT prev_cost = 0;
                for (int c = 0; c < num_children; c++) {
                    int child_idx = fromNode.childIdxs[c];
                    prev_cost += fromNode.enterCosts[child_idx][enum_vec[c]];
                }
                int highest_layer = highest_pin;
                int lowest_layer = lowest_pin;
                for (int e_layer : enum_vec) {
                    highest_layer = max(highest_layer, e_layer);
                    lowest_layer = min(lowest_layer, e_layer);
                }
                for (int layer_idx = 0; layer_idx < layer_cnt; layer_idx++) {
                    int via_bottom = min(lowest_layer, layer_idx);
                    int via_height = max(highest_layer, layer_idx) - via_bottom;
                    db::CostT via_cost =
                        grDatabase.getStackViaCostRelax(gr::GrPoint(via_bottom, fromNode.x, fromNode.y), via_height,wireMap_tpls,viaMap_tpls);
                    db::CostT cost = prev_cost + via_cost;
                    if (cost < fromNode.exitCosts[layer_idx]) {
                        fromNode.exitCosts[layer_idx] = cost;
                        for (int c = 0; c < num_children; c++) {
                            fromNode.exitEnterLayers[fromNode.childIdxs[c]][layer_idx] = enum_vec[c];
                        }
                    }
                }
            }
        }
        // Patterns
        if(debug){
            log() << "refineplacement bf fromNode print..." << std::endl;
            fromNode.print();
            log() << "refineplacement bf toNode print..." << std::endl;
            toNode.print();
        }
        LShape(edge,wireMap_tpls,viaMap_tpls);
        if(debug){
            log() << "refineplacement af fromNode print..." << std::endl;
            fromNode.print();
            log() << "refineplacement af toNode print..." << std::endl;
            toNode.print();
        }
    }
}//end patternRoute


void InitRouteCell::patternRouteMemo() {
    bool debug = false;
    // if(grNet.getName() == database.debug_net) debug = false;
    // add relaxation to router estimator 
    std::vector<std::pair<std::tuple<int,int,int>,double>> wireMap_tpls;
    std::vector<std::pair<std::tuple<int,int,int>,double>> viaMap_tpls;
    for (const auto& guide : grNet.wireRouteGuides) grDatabase.removeWireRelax(guide,wireMap_tpls);
    const auto& viaGuides = grNet.viaRouteGuides;
    for (int g1 = 0; g1 < viaGuides.size(); g1++) {
        for (int g2 = g1 + 1; g2 < viaGuides.size(); g2++) {
            if (abs(viaGuides[g1].layerIdx - viaGuides[g2].layerIdx) != 1) continue;

            auto xIntvl = viaGuides[g1][X].IntersectWith(viaGuides[g2][X]);
            auto yIntvl = viaGuides[g1][Y].IntersectWith(viaGuides[g2][Y]);

            if (xIntvl.IsValid() && yIntvl.IsValid()){
                grDatabase.removeViaRelax({min(viaGuides[g1].layerIdx, viaGuides[g2].layerIdx), xIntvl, yIntvl},viaMap_tpls);
            }
        }
    }


    // memo keep the cost of vias stack in hash table
    std::unordered_map<int,std::vector<db::CostT>> viaStack_CostTb;
    std::unordered_map<int,unordered_map<tuple<int, int>,db::CostT , hash_tuple_cell>> costTb_hash;

    int layer_cnt = database.getLayerNum();  


    if(db::setting.patternRouteMemorization){
        for(auto node_pair : routeNodes ){
            int node_idx = node_pair.first;
            auto node = node_pair.second;
            int lowest_pin = layer_cnt;
            int highest_pin = -1;
            bool steiner_node = false;
            if(node.pinLayers.size() > 0){
                for (int pin_layer : node.pinLayers) {
                    lowest_pin = min(lowest_pin, pin_layer);
                    highest_pin = max(highest_pin, pin_layer);
                }
            }else{
                // lowest_pin = 0;
                // highest_pin = 0;
                steiner_node = true;
            }
             steiner_node = true;
            unordered_map<tuple<int, int>,db::CostT , hash_tuple_cell> costTb_tmp;
            if(!steiner_node){
                for (int layer_idx = 0; layer_idx < layer_cnt; layer_idx++) {
                
                        int lowest_layer = min(lowest_pin, layer_idx);
                        int highest_layer = max(highest_pin, layer_idx);
                        int via_bottom = min(lowest_layer, layer_idx);
                        int via_height = max(highest_layer, layer_idx) - via_bottom;
                        auto via_cost = grDatabase.getStackViaCostRelax(
                                         gr::GrPoint(via_bottom, node.x, node.y),\
                                         via_height,wireMap_tpls,viaMap_tpls);
                        viaStack_CostTb[node_idx].push_back(via_cost);

                        costTb_tmp[std::make_pair(via_bottom,via_height)] = via_cost;

                    

                    // if(debug)

                    //     log() << "node_idx: " << node.idx
                    //         << ", lowest_layer: " << lowest_layer
                    //         << ", highest_layer: " << highest_layer << std::endl;
                    // if(debug)                
                    //     log() << "node_idx: " << node_idx << std::endl;
                    
                    if(debug)
                        log() << "node_idx: " << node.idx
                            << ", b: " << via_bottom
                            << ", x: "  << node.x
                            << ", y: " << node.y 
                            << ", h: " << via_height
                            << ", cc: " << via_cost
                            << std::endl;
                    
                }//end for
            }else{
                for(int lowest_layer = 0;lowest_layer< layer_cnt; lowest_layer++){
                    for (int layer_idx = 0; layer_idx < layer_cnt; layer_idx++) {
                        int highest_layer = layer_idx;
                        int via_bottom = min(lowest_layer, layer_idx);
                        int via_height = max(highest_layer, layer_idx) - via_bottom;
                        // log() << "b: " << via_bottom 
                        //       << ", h: " << via_height << std::endl;
                        auto via_cost = grDatabase.getStackViaCostRelax(
                                         gr::GrPoint(via_bottom, node.x, node.y),\
                                         via_height,wireMap_tpls,viaMap_tpls);
                        costTb_tmp[std::make_pair(via_bottom,via_height)] = via_cost;
                        // viaStack_CostTb[node_idx].push_back(via_cost);

                        // costTb_tmp[std::make_pair(via_bottom,via_height)] = via_cost;

                    

                    // if(debug)

                    //     log() << "node_idx: " << node.idx
                    //         << ", lowest_layer: " << lowest_layer
                    //         << ", highest_layer: " << highest_layer << std::endl;
                    // if(debug)                
                    //     log() << "node_idx: " << node_idx << std::endl;
                    
                    if(debug)
                        log() << "snode_idx: " << node.idx
                            << ", b: " << via_bottom
                            << ", x: "  << node.x
                            << ", y: " << node.y 
                            << ", h: " << via_height
                            << ", cc: " << via_cost
                            << std::endl;
                        
                    }//end for
                }//end for out

                
            }
            costTb_hash[node.idx] = costTb_tmp;
        }
    }//end 



    if(debug)
        log() << "patternRouteRefinePlacement grNet: " << grNet.getName() << std::endl;
    if(debug){
        log() << "grNet: " << grNet.getName() << std::endl;
        log() << "database.getUnitViaCost(): " << database.getUnitViaCost() << std::endl;
        log() << "unitSqrtViaUsage: " << db::setting.unitSqrtViaUsage << std::endl;
        log() << "grDatabase.getLogisticSlope(): " << grDatabase.getLogisticSlope() << std::endl;
        log() << "grDatabase.getUnitViaMultiplier(): " << grDatabase.getUnitViaMultiplier() << std::endl;
    }


    // int layer_cnt = database.getLayerNum();
    for (auto& edge : routeEdges) {
        auto& fromNode = routeNodes[edge.from];
        auto& toNode = routeNodes[edge.to];

        if(debug){
            log() << "fromNode idx: " << fromNode.idx
                  << "[x: " << fromNode.x
                  << ",y: " << fromNode.y << "]" << " -> "
                  << "toNode idx: " << toNode.idx
                  << "[x: " << toNode.x
                  << ",y: " << toNode.y << "]" << std::endl;
        }
        // Initialize fromNode's exitCosts
        fromNode.exitCosts.resize(layer_cnt, numeric_limits<db::CostT>::max());
        int lowest_pin = layer_cnt;
        int highest_pin = -1;
        for (int pin_layer : fromNode.pinLayers) {
            lowest_pin = min(lowest_pin, pin_layer);
            highest_pin = max(highest_pin, pin_layer);
        }
        if (fromNode.childIdxs.size() == 0) {  // no children nodes, only via costs
            for (int layer_idx = 0; layer_idx < layer_cnt; layer_idx++) {
                int lowest_layer = min(lowest_pin, layer_idx);
                int highest_layer = max(highest_pin, layer_idx);

                db::CostT exit_cost = grDatabase.getStackViaCostRelax(gr::GrPoint(lowest_layer, fromNode.x, fromNode.y),
                                                                 highest_layer - lowest_layer,wireMap_tpls,viaMap_tpls);
                fromNode.exitCosts[layer_idx] = exit_cost;
            }
        } else {
            // initialize exitEnterLayers map
            for (int child_idx : fromNode.childIdxs) {
                fromNode.exitEnterLayers[child_idx] = vector<int>(layer_cnt, -1);
            }
            // enumerate all children layer combinations
            int num_children = fromNode.childIdxs.size();
            for (int enum_idx = 0; enum_idx < pow(layer_cnt, num_children); enum_idx++) {
                vector<int> enum_vec;
                int e_idx = enum_idx;
                for (int en = 0; en < num_children; en++) {
                    enum_vec.emplace_back(e_idx % layer_cnt);
                    e_idx /= layer_cnt;
                }
                db::CostT prev_cost = 0;
                for (int c = 0; c < num_children; c++) {
                    int child_idx = fromNode.childIdxs[c];
                    prev_cost += fromNode.enterCosts[child_idx][enum_vec[c]];
                }
                int highest_layer = highest_pin;
                int lowest_layer = lowest_pin;
                for (int e_layer : enum_vec) {
                    highest_layer = max(highest_layer, e_layer);
                    lowest_layer = min(lowest_layer, e_layer);
                }
                if(db::setting.patternRouteMemorization){
                    for (int layer_idx = 0; layer_idx < layer_cnt; layer_idx++) {
                        int via_bottom = min(lowest_layer, layer_idx);
                        int via_height = max(highest_layer, layer_idx) - via_bottom;                        

                        db::CostT via_cost = 
                            // viaStack_CostTb[fromNode.idx][via_height];
                            // grDatabase.getStackViaCost(gr::GrPoint(via_bottom, fromNode.x, fromNode.y), via_height);
                            costTb_hash[fromNode.idx][std::make_pair(via_bottom,via_height)];

                        
                        db::CostT cost = prev_cost + via_cost;
                        
                        if(debug)
                            log() << "node_idx: " << fromNode.idx
                                << ", b: " << via_bottom
                                << ", x: "  << fromNode.x
                                << ", y: " << fromNode.y 
                                << ", h: " << via_height
                                << ", pc: " << prev_cost
                                << ", cc: " << via_cost
                                << ", tc: " << cost << std::endl;
                        // if(debug){
                        //     auto tmp = grDatabase.getStackViaCost(gr::GrPoint(via_bottom, fromNode.x, fromNode.y), via_height,true);
                        //     log() << "grDatabase.getUnitViaMultiplier(): " << grDatabase.getUnitViaMultiplier() << std::endl;
                        //     log() << "grPoint: [l: " << via_bottom
                        //       << ", x: " << fromNode.x
                        //       << ", y: " << fromNode.y 
                        //       << "], via_height: " << via_height 
                        //       << " ,exitCost: " << via_cost << ", prev_cost: " << prev_cost << ", cost: " << cost << std::endl;
                        // }
                            
                        
                        if (cost < fromNode.exitCosts[layer_idx]) {
                            fromNode.exitCosts[layer_idx] = cost;
                            for (int c = 0; c < num_children; c++) {
                                fromNode.exitEnterLayers[fromNode.childIdxs[c]][layer_idx] = enum_vec[c];
                            }
                        }
                    }
                }else{
                    for (int layer_idx = 0; layer_idx < layer_cnt; layer_idx++) {
                        int via_bottom = min(lowest_layer, layer_idx);
                        int via_height = max(highest_layer, layer_idx) - via_bottom;
                        db::CostT via_cost =
                            grDatabase.getStackViaCostRelax(gr::GrPoint(via_bottom, fromNode.x, fromNode.y), via_height,wireMap_tpls,viaMap_tpls);
                        db::CostT cost = prev_cost + via_cost;
                        if (cost < fromNode.exitCosts[layer_idx]) {
                            fromNode.exitCosts[layer_idx] = cost;
                            for (int c = 0; c < num_children; c++) {
                                fromNode.exitEnterLayers[fromNode.childIdxs[c]][layer_idx] = enum_vec[c];
                            }
                        }
                    }// end loop
                }//end if memo
            }
        }
        // Patterns
        if(debug){
            log() << "refineplacement bf fromNode print..." << std::endl;
            fromNode.print();
            log() << "refineplacement bf toNode print..." << std::endl;
            toNode.print();
        }
        LShape(edge,wireMap_tpls,viaMap_tpls);
        if(debug){
            log() << "refineplacement af fromNode print..." << std::endl;
            fromNode.print();
            log() << "refineplacement af toNode print..." << std::endl;
            toNode.print();
        }
    }

    
}//end patternRouteMT


void InitRouteCell::patternRouteMT() {
    bool debug = false;
    // if(grNet.getName() == database.debug_net) debug = false;
    // add relaxation to router estimator 
    std::vector<std::pair<std::tuple<int,int,int>,double>> wireMap_tpls;
    std::vector<std::pair<std::tuple<int,int,int>,double>> viaMap_tpls;
    for (const auto& guide : grNet.wireRouteGuides) grDatabase.removeWireRelax(guide,wireMap_tpls);
    const auto& viaGuides = grNet.viaRouteGuides;
    for (int g1 = 0; g1 < viaGuides.size(); g1++) {
        for (int g2 = g1 + 1; g2 < viaGuides.size(); g2++) {
            if (abs(viaGuides[g1].layerIdx - viaGuides[g2].layerIdx) != 1) continue;

            auto xIntvl = viaGuides[g1][X].IntersectWith(viaGuides[g2][X]);
            auto yIntvl = viaGuides[g1][Y].IntersectWith(viaGuides[g2][Y]);

            if (xIntvl.IsValid() && yIntvl.IsValid()){
                grDatabase.removeViaRelax({min(viaGuides[g1].layerIdx, viaGuides[g2].layerIdx), xIntvl, yIntvl},viaMap_tpls);
            }
        }
    }




    if(debug)
        log() << "patternRouteRefinePlacement grNet: " << grNet.getName() << std::endl;
    if(debug){
        log() << "grNet: " << grNet.getName() << std::endl;
        log() << "database.getUnitViaCost(): " << database.getUnitViaCost() << std::endl;
        log() << "unitSqrtViaUsage: " << db::setting.unitSqrtViaUsage << std::endl;
        log() << "grDatabase.getLogisticSlope(): " << grDatabase.getLogisticSlope() << std::endl;
        log() << "grDatabase.getUnitViaMultiplier(): " << grDatabase.getUnitViaMultiplier() << std::endl;
    }


    int layer_cnt = database.getLayerNum();
    for (auto& edge : routeEdges) {
        auto& fromNode = routeNodes[edge.from];
        auto& toNode = routeNodes[edge.to];

        if(debug){
            log() << "fromNode idx: " << fromNode.idx
                  << "[x: " << fromNode.x
                  << ",y: " << fromNode.y << "]" << " -> "
                  << "toNode idx: " << toNode.idx
                  << "[x: " << toNode.x
                  << ",y: " << toNode.y << "]" << std::endl;
        }
        // Initialize fromNode's exitCosts
        fromNode.exitCosts.resize(layer_cnt, numeric_limits<db::CostT>::max());
        int lowest_pin = layer_cnt;
        int highest_pin = -1;
        for (int pin_layer : fromNode.pinLayers) {
            lowest_pin = min(lowest_pin, pin_layer);
            highest_pin = max(highest_pin, pin_layer);
        }
        if (fromNode.childIdxs.size() == 0) {  // no children nodes, only via costs
            for (int layer_idx = 0; layer_idx < layer_cnt; layer_idx++) {
                int lowest_layer = min(lowest_pin, layer_idx);
                int highest_layer = max(highest_pin, layer_idx);

                db::CostT exit_cost = grDatabase.getStackViaCostRelax(gr::GrPoint(lowest_layer, fromNode.x, fromNode.y),
                                                                 highest_layer - lowest_layer,wireMap_tpls,viaMap_tpls);
                fromNode.exitCosts[layer_idx] = exit_cost;
            }
        } else {
            // initialize exitEnterLayers map
            for (int child_idx : fromNode.childIdxs) {
                fromNode.exitEnterLayers[child_idx] = vector<int>(layer_cnt, -1);
            }
            // enumerate all children layer combinations
            int num_children = fromNode.childIdxs.size();
            #pragma omp parallel
            {
                for (int enum_idx = 0; enum_idx < pow(layer_cnt, num_children); enum_idx++) {
                    vector<int> enum_vec;
                    int e_idx = enum_idx;
                    for (int en = 0; en < num_children; en++) {
                        enum_vec.emplace_back(e_idx % layer_cnt);
                        e_idx /= layer_cnt;
                    }
                    db::CostT prev_cost = 0;
                    for (int c = 0; c < num_children; c++) {
                        int child_idx = fromNode.childIdxs[c];
                        prev_cost += fromNode.enterCosts[child_idx][enum_vec[c]];
                    }
                    int highest_layer = highest_pin;
                    int lowest_layer = lowest_pin;
                    for (int e_layer : enum_vec) {
                        highest_layer = max(highest_layer, e_layer);
                        lowest_layer = min(lowest_layer, e_layer);
                    }
                    #pragma omp for 
                    for (int layer_idx = 0; layer_idx < layer_cnt; layer_idx++) {
                        int via_bottom = min(lowest_layer, layer_idx);
                        int via_height = max(highest_layer, layer_idx) - via_bottom;
                        db::CostT via_cost = 
                            grDatabase.getStackViaCostRelax(gr::GrPoint(via_bottom, fromNode.x, fromNode.y), via_height,wireMap_tpls,viaMap_tpls);
                        db::CostT cost = prev_cost + via_cost;
                        if (cost < fromNode.exitCosts[layer_idx]) {
                            fromNode.exitCosts[layer_idx] = cost;
                            for (int c = 0; c < num_children; c++) {
                                fromNode.exitEnterLayers[fromNode.childIdxs[c]][layer_idx] = enum_vec[c];
                            }
                        }
                    }
                }//end pow(layer_cnt,num_childer
            }
            
        }
        // Patterns
        if(debug){
            log() << "refineplacement bf fromNode print..." << std::endl;
            fromNode.print();
            log() << "refineplacement bf toNode print..." << std::endl;
            toNode.print();
        }
        LShape(edge,wireMap_tpls,viaMap_tpls);
        if(debug){
            log() << "refineplacement af fromNode print..." << std::endl;
            fromNode.print();
            log() << "refineplacement af toNode print..." << std::endl;
            toNode.print();
        }
        
    }
}//end patternRouteMT

db::CostT InitRouteCell::getBufferedWireCost(gr::GrEdge edge) {
    if (wireCostBuffer.find(edge) != wireCostBuffer.end()) {
        return wireCostBuffer[edge];
    } else {
        auto wireCost = grDatabase.getWireCost(edge);
        wireCostBuffer[edge] = wireCost;
        return wireCost;
    }
}

db::CostT InitRouteCell::getBufferedWireCostRelax(gr::GrEdge edge, std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
                                    , std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls) {
    if (wireCostBuffer.find(edge) != wireCostBuffer.end()) {
        return wireCostBuffer[edge];
    } else {
        auto wireCost = grDatabase.getWireCostRelax(edge,wireMap_tpls,viaMap_tpls);
        wireCostBuffer[edge] = wireCost;
        return wireCost;
    }
}

void InitRouteCell::LShape(const RouteEdge& edge, std::vector<std::pair<std::tuple<int,int,int>,double>>& wireMap_tpls
                                    , std::vector<std::pair<std::tuple<int,int,int>,double>>& viaMap_tpls) {
    int layer_cnt =  database.getLayerNum();

    auto& fromNode = routeNodes[edge.from];
    auto& toNode = routeNodes[edge.to];

    vector<db::CostT> enter_cost(layer_cnt, numeric_limits<db::CostT>::max());
    vector<vector<gr::GrPoint>> enter_edges(layer_cnt);

    if (fromNode.x == toNode.x || fromNode.y == toNode.y) {  // has no bending point
        bool dir = (fromNode.x == toNode.x ? X : Y);
        for (int layer = 0; layer < layer_cnt; layer++) {
            if (database.getLayerDir(layer) != dir) continue;
            gr::GrEdge gr_edge(gr::GrPoint(layer, fromNode.x, fromNode.y), gr::GrPoint(layer, toNode.x, toNode.y));
            db::CostT edge_cost = grDatabase.getWireCostRelax(gr_edge,wireMap_tpls,viaMap_tpls);
            db::CostT cost = fromNode.exitCosts[layer] + edge_cost;
            if (cost < enter_cost[layer]) {
                enter_cost[layer] = cost;
                vector<gr::GrPoint> edge_vec;
                edge_vec.emplace_back(gr::GrPoint(layer, fromNode.x, fromNode.y));
                // edge_vec.emplace_back(gr::GrPoint(layer+1, fromNode.x, fromNode.y)); // additional
                // edge_vec.emplace_back(gr::GrPoint(layer, fromNode.x, fromNode.y)); // additional
                edge_vec.emplace_back(gr::GrPoint(layer, toNode.x, toNode.y));
                // edge_vec.emplace_back(gr::GrPoint(layer, toNode.x, toNode.y)); // additional
                // edge_vec.emplace_back(gr::GrPoint(layer, toNode.x, toNode.y)); // additional
                enter_edges[layer] = move(edge_vec);
            }
        }
    } else {  // has a bending point
        for (int to_layer = 0; to_layer < layer_cnt; to_layer++) {
            Dimension to_dir = database.getLayerDir(to_layer);
            int bend_x = -1;  // bending point
            int bend_y = -1;
            if (to_dir == Y) {
                bend_x = fromNode.x;
                bend_y = toNode.y;
            } else {
                bend_x = toNode.x;
                bend_y = fromNode.y;
            }

            // to_layer edge cost
            gr::GrEdge to_gr_edge(gr::GrPoint(to_layer, bend_x, bend_y), gr::GrPoint(to_layer, toNode.x, toNode.y));
            db::CostT to_edge_cost = grDatabase.getWireCostRelax(to_gr_edge,wireMap_tpls,viaMap_tpls);

            db::CostT min_cost = numeric_limits<db::CostT>::max();
            int best_from_layer = -1;
            for (int from_layer = 0; from_layer < layer_cnt; from_layer++) {
                if (database.getLayerDir(from_layer) == to_dir) continue;  // same routing direction

                // from_layer edge cost
                gr::GrEdge from_gr_edge(gr::GrPoint(from_layer, fromNode.x, fromNode.y),
                                        gr::GrPoint(from_layer, bend_x, bend_y));
                db::CostT from_edge_cost = getBufferedWireCostRelax(from_gr_edge,wireMap_tpls,viaMap_tpls);

                db::CostT bend_via_cost = grDatabase.getStackViaCostRelax(
                    gr::GrPoint(min(from_layer, to_layer), bend_x, bend_y), abs(to_layer - from_layer),wireMap_tpls,viaMap_tpls);

                db::CostT cost = fromNode.exitCosts[from_layer] + from_edge_cost + bend_via_cost + to_edge_cost;
                if (cost < min_cost) {
                    min_cost = cost;
                    best_from_layer = from_layer;
                }
            }
            enter_cost[to_layer] = min_cost;
            vector<gr::GrPoint> edge_vec;
            edge_vec.emplace_back(gr::GrPoint(best_from_layer, fromNode.x, fromNode.y));
            bool go_up = (best_from_layer < to_layer);
            for (int cur_layer = best_from_layer;;) {
                edge_vec.emplace_back(gr::GrPoint(cur_layer, bend_x, bend_y));
                if (cur_layer == to_layer) break;
                cur_layer += (go_up ? 1 : -1);
            }
            // edge_vec.emplace_back(gr::GrPoint(best_from_layer, bend_x, bend_y));
            // edge_vec.emplace_back(gr::GrPoint(to_layer, bend_x, bend_y));
            edge_vec.emplace_back(gr::GrPoint(to_layer, toNode.x, toNode.y));
            enter_edges[to_layer] = move(edge_vec);
        }
    }
    toNode.enterCosts[edge.from] = move(enter_cost);
    toNode.enterEdges[edge.from] = move(enter_edges);
}


bool InitRouteCell::buildTopo() {
    std::vector<std::pair<std::tuple<int,int,int>,double>> wireMap_tpls;
    std::vector<std::pair<std::tuple<int,int,int>,double>> viaMap_tpls;
    for (const auto& guide : grNet.wireRouteGuides) grDatabase.removeWireRelax(guide,wireMap_tpls);
    const auto& viaGuides = grNet.viaRouteGuides;
    for (int g1 = 0; g1 < viaGuides.size(); g1++) {
        for (int g2 = g1 + 1; g2 < viaGuides.size(); g2++) {
            if (abs(viaGuides[g1].layerIdx - viaGuides[g2].layerIdx) != 1) continue;

            auto xIntvl = viaGuides[g1][X].IntersectWith(viaGuides[g2][X]);
            auto yIntvl = viaGuides[g1][Y].IntersectWith(viaGuides[g2][Y]);

            if (xIntvl.IsValid() && yIntvl.IsValid()){
                grDatabase.removeViaRelax({min(viaGuides[g1].layerIdx, viaGuides[g2].layerIdx), xIntvl, yIntvl},viaMap_tpls);
            }
        }
    }



    int layer_cnt = database.getLayerNum();
    int rootIdx = -1;
    int root_enter_layer = -1;
    bool debug = false;
    // if(grNet.getName() == database.debug_net) debug = false;
    // build topo tree
    if (routeEdges.size() > 0) {
        rootIdx = routeEdges[routeEdges.size() - 1].to;
        if(debug){
            log() << "rootIdx: " << rootIdx 
                  << "routeEdge size: " <<routeEdges.size() << std::endl;
                 
        }
        RouteNode& rootNode = routeNodes[rootIdx];
        if(debug){
            log() << "rootNode idx: " << rootNode.idx
                  << "[x: " << rootNode.x
                  << ",y: " << rootNode.y << "]" << std::endl;
        }

        int highest_layer = -1;
        int lowest_layer = layer_cnt;
        for (int pin_layer : rootNode.pinLayers) {
            lowest_layer = min(lowest_layer, pin_layer);
            highest_layer = max(highest_layer, pin_layer);
        }

        if(debug)
            log() << "lowest_layer: " << lowest_layer 
                  <<  "highest_layer: " << highest_layer << std::endl;
        // find min cost enter layer
        assert(rootNode.enterCosts.size() == 1);
        auto& enterCosts = rootNode.enterCosts.begin()->second;
        db::CostT min_cost = numeric_limits<db::CostT>::max();
        for (int layer_idx = 0; layer_idx < layer_cnt; layer_idx++) {
            db::CostT cost = enterCosts[layer_idx];

            if(debug)
                log() << "layer_idx: " << layer_idx << ", enterCosts: " << cost << std::endl;

            int via_bottom = min(lowest_layer, layer_idx);
            int via_height = max(highest_layer, layer_idx) - via_bottom;
            db::CostT via_cost =
                grDatabase.getStackViaCostRelax(gr::GrPoint(via_bottom, rootNode.x, rootNode.y), via_height,wireMap_tpls,viaMap_tpls);
            if(debug)
                log() << "via_cost: " << via_cost << ", height: " << via_height << std::endl;
        
            cost += via_cost;  // via cost
            if (cost < min_cost) {
                min_cost = cost;
                root_enter_layer = layer_idx;
                if(debug)
                    log() << "min_cost: " << min_cost <<", root_enter_layer: "
                          << root_enter_layer << std::endl;
            }
        }
    } else {  // local net
        rootIdx = 0;
        root_enter_layer = routeNodes[rootIdx].pinLayers[0];
    }

    // log() << "cost " << min_cost << endl;

    function<shared_ptr<gr::GrSteiner>(int, int)> buildTopo = [&](int root_idx, int exit_layer) {
        if(debug)
            log() << "inside buildTopo..." << std::endl;

        RouteNode& rNode = routeNodes[root_idx];

        if(debug)
            log() << "rNode idx: " << rNode.idx
                  << "[x: " << rNode.x
                  << ",y: " << rNode.y << "]" << std::endl;

        shared_ptr<gr::GrSteiner> topo_root = make_shared<gr::GrSteiner>(gr::GrPoint(exit_layer, rNode.x, rNode.y), -1);

        // build layer topo
        int low_layer = exit_layer;
        int high_layer = exit_layer;

        if(debug)
            log() << "exit layer: " << exit_layer << std::endl;

        for (int child_idx : rNode.childIdxs) {
            if(debug)
                log() << "child_idx: " << child_idx << std::endl;

            int child_enter_layer = exit_layer;  // root has no exit, thus no exitEnterLayers map, we use exit layer to
                                                 // indicate its child's enter layer
            if (rNode.exitEnterLayers.size() > 0) {
                child_enter_layer = rNode.exitEnterLayers[child_idx][exit_layer];
            }
            low_layer = min(low_layer, child_enter_layer);
            high_layer = max(high_layer, child_enter_layer);
        }

        for (int pin_layer : rNode.pinLayers) {
            if (pin_layer != exit_layer) {
                low_layer = min(low_layer, pin_layer);
                high_layer = max(high_layer, pin_layer);
            }
        }

        unordered_map<int, shared_ptr<gr::GrSteiner>> layer2layerTopo;
        layer2layerTopo[exit_layer] = topo_root;
        for (int cur_layer = exit_layer + 1; cur_layer <= high_layer; cur_layer++) {
            layer2layerTopo[cur_layer] = make_shared<gr::GrSteiner>(gr::GrPoint(cur_layer, rNode.x, rNode.y), -1);
            gr::GrSteiner::setParent(layer2layerTopo[cur_layer], layer2layerTopo[cur_layer - 1]);
        }
        for (int cur_layer = exit_layer - 1; cur_layer >= low_layer; cur_layer--) {
            layer2layerTopo[cur_layer] = make_shared<gr::GrSteiner>(gr::GrPoint(cur_layer, rNode.x, rNode.y), -1);
            gr::GrSteiner::setParent(layer2layerTopo[cur_layer], layer2layerTopo[cur_layer + 1]);
        }

        for (int child_idx : rNode.childIdxs) {  // build branches
            int child_enter_layer = exit_layer;
            if (rNode.exitEnterLayers.size() > 0) {
                child_enter_layer = rNode.exitEnterLayers[child_idx][exit_layer];
            }
            auto& enter_edge = rNode.enterEdges[child_idx][child_enter_layer];
            if (enter_edge.size() == 0) {
                log() << "Error: Enter edge missing" << endl;
            }

            int child_exit_layer = enter_edge[0].layerIdx;
            shared_ptr<gr::GrSteiner> child_topo = buildTopo(child_idx, child_exit_layer);
            for (auto& point : enter_edge) {
                shared_ptr<gr::GrSteiner> steiner = make_shared<gr::GrSteiner>(point, -1);
                gr::GrSteiner::setParent(child_topo, steiner);
                child_topo = steiner;
            }
            // if (gr::GrPoint(*child_topo) !=  gr::GrPoint(*topo_root)) {
            gr::GrSteiner::setParent(child_topo, layer2layerTopo[child_enter_layer]);
            // } else {
            //     gr::GrSteiner::setParent(child_topo->children[0], topo_root);
            // }
        }

        return topo_root;
    };

    // dumpTo grNet
    vector<std::shared_ptr<gr::GrSteiner>> gridTopo;
    gridTopo.emplace_back(buildTopo(rootIdx, root_enter_layer));
    gr::GrSteiner::removeRedundancy(gridTopo[0]);

    vector<gr::GrPoint> pin_locs;
    for (auto& kv : routeNodes) {
        auto& node = kv.second;
        for (int pin_layer : node.pinLayers) {
            // log() << "(" << pin_layer << "," << node.x << "," << node.y << ") " << std::endl;
            pin_locs.emplace_back(gr::GrPoint(pin_layer, node.x, node.y));
        }
    }

    // bool debug = false;
    // if(grNet.getName() == database.debug_net) debug = false;

    vector<gr::GrBoxOnLayer> wireRouteGuides;
    vector<gr::GrBoxOnLayer> viaRouteGuides;

    wl_cost= 0;
    cong_cost= 0;
    via_cost= 0;
    via_abs_cost= 0;
    wl_abs_cost= 0;
    path_cost= 0;
    // Note: generate guides by topology
    auto genTopoGuide = [&](std::shared_ptr<gr::GrSteiner> node) {
        bool debug = false;
        if(debug)
            log() << "genTopoGuide ... " << std::endl;
        auto parent = node;
        for (auto child : parent->children) {
            if (parent->layerIdx == child->layerIdx) {
                std::shared_ptr<gr::GrSteiner> lower, upper;
                if ((*parent)[X] < (*child)[X] || (*parent)[Y] < (*child)[Y]) {
                    lower = parent;
                    upper = child;
                } else {
                    lower = child;
                    upper = parent;
                }

                gr::GrEdge edge(gr::GrPoint(lower->layerIdx, lower->x , lower->y ),
                                gr::GrPoint(upper->layerIdx, upper->x, upper->y) );

                DBU edge_cost = grDatabase.getWireCostRelax(edge,wireMap_tpls,viaMap_tpls);

                // give an coefficent to edge cost that are wire. (1->10->100)
                path_cost += edge_cost;
                wl_cost += edge_cost;
                wl_abs_cost += grDatabase.getDist(*parent, *child);

                wireRouteGuides.emplace_back(lower->layerIdx,
                                            utils::IntervalT<int>((*lower)[X], (*upper)[X]),
                                            utils::IntervalT<int>((*lower)[Y], (*upper)[Y]));
                if(debug)
                    log() << "lower: [x: " << lower->x << ", y: " << lower->y 
                        << ", l: " << lower->layerIdx
                        << "], upper: [x: " << upper->x << ", y: " << upper->y 
                        << "l: " << upper->layerIdx << "]"
                        << ", edge_cost: " << edge_cost
                        << std::endl;
            } else {
               
                auto lower_layer = std::min(parent->layerIdx,child->layerIdx);
                auto upper_layer = std::max(parent->layerIdx,child->layerIdx);

                viaRouteGuides.emplace_back(parent->layerIdx,
                                                utils::IntervalT<int>((*parent)[X], (*parent)[X]),
                                                utils::IntervalT<int>((*parent)[Y], (*parent)[Y]));
                viaRouteGuides.emplace_back(child->layerIdx,
                                                utils::IntervalT<int>((*child)[X], (*child)[X]),
                                                utils::IntervalT<int>((*child)[Y], (*child)[Y]));

                DBU via_pt = grDatabase.getStackViaCostRelax(gr::GrPoint(lower_layer, parent->x , parent->y),
                                                std::abs(upper_layer-lower_layer),wireMap_tpls,viaMap_tpls);
                path_cost +=via_pt;
                    
                via_cost +=via_pt;

                via_abs_cost +=1;

                if(debug)
                    log() << "parent: [x: " << parent->x << ", y: " << parent->y 
                        << ", l: " << parent->layerIdx
                        << ", child: [x: " << child->x << ", y: " << child->y 
                        << "l: " << child->layerIdx << "]" 
                        << ", via_pt: " << via_pt
                        << std::endl;
            }
        }
    };


    for (const std::shared_ptr<gr::GrSteiner>& tree : gridTopo) {
        gr::GrSteiner::postOrder(tree,genTopoGuide );
    }

    //-------------------------------
    sliceGuides(wireRouteGuides);
    sliceGuides(viaRouteGuides);

    

    DBU wirelength = 0;
    auto wl_fun = [&](std::shared_ptr<gr::GrSteiner> node) {
        auto parent = node;
        for (auto child : parent->children) {
            if (parent->layerIdx == child->layerIdx) wirelength += grDatabase.getDist(*parent, *child);
        }
    };
    for (const std::shared_ptr<gr::GrSteiner>& tree : gridTopo) {
        gr::GrSteiner::postOrder(tree,wl_fun );
    }
    // log() << "wirelength cell: " << wirelength << std::endl;
    // log() << "path_cost: " << path_cost << std::endl;
    // wl_cost = path_cost;//wirelength;
    // via_cost = viaRouteGuides.size();


    // debug
    // if(grNet.getName() == "net121")
    //     debug=true;
    if(debug)
    {
        log() << "post order " << std::endl;
        log() << endl;
        auto visit = [&](std::shared_ptr<gr::GrSteiner> topo) {
            if (topo->children.size() == 0) return;
            log() << gr::GrPoint(*topo) << "->";
            for(auto& child: topo->children) {
                cout << gr::GrPoint(*child) << " ";
            }
            cout << endl;
        };

        for (const std::shared_ptr<gr::GrSteiner>& tree : gridTopo) {
            gr::GrSteiner::postOrder(tree,visit );
        }
        // grNet.preOrderVisitGridTopo([&](std::shared_ptr<gr::GrSteiner> topo) {
        //     if (topo->children.size() == 0) return;
        //     log() << gr::GrPoint(*topo) << "->";
        //     for(auto& child: topo->children) {
        //         cout << gr::GrPoint(*child) << " ";
        //     }
        //     cout << endl;
        // });
        // log() << "cost " << min_cost << endl;
        log() << endl << "+ + + + + + + + + + + + + + + + + + + + + +" << endl;


        // Note: generate guides by topology
        auto genTopoGuide = [&](std::shared_ptr<gr::GrSteiner> node) {
            auto parent = node;
            for (auto child : parent->children) {
                if (parent->layerIdx == child->layerIdx) {
                    std::shared_ptr<gr::GrSteiner> lower, upper;
                    if ((*parent)[X] < (*child)[X] || (*parent)[Y] < (*child)[Y]) {
                        lower = parent;
                        upper = child;
                    } else {
                        lower = child;
                        upper = parent;
                    }
                    wireRouteGuides.emplace_back(lower->layerIdx,
                                                utils::IntervalT<int>((*lower)[X], (*upper)[X]),
                                                utils::IntervalT<int>((*lower)[Y], (*upper)[Y]));
                } else {
                    viaRouteGuides.emplace_back(parent->layerIdx,
                                                    utils::IntervalT<int>((*parent)[X], (*parent)[X]),
                                                    utils::IntervalT<int>((*parent)[Y], (*parent)[Y]));
                    viaRouteGuides.emplace_back(child->layerIdx,
                                                    utils::IntervalT<int>((*child)[X], (*child)[X]),
                                                    utils::IntervalT<int>((*child)[Y], (*child)[Y]));
                }
            }
        };


        for (const std::shared_ptr<gr::GrSteiner>& tree : gridTopo) {
            gr::GrSteiner::postOrder(tree,genTopoGuide );
        }


        sliceGuides(wireRouteGuides);
        sliceGuides(viaRouteGuides);

        if(debug){
            log() << "wireRouteGuides: " << wireRouteGuides.size() << std::endl;
            log() << "viaRouteGuides: " << viaRouteGuides.size() << std::endl;
            for (const auto& guide : wireRouteGuides) {
                log() << guide.lx() 
                << ", " << guide.ly() 
                << ", " << guide.hx()+1
                << ", " << guide.hy()+1
                << ", " << guide.layerIdx << std::endl; 
            }
            
            for (const auto& guide : viaRouteGuides) {
                log() << guide.lx() 
                << ", " << guide.ly() 
                << ", " << guide.hx()+1
                << ", " << guide.hy()+1
                << ", " << guide.layerIdx << std::endl; 
            }


            log() << "after slice wireRouteGuides: " << wireRouteGuides.size() << std::endl;
            log() << "after slice viaRouteGuides: " << viaRouteGuides.size() << std::endl;
            
        }

        

        DBU wirelength = 0;
        auto wl_fun = [&](std::shared_ptr<gr::GrSteiner> node) {
            auto parent = node;
            for (auto child : parent->children) {
                if (parent->layerIdx == child->layerIdx) wirelength += grDatabase.getDist(*parent, *child);
            }
        };
        for (const std::shared_ptr<gr::GrSteiner>& tree : gridTopo) {
            gr::GrSteiner::postOrder(tree,wl_fun );
        }
        log() << "wirelength cell: " << wirelength << std::endl;

       
    }
   
   
    // if(grNet.getName() == "net121")
    //     log() << "net: " << grNet.getName() 
    //         << ", path cost " << path_cost << std::endl;
    

    int numvio = grDatabase.getVioReportRelax(grNet
                           , wireMap_tpls
                           , viaMap_tpls 
                           , wireRouteGuides
                           , viaRouteGuides);
    // if(numvio > 0)
    //     log() << "grNet: " << grNet.getName() << " has violation! " << std::endl;
          


    // Check Connectivity and Mark Pins (0: no pin; 1: has pin)
    if (!gr::GrSteiner::checkConnectivity(gridTopo[0], pin_locs)) {
        log() << "Error: Connectivity check failed" << endl;
    }

    return numvio>0;
}

void InitRouteCell::sliceGuides(vector<gr::GrBoxOnLayer> &guides, bool mergeAdj) {
    vector<vector<gr::GrBoxOnLayer>> tmpGuides(database.getLayerNum());  // route guides on different layers
    for (auto &guide : guides) {
        tmpGuides[guide.layerIdx].push_back(guide);
    }
    guides.clear();
    for (int layerIdx = 0; layerIdx < database.getLayerNum(); ++layerIdx) {
        gr::GrBoxOnLayer::sliceGrPolygons(tmpGuides[layerIdx], mergeAdj);
        for (auto &guide : tmpGuides[layerIdx]) guides.push_back(guide);
    }
}

void InitRouteCell::plan_fluteOnly() {
    // runFlute();
    runFluteV2();
    
    // bool debug = false;
    // if(grNet.getName() == database.debug_net) debug = false;


    // if(debug){
    //     log() << "after plan_fluteOnly (initRouteCell)..." << std::endl;
    //     for (int i = 0; i < routeEdges.size(); i++) {
    //         log() << "aaaaaaaaaaaaaaaaaaa" << std::endl;
    //         auto& edge = routeEdges[i];
    //         auto& fromNode = routeNodes[edge.from];
    //         auto& toNode = routeNodes[edge.to];
    //         log() << "fromNode idx: " << fromNode.idx
    //             << "[x: " << fromNode.x
    //             << ",y: " << fromNode.y << "]" << " -> "
    //             << "toNode idx: " << toNode.idx
    //             << "[x: " << toNode.x
    //             << ",y: " << toNode.y << "]" << std::endl;
    //     }
        
    // }

    // int s = routeNodes.begin()->first;
    // queue<int> tmp_q;
    // tmp_q.push(s);
    // set<int> vis;
    // while (tmp_q.empty() == false) {
    //     int u = tmp_q.front();
    //     tmp_q.pop();
    //     vis.insert(u);
    //     RouteNode& node = routeNodes[u];
    //     for (int nIdx : node.toConnect) {
    //         if (vis.find(nIdx) != vis.end()) continue;
    //         addUsage2D(node, routeNodes[nIdx], 1);  // update the usage
    //         tmp_q.push(nIdx);
    //     }
    // }
    // costEstimator(routeNodes);
}
// void InitRouteCell::costEstimator(std::map<int, RouteNode>& cur_routeNodes) {
//     // log() << "endge_shift2d ..." << std::endl;
//     // log() << "cur_routeNodes size: " << cur_routeNodes.size() << std::endl;
//     if (cur_routeNodes.size() == 1) return;
//     vector<RouteEdge> edgeset;
//     map<pair<int, int>, set<int>> loc_nodes;  // will be used to merge overlapping nodes
//     int s = cur_routeNodes.begin()->first;
//     queue<int> tmp_q;
//     tmp_q.push(s);
//     set<int> vis;
//     // log() << "grNet: " << grNet.getName()  << std::endl;

//     while (tmp_q.empty() == false) {
//         int u = tmp_q.front();
//         tmp_q.pop();
//         vis.insert(u);

        
//         RouteNode& node = cur_routeNodes[u];
//         loc_nodes[make_pair(node.x, node.y)].insert(u);

//         // log() << "u: " << u << ", node.x: " << node.x << ", node.y: " << node.y << std::endl;

//         for (int nIdx : node.toConnect) {
//             if (vis.find(nIdx) != vis.end()) continue;
//             RouteEdge edge;
//             edge.from = nIdx;
//             edge.to = u;
//             edgeset.emplace_back(edge);
//             // removeUsage2D(node, cur_routeNodes[nIdx], 1);  // Remove this net
//             tmp_q.push(nIdx);
//         }
//     }

//     // log() << "init edge set: " << std::endl;
//     // for(auto e_tmp : edgeset){
//     //     log() << "edge from: " << e_tmp.from << " -> " << e_tmp.to << std::endl;
//     // }


//     auto getStraightCost = [&](int dir, int gridline, utils::IntervalT<int> range) {
//         double cost = 0;
//         for (int cp = range.low; cp < range.high; cp++) {
//             cost += grDatabase.getCost2D(dir, gridline, cp);
//         }
//         return cost;
//     };

//     std::function<double(int, int, int, int)> getCost =
//         [&](int x1, int y1, int x2, int y2) {  // GR point (x1,y1) -> (x2,y2)
//             double ret = 0;
//             int layer_dir_count = 0;
//             if (x1 == x2 || y1 == y2) {  // straight line
//                 int dir = (x1 == x2 ? X : Y);
//                 if (dir == X) {
//                     utils::IntervalT<int> Range(min(y1, y2), max(y1, y2));
//                     ret = getStraightCost(X, x1, Range);
//                 } else {
//                     utils::IntervalT<int> Range(min(x1, x2), max(x1, x2));
//                     ret = getStraightCost(Y, y1, Range);
//                 }
//             } else {  // choose the L with lower cost
//                 int x3 = x1;
//                 int y3 = y2;
//                 int x4 = x2;
//                 int y4 = y1;
//                 ret = min(getCost(x1, y1, x3, y3) + getCost(x2, y2, x3, y3),
//                           getCost(x1, y1, x4, y4) + getCost(x2, y2, x4, y4));
//             }
//             return ret;
//         };


//        std::function<double(int, int, int, int)> getCostVia =
//         [&](int x1, int y1, int x2, int y2) {  // GR point (x1,y1) -> (x2,y2)
//             // double numVias = 0;
//             // double lower_bound = 0;
//             // double upper_bound = 0;
//             // int layer_dir_count = 0;
//             if (x1 == x2 || y1 == y2) {  // straight line
//                 // lower_bound = 2;
//                 // upper_bound = database.getLayerNum() * 2;
//                 return 0;
//             } 
//             // else {  // choose the L with lower cost
//                 // int x3 = x1;
//                 // int y3 = y2;
//                 // int x4 = x2;
//                 // int y4 = y1;
//                 // ret = min(getCost(x1, y1, x3, y3) + getCost(x2, y2, x3, y3),
//                 //           getCost(x1, y1, x4, y4) + getCost(x2, y2, x4, y4));
//                 // lower_bound = 4;
//                 // // for upper bound kind of having multuplier to congestion 
//                 // upper_bound = database.getLayerNum() * 4;
//                 // return 1;
//             // }
//             // the edge is lshape and route nodes are not aligned
//             return 1;
//             // return (lower_bound+upper_bound)/2;
//         };
    

//     auto estimateCostEdgeCong = [&](double& edges_cost) {        
//         for (int i = 0; i < edgeset.size(); i++) {
//             auto& edge = edgeset[i];
//             auto& fromNode = cur_routeNodes[edge.from];
//             auto& toNode = cur_routeNodes[edge.to];
//             edges_cost += getCost(fromNode.x,
//                                 fromNode.y,
//                                 toNode.x,
//                                 toNode.y);
//         }
//     };

//     auto estimateCostEdgeVia = [&](double& edges_cost) {      
//         // log() << "grNet: " << grNet.getName() << std::endl;
//         // log() << "estimate Via cost " << std::endl;

//         // for(auto node : routeNodes){
//         //     log() << "node: [x: " << x
//         // }

//         // if (edgeset.size() == 0) edges_cost = 1;
//         for (int i = 0; i < edgeset.size(); i++) {
//             auto& edge = edgeset[i];
//             auto& fromNode = cur_routeNodes[edge.from];
//             auto& toNode = cur_routeNodes[edge.to];
//             // log() << "fromNode: [x: " << fromNode.x << ", y: " << fromNode.y << "]"
//             //       << ", toNode: [x: " << toNode.x << ", y: " << toNode.y << "]" << std::endl;
//             edges_cost += getCostVia(fromNode.x,
//                                 fromNode.y,
//                                 toNode.x,
//                                 toNode.y);
//         }
//     };


    
//     estimateCostEdgeCong(cong_cost);
//     estimateCostEdgeVia(via_cost);
//     // log() << "cong_cost: " << cong_cost << std::endl;
//     // log() << "via_cost: " << via_cost << std::endl;
    
// }

void InitRouteCell::getRoutingOrder() {

    // log() << "getRoutingOrder net: " << grNet.getName() << std::endl;
    // local net in one gcell
    if (routeNodes.size() == 1) return;
    // pick a degree 1 node as root
    int root = -1;
    for (auto it = routeNodes.begin(); it != routeNodes.end(); it++) {
        if (it->second.degree() == 1) {
            root = it->first;
            break;
        }
    }

    // dfs to decide routing order
    set<int> visited;
    function<void(int)> dfs = [&](int nodeIdx) {
        visited.insert(nodeIdx);
        RouteNode& node = routeNodes[nodeIdx];
        for (int nIdx : node.toConnect) {
            if (visited.find(nIdx) != visited.end()) continue;
            RouteEdge edge;
            edge.from = nIdx;
            edge.to = nodeIdx;
            routeEdges.emplace_back(edge);
            routeNodes[edge.to].childIdxs.emplace_back(edge.from);
            dfs(nIdx);
        }
    };

    if (root != -1) {
        dfs(root);
        reverse(routeEdges.begin(), routeEdges.end());
    } else {
        log() << "Error: Can't find a degree one node (2)" << endl;
    }

    // bool debug = false;
    // if(grNet.getName() == database.debug_net) debug = false;
    if(debug){
        log() << "getRoutingOrder routeEdge: " << std::endl;
        for(auto e_tmp : routeEdges){
            auto& fromNode_tmp = routeNodes[e_tmp.from];
            auto& toNode_tmp = routeNodes[e_tmp.to];
            log() << "edge from: " << e_tmp.from 
                  << "[x: " <<fromNode_tmp.x << ", " << "y: " << fromNode_tmp.y  << "]"
                  << " -> " << e_tmp.to 
                  << "[x: " <<toNode_tmp.x << ", " << "y: " << toNode_tmp.y  << "]"
                  << std::endl;
        }
    }
    
}


void InitRouteCell::getGrPinAccessBoxes(const gr::GCellGrid& gcellGrid) {
    bool debug = false;
    // if(grNet.getName() == "net121") debug = true;
    if(debug) {
        log() << "updatePinAccessBoxes...(initRouteCell grNet)" << std::endl;
    }
    // transform coor to grPoint, construct pinAccessBoxes
    grPinAccessBoxes.clear();
    grPinAccessBoxes.resize(grNet.numOfPins());
    for (int i = 0; i < grNet.numOfPins(); i++) {
        const auto& boxes = pinAccessBoxes[i];
        std::unordered_set<gr::GrPoint> pointSet;
    
        

         DBU smallestVio = std::numeric_limits<DBU>::max();
         vector<const db::BoxOnLayer*> smallestBoxes;
    
        for (const auto& box : boxes) {
            if(debug)
                log() << "i: " << i    
                    << ", box: " << box << std::endl;
            int vio = database.getOvlpFixedMetalArea(box, grNet.dbNet.idx);
            if (vio <= smallestVio) {
                if (vio < smallestVio) smallestBoxes.clear();
                smallestVio = vio;
                smallestBoxes.push_back(&box);
            }
            auto grBox = gcellGrid.rangeSearchGCell(box);//3005436
            if(debug){
                log() << "i: " << i    
                        << ", box: " << box << std::endl;
                log() << "grBox: " << i    
                        << ", lx: " << grBox.lx()
                        << ", ly: " << grBox.ly()
                        << ", hx: " << grBox.hx()
                        << ", hy: " << grBox.hy() << std::endl;

                log() << "vio: " << vio << std::endl;
            }
                

            if (vio == 0) {
                
                
                for (int x = grBox[X].low; x <= grBox[X].high; x++)
                    for (int y = grBox[Y].low; y <= grBox[Y].high; y++) pointSet.emplace(box.layerIdx, x, y);
            }
        }
    //     // all have vio, add those with smallest vio
        if (pointSet.empty()) {
            for (auto box : smallestBoxes) {
                auto grBox = gcellGrid.rangeSearchGCell(*box);
                for (int x = grBox[X].low; x <= grBox[X].high; x++)
                    for (int y = grBox[Y].low; y <= grBox[Y].high; y++) pointSet.emplace(box->layerIdx, x, y);
            }
        }

        for (auto& point : pointSet) {
            grPinAccessBoxes[i].push_back(point);
            if(debug){
                log() << "pt.x: " << point.x
                      << ", pt.y: " << point.y << std::endl;
                log() << "intrval: " << gcellGrid.getCoorIntvl(point,0) << std::endl;
            }
                
        }
    }
}//update pinAccesboxes

// void InitRoute::addUsage2D(RouteNode& u, RouteNode& v, double usage) {
//     if (u.x == v.x) {  // straight edge dir = X
//         int gridline = u.x;
//         utils::IntervalT<int> Range(min(u.y, v.y), max(u.y, v.y));
//         for (int cp = Range.low; cp < Range.high; cp++) {
//             grDatabase.useWire2D(X, gridline, cp, usage);
//         }
//     } else if (u.y == v.y) {
//         int gridline = u.y;
//         utils::IntervalT<int> Range(min(u.x, v.x), max(u.x, v.x));
//         for (int cp = Range.low; cp < Range.high; cp++) {
//             grDatabase.useWire2D(Y, gridline, cp, usage);
//         }
//     } else {  // diagnal edge
//         RouteNode new1(u.x, v.y), new2(v.x, u.y);
//         addUsage2D(new1, u, usage / 2);
//         addUsage2D(new1, v, usage / 2);
//         addUsage2D(new2, u, usage / 2);
//         addUsage2D(new2, v, usage / 2);
//     }
// }

// void InitRoute::removeUsage2D(RouteNode& u, RouteNode& v, double usage) {
//     if (u.x == v.x) {  // straight edge dir = Y
//         int gridline = u.x;
//         utils::IntervalT<int> Range(min(u.y, v.y), max(u.y, v.y));
//         for (int cp = Range.low; cp < Range.high; cp++) {
//             grDatabase.removeUsage2D(X, gridline, cp, usage);
//         }
//     } else if (u.y == v.y) {
//         int gridline = u.y;
//         utils::IntervalT<int> Range(min(u.x, v.x), max(u.x, v.x));
//         for (int cp = Range.low; cp < Range.high; cp++) {
//             grDatabase.removeUsage2D(Y, gridline, cp, usage);
//         }
//     } else {  // diagnal edge
//         RouteNode new1(u.x, v.y), new2(v.x, u.y);
//         removeUsage2D(new1, u, usage / 2);
//         removeUsage2D(new1, v, usage / 2);
//         removeUsage2D(new2, u, usage / 2);
//         removeUsage2D(new2, v, usage / 2);
//     }
// }

// std::map<int, RouteNode>& InitRoute::getRouteNodes() { return routeNodes; }

};