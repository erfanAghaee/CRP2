#include "GenGuide.h"

GuideGeneratorStat guideGenStat;

void GuideGenerator::sliceGuides(vector<gr::GrBoxOnLayer> &guides, bool mergeAdj) {
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

void GuideGenerator::genConnGuides() {
    // keep connectivity
    for (auto &point : grNet.ovlpPoints) {
        grNet.viaRouteGuides.emplace_back(
            point.layerIdx, utils::IntervalT<int>(point[X], point[X]), utils::IntervalT<int>(point[Y], point[Y]));
        if (point.layerIdx + 1 < database.getLayerNum())
            grNet.viaRouteGuides.emplace_back(point.layerIdx + 1,
                                              utils::IntervalT<int>(point[X], point[X]),
                                              utils::IntervalT<int>(point[Y], point[Y]));
        if (point.layerIdx - 1 >= 0)
            grNet.viaRouteGuides.emplace_back(point.layerIdx - 1,
                                              utils::IntervalT<int>(point[X], point[X]),
                                              utils::IntervalT<int>(point[Y], point[Y]));
    }
}

void GuideGenerator::patchPinRegions() {

    bool debug = false;
    if(grNet.getName() == "n_5371"){
        debug = true;
    }

    double patchThresh = 2.0;
    // bool patchPinRegionsMode=true;// 0 is default mode and 1 is advanced mode
    // check if the region allows patching
    auto needPatch = [&](gr::GrPoint point) { return grDatabase.getCellResource(point) < patchThresh; };
    // get surrounding gcells of a point (normally 3 x 3, but will be different at boudary)
    // auto getSurrounding = [&](gr::GrPoint point) {
    //     gr::GrBoxOnLayer box;
    //     box.layerIdx = point.layerIdx;
    //     box.x.low = max(point.x - 1, 0);
    //     box.y.low = max(point.y - 1, 0);
    //     box.x.high = min(point.x + 1, grDatabase.getNumGrPoint(X) - 1);
    //     box.y.high = min(point.y + 1, grDatabase.getNumGrPoint(Y) - 1);

    //     return box;
    // };
    
    auto getSurrounding = [&](gr::GrPoint point) {
        gr::GrBoxOnLayer box;
        
        // bool isCrossBlock=false;
        box.layerIdx = point.layerIdx;
        box.x.low = max(point.x -1, 0);
        box.y.low = max(point.y -1, 0);
        box.x.high = min(point.x+1, grDatabase.getNumGrPoint(X) - 1);
        box.y.high = min(point.y+1, grDatabase.getNumGrPoint(Y) - 1);


        // log() << "box: " << box << std::endl;
        // for(int i = box.x.low; i < box.x.high; i++){
        //     for(int j = box.y.low; j < box.y.high; j++){
        //         log() << "i: " << i << ", j: " << j << std::endl;
        //         auto xl = grDatabase.getCoor(i , X);
        //         auto yl = grDatabase.getCoor(j , Y);
        //         auto xh = grDatabase.getCoor(i+1 , X);
        //         auto yh = grDatabase.getCoor(j+1 , Y);
        //         gr::GrBoxOnLayer box_noblockage(xl,yl,xh,yh);
        //         auto results = database.queryBlockageRTree(box_noblockage,eps);
        //         if(results.size() > 0){
        //             isCrossBlock=true;
        //             break;
        //         }


        //     }//end loop j
        // }//end loop i

        // if(isCrossBlock){
        //     box.x.low = max(point.x - 1 , 0);
        //     box.y.low = max(point.y - 1, 0);
        //     box.x.high = min(point.x + 1, grDatabase.getNumGrPoint(X) - 1);
        //     box.y.high = min(point.y + 1, grDatabase.getNumGrPoint(Y) - 1);
        // }else{
        //     box.x.low = max(point.x - 1 , 0);
        //     box.y.low = max(point.y - 1, 0);
        //     box.x.high = min(point.x + 1, grDatabase.getNumGrPoint(X) - 1);
        //     box.y.high = min(point.y + 1, grDatabase.getNumGrPoint(Y) - 1);
        // }

        

        // auto xl = grDatabase.getCoor(box.x.low , X);
        // auto yl = grDatabase.getCoor(box.y.low , Y);
        // auto xh = grDatabase.getCoor(box.x.high , X);
        // auto yh = grDatabase.getCoor(box.y.high , Y);
        // auto xl = grDatabase.getCoor(guide[Y].low, Y);
        // auto xl = grDatabase.getCoor(guide[X].high + 1, X);
        // auto xl = grDatabase.getCoor(guide[Y].high + 1, Y);

        return box;
    };

    // check if pin access point is in upper layers 
    // avoid generating patch. This is designed to solve
    // not meet guide in detailed routing is ispd19_test5
    auto checkValidPatching = [&](gr::GrNet& grNet) {
        // bool isValid = true;
        for (auto &pbxs : grNet.pinAccessBoxes){
            for (auto &pbx : pbxs) {
                if(pbx.layerIdx >= 1){
                    return false;
                }
            }//end loop 
        }//end loop
        return true;
    };

    if(!db::setting.patchPinRegionsMode){
        bool isValidPatch = true;
        isValidPatch = checkValidPatching(grNet);
        if(isValidPatch){
            for (auto &pbxs : grNet.pinAccessBoxes) {
                for (auto &pbx : pbxs) {
                    bool patched = false;
                    // // check if pbx is covered inside blockage
                    // auto boxTmpGR = gr::GrBoxOnLayer(pbx.layerIdx, {pbx.x}, {pbx.y});
                    // auto xl = grDatabase.getCoor(boxTmpGR[X].low, X);
                    // auto yl = grDatabase.getCoor(boxTmpGR[Y].low, Y);
                    // auto xh = grDatabase.getCoor(boxTmpGR[X].high+1, X);
                    // auto yh = grDatabase.getCoor(boxTmpGR[Y].high+1, Y);
                    // utils::BoxT<DBU> boxTmp(xl,yl,xh,yh); 
                    // if(debug){
                    //     log() << "pbx: " << pbx << std::endl;
                    //     log() << "accesspt: " << boxTmp << std::endl;
                    // }
                    // DBU eps = 10;
                    // auto results = database.queryBlockageRTree(boxTmp,eps);
                    // bool isInsideBlockage = false;
                    // for(auto result : results){
                    //     if(debug){
                    //         log() << "res: " << result << std::endl;
                    //     }
                        
                    //     if(boxTmp.isInside(result)){
                    //         if(debug){
                    //             log() << "is inside"<< std::endl;
                    //         }
                    //         isInsideBlockage = true;
                    //         break;
                            
                    //     }
                    // }
                    // if(isInsideBlockage || (pbx.layerIdx >= 1)){
                    //     continue;
                    // }
                    




                    // patch upper two layers
                    if (pbx.layerIdx < database.getLayerNum() - 2) {
                        guideGenStat.pinRegionPatchCand++;
                        auto bxPlus1 = gr::GrPoint(pbx.layerIdx + 1, pbx.x, pbx.y);
                        auto bxPlus2 = gr::GrPoint(pbx.layerIdx + 2, pbx.x, pbx.y);
                        if (needPatch(bxPlus1) || needPatch(bxPlus2)) {
                            patched = true;
                            guideGenStat.pinRegionPatchNum++;
                            grNet.patchRouteGuides.emplace_back(getSurrounding(bxPlus1));
                            grNet.patchRouteGuides.emplace_back(getSurrounding(bxPlus2));
                        } else {
                            grNet.patchRouteGuides.emplace_back(gr::GrBoxOnLayer(bxPlus1.layerIdx, {pbx.x}, {pbx.y}));
                            grNet.patchRouteGuides.emplace_back(gr::GrBoxOnLayer(bxPlus2.layerIdx, {pbx.x}, {pbx.y}));
                        }
                    }
                    // patch lower two layers
                    if (pbx.layerIdx > 1) {
                        guideGenStat.pinRegionPatchCand++;
                        auto bxMinus1 = gr::GrPoint(pbx.layerIdx - 1, pbx.x, pbx.y);
                        auto bxMinus2 = gr::GrPoint(pbx.layerIdx - 2, pbx.x, pbx.y);
                        if (needPatch(bxMinus1) || needPatch(bxMinus2)) {
                            patched = true;
                            guideGenStat.pinRegionPatchNum++;
                            grNet.patchRouteGuides.emplace_back(getSurrounding(bxMinus1));
                            grNet.patchRouteGuides.emplace_back(getSurrounding(bxMinus2));
                            if(debug){
                                log() << "surronding bxMinus1: " << getSurrounding(bxMinus1) << std::endl;
                                log() << "surronding bxMinus2: " << getSurrounding(bxMinus2) << std::endl;
                            }
                        } else {
                            grNet.patchRouteGuides.emplace_back(gr::GrBoxOnLayer(bxMinus1.layerIdx, {pbx.x}, {pbx.y}));
                            grNet.patchRouteGuides.emplace_back(gr::GrBoxOnLayer(bxMinus2.layerIdx, {pbx.x}, {pbx.y}));

                            if(debug){
                                log() << "bxMinus1: " << gr::GrBoxOnLayer(bxMinus1.layerIdx, {pbx.x}, {pbx.y}) << std::endl;
                                log() << "bxMinus2: " << gr::GrBoxOnLayer(bxMinus2.layerIdx, {pbx.x}, {pbx.y}) << std::endl;
                            }
                        }
                    }
                    // patch original layer
                    if (patched) {
                        grNet.patchRouteGuides.emplace_back(getSurrounding(pbx));
                        if(debug){
                                log() << "getSurrounding(pbx): " << getSurrounding(pbx) << std::endl;
                            }
                    } else {
                        grNet.patchRouteGuides.emplace_back(gr::GrBoxOnLayer(pbx.layerIdx, {pbx.x}, {pbx.y}));
                        if(debug){
                            log() << "pbx: " << gr::GrBoxOnLayer(pbx.layerIdx, {pbx.x}, {pbx.y}) << std::endl;
                        }
                    }
                }
            }
        }//end if isValidPatch
        
    }else{// advancded mode covering all layers
        
        for (auto &pbxs : grNet.pinAccessBoxes) {
            for (auto &pbx : pbxs) {
                guideGenStat.pinRegionPatchCand++;
                for(int layer_idx = 0; layer_idx < database.getLayerNum(); layer_idx++){
                    grNet.patchRouteGuides.emplace_back(gr::GrBoxOnLayer(layer_idx, {pbx.x}, {pbx.y}));
                }
                
            }//end for pbx
        }
                    
    }//end PatchPinRegionsMode
    
}//end patchPinRegions

void GuideGenerator::patchLongSegments() {
    // Note: assuming all guides are single width
    const int patchIntvl = 5;  // min patch interval
    const double patchThresh = 1.0;
    for (const auto &guide : grNet.wireRouteGuides) {
        int layerIdx = guide.layerIdx;
        auto dir = database.getLayerDir(layerIdx);
        int offset = grNet.dbNet.idx % patchIntvl;
        for (int cp = guide[1 - dir].low + offset; cp <= guide[1 - dir].high;) {
            int x = dir == X ? guide[dir].low : cp;
            int y = dir == X ? cp : guide[dir].low;
            if (grDatabase.getCellResource({layerIdx, x, y}) < patchThresh) {
                bool patched = false;
                for (int layer_delta = -1; layer_delta <= 1; layer_delta += 2) {
                    int layer = layerIdx + layer_delta;
                    if (layer < 1 || layer >= database.getLayerNum()) continue;
                    if (grDatabase.getCellResource({layer, x, y}) > patchThresh) {
                        grNet.patchRouteGuides.emplace_back(gr::GrBoxOnLayer(layer, {x}, {y}));
                        guideGenStat.longSegmentPatchNum++;
                        patched = true;
                    }
                }
                if (patched) {
                    cp += patchIntvl;
                } else {
                    cp++;
                }
            } else {
                cp++;
            }
        }
    }
}

void GuideGenerator::patchVioCells() {
    // Note: assuming all guides are single width
    const int queryWidth = 1;

    for (const auto &guide : grNet.wireRouteGuides) {
        int layerIdx = guide.layerIdx;
        auto dir = database.getLayerDir(layerIdx);
        int gridline = guide[dir].low;

        for (int cp = guide[1 - dir].low; cp <= guide[1 - dir].high; cp++) {
            int x = dir == X ? gridline : cp;
            int y = dir == X ? cp : gridline;

            double cellRsrc = grDatabase.getCellResource({layerIdx, x, y});

            if (cellRsrc <= 0) {
                guideGenStat.vioCellNum++;

                gr::GrBoxOnLayer patch(layerIdx, utils::IntervalT<int>(x), utils::IntervalT<int>(y));

                for (int g = gridline - queryWidth; g <= gridline + queryWidth; g++) {
                    if (g < 0 || g >= grDatabase.getNumGrPoint(dir)) continue;

                    double curcellRsrc = dir == X ? grDatabase.getCellResource({layerIdx, g, y})
                                                   : grDatabase.getCellResource({layerIdx, x, g});
                    if (curcellRsrc <= 0) continue;
                    cellRsrc += curcellRsrc;
                    patch[dir].Update(g);
                }

                if (cellRsrc > 0) {
                    const vector<int> layers = {layerIdx + 1, layerIdx - 1};
                    for (auto l : layers) {
                        if (l >= database.getLayerNum() || l <= 1) continue;

                        bool vioFree = true;
                        for (int x = patch.lx(); x <= patch.hx() && vioFree; x++)
                            for (int y = patch.ly(); y <= patch.hy() && vioFree; y++)
                                if (grDatabase.getCellResource({l, x, y}) <= 0) {
                                    vioFree = false;
                                }

                        if (vioFree) {
                            grNet.patchRouteGuides.push_back(patch);
                            grNet.patchRouteGuides.emplace_back(l, patch[X], patch[Y]);
                            guideGenStat.vioCellPatchNum++;
                            break;
                        }
                    }
                }
            }
        }
    }
}

void GuideGenerator::genPatchGuides() {
    grNet.patchRouteGuides.clear();
    if(db::setting.patchPinRegions)
        patchPinRegions();
    if(db::setting.patchLongSegments)
        patchLongSegments();
    if(db::setting.patchVioCells)   
        patchVioCells();
}

void GuideGenerator::genTopoGuides() {
    grNet.wireRouteGuides.clear();
    grNet.viaRouteGuides.clear();

    genConnGuides();

    wl_cost = 0;
    wl_abs_cost = 0;
    // Note: generate guides by topology
    grNet.postOrderVisitGridTopo([&](std::shared_ptr<gr::GrSteiner> node) {
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
                grNet.wireRouteGuides.emplace_back(lower->layerIdx,
                                                   utils::IntervalT<int>((*lower)[X], (*upper)[X]),
                                                   utils::IntervalT<int>((*lower)[Y], (*upper)[Y]));

                gr::GrEdge edge(gr::GrPoint(lower->layerIdx, lower->x , lower->y ),
                                gr::GrPoint(upper->layerIdx, upper->x, upper->y) );

                DBU edge_cost = grDatabase.getWireCost(edge);
                wl_cost += edge_cost;
                wl_abs_cost += grDatabase.getDist(*parent, *child);

            } else {
                grNet.viaRouteGuides.emplace_back(parent->layerIdx,
                                                  utils::IntervalT<int>((*parent)[X], (*parent)[X]),
                                                  utils::IntervalT<int>((*parent)[Y], (*parent)[Y]));
                grNet.viaRouteGuides.emplace_back(child->layerIdx,
                                                  utils::IntervalT<int>((*child)[X], (*child)[X]),
                                                  utils::IntervalT<int>((*child)[Y], (*child)[Y]));
            }
        }
    });

    sliceGuides(grNet.wireRouteGuides);
    sliceGuides(grNet.viaRouteGuides);

    // if(grNet.getName() == "net3136"){
    //     log() << "init viaRouteGuides: " << grNet.viaRouteGuides.size() << std::endl;
    //     log() << "init vias: " << std::endl;
    //     for(auto g : grNet.viaRouteGuides){
    //         log() << g << std::endl;
    //     }
    // }
}
