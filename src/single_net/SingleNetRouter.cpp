#include "SingleNetRouter.h"
#include "InitRoute.h"
#include "MazeRoute.h"

SingleNetRouter::SingleNetRouter(gr::GrNet& grDatabaseNet)
    : grNet(grDatabaseNet), guideGen(grNet), status(db::RouteStatus::SUCC_NORMAL) {
    grDatabase.removeNet(grNet);
    grNet.gridTopo.clear();
}

void SingleNetRouter::finish() {
    utils::timer net_timer;

    guideGen.genTopoGuides();
    grDatabase.useNet(grNet);

    grNet.dbNet.net_timer += net_timer.getTimer();

    // log() << "genTopoGuides grNet: " << grNet.getName()
    //       << ", guideGen wl_cost: " << guideGen.wl_cost
    //       << ", guideGen wl_abs_cost: " << guideGen.wl_abs_cost
    //        << std::endl;
}

void SingleNetRouter::mazeRoute() {
    utils::timer net_timer;

    if (!grNet.needToRoute()) {
        status = db::RouteStatus::SUCC_ONE_PIN;
        
    } else {
        MazeRoute mazeRouter(grNet,FINEGRID);
        mazeRouter.constructGridGraph(guides);
        status = mazeRouter.run();
        // grid graph stream for debugging purpose
        stream_str = mazeRouter.getGridGraphStreamString();
    }

    db::routeStat.increment(db::RouteStage::MAZE, status);

    grNet.dbNet.net_timer += net_timer.getTimer();
}

void SingleNetRouter::initRoutePattern(InitRoute& initRouter) {
    utils::timer net_timer;
    // log() << "initRoutePattern net: " << grNet.getName() << std::endl;
    if (!grNet.needToRoute()) {
        // log() << "needtoRoute failed ..." << std::endl;
        status = db::RouteStatus::SUCC_ONE_PIN;
    } else {
        // log() << "start patternRoute ..." << std::endl;
        initRouter.patternRouteMT();
        initRouter.buildTopo();
        status = db::RouteStatus::SUCC_NORMAL;
    }
    db::routeStat.increment(db::RouteStage::INIT, status);

    grNet.dbNet.net_timer += net_timer.getTimer();
}

void SingleNetRouter::planMazeRoute(const CongestionMap& congMap) {
    utils::timer net_timer;
    // run mazeroute on coarse grid
    const int cellWidth = congMap.getCellWidth();
    const int cellHeight = congMap.getCellHeight();

    gr::GrNet tmpNet = grNet;

    bool debug = false;

    // if(grNet.getName() == "net10214"){
    //     debug = true;
    //     log() << "maze: " << grNet.getName() << std::endl;
    // }    
    MazeRoute mazeRouter(tmpNet,COARSEGRID);
    mazeRouter.constructGridGraph(congMap);
    status = mazeRouter.run();

    stream_coarseGrid_str = mazeRouter.getCoarseGridGraphStreamString();


    // log() << "stream_coarseGrid_str" << stream_coarseGrid_str << std::endl;
    


    if(debug)
        log() << "maze: " << grNet.getName() << ", status: " << status << std::endl;

    // generate guides
    auto getLower = [&](int coor, Dimension dir) {
        if (dir == X)
            return coor * cellWidth;
        else
            return coor * cellHeight;
    };
    auto getUpper = [&](int coor, Dimension dir) {
        if (dir == X)
            return min((coor + 1) * cellWidth, grDatabase.getNumGrPoint(X)) - 1;
        else
            return min((coor + 1) * cellHeight, grDatabase.getNumGrPoint(Y)) - 1;
    };

    

    guides.clear();
    tmpNet.postOrderVisitGridTopo([&](std::shared_ptr<gr::GrSteiner> node) {
        auto parent = node;
        for (auto child : parent->children) {
            // if (tmpNet.getName() == "net6700") printlog(*parent, *child);
            if (parent->layerIdx == child->layerIdx) {
                std::shared_ptr<gr::GrSteiner> lower, upper;
                if ((*parent)[X] < (*child)[X] || (*parent)[Y] < (*child)[Y]) {
                    lower = parent;
                    upper = child;
                } else {
                    lower = child;
                    upper = parent;
                }
                guides.emplace_back(lower->layerIdx,
                                    utils::IntervalT<int>(getLower((*lower)[X], X), getUpper((*upper)[X], X)),
                                    utils::IntervalT<int>(getLower((*lower)[Y], Y), getUpper((*upper)[Y], Y)));
            } else {
                guides.emplace_back(parent->layerIdx,
                                    utils::IntervalT<int>(getLower((*parent)[X], X), getUpper((*parent)[X], X)),
                                    utils::IntervalT<int>(getLower((*parent)[Y], Y), getUpper((*parent)[Y], Y)));
                guides.emplace_back(child->layerIdx,
                                    utils::IntervalT<int>(getLower((*child)[X], X), getUpper((*child)[X], X)),
                                    utils::IntervalT<int>(getLower((*child)[Y], Y), getUpper((*child)[Y], Y)));
            }
        }
    });


    // maintain connectivity
    auto mergedPinAccessBoxes = grNet.getMergedPinAccessBoxes([&](const gr::GrPoint& point) {
        return gr::PointOnLayer(point.layerIdx, point[X] / cellWidth, point[Y] / cellHeight);
    });
    const int neighLayers = 2;
    for (auto& points : mergedPinAccessBoxes) {
        for (auto& point : points) {
            for (int l = point.layerIdx - neighLayers; l <= point.layerIdx + neighLayers; l++) {
                if (l < database.getLayerNum() && l >= 0)
                    guides.emplace_back(l,
                                        utils::IntervalT<int>(getLower(point[X], X), getUpper(point[X], X)),
                                        utils::IntervalT<int>(getLower(point[Y], Y), getUpper(point[Y], Y)));
            }
        }
    }

    grNet.dbNet.net_timer += net_timer.getTimer();
}
