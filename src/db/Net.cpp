#include "Net.h"

#include <fstream>

#include "Setting.h"

namespace db {

NetBase::~NetBase() {}

BoxOnLayer NetBase::getMaxAccessBox(int pinIdx) const {
    DBU maxArea = std::numeric_limits<DBU>::min();
    db::BoxOnLayer bestBox;
    for (const auto& box : pinAccessBoxes[pinIdx]) {
        if (maxArea < box.area()) {
            maxArea = box.area();
            bestBox = box;
        }
    }
    return bestBox;
}

void NetBase::print(ostream& os) const {
    os << "Net " << getName() << " (idx = " << idx << ") with " << numOfPins() << " pins " << std::endl;
    for (int i = 0; i < numOfPins(); ++i) {
        os << "pin " << i << " " << rsynPins[i].getInstanceName() << std::endl;
        for (auto& accessBox : pinAccessBoxes[i]) {
            os << accessBox << std::endl;
        }
    }

    os << routeGuides.size() << " route guides" << std::endl;
    if (routeGuides.size() == gridRouteGuides.size()) {
        for (int i = 0; i < routeGuides.size(); ++i) {
            os << routeGuides[i] << " " << gridRouteGuides[i] << std::endl;
        }
    } else {
        for (auto& routeGuide : routeGuides) {
            os << routeGuide << std::endl;
        }
    }

    os << std::endl;
}

Net::Net(int i, Rsyn::Net net, RsynService& rsynService) {
    idx = i;
    rsynNet = net;
    net_timer = 0;

    // pins
    pinAccessBoxes.reserve(net.getNumPins());
    const Rsyn::Session session;
    const Rsyn::PhysicalDesign& physicalDesign =
        static_cast<Rsyn::PhysicalService*>(session.getService("rsyn.physical"))->getPhysicalDesign();
    const DBU libDBU = physicalDesign.getDatabaseUnits(Rsyn::LIBRARY_DBU);
    for (auto RsynPin : net.allPins()) {
        rsynPins.push_back(RsynPin);
        pinAccessBoxes.emplace_back();
        initPinAccessBoxes(RsynPin, rsynService, pinAccessBoxes.back(), libDBU);
    }
    bool debug = false;
    // if(net.getName() == "net275") debug = false;

    if(debug){
        log() << "initPinAccessBoxes...(Net)" << std::endl;
        for (int i = 0; i < numOfPins(); i++) {
            const auto& boxes = pinAccessBoxes[i];
            for(auto box : boxes){
                log() << "box: " << box << std::endl;
            }
        }
        
    }//end debug 


    // init route guides and detailed routing 
    if(db::setting.hasInputGuide){
            log() << "net: " << net.getName() << std::endl;
            Rsyn::NetGuide rsynGuide = rsynService.routeGuideService->getGuide(net);
            for(auto layer_guide : rsynGuide.allLayerGuides()){
                BoxOnLayer guide_box(layer_guide.getLayer().getRelativeIndex()
                                    ,layer_guide.getBounds().getLower().x
                                    ,layer_guide.getBounds().getLower().y
                                    ,layer_guide.getBounds().getUpper().x
                                    ,layer_guide.getBounds().getUpper().y);
                initRouteGuides.push_back(guide_box);
                // log() << "guide: " << layer_guide.getBounds().getLower().x
                //       << ", "<< layer_guide.getBounds().getLower().y
                //       << ", "<< layer_guide.getBounds().getUpper().x
                //       << ", "<< layer_guide.getBounds().getUpper().y 
                //       << ", "<< layer_guide.getLayer().getRelativeIndex() << std::endl;                
            }//end for 
            // int layerIdx =
            //         rsynService.physicalDesign.getPhysicalLayerByName(segment.clsLayerName).getRelativeIndex();
            // auto wire_descps = net.getWires();
            // for(auto wire : wire_descps){
            //     for(const DefWireSegmentDscp& seg : wire.clsWireSegments ){
            //         std::cout << "layer: " << seg.clsLayerName << std::endl;
            //         std::cout << "clsRoutingPoints size: " << seg.clsRoutingPoints.size() << std::endl;
            //         if(seg.clsRoutingPoints.size() == 2){//wire
            //             for(auto pts : seg.clsRoutingPoints){
            //                 std::cout << "pos.x: " << pts.clsPos.x 
            //                         << ", pos.y: " << pts.clsPos.y << std::endl;
            //                 std::cout << "lx: " << pts.clsRect.getLower().x
            //                         << ", ly: " << pts.clsRect.getLower().y
            //                         << ", hx: " << pts.clsRect.getUpper().x
            //                         << ", hy: " << pts.clsRect.getUpper().y << std::endl;
            //             }

            //         }else if (seg.clsRoutingPoints.size() == 1){//via

            //         }// end if-else seg.clsRoutingPoints.size()
                    
            //     }//end wireSegments
            // }//end wire_descps

            auto wire_descps = net.getWires();
            std::cout << "wire_descps size: " << wire_descps.size() << std::endl;
            for (const DefWireDscp& wire : wire_descps) {
                std::cout << "wire.clsWireSegments size: " << wire.clsWireSegments.size() << std::endl;
                for (const DefWireSegmentDscp& segment : wire.clsWireSegments) {
                    int layerIdx =
                        rsynService.physicalDesign.getPhysicalLayerByName(segment.clsLayerName).getRelativeIndex();
                    const DBU width = segment.clsRoutedWidth;
                    DBUxy pos;
                    DBU ext = 0;
                    std::cout << "segment.clsRoutingPoints: " << segment.clsRoutingPoints.size() << std::endl;
                    for (unsigned i = 0; i != segment.clsRoutingPoints.size(); ++i) {
                        const DefRoutingPointDscp& pt = segment.clsRoutingPoints[i];
                        const DBUxy& nextPos = pt.clsPos;
                        const DBU nextExt = pt.clsHasExtension ? pt.clsExtension : 0;
                        if (i >= 1) {
                            for (unsigned dim = 0; dim != 2; ++dim) {
                                if (pos[dim] == nextPos[dim]) continue;
                                const DBU l = pos[dim] < nextPos[dim] ? pos[dim] - ext : nextPos[dim] - nextExt;
                                const DBU h = pos[dim] < nextPos[dim] ? nextPos[dim] + nextExt : pos[dim] + ext;
                                BoxOnLayer box(layerIdx);
                                box[dim].Set(l, h);
                                box[1 - dim].Set(pos[1 - dim] - width / 2, pos[1 - dim] + width / 2);
                                // fixedMetalVec.emplace_back(box, OBS_NET_IDX);
                                initRouteDR.push_back(box);
                                log() << "wire_box: " << box << std::endl;
                                // ++numSNetObs;
                                break;
                            }
                        }
                        pos = nextPos;
                        ext = nextExt;
                        log() << "pt: " << nextPos << std::endl;
                        log() << "has extension: " <<pt.clsHasExtension << std::endl;
                        log() << "pt.clsHasVia: " << pt.clsHasVia << std::endl;
                        log() << "pt.clsHasRectangle: " << pt.clsHasRectangle << std::endl;
                        if(pt.clsHasRectangle){
                            log() << "rect: " << pt.clsRect << std::endl;
                        }
                        if (!pt.clsHasVia) continue;
                        const Rsyn::PhysicalVia& via = rsynService.physicalDesign.getPhysicalViaByName(pt.clsViaName);
                        const int botLayerIdx = via.getBottomLayer().getRelativeIndex();
                        for (const Rsyn::PhysicalViaGeometry& geo : via.allBottomGeometries()) {
                            Bounds bounds = geo.getBounds();
                            bounds.translate(pos);
                            const BoxOnLayer box(botLayerIdx, getBoxFromRsynBounds(bounds));
                            log() << "via_box: " << box << std::endl;
                            // fixedMetalVec.emplace_back(box, OBS_NET_IDX);
                            // ++numSNetObs;
                        }
                        const int topLayerIdx = via.getTopLayer().getRelativeIndex();
                        for (const Rsyn::PhysicalViaGeometry& geo : via.allTopGeometries()) {
                            Bounds bounds = geo.getBounds();
                            bounds.translate(pos);
                            const BoxOnLayer box(topLayerIdx, getBoxFromRsynBounds(bounds));
                            log() << "via_toplayer_box: " << box << std::endl;
                            // fixedMetalVec.emplace_back(box, OBS_NET_IDX);
                            // ++numSNetObs;
                        }
                        if (via.hasViaRule()) {
                            const utils::PointT<int> numRowCol =
                                via.hasRowCol() ? utils::PointT<int>(via.getNumCols(), via.getNumRows())
                                                : utils::PointT<int>(1, 1);
                            BoxOnLayer botBox(botLayerIdx);
                            BoxOnLayer topBox(topLayerIdx);
                            for (unsigned dimIdx = 0; dimIdx != 2; ++dimIdx) {
                                const Dimension dim = static_cast<Dimension>(dimIdx);
                                const DBU origin = via.hasOrigin() ? pos[dim] + via.getOrigin(dim) : pos[dim];
                                const DBU botOff =
                                    via.hasOffset() ? origin + via.getOffset(Rsyn::BOTTOM_VIA_LEVEL, dim) : origin;
                                const DBU topOff =
                                    via.hasOffset() ? origin + via.getOffset(Rsyn::TOP_VIA_LEVEL, dim) : origin;
                                const DBU length =
                                    (via.getCutSize(dim) * numRowCol[dim] + via.getSpacing(dim) * (numRowCol[dim] - 1)) / 2;
                                const DBU botEnc = length + via.getEnclosure(Rsyn::BOTTOM_VIA_LEVEL, dim);
                                const DBU topEnc = length + via.getEnclosure(Rsyn::TOP_VIA_LEVEL, dim);
                                botBox[dim].Set(botOff - botEnc, botOff + botEnc);
                                topBox[dim].Set(topOff - topEnc, topOff + topEnc);
                            }
                            log() << "botBox: " << botBox << std::endl;
                            log() << "topBox: " << topBox << std::endl;
                            // fixedMetalVec.emplace_back(botBox, OBS_NET_IDX);
                            // fixedMetalVec.emplace_back(topBox, OBS_NET_IDX);
                            // numSNetObs += 2;
                        }
                        if (layerIdx == botLayerIdx)
                            layerIdx = topLayerIdx;
                        else if (layerIdx == topLayerIdx)
                            layerIdx = botLayerIdx;
                        else {
                            // log() << "Error: Special net " << specialNet.getNet().clsName << " via " << pt.clsViaName
                            //     << " on wrong layer " << layerIdx << std::endl;
                            break;
                        }
                    }
                }
            }//end loop wiredescp

            
        }//end if 

        // log() << "init guide boxes net: " << getName() << std::endl;
        // for(auto guide_box : initRouteGuides){
        //     log() << "guide: " << guide_box << std::endl;
        // }

    // checkObidenceGrVsDr(rsynService);


}//end Net constructor 


void Net::initPinAccessBoxes(Rsyn::Pin rsynPin,
                             RsynService& rsynService,
                             vector<BoxOnLayer>& accessBoxes,
                             const DBU libDBU) {
    // PhysicalPort
    if (rsynPin.isPort()) {
        Rsyn::PhysicalPort phPort = rsynService.physicalDesign.getPhysicalPort(rsynPin.getPort());
        getPinAccessBoxes(phPort, accessBoxes);
        return;
    }

    // PhysicalLibraryPin
    Rsyn::PhysicalLibraryPin phLibPin = rsynService.physicalDesign.getPhysicalLibraryPin(rsynPin);

    // PhysicalCell
    Rsyn::Instance instance = rsynPin.getInstance();
    if (instance.getType() != Rsyn::CELL) {
        log() << "Warning: pin is not on a cell " << rsynPin.getNetName() << " " << rsynPin.getInstanceName()
              << std::endl;
        return;
    }
    Rsyn::Cell cell = instance.asCell();
    Rsyn::PhysicalCell phCell = rsynService.physicalDesign.getPhysicalCell(cell);
    Rsyn::PhysicalLibraryCell phLibCell = rsynService.physicalDesign.getPhysicalLibraryCell(cell);
    const DBUxy origin(static_cast<DBU>(std::round(phLibCell.getMacro()->originX() * libDBU)),
                       static_cast<DBU>(std::round(phLibCell.getMacro()->originY() * libDBU)));

    // fill accessBoxes
    getPinAccessBoxes(phLibPin, phCell, accessBoxes, origin);
};

void Net::updatePinAccessBoxes(RsynService& rsynService) {
    auto net = rsynNet;
    pinAccessBoxes.clear();
    // pins
    pinAccessBoxes.reserve(net.getNumPins());
    const Rsyn::Session session;
    const Rsyn::PhysicalDesign& physicalDesign =
        static_cast<Rsyn::PhysicalService*>(session.getService("rsyn.physical"))->getPhysicalDesign();
    const DBU libDBU = physicalDesign.getDatabaseUnits(Rsyn::LIBRARY_DBU);
    for (auto RsynPin : net.allPins()) {
        rsynPins.push_back(RsynPin);
        pinAccessBoxes.emplace_back();
        initPinAccessBoxes(RsynPin, rsynService, pinAccessBoxes.back(), libDBU);
    }
    bool debug = false;
    // if(net.getName() == "net275") debug = false;

    if(debug){
        log() << "updatePinAccessBoxes...(Net)" << std::endl;
        for (int i = 0; i < numOfPins(); i++) {
            const auto& boxes = pinAccessBoxes[i];
            for(auto box : boxes){
                log() << "box: " << box << std::endl;
            }
        }
        
    }
}

void Net::getPinAccessBoxes(Rsyn::PhysicalPort phPort, vector<BoxOnLayer>& accessBoxes) {
    auto displacement = phPort.getPosition();
    auto bounds = phPort.getBounds();
    Bounds dummyCellBounds(displacement, displacement);
    Rsyn::PhysicalTransform transform(dummyCellBounds, phPort.getOrientation());
    bounds.translate(displacement);
    bounds = transform.apply(bounds);
    accessBoxes.emplace_back(phPort.getLayer().getRelativeIndex(), getBoxFromRsynBounds(bounds));
}

void Net::getPinAccessBoxes(Rsyn::PhysicalLibraryPin phLibPin,
                            Rsyn::PhysicalCell phCell,
                            vector<BoxOnLayer>& accessBoxes,
                            const DBUxy& origin) {
    if (!phLibPin.hasPinGeometries()) {
        log() << "Warning: pin of " << phCell.getName() << " has no pinGeometries" << std::endl;
        return;
    }

    const DBUxy displacement = phCell.getPosition() + origin;
    auto transform = phCell.getTransform();
    // for (Rsyn::PhysicalPinGeometry phPinGeo : phLibPin.allPinGeometries()) {
    // TODO: check why multiple PinGeometry on 8t4 inst60849
    auto phPinGeo = phLibPin.allPinGeometries()[0];
    for (Rsyn::PhysicalPinLayer phPinLayer : phPinGeo.allPinLayers()) {
        if (!phPinLayer.hasRectangleBounds()) {
            log() << "Warning: pin has no RectangleBounds" << std::endl;
            continue;
        }
        int layerIdx = phPinLayer.getLayer().getRelativeIndex();
        for (auto bounds : phPinLayer.allBounds()) {
            bounds.translate(displacement);
            bounds = transform.apply(bounds);
            accessBoxes.emplace_back(layerIdx, getBoxFromRsynBounds(bounds));
        }
    }
}//end getPinAccessBoxes



void NetList::init(RsynService& rsynService) {
    bool debug = false;
    nets.clear();
    if (db::setting.dbVerbose >= +db::VerboseLevelT::MIDDLE) {
        log() << "Init NetList ..." << std::endl;
    }
    nets.reserve(rsynService.design.getNumNets());
    int numPins = 0;
    log() << "check init guie solution " << std::endl;
    for (Rsyn::Net net : rsynService.module.allNets()) {
        switch (net.getUse()) {
            case Rsyn::POWER:
                continue;
            case Rsyn::GROUND:
                continue;
            default:
                break;
        }
        nets.emplace_back(nets.size(), net, rsynService);
        numPins += nets.back().pinAccessBoxes.size();
    }
    if (debug) {
        for (const auto& net : nets) {
            for (const auto& accessBoxes : net.pinAccessBoxes) {
                for (const auto& box : accessBoxes) {
                    log() << "net: " << net.getName() << " " << box << std::endl; 
                }
            }
        }
    }//end debug
        

    if (setting.dbVerbose >= +db::VerboseLevelT::MIDDLE) {
        log() << "The number of nets is " << nets.size() << std::endl;
        log() << "The number of pins is " << numPins << std::endl;
        log() << std::endl;
    }
}

}  // namespace db
