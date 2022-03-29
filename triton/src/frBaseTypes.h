/* Authors: Lutong Wang and Bangqi Xu */
/*
 * Copyright (c) 2019, The Regents of the University of California
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE REGENTS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef _FR_BASE_TYPES_H_
#define _FR_BASE_TYPES_H_

#include <vector>
#include <list>
#include <map>
#include <string>
#include <utility>
#include "iostream"

#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>

namespace fr {
  using frLayerNum = int;
  using frCoord = int;
  using frUInt4 = unsigned int;
  using frDist  = double;
  using frString = std::string;
  using frCost = unsigned int;
  using frMIdx = int; // negative value expected 
  template <typename T>
  using frCollection = std::vector<T>;
  template <typename T>
  using frVector     = std::vector<T>;
  template <typename T>
  using frList       = std::list<T>;
  template <typename T>
  using frListIter   = typename std::list<T>::iterator;

  enum frOrientEnum {
      frcR0       = 0, // N
      frcR90      = 1, // W
      frcR180     = 2, // S
      frcR270     = 3, // E
      frcMY       = 4, // FN
      frcMXR90    = 5, // FW
      frcMX       = 6, // FS
      frcMYR90    = 7  // FE
  };
  enum frEndStyleEnum {
      frcTruncateEndStyle = 0, // ext = 0
      frcExtendEndStyle   = 1, // ext = half width
      frcVariableEndStyle = 2  // ext = variable
  };
  enum frPrefRoutingDirEnum {
      frcNotApplicablePrefRoutingDir = 0,
      frcNonePrefRoutingDir          = 1,
      frcHorzPrefRoutingDir          = 2,
      frcVertPrefRoutingDir          = 3 
  };
  enum frBlockObjectEnum {
      frcNet,
      frcTerm,
      frcInst,
      frcVia,
      frcPin,
      frcInstTerm,
      frcRect,
      frcPolygon,
      frcSteiner,
      frcRoute,
      frcPathSeg,
      frcGuide,
      frcBlockage,
      frcLayerBlockage,
      frcBlock,
      frcBoundary,
      frcInstBlockage,
      frcAccessPattern,
      frcMarker,
      frcPatchWire,
      frcAccessPoint,
      frcPinAccess,
      frcCMap,
      frcGCellPattern,
      frcTrackPattern,
      drcNet,
      drcPin,
      drcAccessPattern,
      drcPathSeg,
      drcVia,
      drcMazeMarker,
      drcPatchWire,
      tacTrack,
      tacPin,
      tacPathSeg,
      tacVia,
      gccNet,
      gccPin,
      gccEdge,
      gccRect,
      gccPolygon
  };
  enum class frGuideEnum {
      frcGuideX,
      frcGuideGlobal,
      frcGuideTrunk,
      frcGuideShortConn
  };
  enum class frTermEnum {
    frcNormalTerm,
    frcClockTerm,
    frcPowerTerm,
    frcGroundTerm
  };
  enum class frNetEnum {
    frcNormalNet,
    frcClockNet,
    frcPowerNet,
    frcGroundNet
  };

  enum class frConstraintTypeEnum { // check FlexDR.h fixMode
    frcShortConstraint = 0,
    frcAreaConstraint = 1,
    frcMinWidthConstraint = 2,
    frcSpacingConstraint = 3,
    frcSpacingEndOfLineConstraint = 4,
    frcSpacingEndOfLineParallelEdgeConstraint = 5, // not supported
    frcSpacingTableConstraint = 6, // not supported
    frcSpacingTablePrlConstraint = 7,
    frcSpacingTableTwConstraint = 8,
    frcLef58SpacingTableConstraint = 9, // not supported
    frcLef58CutSpacingTableConstraint = 10, // not supported
    frcLef58CutSpacingTablePrlConstraint = 11, // not supported
    frcLef58CutSpacingTableLayerConstraint = 12, // not supported
    frcLef58CutSpacingConstraint = 13, // not supported
    frcLef58CutSpacingParallelWithinConstraint = 14, // not supported
    frcLef58CutSpacingAdjacentCutsConstraint = 15, // not supported
    frcLef58CutSpacingLayerConstraint = 16, // not supported
    frcCutSpacingConstraint = 17,
    frcMinStepConstraint,
    frcLef58MinStepConstraint,
    frcMinimumcutConstraint,
    frcOffGridConstraint,
    frcMinEnclosedAreaConstraint,
    frcLef58CornerSpacingConstraint, // not supported
    frcLef58CornerSpacingConcaveCornerConstraint, // not supported
    frcLef58CornerSpacingConvexCornerConstraint, // not supported
    frcLef58CornerSpacingSpacingConstraint, // not supported
    frcLef58CornerSpacingSpacing1DConstraint, // not supported
    frcLef58CornerSpacingSpacing2DConstraint, // not supported
    frcLef58SpacingEndOfLineConstraint, // not supported
    frcLef58SpacingEndOfLineWithinConstraint, // not supported
    frcLef58SpacingEndOfLineWithinEndToEndConstraint, // not supported
    frcLef58SpacingEndOfLineWithinParallelEdgeConstraint, // not supported
    frcLef58SpacingEndOfLineWithinMaxMinLengthConstraint, // not supported
    frcLef58CutClassConstraint, // not supported
    frcNonSufficientMetalConstraint,
    frcSpacingSamenetConstraint,
    frcLef58RightWayOnGridOnlyConstraint,
    frcLef58RectOnlyConstraint,
    frcRecheckConstraint
  };

  enum class frCornerTypeEnum {
    UNKNOWN,
    CONCAVE,
    CONVEX
  };

  enum class frCornerDirEnum {
    UNKNOWN,
    NE,
    SE,
    SW,
    NW
  };

  enum class frMinimumcutConnectionEnum {
    UNKNOWN = -1,
    FROMABOVE = 0,
    FROMBELOW = 1
  };

  enum class frMinstepTypeEnum {
    UNKNOWN = -1,
    INSIDECORNER = 0,
    OUTSIDECORNER = 1,
    STEP = 2
  };

  #define OPPOSITEDIR 7 // used in FlexGC_main.cpp
  enum class frDirEnum { UNKNOWN = 0, D = 1, S = 2, W = 3, E = 4, N = 5, U = 6 };

  enum class frLayerTypeEnum {
    CUT,
    ROUTING,
    IMPLANT,
    MASTERSLICE
  };


  enum class AccessPointTypeEnum {
    Ideal,
    Good,
    Offgrid,
    None
  };

  enum class MacroClassEnum {
    UNKNOWN,
    CORE,
    CORE_TIEHIGH,
    CORE_TIELOW,
    CORE_WELLTAP,
    CORE_SPACER,
    CORE_ANTENNACELL,
    COVER,
    ENDCAP_PRE,
    BLOCK,
    RING, // ispd19
    PAD, // ispd19
    PAD_POWER, // ispd19
    PAD_SPACER,
    ENDCAP_BOTTOMLEFT // ispd19
  };

  // note: FlexPA hardcoded the cost, don't change here
  enum class frAccessPointEnum {
    frcOnGridAP   = 0,
    frcHalfGridAP = 1,
    frcCenterAP   = 2,
    frcEncOptAP   = 3
  };

  namespace bg  = boost::geometry;

  typedef bg::model::d2::point_xy<frCoord, bg::cs::cartesian>  point_t;
  typedef bg::model::box<point_t>                              box_t;
  typedef bg::model::segment<point_t>                          segment_t;

  class frBox;

  template <typename T>
  using rq_box_value_t = std::pair<frBox, T>;
}

#endif
