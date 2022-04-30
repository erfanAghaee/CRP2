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

#include <iostream>
#include "global.h"
#include "FlexRoute.h"
#include "io/io.h"
#include "pa/FlexPA.h"
#include "ta/FlexTA.h"
#include "dr/FlexDR.h"
#include "gc/FlexGC.h"
#include "rp/FlexRP.h"
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace fr;

void FlexRoute::init() {
  if (DBPROCESSNODE == "GF14_13M_3Mx_2Cx_4Kx_2Hx_2Gx_LB") {
    VIAINPIN_BOTTOMLAYERNUM = 2;
    VIAINPIN_TOPLAYERNUM = 2;
    USENONPREFTRACKS = false;
    BOTTOM_ROUTING_LAYER = 4;
    TOP_ROUTING_LAYER = 18;
    ENABLE_VIA_GEN = false;
  }

  io::Parser parser(getDesign());
  parser.readLefDef();
  parser.readGuide();
  parser.postProcess();
  FlexPA pa(getDesign());
  pa.main();
  parser.postProcessGuide();
}

void FlexRoute::prep() {
  FlexRP rp(getDesign(), getDesign()->getTech());
  rp.main();
}

void FlexRoute::ta() {
  FlexTA ta(getDesign());
  ta.main();
  io::Writer writer(getDesign());
  writer.writeFromTA();
}

void FlexRoute::dr() {
  bool debug_erfan = false;
  if(debug_erfan) std::cout << "start DR erfan" << std::endl;
  FlexDR dr(getDesign());
  dr.main();
  dr.getReport();
  dr.getReportNets();
  // // dr.getCongestion();
  // dr.outOffGuideReport();
  
}

void FlexRoute::endFR() {
  io::Writer writer(getDesign());
  writer.writeFromDR();
}

void FlexRoute::checkGuideVsDR(){
  std::cout << "checkGuideVsDR" << std::endl;
} 



int FlexRoute::main() {
  bool debug_erfan = true;

  // get filter nets
  // filter_nets_set
  std::vector<std::string> filter_nets;
  boost::split(filter_nets, filter_nets_name, boost::is_any_of(","));

  for(auto tmp : filter_nets){
    if(tmp != "-1")
      filter_nets_set.insert(tmp);  
  }
  // 

  
  init();
  // if(debug_erfan) checkGuideVsDR();
  if(debug_erfan) {
    std::cout<<"DR_CRPFixedNetHelper: " << DR_CRPFixedNetHelper << std::endl;
    std::cout<<"defFixedNets: " << defFixedNets << std::endl;
  }

  bool debug = false;
  if(debug){
    FlexDR dr(getDesign());
    dr.logNets();
    dr.logCells();
  }else{
    prep();
    ta();
    // if(DR_CRPFixedNetHelper != "1" ){
    dr();
    // }
    
    endFR();
  }
  

  return 0;
}

