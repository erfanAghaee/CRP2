#!/usr/bin/env python

import argparse
import os
from pydoc import describe
import sys
import pandas as pd

'''
Python Project Source File imports 
'''

import_path = os.path.abspath(os.path.join(os.path.join(__file__, ".."), ".."))
sys.path.insert(0, import_path)


from backend.cell import *
from backend.net import *
from backend.vio import *
from backend.drc import *
from backend.db import *
from backend.coef import *
from backend.astar import *
from backend.score import *
from backend.patternroute import *
from backend.congestionEdge import *
from backend.fixedMetals import *
from backend.utils import *
from backend.gcell import *
from backend.pltcairo import *
import matplotlib.pyplot as plt

def plotBox(plt_obj,boxs,color,alpha,type):
    try:
        if len(boxs[0]) != 4:
            return
    except:
        return

    xls = [box[XL] for box in boxs]
    yls = [box[YL] for box in boxs]
    xhs = [box[XH] for box in boxs]
    yhs = [box[YH] for box in boxs]

    ws = [np.abs(xhs[i]-xls[i]) for i in np.arange(len(xls))]
    hs = [np.abs(yhs[i]-yls[i]) for i in np.arange(len(yls))]


    plt_obj.run(xls,yls,ws,hs,color,alpha)

def mainAnimation(args):
    # specific_net = "net60637"
    specific_net = "net29885"
    cell_name = "inst8879"

    # iteration of gr
    for l in np.arange(9):
        for iter_gr in np.arange(4):
            for iter_dr in np.arange(4):

                if not(iter_gr == 0 and (l==0) \
                    and (iter_dr == 0)):
                    continue

                db = getDB(args,iter_gr=iter_gr,iter_dr=iter_dr)

                cells_candidate_df = db["cellcandidate"]
                cell_grps = cells_candidate_df.groupby(["cell_name"])
                

                for cell_name in cell_grps.groups:

                    if(len(cell_grps.get_group(cell_name)) <= 1 ):
                        continue

                    die_df = db["die"]
                    net_DRGuide_obj = Net(db,"netDRGuide")
                    fixedMetals_obj = FixedMetals(db)
                    cell_obj = Cell(db,"cell")
                    cell_candidate_obj = Cell(db,"cellcandidate")
                    gcell_obj = GCell(db)
                    net_obj = Net(db,"net")
                    vio_obj = Vio(db)
                    drnet_obj = Net(db,"drnet")
                    drc_obj = DRC(db)
                    congestion_obj = CongestionEdge(db)
                    patternroute_obj = PatternRoute(db,"patternroute")
                    score_obj = Score(db,"score")
                    # 81.5
                    # window = [315.30,574.2,333.2,589.5]
                    # window = [81.5,574.2,333.2,589.5]
                    # window = [x*2000 for x in window]
                    # window = [642000,1164000,1371000,178000]
                    # window = [417,582,445,589]
                    # window = [417,562,660,589]
                    # window = [x*2000 for x in window]
                    # window = [684,397,689,399]
                    # window = [x*2000 for x in window]


                    window = [
                          die_df["die_xl"].values[0]
                        , die_df["die_yl"].values[0]
                        , die_df["die_xh"].values[0]
                        , die_df["die_yh"].values[0]
                    ]


                    plt_obj = PltCairo()
                    surface = plt_obj.init(window)

                    cell_obj.getWindow(window,plt_obj,text=True)
                    cell_candidate_obj.getWindow(window,plt_obj,cell_name=cell_name)
                    
                    # net_obj.getWindow(window,plt_obj,l=l)
                    # vio_obj.getWindow(window,plt_obj,l=l)
                    # fixedMetals_obj.getWindow(window,plt_obj,l=l)
                    # drc_obj.getWindow(window,plt_obj,l=l)
                    # # # net_DRGuide_obj.getWindow(window,plt_obj,(0.8,0.1,0.1),0.8,net_name=specific_net)
                    # drnet_obj.getWindow(window,plt_obj,l=l)
                    # congestion_obj.getWindow(window,plt_obj,l=l,text=True)
                    # patternroute_obj.getWindow(window,plt_obj,(0,1,0),0.5,net_name=specific_net,l=l)
                    # gcell_obj.getWindow(window,plt_obj,(1,0,0),0.001)
                    # score_obj.getWindow(window,plt_obj,(1,0,0),1,net_name=specific_net,l=l)
                    
                

                    surface.write_to_png(args.dir+ "imgs/" + cell_name  
                        +".gr."+str(iter_gr)
                        +".dr."+str(iter_dr)
                        +".l."+str(l)+".png")

            #     break
            # break
    # OFGW = score_obj.getOutOfGuideTotal("wire")
    # OFGV = score_obj.getOutOfGuideTotal("via")
    
    # www = score_obj.getWrongWayWiring()

    # oftw = score_obj.getOffTrackTotal("wire")
    # oftv = score_obj.getOffTrackTotal("via")

    # short = score_obj.getShortAreaTotal()
    # minArea = score_obj.getMinAreaTotal()

    # it is correctly calculated 
    # wl,vias = score_obj.getWirelengthViasTotal()    
    # print("wl(wirelength):",wl)
    # print("vias:",vias)


    # print("OFGW(Out of guide wirelength):",OFGW)
    # print("OFGV(Out of guide Vias):",OFGV)
    # print("WWW (Wrong Way Wiring):",www)
    # print("oftw(off track wiring):",oftw)
    # print("oftv(off track via insertion):",oftv)
    # print("shortArea: ",short)
    # print("# MinArea: ",minArea)
    

        
def mainCoef(args): 
    coef_obj = Coef(args,num_iter=4)
    coef_obj.drawLogisticSlope()
    coef_obj.drawViaCost()
    coef_obj.drawWireCost()

    
def mainAstar(args):
    # specific_net = "net60637"
    specific_net = "net1237"

    directory = args.dir+"imgs/"+"astar"
    if not os.path.exists(directory):
        os.mkdir(directory)

   

    # iteration of gr
    for i in np.arange(4):
        for l in np.arange(8):

            # if i != 3:
            #     continue

            db = getDB(args,iter_gr=i,iter_dr=0)

            die_df = db["die"]
            net_DRGuide_obj = Net(db,"netDRGuide")
            fixedMetals_obj = FixedMetals(db)
            cell_obj = Cell(db)
            gcell_obj = GCell(db)
            net_obj = Net(db,"net")
            vio_obj = Vio(db)
            drnet_obj = Net(db,"drnet")
            drc_obj = DRC(db)
            congestion_obj = CongestionEdge(db)
            patternroute_obj = PatternRoute(db,"patternroute")
            astar_obj = Astar(db,"astar")
            astar_coarsegrid_obj = Astar(db,"astarcoarsegrid")
            # 81.5
            # window = [315.30,574.2,333.2,589.5]
            # window = [81.5,574.2,333.2,589.5]
            # window = [x*2000 for x in window]
            # window = [642000,1164000,1371000,178000]
            # window = [417,582,445,589]
            # window = [417,562,660,589]
            # window = [x*2000 for x in window]


            window = [
                  die_df["die_xl"].values[0]
                , die_df["die_yl"].values[0]
                , die_df["die_xh"].values[0] 
                , die_df["die_yh"].values[0] 
            ]


            plt_obj = PltCairo()
            surface = plt_obj.init(window)

            cell_obj.getWindow(window,plt_obj,(0,0,1),1,text=False)
            
            # net_obj.getWindow(window,plt_obj,(0,1,0),0.5,net_name=specific_net,l=l)
            net_obj.getWindow(window,plt_obj,(0,1,0),0.5,l=l)
            vio_obj.getWindow(window,plt_obj,(0.6,0.1,0.6),0.8,l=l)
            fixedMetals_obj.getWindow(window,plt_obj,(0,0,0),1,l=0)
            drc_obj.getWindow(window,plt_obj,(1,0,0),1,l=l)
            # # # net_DRGuide_obj.getWindow(window,plt_obj,(0.8,0.1,0.1),0.8,net_name=specific_net)
            # # # drnet_obj.getWindow(window,plt_obj,(0,0,1),1,net_name=specific_net)
            # congestion_obj.getWindow(window,plt_obj,(0.5,0.1,0.7),0.1,l=l)
            # # # # patternroute_obj.getWindow(window,plt_obj,(0,1,0),0.5,net_name=specific_net,l=l)
            # gcell_obj.getWindow(window,plt_obj,(1,0,0),0.001)
            # astar_obj.getWindow(window,plt_obj,(0,1,1),0.4,net_name=specific_net,l=l)
            # astar_coarsegrid_obj.getWindow(window,plt_obj,(0.8,0.5,0.3),0.2,net_name=specific_net,l=l)

            surface.write_to_png(args.dir+ "imgs/astar/" +"gridgraph." +specific_net+"."
                +str(i)+ ".l." + str(l)+".png")


def mainCongestion(args):
    
    specific_net = "net1237"

    directory = args.dir+"imgs/"+"congestion"
    if not os.path.exists(directory):
        os.mkdir(directory)



    # iteration of gr
    for i in np.arange(4):
        if i!=3:
            continue
        

        for j in np.arange(1):

            db = getDB(args,iter_gr=i,iter_dr=j)

            die_df = db["die"]
            net_DRGuide_obj = Net(db,"netDRGuide")
            fixedMetals_obj = FixedMetals(db)
            cell_obj = Cell(db)
            gcell_obj = GCell(db)
            net_obj = Net(db,"net")
            vio_obj = Vio(db)
            drnet_obj = Net(db,"drnet")
            drc_obj = DRC(db)
            congestion_obj = CongestionEdge(db)
            patternroute_obj = PatternRoute(db,"patternroute")
            astar_obj = Astar(db,"astar")
            astar_coarsegrid_obj = Astar(db,"astarcoarsegrid")
            # 81.5
            # window = [315.30,574.2,333.2,589.5]
            # window = [81.5,574.2,333.2,589.5]
            # window = [x*2000 for x in window]
            # window = [642000,1164000,1371000,178000]
            # window = [417,582,445,589]
            # window = [417,562,660,589]
            # window = [x*2000 for x in window]


            window = [
                    die_df["die_xl"].values[0]
                , die_df["die_yl"].values[0]
                , die_df["die_xh"].values[0] 
                , die_df["die_yh"].values[0] 
            ]


            plt_obj = PltCairo()
            surface = plt_obj.init(window)

            # cell_obj.getWindow(window,plt_obj,(0,0,1),1,text=False)
            
            # # net_obj.getWindow(window,plt_obj,(0,1,0),0.5,net_name=specific_net,l=l)
            # net_obj.getWindow(window,plt_obj,(0,1,0),0.5,l=l)
            # vio_obj.getWindow(window,plt_obj,(0.6,0.1,0.6),0.8,l=l)
            # fixedMetals_obj.getWindow(window,plt_obj,(0,0,0),1,l=0)
            drc_obj.getStatistics(window)
            # # # net_DRGuide_obj.getWindow(window,plt_obj,(0.8,0.1,0.1),0.8,net_name=specific_net)
            # # # drnet_obj.getWindow(window,plt_obj,(0,0,1),1,net_name=specific_net)
            congestion_obj.getStatistics(window)
            # # # # patternroute_obj.getWindow(window,plt_obj,(0,1,0),0.5,net_name=specific_net,l=l)
            # gcell_obj.getWindow(window,plt_obj,(1,0,0),0.001)
            # astar_obj.getWindow(window,plt_obj,(0,1,1),0.4,net_name=specific_net,l=l)
            # astar_coarsegrid_obj.getWindow(window,plt_obj,(0.8,0.5,0.3),0.2,net_name=specific_net,l=l)

            # surface.write_to_png(args.dir+ "imgs/astar/" +"gridgraph." +specific_net+"."
            #     +str(i)+ ".l." + str(l)+".png")
            plt.xlabel("layers")
            plt.ylabel("Norms")
            plt.legend()
            plt.savefig(args.dir+ "imgs/congestion/congestionMap.via.png")


       
        
def debugCUGRISPD2019Test5(args):
    # specific_net = "net60637"
    specific_net = "net29885"
    cell_name = "inst8879"

    # iteration of gr
    for l in np.arange(6):
        for iter_gr in np.arange(4):

            if(iter_gr != 0):
                continue
            
            db = getDB(args,iter_gr=iter_gr,iter_dr=0)

            
                

            

            die_df = db["die"]
            fixedMetals_obj = FixedMetals(db)
            cell_obj = Cell(db,"cell")
            
            gcell_obj = GCell(db)
            net_obj = Net(db,"net")
            vio_obj = Vio(db)
            gcell_obj = GCell(db)
            
            congestion_obj = CongestionEdge(db)
            # patternroute_obj = PatternRoute(db,"patternroute")
            # score_obj = Score(db,"score")
               


            window = [
                    die_df["die_xl"].values[0]
                , die_df["die_yl"].values[0]
                , die_df["die_xh"].values[0]
                , die_df["die_yh"].values[0]
            ]


            plt_obj = PltCairo()
            surface = plt_obj.init(window)

            cell_obj.getWindow(window,plt_obj,text=False)
            # gcell_obj.getWindow(window,plt_obj)

            # 
                    
                    
            
                    # vio_obj.getWindow(window,plt_obj,l=l)
            fixedMetals_obj.getWindow(window,plt_obj,l=l)
            net_obj.getWindow(window,plt_obj,l=l)
                    # drc_obj.getWindow(window,plt_obj,l=l)
                    # # # net_DRGuide_obj.getWindow(window,plt_obj,(0.8,0.1,0.1),0.8,net_name=specific_net)
                    # drnet_obj.getWindow(window,plt_obj,l=l)
                    # congestion_obj.getWindow(window,plt_obj,l=l,text=True)
                    # patternroute_obj.getWindow(window,plt_obj,(0,1,0),0.5,net_name=specific_net,l=l)
                    # gcell_obj.getWindow(window,plt_obj,(1,0,0),0.001)
                    # score_obj.getWindow(window,plt_obj,(1,0,0),1,net_name=specific_net,l=l)
                    
                

            surface.write_to_png(args.dir+ "imgs/" + cell_name  
                +".gr."+str(iter_gr)
                +".dr."+str(0)
                +".l."+str(l)+".png")


    # print how net connects
    surface = plt_obj.init(window)
    cell_obj.getWindow(window,plt_obj,text=False)
    fixedMetals_obj.getWindow(window,plt_obj)
    net_obj.getWindow(window,plt_obj)
    surface.write_to_png(args.dir+ "imgs/" +"net.png")
            # break

            #     break
            # break
    # OFGW = score_obj.getOutOfGuideTotal("wire")
    # OFGV = score_obj.getOutOfGuideTotal("via")
    
    # www = score_obj.getWrongWayWiring()

    # oftw = score_obj.getOffTrackTotal("wire")
    # oftv = score_obj.getOffTrackTotal("via")

    # short = score_obj.getShortAreaTotal()
    # minArea = score_obj.getMinAreaTotal()

    # it is correctly calculated 
    # wl,vias = score_obj.getWirelengthViasTotal()    
    # print("wl(wirelength):",wl)
    # print("vias:",vias)


    # print("OFGW(Out of guide wirelength):",OFGW)
    # print("OFGV(Out of guide Vias):",OFGV)
    # print("WWW (Wrong Way Wiring):",www)
    # print("oftw(off track wiring):",oftw)
    # print("oftv(off track via insertion):",oftv)
    # print("shortArea: ",short)
    # print("# MinArea: ",minArea)
def debugGCELL(args):
    # specific_net = "net60637"
    specific_net ="n_5371"
    # specific_net =  "n_5371"
    cell_name = "inst8879"


    # iteration of gr
    print(args.dir,args.bench)
    
            
    db = getDB(args,iter_gr=4,iter_dr=0)

            
                

            

    die_df = db["die"]
    fixedMetals_obj = FixedMetals(db)
    cell_obj = Cell(db,"cell")
    net_obj = Net(db,"net")
    net_DRGuide_obj = Net(db,"netDRGuide")
            
    gcell_obj = GCell(db)

    # window = [
    #     112*1000,
    #     900*1000,
    #     125*1000,
    #     908*1000
    # ]

    window = [
            die_df["die_xl"].values[0]-6000
        , die_df["die_yl"].values[0]-6000
        , die_df["die_xh"].values[0]+6000
        , die_df["die_yh"].values[0]+6000
    ]


    plt_obj = PltCairo()
    surface = plt_obj.init(window)
    

    cell_obj.getWindow(window,plt_obj,text=False)
    fixedMetals_obj.getWindow(window,plt_obj)
    # gcell_obj.getWindow(window,plt_obj,l=0)
    # net_obj.getWindow(window,plt_obj,net_name=specific_net)
    # net_DRGuide_obj.getWindow(window,plt_obj,net_name=specific_net)
    
            
                    # drc_obj.getWindow(window,plt_obj,l=l)
                    # # # net_DRGuide_obj.getWindow(window,plt_obj,(0.8,0.1,0.1),0.8,net_name=specific_net)
                    # drnet_obj.getWindow(window,plt_obj,l=l)
                    # congestion_obj.getWindow(window,plt_obj,l=l,text=True)
                    # patternroute_obj.getWindow(window,plt_obj,(0,1,0),0.5,net_name=specific_net,l=l)
                    # gcell_obj.getWindow(window,plt_obj,(1,0,0),0.001)
                    # score_obj.getWindow(window,plt_obj,(1,0,0),1,net_name=specific_net,l=l)
                    
                

    surface.write_to_png(args.dir+ "imgs/GCELL.png")


    

def drawBenchmarks(args):
    # specific_net = "net60637"
    specific_net = "net29885"
    cell_name = "inst8879"

    # iteration of gr
    
    
            
    db = getDB(args,iter_gr=0,iter_dr=0)

    die_df = db["die"]
    fixedMetals_obj = FixedMetals(db)
    cell_obj = Cell(db,"cell")
    
    gcell_obj = GCell(db)
    net_obj = Net(db,"net")
    vio_obj = Vio(db)
    
    congestion_obj = CongestionEdge(db)
    # patternroute_obj = PatternRoute(db,"patternroute")
    # score_obj = Score(db,"score")
        


    window = [
            die_df["die_xl"].values[0]
        , die_df["die_yl"].values[0]
        , die_df["die_xh"].values[0]
        , die_df["die_yh"].values[0]
    ]


    plt_obj = PltCairo()
    surface = plt_obj.init(window)

    cell_obj.getWindow(window,plt_obj,text=False)
            
            
    
            # vio_obj.getWindow(window,plt_obj,l=l)
    fixedMetals_obj.getWindow(window,plt_obj)
    # net_obj.getWindow(window,plt_obj,l=l)
            # drc_obj.getWindow(window,plt_obj,l=l)
            # # # net_DRGuide_obj.getWindow(window,plt_obj,(0.8,0.1,0.1),0.8,net_name=specific_net)
            # drnet_obj.getWindow(window,plt_obj,l=l)
            # congestion_obj.getWindow(window,plt_obj,l=l,text=True)
            # patternroute_obj.getWindow(window,plt_obj,(0,1,0),0.5,net_name=specific_net,l=l)
            # gcell_obj.getWindow(window,plt_obj,(1,0,0),0.001)
            # score_obj.getWindow(window,plt_obj,(1,0,0),1,net_name=specific_net,l=l)
                    
                

    surface.write_to_png(args.dir+ "imgs/" + cell_name  
        +".gr."+str(0)
        +".dr."+str(0)
        +".l."+str(0)+".png")
