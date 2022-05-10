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
from backend.congestionEdge import *
from backend.fixedMetals import *
from backend.outofguide import *
from backend.utils import *
from backend.gcell import *
from backend.pltcairo import *
from backend.animation import *
from backend.score import *


def drawCRP(args):
    specific_net = "net60637"

    print(args.dir)
    db = {}
    db["cell"] = pd.read_csv(args.dir + "cells/" + args.bench + ".cell.0.csv", \
        dtype={"x":'float64',"y":'float64',"w":'float64',"h":"float64"})
    db["die"] = pd.read_csv(args.dir + "misc/" + args.bench + ".die.csv",\
        dtype={"die_xl":'float64',"die_yl":'float64',"die_xh":'float64',"die_yh":"float64"} )

    db["fixedMetals"] = pd.read_csv(args.dir + "misc/" + args.bench + ".fixedMetals.csv",\
        dtype={"xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64"} )


    db["args"] = args
    db["gcell"] = pd.read_csv(args.dir + "misc/" + args.bench + ".gcell.csv", \
        dtype={"l":"float64","x":'float64',"y":'float64',"w":'float64',"h":"float64"})

    db["net"] = pd.read_csv(args.dir + "nets/" + args.bench + ".net.0.csv", \
        dtype={"x":'float64',"y":'float64',"w":'float64',"h":"float64"})
    
    db["netDRGuide"] = pd.read_csv(args.dir + "DR/" + args.bench + ".dr.guide.csv", \
        dtype={"x":'float64',"y":'float64',"w":'float64',"h":"float64"})

    die_df = db["die"]

    net_DRGuide_obj = Net(db,"netDRGuide")
    


    fixedMetals_obj = FixedMetals(db)


    # iteration of gr
    for i in np.arange(5):

        if i != 4:
            continue

        


        db["net"] = pd.read_csv(args.dir + "nets/" + args.bench + ".net."+str(i)+".csv", \
            dtype={"l": "float64","xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64"})
        

        db["vio"] = pd.read_csv(args.dir + "vios/" + args.bench + ".vio."+str(i)+".csv", \
            dtype={"l": "float64","xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64"})



        db["congestion"] = pd.read_csv(args.dir + "congestions/" + args.bench + ".congestion.edge."+str(i)+".csv", \
            dtype={"l": "float64","xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64",\
                "wireUsage":"float64","fixedUsage":"float64",\
                    "viaUsage":"float64","numTracks":"float64",\
                        "overflow":"float64"})

        cell_obj = Cell(db)
        # cell_obj.run()

        gcell_obj = GCell(db)
        # gcell_obj.run()
        net_obj = Net(db,"net")
        

        vio_obj = Vio(db)

        congestion_obj = CongestionEdge(db)


        


        # 81.5
        # window = [315.30,574.2,333.2,589.5]
        # window = [81.5,574.2,333.2,589.5]
        # window = [x*2000 for x in window]
        # window = [642000,1164000,1371000,178000]
        # window = [417,582,445,589]
        window = [417,562,660,589]
        window = [x*2000 for x in window]


        # window = [
        #           die_df["die_xl"].values[0]
        #         , die_df["die_yl"].values[0]
        #         , die_df["die_xh"].values[0]
        #         , die_df["die_yh"].values[0]
        #     ]

        
        



        # surface.write_to_png(args.dir+ "imgs/" +"net.window.png")  # Output to PNG
        l = 3
        # iteration of dr
        for j in np.arange(15):
            plt_obj = PltCairo()
            surface = plt_obj.init(window)
            # cell_obj.getWindow(window,plt_obj,(0,0,1),1,text=False)
            # gcell_obj.getWindow(window,plt_obj,(1,0,0),0.001)

            
            net_obj.getWindow(window,plt_obj,(0,1,0),0.8,net_name=specific_net)
            # net_obj.getWindow(window,plt_obj,(0,1,0),0.8,l=l)
            
            vio_obj.getWindow(window,plt_obj,(1,0,0),0.8,l=l)
            fixedMetals_obj.getWindow(window,plt_obj,(0,0,0),1)
            # net_DRGuide_obj.getWindow(specific_net,window,plt_obj,(1,0,0),0.5)
            
            db["drc"] = pd.read_csv(args.dir + "drc/" + args.bench + ".dr.drc."+str(j)+".csv", \
                dtype={"l": "float64","xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64"})
            db["drnet"] = pd.read_csv(args.dir + "drnets/" + args.bench + ".dr.net."+str(j)+".csv", \
                dtype={"l": "float64","xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64"})
            drnet_obj = Net(db,"drnet")
            drc_obj = DRC(db)

            # drc_obj.getWindow(window,plt_obj,(0.9,1,0.1),1,l=l)
            drc_obj.getWindow(window,plt_obj,(1,0,0),1)
            
            drnet_obj.getWindow(window,plt_obj,(0,0,1),1)
            congestion_obj.getWindow(window,plt_obj,(0.5,0.1,0.7),0.1,l=l)
            surface.write_to_png(args.dir+ "imgs/" +"net.gr."+str(i)+".dr."+str(j)+".png")  # Output to PNG
    # surface.write_to_png(args.dir+ "imgs/" +"net.window."+str(i)+".png")  # Output to PNG



def drawGuideCompare(args):
    # specific_net = "net60637"
    specific_net = "net28713"

    print(args.dir)
    db = {}
    db["cell"] = pd.read_csv(args.dir + "cells/" + args.bench + ".cell.0.csv", \
        dtype={"x":'float64',"y":'float64',"w":'float64',"h":"float64"})
    db["die"] = pd.read_csv(args.dir + "misc/" + args.bench + ".die.csv",\
        dtype={"die_xl":'float64',"die_yl":'float64',"die_xh":'float64',"die_yh":"float64"} )

    db["fixedMetals"] = pd.read_csv(args.dir + "misc/" + args.bench + ".fixedMetals.csv",\
        dtype={"xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64"} )


    db["args"] = args
    db["gcell"] = pd.read_csv(args.dir + "misc/" + args.bench + ".gcell.csv", \
        dtype={"l":"float64","x":'float64',"y":'float64',"w":'float64',"h":"float64"})

    db["net"] = pd.read_csv(args.dir + "nets/" + args.bench + ".net.0.csv", \
        dtype={"xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64"})
    
    db["netDRGuide"] = pd.read_csv(args.dir + "DR/" + args.bench + ".dr.guide.csv", \
        dtype={"xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64"})

    die_df = db["die"]

    net_DRGuide_obj = Net(db,"netDRGuide")
    


    fixedMetals_obj = FixedMetals(db)

    i=4
    # compareTwoRoute(db["net"],db["netDRGuide"]) 
    # return


    db["net"] = pd.read_csv(args.dir + "nets/" + args.bench + ".net."+str(i)+".csv", \
        dtype={"l": "float64","xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64"})

    db["vio"] = pd.read_csv(args.dir + "vios/" + args.bench + ".vio."+str(i)+".csv", \
        dtype={"l": "float64","xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64"})

    cell_obj = Cell(db)
        # cell_obj.run()

    gcell_obj = GCell(db)
        # gcell_obj.run()
    net_obj = Net(db,"net")
    net_obj_preprocessing = Net(db,"netDRGuide")
        

    vio_obj = Vio(db)

        


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
    # gcell_obj.getWindow(window,plt_obj,(0,0,0),0.001)
    net_obj.getWindow(specific_net,window,plt_obj,(0,1,0),0.8)
    surface.write_to_png(args.dir+ "imgs/" +"netCUGR.png")  # Output to PNG
        # vio_obj.getWindow(window,plt_obj,(1,0,0),1)
        # fixedMetals_obj.getWindow(window,plt_obj,(0,0,0),1)

    plt_obj = PltCairo()
    surface = plt_obj.init(window)
    cell_obj.getWindow(window,plt_obj,(0,0,1),1,text=False)
    # gcell_obj.getWindow(window,plt_obj,(1,0,0),0.001)
    net_obj_preprocessing.getWindow(specific_net,window,plt_obj,(1,0,0),0.8)
    surface.write_to_png(args.dir+ "imgs/" +"netDRPreprocessing.png")  # Output to PNG
       

def interval(row,wind):
    return pd.Interval(row.xl,row.xh).overlaps(pd.Interval(wind[0], wind[1]))

def menu(args):
    debugGCELL(args)
    # drawBenchmarks(args)
    # debugCUGRISPD2019Test5(args)
    # mainAnimation(args)
    # drawCRP(args)
    # drawGuideCompare(args)
    # mainCoef(args)
    # mainAstar(args)
    # mainCongestion(args)
    # outofguideInit(args)
    
    # mainScore(args)
    # outofguideModuleTest(args)
    
    
    
    
    

    


if __name__ == "__main__":
    print("Python program starts...")
    parser = argparse.ArgumentParser("python main.py", \
                            description="A program to excute for debugging purposes.", \
                            add_help=False, formatter_class=argparse.RawTextHelpFormatter)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")

    required.add_argument("-dir", "--dir", type=str, required=True, \
                            help="\t path to read benchmarks.", metavar="\b")
    required.add_argument("-bench", "--bench", type=str, required=True, \
                            help="\t path to read benchmarks.", metavar="\b")
    # optional.add_argument("-debug_output", "--DebugOutput", action="store_true", \
    #                         help="\t Store debug output json file.")
    args   = parser.parse_args()
    menu(args)