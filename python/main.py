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
from backend.gcell import *
from backend.pltcairo import *


def drawVio(plt_obj,color,alpha):
    # rect = [ 382.5, 583.55, 385.5*2, 585.05*2]
    rect = [639000,1149100,645000,1152100]
    w = np.abs(rect[XL]-rect[XH])
    h = np.abs(rect[YL]-rect[YH])
    rect = [rect[XL],rect[YH],w,h]
    # rect = [x*2000.0 for x in rect]
    xs=[]; xs.append(rect[XL])
    ys=[]; ys.append(rect[YL])
    ws=[]; ws.append(w)
    hs=[]; hs.append(h)

    plt_obj.run(xs,ys,\
                ws,hs,color,alpha)


def menu(args):
    print(args.dir)
    db = {}
    db["cell"] = pd.read_csv(args.dir + "cells/" + args.bench + ".cell.0.csv", \
        dtype={"x":'float64',"y":'float64',"w":'float64',"h":"float64"})
    db["die"] = pd.read_csv(args.dir + "misc/" + args.bench + ".die.csv",\
        dtype={"die_xl":'float64',"die_yl":'float64',"die_xh":'float64',"die_yh":"float64"} )
    db["args"] = args
    db["gcell"] = pd.read_csv(args.dir + "misc/" + args.bench + ".gcell.csv", \
        dtype={"l":"float64","x":'float64',"y":'float64',"w":'float64',"h":"float64"})

    db["net"] = pd.read_csv(args.dir + "nets/" + args.bench + ".net.0.csv", \
        dtype={"x":'float64',"y":'float64',"w":'float64',"h":"float64"})

    die_df = db["die"]

    for i in np.arange(5):

        # if i != 4:
        #     continue


        db["net"] = pd.read_csv(args.dir + "nets/" + args.bench + ".net."+str(i)+".csv", \
            dtype={"l": "float64","xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64"})

        db["vio"] = pd.read_csv(args.dir + "vios/" + args.bench + ".vio."+str(i)+".csv", \
            dtype={"l": "float64","xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64"})

        cell_obj = Cell(db)
        # cell_obj.run()

        gcell_obj = GCell(db)
        # gcell_obj.run()
        net_obj = Net(db)

        vio_obj = Vio(db)

        


        # 81.5
        # window = [315.30,574.2,333.2,589.5]
        # window = [81.5,574.2,333.2,589.5]
        # window = [x*2000 for x in window]
        # window = [642000,1164000,1371000,178000]
        # window = [417,582,445,589]
        window = [417,482,660,589]
        window = [x*2000 for x in window]


        # window = [
        #           die_df["die_xl"].values[0]
        #         , die_df["die_yl"].values[0]
        #         , die_df["die_xh"].values[0]
        #         , die_df["die_yh"].values[0]
        #     ]

        plt_obj = PltCairo()
        surface = plt_obj.init(window)
        cell_obj.getWindow(window,plt_obj,(0,0,1),1,text=False)
        gcell_obj.getWindow(window,plt_obj,(1,0,0),0.001)
        net_obj.getWindow("net60637",window,plt_obj,(0,1,0),0.8)
        vio_obj.getWindow(window,plt_obj,(1,0,0),1)
        
        # drawVio(plt_obj,(1,0,0),0.9)
        # net_obj.getNet("net2110",window,plt_obj,(0,1,0),1)


        # surface.write_to_png(args.dir+ "imgs/" +"net.window.png")  # Output to PNG
        surface.write_to_png(args.dir+ "imgs/" +"net.window."+str(i)+".png")  # Output to PNG


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