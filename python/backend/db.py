

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




def getDB(args,iter_gr,iter_dr):
    print(args.dir)
    db = {}
    try:
        db["cell"] = pd.read_csv(args.dir + "cells/" + args.bench + ".cell." \
            +str(iter_gr)+".csv", \
            dtype={"xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64"})
    except:
        print("No cell input!")

    try:
        db["cellcandidate"] = pd.read_csv(args.dir + "placement/" + args.bench + ".cells.candidates." \
            +str(iter_gr)+".csv", \
            dtype={"xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64","cost":"float64"})
    except:
        print("No cell candidate input!")
    
    try:
        db["cellmoved"] = pd.read_csv(args.dir + "placement/" + args.bench + ".cells.moved." \
            +str(iter_gr)+".csv", \
            dtype={"xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64",\
                "newXl":'float64',"newYl":'float64',"newXh":'float64',"newYh":"float64"})
    except:
        print("No cell moved input!")

    try:
        db["die"] = pd.read_csv(args.dir + "misc/" + args.bench + ".die.csv",\
        dtype={"die_xl":'float64',"die_yl":'float64',"die_xh":'float64',"die_yh":"float64"} )
    except:
        print("No die input!")

    try:
        db["fixedMetals"] = pd.read_csv(args.dir + "fixedmetals/" + args.bench + ".fixedMetals."\
        +str(iter_gr)+".csv",\
        dtype={"xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64"} )
    except:
        print("No fixedMetals input!")

    try:
        db["args"] = args
    except:
        print("No args input!")

    try:
        db["gcell"] = pd.read_csv(args.dir + "misc/" + args.bench + ".gcell.csv", \
        dtype={"l":"float64","x":'float64',"y":'float64',"w":'float64',"h":"float64"})
    except:
        print("No gcell input!")

    try:
        db["net"] = pd.read_csv(args.dir + "nets/" + args.bench + ".net."\
        +str(iter_gr)+".csv", \
        dtype={"x":'float64',"y":'float64',"w":'float64',"h":"float64"})
    except:
        print("No net input!")
    
    try:
        db["netDRGuide"] = pd.read_csv(args.dir + "DR/" + args.bench + ".dr.guide.csv", \
        dtype={"x":'float64',"":'float64',"w":'float64',"h":"float64"})
    except:
        print("No netDRGuide input!")

    try:
        db["vio"] = pd.read_csv(args.dir + "vios/" + args.bench + ".vio."\
        +str(iter_gr)+".csv", \
        dtype={"l": "float64","xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64"})
    except:
        print("No vio input!")


    try:
        db["congestion"] = pd.read_csv(args.dir + "congestions/" + args.bench + ".congestion.edge."\
        +str(iter_gr)+".csv", \
        dtype={"l": "float64","xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64",\
            "wireUsage":"float64","fixedUsage":"float64",\
                "viaUsage":"float64","numTracks":"float64",\
                    "overflow":"float64"})
    except:
        print("No congestion input!")

    try:
        db["drc"] = pd.read_csv(args.dir + "drc/" + args.bench + ".dr.drc."\
        +str(iter_dr)+".csv", \
        dtype={"l": "float64","xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64"})
    except:
        print("No drc input!")

    try:
        db["drnet"] = pd.read_csv(args.dir + "drnets/" + args.bench + ".dr.net."\
        +str(iter_dr)+".csv", \
        dtype={"l": "float64","xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64"})
    except:
        print("No drnet input!")
    

    try:
        db["patternroute"] = pd.read_csv(args.dir + "patternroute/" + args.bench + ".patternroute."\
        +str(iter_gr)+".csv", \
        dtype={"xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64",\
            "cost":"str"})
    except:
        print("No patternRoute input!")


    
    try:
        db["coef"] = pd.read_csv(args.dir + "coef/" + args.bench + ".coef."\
        +str(iter_gr)+".csv",\
            dtype={"name":"str","value":"float64"})
    except:
        print("No coef input!")

    try:
        db["astar"] = pd.read_csv(args.dir + "astar/" + args.bench + ".astar."\
        +str(iter_gr)+".csv",\
            dtype={"idx":"float64","net_name":"str",\
                "l":"float64","xl":"float64","yl":"float64","xh":"float64",\
                "hy":"float64","eBackward":"float64","eForward":"float64",\
                "eDown":"float64","eUp":"float64","eBackwardCost":"float64",\
                "eForwardCost":"float64","eDownCost":"float64","eUpCost":"float64"})
    except:
        print("No astar input!")

    try:
        db["astarcoarsegrid"] = pd.read_csv(args.dir + "astar/" + args.bench + ".astar.coarsegrid."\
        +str(iter_gr)+".csv",\
            dtype={"idx":"float64","net_name":"str",\
                "l":"float64","xl":"float64","yl":"float64","xh":"float64",\
                "hy":"float64","eBackward":"float64","eForward":"float64",\
                "eDown":"float64","eUp":"float64","eBackwardCost":"float64",\
                "eForwardCost":"float64","eDownCost":"float64","eUpCost":"float64"})
    except:
        print("No astar input!")


    try:
        db["legalizerBoard"] = pd.read_csv(args.dir + "legalizer/" + args.bench + ".legalizer.board.csv", \
            dtype={"xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64","cost":"float64"})
    except:
        print("No legalizerBoard input!")

    try:
        db["legalizer"] = pd.read_csv(args.dir + "legalizer/" + args.bench + ".legalizer.csv", \
            dtype={"xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64","cost":"float64"})
    except:
        print("No legalizer input!")

    try:
        db["legalizerWeight"] = pd.read_csv(args.dir + "legalizer/" + args.bench + ".legalizer.weights.csv", \
            dtype={"weight_idx":"float64",\
                "xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64","cost":"float64"})
    except:
        print("No legalizer weights input!")

    try:
        db["legalizerOverlaps"] = pd.read_csv(args.dir + "legalizer/" + args.bench + ".legalizer.overlaps.csv", \
            dtype={"weight_idx1":"float64",\
                "xl1":'float64',"yl1":'float64',"xh1":'float64',"yh1":"float64",\
                "weight_idx2":"float64",\
                "xl2":'float64',"yl2":'float64',"xh2":'float64',"yh2":"float64"})
    except:
        print("No legalizer overlaps input!")


    try:
        db["legalizerSolution"] = pd.read_csv(args.dir + "legalizer/" + args.bench + ".legalizer.solution.csv", \
            dtype={"weight_idx":"float64",\
                "xl":'float64',"yl":'float64',"xh":'float64',"yh":"float64"})
    except:
        print("No legalizer weights input!")



    return db

    