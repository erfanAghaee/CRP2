

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
            dtype={"x":'float64',"y":'float64',"w":'float64',"h":"float64"})
    except:
        print("No cell input!")

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


    return db

    