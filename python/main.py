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
from backend.gcell import *


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

    
    cell_obj = Cell(db)
    cell_obj.run()

    gcell_obj = GCell(db)
    gcell_obj.run()



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