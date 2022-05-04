import pandas as pd
import numpy as np
import xlsxwriter
import argparse


import os
import sys


import_path = os.path.abspath(os.path.join(os.path.join(__file__, ".."), ".."))
sys.path.insert(0, import_path)

from backend.raw_score import RawScore
from backend.cost_funcion import CostFunction
from backend.weights import Weight


flow_naming = {
    "PL: Contest; GR: Contest; DR: DRCU" : "Con+Con+DRCU",
    "PL: CRP(10); GR: CRP(10); DR: Triton (path_cost; PR last iteration of GR is applied)" : "Con+CRP10(PathCost)+Triton",
    "PL: CRP(1); GR: CRP(1); DR: Triton (path_cost; PR last iteration of GR is applied)" : "Con+CRP1(PathCost)+Triton",
    "PL: Contest; GR: CUGR; DR: Triton" : "Con+CUGR+Triton",
    "PL: Contest; GR: Contest; DR: Triton" : "Con+Con+Triton",
    "PL: Eh?Placer; GR: CUGR; DR: Triton" : "Con+EhPlacer+CUGR+Triton",
    "PL: DreamPlace; GR: CUGR; DR: Triton" : "Con+DreamPlacer+CUGR+Triton",
    "PL: ILP-Based; GR: ILP-based; DR: Triton" : "Con+ILP-based+Triton",
    "PL: CRP(10); GR: CRP(10); DR: Triton (wl_abs_cost; PR first iteration of GR is applied)" : "Con+CRP10(WLCost)+Triton",
    "PL: Contest; GR: Contest; DR: LuLuRoute" : "Con+Con+LuLu",
    "PL: Contest; GR: CRP-FPM; DR: Triton (Improve UnitShortCost-Pitch of the layers)" : "Con+CRP-FPM+Triton",
    "PL: Contest; GR: CUGR_ORIG; DR: Triton" : "Con+CUGR_ORG+Triton"
}




def raw_scores(args):
    ispd18 = pd.read_csv(args.dir + args.bench+"_scores.csv")
    benchs = ispd18.groupby(["benchmark"])
    writer = pd.ExcelWriter('ispd18_scores.xlsx', engine='xlsxwriter')

    all_tables = ['wl_gr', 'vias_gr', 'time_gr', 'mem_gr',
        'wl_dr', 'vias_dr', 'OFGW', 'OFGV', 'OFTW', 'OFTV', 'WWW', 'short_area',
        'min_area', 'spacing', 'score', 'time_dr', 'mem_dr']
    
    raw_score = RawScore(flow_naming=flow_naming)

    for col in all_tables:
        print(col)
        raw_score.run(benchs=benchs,metric=col,writer=writer)    
        
    writer.save()

def cost_function_behavior(args):
    cost_obj = CostFunction(args)
    cost_obj.run()

def weight_table(args):
    weigth = Weight(args)
    weigth.run()
def menu(args):
    raw_scores(args)
    # cost_function_behavior(args)
    # weight_table(args)
    


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
