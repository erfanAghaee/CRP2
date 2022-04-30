from functools import total_ordering
from genericpath import exists
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import sys



'''
Python Project Source File imports 
'''

import_path = os.path.abspath(os.path.join(os.path.join(__file__, ".."), ".."))
sys.path.insert(0, import_path)


table = {
    "run_exp144.sh" : "run_exp40.sh",
    "run_exp145.sh" : "run_exp88.sh",
    "run_exp146.sh" : "run_exp90.sh",
    "run_exp147.sh" : "run_exp91.sh"
}
BL="run_exp145.sh"
EH="run_exp144.sh"
CRP1="run_exp146.sh"
CRP10="run_exp147.sh"


scenarios = {
    BL    : "BL",
    EH    : "Ehplacer",
    CRP1  : "CRP1",
    CRP10 : "CRP10"
}

exist ={
    BL    : False,
    EH    : False,
    CRP1  : False,
    CRP10 : False
}


dir = "../reports/"

def getDiff(x,y):
    return abs(x-y)

def getMovs(base,new):
    result = base.merge(new,on="cell_name")
    result["diff_x"] = result.apply(lambda row: abs(row.xl_x -row.xl_y),axis=1)
    result["diff_y"] = result.apply(lambda row: abs(row.yl_x -row.yl_y),axis=1)
    result["moved"] = result.apply(lambda row: (row.diff_x > 0 and row.diff_y > 0),axis=1)
    # result["diff_y"] = result.apply(lambda x: )
    return result

def getWL(xl,yl,xh,yh):
    return abs(xl-xh) + abs(yl-yh)


def getNets(tb_nets):
    num_layers = max(tb_nets[BL]["l"].values) + 1
    print(num_layers)

    df = pd.DataFrame()
    df["l"] = range(num_layers)

    # print(df)

    cols = [BL,EH,CRP1,CRP10]
    # cols = [BL]
    # cols = []
    # print("cols")
    # for col in cols_total:
    #     print(col)
    #     if(exist[col]):
    #         cols.append(col)



    
    for col in cols:
        if(exist[col]):
            df_tmp = tb_nets[col]
            vals = []
        
            for l in range(num_layers):
                # df_tmp = df_tmp.loc[(df_tmp.l == l) & (df_tmp.type == "wire")]
                df_filter = df_tmp.loc[(df_tmp.l == l) ]
                
                wl = df_filter.apply(lambda row: getWL(row["xl"],row["yl"],row["xh"],row["yh"]) ,axis=1)

                # print("layer: ", l , wl.values)
                # vals.append(len(df_tmp.loc[(df_tmp.l == l) & (df_tmp.type == "wire")]))
                vals.append(np.sum(wl.values))
            #     print(np.sum(wl.values))

            df[col] = vals
        else:
            df[col] = -1

    return df 
    #     print("new l: ", l , ", numWires: ", len(new.loc[(new.l == l) & (new.type == "wire")]))


def collectData(bench):
    tbs_cells = dict()
    tbs_nets = dict()
    for table_key in table.keys():
        print(table_key)
        try:
            file = dir + table_key + "/" + bench + "/drnets/" + bench + ".dr.cells.0.csv"
            tbs_cells[table_key] = pd.read_csv(file,dtype={"xl":"float64","yl":"float64","xh":"float64","yh":"float64"})
            file = dir + table_key + "/" + bench + "/drnets/" + bench + ".dr.net.0.csv"
            tbs_nets[table_key] = pd.read_csv(file,dtype={"xl":"float64","yl":"float64","xh":"float64","yh":"float64",\
                "l":"int"})
            exist[table_key]=True
        except:
            print(table_key," doesn't exist!")

    return tbs_cells,tbs_nets


def getAvg(vals):
    num=0; total=0
    for val in vals:
        if(val > 0):
            total += val
            num+=1

    return total/num


    


def logMovs(tbs_cells,scenario):
    d = {"avg_x":[0],"avg_y":[0],"max_x":[0],"max_y":[0],"moved":[0],"total":[0],}
    df = pd.DataFrame(d)

    res = getMovs(tbs_cells[BL],tbs_cells[scenario])

    df["avg_x"] = getAvg(res["diff_x"].values)
    df["avg_y"] = getAvg(res["diff_y"].values)
    df["max_x"] = np.max(res["diff_x"].values)
    df["max_y"] = np.max(res["diff_y"].values)
    df["moved"] = np.sum(res["moved"].values)
    df["total"] = np.sum(len(res["moved"]))
    # print(df)
    df.to_csv(args.dir + "mov."+scenarios[scenario]+".csv", index=False)



    # print("ehplacer moved: ",np.sum(res["moved"]), ", total: ",len(res))

    # res = getMovs(tbs_cells[BL],tbs_cells[CRP1])
    # print("CRP1 moved: ",np.sum(res["moved"]), ", total: ",len(res))

    # res = getMovs(tbs_cells[BL],tbs_cells[CRP10])
    # print("CRP10 moved: ",np.sum(res["moved"]), ", total: ",len(res))


def menu(args):
    benchs = ["ispd18_test1","ispd18_test2","ispd18_test3","ispd18_test4","ispd18_test5","ispd18_test6",\
       "ispd18_test7","ispd18_test8","ispd18_test9","ispd18_test10"]

    if not os.path.exists(args.dir):
        os.mkdir(args.dir)

    tbs_cells,tbs_nets = collectData(args.bench)

    res = getNets(tbs_nets)
    

    if(exist[EH]):
        logMovs(tbs_cells,EH)
    if(exist[CRP1]):
        logMovs(tbs_cells,CRP1)
    if(exist[CRP10]):
        logMovs(tbs_cells,CRP10)
    
    res = res.rename(columns={"run_exp145.sh":"BL","run_exp144.sh":"Eh?Placer",\
         "run_exp146.sh":"CRP1", "run_exp147.sh":"CRP10"})

    res.to_csv(args.dir + "congestion.csv", index=False)

    
    


    


if __name__ == "__main__":
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
    
    
