

# from email.policy import default
import os
import sys
import pandas as pd
import numpy as np


import_path = os.path.abspath(os.path.join(os.path.join(__file__, ".."), ".."))
sys.path.insert(0, import_path)

from backend.pltcairo import *
from backend.utils import *
from backend.param import *

class Net:
    def __init__(self,db,type):
        self.db = db
        self.type = type

    def getWindow(self,window,plt_obj,net_name="-1",l=-1):
        if not(self.type in self.db):
            return
            
        db = self.db
        die_df = db["die"]
        
        net_df = db[self.type]
        # args = db["args"]

        

        # net_filter = net_df.loc[net_df.apply(lambda row: getIntervals(row.xl,row.xh,window[XL],window[XH]),axis=1)]
        # net_filter = net_filter.loc[net_filter.apply(lambda row: getIntervals(row.yl,row.yh,window[YL],window[YH]),axis=1)]
        

        
        net_filter = net_df.loc[ (net_df.xl >= window[XL] )& (net_df.xh <= window[XH])]
        net_filter = net_filter.loc[ (net_df.yl >= window[YL] )& (net_df.yh <= window[YH])]
       
        if(net_name != "-1"):
            net_filter = net_df.loc[ net_df.net_name == net_name]
        # only second layer
        if (l != -1):
            net_filter = net_filter.loc[net_filter.l == l]

        

        xls = net_filter.xl.values
        yls = net_filter.yl.values
        xhs = net_filter.xh.values
        yhs = net_filter.yh.values

        w = 0
        if(self.type == "drnet"):
            w = 100
        

        ws = [np.abs(xhs[i]-xls[i]+w) for i in np.arange(len(xls))]
        hs = [np.abs(yhs[i]-yls[i]+w) for i in np.arange(len(yls))]

        if(self.type == "net"):
            colors =[(0,1,0) for i in range(len(xls))]
            alphas =[0.5 for i in range(len(xls))]
            plt_obj.run(net_filter.xl.values,net_filter.yl.values,\
                        ws,hs,colors,alphas)
        elif(self.type =="drnet"):
            colors =[(0,0,1) for i in range(len(xls))]
            alphas =[1 for i in range(len(xls))]
            plt_obj.run(net_filter.xl.values,net_filter.yl.values,\
                        ws,hs,colors,alphas)
        elif(self.type =="netDRGuide"):
            colors =[(1,0,0) for i in range(len(xls))]
            alphas =[0.5 for i in range(len(xls))]
            plt_obj.run(net_filter.xl.values,net_filter.yl.values,\
                        ws,hs,colors,alphas)


        