

import os
import sys
import pandas as pd
import numpy as np


import_path = os.path.abspath(os.path.join(os.path.join(__file__, ".."), ".."))
sys.path.insert(0, import_path)

from backend.pltcairo import *
from backend.utils import *
from backend.param import *

class DRC:
    def __init__(self,db):
        self.db = db


    def run(self):
        pass
        # db =  self.db
        # print(cells_df)
        # self.pltAllCells()
        # window = [320.3,574.2,333.2,579.0]
        # window = [x*2000 for x in window]
        # self.pltWindow(window)

    def getWindow(self,window,plt_obj,color,alpha,l=-1):
        if not("drc" in self.db):
            return

        db = self.db
        die_df = db["die"]
        drc_df = db["drc"]
        args = db["args"]
        
        # surface = plt_obj.init(window)
        drc_filter = drc_df.loc[drc_df.apply(lambda row: getIntervals(row.xl,row.xh,window[XL],window[XH]),axis=1)]
        drc_filter = drc_filter.loc[drc_filter.apply(lambda row: getIntervals(row.yl,row.yh,window[YL],window[YH]),axis=1)]
       
        drc_filter = drc_df.loc[ (drc_df.xl >= window[XL] )& (drc_df.xh <= window[XH])]
        drc_filter = drc_filter.loc[ (drc_df.yl >= window[YL] )& (drc_df.yh <= window[YH])]
        if(l != -1):
            drc_filter = drc_filter.loc[ drc_df.l == l]

        xls = drc_filter.xl.values
        yls = drc_filter.yl.values
        xhs = drc_filter.xh.values
        yhs = drc_filter.yh.values

        ws = [np.abs(xhs[i]-xls[i]) for i in np.arange(len(xls))]
        hs = [np.abs(yhs[i]-yls[i]) for i in np.arange(len(yls))]


        plt_obj.run(xls,yls,\
                    ws,hs,color,alpha)


    

    def getStatistics(self,window):
        if not("drc" in self.db):
            return

        db = self.db
        die_df = db["die"]
        drc_df = db["drc"]
        args = db["args"]
        
        # print(congestion_df)
        grps = drc_df.groupby(['l'])
        nums = drc_df.groupby(['l']).size()
        norm = [float(i)/max(nums) for i in nums]
        print(nums)

        layers = grps.groups.keys()

        plt.plot(layers,norm,label="drc.dr")
        
        # plt.plot(grps)
        

