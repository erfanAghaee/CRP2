

import os
import sys
from tkinter import font
import pandas as pd
import numpy as np


import_path = os.path.abspath(os.path.join(os.path.join(__file__, ".."), ".."))
sys.path.insert(0, import_path)

from backend.pltcairo import *
from backend.utils import *
from backend.param import *

class Cell:
    def __init__(self,db,type_):
        self.db = db
        self.type = type_

    def getWindow(self,window,plt_obj,cell_name = "-1",text=False):
        if not(self.type in self.db):
            return

        db = self.db
        die_df = db["die"]
        cell_df = db[self.type]
        args = db["args"]
        
        # surface = plt_obj.init(window)

        # cell_filter = cell_df.loc[cell_df.apply(lambda row: getIntervals(row.x,row.x+row.w,window[XL],window[XH]),axis=1)]
        # cell_filter = cell_filter.loc[cell_filter.apply(lambda row: getIntervals(row.y,row.y+row.h,window[YL],window[YH]),axis=1)]

        cell_filter = cell_df.loc[ (cell_df.xl >= window[XL] )& (cell_df.xh <= window[XH])]
        cell_filter = cell_filter.loc[ (cell_filter.yl >= window[YL] )& (cell_filter.yh <= window[YH])]
        
        if(cell_name != "-1"):
            cell_filter = cell_filter.loc[ cell_filter.cell_name == cell_name]
        # cell_filter = cell_filter.loc[ (cell_df.y >= window[YL] )& (cell_df.y <= window[YH])]
        
        txts = cell_filter.cell_name.values

        xls = cell_filter.xl.values
        yls = cell_filter.yl.values
        xhs = cell_filter.xh.values
        yhs = cell_filter.yh.values
        ws = [np.abs(xhs[i]-xls[i]) for i in np.arange(len(xls))]
        hs = [np.abs(yhs[i]-yls[i]) for i in np.arange(len(yls))]

        if(self.type == "cell"):
            colors =[(0,0,1) for i in range(len(xls))]
            alphas =[1 for i in range(len(xls))]

            plt_obj.run(xls,yls,\
                        ws,hs,colors,alphas)

            if(text):
                plt_obj.drawText(txts,font=12)
        
        elif(self.type == "cellcandidate"):
            colors =[(1,0,0) for i in range(len(xls))]
            alphas =[0.6 for i in range(len(xls))]
            plt_obj.run(xls,yls,\
            ws,hs,colors,alphas)






    # def getWindow(self,window,plt_obj,color,alpha,text=False):
    #     if not("cell" in self.db):
    #         return

    #     db = self.db
    #     die_df = db["die"]
    #     cell_df = db["cell"]
    #     args = db["args"]
        
    #     # surface = plt_obj.init(window)

    #     # cell_filter = cell_df.loc[cell_df.apply(lambda row: getIntervals(row.x,row.x+row.w,window[XL],window[XH]),axis=1)]
    #     # cell_filter = cell_filter.loc[cell_filter.apply(lambda row: getIntervals(row.y,row.y+row.h,window[YL],window[YH]),axis=1)]
        
    #     cell_filter = cell_df.loc[ (cell_df.x >= window[XL] )& (cell_df.x <= window[XH])]
        
    #     cell_filter = cell_filter.loc[ (cell_df.y >= window[YL] )& (cell_df.y <= window[YH])]
        
    #     txts = cell_filter.cell_name.values

    #     plt_obj.run(cell_filter.x.values,cell_filter.y.values,\
    #                 cell_filter.w.values,cell_filter.h.values,color,alpha)

    #     if(text):
    #         plt_obj.drawText(txts,font=42)
        
    #     # return surface

    