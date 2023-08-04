

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

class Legalizer:
    def __init__(self,db,type_):
        self.db = db
        self.type = type_

    def getWindowSize(self):
        eps = 1000
        db = self.db
        die_df = db["die"]
        legalizer_df = db[self.type]
        args = db["args"]
        xl = np.min(legalizer_df["xl"].values)
        yl = np.min(legalizer_df["yl"].values)
        xh = np.max(legalizer_df["xh"].values)
        yh = np.max(legalizer_df["yh"].values)

        return [xl-eps,yl-eps,xh+eps,yh+eps]

    def getWindowBoard(self,window,plt_obj,text=False):
        db = self.db
        die_df = db["die"]
        legalizer_df = db[self.type]
        args = db["args"]
        
        # surface = plt_obj.init(window)

        # cell_filter = cell_df.loc[cell_df.apply(lambda row: getIntervals(row.x,row.x+row.w,window[XL],window[XH]),axis=1)]
        # cell_filter = cell_filter.loc[cell_filter.apply(lambda row: getIntervals(row.y,row.y+row.h,window[YL],window[YH]),axis=1)]

        legalizer_filter = legalizer_df.loc[ (legalizer_df.xl >= window[XL] )& (legalizer_df.xh <= window[XH])]
        legalizer_filter = legalizer_filter.loc[ (legalizer_filter.yl >= window[YL] )& (legalizer_filter.yh <= window[YH])]
        
        
        # cell_filter = cell_filter.loc[ (cell_df.y >= window[YL] )& (cell_df.y <= window[YH])]
        
        txts = [str(cost) for cost in legalizer_filter.cost.values]

        xls = legalizer_filter.xl.values
        yls = legalizer_filter.yl.values
        xhs = legalizer_filter.xh.values
        yhs = legalizer_filter.yh.values
        ws = [np.abs(xhs[i]-xls[i]) for i in np.arange(len(xls))]
        hs = [np.abs(yhs[i]-yls[i]) for i in np.arange(len(yls))]

        
        colors =[(0,0,1) for i in range(len(xls))]
        alphas =[0.2 for i in range(len(xls))]

        plt_obj.run(xls,yls,\
                    ws,hs,colors,alphas)

        if(text):
            plt_obj.drawText(txts,font=45)

    def getWindowLegalizer(self,window,plt_obj,text=False):
        db = self.db
        die_df = db["die"]
        legalizer_df = db[self.type]
        args = db["args"]
        
        # surface = plt_obj.init(window)

        # cell_filter = cell_df.loc[cell_df.apply(lambda row: getIntervals(row.x,row.x+row.w,window[XL],window[XH]),axis=1)]
        # cell_filter = cell_filter.loc[cell_filter.apply(lambda row: getIntervals(row.y,row.y+row.h,window[YL],window[YH]),axis=1)]

        legalizer_filter = legalizer_df.loc[ (legalizer_df.xl >= window[XL] )& (legalizer_df.xh <= window[XH])]
        legalizer_filter = legalizer_filter.loc[ (legalizer_filter.yl >= window[YL] )& (legalizer_filter.yh <= window[YH])]
        
        
        # cell_filter = cell_filter.loc[ (cell_df.y >= window[YL] )& (cell_df.y <= window[YH])]
        
        txts = [str(cost) for cost in legalizer_filter.cost.values]

        xls = legalizer_filter.xl.values
        yls = legalizer_filter.yl.values
        xhs = legalizer_filter.xh.values
        yhs = legalizer_filter.yh.values
        ws = [np.abs(xhs[i]-xls[i]) for i in np.arange(len(xls))]
        hs = [np.abs(yhs[i]-yls[i]) for i in np.arange(len(yls))]

        
        colors =[(1,0,1) for i in range(len(xls))]
        alphas =[0.1 for i in range(len(xls))]

        plt_obj.run(xls,yls,\
                    ws,hs,colors,alphas)

        if(text):
            plt_obj.drawText(txts,font=45)

    def getWindowLegalizerSol(self,window,plt_obj,text=False):
        db = self.db
        die_df = db["die"]
        legalizer_df = db[self.type]
        args = db["args"]
        
        # surface = plt_obj.init(window)

        # cell_filter = cell_df.loc[cell_df.apply(lambda row: getIntervals(row.x,row.x+row.w,window[XL],window[XH]),axis=1)]
        # cell_filter = cell_filter.loc[cell_filter.apply(lambda row: getIntervals(row.y,row.y+row.h,window[YL],window[YH]),axis=1)]

        legalizer_filter = legalizer_df.loc[ (legalizer_df.xl >= window[XL] )& (legalizer_df.xh <= window[XH])]
        legalizer_filter = legalizer_filter.loc[ (legalizer_filter.yl >= window[YL] )& (legalizer_filter.yh <= window[YH])]
        
        
        # cell_filter = cell_filter.loc[ (cell_df.y >= window[YL] )& (cell_df.y <= window[YH])]
        
        txts = [name for name in legalizer_filter.cell_name.values]

        xls = legalizer_filter.xl.values
        yls = legalizer_filter.yl.values
        xhs = legalizer_filter.xh.values
        yhs = legalizer_filter.yh.values
        ws = [np.abs(xhs[i]-xls[i]) for i in np.arange(len(xls))]
        hs = [np.abs(yhs[i]-yls[i]) for i in np.arange(len(yls))]

        
        colors =[(1,0,1) for i in range(len(xls))]
        alphas =[0.1 for i in range(len(xls))]

        plt_obj.run(xls,yls,\
                    ws,hs,colors,alphas)

        if(text):
            plt_obj.drawText(txts,font=45)
    

    def getWindow(self,window,plt_obj,text=False):
        if not(self.type in self.db):
            return

        if(self.type == "legalizerBoard"):
            self.getWindowBoard(window,plt_obj,text=text)
        elif(self.type == "legalizer"):
            self.getWindowLegalizer(window,plt_obj,text=text)
        elif(self.type == "legalizerSolution"):
            self.getWindowLegalizerSol(window,plt_obj,text=text)


        
        







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

    