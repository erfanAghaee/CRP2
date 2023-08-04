

import os
import sys
import pandas as pd
import numpy as np


import_path = os.path.abspath(os.path.join(os.path.join(__file__, ".."), ".."))
sys.path.insert(0, import_path)

from backend.pltcairo import *
from backend.utils import *
from backend.param import *

class GCell:
    def __init__(self,db):
        self.db = db


    def run(self):
        # db =  self.db
        # print(cells_df)
        self.pltAllGCells()
        window = [320.3,574.2,333.2,579.0]
        window = [x*2000 for x in window]
        self.pltWindow(window)

    def getWindow(self,window,plt_obj,l=-1):
        if not("gcell" in self.db):
            return

        db = self.db
        die_df = db["die"]
        cell_df = db["gcell"]
        # args = db["args"]
        
        # surface = plt_obj.init(window)

        # cell_filter = cell_df.loc[cell_df.apply(lambda row: getIntervals(row.x,row.x + row.w,window[XL],window[XH]),axis=1)]
        # cell_filter = cell_filter.loc[cell_filter.apply(lambda row: getIntervals(row.y,row.y+row.h,window[YL],window[YH]),axis=1)]

        
        cell_filter = cell_df.loc[ (cell_df.x >= window[XL] )& (cell_df.x <= window[XH])]
        
        cell_filter = cell_filter.loc[ (cell_df.y >= window[YL] )& (cell_df.y <= window[YH])]
        if(l!= -1):
            cell_filter = cell_filter.loc[cell_filter.l == l]

        xs_txt = cell_filter.gcellX.values
        ys_txt = cell_filter.gcellY.values
        
        txts = [ "("+str(xs_txt[i]) +","+ str(ys_txt[i]) + ")" \
            for i in np.arange(len(cell_filter.gcellX.values))]

        colors =[(1,0,0) for i in range(len(cell_filter.x.values))]
        alphas =[0.09 for i in range(len(cell_filter.x.values))]

        plt_obj.run(cell_filter.x.values,cell_filter.y.values,\
                    cell_filter.w.values,cell_filter.h.values,colors,alphas)

        # plt_obj.drawText(txts,font=72,up=True)
        
        # return surface



    def pltWindow(self,window):
        db = self.db
        plt_obj = PltCairo()
        die_df = db["die"]
        cell_df = db["gcell"]
        args = db["args"]
        
        surface = plt_obj.init(window)

        
        cell_filter = cell_df.loc[ (cell_df.x >= window[XL] )& (cell_df.x <= window[XH])]
        
        cell_filter = cell_filter.loc[ (cell_df.y >= window[YL] )& (cell_df.y <= window[YH])]
        
        # txts = cell_filter.cell_name.values

        plt_obj.run(cell_filter.x.values,cell_filter.y.values,\
                    cell_filter.w.values,cell_filter.h.values)

        # plt_obj.drawText(txts)
        
        surface.write_to_png(args.dir+ "imgs/" +"gcell.window.png")  # Output to PNG




    def pltAllGCells(self):
        db = self.db
        plt_obj = PltCairo()
        die_df = db["die"]
        gcell_df = db["gcell"]
        args = db["args"]
        window = [0,0,0,0]
        if die_df is not np.nan:
            window = [
                  die_df["die_xl"].values[0]
                , die_df["die_yl"].values[0]
                , die_df["die_xh"].values[0]
                , die_df["die_yh"].values[0]
            ]
        else:
            print("invalid window box to plot in cell class!")
            sys.exit()
        surface = plt_obj.init(window)
        plt_obj.run(gcell_df.x.values,gcell_df.y.values,\
                    gcell_df.w.values,gcell_df.h.values)
        
        surface.write_to_png(args.dir+ "imgs/" +"gcells.png")  # Output to PNG
        


