

import os
import sys
import pandas as pd
import numpy as np


import_path = os.path.abspath(os.path.join(os.path.join(__file__, ".."), ".."))
sys.path.insert(0, import_path)

from backend.pltcairo import *
from backend.param import *

class Vio:
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
        db = self.db
        die_df = db["die"]
        vio_df = db["vio"]
        args = db["args"]
        
        # surface = plt_obj.init(window)
       
        vio_filter = vio_df.loc[ (vio_df.xl >= window[XL] )& (vio_df.xh <= window[XH])]
        vio_filter = vio_filter.loc[ (vio_df.yl >= window[YL] )& (vio_df.yh <= window[YH])]
        if(l != -1):
            vio_filter = vio_filter.loc[ vio_df.l == l]

        xls = vio_filter.xl.values
        yls = vio_filter.yl.values
        xhs = vio_filter.xh.values
        yhs = vio_filter.yh.values

        ws = [np.abs(xhs[i]-xls[i]) for i in np.arange(len(xls))]
        hs = [np.abs(yhs[i]-yls[i]) for i in np.arange(len(yls))]


        plt_obj.run(vio_filter.xl.values,vio_filter.yl.values,\
                    ws,hs,color,alpha)

        # plt_obj.drawText(txts)
        
        # return surface

    # def pltWindow(self,window):
    #     db = self.db
    #     plt_obj = PltCairo()
    #     die_df = db["die"]
    #     cell_df = db["cell"]
    #     args = db["args"]
        
    #     surface = plt_obj.init(window)

        
    #     cell_filter = cell_df.loc[ (cell_df.x >= window[XL] )& (cell_df.x <= window[XH])]
        
    #     cell_filter = cell_filter.loc[ (cell_df.y >= window[YL] )& (cell_df.y <= window[YH])]
        
    #     txts = cell_filter.cell_name.values

    #     plt_obj.run(cell_filter.x.values,cell_filter.y.values,\
    #                 cell_filter.w.values,cell_filter.h.values)

    #     plt_obj.drawText(txts)
        
    #     surface.write_to_png(args.dir+ "imgs/" +"cells.window.png")  # Output to PNG


    # def pltAllCells(self):
    #     db = self.db
    #     plt_obj = PltCairo()
    #     die_df = db["die"]
    #     cell_df = db["cell"]
    #     args = db["args"]
    #     window = [0,0,0,0]
    #     if die_df is not np.nan:
    #         window = [
    #               die_df["die_xl"].values[0]
    #             , die_df["die_yl"].values[0]
    #             , die_df["die_xh"].values[0]
    #             , die_df["die_yh"].values[0]
    #         ]
    #     else:
    #         print("invalid window box to plot in cell class!")
    #         sys.exit()
    #     surface = plt_obj.init(window)
    #     plt_obj.run(cell_df.x.values,cell_df.y.values,\
    #                 cell_df.w.values,cell_df.h.values)
        
    #     surface.write_to_png(args.dir+ "imgs/" +"cells.all.png")  # Output to PNG
        


