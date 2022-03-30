
import os
import sys 
import pandas as pd
import numpy as np


import_path = os.path.abspath(os.path.join(os.path.join(__file__, ".."), ".."))
sys.path.insert(0, import_path)

from backend.pltcairo import *

class Cell:
    def __init__(self,db):
        self.db = db


    def run(self):
        # db =  self.db
        # print(cells_df)
        self.pltAllCells()

    def pltAllCells(self):
        db = self.db
        plt_obj = PltCairo()
        die_df = db["die"]
        cell_df = db["cell"]
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
        plt_obj.run(cell_df.x.values,cell_df.y.values,\
                    cell_df.w.values,cell_df.h.values)
        
        surface.write_to_png(args.dir+ "imgs/" +"cells.png")  # Output to PNG
        pass


