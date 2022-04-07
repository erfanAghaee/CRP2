

import os
import sys
import pandas as pd
import numpy as np


import_path = os.path.abspath(os.path.join(os.path.join(__file__, ".."), ".."))
sys.path.insert(0, import_path)

from backend.pltcairo import *
from backend.param import *

class Cell:
    def __init__(self,db):
        self.db = db


    def getWindow(self,window,plt_obj,color,alpha,text=False):
        if not("cell" in self.db):
            return

        db = self.db
        die_df = db["die"]
        cell_df = db["cell"]
        args = db["args"]
        
        # surface = plt_obj.init(window)

        
        cell_filter = cell_df.loc[ (cell_df.x >= window[XL] )& (cell_df.x <= window[XH])]
        
        cell_filter = cell_filter.loc[ (cell_df.y >= window[YL] )& (cell_df.y <= window[YH])]
        
        txts = cell_filter.cell_name.values

        plt_obj.run(cell_filter.x.values,cell_filter.y.values,\
                    cell_filter.w.values,cell_filter.h.values,color,alpha)

        if(text):
            plt_obj.drawText(txts)
        
        # return surface

    