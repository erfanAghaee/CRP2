

import os
import sys
import pandas as pd
import numpy as np


import_path = os.path.abspath(os.path.join(os.path.join(__file__, ".."), ".."))
sys.path.insert(0, import_path)

from backend.pltcairo import *
from backend.utils import *
from backend.param import *
import matplotlib.pyplot as plt

class CongestionEdge:
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

    def getWindow(self,window,plt_obj,l=-1,text=False):
        if not("congestion" in self.db):
            return

        db = self.db
        die_df = db["die"]
        congestion_df = db["congestion"]
        args = db["args"]
        
        # surface = plt_obj.init(window)

        
       
        congestion_filter = congestion_df.loc[congestion_df.apply(lambda row: getIntervals(row.xl,row.xh,window[XL],window[XH]),axis=1)]
        congestion_filter = congestion_filter.loc[congestion_filter.apply(lambda row: getIntervals(row.yl,row.yh,window[YL],window[YH]),axis=1)]

        # congestion_filter = congestion_df.loc[ (congestion_df.xl >= window[XL] )& (congestion_df.xh <= window[XH])]
        # congestion_filter = congestion_filter.loc[ (congestion_df.yl >= window[YL] )& (congestion_df.yh <= window[Y

        # i1 = pd.Interval(0, 2)
        # i2 = pd.Interval(1, 3)
        # print(i1.overlaps(i2))


        if(l != -1):
            congestion_filter = congestion_filter.loc[ congestion_df.l == l]


        txts_wireUsage = congestion_filter.wireUsage.values
        txts_fixedUsage = congestion_filter.fixedUsage.values
        txts_viaUsage = congestion_filter.viaUsage.values
        txts_numTracks = congestion_filter.numTracks.values
        txts = []
        
        if(text):
            for i in np.arange(len(txts_wireUsage)):
                # txt = "w:{},f:{}\nv:{},t:{}\nr:{}".format(txts_wireUsage[i],\
                #     txts_fixedUsage[i],\
                #     txts_viaUsage[i],\
                #     txts_numTracks[i])
                # txt = "{:.2f}".format(txts_wireUsage[i]+txts_fixedUsage[i]+txts_viaUsage[i]-txts_numTracks[i])
                # txt = "{:.1f}".format(txts_wireUsage[i])
                txt = "{:.1f}".format(txts_fixedUsage[i])
                
                # txt2 = "{:.1f}".format(txts_viaUsage[i])
                # txt3 = "{:.1f}".format(txts_numTracks[i])
                # txt = txt3
                # txt = txt0 + "|"+txt1 + "|"+txt2 + "|"+txt3 
                # # txt =txt1 +"|"+txt3 
                # print("txt0:",txt0)
                    
                txts.append(txt)
        
        # txts = [str(s) for s in np.arange(len(txts_wireUsage))]

        # txts = [txts[i]+"/"+str(txts_viaUsage[i]) \
        #     for i in np.arange(len(txts_viaUsage))]

        xls = congestion_filter.xl.values
        yls = congestion_filter.yl.values
        xhs = congestion_filter.xh.values
        yhs = congestion_filter.yh.values

        ws = [np.abs(xhs[i]-xls[i]) for i in np.arange(len(xls))]
        hs = [np.abs(yhs[i]-yls[i]) for i in np.arange(len(yls))]

        colors =[(i,0,0) for i in congestion_filter.wireUsage.values]
        alphas =[0.1 for i in range(len(xls))]


        plt_obj.run(xls,yls,\
                    ws,hs,colors,alphas)

        if(text):
            plt_obj.drawText(txts,color=(1,1,1),font=10)
        
    def getStatistics(self,window):
        if not("congestion" in self.db):
            return

        db = self.db
        die_df = db["die"]
        congestion_df = db["congestion"]
        args = db["args"]
        
        # print(congestion_df)
        grps = congestion_df.groupby(['l'])
        avgs = grps["viaUsage"].mean()
        norms = [float(i)/max(avgs) for i in avgs.values]
        print(avgs)

        layers = grps.groups.keys()

        plt.plot(layers,norms,label="congestion.gr")
        
        
        
 
        