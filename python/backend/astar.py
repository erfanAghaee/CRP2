

# from email.policy import default
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 


import_path = os.path.abspath(os.path.join(os.path.join(__file__, ".."), ".."))
sys.path.insert(0, import_path)

from backend.pltcairo import *
from backend.utils import *
from backend.param import *

class Astar:
    def __init__(self,db,type):
        self.db = db
        self.type = type

    def getWindow(self,window,plt_obj,color,alpha,net_name="-1",l=-1):
        # idx,net_name,l,xl,yl,xh,hy,
        # eBackward,eForward,eDown,eUp,eBackwardCost,eForwardCost,eDownCost,eUpCost
        if not(self.type in self.db):
            return
            
        db = self.db
        die_df = db["die"]
        
        net_df = db[self.type]
        args = db["args"]

        
        net_filter = net_df.loc[net_df.apply(lambda row: getIntervals(row.xl,row.xh,window[XL],window[XH]),axis=1)]
        net_filter = net_filter.loc[net_filter.apply(lambda row: getIntervals(row.yl,row.yh,window[YL],window[YH]),axis=1)]

        # net_filter = net_df.loc[ (net_df.xl >= window[XL] )& (net_df.xh <= window[XH])]
        # net_filter = net_filter.loc[ (net_df.yl >= window[YL] )& (net_df.yh <= window[YH])]        
       
        if(net_name != "-1"):
            net_filter = net_filter.loc[ net_filter.net_name == net_name]
                
        if (l != -1):
            net_filter = net_filter.loc[net_filter.l == l]

        net_filter = net_filter.set_index(net_filter.idx)

        net_filter.index = net_filter.index.astype(int)


        # print(net_filter)
        # print(net_filter)
        xls = []
        yls = []
        xhs = []
        yhs = []
        
        txts=[]

        for i in np.arange(len(net_filter)):
            node_1 = net_filter.iloc[i]

            xl = node_1.xl
            yl = node_1.yl
            xh = node_1.xh
            yh = node_1.yh
            
            if(node_1.eBackward != -1):
                # pass
                node_backward = net_filter.loc[int(node_1.eBackward)]
                xl_backward = node_backward.xl
                yl_backward = node_backward.yl
                xh_backward = node_backward.xh
                yh_backward = node_backward.yh


                xls.append(np.min([xl,xl_backward,xh,xh_backward]))
                yls.append(np.min([yl,yl_backward,yh,yh_backward]))
                xhs.append(np.max([xl,xl_backward,xh,xh_backward]))
                yhs.append(np.max([yl,yl_backward,yh,yh_backward]))
                txts.append("{:.2f}".format(node_1.eBackwardCost))
            
            if(node_1.eForward != -1):
                # pass
                node_forward = net_filter.loc[int(node_1.eForward)]
                xl_forward = node_forward.xl
                yl_forward = node_forward.yl
                xh_forward = node_forward.xh
                yh_forward = node_forward.yh


                xls.append(np.min([xl,xl_forward,xh,xh_forward]))
                yls.append(np.min([yl,yl_forward,yh,yh_forward]))
                xhs.append(np.max([xl,xl_forward,xh,xh_forward]))
                yhs.append(np.max([yl,yl_forward,yh,yh_forward]))

                txts.append("{:.2f}".format(node_1.eForwardCost))


            # xl = 
            # yl =
            # xh = 
            # yh = 
            

        # xls = net_filter.xl.values
        # yls = net_filter.yl.values
        # xhs = net_filter.xh.values
        # yhs = net_filter.yh.values



        w = 0
        if(self.type == "drnet"):
            w = 100
        
        # to show the level of edge in the gridgraph
        # shrink_grid_edge=1500
        # midx = ((xls[i]+xhs[i])/2.0)
        # midy = ((yls[i]+yhs[i])/2.0)

        # (xls[i] + ((xls[i]+xhs[i])/2.0))/2.0
        # xls = [xls[i] + shrink_grid_edge*(abs(xls[i]-xhs[i])) for i in np.arange(len(xls))]
        # yls = [yls[i] + shrink_grid_edge*(abs(yls[i]-yhs[i])) for i in np.arange(len(yls))]
        # xhs = [xhs[i] - shrink_grid_edge*(abs(xls[i]-xhs[i])) for i in np.arange(len(xhs))]
        # yhs = [yhs[i] - shrink_grid_edge*(abs(yls[i]-yhs[i])) for i in np.arange(len(yhs))]
        xls = [(xls[i] + ((xls[i]+xhs[i])/2.0))/2 for i in np.arange(len(xls))]
        yls = [(yls[i] + ((yls[i]+yhs[i])/2.0))/2 for i in np.arange(len(yls))]
        xhs = [(xhs[i] + ((xls[i]+xhs[i])/2.0))/2 for i in np.arange(len(xhs))]
        yhs = [(yhs[i] + ((yls[i]+yhs[i])/2.0))/2  for i in np.arange(len(yhs))]


        ws = [np.abs(xhs[i]-xls[i]+w) for i in np.arange(len(xls))]
        hs = [np.abs(yhs[i]-yls[i]+w) for i in np.arange(len(yls))]


        # txts = net_filter.cost.values

        plt_obj.run(xls,yls,\
                    ws,hs,color,alpha)

        
        plt_obj.drawText(txts,font=72)

        