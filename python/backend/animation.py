#!/usr/bin/env python

import argparse
import os
from pydoc import describe
import sys
import pandas as pd

'''
Python Project Source File imports 
'''

import_path = os.path.abspath(os.path.join(os.path.join(__file__, ".."), ".."))
sys.path.insert(0, import_path)


from backend.cell import *
from backend.net import *
from backend.vio import *
from backend.drc import *
from backend.db import *
from backend.congestionEdge import *
from backend.fixedMetals import *
from backend.utils import *
from backend.gcell import *
from backend.pltcairo import *



def mainAnimation(args):
    specific_net = "net60637"

   

    # iteration of gr
    l = 2
    for i in np.arange(5):

        # if i != 3:
        #     continue

        db = getDB(args,iter_gr=i,iter_dr=0)

        die_df = db["die"]
        net_DRGuide_obj = Net(db,"netDRGuide")
        fixedMetals_obj = FixedMetals(db)
        cell_obj = Cell(db)
        gcell_obj = GCell(db)
        net_obj = Net(db,"net")
        vio_obj = Vio(db)
        drnet_obj = Net(db,"drnet")
        drc_obj = DRC(db)
        congestion_obj = CongestionEdge(db)

        # 81.5
        # window = [315.30,574.2,333.2,589.5]
        # window = [81.5,574.2,333.2,589.5]
        # window = [x*2000 for x in window]
        # window = [642000,1164000,1371000,178000]
        # window = [417,582,445,589]
        window = [417,562,660,589]
        window = [x*2000 for x in window]


        # window = [
        #         die_df["die_xl"].values[0]
        #     , die_df["die_yl"].values[0]
        #     , die_df["die_xh"].values[0]
        #     , die_df["die_yh"].values[0]
        # ]


        plt_obj = PltCairo()
        surface = plt_obj.init(window)

        cell_obj.getWindow(window,plt_obj,(0,0,1),1,text=False)
        gcell_obj.getWindow(window,plt_obj,(1,0,0),0.001)
        net_obj.getWindow(window,plt_obj,(0,1,0),0.8,net_name=specific_net)
        vio_obj.getWindow(window,plt_obj,(0.6,0.1,0.6),0.8,l=l)
        fixedMetals_obj.getWindow(window,plt_obj,(0,0,0),1)
        drc_obj.getWindow(window,plt_obj,(1,0,0),1)
        net_DRGuide_obj.getWindow(window,plt_obj,(0.8,0.1,0.1),0.8,net_name=specific_net)
        drnet_obj.getWindow(window,plt_obj,(0,0,1),1,net_name=specific_net)
        congestion_obj.getWindow(window,plt_obj,(0.5,0.1,0.7),0.1,l=l)

        surface.write_to_png(args.dir+ "imgs/" +"net.gr."+str(i)+".png")
        

        
        



       
        
