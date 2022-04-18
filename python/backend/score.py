#!/usr/bin/env python

import argparse
# from ast import pattern
import os
# from pydoc import describe
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
from backend.coef import *
from backend.astar import *
from backend.patternroute import *
from backend.congestionEdge import *
from backend.fixedMetals import *
from backend.utils import *
from backend.gcell import *
from backend.pltcairo import *
# from backend.outofguide import *
import matplotlib.pyplot as plt
from pandarallel import pandarallel
from joblib import Parallel, delayed
import multiprocessing

class Score:
    def __init__(self,db,type):
        self.db = db
        self.type = type

    def getOutOfGuide(self,dr,gr):
        dr_res = copy.deepcopy(dr)
        for i in range(len(gr)):
            dr_queue = []
            
            for j in range(len(dr_res)):
                dr_queue.append(dr_res[j])
            dr_cropped = []
            
            for j in range(len(dr_queue)):
                # dr_tmp = getBoxDifference(dr_res[j],gr[i])
                
                dr_tmp = getBoxDifference(dr_queue[j],gr[i])
                
                for k in range(len(dr_tmp)):
                    if(len(dr_tmp[k]) > 0):
                        dr_cropped.append(dr_tmp[k])
            dr_res = dr_cropped

        return dr_res


    def processOutOfGuide(self,i):
        grp = self.grp_dict[i]
        grpby_drnet = self.grpby_drnets.get_group(grp)
        grpby_grnet = self.grpby_grnets.get_group(grp)

        outofguides = []
        for l in range(9):
            grpby_drnets_layer = grpby_drnet.loc[grpby_drnet.l == l]
            grpby_grnets_layer = grpby_grnet.loc[grpby_grnet.l == l]

            dr_xls = grpby_drnets_layer.xl.values
            dr_yls = grpby_drnets_layer.yl.values
            dr_xhs = grpby_drnets_layer.xh.values
            dr_yhs = grpby_drnets_layer.yh.values

            gr_xls = grpby_grnets_layer.xl.values
            gr_yls = grpby_grnets_layer.yl.values
            gr_xhs = grpby_grnets_layer.xh.values
            gr_yhs = grpby_grnets_layer.yh.values
            dr = [[dr_xls[i],dr_yls[i],dr_xhs[i],dr_yhs[i]] for i in range(len(dr_xls))]
            gr = [[gr_xls[i],gr_yls[i],gr_xhs[i],gr_yhs[i]] for i in range(len(gr_xls))]

            outofguide = self.getOutOfGuide(dr=dr,gr=gr)
            # print(outofguide)
            for box in outofguide:
                if(len(box) > 0):
                    outofguides.append(box)

            # break
        return outofguides
        


    def getOutOfGuideTotal(self,type="wire"):
        pandarallel.initialize()

        db = self.db
        
        drnet_df = db["drnet"]
        grnet_df = db["net"]

        # drnet_df = drnet_df.loc[drnet_df.net_name=="net446"]
        # grnet_df = grnet_df.loc[grnet_df.net_name=="net446"]
        if(type=="wire"):
            drnet_df = drnet_df.loc[(drnet_df.type == "wire") | (drnet_df.type == "patch") ]
        elif(type=="via"):
            drnet_df = drnet_df.loc[(drnet_df.type == "via") ]
        # grnet_df = grnet_df.loc[grnet_df.type == "wire"]

        
        outofguides = []
        self.grpby_drnets = drnet_df.groupby(["net_name"])
        self.grpby_grnets = grnet_df.groupby(["net_name"])

        i = 0
        self.grp_dict = {}
        nets = self.grpby_drnets.groups
        for grp in self.grpby_drnets.groups:
            self.grp_dict[i] = grp
            i +=1

        # for grp in grpby_drnets.groups:
        # Parallel(n_jobs=-1)(delayed(self.process)(i) for i in range(len(self.grpby_drnets.groups)))

        pool_obj = multiprocessing.Pool()
        outofguides = pool_obj.map(self.processOutOfGuide,range(0,len(self.grpby_drnets.groups)))
        # outofguides = pool_obj.map(self.process,range(0,2))
        # print(outofguides)

        outofguidesTotal = []
        for guide_list in outofguides:
            for guide in guide_list:
                if(len(guide) > 0):
                    outofguidesTotal.append(guide)

        # print(outofguidesTotal)
        wls = [getWirelength(box) for box in outofguidesTotal]
        
        
        if type=="via":
            return len(outofguidesTotal)

        return np.sum(wls)


    def getWrongWayWiring(self):
        pandarallel.initialize()
        # layers direction
        layers_dir=[V,H,V,H,V,H,V,H,V]
        db = self.db
        die_df = db["die"]
        
        # net_df = db[self.type]
        drnet_df = db["drnet"]  

        total_www = 0
        for l in range(9):
            drnet_filter = drnet_df.loc[drnet_df.l==l]
            drnet_filter = drnet_filter.loc[drnet_filter.type=="wire"]

            res = drnet_filter.parallel_apply(lambda row: getWirelengthWWW([row["xl"],row["yl"],row["xh"],row["yh"]],layers_dir[l]) ,axis=1 )    
            total_www += np.sum(res)

        return total_www


    def getWirelengthViasTotal(self,net_name="-1",l=-1):
        pandarallel.initialize()

        db = self.db
        drnet_df = db["drnet"]  
        drnet_filter = drnet_df.loc[drnet_df.type=="wire"]
        # res = drnet_filter.apply(lambda row: abs(row["xl"]-row["xh"])+abs(row["yl"]-row["yh"]),axis=1 )
        res = drnet_filter.parallel_apply(lambda row: getWirelength([row["xl"],row["yl"],row["xh"],row["yh"]]) ,axis=1 )
        wl = np.sum(res)
        drnet_filter = drnet_df.loc[drnet_df.type=="via"]        
        num_vias = len(drnet_filter)/2

        return wl,num_vias


    def getOffTrackTotal(self):
        #  start,  numTracks, steps
        # TRACKS Y 200 DO 2946 STEP 400 LAYER Metal9 ;
        # TRACKS X 500 DO 4363 STEP 400 LAYER Metal9 ;
        # TRACKS X 500 DO 4363 STEP 400 LAYER Metal8 ;
        # TRACKS Y 200 DO 2946 STEP 400 LAYER Metal8 ;
        # TRACKS Y 200 DO 3928 STEP 300 LAYER Metal7 ;
        # TRACKS X 400 DO 5818 STEP 300 LAYER Metal7 ;
        # TRACKS X 400 DO 5818 STEP 300 LAYER Metal6 ;
        # TRACKS Y 200 DO 3928 STEP 300 LAYER Metal6 ;
        # TRACKS Y 200 DO 5891 STEP 200 LAYER Metal5 ;
        # TRACKS X 100 DO 8728 STEP 200 LAYER Metal5 ;
        # TRACKS X 100 DO 8728 STEP 200 LAYER Metal4 ;
        # TRACKS Y 200 DO 5891 STEP 200 LAYER Metal4 ;
        # TRACKS Y 200 DO 5891 STEP 200 LAYER Metal3 ;
        # TRACKS X 100 DO 8728 STEP 200 LAYER Metal3 ;
        # TRACKS X 100 DO 8728 STEP 200 LAYER Metal2 ;
        # TRACKS Y 200 DO 5891 STEP 200 LAYER Metal2 ;
        # TRACKS Y 200 DO 5891 STEP 200 LAYER Metal1 ;
        # TRACKS X 100 DO 8728 STEP 200 LAYER Metal1 ;
        
        # csv
        data = {
            'dir'  : ['X', 'Y', 'X', 'Y','X', 'Y', 'X', 'Y','X', 'Y', 'X', 'Y','X', 'Y', 'X', 'Y','X','Y'], 
            'start': [100,200,200,100,100,200,200,100,100,200,200,400,400,200,200,500,500,200],
            'numTracks': [8728,5891,5891,8728,8728,5891,5891,8728,8728,5891,3928,5818,5818,3928,2946,4363,4363,2946],
            'step':[200,200,200,200,200,200,200,200,200,200,300,300,300,300,400,400,400,400],
            'l':[0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8]
        }  
            
        track_pattern_df = pd.DataFrame(data)
        pandarallel.initialize()
        # layers direction
        layers_dir=[V,H,V,H,V,H,V,H,V]
        db = self.db
        die_df = db["die"]
        
        # net_df = db[self.type]
        drnet_df = db["drnet"]  

        total_offtrack = 0
        for l in range(9):
            # if l != 0:
            #     continue

            track_pattern_filter = track_pattern_df.loc[track_pattern_df.l == l]

            drnet_filter = drnet_df.loc[drnet_df.l==l]
            drnet_filter = drnet_filter.loc[drnet_filter.type=="wire"]
            # drnet_filter = drnet_filter.loc[drnet_filter.net_name=="net446"]
            # print(drnet_filter)
            # res = drnet_filter.parallel_apply(lambda row: \
            #     getOffTrack([row["xl"],row["yl"],row["xh"],row["yh"]],layers_dir[l],track_pattern_filter) ,axis=1 )    
            res = drnet_filter.parallel_apply(lambda row: \
                getOffTrack([row["xl"],row["yl"],row["xh"],row["yh"]],layers_dir[l],track_pattern_filter) ,axis=1 )    
            # print("l:",l,res)
            total_offtrack += np.sum(res)
        # l = 0
        # track_pattern_filter = track_pattern_df.loc[track_pattern_df.l == l]
        # getOffTrack([381900.0, 861400.0, 381900.0, 862400.0],l,track_pattern_filter)
        return total_offtrack
        
        
        



        

        
        
        
        
        
        
        
        
        
        
        
            
            
            
            
            
            
        track_pattern_df = pd.DataFrame(data) 


    def getWindow(self,window,plt_obj,color,alpha,net_name="-1",l=-1):
        # if not(self.type in self.db):
        #     return
            
        db = self.db
        die_df = db["die"]
        
        # net_df = db[self.type]
        drnet_df = db["drnet"]
        grnet_df = db["net"]
        args = db["args"]

        # net_filter = net_df.loc[net_df.apply(lambda row: getIntervals(row.xl,row.xh,window[XL],window[XH]),axis=1)]
        # net_filter = net_filter.loc[net_filter.apply(lambda row: getIntervals(row.yl,row.yh,window[YL],window[YH]),axis=1)]
        

        drnet_filter = drnet_df.loc[ (drnet_df.xl >= window[XL] )& (drnet_df.xh <= window[XH])]
        drnet_filter = drnet_filter.loc[ (drnet_filter.yl >= window[YL] )& (drnet_filter.yh <= window[YH])]

        grnet_filter = grnet_df.loc[ (grnet_df.xl >= window[XL] )& (grnet_df.xh <= window[XH])]
        grnet_filter = grnet_filter.loc[ (grnet_filter.yl >= window[YL] )& (grnet_filter.yh <= window[YH])]
       
        if(net_name != "-1"):
            drnet_filter = drnet_filter.loc[ drnet_filter.net_name == net_name]
            grnet_filter = grnet_filter.loc[ grnet_filter.net_name == net_name]
        # only second layer
        if (l != -1):
            drnet_filter = drnet_filter.loc[drnet_filter.l == l]
            grnet_filter = grnet_filter.loc[grnet_filter.l == l]

        
        dr_xls = drnet_filter.xl.values
        dr_yls = drnet_filter.yl.values
        dr_xhs = drnet_filter.xh.values
        dr_yhs = drnet_filter.yh.values

        gr_xls = grnet_filter.xl.values
        gr_yls = grnet_filter.yl.values
        gr_xhs = grnet_filter.xh.values
        gr_yhs = grnet_filter.yh.values

        dr = [[dr_xls[i],dr_yls[i],dr_xhs[i],dr_yhs[i]] for i in range(len(dr_xls))]
        gr = [[gr_xls[i],gr_yls[i],gr_xhs[i],gr_yhs[i]] for i in range(len(gr_xls))]

        outofguide = self.getOutOfGuide(dr=dr,gr=gr)
        print("l: ",l,", outofguide: ",outofguide)

        if not(len(outofguide) > 0):
           return 
            
        xls = [x[XL] for x in outofguide ]
        yls = [x[YL] for x in outofguide ]
        xhs = [x[XH] for x in outofguide ]
        yhs = [x[YH] for x in outofguide ]

        # w = 0
        # if(self.type == "drnet"):
        w = 100
        

        ws = [np.abs(xhs[i]-xls[i]+w) for i in np.arange(len(xls))]
        hs = [np.abs(yhs[i]-yls[i]+w) for i in np.arange(len(yls))]


        plt_obj.run(xls,yls,\
                    ws,hs,color,alpha)




            



            


            


    