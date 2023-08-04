import os
import sys
import pandas as pd
import numpy as np


import_path = os.path.abspath(os.path.join(os.path.join(__file__, ".."), ".."))
sys.path.insert(0, import_path)


from backend.param import *

def getIntervals(l,h,wl,wh):
    return pd.Interval(l,h).overlaps(pd.Interval(wl, wh))


def getArea(row):
    return abs(row.xh-row.xl)+abs(row.yh-row.yl)

def compareTwoRoute(route_1,route_2):
    print(len(route_1))
    print(len(route_2))
    r1_gp = route_1.groupby("net_name")
    r2_gp = route_2.groupby("net_name")
    print(len(r1_gp))
    print(len(r2_gp))
    nets = set()
    for grp in r1_gp.groups:
        r1_df = r1_gp.get_group(grp)
        r2_df = r2_gp.get_group(grp)

        r1_df["area"] = r1_df.apply(lambda row: getArea(row),axis=1)
        r2_df["area"] = r2_df.apply(lambda row: getArea(row),axis=1)

        print(r1_df)
        # if(len(r1_df) != len(r2_df)):
        if(r1_df["area"].sum() != \
           r2_df["area"].sum()):
            nets.add(grp)
            print("not equal")


        break
        

        # if(r1_df[["l","xl","yl","xh","yh"]].equals(r2_df[["l","xl","yl","xh","yh"]])):
        #     pass
        # else:
        #     nets.add(grp)
        
    print(len(nets))
    i = 0
    for net in nets:
        print(net)
        if i == 2:
            break
        i += 1
    return nets

# can be used for get out of guide
def getBoxDifference(A,B):
    
    if ((A[XH] <= B[XL]) or (A[XL] >= B[XH]) or (A[YH] <= B[YL]) or (A[YL] >= B[YH])): 
        # print("no intersection")
        return [[A[XL],A[YL],A[XH],A[YH]]] # no intersection
    else:
        if((A[XL] >= B[XL]) and (A[YL] >= B[YL]) and (A[XH] <= B[XH]) and (A[YH] <= B[YH]) ):
            return [[]] # within and detailed router is inside the box

        # 4-direction
        elif ((A[XL] < B[XL]) and (A[YL] >= B[YL]) and (A[XH] <= B[XH]) and (A[YH] <= B[YH])): # left 
            return  [[A[XL],A[YL],B[XL],A[YH]]]
        elif ((A[XL] >= B[XL]) and (A[YL] >= B[YL]) and (A[XH] > B[XH]) and (A[YH] <= B[YH])): # right
            return  [[B[XH],A[YL],A[XH],A[YH]]]
        elif ((A[XL] >= B[XL]) and (A[YL] >= B[YL]) and (A[XH] <= B[XH]) and (A[YH] > B[YH])): # top
            return  [[A[XL],B[YH],A[XH],A[YH]]]
        elif ((A[XL] >= B[XL]) and (A[YL] < B[YL]) and (A[XH] <= B[XH]) and (A[YH] <= B[YH])): # bottom 
            return  [[A[XL],A[YL],A[XH],B[YL]]]

        # 2-cross
        elif ((A[XL] < B[XL]) and (A[YL] >= B[YL]) and (A[XH] > B[XH]) and (A[YH] <= B[YH])): # left-right
            return  [[A[XL],A[YL],B[XL],A[YH]],[B[XH],A[YL],A[XH],A[YH]]]
        elif ((A[XL] >= B[XL]) and (A[YL] < B[YL]) and (A[XH] <= B[XH]) and (A[YH] > B[YH])): # top-bottom
            return  [[A[XL],A[YL],A[XH],B[YL]],[A[XL],B[YH],A[XH],A[YH]]]


        # 4-corners
        elif ((A[XL] < B[XL]) and (A[YL] < B[YL]) and (A[XH] <= B[XH]) and (A[YH] <= B[YH])): # bottom-left
            return  [[A[XL],A[YL],B[XL],A[YH]],[B[XL],A[YL],A[XH],B[YL]]]
        elif ((A[XL] >= B[XL]) and (A[YL] < B[YL]) and (A[XH] > B[XH]) and (A[YH] <= B[YH])): # bottom-right
            return  [[B[XH],A[YL],A[XH],A[YH]],[A[XL],A[YL],B[XH],B[YL]]]
        elif ((A[XL] < B[XL]) and (A[YL] >= B[YL]) and (A[XH] <= B[XH]) and (A[YH] > B[YH])): # top-left
            return  [[A[XL],A[YL],B[XL],A[YH]],[B[XL],B[YH],A[XH],A[YH]]]
        elif ((A[XL] >= B[XL]) and (A[YL] >= B[YL]) and (A[XH] > B[XH]) and (A[YH] > B[YH])): # top-right
            return  [[B[XH],A[YL],A[XH],A[YH]],[A[XL],B[YH],B[XH],A[YH]]]


        # 4-corner crossing
        elif ((A[XL] < B[XL]) and (A[YL] >= B[YL]) and (A[XH] > B[XH]) and (A[YH] > B[YH])): # top-left-right-corner
            return  [[A[XL],A[YL],B[XL],A[YH]],[B[XL],B[YH],B[XH],A[YH]],[B[XH],A[YL],A[XH],A[YH]]]
        elif ((A[XL] < B[XL]) and (A[YL] < B[YL]) and (A[XH] > B[XH]) and (A[YH] <= B[YH])): # bottom-left-right-corner
            return  [[A[XL],A[YL],B[XL],A[YH]],[B[XL],A[YL],B[XH],B[YL]],[B[XH],A[YL],A[XH],A[YH]]]
        elif ((A[XL] < B[XL]) and (A[YL] < B[YL]) and (A[XH] <= B[XH]) and (A[YH] > B[YH])): # left-top-bottom-corner
            return  [[A[XL],A[YL],A[XH],B[YL]],[A[XL],B[YL],B[XL],B[YH]],[A[XL],B[YH],A[XH],A[YH]]]
        elif ((A[XL] >= B[XL]) and (A[YL] < B[YL]) and (A[XH] > B[XH]) and (A[YH] > B[YH])): # right-top-bottom-corner
            return  [[A[XL],A[YL],A[XH],B[YL]],[B[XH],B[YL],A[XH],B[YH]],[A[XL],B[YH],A[XH],A[YH]]]


        # wrapped
        elif ((A[XL] < B[XL]) and (A[YL] < B[YL]) and (A[XH] > B[XH]) and (A[YH] > B[YH])): # wrap
            return  [[A[XL],A[YL],A[XH],B[YL]],[B[XH],B[YL],A[XH],B[YH]],[A[XL],B[YH],A[XH],A[YH]],[A[XL],B[YL],B[XL],B[YH]]]


def getBoxIntersection(A,B):
    
    if ((A[XH] <= B[XL]) or (A[XL] >= B[XH]) or (A[YH] <= B[YL]) or (A[YL] >= B[YH])): 
        # print("no intersection")
        return [[]] # no intersection
    else:
        if((A[XL] >= B[XL]) and (A[YL] >= B[YL]) and (A[XH] <= B[XH]) and (A[YH] <= B[YH]) ):
            return [[A[XL],A[YL],A[XH],A[YH]]] # within and detailed router is inside the box

        # 4-direction
        elif ((A[XL] < B[XL]) and (A[YL] >= B[YL]) and (A[XH] <= B[XH]) and (A[YH] <= B[YH])): # left 
            return  [[B[XL],A[YL],A[XH],A[YH]]]
        elif ((A[XL] >= B[XL]) and (A[YL] >= B[YL]) and (A[XH] > B[XH]) and (A[YH] <= B[YH])): # right
            return  [[A[XL],A[YL],B[XH],A[YH]]]
        elif ((A[XL] >= B[XL]) and (A[YL] >= B[YL]) and (A[XH] <= B[XH]) and (A[YH] > B[YH])): # top
            return  [[A[XL],A[YL],A[XH],B[YH]]]
        elif ((A[XL] >= B[XL]) and (A[YL] < B[YL]) and (A[XH] <= B[XH]) and (A[YH] <= B[YH])): # bottom 
            return  [[A[XL],B[YL],A[XH],A[YH]]]

        # 2-cross
        elif ((A[XL] < B[XL]) and (A[YL] >= B[YL]) and (A[XH] > B[XH]) and (A[YH] <= B[YH])): # left-right
            return  [[B[XL],A[YL],B[XH],A[YH]]]
        elif ((A[XL] >= B[XL]) and (A[YL] < B[YL]) and (A[XH] <= B[XH]) and (A[YH] > B[YH])): # top-bottom
            return  [[A[XL],B[YL],A[XH],B[YH]]]


        # 4-corners
        elif ((A[XL] < B[XL]) and (A[YL] < B[YL]) and (A[XH] <= B[XH]) and (A[YH] <= B[YH])): # bottom-left
            return  [[B[XL],B[YL],A[XH],A[YH]]]
        elif ((A[XL] >= B[XL]) and (A[YL] < B[YL]) and (A[XH] > B[XH]) and (A[YH] <= B[YH])): # bottom-right
            return  [[A[XL],B[YL],B[XH],A[YH]]]
        elif ((A[XL] < B[XL]) and (A[YL] >= B[YL]) and (A[XH] <= B[XH]) and (A[YH] > B[YH])): # top-left
            return  [[B[XL],A[YL],A[XH],B[YH]]]
        elif ((A[XL] >= B[XL]) and (A[YL] >= B[YL]) and (A[XH] > B[XH]) and (A[YH] > B[YH])): # top-right
            return  [[A[XL],A[YL],B[XH],B[YH]]]


        # 4-corner crossing
        elif ((A[XL] < B[XL]) and (A[YL] >= B[YL]) and (A[XH] > B[XH]) and (A[YH] > B[YH])): # top-left-right-corner
            return  [[B[XL],A[YL],B[XH],B[YH]]]
        elif ((A[XL] < B[XL]) and (A[YL] < B[YL]) and (A[XH] > B[XH]) and (A[YH] <= B[YH])): # bottom-left-right-corner
            return  [[B[XL],B[YL],B[XH],A[YH]]]
        elif ((A[XL] < B[XL]) and (A[YL] < B[YL]) and (A[XH] <= B[XH]) and (A[YH] > B[YH])): # left-top-bottom-corner
            return  [[B[XL],B[YL],A[XH],B[YH]]]
        elif ((A[XL] >= B[XL]) and (A[YL] < B[YL]) and (A[XH] > B[XH]) and (A[YH] > B[YH])): # right-top-bottom-corner
            return  [[A[XL],B[YL],B[XH],B[YH]]]


        # wrapped
        elif ((A[XL] < B[XL]) and (A[YL] < B[YL]) and (A[XH] > B[XH]) and (A[YH] > B[YH])): # wrap
            return  [[B[XL],B[YL],B[XH],B[YH]]]


def getWirelength(box):
    return abs(box[XL]-box[XH])+abs(box[YL]-box[YH])

def getDirection(box):
    return (H if abs(box[XL]-box[XH]) > abs(box[YL]-box[YH]) else V)

def getWirelengthWWW(box,dir):
    if getDirection(box) != dir:
        return getWirelength(box)
    return 0

def getOffTrack(box,dir,patternTrack,l):
    debug = False
    dir = getDirection(box)
    if(debug):
        print("dir:",dir)

    if dir == V:
        center = (box[XL] + box[XH])/2.0
    else:
        center = (box[YL] + box[YH])/2.0

    if(debug):
        print(center)
    
    if dir == V:
        patternTrack = patternTrack.loc[patternTrack.dir  == 'X']
    else:
        patternTrack = patternTrack.loc[patternTrack.dir  == 'Y']

    if(debug):
        print(patternTrack)

    start = patternTrack.start.values[0]
    step = patternTrack.step.values[0]

    if(debug):
        print("start: ",start)
        print("step: ",step)

    a = float((center - start))/float(step)
    new_center = int(a)*step + start

    if(debug):
        print("old center: ",center,", new_center: ",new_center)

    if center != new_center:
        if(debug):
            print("old center: ",center,", new_center: ",new_center,", box: ",box, ", l: ",l)
        # # print(new_center,center)
        # print("wl:", getWirelength(box))
        return getWirelength(box)
    return 0


def getArea(box,width):
    x = abs(box[XH]-box[XL]) 
    y = abs(box[YH]-box[YL]) 
    # if( x == 0):
    x += (width*2000)
    # if( y == 0):
    y += (width*2000)

    return x/2000.0*y/2000.0

