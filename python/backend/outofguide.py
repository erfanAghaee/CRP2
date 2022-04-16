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


from backend.pltcairo import *
from backend.param import *
from rtree import index 


def plot(plt_obj,boxs,color,alpha,type):
    xls = [box[XL] for box in boxs]
    yls = [box[YL] for box in boxs]
    xhs = [box[XH] for box in boxs]
    yhs = [box[YH] for box in boxs]

    ws = [np.abs(xhs[i]-xls[i]) for i in np.arange(len(xls))]
    hs = [np.abs(yhs[i]-yls[i]) for i in np.arange(len(yls))]


    plt_obj.run(xls,yls,ws,hs,color,alpha)

def initRtree(boxs,idx):

    i = 0
    for box in boxs:
        idx.insert(i, (box[XL], box[YL], box[XH], box[YH]))
        i +=1
def getUnionBox(boxs):
    xls = [min(box[XL],box[XH]) for box in boxs]
    yls = [min(box[YL],box[YL]) for box in boxs]
    xhs = [max(box[XH],box[XH]) for box in boxs]
    yhs = [max(box[YH],box[YH]) for box in boxs]
    return [[min(xls),min(yls),max(xhs),max(yhs)]]


# def contain(rectA,rectB):
#     wA = abs(rectA[XL]-rectA[XH])
#     hA = abs(rectA[YL]-rectA[YH])
#     wB = abs(rectB[XL]-rectB[XH])
#     hB = abs(rectB[YL]-rectB[YH])
#     return (rectB[XL] >= rectA[XL]) and 
#            (rectB[YL] >= rectA[YL]) and
#            (rectB[XL] + wB <= rectA[XL] + wA) and
#            (rectB[YL] + hB <= rectA[YL] + hA)

# # https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Rectangle_difference
# def intersect(rectA,rectB):
#     wA = abs(rectA[XL]-rectA[XH])
#     hA = abs(rectA[YL]-rectA[YH])
#     wB = abs(rectB[XL]-rectB[XH])
#     hB = abs(rectB[YL]-rectB[YH])
#     return not((rectB[XL] + wB <= rectA[XL]) or
#              (rectB[YL] + hB <= rectA[YL]) or
#              (rectB[XL] >= rectA[XL] + wA) or
#              (rectB[YL] >= rectA[YL] + hA))



def checkoutofguide(dr,gr):
    idx_gr = index.Index()
    # idx_dr = index.Index()

    # initRtree(dr,idx_dr)
    initRtree(gr,idx_gr)

    for dr_box in dr:
        ovrlp_boxs = [gr[i] for i in list(idx_gr.intersection((dr_box[XL],dr_box[YL],dr_box[XH],dr_box[YH])))]

        unionbox = getUnionBox(ovrlp_boxs)

        for box in ovrlp_boxs:
            print(box,dr)
            print("intersected: ",intersect(dr[0],box))
        
        print(ovrlp_boxs)

def initGCells(window):
    gcells = []
    for i in np.arange(window[XH]):
        for j in np.arange(window[YH]):
            # print([i,j,i+0.91,j+1])
            gcells.append([i,j,i+0.9,j+0.9])
    return gcells

def bitWiseCompare(A,B):
    if (A < B):
        return 1
    return 0




def getBoxDifference(A,B):
    
    if ((A[XH] <= B[XL]) or (A[XL] >= B[XH]) or (A[YH] <= B[YL]) or (A[YL] >= B[YH])): 
        print("no intersection")
        return [A[XL],A[YL],A[XH],A[YH]] # no intersection
    else:
        if((A[XL] >= B[XL]) and (A[YL] >= B[YL]) and (A[XH] <= B[XH]) and (A[YH] <= B[YH]) ):
            print("within")
            return [] # within and detailed router is inside the box

        # 4-direction
        elif ((A[XL] < B[XL]) and (A[YL] >= B[YL]) and (A[XH] <= B[XH]) and (A[YH] <= B[YH])): # left 
            return  [A[XL],A[YL],B[XL],A[YH]]
        elif ((A[XL] >= B[XL]) and (A[YL] >= B[YL]) and (A[XH] > B[XH]) and (A[YH] <= B[YH])): # right
            return  [B[XH],A[YL],A[XH],A[YH]]
        elif ((A[XL] >= B[XL]) and (A[YL] >= B[YL]) and (A[XH] <= B[XH]) and (A[YH] > B[YH])): # top
            return  [A[XL],B[YH],A[XH],A[YH]]
        elif ((A[XL] >= B[XL]) and (A[YL] < B[YL]) and (A[XH] <= B[XH]) and (A[YH] <= B[YH])): # bottom 
            return  [A[XL],A[YL],A[XH],B[YL]]

        # 2-cross
        elif ((A[XL] < B[XL]) and (A[YL] >= B[YL]) and (A[XH] > B[XH]) and (A[YH] <= B[YH])): # left-right
            return  [[A[XL],A[YL],B[XL],A[YH]],[B[XH],A[YL],A[XH],A[YH]]]
        elif ((A[XL] >= B[XL]) and (A[YL] < B[YL]) and (A[XH] <= B[XH]) and (A[YH] > B[YH])): # top-bottom
            return  [[A[XL],A[YL],A[XH],B[YL]],[A[XL],B[YH],A[XH],A[YH]]]


    #     # 4-corners
    #     elif ((A[XL] < B[XL]) and (A[YL] < B[YL]) and (A[XH] <= B[XH]) and (A[YH] <= B[YH])): # bottom-left
    #         return  [[A[XL],A[YL],B[XL],A[YH]],[B[XL],A[YL],A[XH],B[YL]]]
    #     elif ((A[XL] >= B[XL]) and (A[YL] < B[YL]) and (A[XH] > B[XH]) and (A[YH] <= B[YH])): # bottom-right
    #         return  [[B[XH],A[YL],A[XH],A[YH]],[A[XL],A[YL],B[XH],B[YL]]]
    #     elif ((A[XL] < B[XL]) and (A[YL] >= B[YL]) and (A[XH] <= B[XH]) and (A[YH] > B[YH])): # top-left
    #         return  [[A[XL],A[YL],B[XL],A[YH]],[B[XL],B[YH],A[XH],A[YH]]]
    #     elif ((A[XL] >= B[XL]) and (A[YL] >= B[YL]) and (A[XH] > B[XH]) and (A[YH] > B[YH])): # top-right
    #         return  [[B[XH],A[YL],A[XH],A[YH]],[A[XL],B[YH],B[XH],A[YH]]]


def draw(args,dr,gr,name):
    window = [
        0,0,20,20
    ]

    gcells = initGCells(window)


    plt_obj = PltCairo()
    surface = plt_obj.init(window)

    
    plot(plt_obj,dr,(0,1,0),1,"dr")
    plot(plt_obj,gr,(0,0,1),0.2,"gr")
    plot(plt_obj,gcells,(1,0,0),0.1,"gcells")

    # checkoutofguide(dr,gr)

    surface.write_to_png(args.dir+ name +".png")

def testNointersection(args):
    dr = [[1,5,4,6]]
    gr = [[5,1,15,10]]

    print("noIntersection:",getBoxDifference(dr[0],gr[0]))
    
    draw(args,dr,gr,"noIntersection")

def testWithin(args):
    dr = [[5,5,10,6]]
    gr = [[5,1,15,10]]

    print("within: ",getBoxDifference(dr[0],gr[0]))
    
    draw(args,dr,gr,"within")
   
def testLeft(args):
    dr = [[3,5,10,6]]
    gr = [[5,1,15,10]]

    print("left: ",getBoxDifference(dr[0],gr[0]))
    
    draw(args,dr,gr,"left")

def testRight(args):
    dr = [[5,5,16,6]]
    gr = [[5,1,15,10]]

    print("right: ",getBoxDifference(dr[0],gr[0]))
    
    draw(args,dr,gr,"right")

def testBottom(args):
    dr = [[5,0,8,6]]
    gr = [[5,1,15,10]]

    print("bottom: ",getBoxDifference(dr[0],gr[0]))
    
    draw(args,dr,gr,"bottom")

def testTop(args):
    dr = [[5,8,8,12]]
    gr = [[5,1,15,10]]

    print("top: ",getBoxDifference(dr[0],gr[0]))
    
    draw(args,dr,gr,"top")

def testLeftToRight(args):
    dr = [[5,8,8,12]]
    gr = [[5,1,15,10]]

    print("top: ",getBoxDifference(dr[0],gr[0]))
    
    draw(args,dr,gr,"top")



def outofguideInit(args):
    # dr= [[1,5,10,6]]
    # gr= [[5,1,15,10],[3,1,4,10],[3,7,10,8]]

    # no intersection
    testNointersection(args)
    testWithin(args)
    testLeft(args)
    testRight(args)
    testBottom(args)
    testTop(args)