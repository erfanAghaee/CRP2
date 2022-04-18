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
import copy


def plot(plt_obj,boxs,color,alpha,type):
    try:
        if len(boxs[0]) != 4:
            return
    except:
        return

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

        # unionbox = getUnionBox(ovrlp_boxs)

        for box in ovrlp_boxs:
            print(box,dr)
            # print("intersected: ",intersect(dr[0],box))
        
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




def draw(args,dr,gr,outofguide,name,window=[0,0,20,20]):
    # window = [
    #     0,0,20,20
    # ]
    xs=[];ys=[]
    for i in range(len(dr)):
        xs.append(dr[i][XL]);xs.append(dr[i][XH])
    for i in range(len(gr)):
        xs.append(gr[i][XL]);xs.append(gr[i][XH])
    for i in range(len(dr)):
        ys.append(dr[i][YL]);ys.append(dr[i][YH])
    for i in range(len(gr)):
        ys.append(gr[i][YL]);ys.append(gr[i][YH])

    window = [np.min(xs),np.min(ys),np.max(xs),np.max(ys)]
    

    # gcells = initGCells(window)
    print(gr)
    print(dr)
    print(window)
    


    plt_obj = PltCairo()
    surface = plt_obj.init(window)

    
    plot(plt_obj,dr,(0,1,0),1,"dr")
    plot(plt_obj,gr,(0,0,1),0.2,"gr")
    # plot(plt_obj,gcells,(1,0,0),0.1,"gcells")
    plot(plt_obj,outofguide,(1,0,0.3),1,"outofguide")
    

    # checkoutofguide(dr,gr)

    directory = args.dir+"imgs/"+"outofguides"
    if not os.path.exists(directory):
        os.mkdir(directory)

    surface.write_to_png(directory + "/"+name +".png")

def testNointersection(args):
    dr = [[1,5,4,6]]
    gr = [[5,1,15,10]]

    outofguide = getBoxDifference(dr[0],gr[0])
    
    draw(args,dr,gr,outofguide,"noIntersection")

def testWithin(args):
    dr = [[5,5,10,6]]
    gr = [[5,1,15,10]]

    outofguide = getBoxDifference(dr[0],gr[0])
    
    draw(args,dr,gr,outofguide,"within")
   
def testLeft(args):
    dr = [[3,5,10,6]]
    gr = [[5,1,15,10]]

    outofguide = getBoxDifference(dr[0],gr[0])
    
    draw(args,dr,gr,outofguide,"left")

def testRight(args):
    dr = [[5,5,16,6]]
    gr = [[5,1,15,10]]

    outofguide = getBoxDifference(dr[0],gr[0])
    
    draw(args,dr,gr,outofguide,"right")

def testBottom(args):
    dr = [[5,0,8,6]]
    gr = [[5,1,15,10]]

    outofguide = getBoxDifference(dr[0],gr[0])
    
    draw(args,dr,gr,outofguide,"bottom")

def testTop(args):
    dr = [[5,8,8,12]]
    gr = [[5,1,15,10]]

    outofguide = getBoxDifference(dr[0],gr[0])
    
    draw(args,dr,gr,outofguide,"top")

def testLeftToRight(args):
    dr = [[5,8,8,12]]
    gr = [[5,1,15,10]]

    outofguide = getBoxDifference(dr[0],gr[0])
    
    draw(args,dr,gr,outofguide,"top")

def testLefttoRight(args):
    dr = [[4,8,16,10]]
    gr = [[5,1,15,10]]

    outofguide = getBoxDifference(dr[0],gr[0])
    
    draw(args,dr,gr,outofguide,"leftright")
def testUptoDown(args):
    dr = [[5,0,8,12]]
    gr = [[5,1,15,10]]

    outofguide = getBoxDifference(dr[0],gr[0])
    
    draw(args,dr,gr,outofguide,"updown")

def testBottomLeft(args):
    dr = [[4,0,8,2]]
    gr = [[5,1,15,10]]

    outofguide = getBoxDifference(dr[0],gr[0])
    
    draw(args,dr,gr,outofguide,"bottomleft")
def testBottomRight(args):
    dr = [[13,0,16,2]]
    gr = [[5,1,15,10]]

    outofguide = getBoxDifference(dr[0],gr[0])
    
    draw(args,dr,gr,outofguide,"bottomRight")
def testTopRight(args):
    dr = [[13,9,16,12]]
    gr = [[5,1,15,10]]

    outofguide = getBoxDifference(dr[0],gr[0])
    
    draw(args,dr,gr,outofguide,"topRight")

def testTopLeft(args):
    dr = [[4,9,6,12]]
    gr = [[5,1,15,10]]

    outofguide = getBoxDifference(dr[0],gr[0])
    
    draw(args,dr,gr,outofguide,"topLeft")

def testTopLeftRightCorner(args):
    dr = [[4,9,16,12]]
    gr = [[5,1,15,10]]

    outofguide = getBoxDifference(dr[0],gr[0])
    
    draw(args,dr,gr,outofguide,"topLeftRightCorner")

def testLeftTopBottomCorner(args):
    dr = [[4,0,6,12]]
    gr = [[5,1,15,10]]

    outofguide = getBoxDifference(dr[0],gr[0])
    
    draw(args,dr,gr,outofguide,"leftTopBottomCorner")

def testRightTopBottomCorner(args):
    dr = [[14,0,16,12]]
    gr = [[5,1,15,10]]

    outofguide = getBoxDifference(dr[0],gr[0])
    
    draw(args,dr,gr,outofguide,"rightTopBottomCorner")

def testBottomLeftRightCorner(args):
    dr = [[4,0,16,2]]
    gr = [[5,1,15,10]]

    outofguide = getBoxDifference(dr[0],gr[0])
    
    draw(args,dr,gr,outofguide,"bottomLeftRightCorner")


def testWrap(args):
    dr = [[4,0,16,12]]
    gr = [[5,1,15,10]]

    outofguide = getBoxDifference(dr[0],gr[0])
    
    draw(args,dr,gr,outofguide,"wrap")


def testRandomBoxs(args,n):
    rng = 10
    for i in range(n):
        dr_xl = np.random.randint(0,rng); dr_yl = np.random.randint(0,rng);
        dr_w = np.random.randint(0,rng); dr_h = np.random.randint(0,rng);

        print(dr_xl,dr_yl,dr_xl+dr_w,dr_yl+dr_h)

        # gr_xl = np.random.randint(0,rng); gr_yl = np.random.randint(0,rng);
        # gr_w = np.random.randint(0,rng); gr_h = np.random.randint(0,rng);
        gr = [[5,1,15,10]]

        dr = [[dr_xl,dr_yl,dr_xl+dr_w,dr_yl+dr_h]]
        # gr = [[gr_xl,gr_yl,gr_xl+gr_w,gr_yl+gr_h]]

        outofguide = getBoxDifference(dr[0],gr[0])
    
        draw(args,dr,gr,outofguide,"test"+str(i))

def testMultiBlock(args):
    gr = [[3,1,4,10],[5,1,8,10]]
    dr = [[1,3,9,4],[8,4,9,5],[7,4,9,5]]



    dr=[
        [710350.0,   525290.0,  710450.0,   525510.0],
        [466790.0,  1130750.0,  467010.0,  1130850.0],
        [587450.0,  1170690.0,  587550.0,  1170910.0],
        [578050.0,  1164290.0,  578150.0,  1164510.0],
        [671850.0,  1170690.0,  671950.0,  1170910.0],
        [438590.0,   667750.0,  438810.0,   667850.0],
        [724650.0,  1161090.0,  724750.0,  1161310.0],
        [824650.0,  1159490.0,  824750.0,  1159710.0],
        [852650.0,  1164290.0,  852750.0,  1164510.0],
        [845050.0,   493890.0,  845150.0,   494110.0]
    ]


    gr= [
        [843000,   486100,  846000,   498100],
        [708000,   525100,  711000,   528100],
        [438000,   666100,  441000,   669100],
        [465000,  1128100,  468000,  1131100],
        [723000,  1158100,  726000,  1161100],
        [822000,  1158100,  825000,  1161100],
        [576000,  1164100,  579000,  1167100],
        [852000,  1164100,  855000,  1167100],
        [585000,  1170100,  588000,  1173100],
        [669000,  1170100,  672000,  1173100]
    ]

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

    
    draw(args,dr,gr,dr_res,"multi")


def outofguideModuleTest(args):
    # dr= [[1,5,10,6]]
    # gr= [[5,1,15,10],[3,1,4,10],[3,7,10,8]]

    # no intersection
    # testNointersection(args)
    # testWithin(args)
    # testLeft(args)
    # testRight(args)
    # testBottom(args)
    # testTop(args)
    # testLefttoRight(args)
    # testUptoDown(args)
    # testBottomLeft(args)
    # testBottomRight(args)
    # testTopRight(args)
    # testTopLeft(args)
    # testTopLeftRightCorner(args)
    # testLeftTopBottomCorner(args)
    # testRightTopBottomCorner(args)
    # testBottomLeftRightCorner(args)
    # testWrap(args)
    # testRandomBoxs(args,10)
    testMultiBlock(args)
    pass



def getOutOfGuides(args,dr,gr):
    # gr = [[3,1,4,10],[5,1,8,10]]
    # dr = [[1,3,9,4],[8,4,9,5],[7,4,9,5]]
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