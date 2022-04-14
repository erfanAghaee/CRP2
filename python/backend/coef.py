#!/usr/bin/env python

import argparse
from cProfile import label

import os
from pydoc import describe
import sys
import pandas as pd
import numpy as np

'''
Python Project Source File imports 
'''

import_path = os.path.abspath(os.path.join(os.path.join(__file__, ".."), ".."))
sys.path.insert(0, import_path)


from backend.db import *
from backend.param import *
from backend.utils import *
import matplotlib.pyplot as plt



class Coef:
    def __init__(self,args,num_iter):
        self.args = args
        self.init(num_iter)
        self.num_rows = int(num_iter/2)
        self.num_cols = int(num_iter/2)
        plt.rcParams["figure.figsize"] = (16,6)

    
    def init(self,num_iter):
        
        args = self.args
        directory = args.dir+"imgs/"+"coef"
        if not os.path.exists(directory):
            os.mkdir(directory)
        df = pd.DataFrame({})
        for i in np.arange(num_iter):
            db = getDB(args,iter_gr=i,iter_dr=0)
            coef_df = db["coef"]

            if(i == 0):
                df["name"] = coef_df["name"].values
            df["iter"+str(i)] = coef_df["value"].values
        df.to_csv(args.dir+ "imgs/coef/"+"coef.csv",index=False)

    def getCoef(self):
        args = self.args
        df = pd.read_csv(args.dir+ "imgs/coef/"+"coef.csv")
        df = df.set_index(["name"])
        df = df.T
        return df


    def drawLogisticSlope(self):
        args = self.args

        df = self.getCoef()

        num_rows=self.num_rows
        num_cols=self.num_cols

        fig, axes = plt.subplots(num_rows, num_cols, sharex=True, sharey=True)
        
        
        slopes = df["logisticSlope"].values
        x = np.linspace(-10,10,100)
        
        for i, ax in enumerate(axes.flat):
            s = slopes[i]
            # logistic function
            z = 1/(1 + np.exp(-1*s*(x))) 
            ax.plot(x,z)
            ax.set_title(f'iter: {i}, Slope {s}')

        # df["logisticSlope"].plot()
        plt.savefig(args.dir+ "imgs/coef/"+"coef_logisticSlope.png")

    
    
    def drawViaCost(self):
        args = self.args
        df = self.getCoef()
        num_rows=self.num_rows
        num_cols=self.num_cols

        # fig, axes = plt.subplots(num_rows, num_cols, sharex=True, sharey=True)
        fig, axes = plt.subplots(num_rows, num_cols, sharex=True, sharey=False)

        slopes = df["logisticSlope"].values
        unitViaCosts = df["unitViaCost"].values
        unitViaMultipliers = df["unitViaMultiplier"].values

        x = np.linspace(-10,10,100)

        for i, ax in enumerate(axes.flat):
            s = slopes[i]
            unitViaCost = unitViaCosts[i]
            unitViaMultiplier = unitViaMultipliers[i]
            # viaCost
            # database.getUnitViaCost() * (unitViaMultiplier + getViaShortCost(via,debug));
            z = 1/(1 + np.exp(-1*s*(x)))
            viaCost = unitViaCost*(unitViaMultiplier + z)


            # z = 1/(1 + np.exp(-1*s*(x)))
            ax.plot(x,viaCost)
            ax.set_xlabel(f'{subfigure_label[i]}')    
            ax.set_title(f'iter: {i}, unitViaCost: {unitViaCost},unitViaMultiplier: {unitViaMultiplier} , Slope: {s}')

        # df["logisticSlope"].plot()
        fig.suptitle("viaCost = unitViaCost*(unitViaMultiplier + logsticFunction)", fontsize=16)
        plt.savefig(args.dir+ "imgs/coef/"+"coef_viaCost.png")
    
    def drawWireCost(self):
        args = self.args
        df = self.getCoef()
        num_rows=self.num_rows
        num_cols=self.num_cols

        # fig, axes = plt.subplots(num_rows, num_cols, sharex=True, sharey=True)
        fig, axes = plt.subplots(num_rows, num_cols, sharex=True, sharey=False)

        slopes = df["logisticSlope"].values
        unitViaCosts = df["unitViaCost"].values
        unitViaMultipliers = df["unitViaMultiplier"].values
        

        x = np.linspace(-10,10,100)

        for i, ax in enumerate(axes.flat):
            s = slopes[i]
            unitViaCost = unitViaCosts[i]
            unitViaMultiplier = unitViaMultipliers[i]
            # viaCost
            # database.getUnitViaCost() * (unitViaMultiplier + getViaShortCost(via,debug));
            z = 1/(1 + np.exp(-1*s*(x)))
            viaCost = unitViaCost*(unitViaMultiplier + z)

            # wirecost
            #getWireDistCost(edge) + getWireShortCost(edge,debug)
            for l in np.arange(9):
                unitShortVioCostDiscounteds = \
                    df["unitShortVioCostDiscounted"+str(l)].values
                unitShortVioCostDiscounted = unitShortVioCostDiscounteds[i]
                wireCost = 3000 + (z * unitShortVioCostDiscounted  )
                # z = 1/(1 + np.exp(-1*s*(x)))
                ax.plot(x,wireCost,label="layer"+str(l))  
                ax.set_xlabel(f'{subfigure_label[i]}')
            plt.legend()
            
            ax.set_title(f'iter: {i}, Slope: {s}')

        # df["logisticSlope"].plot()
        fig.suptitle("wireCost = getWireDistCost(edge) + getWireShortCost(edge)", fontsize=16)
        plt.savefig(args.dir+ "imgs/coef/"+"coef_wireCost.png")