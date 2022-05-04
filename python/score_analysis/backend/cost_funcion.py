from cProfile import label
import numpy as np
import argparse

import matplotlib.pyplot as plt


class CostFunction:
    def __init__(self,args):
        self.args = args

    def run(self):
        self.runPenaltyCahnges()
        
    def runPenaltyCahnges(self):
        args = self.args
        # demand function changes
        d = np.linspace(1,10)
        # capacity number of available tracks
        c = 5
        # GCell size 3000
        exp_len = 3000
        # how it changes in cugr
        slopes = [1,2,4,8,100]

        bs = [0]

        for b in bs:
            for s in slopes:
                label = "slope: " + str(s)
                y = exp_len * (1/(1 + np.exp(-1*s*(d-(c-b)))))
                plt.plot(d,y,label=label)
            
        plt.legend()
        
        plt.savefig("effect_of_slope_on_congestion_penalty.png")
