import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from termcolor import colored
from pathlib import Path
import copy
import cairocffi as cairo
import math
import time
import collections 
import seaborn as sns

import_path = os.path.abspath(os.path.join(os.path.join(__file__, ".."), ".."))
sys.path.insert(0, import_path)

from backend.param import *


class PltCairo:         

    def __init__(self):
        pass


    def bin_xl(self,id_x):
        """
        @param id_x horizontal index 
        @return bin xl
        """
        return self.xl + id_x * self.bin_size_x

    def bin_xh(self,id_x):
        """
        @param id_x horizontal index 
        @return bin xh
        """
        return min(self.bin_xl(id_x) + self.bin_size_x, self.xh)

    def bin_yl(self,id_y):
        """
        @param id_y vertical index 
        @return bin yl
        """
        return self.yl + id_y * self.bin_size_y

    def bin_yh(self,id_y):
        """
        @param id_y vertical index 
        @return bin yh
        """
        return min(self.bin_yl(id_y) + self.bin_size_y, self.yh)

    def normalize_x(self,xx):
        return ((xx - (self.layout_xl - self.padding * self.bin_size_x)) / (
            self.layout_xh - self.layout_xl + self.padding * 2 * self.bin_size_x)) * self.width

    def normalize_y(self,xx):
        return (xx - (self.layout_yl - self.padding * self.bin_size_y)) / (
            self.layout_yh - self.layout_yl + self.padding * 2 * self.bin_size_y) * self.height
       

    def draw_rect(self,x1, y1, x2, y2):
        self.ctx.move_to(x1, y1)
        self.ctx.line_to(x1, y2)
        self.ctx.line_to(x2, y2)
        self.ctx.line_to(x2, y1)
        self.ctx.close_path()
        self.ctx.stroke()


    def drawLayout(self):
        self.ctx.set_source_rgb(0.8, 0.8, 0.8)
        # 225, 227, 225
        self.draw_layout_xl = self.normalize_x(self.layout_xl - self.padding * self.bin_size_x)
        self.draw_layout_yl = self.normalize_y(self.layout_yl - self.padding * self.bin_size_y)
        self.draw_layout_xh = self.normalize_x(self.layout_xh + self.padding * self.bin_size_x)
        self.draw_layout_yh = self.normalize_y(self.layout_yh + self.padding * self.bin_size_y)
        self.ctx.rectangle(self.draw_layout_xl, self.draw_layout_yl, self.draw_layout_xh,
                        self.draw_layout_yh)
        self.ctx.fill()
        self.ctx.set_line_width(self.line_width)
        self.ctx.set_source_rgba(0.1, 0.1, 0.1, alpha=0.8)
        self.ctx.move_to(self.normalize_x(self.xl), self.normalize_y(self.yl))
        self.ctx.line_to(self.normalize_x(self.xl), self.normalize_y(self.yh))
        self.ctx.line_to(self.normalize_x(self.xh), self.normalize_y(self.yh))
        self.ctx.line_to(self.normalize_x(self.xh), self.normalize_y(self.yl))
        self.ctx.close_path()
        self.ctx.stroke()

    def init(self,widow):
        self.xl = int(widow[XL])
        self.yl = int(widow[YL])
        self.xh = int(widow[XH])
        self.yh = int(widow[YH])

        
        if self.xh - self.xl < self.yh - self.yl:
            self.height = 800 
            self.width = round(self.height * (self.xh - self.xl) / (self.yh - self.yl))
        else:
            self.width = 800 
            self.height = round(self.width * (self.yh - self.yl) / (self.xh -self.xl))
        # self.width = 800
        # self.height = 800
        self.line_width = 0.1

        self.layout_xl = self.xl
        self.layout_yl = self.yl
        self.layout_xh = self.xh
        self.layout_yh = self.yh

        self.padding = 0

        self.bin_size_x = 1
        self.bin_size_y = 1

        # init surface 
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, self.width, self.height)
        self.ctx = cairo.Context(surface)

        self.drawLayout()

        return surface

    
    def drawText(self,texts):
        self.ctx.set_font_size(24)
        self.ctx.select_font_face("monospace", cairo.FONT_SLANT_NORMAL,
                                cairo.FONT_WEIGHT_NORMAL)

        for i in np.arange(len(self.node_xl)):
            x = self.node_xl[i]
            y = self.node_yh[i]
            self.ctx.move_to(x,y)
            self.ctx.set_source_rgba(1, 0, 0, alpha=1)
            self.ctx.show_text(texts[i])

    # x,y coordinations and width and height of boxes
    def run(self,xs,ys,ws,hs):
        self.x= np.array(xs)
        self.y= np.array(ys)

        self.node_size_x=np.array(ws)
        self.node_size_y=np.array(hs)
       
        self.draw()

    def draw(self):
            # self.ctx.set_font_size(16)
        # self.ctx.select_font_face("monospace", cairo.FONT_SLANT_NORMAL,
        #                         cairo.FONT_WEIGHT_NORMAL)

        
        self.node_xl = self.x
        self.node_yl = self.layout_yl + self.layout_yh - (self.y + self.node_size_y[0:len(self.y)]
                                            )  # flip y
        self.node_xh = self.node_xl + self.node_size_x[0:len(self.x)]
        self.node_yh = self.layout_yl + self.layout_yh - self.y  # flip y        

        self.node_xl = self.normalize_x(self.node_xl)
        self.node_yl = self.normalize_y(self.node_yl)
        self.node_xh = self.normalize_x(self.node_xh)
        self.node_yh = self.normalize_y(self.node_yh)

        
        # self.ctx.set_line_width(self.line_width)
        # self.ctx.set_line_width(self.line_width *5)
        # self.ctx.set_source_rgba(0, 0, 0.8, alpha=0.5)  # Solid color
        # self.ctx.set_source_rgba(color[0], color[0], color[1], alpha=0.5)  # Solid color
        # self.ctx.set_source_rgba(color[0], color[0], color[1])  # Solid color
        for i in range(len(self.x)):
            self.ctx.rectangle(self.node_xl[i], self.node_yl[i], \
                            self.node_xh[i] - self.node_xl[i],
                            self.node_yh[i] -
                            self.node_yl[i])  # Rectangle(xl, yl, w, h)
            self.ctx.fill()

        self.ctx.set_source_rgba(0.8, 0.8, 0.8, alpha=0.8)  # Solid color
        for i in range(len(self.x)):
            self.draw_rect(self.node_xl[i], self.node_yl[i], self.node_xh[i], self.node_yh[i])
