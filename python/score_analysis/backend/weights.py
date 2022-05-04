import numpy as np
import argparse

import matplotlib.pyplot as plt


class Weight:
    def __init__(self,args):
        self.args = args

    def run(self):
        print("hi")
        pass
       
       
    # def writeToExcel(self,writer,df,name,caption,row):
    #     df.to_excel(writer, sheet_name=name, index=False,startrow=row)

    #     workbook  = writer.book
    #     worksheet = writer.sheets[name]
    #     worksheet.write('A1', caption)

    #     my_format = workbook.add_format()
    #     my_format.set_align('right')

    #     worksheet.set_column('C:D', None, my_format)
    #     worksheet.set_column('G:H', None, my_format)


    #     self.fitCols(writer,df_even,name)
    #     self.fitCols(writer,df_odd,name)