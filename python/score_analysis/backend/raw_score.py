
import argparse
import os
from pydoc import describe
import sys
import pandas as pd
import numpy as np

import_path = os.path.abspath(os.path.join(os.path.join(__file__, ".."), ".."))
sys.path.insert(0, import_path)


class RawScore:
    def __init__(self,flow_naming):
        self.flow_naming = flow_naming

    def run(self,benchs,metric,writer):
        df_score = self.getMetric(benchs=benchs,metric=metric)
        self.changeNamingFlow(df_score)
        df_score_odd, df_score_even = self.getOddEven(df_score=df_score)
        self.writeToExcelAll(writer,df_score_odd,df_score_even,name="ispd18_"+metric,caption="Ranking " + metric)


    def getGap(self,x):
        # print(len(x))
        gaps = []
        gaps.append(0)
        for i in range(1,len(x)):
            # print(x)
            if((x[i] != "dnf") and (x[i] != "na")):
                expr = (float(x[i]) - float(x[0]))
                gaps.append(expr)
            else:
                if(x[i] == "dnf"):
                    gaps.append("dnf")
                elif(x[i] == "na"): 
                    gaps.append("na")
        
        return gaps

    def getNumericForSort(self,vals):
        vals_float = []
        for val in vals:
            if(val == "dnf"):
                vals_float.append(-1)
            elif(val == "na"):
                vals_float.append(-2)
            else:
                vals_float.append(float(val))

        max_num = np.max(vals_float)
        for i in range(len(vals_float)):
            if(vals_float[i] == -1):
                vals_float[i] = max_num+1
            elif(vals_float[i] == -2):
                vals_float[i] = max_num+2
        return vals_float

    def getFormat(self,val,metric):
        if(val == "dnf"):
            return "dnf"
        elif(val == "na"): 
            return "na"
        else:
            if(metric == "time_gr" or metric == "time_dr"):
                return "{:.2f}".format(float(val))
            else:
                return float(val)
        
        
    def getMetric(self,benchs,metric):
        dfs = []
        for grp in benchs.groups:
    
            df = benchs.get_group(grp)
            df['score_num'] = self.getNumericForSort(df[metric].values)
            # df.apply(lambda row: getNumeric(row[metric]),axis=1)          
            df = df.sort_values(by=["score_num"])
            df_new = df[["benchmark","flow",metric]].reset_index(drop=True)
            # df_new[metric] = df.apply(lambda row: getFormat(row[metric],metric),axis=1)
            df_new["gap"] = self.getGap(df_new[metric].values)
            df_new["rank"] = np.arange(1,len(df_new)+1,dtype=int)
            df_new = pd.DataFrame(df_new,columns = ["benchmark","rank","flow",metric,"gap"])
            dfs.append(df_new)


        df = pd.concat(dfs)
        a = df["benchmark"].values
        new_idx = [int(i.split("ispd18_test")[-1]) for i in a]
        df["index"] =new_idx
        df = df.sort_values(by=["index","rank"])
        df = df.drop(["index"],axis=1)
        return df


    def changeNamingFlowFunc(self,row):
        return self.flow_naming[row]

    def changeNamingFlow(self,df):
        df["flow"] = df.apply(lambda row: self.changeNamingFlowFunc(row["flow"]),axis=1)

    def getOddEven(self,df_score):
        df_score_odd = df_score.loc[(df_score.benchmark=="ispd18_test1") | \
        (df_score.benchmark=="ispd18_test3") | \
        (df_score.benchmark=="ispd18_test5") | \
        (df_score.benchmark=="ispd18_test7") | \
        (df_score.benchmark=="ispd18_test9") ]
        df_score_even = df_score.loc[(df_score.benchmark=="ispd18_test2") | \
            (df_score.benchmark=="ispd18_test4") | \
            (df_score.benchmark=="ispd18_test6") | \
            (df_score.benchmark=="ispd18_test8") | \
            (df_score.benchmark=="ispd18_test10") ]
        return df_score_odd,df_score_even

    def fitCols(self,writer,df,name):
        for column in df:
            column_width = max(df[column].astype(str).map(len).max(), len(column))
            col_idx = df.columns.get_loc(column)
            writer.sheets[name].set_column(col_idx, col_idx, column_width)


    def writeToExcel(self,writer,df_odd,df_even,name,caption,row):
        df_odd.to_excel(writer, sheet_name=name, index=False,startrow=row)
        df_even.to_excel(writer, sheet_name=name, index=False,startrow=row,startcol=len(df_odd.columns))
        workbook  = writer.book
        worksheet = writer.sheets[name]
        worksheet.write('A1', caption)

        my_format = workbook.add_format()
        my_format.set_align('right')

        worksheet.set_column('C:D', None, my_format)
        worksheet.set_column('G:H', None, my_format)


        self.fitCols(writer,df_even,name)
        self.fitCols(writer,df_odd,name)
        

    def writeToExcelAll(self,writer,df_odd,df_even,name,caption):
        df_odd_grps = df_odd.groupby(["benchmark"])
        df_even_grps = df_even.groupby(["benchmark"])
        benchs_odd = ["ispd18_test1","ispd18_test3","ispd18_test5","ispd18_test7","ispd18_test9"]
        benchs_even = ["ispd18_test2","ispd18_test4","ispd18_test6","ispd18_test8","ispd18_test10"]
        
        row = 2
        rows = []
        for i in range(len(benchs_even)):
            df_odd_tmp = df_odd_grps.get_group(benchs_odd[i])
            df_even_tmp = df_even_grps.get_group(benchs_even[i])
            df_even_tmp = df_even_tmp.drop(["benchmark"],axis=1)
            df_odd_tmp = df_odd_tmp.drop(["benchmark"],axis=1)
            self.writeToExcel(writer,df_odd_tmp,df_even_tmp,name,caption,row)
            rows.append(row)
            row = row + len(df_odd_tmp)+2

        # print(rows)
        
        workbook  = writer.book
        worksheet = writer.sheets[name]
        merge_format = workbook.add_format({
            'align': 'center',
            'valign': 'vcenter'})
        bold_format = workbook.add_format({'bold': True})

        
        benchmarks= ["ispd18_test1","ispd18_test2","ispd18_test3","ispd18_test4","ispd18_test5","ispd18_test6",\
            "ispd18_test7","ispd18_test8","ispd18_test9","ispd18_test10"]
        
        mergs_r1_c1="A";mergs_r1_c2="D"
        mergs_r2_c1="E";mergs_r2_c2="H"
        i = 0
        for b in range(int(len(benchmarks)/2)):
            worksheet.merge_range(mergs_r1_c1+str(rows[i])+":"+mergs_r1_c2+str(rows[i]), benchs_odd[b], merge_format)
            worksheet.merge_range(mergs_r2_c1+str(rows[i])+":"+mergs_r2_c2+str(rows[i]), benchs_even[b], merge_format)
            # print(mergs_r1_c1+str(rows[i])+":"+mergs_r1_c2+str(rows[i]))
            # print(mergs_r2_c1+str(rows[i])+":"+mergs_r2_c2+str(rows[i]))
            i = i + 1#len(df_odd) + 2
            

        # worksheet.merge_range('A2:D2', 'ispd18_test1', merge_format)
        # worksheet.merge_range('E2:H2', 'ispd18_test2', merge_format)
        # worksheet.merge_range('A14:D14', 'ispd18_test3', merge_format)
        # worksheet.merge_range('E14:H14', 'ispd18_test4', merge_format)
        # worksheet.merge_range('A26:D26', 'ispd18_test5', merge_format)
        # worksheet.merge_range('E26:H26', 'ispd18_test6', merge_format)
        # worksheet.merge_range('A38:D38', 'ispd18_test7', merge_format)
        # worksheet.merge_range('E38:H38', 'ispd18_test8', merge_format)
        # worksheet.merge_range('A50:D50', 'ispd18_test9', merge_format)
        # worksheet.merge_range('E50:H50', 'ispd18_test10', merge_format)
        