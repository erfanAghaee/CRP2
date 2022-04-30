import pandas as pd
import numpy as np
import xlsxwriter


flow_naming = {
    "PL: Contest; GR: Contest; DR: DRCU" : "Con+Con+DRCU",
    "PL: CRP(10); GR: CRP(10); DR: Triton (path_cost; PR last iteration of GR is applied)" : "Con+CRP10(PathCost)+Triton",
    "PL: CRP(1); GR: CRP(1); DR: Triton (path_cost; PR last iteration of GR is applied)" : "Con+CRP1(PathCost)+Triton",
    "PL: Contest; GR: CUGR; DR: Triton" : "Con+CUGR+Triton",
    "PL: Contest; GR: Contest; DR: Triton" : "Con+Con+Triton",
    "PL: Eh?Placer; GR: CUGR; DR: Triton" : "Con+EhPlacer+CUGR+Triton",
    "PL: DreamPlace; GR: CUGR; DR: Triton" : "Con+DreamPlacer+CUGR+Triton",
    "PL: ILP-Based; GR: ILP-based; DR: Triton" : "Con+ILP-based+Triton",
    "PL: CRP(10); GR: CRP(10); DR: Triton (wl_abs_cost; PR first iteration of GR is applied)" : "Con+CRP10(WLCost)+Triton",
    "PL: Contest; GR: Contest; DR: LuLuRoute" : "Con+Con+LuLu"
}




def getGap(x):
    
    print(x)
    # print(len(x))
    gaps = []
    gaps.append(0)
    for i in range(1,len(x)):
        # print(x)
        if((x[i] != "dnf") and (x[i] != "na")):
            expr = (int(x[i]) - int(x[0]))
            gaps.append(expr)
        else:
            if(x[i] != "dnf"):
                gaps.append("dnf")
            elif(x[i] != "na"): 
                gaps.append("na")

    return gaps

def getNumeric(val):
    if((val == "dnf") or ( val == "na")):
        return float('inf')
    else:
        return int(val)
    
    
def getMetric(benchs,metric):
    dfs = []
    for grp in benchs.groups:
    
        df = benchs.get_group(grp)
        # print(grp)
            # print( df)
        df['score_num'] = df.apply(lambda row: getNumeric(row[metric]),axis=1)
            
        df = df.sort_values(by=["score_num"])

        df_new = df[["benchmark","flow",metric]].reset_index(drop=True)
        print("len: ",len(df_new),", lenGetGaps: ",len(getGap(df_new[metric].values)))
        df_new["gap"] = getGap(df_new[metric].values)
        df_new["rank"] = np.arange(1,len(df_new)+1,dtype=int)

        df_new = pd.DataFrame(df_new,columns = ["benchmark","rank","flow",metric,"gap"])
        # 
        dfs.append(df_new)


    df = pd.concat(dfs)
    a = df["benchmark"].values
    new_idx = [int(i.split("ispd18_test")[-1]) for i in a]
    df["index"] =new_idx
    df = df.sort_values(by=["index","rank"])
    df = df.drop(["index"],axis=1)
    return df


def changeNamingFlowFunc(row):
    return flow_naming[row]

def changeNamingFlow(df):
    df["flow"] = df.apply(lambda row: changeNamingFlowFunc(row["flow"]),axis=1)

def getOddEven(df_score):
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

def fitCols(writer,df,name):
    for column in df:
        column_width = max(df[column].astype(str).map(len).max(), len(column))
        col_idx = df.columns.get_loc(column)
        writer.sheets[name].set_column(col_idx, col_idx, column_width)


def writeToExcel(writer,df_odd,df_even,name,caption,row):
    df_odd.to_excel(writer, sheet_name=name, index=False,startrow=row)
    df_even.to_excel(writer, sheet_name=name, index=False,startrow=row,startcol=len(df_odd.columns))
    workbook  = writer.book
    worksheet = writer.sheets[name]
    worksheet.write('A1', caption)

    my_format = workbook.add_format()
    my_format.set_align('right')

    worksheet.set_column('D:E', None, my_format)
    worksheet.set_column('I:J', None, my_format)


    fitCols(writer,df_even,name)
    fitCols(writer,df_odd,name)
    

def writeToExcelAll(writer,df_odd,df_even,name,caption):
    df_odd_grps = df_odd.groupby(["benchmark"])
    df_even_grps = df_even.groupby(["benchmark"])
    benchs_odd = ["ispd18_test1","ispd18_test3","ispd18_test5","ispd18_test7","ispd18_test9"]
    benchs_even = ["ispd18_test2","ispd18_test4","ispd18_test6","ispd18_test8","ispd18_test10"]
    
    row = 2
    for i in range(len(benchs_even)):
        df_odd_tmp = df_odd_grps.get_group(benchs_odd[i])
        df_even_tmp = df_even_grps.get_group(benchs_even[i])
        writeToExcel(writer,df_odd_tmp,df_even_tmp,name,caption,row)


        row = row + len(df_odd_tmp)+2
    
    workbook  = writer.book
    worksheet = writer.sheets[name]
    merge_format = workbook.add_format({
        'align': 'center',
        'valign': 'vcenter'})

    worksheet.merge_range('B2:E2', 'ispd18_test1', merge_format)
    worksheet.merge_range('G2:J2', 'ispd18_test2', merge_format)
    worksheet.merge_range('B14:E14', 'ispd18_test3', merge_format)
    worksheet.merge_range('G14:J14', 'ispd18_test4', merge_format)
    worksheet.merge_range('B26:E26', 'ispd18_test5', merge_format)
    worksheet.merge_range('G26:J26', 'ispd18_test6', merge_format)
    worksheet.merge_range('B38:E38', 'ispd18_test7', merge_format)
    worksheet.merge_range('G38:J38', 'ispd18_test8', merge_format)
    worksheet.merge_range('B50:E50', 'ispd18_test9', merge_format)
    worksheet.merge_range('G50:J50', 'ispd18_test10', merge_format)
    


def menu():
    ispd18 = pd.read_csv("ispd18_scores.csv")
    benchs = ispd18.groupby(["benchmark"])
    writer = pd.ExcelWriter('ispd18_scores.xlsx', engine='xlsxwriter')

    all_tables = ['wl_gr', 'vias_gr', 'time_gr', 'mem_gr',
        'wl_dr', 'vias_dr', 'OFGW', 'OFGV', 'OFTW', 'OFTV', 'WWW', 'short_area',
        'min_area', 'spacing', 'score', 'time_dr', 'time_mem']


    for col in all_tables:
        print(col)
        df_score = getMetric(benchs=benchs,metric=col)
        changeNamingFlow(df_score)
        df_score_odd, df_score_even = getOddEven(df_score=df_score)

        writeToExcelAll(writer,df_score_odd,df_score_even,name="ispd18_"+col,caption="Ranking " + col)    
        
    writer.save()


if __name__ == "__main__":
    menu()
