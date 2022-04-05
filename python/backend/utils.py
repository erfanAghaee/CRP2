import pandas as pd

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