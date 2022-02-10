import pandas as pd
df = pd.read_csv("yair.csv")


def parseNetglyc(file_path, df):
    glyc_dict = {}
    with open(file_path, 'r') as f:
        for line in f.readlines():
            line = line.split()
            glyc_dict[line[2]] = line[3]
    
    potential_list = []
    for pep in df["peptide"]:
        for key,value in glyc_dict.items():
            if key in pep:
                potential_list.append(value)
    print(potential_list)


        

    df["Potential"] = potential_list   
    return df

df  = parseNetglyc("netGlycResult.txt", df)