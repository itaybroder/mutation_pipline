import pandas as pd

def parseNetglyc(file_path):
    cols = ["seq", "position", "pep", "potential"]
    rows = []
    with open(file_path, 'r') as file_:
        for line in file_:
            if(line.startswith("seq")):
                elements = line.split()
                seq = elements[0]
                pos = elements[1]
                pep = elements[2]
                potential = elements[3]
                rows.append([seq, pos, pep, potential])
    return pd.DataFrame(rows, columns=cols)


df = parseNetglyc("../data/output_files/glyc_output.txt")
print("HEllo")






