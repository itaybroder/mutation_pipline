import pandas as pd

name_map = {"seq1": "Alpha B.1.1.7",
            "seq2": "Beta B.1.351",
            "seq3": "Gamma P.1",
            "seq4": "Kappa B.1.617.1",
            "seq5": "Delta B.1.617.2",
            "seq6": "Lambda C.37",
            "seq7": "Mu B.1.621",
            "seq8": "Omicron BA.1",
            "seq9": "Omicron BA.2",
            "seq10": "AY.4.2",
            "seq11": "C.1.2",
            "seq12": "Eta B.1.525",
            "seq13": "Iota B.1.526",
            "seq14": "original"}

def parseNetglyc(file_path):
    cols = ["varient", "position", "peptide", "Potential"]
    rows = []
    with open(file_path, 'r') as file_:
        for line in file_:
            if(line.startswith("seq")):
                elements = line.split()
                seq = name_map[elements[0]]
                pos = elements[1]
                pep = elements[2]
                potential = elements[3]
                rows.append([seq, pos, pep, potential])
    return pd.DataFrame(rows, columns=cols)


# df = parseNetglyc("../data/output_files/glyc_output.txt")
# print("HEllo")






