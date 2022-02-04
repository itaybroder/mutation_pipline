from os import stat
import pandas as pd
from tools import PROTEIN, split_for_all_peptides

#getting a list of all mutations from the file
mutation_list = []
with open("mutation.txt", "r") as f:
    lines = f.readlines()

for line in lines:
    mutation_list.append(line.replace(",", "").replace("'", "").replace("\n", "").replace(" ", ""))
#----------------------------------------------------------------------------------------------

#getting all the peptides from the original protein
peptides = split_for_all_peptides(PROTEIN)
#---------------------------------------------------------------------------------------------

cols = ["mutation poition", "mutation", "start_pos", "end_pos", "original_peptide", "after_mutation"]
rows = []

for peptide in peptides:
    for mut in mutation_list:
        mut_pos = int(mut[1:-1])
        new_aa = mut[-1]
        start_pos = int(peptide[1])
        end_pos = int(peptide[2])
        if(mut_pos >= start_pos and  end_pos >= mut_pos):
            relative_index = mut_pos-start_pos-1
            new_peptide = peptide[0]
            new_peptide = list(new_peptide)
            new_peptide[relative_index] = new_aa 
            new_peptide = "".join(new_peptide)
            rows.append([mut_pos, new_aa, start_pos, end_pos, peptide[0], new_peptide])
    
df = pd.DataFrame(rows, columns=cols)
df.to_csv("mutation_peptide_df.csv")
    
