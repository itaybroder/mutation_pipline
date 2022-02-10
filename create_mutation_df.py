from os import stat
import pandas as pd
from tools import PROTEIN, split_for_all_peptides, json_txt_to_dict,create_varients_dict


#setting the main varients protiens dict
mutation_dict = json_txt_to_dict("data/vocs.txt")
protien_varients_dict = create_varients_dict(mutation_dict)
print(protien_varients_dict["Alpha B.1.1.7"]["protein"])



#getting all the peptides from the original protein
peptides = split_for_all_peptides(PROTEIN)
#---------------------------------------------------------------------------------------------

#creating the first mutation dataframe
cols = ["mutation poition", "mutation", "start_pos", "end_pos", "original_peptide", "after_mutation", "Vocs"]
rows = []

for peptide in peptides:
    for mut in mutation_list:
        if("-" in mut):
            continue
        mut_pos = int(mut[1:-1]) 
        new_aa = mut[-1]
        start_pos = int(peptide[1])
        end_pos = int(peptide[2])
        if(mut_pos >= start_pos and  end_pos > mut_pos):
            relative_index = mut_pos-start_pos 
            new_peptide = peptide[0]
            new_peptide = list(new_peptide)
            new_peptide[relative_index] = new_aa 
            new_peptide = "".join(new_peptide)
            rows.append([mut_pos, new_aa, start_pos, end_pos, peptide[0], new_peptide, is_vocs(mut, voc_mutation_list)])
    
df = pd.DataFrame(rows, columns=cols)
df.to_csv("./data/mutation_peptide_df.csv")
#------------------------------------------------------------------------------------------------------------------