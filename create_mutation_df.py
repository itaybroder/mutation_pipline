from os import stat
import pandas as pd

import tools
from tools import PROTEIN, split_for_all_peptides, json_txt_to_dict,create_varients_dict, create_protien_fasta_file, create_peptides_input_file
from predictors.netmhcpan import feed_to_NetMHCPan, create_dataframe_from_netmhcpan
from predictors.netchop import create_dataframe_from_netchop
import importlib

importlib.reload(tools)

#prapring the main varients protiens dict
mutation_dict = json_txt_to_dict("data/vocs.txt")
protien_varients_dict = create_varients_dict(mutation_dict)
#-----------------------------------------------------------------------------------------------------------------------

#the main pipeline
def main_pipline(protien_varients_dict):
    #creating the base df
    cols = ["peptide", "start_pos", "end_pos", "varient"]
    rows = []
    for varient in protien_varients_dict:
        for pep in protien_varients_dict[varient]["peptides"]:
            rows.append([pep[0], pep[1], pep[2], varient])

    base_df = pd.DataFrame(rows, columns=cols)

    #creating the input files for the predictors
    create_peptides_input_file(base_df)
    create_protien_fasta_file(protien_varients_dict)

    #run the predictors
    feed_to_NetMHCPan()

    netmhcpan_df = create_dataframe_from_netmhcpan()
#-----------------------------------------------------------------------------------------------------------------------

main_pipline(protien_varients_dict)