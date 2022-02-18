from os import stat
import pandas as pd
import tools
from tools import PROTEIN, split_for_all_peptides, json_txt_to_dict, create_varients_dict, create_protien_fasta_file, \
    create_peptides_input_file, merge_netmhcpan, prapere_file_to_glyc, merge_base_df_with_glyc
from predictors.netmhcpan import feed_to_NetMHCPan, create_dataframe_from_netmhcpan
from predictors.netchop import create_dataframe_from_netchop, feed_to_Netchop
from predictors.glyc import parseNetglyc
import importlib

# prapring the main varients protiens dict
mutation_dict = json_txt_to_dict("data/vocs.txt")
protien_varients_dict = create_varients_dict(mutation_dict)
protien_varients_dict = prapere_file_to_glyc(protien_varients_dict)
# -----------------------------------------------------------------------------------------------------------------------

# the main pipeline
def main_pipline(protien_varients_dict):
    # creating the base df
    cols = ["peptide", "start_pos", "end_pos", "varient", "Chopped"]
    rows = []
    for varient in protien_varients_dict:
        for pep in protien_varients_dict[varient]["peptides"]:
            rows.append([pep[0], pep[1], pep[2], varient, False])

    base_df = pd.DataFrame(rows, columns=cols)
    print(base_df.head())
    # # creating the input files for the predictors
    # create_peptides_input_file(base_df)
    # create_protien_fasta_file(protien_varients_dict)
    
    # run the predictors
    feed_to_NetMHCPan()
    feed_to_Netchop()
    feed_to_Netchop()

    # base_df = create_dataframe_from_netchop(base_df, protien_varients_dict)
    # netmhcpan_df = create_dataframe_from_netmhcpan()
    # base_df = merge_netmhcpan(netmhcpan_df, base_df)
    # glyc_df = parseNetglyc("data/output_files/glyc_output.txt")
    # base_df = merge_base_df_with_glyc(base_df, glyc_df)
    # base_df.to_csv("data/output_files/final_base_df.csv")
# -----------------------------------------------------------------------------------------------------------------------


main_pipline(protien_varients_dict)

# base_df = pd.read_csv("data/output_files/final_base_df.csv")

# base_df.loc[base_df.varient == "original", "end_pos"] -=1
# base_df.loc[base_df.varient == "original", "start_pos"] +=1

# base_df.to_csv("data/output_files/final_base_df.csv")