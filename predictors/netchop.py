import os
import subprocess
from predictors.config import PATH_TO_NETCHOP, INPUT_DIR, OUTPUT_DIR

def feed_to_Netchop():
    input_file_path = INPUT_DIR + "/protiens.fasta"
    output_file_path=OUTPUT_DIR + "/netchop_output.txt"
    path_to_tool=PATH_TO_NETCHOP
    command = path_to_tool + "/netchop " + input_file_path + " > " + output_file_path
    print(command + " :)")
    subprocess.check_output('%s' % command, shell=True)


def seq_map(mutation_dict, varient_name):
    for varient in mutation_dict:
        if(varient_name in varient):
            return varient


#creating a dataframe from the tool
def create_dataframe_from_netchop(base_df, protien_varients_dict):
    output_file = "data/output_files/netchop_output.txt"
    with open(output_file, "r+") as file:
        omicronflag=False
        omicron2flag = False
        for line in file:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            elif line.startswith("-"):
                continue
            elif line.startswith("N"):
                continue
            elif line.startswith("p"):
                continue
            elif line.startswith(" pos  AA  C      score      Ident"):
                if(omicronflag):
                    omicron2flag = True
                if(omicron2flag):
                    omicron2flag = False
                continue
            elif line == "":
                continue
            else:
                
                line = line.split()
                end_pos = int(line[0])
                chopped = (line[2] == 'S')
                curr_seq = seq_map(protien_varients_dict, line[4])
                if(curr_seq == "Omicron BA.1"):
                    omicronflag=True
                if(omicron2flag == True):
                    curr_seq = "Omicron BA.2"
                base_df.loc[(base_df["end_pos"] == end_pos) & (base_df["varient"] == curr_seq), "Chopped"] = chopped



        return base_df