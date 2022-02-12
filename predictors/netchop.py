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



#creating a dataframe from the tool
def create_dataframe_from_netchop(output_file):
    rows = []
    directory=dir
    for txt in os.listdir(directory):
        if txt.endswith("out.txt"):
            with  open(directory+"/"+txt, "r+") as file:
                lis=[]
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
                        continue
                    elif line == "":
                        continue
                    else:
                        line = line.split()
                        pos = int(line[0])
                        choped = (line[2] == 'S')
                        curr_seq = line[4]
                        if(curr_seq == seq):
                            a = df["end_pos"]
                            df.loc[df['end_pos'] == pos, 'Chopped'] = choped
                        lis.append(line)
                
                return df