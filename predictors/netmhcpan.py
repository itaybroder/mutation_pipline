import subprocess
import os
import pandas as pd
from config import PATH_TO_NETMHCPAN, OUTPUT_DIR

#uses the tool
def feed_to_NetMHCPan(input_file_path, output_file_path, pep_length):
    path_to_tool = PATH_TO_NETMHCPAN
    if os.path.exists(path_to_tool):
        print('found all paths  :) ')
    else:
        print("coudn't find ", " ", path_to_tool)

    if os.path.exists(output_file_path):
        print('found all paths  :) ')
    else:
        print("coudn't find ", " ", output_file_path)
    if os.path.exists(input_file_path):
         print('found all paths  :) ')
    else:
        print("coudn't find ", " ", input_file_path)


    HLA_str = 'HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A24:02,HLA-A26:01,HLA-B07:02,HLA-B08:01,HLA-B27:05,HLA-B39:01,' \
              'HLA-B40:01,HLA-B58:01,HLA-B15:01'
    for length in pep_length:
        command = path_to_tool + "/netMHCpan -p " + input_file_path + " -xls " + " -l " + str(length) + " -a " + HLA_str + " -s " + " >" + output_file_path
    print(command)
    subprocess.check_output('%s' % command, shell=True)

#creting a pandas DataFrame from the output file of the tool

def create_dataframe_from_netmhcpan(netmhcpan_output_file):
    with open(netmhcpan_output_file, 'r') as f:
        output_file_lines = f.readlines()
        lis=[]
        i = 0
        while(i<len(output_file_lines)):
            
            while(i<len(output_file_lines)):
                if(output_file_lines[i].startswith(" Pos")):
                    i+=2
                    break
                i+=1
            while(i<len(output_file_lines) and (not output_file_lines[i].startswith("-"))):
                line = output_file_lines[i].split()
                mhc_type = line[1]
                peptide = line[2]
                rank = line[12]
                row = [mhc_type, peptide, rank]
                print(row)
                lis.append(row)
                i+=1
        
        cols = ["mhc_type", "peptide", "rank"]
        mhc_frame = pd.DataFrame(lis, columns=cols)
        mhc_frame.to_csv(OUTPUT_DIR + "/mhc_resualt.csv")



