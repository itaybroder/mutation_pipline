import json
from predictors.config import INPUT_DIR
#splitting RNA to all the possible peptides by certain k
def split_by_k(RNA, k, pos_list):
    peptides = []
    pos_counter = 1
    for s in RNA[:len(RNA)-k]:
        peptides.append([RNA[pos_counter-1: pos_counter+k-1], pos_list[pos_counter-1], pos_list[pos_counter+k-2]])
        pos_counter+=1
    return peptides

#returns a list with all peps
def split_for_all_peptides(str, pos_list):
    return split_by_k(str, 8, pos_list) + split_by_k(str, 9, pos_list) + split_by_k(str, 10, pos_list)

def split_by_k_r(RNA, k):
    peptides = []
    pos_counter = 1
    for s in RNA[:len(RNA)-k]:
        peptides.append([RNA[pos_counter-1: pos_counter+k-1], pos_counter-1, pos_counter+k])
        pos_counter+=1
    return peptides

#returns a list with all peps
def split_for_all_peptides_r(str):
    return split_by_k_r(str, 8) + split_by_k_r(str, 9) + split_by_k_r(str, 10)


def is_vocs(mutation, voc_list):
    for m in voc_list:
        if(m == mutation):
            return True

    return False

def json_txt_to_dict(file_path):
    d = json.load(open(file_path))
    return d

def create_varients_protien(varient_name, mutation_list, protien, pos_list):
    for mut in mutation_list:
        mut_pos = int(mut[1:-1]) 
        new_aa = mut[-1]
        if(new_aa == "-"):
            new_aa = ""
            pos_list.pop(mut_pos-1)
        new_protien = protien
        new_protien = list(new_protien)
        new_protien[mut_pos] = new_aa 
        new_protien = "".join(new_protien)
        return [new_protien, pos_list]


def create_varients_dict(varients_dict):
    new_dict = {}
    
    for varient in varients_dict:
        pos_list = list(range(1, len(PROTEIN)+1))
        new_protien, pos_list = create_varients_protien(varient, varients_dict[varient], PROTEIN, pos_list)
        new_dict[varient] = {'protein': new_protien, 'peptides':split_for_all_peptides(new_protien, pos_list), 'pos_list': pos_list}

    new_dict["original"] = {'protein': PROTEIN, 'peptides': split_for_all_peptides_r(PROTEIN),'pos_list': list(range(1, len(PROTEIN)+1))}
    return new_dict

def create_protien_fasta_file(varient_dict):
    input_file = open("data/input_files/protiens.fasta", "w")
    for varient in varient_dict:
        input_file.write(">" + varient)
        input_file.write("\n")
        input_file.write(varient_dict[varient]["protein"])
        input_file.write("\n")
    input_file.close()



def create_peptides_input_file(base_df):
    peptides = list(base_df["peptide"].unique())
    input_file = open("data/input_files/all_peptides.txt", "w")
    for pep in peptides:
        input_file.write(pep)
        input_file.write("\n")
    input_file.close()

def merge_netmhcpan(mhcpan_df, base_df):
    MHC_TYPES = list(set(mhcpan_df["mhc_type"]))

    for mhc_type in MHC_TYPES[::-1]:
        base_df = base_df.merge(mhcpan_df[mhcpan_df["mhc_type"] == mhc_type][["peptide", "rank"]].rename(
            columns={"rank": mhc_type + "_rank"}), left_on="peptide", right_on="peptide")
    return base_df



def prapere_file_to_glyc(varient_dict):
    input_file = open("data/input_files/glyc_protiens.fasta", "w")
    counter = 1
    for varient in varient_dict:
        input_file.write(">seq" + str(counter))
        input_file.write("\n")
        input_file.write(varient_dict[varient]["protein"])
        input_file.write("\n")
        varient_dict[varient]["alt_name"] = "seq" + str(counter)
        counter += 1
    input_file.close()
    return varient_dict
PROTEIN = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"
