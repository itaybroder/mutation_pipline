import json

#splitting RNA to all the possible peptides by certain k
def split_by_k(RNA, k, pos_list):
    peptides = []
    pos_counter = 1
    for s in RNA[:len(RNA)-k]:
        peptides.append([RNA[pos_counter-1: pos_counter+k-1], pos_counter, pos_counter+k])
        pos_counter+=1
    return peptides

#returns a list with all peps
def split_for_all_peptides(str):
    return split_by_k(str, 8) + split_by_k(str, 9) + split_by_k(str, 10)

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
        new_dict[varient] = {'protein': new_protien, 'peptides':split_for_all_peptides(new_protien), 'pos_list': pos_list}
    return new_dict

PROTEIN = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"
