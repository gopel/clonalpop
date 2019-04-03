# Importations
import os

def concatenation (folder) :
    list_files = os.listdir(folder)
    new_list = []
    for element in list_files :
        if "fastq" in element :
            new_list.append(element)
    list_files = new_list
    if ".gz" in list_files[1]:
        input = "cat "
        for element in list_files :
            input += folder + "/" + element + " "
        input += " > " + folder + "/concatenated_final.fastq.gz"
        os.system(input)
    else :
        with open(folder + "/concatenated_final.fastq", "w") as outfile:
            l = [folder + "/" + element for element in list_files]
            for fname in l:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
    return



# Illumina samples
def assembler_spades_illumina(input_path, output_path, list_spades_illumina,fichier):
    for element in list_spades_illumina :
        dir = os.listdir(input_path + "/" + element)
        new_dir = []
        for subfile in dir:
            if "fastq" in subfile:
                new_dir.append(subfile)
        if "_R1" in new_dir[0] or "_R2" in new_dir[0]:
            pre_dir_1 = dir
            pre_dir_2 = []
            dir = [0, 0]
            for pre_element in pre_dir_1:
                if 'fastq' in pre_element:
                    pre_dir_2.append(pre_element)
            for pre_element in pre_dir_2:
                if 'R1' in pre_element:
                    dir[0] = pre_element
                if 'R2' in pre_element:
                    dir[1] = pre_element
            os.system("spades.py --pe1-1  " + input_path + "/" + element + "/" + new_dir[0]
                      + " --pe1-2 " + input_path + "/" + element + "/" + new_dir[1]
                      + " -o " + output_path + "/" + element + "/Spades")
        else:
            if len(new_dir) == 1:
                os.system("spades.py --s1 " + input_path + "/" + element + "/" + new_dir[0]
                          + " -o " + output_path + "/" + element + "/Spades")
            else:
                concatenation(input_path + "/" + element)
                file = ""
                for mini_file in os.listdir(input_path + "/" + element):
                    if "concatenated_final" in mini_file:
                        file+= mini_file
                os.system("spades.py --s1 " + input_path + "/" + element + "/" + file
                          + " -o " + output_path + "/" + element + "/Spades")
        fichier.write("\n       " + element + " has been assembled")
    return


# Illumina unmatched samples      A CORRIGER
def assembler_spades_unmatched_illumina(output_path, acces_dossier_compare, reference_file):
    for element in acces_dossier_compare:
        dir = os.listdir(output_path + "/" + element + "/Breseq/" + reference_file + "/data")
        new_dir = []
        for subfile in dir:
            if "fastq" in subfile:
                new_dir.append(subfile)
        if len(new_dir) == 1:
            os.system("spades.py --s1 " + output_path + "/" + element + "/Breseq/" + reference_file + "/data/" + new_dir[0]
                      + " -o " + output_path + "/" + element + "/Spades_unmatched/" + reference_file)
        else:
            concatenation(output_path + "/" + element + "/Breseq/" + reference_file + "/data")
            file = ""
            for mini_file in os.listdir(output_path + "/" + element + "/Breseq/" + reference_file + "/data"):
                if "R1" in mini_file:
                    file+= mini_file
                if "R2" in mini_file:
                    file += mini_file
            os.system("spades.py --s1 " + output_path + "/" + element + "/Breseq/" + reference_file + "/data/concatenated_final.fastq"""
                      + " -o " + output_path + "/" + element + "/Spades_unmatched/" + reference_file)
    return


# Nanopore samples
def assembler_nanopore(input_path, output_path, list_nanopore, fichier) :
    for element in list_nanopore :
        dir = os.listdir(input_path + "/" + element)
        if len(dir) > 1 :
            if "concatenated_final.fastq" not in dir :
                concatenation(input_path + "/" + element)
        new_element = element[0:-9] + "_illumina"
        new_dir = os.listdir(input_path + "/" + new_element)
        new_dir_2 = ['','']
        for element_dir in new_dir :
            if "fastq" in element_dir :
                if "_R1" in element_dir :
                    new_dir_2[0] = element_dir
                elif "_R2" in element_dir :
                    new_dir_2[1] = element_dir
        os.system("unicycler "
                  + " -1 " + input_path + "/" + new_element + "/" + new_dir_2[0]
                  + " -2 " + input_path + "/" + new_element + "/" + new_dir_2[1]
                  + " -l " + input_path + "/" + element + "/concatenated_final.fastq "
                  + " -o " + output_path + "/" + element + "/Unicycler")
        fichier.write("\n       " + element + " has been assembled")
    return
