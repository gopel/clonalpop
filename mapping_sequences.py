# -*-coding:Latin-1 -*
# Importations
import os

def concatenation (folder) :
    list_files = os.listdir(folder)
    with open(folder + "/concatenated.fastq","w") as outfile :
        l = [folder + "/" + element for element in list_files]
        for fname in l:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
    return


# illumina + nanopore
def breseq(input_path, output_path, input_nanopore, input_illumina, fichier): # deux premiers = pour la pipeline    # deux suivants = pour breseq
    for reference_file in input_nanopore :
        fichier.write("\n   Mapping on "+ reference_file)
        for element in input_illumina :
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
                os.system("breseq -r " + output_path + "/" + reference_file +"/Prokka/" + reference_file + ".gbk "
                          + " -n " + element + "_mapped "
                          + input_path + "/" + element + "/" + new_dir[0] + " "
                          + input_path + "/" + element + "/" + new_dir[1]
                          + " -o " + output_path + "/" + element + "/Breseq/" + reference_file)
            else :
                if len(new_dir) == 1 :
                    os.system("breseq -r " + output_path + "/" + reference_file +"/Prokka/" + reference_file + ".gbk "
                              + " -n " + element + "_mapped "
                              + input_path + "/" + element + "/concatenated.fastq "
                              + " -o " + output_path + "/" + element + "/Breseq/" + reference_file)
                else :
                    pre_dir = os.listdir(input_path + "/" + element)
                    dir = []
                    for file in pre_dir:
                        if "concatenated" in file:
                            dir.append(file)
                    os.system("breseq -r " + output_path + "/" + reference_file + "/Prokka/" + reference_file + ".gbk "
                              + " -n " + element + "_mapped "
                              + input_path + "/" + element + "/" + new_dir[0]
                              + " -o " + output_path + "/" + element + "/Breseq/" + reference_file)
            fichier.write("\n       " + element + " has been mapped on " + reference_file)
    return



# illumina + reference_gene
def breseq_illumina_reference_gene(input_path, output_path, input_illumina, input_ref_gene, fichier) : #list_illumina     reference_gene
    for reference_file in input_ref_gene[0] :
        fichier.write("\n   Mapping on " + reference_file)
        for element in input_illumina :
            dir = os.listdir(input_path + "/" + element)
            new_dir = []
            for subfile in dir :
                if "fastq" in subfile :
                    new_dir.append(subfile)
            if "_R1" in new_dir[0] or "_R2" in new_dir[0]  :
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
                os.system("breseq -r " + output_path + "/" + reference_file +"/Prokka/" + reference_file + ".gbk "
                          + " -n " + element + "_mapped "
                          + input_path + "/" + element + "/" + new_dir[0] + " "
                          + input_path + "/" + element + "/" + new_dir[1]
                          + " -o " + output_path + "/" + element + "/Breseq/" + reference_file)
            else :
                if len(new_dir) == 1 :
                    os.system("breseq -r " + output_path + "/" + reference_file +"/Prokka/" + reference_file + ".gbk "
                              + " -n " + element + "_mapped "
                              + input_path + "/" + element + "/concatenated.fastq "
                              + " -o " + output_path + "/" + element + "/Breseq/" + reference_file)
                else :
                    pre_dir = os.listdir(input_path + "/" + element)
                    dir = []
                    for file in pre_dir:
                        if "concatenated" in file:
                            dir.append(file)
                    os.system("breseq -r " + output_path + "/" + reference_file + "/Prokka/" + reference_file + ".gbk "
                              + " -n " + element + "_mapped "
                              + input_path + "/" + element + "/" + new_dir[0]
                              + " -o " + output_path + "/" + element + "/Breseq/" + reference_file)
                fichier.write("\n       " + element + " has been mapped on " + reference_file)
    for reference_file in input_ref_gene[1] :
        fichier.write("\n   Mapping on " + reference_file)
        for element in input_illumina :
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
                os.system("breseq -r " + output_path + "/" + reference_file +"/Prokka/" + reference_file + ".gbk "
                          + " -n " + element + "_mapped "
                          + input_path + "/" + element + "/" + new_dir[0] + " "
                          + input_path + "/" + element + "/" + new_dir[1] + " "
                          + " -o " + output_path + "/" + element + "/Breseq/" + reference_file)
            else :
                if len(new_dir) == 1 :
                    os.system("breseq -r " + output_path + "/" + reference_file +"/Prokka/" + reference_file + ".gbk "
                              + " -n " + element + "_mapped "
                              + input_path + "/" + element + "/concatenated.fastq "
                              + " -o " + output_path + "/" + element + "/Breseq/" + reference_file)
                else :
                    pre_dir = os.listdir(input_path + "/" + element)
                    dir = []
                    for file in pre_dir:
                        if "concatenated" in file:
                            dir.append(file)
                    os.system("breseq -r " + output_path + "/" + reference_file + "/Prokka/" + reference_file + ".gbk "
                              + " -n " + element + "_mapped "
                              + input_path + "/" + element + "/" + new_dir[0]
                              + " -o " + output_path + "/" + element + "/Breseq/" + reference_file)
            fichier.write("\n       " + element + " has been mapped on " + reference_file)
    for reference_file in input_ref_gene[2] :
        fichier.write("\n   Mapping on " + reference_file)
        for element in input_illumina :
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
                os.system("breseq -r " + output_path + "/" + reference_file +"/Prokka/" + reference_file + ".gbk "
                          + " -n " + element + "_mapped "
                          + input_path + "/" + element + "/" +  new_dir[0] + " "
                          + input_path + "/" + element + "/" +  new_dir[1]
                          + " -o " + output_path + "/" + element + "/Breseq/" + reference_file)
            else :
                if len(new_dir) == 1 :
                    os.system("breseq -r " + output_path + "/" + reference_file +"/Prokka/" + reference_file + ".gbk "
                              + " -n " + element + "_mapped "
                              + input_path + "/" + element + "/concatenated.fastq "
                              + " -o " + output_path + "/" + element + "/Breseq/" + reference_file)
                else :
                    pre_dir = os.listdir(input_path + "/" + element)
                    dir = []
                    for file in pre_dir:
                        if "concatenated" in file:
                            dir.append(file)
                    os.system("breseq -r " + output_path + "/" + reference_file + "/Prokka/" + reference_file + ".gbk "
                              + " -n " + element + "_mapped " + " "
                              + input_path + "/" + element + "/" + new_dir[0]
                              + " -o " + output_path + "/" + element + "/Breseq/" + reference_file)
            fichier.write("\n       " + element + " has been mapped on " + reference_file)
    return


# Illumina + reference_prot
def breseq_illumina_reference_prot(input_path, output_path, input_illumina, input_ref_prot, fichier) : # list_illumina  reference_prot
    for reference_file in input_ref_prot :
        fichier.write("\n   Mapping on " + reference_file)
        for element in input_illumina :
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
                os.system("breseq -r " + input_path + "/" + reference_file
                          + " -n " + element + "_mapped "
                          + input_path + "/" + element + "/" + new_dir[0] + " "
                          + input_path + "/" + element + "/" + new_dir[1]
                          + " -o " + output_path + "/" + element + "/Breseq/" + reference_file)
            else:
                if len(new_dir) == 1:
                    os.system("breseq -r " + input_path + "/" + reference_file
                              + " -n " + element + "_mapped "
                              + input_path + "/" + element + "/concatenated.fastq "
                              + " -o " + output_path + "/" + element + "/Breseq/" + reference_file)
                else:
                    pre_dir = os.listdir(input_path + "/" + element)
                    dir = []
                    for file in pre_dir:
                        if "concatenated" in file:
                            dir.append(file)
                    os.system("breseq -r " + input_path + "/" + reference_file
                              + " -n " + element + "_mapped "
                              + input_path + "/" + element + "/" + new_dir[0]
                              + " -o " + output_path + "/" + element + "/Breseq/" + reference_file)
            fichier.write("\n       " + element + " has been mapped on " + reference_file)
    return


### Using Gdtool compare to compare all the samples' mutations

def gdtools_compare_genic (input_path, output_path, reference_genes, comparing_genes):
    for reference_file in reference_genes :
        input_phrase = "gdtools COMPARE -f tsv -o " + output_path + "/" + reference_file + "/comparison.tsv " + " -r " + output_path + "/" + reference_file + "/Prokka/" + reference_file + ".gbk "
        input_phrase2 = "gdtools COMPARE -f html -o " + output_path + "/" + reference_file + "/comparison.html " + " -r " + output_path + "/" + reference_file + "/Prokka/" + reference_file + ".gbk "
        for element in comparing_genes:
            input_phrase += output_path + "/" + element + "/Breseq/" + reference_file + "/output/output.gd "
            input_phrase2 += output_path + "/" + element + "/Breseq/" + reference_file + "/output/output.gd "
        os.system(input_phrase)
        os.system(input_phrase2)
        return

def gdtools_compare_proteic(input_path, output_path, references_gbk, comparing_genes):
    for reference_file in references_gbk:
        input_phrase = "gdtools COMPARE -f tsv -o " + output_path + "/" + reference_file + "/comparison.tsv " + " -r " + input_path + "/" + reference_file + " "
        input_phrase2 = "gdtools COMPARE -f html -o " + output_path + "/" + reference_file + "/comparison.html " + " -r " + input_path + "/" + reference_file + " "
        for element in comparing_genes:
            input_phrase += output_path + "/" + element + "/Breseq/" + reference_file + "/output/output.gd "
            input_phrase2 += output_path + "/" + element + "/Breseq/" + reference_file + "/output/output.gd "
        os.system(input_phrase)
        os.system(input_phrase2)
    return

def gdtools_compare_total (input_path, output_path, reference_genes, references_gbk, comparing_genes):
    gdtools_compare_genic(input_path, output_path, reference_genes, comparing_genes)
    gdtools_compare_proteic(input_path, output_path, references_gbk, comparing_genes)
    return

'''
Valid output formats:
  HTML    Descriptive table viewable in a web browser
  GD      GenomeDiff with added annotation of mutations
  TSV     Tab-separated values file suitable for input into R or Excel
  PHYLIP  Alignment file suitable for input into PHYLIP
  JSON  	JavaScript object notation file suitable for parsing'''

def extraction_gdtools(file) :
    compare_gdtool_file = open(file, "r")
    compare_gdtool = compare_gdtool_file.read()
    phrases = compare_gdtool.split('\n')
    new_phrases = []
    for phrase in phrases :
        cuted_phrase = phrase.split('	')
        new_phrases.append(cuted_phrase)
    return new_phrases[0:-1]


### Using Breseq for the phylogenetic Tree (phylip document)

def gdtools_phylip_genic (input_path, output_path, reference_genes, comparing_genes):
    for reference_file in reference_genes :
        input_phrase = "breseq gdtools COMPARE -f phylip -o " + output_path + "/" + reference_file + "/comparison.phylip " + " -r " + output_path + "/" + reference_file + "/Prokka/" + reference_file + ".gbk "
        for element in comparing_genes:
            input_phrase += output_path + "/" + element + "/Breseq/" + reference_file + "/output/output.gd "
        os.system(input_phrase)
        return

def gdtools_phylip_proteic(input_path, output_path, references_gbk, comparing_genes):
    for reference_file in references_gbk:
        input_phrase = "breseq gdtools COMPARE -f phylip -o " + output_path + "/" + reference_file + "/comparison.phylip " + " -r " + input_path + "/" + reference_file + " "
        for element in comparing_genes:
            input_phrase += output_path + "/" + element + "/Breseq/" + reference_file + "/output/output.gd "
        os.system(input_phrase)
    return

def gdtools_phylip_total (input_path, output_path, reference_genes, references_gbk, comparing_genes):
    gdtools_phylip_genic(input_path, output_path, reference_genes, comparing_genes)
    gdtools_phylip_proteic(input_path, output_path, references_gbk, comparing_genes)
    return


