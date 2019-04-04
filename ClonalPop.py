#!/usr/bin/env python # -*-coding:Latin-1 -*
import os
import sys
import assembling_sequences
import annotating_sequences
import mapping_sequences
import qualite_sequencage_unicycler as qsu
import qualite_sequence as qsc
import abricate
import gestion_excel


### Activation and options
input_path = sys.argv[1]
output_path = sys.argv[2]

os.makedirs(output_path, exist_ok=True)

### Creating the simplified log
fichier = open(output_path + "/simplified_log.txt", "w")
fichier.write("Pipeline started")


### Getting the files
files_and_folders = os.listdir(input_path)

list_illumina = []
list_nanopore = []
list_iontorrent = []
reference_gene = []
reference_prot = []
mini_list_0 = []
mini_list_1 = []
mini_list_2 = []
for element in files_and_folders :
    if "illumina" in element :
        list_illumina.append(element)
    elif "nanopore" in element :
        list_nanopore.append(element)
    elif "reference" in element :
        if 'fasta' or 'fas' or 'fa' in element :
            if 'fasta' in element :
                mini_list_0.append(element[0:-6])
            elif 'fas' in element :
                mini_list_1.append(element[0:-4])
            elif 'fa' in element :
                mini_list_2.append(element[0:-3])
        if '.gbk' in element :
            reference_prot.append(element)
        elif '.gb' in element :
            reference_prot.append(element + "k")
reference_gene.append(mini_list_0)
reference_gene.append(mini_list_1)
reference_gene.append(mini_list_2)

list_illumina.sort()
list_nanopore.sort()
list_iontorrent.sort()
# no reference_gene.sort() because the list contains 3 mini_lists : [ [fasta] [fas] [fa]]
reference_prot.sort()

all_samples = [list_illumina,
               list_nanopore,
               list_iontorrent,
               reference_gene,
               reference_prot]

#### Creating a sub_folder for each file in the final folder
for element in list_illumina :
    os.makedirs(output_path + "/" + element, exist_ok=True)
for element in list_nanopore :
    os.makedirs(output_path + "/" + element, exist_ok=True)
for element in list_iontorrent :
    os.makedirs(output_path + "/" + element, exist_ok=True)
for sub_list in reference_gene :
    for element in sub_list :
        os.makedirs(output_path + "/" + element, exist_ok=True)
for element in reference_prot :
    os.makedirs(output_path + "/" + element + "_prot", exist_ok=True)

all_samples_extension_recup = os.listdir(output_path)

list_illumina_extension = []
list_nanopore_extension = []
list_iontorrent_extension = []
list_reference_gene_extension = []
list_reference_prot_extension = []

for element in all_samples_extension_recup:
    if "illumina" in element :
        list_illumina_extension.append(element)
    elif "nanopore" in element :
        list_nanopore_extension.append(element)
    elif "iontorrent" in element :
        list_iontorrent_extension.append(element)
    elif "reference_gene" in element :
        list_reference_gene_extension.append(element)
    elif "reference_gene" in element :
        list_reference_prot_extension.append(element)

list_illumina_extension.sort()
list_nanopore_extension.sort()
list_iontorrent_extension.sort()
list_reference_gene_extension.sort()
list_reference_prot_extension.sort()
all_samples_extension = [list_illumina_extension,
                         list_nanopore_extension,
                         list_iontorrent_extension,
                         list_reference_gene_extension,
                         list_reference_prot_extension]

# Merging all the lists in two : reference and sample
acces_dossier_compare = [element for element in list_illumina]
acces_dossier_reference = [element for element in list_nanopore]
for sub_list in reference_gene :
    for element in sub_list :
        acces_dossier_reference.append(element)
for element in reference_prot:
    acces_dossier_reference.append(element)


# Unzipping all the illumina files
fichier.write("\n\nUnzipping of the potential .gz files")
for element in list_illumina :
    dir = os.listdir(input_path + "/" + element)
    for sub_element in dir :
        if ".gz" in sub_element:
            os.system("gunzip " + input_path + "/" + element + "/" + sub_element)
            fichier.write("\n   " + element + "has been unzipped")
fichier.write("\nAll the files are unzipped\n")

# Turning _fastq into .fastq
for element in list_illumina :
    dir = os.listdir(input_path + "/" + element)
    for sub_element in dir :
        if "_fastq" in sub_element:
            stringl = input_path + "/" + element + "/" + sub_element
            stringl2 = stringl[0:-6] + ".fastq"
            os.rename(stringl, stringl2)

### Assembling the files
fichier.write("\nStarting assembly")
## Illumina
fichier.write("\n   Assembling the Illumina sequences (SPAdes)")
assembling_sequences.assembler_spades_illumina(input_path, output_path, list_illumina, fichier)
fichier.write("\n   All the Illumina sequences have been assembled")
## Nanopore
fichier.write("\n   Assembling the Nanopore sequences (Unicycler)")
assembling_sequences.assembler_nanopore(input_path, output_path, list_nanopore, fichier)
fichier.write("\n   All the Nanopore sequences have been assembled")
fichier.write("\nEnding assembly\n")

### Annotating the files
fichier.write("\nStarting annotation (Prokka)")
fichier.write("\n   Annotating the Illumina sequences")
annotating_sequences.prokka (all_samples_extension[0], input_path, output_path, fichier)
fichier.write("\n   All the Illumina sequences have been annotated")
annotating_sequences.prokka(all_samples_extension[2], input_path, output_path, fichier)
fichier.write("\n   Annotating the Nanopore sequences")
annotating_sequences.prokka_unicycler(list_nanopore, input_path, output_path, fichier)
fichier.write("\n   All the Nanopore sequences have been annotated")
fichier.write("\n   Annotating the reference sequences")
annotating_sequences.prokka_reference (all_samples[3], input_path, output_path, fichier)
fichier.write("\n   All the reference sequences have been annotated")
fichier.write("\nEnding annotation\n")


### Gestion of the reference sequences
reference_sequences = []
for element in all_samples_extension_recup :
    if "nanopore" or "reference" in element :
        reference_sequences.append(element)
reference_sequences.sort()
del reference_sequences[reference_sequences.index('simplified_log.txt')]


### Mapping the sequences
fichier.write("\nStarting mapping (Breseq)")
reference_non_prokka = reference_prot
mapping_sequences.breseq(input_path, output_path, list_nanopore,list_illumina, fichier)
mapping_sequences.breseq_illumina_reference_gene(input_path, output_path,
                                                 list_illumina, reference_gene, fichier)
mapping_sequences.breseq_illumina_reference_prot(input_path, output_path,
                                                 list_illumina, reference_prot, fichier)
fichier.write("\nEnding mapping\n")

# Making the comparison.tsv file
fichier.write("\nComparing the genomes")
mapping_sequences.gdtools_compare_total (input_path, output_path,
                                         acces_dossier_reference, reference_prot,
                                         acces_dossier_compare)
reference_number = len(acces_dossier_reference) + len(reference_prot)
fichier.write("\nCreating the .phylip alignment\n")
# Makin the Phylip file
mapping_sequences.gdtools_phylip_total (input_path,
                                        output_path,
                                        acces_dossier_reference,
                                        reference_prot,
                                        acces_dossier_compare)

'''
### Abricate
fichier.write("\nStarting content analysis (Abricate)")
for element in acces_dossier_compare:
    fichier.write("\n   Analysis on " + element)
    abricate.abricate(output_path, element)
protein_string_AntibioRes_CARD, \
protein_string_Virulence_ECVF, \
protein_string_Virulence_VFDB, \
protein_string_Plasmids_PlasmidFinder, \
protein_string_AntimicRes_ResFinder, \
protein_string_AntimicRes_NCBI \
    = abricate.extracting_data_for_protein_index\
    (acces_dossier_compare, output_path)
fichier.write("\nEnding content analysis\n")
'''

### Making the Excel file
fichier.write("\nCreating the Excel output\n")
import xlsxwriter
import string
letters_excel = list(string.ascii_uppercase)
letters_excel.extend([i+b for i in letters_excel for b in letters_excel])

wb = xlsxwriter.Workbook(output_path + "/Output3_wesh.xlsx")

style = wb.add_format({'bold': True, 'font_name': 'Calibri',
                       'font_size' : 17.5, 'bottom':True, 'top': True,
                       'left': True, 'right': True, 'bg_color': '#2FCCFF', 'text_wrap': True})
style_database = wb.add_format({'bold': True,
                                'font_name': 'Calibri', 'font_size' : 15,
                                'bottom':True, 'top': True, 'left': True,
                                'right': True, 'bg_color': '#3EFEFF',
                                'text_wrap': True})
style_gene_index = wb.add_format({'bold': False,
                                  'font_name': 'Calibri',
                                  'font_size' : 12.5,
                                  'text_wrap': True})
style_cells = wb.add_format({'bold': False, 'font_name': 'Calibri',
                             'font_size' : 12.5, 'text_wrap': True})
style_cell_headline = wb.add_format({'bold': True, 'font_name': 'Calibri',
                                     'font_size' : 12.5, 'bottom':True,
                                     'top': True, 'left': True,
                                     'right': True, 'bg_color': '#2FCCFF',
                                     'text_wrap': True})
style_citation_pipeline = wb.add_format({'bold': True, 'font_name': 'Calibri',
                                         'font_size' : 12.5, 'bottom':True, 'top': True,
                                         'left': True, 'right': True,
                                         'bg_color': '#FFFF01', 'text_wrap': True})
style.set_align('center')
style.set_align('vcenter')
style_database.set_align('vcenter')
style_database.set_align('center')
style_gene_index.set_align('vcenter')
style_cells.set_align('center')
style_cells.set_align('vcenter')
style_cell_headline.set_align('center')
style_cell_headline.set_align('vcenter')
style_citation_pipeline.set_align('center')
style_citation_pipeline.set_align('vcenter')


## Defining the xl quality file for Illumina sequences
feuil1 = wb.add_worksheet('Illumina quality')
gestion_excel.illumina_sheet_head(feuil1, style)
gestion_excel.illumina_sheet_corps(output_path,
                                   acces_dossier_compare,
                                   feuil1,
                                   style_cells,
                                   letters_excel)


## Defining the xl quality file for Nanopore (reference) sequences
feuil1 = wb.add_worksheet('Nanopore quality')
gestion_excel.nanopore_sheet_head(feuil1, style)
gestion_excel.nanopore_sheet_corps(output_path, list_nanopore, feuil1, style_cells, letters_excel)


## Mutational comparison (Breseq)
gestion_excel.display_comparison(output_path, acces_dossier_reference, wb, style, style_gene_index)

for reference_file in acces_dossier_reference:
    assembling_sequences.assembler_spades_unmatched_illumina(output_path,
                                                             acces_dossier_compare,
                                                             reference_file)

'''
### Abricate general
feuil4 = wb.add_worksheet('Genic overview')
gestion_excel.abricate_overview_head(feuil4, style)
gestion_excel.abricate_overview_corps(feuil4, style_cells, acces_dossier_compare, output_path)


### Abricate across samples
protein_string_AntibioRes_CARD, \
protein_string_Virulence_ECVF, \
protein_string_Virulence_VFDB, \
protein_string_Plasmids_PlasmidFinder, \
protein_string_AntimicRes_ResFinder, \
protein_string_AntimicRes_NCBI \
    = abricate.extracting_data_for_protein_index(acces_dossier_compare, output_path)
gestion_excel.abricate_across_sample_excel(output_path, feuil4, style_cell_headline, style_cells, wb)

### Gene index
feuilx = wb.add_worksheet('Gene index')
protein_list_AntibioRes_CARD, protein_list_Virulence_ECVF, \
protein_list_Virulence_VFDB, protein_list_Plasmids_PlasmidFinder, \
protein_list_AntimicRes_ResFinder,protein_list_AntimicRes_NCBI \
    = abricate.extracting_data_for_protein_index_2(output_path, acces_dossier_compare)

gestion_excel.gene_index(feuilx,
                         protein_list_AntibioRes_CARD,
                         protein_list_Virulence_ECVF,
                         protein_list_Virulence_VFDB,
                         protein_list_Plasmids_PlasmidFinder,
                         protein_list_AntimicRes_ResFinder,
                         protein_list_AntimicRes_NCBI, style,
                         style_cells, style_database, style_gene_index, letters_excel)
'''

### Phylogenetic tree(s)
feuilx = wb.add_worksheet('Phylo_tree')
feuilx.set_row(0, 60)
for reference_file in acces_dossier_reference:
    os.system("python mini_phylo_tree.py "
              + output_path + " " + reference_file
              + " > " + output_path + "/" + reference_file + "/phylo_tree.txt")


gestion_excel.display_tree(output_path,
                           acces_dossier_reference,
                           style, style_database,
                           style_gene_index,
                           feuilx, letters_excel)


### Citations
feuilx = wb.add_worksheet('Citations')
gestion_excel.display_citations(feuilx,
                                style_gene_index,
                                style_citation_pipeline,
                                style_database, style)


### Suppressing every byproduct
fichier.write("\nSuppressing unwanted byproducts\n")
for element in acces_dossier_compare:
    if 'spades' not in sys.argv:
        fichier.write("\n   Removing SPAdes byproducts")
        os.remove(output_path + "/" + element + "/Spades")
    if 'prokka' not in sys.argv:
        fichier.write("\n   Removing Prokka byproducts")
        os.remove(output_path + "/" + element + "/Prokka")
    if 'breseq' not in sys.argv:
        fichier.write("\n   Removing Breseq byproducts")
        os.remove(output_path + "/" + element + "/Breseq")
    #if 'abricate' not in sys.argv:
    #    fichier.write("\n   Removing ABRicate byproducts")
    #    os.remove(output_path + "/" + element + "/Abricate")
    if 'quality' not in sys.argv:
        fichier.write("\n   Removing Quality Graphs byproducts")
        os.remove(output_path + "/" + element + "/Quality_check")

wb.close()
'''
# a la fin runner abricate sur les fasta unmatched et afficher un g�nic overview des g�nesnon mapp�s

for reference_file in acces_dossier_reference :
    for element in acces_dossier_compare:
        print(output_path + "/" + element + "/Breseq/" + reference_file + "/data/" + element + ".unmatched.fastq")
        mini_dir = os.listdir(output_path + "/" + element + "/Breseq/" + reference_file + "/data")
        assembling_sequences.assembler_spades_unmatched_illumina(output_path, acces_dossier_compare, reference_file)
'''
fichier.write("End of the pipeline")
fichier.close()