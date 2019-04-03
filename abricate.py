# -*-coding:Latin-1 -*
import os


def abricate(output_path, element) :
    '''os.makedirs(output_path + "/" + element + "/Abricate", exist_ok=True)'''
    os.system("docker run replikation/abricate --db card " + output_path + "/" + element + "/Prokka/" + element + ".fna > " + output_path + "/" + element +  "/Abricate/" +  element + "_AntibioRes_CARD.txt")
    os.system("abricate --db resfinder " + output_path + "/" + element + "/Prokka/" + element + ".fna > " + output_path + "/" + element +  "/Abricate/" +  element + "_AntimicRes_ResFinder.txt")
    os.system("abricate --db ncbi "  + output_path + "/" + element + "/Prokka/" + element + ".fna > " + output_path + "/" + element +  "/Abricate/" +  element + "_AntimicRes_NCBI.txt")
    os.system("abricate --db ecoli_vf "  + output_path + "/" + element + "/Prokka/" + element + ".fna > " + output_path + "/" + element +  "/Abricate/" +  element +  "_Virulence_ECVF.txt")
    os.system("abricate --db vfdb "  + output_path + "/" + element + "/Prokka/" + element + ".fna > " + output_path + "/" + element +  "/Abricate/" +  element +  "_Virulence_VFDB.txt")
    os.system("abricate --db plasmidfinder "  + output_path + "/" + element + "/Prokka/" + element + ".fna > " + output_path + "/" + element +  "/Abricate/" +  element + "_Plasmids_PlasmidFinder.txt")


#abricate(letter_illumina)
# Faire apparaitre les gines puis faireune matrice d'analyse presence absence
# 1 gros tableau bourrin avec toutes les infos comme prevu au debut
# Plusieurs slides apres av
def bacteria_resistance(file) :
    fichier = open(file, "r")
    contenu = fichier.read()
    phrases = []
    list_protein_result = []
    for phrase in contenu.split('\n') :
        phrases.append(phrase)
    new_phrases = phrases [13:-1]
    result = ""
    gene_result =""
    protein_result = []
    for phrase in new_phrases :
        #print(phrase)
        mini_phrase = phrase.split()
        #print(mini_phrase)
        #new_phrases = new_phrases [13:]
        localisation = mini_phrase[1]
        gene = mini_phrase[4]
        coverage = float(mini_phrase[8])
        identity = float(mini_phrase[9])
        trust_coefficient = str(round(coverage*identity/10000,2))
        product = ""
        for k in range(12, len(mini_phrase)):
            product += mini_phrase[k] + " " #str(mini_phrase[12:])
        #mini_result = "Gene: "+ gene + ", Protein: " + product + " (" + trust_coefficient + ") \n "
        mini_protein_result = "Gene: " + gene + ", Protein: " + product + " \n "
        #mini_gene_result = [ gene +  "(" + trust_coefficient + " \n "]
        mini_gene_result = gene + "\n "
        gene_result += mini_gene_result
        mini_list_protein_result = [gene, product]
        list_protein_result.append(mini_list_protein_result)
        protein_result.append(mini_protein_result)
        #gene_result.append(mini_gene_result)
        #print(list_protein_result)
    fichier.close()
    return (gene_result, protein_result, list_protein_result)


def bacteria_virulence_ECVF(file) :
    fichier = open(file, "r")
    contenu = fichier.read()
    phrases = []
    list_protein_result = []
    for phrase in contenu.split('\n') :
        phrases.append(phrase)
    new_phrases = phrases [13:-1]
    result = ""
    gene_result = ''#[]
    protein_result = []
    for phrase in new_phrases :
        #print(phrase)
        mini_phrase = phrase.split()
        #print(mini_phrase)
        #new_phrases = new_phrases [13:]
        localisation = mini_phrase[1]
        gene = mini_phrase[4]
        coverage = float(mini_phrase[8])
        identity = float(mini_phrase[9])
        trust_coefficient = str(round(coverage*identity/10000,2))
        product = ""
        for k in range(12, len(mini_phrase)):
            product += mini_phrase[k] + " " #str(mini_phrase[12:])
        #mini_result = "Gene: "+ gene + ", Protein: " + product + " (" + trust_coefficient + ") \n "
        mini_protein_result = "Gene: " + gene + ", Protein: " + product + " \n "
        mini_gene_result = gene + "\n "
        gene_result += mini_gene_result
        mini_list_protein_result = [gene, product]
        list_protein_result.append(mini_list_protein_result)
        protein_result.append(mini_protein_result)
        #gene_result.append(mini_gene_result)
        #print(list_protein_result)
    fichier.close()
    return (gene_result, protein_result, list_protein_result)


def bacteria_virulence_VDFB(file) :
    fichier = open(file, "r")
    contenu = fichier.read()
    phrases = []
    list_protein_result= []
    for phrase in contenu.split('\n') :
        phrases.append(phrase)
    new_phrases = phrases [13:-1]
    result = ""
    gene_result =""
    protein_result = []
    for phrase in new_phrases :
        #print(phrase)
        mini_phrase = phrase.split()
        #print(mini_phrase)
        #new_phrases = new_phrases [13:]
        localisation = mini_phrase[1]
        gene = mini_phrase[4]
        coverage = float(mini_phrase[8])
        identity = float(mini_phrase[9])
        trust_coefficient = str(round(coverage*identity/10000,2))
        product = ""
        for k in range(12, len(mini_phrase)):
            product += mini_phrase[k] + " " #str(mini_phrase[12:])
        #mini_result = "Gene: "+ gene + ", Protein: " + product + " (" + trust_coefficient + ") \n "
        mini_protein_result = "Gene: " + gene + ", Protein: " + product + " \n "
        mini_gene_result = gene + "\n "
        gene_result += mini_gene_result
        mini_list_protein_result = [gene, product]
        list_protein_result.append(mini_list_protein_result)
        protein_result.append(mini_protein_result)
        #gene_result.append(mini_gene_result)
        #print(list_protein_result)
    fichier.close()
    return (gene_result, protein_result, list_protein_result)


def bacteria_PlasmidFinder(file) :
    fichier = open(file, "r")
    contenu = fichier.read()
    phrases = []
    for phrase in contenu.split('\n') :
        element = phrase.split('\t')
        for mini_element in element :
            phrases.append(mini_element)
    #print(phrases)
    new_phrases = phrases [13:-1]
    n_phrases = int(len(new_phrases)/13)
    result = ""
    gene_result = ""
    protein_result = []
    list_protein_result= []
    #print(new_phrases)
    for k in range(n_phrases) :
        mini_phrase = new_phrases[0:13]
        new_phrases = new_phrases[13:]
        #print(mini_phrase)
        #new_phrases = new_phrases [13:]
        localisation = mini_phrase[1]
        gene = mini_phrase[4]
        coverage = float(mini_phrase[8])
        identity = float(mini_phrase[9])
        trust_coefficient = str(round(coverage*identity/10000,2))
        #product_1 = mini_phrase[12]
        #product_2 = mini_phrase[13:]
        #product = str(product_1) + "("
        #for mini_product in product_2 :
        #    product += str(mini_product) + " "
        #product+= ")"
        product = ""
        for k in range(12, len(mini_phrase)):
            product += mini_phrase[k] + " " #str(mini_phrase[12:])
        mini_protein_result = "Gene: " + gene + ", Protein: " + product + " \n "
        mini_gene_result = gene + "\n "
        gene_result += mini_gene_result
        mini_list_protein_result = [gene, product]
        list_protein_result.append(mini_list_protein_result)
        protein_result.append(mini_protein_result)
        #gene_result.append(mini_gene_result)
        #print(list_protein_result)
    fichier.close()
    return (gene_result, protein_result, list_protein_result)


def bacteria_AntimicRes_ResFinder(file) :
    fichier = open(file, "r")
    contenu = fichier.read()
    phrases = []
    for phrase in contenu.split('\n') :
        element = phrase.split('\t')
        for mini_element in element :
            phrases.append(mini_element)
    #print(phrases)
    new_phrases = phrases [13:-1]
    n_phrases = int(len(new_phrases)/13)
    result = ""
    gene_result = ""
    protein_result = []
    list_protein_result =[]
    #print(new_phrases)
    for k in range(n_phrases) :
        mini_phrase = new_phrases[0:13]
        new_phrases = new_phrases[13:]
        #print(mini_phrase)
        #new_phrases = new_phrases [13:]
        localisation = mini_phrase[1]
        gene = mini_phrase[4]
        coverage = float(mini_phrase[8])
        identity = float(mini_phrase[9])
        trust_coefficient = str(round(coverage*identity/10000,2))
        #product_1 = mini_phrase[12]
        #product_2 = mini_phrase[13:]
        #product = str(product_1) + "("
        #for mini_product in product_2 :
        #    product += str(mini_product) + " "
        #product+= ")"
        product = ""
        for k in range(12, len(mini_phrase)):
            product += mini_phrase[k] + " " #str(mini_phrase[12:])
        mini_protein_result = "Gene: " + gene + ", Protein: " + product + " \n "
        mini_gene_result = gene + "\n "
        gene_result += mini_gene_result
        mini_list_protein_result = [gene, product]
        list_protein_result.append(mini_list_protein_result)
        protein_result.append(mini_protein_result)
        #gene_result.append(mini_gene_result)
        #print(list_protein_result)
    fichier.close()
    return (gene_result, protein_result, list_protein_result)


def bacteria_AntimicRes_NCBI(file) :
    fichier = open(file, "r")
    contenu = fichier.read()
    phrases = []
    for phrase in contenu.split('\n') :
        element = phrase.split('\t')
        for mini_element in element :
            phrases.append(mini_element)
    #print(phrases)
    new_phrases = phrases [13:-1]
    n_phrases = int(len(new_phrases)/13)
    result = ""
    gene_result = ""
    protein_result = []
    list_protein_result = []
    #print(new_phrases)
    for k in range(n_phrases) :
        mini_phrase = new_phrases[0:13]
        new_phrases = new_phrases[13:]
        #print(mini_phrase)
        #new_phrases = new_phrases [13:]
        localisation = mini_phrase[1]
        gene = mini_phrase[4]
        coverage = float(mini_phrase[8])
        identity = float(mini_phrase[9])
        trust_coefficient = str(round(coverage*identity/10000,2))
        #product_1 = mini_phrase[12]
        #product_2 = mini_phrase[13:]
        #product = str(product_1) + "("
        #for mini_product in product_2 :
        #    product += str(mini_product) + " "
        #product+= ")"
        product = ""
        for k in range(12, len(mini_phrase)):
            product += mini_phrase[k] + " " #str(mini_phrase[12:])
        mini_protein_result = "Gene: " + gene + ", Protein: " + product + " \n "
        mini_gene_result = gene + "\n "
        gene_result += mini_gene_result
        mini_list_protein_result = [gene, product]
        list_protein_result.append(mini_list_protein_result)
        protein_result.append(mini_protein_result)
        #gene_result.append(mini_gene_result)
        #print(list_protein_result)
    fichier.close()
    return (gene_result, protein_result, list_protein_result)

## REGLER CA AUSSI
def extracting_everything_abricate(output_path, acces_dossier_compare) :
    gene_list = [['resistance','virulence_ECVF','virulence_VDFB','PlasmidFinder','AntimicRes_ResFinder','AntimicRes_NCBI']]
    protein_list = [['resistance', 'virulence_ECVF', 'virulence_VDFB', 'PlasmidFinder', 'AntimicRes_ResFinder','AntimicRes_NCBI']]
    for element in acces_dossier_compare:
        mini_gene_list = []
        output_path + "/" + element + "/Abricate/" + element
        mini_gene_list.append(bacteria_resistance(output_path + "/" + element + "/Abricate/" + element + "_AntibioRes_CARD.txt")[0])
        mini_gene_list.append(bacteria_virulence_ECVF(output_path + "/" + element + "/Abricate/" + element + "_Virulence_ECVF.txt")[0])
        mini_gene_list.append(bacteria_virulence_VDFB(output_path + "/" + element + "/Abricate/" + element + "_Virulence_VFDB.txt")[0])
        mini_gene_list.append(bacteria_PlasmidFinder(output_path + "/" + element + "/Abricate/" + element + "_Plasmids_PlasmidFinder.txt")[0])
        mini_gene_list.append(bacteria_AntimicRes_ResFinder(output_path + "/" + element + "/Abricate/" + element + "_AntimicRes_ResFinder.txt")[0])
        mini_gene_list.append(bacteria_AntimicRes_NCBI(output_path + "/" + element + "/Abricate/" + element + "_AntimicRes_NCBI.txt")[0])
        gene_list.append(mini_gene_list)
    return gene_list

'''
# Feuille generale
#extracting_everything_abricate()

total_abricate = abricate.extracting_everything_abricate()
total_abricate = total_abricate[1:]

for element in letter_illumina:
    sample_ID = ID + element
    feuil3.write(k, 0, sample_ID, style_cells)
    for element in total_abricate:
        print(element)
        #feuil3.write(k, 1, element[0], style_cells)
    k+=1

'''

#print(extracting_everything_abricate())
# Liste avec tous les echantillons, dans chaque sous-liste 6 listes (on garde les 6)
# Liste interessante : list_protein_result, chaque element = une liste de deux elements (on garde la 2 (troisieme)
    # [gene, product] peut etre garder juste gene et ajouter '\n' à chaque fois (on garde le 0 (premier))


### Extraire des donnes pour lindex des proteines (termine)
def extracting_data_for_protein_index(acces_dossier_compare, output_path):
    gene_list = []
    protein_list = []
    for element in acces_dossier_compare:
        mini_protein_list = []
        mini_protein_list.append(bacteria_resistance(output_path + "/" + element + "/Abricate/" + element +"_AntibioRes_CARD.txt")[1])
        mini_protein_list.append(
            bacteria_resistance(output_path + "/" + element + "/Abricate/" + element + "_Virulence_ECVF.txt")[1])
        mini_protein_list.append(
            bacteria_resistance(output_path + "/" + element + "/Abricate/" + element + "_Virulence_VFDB.txt")[1])
        mini_protein_list.append(
            bacteria_resistance(output_path + "/" + element + "/Abricate/" + element + "_Plasmids_PlasmidFinder.txt")[1])
        mini_protein_list.append(
            bacteria_resistance(output_path + "/" + element + "/Abricate/" + element + "_AntimicRes_ResFinder.txt")[1])
        mini_protein_list.append(
            bacteria_resistance(output_path + "/" + element + "/Abricate/" + element + "_AntimicRes_NCBI.txt")[1])
        #mini_protein_list.append(bacteria_resistance("/Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(
        #    child) + "/" + sample_ID + "/Abricate/" + sample_ID + "_AntibioRes_CARD.txt")[1])
        #mini_protein_list.append(bacteria_virulence_ECVF(
        #    "/Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(
        #        child) + "/" + sample_ID + "/Abricate/" + sample_ID + "_Virulence_ECVF.txt")[1])
        #mini_protein_list.append(bacteria_virulence_VDFB(
        #    "/Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(
        #        child) + "/" + sample_ID + "/Abricate/" + sample_ID + "_Virulence_VFDB.txt")[1])
        #mini_protein_list.append(bacteria_PlasmidFinder(
        #    "/Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(
        #        child) + "/" + sample_ID + "/Abricate/" + sample_ID + "_Plasmids_PlasmidFinder.txt")[1])
        #mini_protein_list.append(bacteria_AntimicRes_ResFinder(
        #    "/Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(
        #        child) + "/" + sample_ID + "/Abricate/" + sample_ID + "_AntimicRes_ResFinder.txt")[1])
        #mini_protein_list.append(bacteria_AntimicRes_NCBI(
        #    "/Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(
        #        child) + "/" + sample_ID + "/Abricate/" + sample_ID + "_AntimicRes_NCBI.txt")[1])
        protein_list.append(mini_protein_list)
    #print(protein_list)
    protein_list_AntibioRes_CARD = []
    protein_list_Virulence_ECVF = []
    protein_list_Virulence_VFDB = []
    protein_list_Plasmids_PlasmidFinder = []
    protein_list_AntimicRes_ResFinder = []
    protein_list_AntimicRes_NCBI = []
    for element in protein_list :
        protein_list_AntibioRes_CARD.append(element[0])
        protein_list_Virulence_ECVF.append(element[1])
        protein_list_Virulence_VFDB.append(element[2])
        protein_list_Plasmids_PlasmidFinder.append(element[3])
        protein_list_AntimicRes_ResFinder.append(element[4])
        protein_list_AntimicRes_NCBI.append(element[5])
    #print(protein_list_AntibioRes_CARD)
    protein_string_AntibioRes_CARD =""
    protein_string_Virulence_ECVF = ""
    protein_string_Virulence_VFDB = ""
    protein_string_Plasmids_PlasmidFinder = ""
    protein_string_AntimicRes_ResFinder = ""
    protein_string_AntimicRes_NCBI = ""
    for element in protein_list_AntibioRes_CARD :
        for sous_element in element :
            if sous_element not in protein_string_AntibioRes_CARD :
                protein_string_AntibioRes_CARD += sous_element
    for element in protein_list_Virulence_ECVF :
        for sous_element in element :
            if sous_element not in protein_string_Virulence_ECVF :
                protein_string_Virulence_ECVF += sous_element
    for element in protein_list_Virulence_VFDB :
        for sous_element in element :
            if sous_element not in protein_string_Virulence_VFDB :
                protein_string_Virulence_VFDB += sous_element
    for element in protein_list_Plasmids_PlasmidFinder :
        for sous_element in element :
            if sous_element not in protein_string_Plasmids_PlasmidFinder :
                protein_string_Plasmids_PlasmidFinder += sous_element
    for element in protein_list_AntimicRes_ResFinder :
        for sous_element in element :
            if sous_element not in protein_string_AntimicRes_ResFinder :
                protein_string_AntimicRes_ResFinder += sous_element
    for element in protein_list_AntimicRes_NCBI :
        for sous_element in element :
            if sous_element not in protein_string_AntimicRes_NCBI :
                protein_string_AntimicRes_NCBI += sous_element
    return protein_string_AntibioRes_CARD, protein_string_Virulence_ECVF, protein_string_Virulence_VFDB, protein_string_Plasmids_PlasmidFinder, protein_string_AntimicRes_ResFinder, protein_string_AntimicRes_NCBI
#print(extracting_data_for_protein_index()[2])


def extracting_data_for_protein_index_2(output_path, acces_dossier_compare) :
    gene_list = []
    protein_list = []
    for element in acces_dossier_compare:
        mini_protein_list = []
        mini_protein_list.append(bacteria_resistance(output_path + "/" + element + "/Abricate/" + element + "_AntibioRes_CARD.txt")[2])
        mini_protein_list.append(bacteria_virulence_ECVF(output_path + "/" + element + "/Abricate/" + element + "_Virulence_ECVF.txt")[2])
        mini_protein_list.append(bacteria_virulence_VDFB(
            output_path + "/" + element + "/Abricate/" + element + "_Virulence_VFDB.txt")[2])
        mini_protein_list.append(bacteria_PlasmidFinder(
            output_path + "/" + element + "/Abricate/" + element + "_Plasmids_PlasmidFinder.txt")[2])
        mini_protein_list.append(bacteria_AntimicRes_ResFinder(
            output_path + "/" + element + "/Abricate/" + element + "_AntimicRes_ResFinder.txt")[2])
        mini_protein_list.append(bacteria_AntimicRes_NCBI(
            output_path + "/" + element + "/Abricate/" + element + "_AntimicRes_NCBI.txt")[2])
        protein_list.append(mini_protein_list)
    #print(protein_list)
    protein_list_AntibioRes_CARD = []
    protein_list_Virulence_ECVF = []
    protein_list_Virulence_VFDB = []
    protein_list_Plasmids_PlasmidFinder = []
    protein_list_AntimicRes_ResFinder = []
    protein_list_AntimicRes_NCBI = []
    for element in protein_list:
        for sous_element in element[0] :
            if sous_element not in protein_list_AntibioRes_CARD :
                protein_list_AntibioRes_CARD.append(sous_element)
        for sous_element in element[1]:
            if sous_element not in protein_list_Virulence_ECVF :
                protein_list_Virulence_ECVF.append(sous_element)
        for sous_element in element[2]:
            if sous_element not in protein_list_Virulence_VFDB :
                protein_list_Virulence_VFDB.append(sous_element )
        for sous_element in element[3]:
            if sous_element not in protein_list_Plasmids_PlasmidFinder:
                protein_list_Plasmids_PlasmidFinder.append(sous_element )
        for sous_element in element[4]:
            if sous_element not in protein_list_AntimicRes_ResFinder :
                protein_list_AntimicRes_ResFinder.append(sous_element)
        for sous_element in element[5]:
            if sous_element not in protein_list_AntimicRes_NCBI:
                protein_list_AntimicRes_NCBI.append(sous_element )
    return protein_list_AntibioRes_CARD , protein_list_Virulence_ECVF , protein_list_Virulence_VFDB, protein_list_Plasmids_PlasmidFinder,protein_list_AntimicRes_ResFinder,protein_list_AntimicRes_NCBI
'''
protein_list_AntibioRes_CARD , protein_list_Virulence_ECVF , protein_list_Virulence_VFDB, protein_list_Plasmids_PlasmidFinder,protein_list_AntimicRes_ResFinder,protein_list_AntimicRes_NCBI = extracting_data_for_protein_index_2()



### Combining reports across samples

def abricate_report_across_samples(letter_illumina) :
#    for element in letter_illumina:
#        sample_ID = ID + element
#        try:
#            os.mkdir("/Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/" + sample_ID + "/Abricate/")
#        except OSError:
#            pass
#        os.system("abricate --db card /Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/" + sample_ID + "/Prokka/" + sample_ID + "_illumina_prokka.fna > /Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/" + sample_ID + "/Abricate/" +  sample_ID + "_AntibioRes_CARD.tab")
#        os.system("abricate --db resfinder /Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/" + sample_ID + "/Prokka/" + sample_ID + "_illumina_prokka.fna > /Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/" + sample_ID + "/Abricate/" + sample_ID + "_AntimicRes_ResFinder.tab")
#        os.system("abricate --db ncbi /Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/" + sample_ID + "/Prokka/" + sample_ID + "_illumina_prokka.fna > /Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/" + sample_ID + "/Abricate/" + sample_ID + "_AntimicRes_NCBI.tab")
#        os.system("abricate --db ecoli_vf /Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/" + sample_ID + "/Prokka/" + sample_ID + "_illumina_prokka.fna > /Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/" + sample_ID + "/Abricate/" + sample_ID + "_Virulence_ECVF.tab")
#        os.system("abricate --db vfdb /Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/" + sample_ID + "/Prokka/" + sample_ID + "_illumina_prokka.fna > /Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/" + sample_ID + "/Abricate/" + sample_ID + "_Virulence_VFDB.tab")
#        os.system("abricate --db plasmidfinder /Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/" + sample_ID + "/Prokka/" + sample_ID + "_illumina_prokka.fna > /Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/" + sample_ID + "/Abricate/" + sample_ID + "_Plasmids_PlasmidFinder.tab")
    sentence_AntibioRes_CARD =''
    sentence_AntimicRes_ResFinder =''
    sentence_AntimicRes_NCBI =''
    sentence_Virulence_ECVF =''
    sentence_Virulence_VFDB = ''
    sentence_Plasmids_PlasmidFinder =''
    for element in letter_illumina:
        sample_ID = ID + element
        sentence_AntibioRes_CARD += "/Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/" + sample_ID + "/Abricate/" +  sample_ID + "_AntibioRes_CARD.tab "
        sentence_AntimicRes_ResFinder += "/Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/" + sample_ID + "/Abricate/" +  sample_ID + "_AntimicRes_ResFinder.tab "
        sentence_AntimicRes_NCBI += "/Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/" + sample_ID + "/Abricate/" +  sample_ID + "_AntimicRes_NCBI.tab "
        sentence_Virulence_ECVF += "/Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/" + sample_ID + "/Abricate/" +  sample_ID + "_Virulence_ECVF.tab "
        sentence_Virulence_VFDB  += "/Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/" + sample_ID + "/Abricate/" +  sample_ID + "_Virulence_VFDB.tab "
        sentence_Plasmids_PlasmidFinder += "/Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/" + sample_ID + "/Abricate/" +  sample_ID + "_Plasmids_PlasmidFinder.tab "
    try:
        os.mkdir("/Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/Abricate_Reports_across_samples")
    except OSError:
        pass
    os.system("abricate --summary " + sentence_AntibioRes_CARD +" > /Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/Abricate_Reports_across_samples/AntibioRes_CARD_report_samples.txt")
    os.system("abricate --summary " + sentence_AntimicRes_ResFinder + " > /Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/Abricate_Reports_across_samples/AntimicRes_ResFinder_report_samples.txt")
    os.system("abricate --summary " + sentence_AntimicRes_NCBI + " > /Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/Abricate_Reports_across_samples/AntimicRes_NCBI_report_samples.txt")
    os.system("abricate --summary " + sentence_Virulence_ECVF + " > /Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/Abricate_Reports_across_samples/Virulence_ECVF_report_samples.txt")
    os.system("abricate --summary " + sentence_Virulence_VFDB + " > /Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/Abricate_Reports_across_samples/Virulence_VFDB_report_samples.txt")
    os.system("abricate --summary " + sentence_Plasmids_PlasmidFinder + " > /Users/Yanis/Desktop/Projet_de_recherche/Genomes/enfant_" + str(child) + "/Abricate_Reports_across_samples/Plasmids_PlasmidFinder_report_samples.txt")

abricate_report_across_samples(letter_illumina)'''


# Extracting phrases report samples

def extracting_report_samples(file):
    fichier = open(file, "r")
    contenu = fichier.read()
    phrases = []
    for phrase in contenu.split('\n') :
        phrases.append(phrase)
    #print(phrases)
    new_phrases = []
    for element in phrases :
        new_phrases.append(element.split('\t'))
    n = len(new_phrases[0])
    return(new_phrases, n)