# -*-coding:Latin-1 -*
import re
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt

### Extracting data

## File --> interesting lines
def phrases_interessantes(contenu):
    phrases_interessantes = []
    for phrase in contenu.split():
        if ">" in phrase:
            phrases_interessantes.append(phrase)
    return (phrases_interessantes)

### Assessing the sequence' quality

## On the number of contigs

def nb_contig(contenu):
    nb_contigs = contenu.count(">")
    return nb_contigs

def quali_nb_contig(contenu, nb_contigs):
    nb_contigs = nb_contig(contenu)
    if 0 < nb_contigs < 80:
        r = "Excellent sequencing quality."
    elif 80 < nb_contigs < 170:
        r = "Good sequencing quality."
    else:
        r ="Medium sequencing quality"
    phrase = str(nb_contigs) + " contigs. \n" + r
    return phrase

## On the circularity of the contigs

def is_circular(phrases, contenu) :
    nb_contigs = contenu.count(">")
    nb_circular = contenu.count("circular=true")
    if nb_circular == nb_contigs :
            if nb_contigs == 1 :
                phrase = "The contig is circular."
            else :
                phrase = "All of the " + str(nb_contigs) + " contigs are circular."
            return phrase
    else :
        if nb_contig == 1 :
            phrase = "The contig is not circular"
            return phrase
        else :
            phrase_1 = "All the contigs are not circular"
            circular_list = re.findall(r"circular=(?P<circular>\w+)", contenu)
            circular_ones = []
            non_circular_ones = []
            for k in range(len(circular_list)) :
                if circular_list[k] == "true" :
                    circular_ones.append(k)
                else :
                    non_circular_ones.append(k)
            circular_ones = str(circular_ones)
            limit_circular = len(circular_ones) -1
            non_circular_ones = str(non_circular_ones)
            limit_non_circular = len(non_circular_ones) -1
            phrase_2 = "\nCircular contigs : " + str(circular_ones[1:limit_circular]) +"\nNon circular contigs : "+ str(non_circular_ones[1:limit_non_circular])
            phrase = phrase_1 + phrase_2
            return phrase

## On the size of the genome

def genome_size(contenu):
    length_list = re.findall(r"length=(?P<circular>\d+)", contenu)
    length_list = [int(i) for i in length_list]
    genome_size = sum([int(i) for i in length_list])
    if genome_size < 4000000 or genome_size > 6000000:
        r =("Genome size not satisfying")
    else:
        r = ("Satisfying genome size")
    phrase = str(round(genome_size/1000000,2)) + " Mb \n" + r
    return phrase, length_list, genome_size

## Getting the median contig N_50

def N_50f(length_list):
    sorted_length = sorted(length_list)
    compteur_contig = 0
    index = 0
    while compteur_contig < sum(length_list) / 2:
        compteur_contig += sorted_length[index]
        index += 1
    index_N50 = index - 1
    contig_N50 = sorted_length[index_N50]
    index_N50 = length_list.index(contig_N50) + 1
    return index_N50, contig_N50

### Graphical analysis

## Depth

def depth(contenu, output_path, element):
    depth_list = re.findall(r"depth=(?P<circular>....)", contenu)
    depth_list = [float(element) for element in depth_list]
    len_depth = len(depth_list)
    contig = [int(i + 1) for i in range(len_depth)]
    if len(depth_list) <= 5:
        phrase_depth = ""
        for element in contig:
            phrase_depth += "\nContig " + str(element) + " : " + str(depth_list[element - 1])
        return phrase_depth
    else:
        graph_depth = plt.figure(figsize=(6,5))
        plt.plot(contig, depth_list, 'g', marker='', linewidth=1)
        plt.savefig(output_path + "/" + element + "/Quality_check/depth")
        plt.close()
        return (graph_depth)

## Length of contigs

# Plot
def length_contig(contenu, output_path, element):
    length_list = re.findall(r"length=(?P<circular>\d+)", contenu)
    length_list = [float(element) for element in length_list]
    len_depth = len(length_list)
    contig = [int(i+1) for i in range(len_depth)]
    if len(length_list) <= 5 :
        phrase_length = ""
        for element in contig :
            phrase_length += "\nContig " + str(element) + " : " + str(length_list[element-1])
        return phrase_length
    else :
        sequence_figure = plt.figure(figsize=(6,5))
        plt.plot(contig, length_list, 'g', marker='', linewidth=1)
        plt.savefig(output_path + "/" + element + "/Quality_check/length_contig")
        plt.close()
        return(sequence_figure)

# Histogram of length repartition

def hist_length_contig(length_list, output_path, element, contenu, contig_N50):
    nb_contigs = contenu.count(">")
    if nb_contigs == 1:
        return ("There is only one contig so no histogram to print")
    else:
        hist_length_contig_h = plt.figure(figsize=(6,5))
        plt.hist(length_list, bins=np.linspace(0, max(length_list) - 1, 10), normed='True')
        plt.axvline(x=contig_N50, color='red')
        plt.savefig(output_path + "/" + element + "/Quality_check/length_repartition")
        plt.close()
        return (hist_length_contig_h)