# -*-coding:Latin-1 -*
### Importations
import re
import numpy as np
import matplotlib.pyplot as plt

### Extracting data

# Number of contigs
def nb_contig(contenu):
    nb_contigs = contenu.count(">")
    return nb_contigs


# File --> list of interesting phrases (those containing the infos)
def phrases_interessantes(contenu):
    phrases_interessantes = []
    for phrase in contenu.split():
        if "NODE" in phrase:
            phrases_interessantes.append(phrase)
    return (phrases_interessantes)

# Extraction of data keeping only the contigs of more than 1000bp)
def extraction(contenu):
    NODE = []
    length = []
    coverage = []
    for k in phrases_interessantes(contenu):  # no range because we want strings
        b = re.split("_", k) # Splitting phrase in each _
        if int(b[3]) > 1000 :
            NODE.append(int(b[1]))
            length.append(int(b[3]))
            coverage.append(round(float(b[5]), 0))  # round to the closest int
    contig = [i + 1 for i in range(len(coverage))]
    return contig, NODE, length, coverage

def nb_contig_sup_1000(contig):
    nb_contigs = len(contig)
    return nb_contigs

### Analyzing data
'''We have :
- nb_contig : number of contig in the sequence
- contig : list of numbers from 1 to nb_contig
- NODE : all the NODEs IDs
- length : length of all the contigs in the sequence
- coverage : coverage of all the contigs in the sequence
'''

# Sequencing quality
def quali_nb_contig(contenu, contig):
    nb_contigs = nb_contig(contenu)
    nb_contig_sup_1000_v = nb_contig_sup_1000(contig)
    if 0 < nb_contigs < 80:
        r = "Excellent sequencing quality."
    elif 80 < nb_contigs < 170:
        r = "Good sequencing quality."
    else:
        r ="Medium sequencing quality"
    phrase = str(nb_contigs) + " contigs. \n" + r + "\n "+ str(nb_contig_sup_1000_v) + " contigs > 1kb"
    return phrase

# Genome size
def genome_size(length):
    new_length = [float(element) for element in length]
    genome_size = sum(length)
    if genome_size < 4000000 or genome_size > 6000000:
        r =("Genome size not satisfying")
    else:
        r = ("Satisfying genome size")
    phrase = str(round(genome_size/1000000,2)) + " Mb \n" + r
    return phrase

# Median contig N50
def N_50f(length,NODE):
    sorted_length = sorted(length)
    compteur_contig = 0
    index = 0
    while compteur_contig < sum(sorted_length) / 2:
        compteur_contig += sorted_length[index]
        index += 1
    index_N50 = index - 1
    contig_N50 = sorted_length[index_N50]
    index_N50 = length.index(contig_N50)
    NODE_N50 = NODE[index_N50]
    return index_N50, contig_N50, NODE_N50

# Coverage
def hist_cover_f(coverage, output_path, element):
    hist_cover = plt.figure(figsize=(6,5), facecolor='w', edgecolor='k')
    plt.hist(coverage, bins=np.linspace(0, 250, 25 + 1))
    plt.savefig(output_path + "/" + element + "/Quality_check/hist_cover")
    plt.close()
    return(hist_cover)

# Length of the contigs
def graph_length(NODE, length, output_path, element):
    graph_length = plt.figure(figsize=(6,5))
    plt.plot(NODE, length, 'g', linewidth=0.8, marker='')
    plt.savefig(output_path + "/" + element + "/Quality_check/graph_length")
    plt.close()
    return(graph_length)

# Coverage of the contigs
def graph_cover_f(NODE, coverage, output_path, element):
    graph_cover = plt.figure(figsize=(6,5))
    plt.plot(NODE, coverage, "b", linewidth=0.8, marker='')
    plt.savefig(output_path + "/" + element + "/Quality_check/graph_cover")
    plt.close()
    return (graph_cover)

# Length repartition
def hist_contig_f(length,NODE, output_path, element):
    hist_contig = plt.figure(figsize=(6,5))
    plt.hist(length, bins=np.linspace(0, max(length), 100 + 1))
    plt.axvline(x= N_50f(length,NODE)[1], color='red')
    plt.savefig(output_path + "/" + element + "/Quality_check/hist_contig")
    plt.close()
    return(hist_contig)

# Coverage by length
def length_coverage_f(length, coverage, output_path, element):
    length_coverage_v = plt.figure(figsize=(6,5))
    plt.plot(length, coverage, 'g', linewidth=0.0, marker='.')
    plt.loglog()
    plt.savefig(output_path + "/" + element + "/Quality_check/length_coverage")
    plt.close()
    return (length_coverage_v)
