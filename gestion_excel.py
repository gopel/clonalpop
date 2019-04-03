# -*-coding:Latin-1 -*
import os
import qualite_sequence as qsc
import qualite_sequencage_unicycler as qsu
import mapping_sequences
import abricate
import string

### Defining the xl quality file for Illumina sequences

def illumina_sheet_head(feuil1, style):
    feuil1.set_default_row(100)
    feuil1.set_row(0, 60)
    feuil1.set_column(0, 8, 25)
    feuil1.write('A1', 'Sample', style)
    feuil1.write('B1', 'Sequencing quality', style)
    feuil1.write('C1', 'Genome size (Mb)', style)
    feuil1.write('D1', 'N_50', style)
    feuil1.write('E1', 'Coverage by contig', style)
    feuil1.write('F1', 'Coverage repartition', style)
    feuil1.write('G1', 'Length of contigs', style)
    feuil1.write('H1', 'Length repartition', style)
    feuil1.write('I1', 'Coverage by length', style)
    return


def illumina_sheet_corps(output_path, acces_dossier_compare, feuil1, style_cells, letters_excel):
    line = 1
    column = 0
    for element in acces_dossier_compare:
        os.makedirs(output_path + "/" + element + "/Quality_check", exist_ok=True)
        fichier = open(output_path + "/" + element + "/Spades/contigs.fasta")
        contenu = fichier.read()
        contig, NODE, length, coverage = qsc.extraction(contenu)
        phrases = qsc.phrases_interessantes(contenu)
        nb_contig = qsc.nb_contig(contenu)
        phrase_genome_size = qsc.genome_size(length)
        index_N50, contig_N50, NODE_N50 = qsc.N_50f(length, NODE)
        phrase_N50 = "NODE N_50 = " + str(NODE_N50) + ". \nLength N_50 = " + str(int(contig_N50 / 1000)) + "kb."
        feuil1.write(letters_excel[0] + str(line + 1), str(element), style_cells)
        feuil1.write(letters_excel[1] + str(line + 1), qsc.quali_nb_contig(contenu, contig), style_cells)
        feuil1.write(letters_excel[2] + str(line + 1), phrase_genome_size, style_cells)
        feuil1.write(letters_excel[3] + str(line + 1), phrase_N50, style_cells)
        qsc.hist_cover_f(coverage, output_path, element)
        qsc.graph_length(NODE, length, output_path, element)
        qsc.graph_cover_f(NODE, coverage, output_path, element)
        qsc.hist_contig_f(length, NODE, output_path, element)
        qsc.length_coverage_f(length, coverage, output_path, element)
        feuil1.insert_image(letters_excel[4] + str(line + 1), output_path + "/" + element + '/Quality_check/graph_cover.png',
                               {'x_offset': 2, 'y_offset': 2, 'x_scale': 0.3, 'y_scale': 0.3})
        feuil1.insert_image(letters_excel[5] + str(line + 1),
                            output_path + "/" + element + '/Quality_check/hist_cover.png',
                            {'x_offset': 2, 'y_offset': 2, 'x_scale': 0.3, 'y_scale': 0.3})
        feuil1.insert_image(letters_excel[6] + str(line + 1),
                            output_path + "/" + element + '/Quality_check/graph_length.png',
                            {'x_offset': 2, 'y_offset': 2, 'x_scale': 0.3, 'y_scale': 0.3})
        feuil1.insert_image(letters_excel[7] + str(line + 1),
                            output_path + "/" + element + '/Quality_check/hist_contig.png',
                            {'x_offset': 2, 'y_offset': 2, 'x_scale': 0.3, 'y_scale': 0.3})
        feuil1.insert_image(letters_excel[8] + str(line + 1),
                            output_path + "/" + element + '/Quality_check/length_coverage.png',
                            {'x_offset': 2, 'y_offset': 2, 'x_scale': 0.3, 'y_scale': 0.3})
        fichier.close()
        line += 1
    return


### Defining the xl quality file for Nanopore (reference) sequences
def nanopore_sheet_head(feuil1, style):
    feuil1.set_row(0, 60)
    feuil1.set_default_row(100)
    feuil1.set_column(0, 0, 20)
    feuil1.set_column(1, 8, 25)
    feuil1.write(0, 0, 'Sample', style)
    feuil1.write(0, 1, 'Sequencing quality', style)
    feuil1.write(0, 2, 'Genome size (Mb)', style)
    feuil1.write(0, 3, 'N_50', style)
    feuil1.write(0, 4, 'Genome Circularity', style)
    feuil1.write(0, 5, 'Genome Depth', style)
    feuil1.write(0, 6, 'Length of contigs', style)
    feuil1.write(0, 7, 'Length repartition', style)
    return


def nanopore_sheet_corps(output_path, list_nanopore, feuil1, style_cells, letters_excel):
    line = 1
    column = 0
    for element in list_nanopore:
        os.makedirs(output_path + "/" + element + "/Quality_check", exist_ok=True)
        fichier = open(output_path + "/" + element + "/Unicycler/assembly.fasta")
        contenu = fichier.read()
        phrases = qsu.phrases_interessantes(contenu)
        nb_contig = qsu.nb_contig(contenu)
        phrase_genome_size, length_list, total_length = qsu.genome_size(contenu)
        index_N50, contig_N50 = qsu.N_50f(length_list)
        phrase_N50 = "Index N_50 = " + str(index_N50) + ". \nLength N_50 = " + str(contig_N50)
        qsu.hist_length_contig(length_list, output_path, element, contenu, contig_N50)
        feuil1.write(letters_excel[0] + str(line + 1), str(element), style_cells)
        feuil1.write(letters_excel[1] + str(line + 1), qsu.quali_nb_contig(contenu, nb_contig), style_cells)
        feuil1.write(letters_excel[2] + str(line + 1), qsu.genome_size(contenu)[0], style_cells)
        feuil1.write(letters_excel[3] + str(line + 1), phrase_N50, style_cells)
        feuil1.write(letters_excel[4] + str(line + 1), qsu.is_circular(phrases, contenu), style_cells)
        if nb_contig <= 5:
            feuil1.write(letters_excel[5] + str(line + 1), qsu.depth(contenu, output_path, element), style_cells)
            feuil1.write(letters_excel[6] + str(line + 1), qsu.length_contig(contenu, output_path, element), style_cells)
        else:
            feuil1.insert_image(letters_excel[5] + str(line + 1),
                                output_path + "/" + element + "/Quality_check/depth.png",
                                {'x_offset': 2, 'y_offset': 2, 'x_scale': 0.3, 'y_scale': 0.3})
            feuil1.insert_image(letters_excel[6] + str(line + 1),
                                output_path + "/" + element + "/Quality_check/length_contig.png",
                                {'x_offset': 2, 'y_offset': 2, 'x_scale': 0.3, 'y_scale': 0.3})
        feuil1.insert_image(letters_excel[7] + str(line + 1),
                            output_path + "/" + element + "/Quality_check/length_repartition.png",
                            {'x_offset': 2, 'y_offset': 2, 'x_scale': 0.3, 'y_scale': 0.3})
        fichier.close()
        line += 1
    return


### Mutational comparison (Breseq)
def display_comparison(output_path, acces_dossier_reference, wb, style, style_gene_index):
    k=0
    for reference_file in acces_dossier_reference :
        k+=1

        feuil2 = wb.add_worksheet('Mutational comparison_' + str(k))
        feuil2.set_row(0, 60)

        for column in [1, 2, 32]:
            feuil2.set_column(column, column, 30)
        for column in [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16, 18, 22, 23, 24, 25, 26, 27, 28, 29, 31, 33, 34, 35]:
            feuil2.set_column(column, column, 10)
        for gene in [14, 15, 16, 30]:
            feuil2.set_column(column, column, 20)
        lists_compare = mapping_sequences.extraction_gdtools(output_path + "/" + reference_file + "/comparison.tsv")
        list2 = []
        numeros = [i for i in range(0,36)]
        for num_mini_list in range(len(lists_compare)):
            mini_list_2 = []
            for k in [9,10,11,12]:
                mini_list_2.append(lists_compare[num_mini_list][k])
            for k in numeros[0:9]:
                mini_list_2.append(lists_compare[num_mini_list][k])
            for k in numeros[13:]:
                mini_list_2.append(lists_compare[num_mini_list][k])
            list2.append(mini_list_2)
        lists_compare = list2
        for k in range(len(lists_compare)) : # k est la ligne
            indice = 0
            if k == 0 :
                for element in lists_compare[k]:
                    feuil2.write(k, indice, element, style)
                    indice+=1
            else :
                for element in lists_compare[k]:
                    feuil2.write(k, indice, element, style_gene_index)
                    indice+=1


### Abricate general
def abricate_overview_head(feuil4, style):
    feuil4.set_row(0, 60)
    feuil4.set_column(1, 1000, 30)
    feuil4.set_column(0, 0, 20)
    feuil4.write('A1', 'Sample', style)
    feuil4.write('B1', 'Antibiotic Resistance', style)
    feuil4.write('C1', 'Antimicrobial Resistance (ResFinder)', style)
    feuil4.write('D1', 'Antimicrobial Resistance (NCBI)', style)
    feuil4.write('E1', 'Virulence (general)', style)
    feuil4.write('F1', 'Virulence (E.coli)', style)
    feuil4.write('G1', 'Plasmids', style)
    return


def abricate_overview_corps(feuil4, style_cells, acces_dossier_compare, output_path, letters_excel):
    total_abricate = abricate.extracting_everything_abricate(output_path, acces_dossier_compare)
    total_abricate = total_abricate[1:]
    line = 1
    for compare in acces_dossier_compare:
        feuil4.write(letters_excel[0] + str(line + 1), compare, style_cells)
        # feuil3.row(k).height_mismatch = True
        # feuil3.row(k).height = 1000 * 20
        for element in total_abricate:
            feuil4.write(letters_excel[1] + str(line + 1), element[0], style_cells)  # element[0][0]
            feuil4.write(letters_excel[2] + str(line + 1), element[4], style_cells)
            feuil4.write(letters_excel[3] + str(line + 1), element[5], style_cells)
            feuil4.write(letters_excel[4] + str(line + 1), element[2], style_cells)
            feuil4.write(letters_excel[5] + str(line + 1), element[1], style_cells)
            feuil4.write(letters_excel[6] + str(line + 1), element[3], style_cells)
        line += 1
    return


### Abricate across samples

#protein_string_AntibioRes_CARD, protein_string_Virulence_ECVF, protein_string_Virulence_VFDB, protein_string_Plasmids_PlasmidFinder, protein_string_AntimicRes_ResFinder, protein_string_AntimicRes_NCBI = abricate.extracting_data_for_protein_index(acces_dossier_compare, output_path)

def displaying_abricate_across_samples(file, element, short_element, feuil4, style_cell_headline, style_cells, wb):
    feuil4 = wb.add_worksheet(short_element)
    feuil4.set_row(0, 60)
    feuil4.set_column(1, 1000, 10)
    feuil4.set_column(0, 0, 20)
    abricate_list = abricate.extracting_report_samples(file)
    abricate_list_0 = abricate_list[0][0:-1]
    abricate_list_1 = abricate_list[1]
    for k in range(len(abricate_list_0)):
        if k == 0:
            for element in abricate_list_0[k]:
                indice = abricate_list_0[k].index(element)
                feuil4.write(k, indice, element, style_cell_headline)
        else:
            line = abricate_list_0[k]
            line = line[1:]
            new_line = []
            new_line.append('XXX')     # acces_dossier_compare[k]
            for element in line:
                new_line.append(element)
            for r in range(len(new_line)):
                feuil4.write(k, r, new_line[r], style_cells)


def abricate_across_sample_excel(output_path, feuil4, style_cell_headline, style_cells, wb):
    abricate_compare = os.listdir(output_path + "/Abricate_Reports_across_samples")
    abricate_compare2 = [element[0:-19] for element in abricate_compare]
    indice = 0
    for element in abricate_compare :
        file = output_path + "/Abricate_Reports_across_samples/" + element
        short_element = abricate_compare2[indice]
        displaying_abricate_across_samples(file, element, short_element, feuil4, style_cell_headline, style_cells, wb)
        indice +=1


### Gene index
def gene_index(feuilx, protein_list_AntibioRes_CARD , protein_list_Virulence_ECVF , protein_list_Virulence_VFDB, protein_list_Plasmids_PlasmidFinder, protein_list_AntimicRes_ResFinder,protein_list_AntimicRes_NCBI, style, style_cells, style_database, style_gene_index, letters_excel):
    feuilx.set_row(0, 60)
    feuilx.set_column(1, 1, 170)
    feuilx.set_column(0, 0, 30)
    #feuilx.col(0).width = 7000

    k=0
    feuilx.write(letters_excel[0] + str(k + 1), 'Protein', style)
    feuilx.write(letters_excel[1] + str(k + 1), 'Characteristics', style)


    k+=3
    feuilx.write(letters_excel[0] + str(k + 1), 'Database : ResFinder' , style_database)
    k+=1
    for element in protein_list_AntimicRes_ResFinder :
        feuilx.write(letters_excel[0] + str(k + 1), element[0], style_cells)
        feuilx.write(letters_excel[1] + str(k + 1), element[1], style_gene_index)
        k+=1


    k+=3
    feuilx.write(letters_excel[0] + str(k + 1), 'Database : NCBI' , style_database)
    k+=1
    for element in  protein_list_AntimicRes_NCBI  :
        feuilx.write(letters_excel[0] + str(k + 1), element[0], style_cells)
        feuilx.write(letters_excel[1] + str(k + 1), element[1], style_gene_index)
        k+=1


    k+= 3
    feuilx.write(letters_excel[0] + str(k + 1), 'Database : CARD', style_database)
    k+=1
    for element in protein_list_AntibioRes_CARD :
        feuilx.write(letters_excel[0] + str(k + 1), element[0], style_cells)
        feuilx.write(letters_excel[1] + str(k + 1), element[1], style_gene_index)
        k+=1


    k+= 3
    feuilx.write(letters_excel[0] + str(k + 1), 'Database : VFDB', style_database)
    k+=1
    for element in protein_list_Virulence_VFDB :
        feuilx.write(letters_excel[0] + str(k + 1), element[0], style_cells)
        feuilx.write(letters_excel[1] + str(k + 1), element[1], style_gene_index)
        k+=1


    k+= 3
    feuilx.write(letters_excel[0] + str(k + 1), 'Database : ECVF', style_database)
    k+=1
    for element in protein_list_Virulence_ECVF :
        feuilx.write(letters_excel[0] + str(k + 1), element[0], style_cells)
        feuilx.write(letters_excel[1] + str(k + 1), element[1], style_gene_index)
        k+=1


    k+=3
    feuilx.write(letters_excel[0] + str(k + 1), 'Database : PlasmidFinder' , style_database)
    k+=1
    for element in protein_list_Plasmids_PlasmidFinder :
        feuilx.write(letters_excel[0] + str(k + 1), element[0], style_cells)
        feuilx.write(letters_excel[1] + str(k + 1), element[1], style_gene_index)
        k+=1
    return


### Phylogenetic tree(s)
def display_tree(output_path, acces_dossier_reference, style, style_database, style_gene_index, feuilx, letters_excel):
    feuilx.set_row(0, 60)
    feuilx.set_column(1, 1, 170)
    feuilx.set_column(0, 0, 30)
    feuilx.write('A1', 'Sample', style)
    feuilx.write('B1', 'Tree', style)
    line=1
    for reference_file in acces_dossier_reference:
        file = open(output_path + "/" + reference_file + "/phylo_tree.txt")
        tree = file.read()
        feuilx.write(letters_excel[0] + str(line + 1), reference_file, style_database)
        feuilx.write(letters_excel[1] + str(line + 1), tree, style_gene_index)
        line+=1


### Citations
def display_citations(feuilx, style_gene_index, style_citation_pipeline, style_database, style):
    feuilx.set_column(1, 1, 170)
    feuilx.set_column(0, 0, 30)

    feuilx.write(0, 1, 'XXX Phrase de demande de citation de l\'article pour cette pipeline', style_gene_index)

    feuilx.write(2, 0, 'This pipeline', style_citation_pipeline)
    feuilx.write(2, 1, 'Article de la pipeline', style_citation_pipeline)


    feuilx.write(4, 1, 'Here are the softwares used by this pipeline:', style_gene_index)

    feuilx.set_row(5, 30)

    feuilx.write(6, 0, 'CORE \nSOFTWARES', style)
    feuilx.write(7, 0, 'SPADES', style_database)
    feuilx.write(7, 1, '@article{spades,\n Author = {Anton Bankevich, Sergey Nurk, Dmitry Antipov, Alexey A. Gurevich, Mikhail Dvorkin, Alexander S. Kulikov, Valery M. Lesin, Sergey I. Nikolenko, Son Pham, Andrey D. Prjibelski, Alexey V. Pyshkin, Alexander V. Sirotkin, Nikolay Vyahhi, Glenn Tesler, Max A. Alekseyev, and Pavel A. Pevzner}, \nJournal = {Journal of computational biology}, \nTitle = {SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing}, \nYear = {2012}}', style_gene_index)

    feuilx.write(8, 0, 'BAYESHAMMER', style_database)
    feuilx.write(8, 1, '@article{bayeshammer, \nAuthor = {Sergey I. Nikolenko, Anton I. Korobeynikov, and Max A. Alekseyev}, \nJournal = {BMC Genomics},\nTitle = {BayesHammer: Bayesian clustering for error correction in single-cell sequencing},\nYear = {2012}}', style_gene_index)

    feuilx.write(9, 0, 'UNICYCLER', style_database)
    feuilx.write(9, 1, '@article{unicycler_assemble,\nAuthor = {Ryan R. Wick, Louise M. Judd, Claire L. Gorrie, Kathryn E. Holt},\nJournal = {PLOS Computational Biology},\nTitle = {Unicycler: Resolving bacterial genome assemblies from short and long sequencing reads},\nYear = {2017}}', style_gene_index)

    feuilx.write(10, 0, 'PROKKA', style_database)
    feuilx.write(10, 1, '@article{prokka,\nAuthor = {Seemann T.},\nJournal = {Bioinformatics},\nTitle = {Prokka: rapid prokaryotic genome annotation},\nYear = {2014}}', style_gene_index)

    feuilx.write(11, 0, 'PROKKA DATABASES', style_database)
    feuilx.write(11, 1, '@article{infernal,\nAuthor = {Kolbe DL and Eddy SR.},\nJournal = {Bioinformatics},\nTitle = {Fast filtering for RNA homology search},\nYear = {2011}} \n ~ \n@article{signalp4,\nAuthor = {Petersen TN, et al. },\nJournal = {Nature Methods},\nTitle = {SignalP 4.0: discriminating signal peptides from transmembrane regions},\nYear = {2011}} \n ~ \n@article{aragorn,\nAuthor = {Laslett D. and  Canback B. },\nJournal = {Nucleic Acids Research},\nTitle = {ARAGORN, a program to detect tRNA genes and tmRNA genes in nucleotide sequences},\nYear = {2004}}\n ~ \n@article{rnahammer,\nAuthor = {Lagesen K, et al.},\nJournal = {Nucleic Acids Research},\nTitle = {RNAmmer: consistent and rapid annotation of ribosomal RNA genes},\nYear = {2007}} \n ~ \n@article{prodigal,\nAuthor = {Hyatt D, et al.},\nJournal = {BMC Bioinformatics},\nTitle = {Prodigal: prokaryotic gene recognition and translation initiation site identification},\nYear = {2010}}', style_gene_index)

    feuilx.write(12, 0, 'BRESEQ', style_database)
    feuilx.write(12, 1,
                 '@article{breseq,\nAuthor = {Deatherage, D.E., Barrick, J.E.},\nJournal = {Methods in Molecular Biology},\nTitle = {Identification of mutations in laboratory-evolved microbes from next-generation sequencing data using breseq},\nYear = {2014}}',
                 style_gene_index)

    feuilx.write(13, 0, 'ABRICATE', style_database)
    feuilx.write(13, 1,
                 '@url{abricate,\nAuthor = {Torsten Seemann},\nDate-Added = {2019-01-03 01:21:48 +0100},\nDate-Modified = {2019-01-03 01:22:30 +0100},\nTitle = {ABRicate},\nUrldate = {https://github.com/tseemann/abricate}}',
                 style_gene_index)

    feuilx.set_row(14, 30)

    feuilx.write(15, 0, 'BACTERIAL\nANALYSES', style)
    feuilx.write(16, 0, 'RESFINDER', style_database)
    feuilx.write(16, 1, '@article{resfinder,\nAuthor = {Zankari E, Hasman H, Cosentino S, Vestergaard M, R"',
                 style_gene_index)

    feuilx.write(17, 0, 'NCBI\nANTIMICRES', style_database)
    feuilx.write(17, 1,
                 '@url{NCBIdatabase,\nAuthor = {NBCI},\nDate-Added = {2019-01-03 01:14:28 +0100},\nDate-Modified = {2019-01-03 01:15:00 +0100},\nTitle = {Bacterial Antimicrobial Resistance Reference Gene Database},\nUrldate = {https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313047}}',
                 style_gene_index)

    feuilx.write(18, 0, 'CARD DATABASE', style_database)
    feuilx.write(18, 1,
                 '@article{card,\nAuthor = {Jia B, Raphenya AR, Alcock B, Waglechner N, Guo P, Tsang KK, Lago BA, Dave BM, Pereira S, Sharma AN, Doshi S, Courtot M, Lo R, Williams LE, Frye JG, Elsayegh T, Sardar D, Westman EL, Pawlowski AC, Johnson TA, Brinkman FS, Wright GD, McArthur AG}, \nJournal = {Nucleic Acids Research}, \nTitle = {CARD 2017: expansion and model-centric curation of the comprehensive antibiotic resistance database}, \nYear = {2017}}',
                 style_gene_index)

    feuilx.write(19, 0, 'VFDB', style_database)
    feuilx.write(19, 1,
                 '@article{vfdb,\nAuthor = {Chen L, Yang J, Yu J, Yao Z, Sun L, Shen Y, Jin Q},\nDate-Added = {2019-01-03 01:08:35 +0100},\nDate-Modified = {2019-01-03 01:09:08 +0100},\nJournal = {Nucleic Acids Research},\nTitle = {VFDB: a reference database for bacterial virulence factors},\nYear = {2005}}',
                 style_gene_index)

    feuilx.write(20, 0, 'ECVF DATABASE', style_database)
    feuilx.write(20, 1,
                 '@url{ecvf,\nDate-Added = {2019-01-03 01:15:47 +0100},\nDate-Modified = {2019-01-03 01:16:52 +0100},\nTitle = {Escherichia coli virulence factors},\nUrldate = {https://github.com/phac-nml/ecoli_vf}}',
                 style_gene_index)

    feuilx.write(21, 0, 'PLASMIDFINDER', style_database)
    feuilx.write(21, 1,
                 '@article{plasmidfinder,\nAuthor = {Carattoli A, Zankari E, Garcia-Fernandez A, Voldby Larsen M, Lund O, Villa L, Aarestrup FM, Hasman H},\nJournal = {Antimicrobial Agents and Chemotherapy},\nTitle = {PlasmidFinder and pMLST: in silico detection and typing of plasmids},\nYear = {2014}}',
                 style_gene_index)