from Bio.Blast.Applications import NcbiblastpCommandline
from io import StringIO
#from StringIO import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from datetime import datetime
import time



#/Users/Yanis/Desktop/Projet_de_recherche/Output_folder/162106A_nanopore/Unicycler/assembly.fasta
#/Users/Yanis/Desktop/Projet_de_recherche/Output_folder/B_illumina/Spades/contigs.fasta

# # Create two sequence files
# seq1 = SeqRecord(Seq("FQTWEEFSRAAEKLYLADPMKVRVVLKYRHDGNLCIKVTDDLVCLVYRTDQAQDVKKIEKF"),
#                    id="seq1")
# seq2 =SeqRecord(Seq("FQTWEEFSRAEKLYLADPMKVRVVLRYRHVGNLCIKVTDDLICLVYRTDQAQDVKKIEKF"),
#                    id="seq2")
# SeqIO.write(seq1, "seq1.fasta", "fasta")
# SeqIO.write(seq2, "seq2.fasta", "fasta")

t1 = datetime.now().time()
start_time = time.time()


#blast_result_record = NCBIXML.read(StringIO("/Users/Yanis/Desktop/Projet_de_recherche/Output_folder/opuntia.xml"))
#print(blast_result_record)
#for alignment in blast_result_record.alignments:
#    for hsp in alignment.hsps:
#        print ('****Alignment****')
#        print ('sequence:', alignment.title)
#        print ('length:', alignment.length)
#        print ('e value:', hsp.expect)
#        print (hsp.query)
#        print (hsp.match)
#        print (hsp.sbjct)
#
#




#print(t1)

'''from Bio.Blast.Applications import NcbiblastxCommandline
blastx_cline = NcbiblastxCommandline(
    query="/Users/Yanis/Desktop/Projet_de_recherche/Output_folder/B_illumina/Spades_unmatched/162106A_nanopore/contigs.fasta",
    subject="/Users/Yanis/Desktop/Projet_de_recherche/Output_folder/162106A_nanopore/Unicycler/assembly.fasta",
    evalue=0.001,outfmt=5, out="/Users/Yanis/Desktop/Projet_de_recherche/Output_folder/opuntia_le_retour.txt")
blastx_cline
print(blastx_cline)
stdout, stderr = blastx_cline()
print(stdout, stderr)
print('Done')'''

'''from Bio import SearchIO
resultHandle = "/Users/Yanis/Desktop/Projet_de_recherche/Output_folder/opuntia_le_retour.xml"
blast_result = SearchIO.parse(resultHandle, 'blast-xml')
for qresult in blast_result:
    print("%s %s" % (qresult.id, qresult.description))'''


#print(blast_result[0][0])
#start = blast_result[0][0].hit_start
#end = blast_result[0][0].hit_end


handle = open("/Users/Yanis/Desktop/Projet_de_recherche/Output_folder/opuntia_le_retour.xml")
file = handle.read()

blast_result_record = NCBIXML.parse(StringIO(file))
#print(file)
#for element in NCBIXML.parse(file):
#    print(element)

'''for blast_record in blast_records:
for alignment in *blast_record.alignments*:
for hsp in alignment.hsps:
if hsp.expect < 0.04:
print '****Alignment****'
print 'sequence:', alignment.title
print 'length:', alignment.length
print 'e value:', hsp.expect
print hsp.query[0:75] + '...'
print hsp.match[0:75] + '...'
print hsp.sbjct[0:75] + '...'''''


for blast_record in blast_result_record:
    #print(blast_record)
    print(blast_record)
    print('Alignement', blast_record.alignments)
    for alignment in blast_record.alignments:
        print(alignment)
        for hsp in alignment.hsps:
            print(hsp)
            print ('****Alignment****')
            print ('sequence:', alignment.title)
            print ('length:', alignment.length)
            print ('e value:', hsp.expect)
            print (hsp.query)
            print (hsp.match)
            print (hsp.sbjct)

#t2 = datetime.now().time()
#print(t2)
print("--- %s seconds ---" % (time.time() - start_time))

''''# Run BLAST and parse the output as XML
output = NcbiblastpCommandline(query="/Users/Yanis/Desktop/Projet_de_recherche/Output_folder/162106A_nanopore/Unicycler/assembly.fasta", subject="/Users/Yanis/Desktop/Projet_de_recherche/Output_folder/B_illumina/Spades/contigs.fasta",
outfmt=5)()[0]
print('Done 2')
# query="opuntia.fasta", db="nr", evalue=0.001,
# ...                                      outfmt=5, out="opuntia.xml"
blast_result_record = NCBIXML.read(StringIO(output))
print(blast_result_record)

'''




'''### Test total :

from Bio.Blast.Applications import NcbiblastxCommandline
blastx_cline = NcbiblastxCommandline(
    query="/Users/Yanis/Desktop/Projet_de_recherche/Output_folder/A_illumina/Spades/contigs.fasta",
    subject="/Users/Yanis/Desktop/Projet_de_recherche/Output_folder/162106A_nanopore/Unicycler/assembly.fasta",
    evalue=0.001,outfmt=5, out="/Users/Yanis/Desktop/Projet_de_recherche/Output_folder/opuntia_le_retour_du_duel.xml")
blastx_cline
print(blastx_cline)
stdout, stderr = blastx_cline()
print(stdout, stderr)
print('Done')

handle = open("/Users/Yanis/Desktop/Projet_de_recherche/Output_folder/opuntia_le_retour_du_duel.xml")
file = handle.read()
blast_result_record = NCBIXML.parse(StringIO(file))

for blast_record in blast_result_record:
    #print(blast_record)
    print('Alignement',blast_record.alignments)
    for alignment in blast_record.alignments:
        print(alignment)
        for hsp in alignment.hsps:
            print ('****Alignment****')
            print ('sequence:', alignment.title)
            print ('length:', alignment.length)
            print ('e value:', hsp.expect)
            print (hsp.query)
            print (hsp.match)
            print (hsp.sbjct)'''