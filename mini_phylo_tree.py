#!/usr/bin/env python

import sys
#print(sys.argv)

loose = sys.argv
#print(loose)

output_path = loose[1]
reference_file = loose[2]

print(output_path)
print(reference_file)

#import re
#import os
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.Consensus import *
#from whole_pipeline_pu import output_path, acces_dossier_reference
#import matplotlib.pyplot as plt



aln = AlignIO.read(output_path + "/" + reference_file + "/comparison.phylip", 'phylip')
print(aln)
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)
constructor = DistanceTreeConstructor(calculator, 'nj')
tree = constructor.build_tree(aln)
tree.ladderize()  # Flip branches so deeper clades are displayed at top
Phylo.draw_ascii(tree)
