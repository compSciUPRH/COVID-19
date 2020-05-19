#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
COVIDsequenceAnalyss.py

Script that selects, analyzes DNA sequences using informative sequence sites
selected by their high entropy among the population. Based on a paper by 
Z. Zhao, B. A. Sokhansanj, and G. Rosen, “Characterizing geographical and
temporal dynamics of novel coronavirus sars-cov-2 using informative subtype
markers,” bioRxiv, 2020.

Created on Sat May 16 12:39:30 2020

@author: José O. Sotero Esteva
    Universidad de Puerto Rico en Humacao
    jose.sotero@upr.edu
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License version 3 as published by
    the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program (gpl.txt).  If not, see <http://www.gnu.org/licenses/>.

    Acknowledgements: The main funding source for this project has been provided
    by the UPR-Penn Partnership for Research and Education in Materials program, 
    USA National Science Foundation grant number NSF-DMR-1523463,. 

"""

import Bio.SeqIO
import random
import itertools
import math

# constants
SEQS_PER_PLACE  = 13  # the ammount of sequences from Puerto Riico at this moment
UMBRAL_ENTROPIA = 0.4

covid = Bio.SeqIO.parse("GISAID/gisaid_hcov-19_2020_05_16_17.fasta", "fasta")
#covid = Bio.SeqIO.parse("selectedSequences.fasta", "fasta")

sequences_by_country = {"USA/CA":[], "USA/NY":[], "Wuhan":[], "Italy":[], 
                   'USA/FL':[], 'Puerto Rico':[]}

#============================================================================
# read all sequences
#============================================================================
for virus in covid:
    for lugar in sequences_by_country:
        if lugar in virus.description and len(virus.seq) > 29000:
            # fix IDs
            virus.id = virus.description.split('|')[1]
            virus.description = virus.description.replace('Puerto Rico', 'Puerto_Rico')
            sequences_by_country[lugar].append(virus)
            break

#============================================================================
# reduce lists sizes to SEQS_PER_PLACE
#============================================================================
for lugar in sequences_by_country:
    n = len(sequences_by_country[lugar])
    for i in range(n-SEQS_PER_PLACE):
        del sequences_by_country[lugar][ random.randint(0,n-i-1)]


#============================================================================
# Align sequences
#============================================================================

secuencias = itertools.chain.from_iterable(sequences_by_country.values())
count = Bio.SeqIO.write(secuencias, "selectedSequences.fasta", "fasta")

''' Actually done manually in a server
from Bio.Align.Applications import ClustalwCommandline
clustalw_cline = ClustalwCommandline("clustalw-mpi", infile="selectedSequences.fasta", outfile='alignedSequences.clw')
stdout, stderr = clustalw_cline()
'''

alignments = Bio.AlignIO.read("alignedSequences.clw", "clustal")
#covid = list(Bio.SeqIO.parse("alignedSequences.fasta", "fasta"))

#============================================================================
# Compute bases frequencies (pk)
#============================================================================
P = []
for i in range(len(alignments[0].seq)):
    # conteo de bases
    pk = {}
    for virus in alignments:
        if virus.seq[i] in 'ACGT':
            if virus.seq[i] in pk:
                pk[virus.seq[i]] += 1
            else:
                pk[virus.seq[i]] = 1
            
    # cómpuito de frecuencias
    P.append({})
    for k in pk:
        P[-1][k] = pk[k] / sum(pk.values())

#============================================================================
# Masked entropy function
#============================================================================
H=[-sum(map(lambda x: x * math.log2(x), P[i].values())) for i in range(len(alignments[0].seq))]

#============================================================================
# draw entropy function
#============================================================================
import matplotlib.pyplot as plt
plt.plot(H, marker='o', markersize=3, color='k')
plt.title('Masked entropy function COVID-19')
plt.ylabel('Masked entropy')
plt.xlabel('Position')
plt.savefig("entropy.png")

#============================================================================
# Select Informative Subtype Markers
#============================================================================
ISM = []
for i in range(len(alignments[0].seq)):
    if H[i] >= UMBRAL_ENTROPIA:
        ISM.append(i)

#============================================================================
# print ISM sites table
#============================================================================
print("{:^6}  {:^20}  {:^13}".format('Site', 'Nucleotide Position', 'Entropy'))
for i in range(len(ISM)):
    print("{:^6}  {:^20}  {:^13.10}".format(i+1, ISM[i], H[ISM[i]]))

#============================================================================
# correct errors in ISMs
#============================================================================
import numpy as np
seq_ISM = {virus.id:np.array([virus.seq[i] for i in ISM]) for virus in alignments }
corrected_seq = seq_ISM.copy()
symbol = {'A':'A', 'C':'C', 'G':'G', 'T':'T', 'AT':'W', 'CG':'S', 'AC':'M',
          'GT':'K', 'AG':'R', 'CT':'Y', 'CGT':'B', 'AGT':'D', 'ACT':'H',
          'ACG':'V', 'ACGT':'N'}

prev_N = sum([sum(seq_ISM[sid] == 'N') for sid in seq_ISM])  # counts Ns
print("Amount of ambiguous bases N before error correction: ", prev_N)
for sid in seq_ISM:
    if 'N' in seq_ISM[sid]:
        # get identical ISMs
        identical = []
        for id_ref in seq_ISM:
            if (not 'N' in seq_ISM[id_ref]) and np.logical_or(seq_ISM[sid] == 'N', seq_ISM[sid] == seq_ISM[id_ref]).all():
                identical.append(seq_ISM[id_ref])
                #print(seq_ISM[id_ref])
        identical = np.array(identical)  # a matrix with ISMs as rows
        
        # creates string from unique bases in each position and gets corrsponding symbol
        corrected_seq[sid] = np.array([symbol[''.join(sorted(list(set(identical[:,i]))))] 
                                    for i in range(len(identical[0]))])
        corrected_seq[sid] = ''.join(corrected_seq[sid])
    else:
        corrected_seq[sid] = ''.join(seq_ISM[sid])

pos_N = sum([corrected_seq[sid].count('N') for sid in corrected_seq])  # counts Ns
print("Amount of ambiguous bases N after error correction: ", pos_N)

#============================================================================
# count unique ISMs
#============================================================================
unique_ISMs = set(corrected_seq.values())
ISM_counts = []
for ism in unique_ISMs:
    ISM_counts.append((list(corrected_seq.values()).count(ism), ism))
ISM_counts.sort(reverse=True)

#============================================================================
# plot bar chart of unique ISMs
#============================================================================
y, x = list(zip(*ISM_counts))
plt.clf()
plt.bar(range(len(y)), y)
plt.xticks(range(len(y)), x,rotation=90)
plt.ylabel('Counts')
plt.tight_layout () 
plt.savefig("ISMcount.png")

#============================================================================
# print descriptions of the secuences used
#============================================================================
covid = Bio.SeqIO.parse("GISAID/gisaid_hcov-19_2020_05_16_17.fasta", "fasta")
idToDescription = {v.description.split('|')[1]:v.description for v in covid}
print("ISM             Sequence")
places_ids = {}
places = list(sequences_by_country.keys())
for place in places:
    print("{:10}  &  ".format(place), end='')
    places_ids[place] = [virus.id for virus in sequences_by_country[place] ]
print()
for i in range(len(places_ids['Puerto_Rico'])):
    for place in places:
        print("{:10}  &  ".format(places_ids[place][i]), end='')
    print()
    '''
        for i in range(len(ISM)):
            print(virus.seq[ISM[i]], end='')
        print("  ", idToDescription[virus.id])
    '''

#============================================================================
# cout ISMs in each place
#============================================================================
isms_counts = {}
for place in sequences_by_country:
    place_ids = [virus.id for virus in sequences_by_country[place] ]
    place_isms_counts = {}
    for sid in corrected_seq:
        if sid in place_ids:
            if corrected_seq[sid] in place_isms_counts:
                place_isms_counts[corrected_seq[sid]] += 1
            else:
                place_isms_counts[corrected_seq[sid]]  = 1
    isms_counts[place] = place_isms_counts

#============================================================================
# plot pie charts
#============================================================================
from matplotlib import cm
viridis = cm.get_cmap('tab20', len(unique_ISMs))
vcolors = viridis(range(len(unique_ISMs)))
color_dict = {ism:vcolors[i] for ism, i in zip(unique_ISMs, range(len(unique_ISMs)))}

sp = 1
plt.clf()
plt.figure(num=None, figsize=(7, 6), dpi=80, facecolor='w', edgecolor='k')
for place in isms_counts:
    plt.subplot(3,2,sp)
    place_colors = [color_dict[ism] for ism in isms_counts[place]]
    plt.pie(isms_counts[place].values(), colors=place_colors,textprops=dict(color="w", weight="bold"))  # , autopct='%1.1f%%'
    plt.legend(list(isms_counts[place].keys()), loc='upper right', bbox_to_anchor=(2.3, 0.9), prop={'size': 9})
    plt.title(place)
    sp += 1
plt.subplots_adjust(wspace = 2.0)
plt.suptitle('COVID-19 subtypes distributons')
plt.savefig("COVIDpie.png")
