#!/usr/bin/env python3
# Yizhu Lin 02/21/2023
# for isoseq isoform & mutation analysis

import sys
import os
from isoseqUtil import *
import pandas as pd
import numpy as np

# input is geneid ENSGxxxx
gene = sys.argv[1]
mydir = sys.argv[2]
fadir = sys.argv[3]

# read input files file structure assume .sam, .intersect, .fa, in the same folder
sam = pd.read_table(os.path.join(mydir,'%s.sam' % gene),header=None, names=['read_id', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'SEQ', 'UMI'],
                    usecols=list(range(0,6))+[9, 11], comment='@')
bed = pd.read_table(os.path.join(mydir,'%s.intersect' % gene), header=None, usecols= [0,1,2,3,5,6,7,9,10,11,13,14,17,18,19,20],
                   names=['chrom','chromStart','chromEnd','read_id2','strand', 'thickStart', 'thickEnd', 'blockCount',
                         'blockSizes', 'blockStarts','annoStart', 'annoEnd', 'geneid', 'SYMBOL', 'ENTREZID', 'overlapLen'])
fa, geneStrand, geneStart, geneEnd = read_fa(os.path.join(fadir,'%s.fasta' % gene))

# concat sam and bed
df = pd.concat([sam, bed], axis=1)
df2 = df.apply(row_parse_cigar, axis=1)
df2 = pd.DataFrame.from_dict(dict(zip(df2.index, df2.values))).T
df2.columns = ['aln_span', 'errorCount', 'mutNts','mutLocs','mutCount', 'padDownstream','outputGeneStr']

df3 = df.join(df2[['aln_span', 'errorCount', 'mutNts','mutLocs','mutCount', 'padDownstream']])

# clean UMIs, keep sequence with least mutations
df3 = df3.sort_values(['RNAME', 'UMI', 'mutCount'])
df3 = df3.drop_duplicates(subset=['UMI'], keep='first')
print("Reads aligned to %s after remove UMI dups: %d" % (gene, len(df3)))

# get mut types for each read
df4 = df3.apply(row_get_mut_type, axis=1,  fa=fa, geneStrand=geneStrand, geneStart=geneStart)
df4 = pd.DataFrame.from_dict(dict(zip(df4.index, df4.values))).T
df4.columns = ['refNts', 'muts', 'AtoG_count', 'AtoGs','CtoU_count', 'CtoUs','dualMuts']
df5 = df3.join(df4)
# save result in tab
df5 = df5[["read_id","FLAG","RNAME","POS","errorCount","mutCount","AtoG_count","CtoU_count",
    "refNts","mutNts","mutLocs","padDownstream","muts","AtoGs","CtoUs","dualMuts",
    "MAPQ","strand","blockCount","blockSizes","blockStarts","geneid", "SYMBOL","ENTREZID",
    "overlapLen","aln_span"]]
df5 = df5.replace(r'^\s*$', np.nan, regex=True)
df5.to_csv(os.path.join(mydir,'%s_result.tab' % gene), sep = '\t', na_rep='NA', index=False)
