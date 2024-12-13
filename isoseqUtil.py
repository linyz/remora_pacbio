#!/usr/bin/env python3

# Yizhu Lin 10/15/2021
# for isoseq isoform & mutation analysis
import re
def rev_comp(seq):
    '''rev-comp of a seqeunce'''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

def parse_cigar(cigar, seq, strand):
    """parse cigar str, output 2 arrays, one for mutations (X), one for seq occupancy (MDN=X)
    correct strand, output should be 5' to 3'  """
    cigarList = re.split('(\d+)',cigar)
    cigarLocs = [int(x) for x in cigarList[1::2]]
    cigarSTRs = cigarList[2::2]
    if strand == '-':
        cigarLocs = cigarLocs[::-1]
        cigarSTRs = cigarSTRs[::-1]
        seq = rev_comp(seq)
    outputStr = ''
    errorCount = 0
    seqLoc = 0
    mutNts = ''
    mutLocs = []
    mutCount = 0
    aln_span = 0
    for i, cigarChar in enumerate(cigarSTRs):
        if cigarChar in 'MD=':
            outputStr += cigarChar*cigarLocs[i]
            aln_span += cigarLocs[i]
        if cigarChar == 'N':
            outputStr = outputStr+'{0:0>8}'.format(cigarLocs[i])+'N'
            aln_span += cigarLocs[i]
        if cigarChar in 'IDX':
            errorCount += cigarLocs[i]
        if cigarChar in 'MIS=':
            seqLoc += cigarLocs[i]
        if cigarChar == 'X':
            for j in range(cigarLocs[i]):
                mutCount += 1   
                mutNts += seq[seqLoc]
                aln_span += 1
                mutLocs.append(aln_span)   
                seqLoc += 1
                
        else:
            continue

    return outputStr, errorCount, mutNts,mutLocs,mutCount,aln_span

def row_parse_cigar(row):
    """ parse row (after join sam and intersect), get mutation info from CIGAR"""
    cigar = row['CIGAR']
    seq = row['SEQ']
    strand = row['strand'] # this is read strand (mapped to genome)
    geneStart = row['annoStart'] 
    geneEnd = row['annoEnd']
    outputStr, errorCount, mutNts,mutLocs,mutCount,aln_span = parse_cigar(cigar, seq, strand)
    alnStart = row['POS']
    alnEnd = row['POS'] + aln_span
    
    # pad only 3'
    if strand == '+':
        padDownstream = int(geneEnd - alnEnd)
    else:
        padDownstream = int(alnStart - geneStart)
    if padDownstream >= 0:
        outputGeneStr = outputStr + 'V'*padDownstream
    else:
        outputGeneStr = outputStr[0:padDownstream]
        # mutLoc to chr loc
    mut_locs_in_chrom = []
    if strand == '+':
        for loc in mutLocs:
            mut_locs_in_chrom.append(alnStart+loc-1)
    else:
        for loc in mutLocs:
            mut_locs_in_chrom.append(alnEnd-loc)       
    return [aln_span, errorCount, mutNts, '_'.join([str(x) for x in mut_locs_in_chrom]), mutCount, padDownstream, outputGeneStr]

def parse_gff_info(row, info_col = 'info'):
    """parse gff input annotation column"""
    info = row[info_col]
    [gene_id, transcript_id] = info.split(";")[0:2]
    gene_id = gene_id.strip().lstrip('gene_id').strip().strip('"')
    transcript_id = transcript_id.strip().lstrip('transcript_id').strip().strip('"')
    exon_id = ''
    if row['type'] == 'exon':
        exon_id = gene_id+'_'+str(row['start'])+'_'+str(row['end'])
    return [gene_id, transcript_id, exon_id]

def parse_gencode_gff(row, info_col = 'annoInfo', sep='='):
    """parse gff input"""
    # sep "=" if gff, " " if gtf

    info = row[info_col]
    items = info.split(';')
    ID = ''
    gene_id = ''
    gene_name = ''
    for item in items:
        item = item.strip()
        if item.startswith('ID'):
            ID = item.split(sep)[1]
        elif item.startswith('gene_id'):
            gene_id = item.split(sep)[1]
        elif item.startswith('gene_name'):
            gene_name = item.split(sep)[-1]
    return  [ID, gene_id, gene_name]

def exon_dist(exonA, exonB):
    """calculate distance between 2 exons, by start / end locations"""
    [exonA_start, exonA_end] = exonA.split('_')[1:2]
    [exonB_start, exonB_end] = exonB.split('_')[1:2]
    dist = (int(exonA_start) - int(exonB_start))**2 + (int(exonA_end) - int(exonB_end))**2 
    return dist

def anno_to_exons(row):
    """exons orientation and locations"""
    blockSizes = [int(x) for x in row['blockSizes'].split(',')[0:-1]]
    blockStarts = [int(x) for x in row['blockStarts'].split(',')[0:-1]]
    # return exon from 5' to 3'
    exons = []
    if row['strand'] == '-':
        blockSizes = blockSizes[::-1]
        blockStarts = blockStarts[::-1]
    for i, value in enumerate(blockSizes):
        exonStart = int(row['chromStart'])+blockStarts[i]
        exonEnd = exonStart + value
        exons.append('_'.join([row['chrom'], str(exonStart), str(exonEnd)]))
    return exons
    
def read_fa(filename):
    """read fasta file (one entry per file)"""
    f = open(filename, "r")
    lines = f.readlines()
    [geneid,SYMBOL,chrom,geneStrand, grange] = lines[0].lstrip('>').split(':')
    [geneStart, geneEnd] = grange.split('-')
    seq = ''.join(x.strip() for x in lines[1:])
    return seq, geneStrand, geneStart, geneEnd

def get_ref_nt(fa, geneStrand,geneStart, loc):
    """return nt in reference sequence"""
    locG = int(loc) - int(geneStart)
    if locG >= len(fa) or locG <0:
        return 'O' # out of gene annotation boundary
    else:
        if geneStrand == '+':
            return fa[locG]
        else:
            return fa[::-1][locG]

def get_mut_type(mutLoc, mutNt, fa, geneStrand, geneStart):
    """get refs and mut types for all point muts"""
    refNt = get_ref_nt(fa, geneStrand, geneStart, mutLoc)
    return refNt, '%s>%s' % (refNt, mutNt)

def row_get_mut_type(row, fa, geneStrand, geneStart):
    """get all mut for a read (a row)"""
    mutLocs = row['mutLocs'].split('_')
    mutNts = row['mutNts']
    l = len(mutNts)
    refNts = ''
    muts = []
    AtoG_count = 0
    AtoGs = []
    CtoU_count = 0
    CtoUs = []
    dualMuts = []
    for i in range(l):
        mutLoc = mutLocs[i]
        mutNt = mutNts[i]
        refNt, mutType = get_mut_type(mutLoc, mutNt, fa, geneStrand, geneStart)
        refNts = refNts + refNt
        mut = mutLoc+':'+mutType
        muts.append(mut)
        if mutType == 'C>T':
            CtoU_count += 1
            CtoUs.append(mutLoc)
            dualMuts.append(mut)
        elif mutType == 'A>G':
            AtoG_count += 1
            AtoGs.append(mutLoc)
            dualMuts.append(mut)
    return refNts, '_'.join(muts), str(AtoG_count), '_'.join(AtoGs), str(CtoU_count), '_'.join(CtoUs), '_'.join(dualMuts)
