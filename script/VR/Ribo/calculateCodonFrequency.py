#!/usr/bin/env python3

## This script is used to caclulate the translation start and stop
## positions in a transcript given a GFF file.

## Script from René Dreos


import sys
import math
import argparse

def parserFunction():
    # parse the command line arguments:
    parser = argparse.ArgumentParser(description='Given the SGA file genearated by'
                                     ' makeCodonFile.py, it calculates the codon '
                                     'frequencies for each transcript.')
    parser.add_argument('-i','--input',
                        help='SGA file with codon positions.',
                        nargs=1,
                        type=str)
    parser.add_argument('-n','--norm',
                        help='Report normalisation factor, otherwise print density',
                        action="store_true")
    parser.add_argument('-d','--debug',
                        help='Print verbose output for debugging',
                        action="store_true")
    return parser.parse_args()

def initializeCodon(step):
    c = ["A", "C", "G", "T"]
    x = dict()
    if step == 3:
        for I in c:
            for K in c:
                for J in c:
                    x[I+K+J] = 0
    else:
        for I in c:
            for K in c:
                for J in c:
                    for L in c:
                        for M in c:
                            for N in c:
                                x[I+K+J+L+M+N] = 0
    # x = {
    #     "TTT":0,
    #     "TCT":0,
    #     "TAA":0,
    #     "TAG":0,
    #     "TAT":0,
    #     "TGA":0,
    #     "TGT":0,
    #     "TTC":0,
    #     "TCC":0,
    #     "TAC":0,
    #     "TGC":0,   
    #     "TTA":0,
    #     "TCA":0,
    #     "TTG":0,
    #     "TCG":0,
    #     "TGG":0,
    #     "CTT":0,
    #     "CCT":0,
    #     "CAT":0,
    #     "CGT":0,
    #     "CTC":0,
    #     "CCC":0,
    #     "CAC":0,
    #     "CGC":0,
    #     "CTA":0,
    #     "CCA":0,
    #     "CAA":0,
    #     "CGA":0,
    #     "CTG":0,
    #     "CCG":0,
    #     "CAG":0,
    #     "CGG":0, 
    #     "ATT":0,
    #     "ACT":0,
    #     "AAT":0,
    #     "AGT":0,
    #     "ATC":0,
    #     "ACC":0,
    #     "AAC":0,
    #     "AGC":0,
    #     "ATA":0,
    #     "ACA":0,
    #     "AAA":0,
    #     "AGA":0,
    #     "ATG":0,
    #     "ACG":0,
    #     "AAG":0,
    #     "AGG":0, 
    #     "GTT":0,
    #     "GCT":0,
    #     "GAT":0,
    #     "GGT":0,
    #     "GTC":0,
    #     "GCC":0,
    #     "GAC":0,
    #     "GGC":0,
    #     "GTA":0,
    #     "GCA":0,
    #     "GAA":0,
    #     "GGA":0,
    #     "GTG":0,
    #     "GCG":0,
    #     "GAG":0,
    #     "GGG":0 
    # }
    return(x)
    

def parseSga(inf, debug):
    allCodonCounts = dict()
    tId            = ''
    cdsLength      = dict()
    begPos         = dict()
    endPos         = dict()
    Ncod           = dict()
    for line in inf:
        fields = line.split('\t')
        codon = fields[1]
        lcodon = len(codon)
        if debug: print('Codon length:', lcodon)
        if (tId == ''):
            pos = int(fields[2])
            tId = fields[0]
            begPos[tId] = pos
            Ncod[tId] = 1
            transCodonCounts = initializeCodon(lcodon)
            transCodonCounts[codon] = 1
            if debug:
                print('Starting reading SGA')
                print('Transcript ID:', tId)
                print('Codon:', codon)
                print('CDS start position:', pos)
        elif(tId == fields[0]):
            pos = int(fields[2])
            transCodonCounts[codon] += 1
            Ncod[tId] += 1
            if debug: print(codon, pos, '=', transCodonCounts[codon])
        else:
            endPos[tId] = pos
            allCodonCounts[tId] = transCodonCounts
            if debug:
                print('End position for', tId, '=', endPos[tId])
                print('Number of codons for', tId, '=', Ncod[tId])
                print(allCodonCounts[tId])
            pos = int(fields[2])
            transCodonCounts = initializeCodon(lcodon)
            transCodonCounts[codon] = 1
            tId = fields[0]
            begPos[tId] = pos
            Ncod[tId] = 1
            if debug:
                print('Change transcript to', tId)
                print('Codon:', codon)
                print('CDS start position:', begPos[tId])
    endPos[tId] = pos
    allCodonCounts[tId] = transCodonCounts
    Ncod[tId] += 1
    ## return(allCodonCounts, begPos, endPos, Ncod)
    return(allCodonCounts, Ncod)

def rescale(x):
    xN = x
    if sum(x) > 0:
        xmax = max(x)
        xmin = min(x)
        xAv = sum(x) / len(x)
        for i in range(len(x)):
            xN[i] = math.ceil(((x[i] - xmin)/(xmax - xmin))*999+1)
    return(xN)
        
def normalize(x):
    xN = x
    xAv = sum(x) / len(x)
    for i in range(len(x)):
        if x[i] > 0:
            xN[i] = xAv / x[i]
        else:
            xN[i] = 0
    return(rescale(xN))

def codonDist(cCounts, Ncod, norm):
    tIds = sorted(Ncod)
    codons = sorted(cCounts[tIds[0]])
    cDist = dict()
    cNorm = dict()
    for T in tIds:
        for C in codons:
            fcod = 1*codonCounts[T][C]/Ncod[T]
            ## print(T, C, fcod, sep="\t")
            if C in cDist:
                cDist[C].append(fcod)
            else:
                cDist[C] = []
                cDist[C].append(fcod)
    for C in codons:
        cmin = round(min(cDist[C]), 4)
        cmax = round(max(cDist[C]), 4)
        cave = round(sum(cDist[C]) / len(cDist[C]), 4)
        if norm:
            cNorm[C] = normalize(cDist[C])
        else:
            cNorm[C] = cDist[C]
            ## print('Codon', C, ': min =', cmin, 'max =', cmax, 'average =', cave)
    return cNorm

######################################################################
### main start here
if __name__ == '__main__':
    args  = parserFunction()
    debug = args.debug
    norm  = args.norm
    if args.input:
        inf = open(args.input[0])
    else:
        inf = sys.stdin

    codonCounts, Ncod = parseSga(inf, debug)
    codonNorm = codonDist(codonCounts, Ncod, norm)
    tIds = sorted(Ncod)
    codons = sorted(codonCounts[tIds[0]])
    for i in range(len(tIds)):
        for C in codons:
            print(tIds[i], C, codonNorm[C][i], sep='\t')
