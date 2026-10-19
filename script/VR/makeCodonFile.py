#!/usr/bin/env python3

## This script is used to caclulate the translation start and stop
## positions in a transcript given a GFF file.

## Script from René Dreos


import HTSeq
import collections
import argparse

def parserFunction():
    # parse the command line arguments:
    parser = argparse.ArgumentParser(description='Given a GTF and a'
                                     ' FASTA file the script '
                                     'outputs each codon coordinates'
                                     ' for all protein coding transcripts'
                                     ' in SGA format. This is used to make'
                                     ' RUST plots.')
    parser.add_argument('-g','--gtf',
                        help='GTF file with genomic regions.',
                        nargs=1,
                        type=str)
    parser.add_argument('-f','--fasta',
                        help='FASTA file with transcripts sequences.',
                        nargs=1,
                        type=str)
    parser.add_argument('--dicodons',
                        help='Get dicodons instead of codons.',
                        action='store_true')
    parser.add_argument('--debug',
                        help='Print verbose output for debugging',
                        action="store_true")
    return parser.parse_args()


def getCDS(gtfFile, debug):
    gtf_file = HTSeq.GFF_Reader( gtfFile, end_included=True )
    intervals = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
    if debug: print('Parsing GTF file')
    tPos = dict()
    cdsLength = collections.Counter( )
    exons = dict()
    transcriptLength = collections.Counter( )
    cdsStart = dict()
    geneId = dict()
    geneBiotype = dict()
    trnsBiotype = dict()
    for feature in gtf_file:
        if len(str(feature.iv.chrom)) > 2:
            continue
        if feature.type == 'transcript':
            transcript_id = feature.attr['transcript_id']
            tPos[transcript_id] = feature.iv
            geneId[transcript_id] = feature.attr['gene_id']
            geneBiotype[transcript_id] = feature.attr['gene_biotype']
            trnsBiotype[transcript_id] = feature.attr['transcript_biotype']
            cdsLength[transcript_id] = 0
            if debug:
                print('Transcript id:', transcript_id)
                print('  GID:', geneId[transcript_id])
                print('  GBT:', geneBiotype[transcript_id])
                print('  TBT:', trnsBiotype[transcript_id])
                print('  CHR:', str(feature.iv.chrom))
        if feature.type == 'exon':
            transcript_id = feature.attr['transcript_id']
            exons[transcript_id] = feature.iv
            transcriptLength[transcript_id] += feature.iv.length
            if debug: print('  Cumulative Transcript length: ', transcriptLength[transcript_id])
        if feature.type == 'CDS':
            transcript_id = feature.attr['transcript_id']
            cdsLength[transcript_id] += feature.iv.length
            if debug:
                print('  CDS length:', feature.iv.length)
                print('  Cumulative CDS length:', cdsLength[transcript_id])
        if feature.type == 'start_codon':
            transcript_id = feature.attr['transcript_id']
            strand = feature.iv.strand
            if transcript_id in cdsStart:
                continue
            else:
                if strand == '+':
                    start = feature.iv.start
                    tStart = tPos[transcript_id].start
                    # cdsStart[transcript_id] = start - tStart
                    cdsStart[transcript_id] = transcriptLength[transcript_id] - (exons[transcript_id].end - start - 1)
                else:
                    start = feature.iv.end
                    tStart = tPos[transcript_id].end
                    # cdsStart[transcript_id] = tStart - start
                    cdsStart[transcript_id] = transcriptLength[transcript_id] - (start - exons[transcript_id].start - 1)
                if debug:
                    print('  Strand: ', strand)
                    print('  CDS start: ', start)
                    print('  Transcript Start: ', tStart)
                    print('  CDS Trns Start: ', cdsStart[transcript_id])
    return(geneId, cdsStart, cdsLength)

def getFASTA(fasta_file, debug):
    if debug: print('Parsing FASTA file')
    sequences = dict()
    for sequence in HTSeq.FastaReader( fasta_file ):
        if debug: print("Sequence '%s' has length %d." % ( sequence.name, len(sequence) ))
        ## names = sequence.name.split(' ')
        ## if debug: print("Sequence new name: '%s'" % ( names[0] ))
        sequences[sequence.name] = sequence
    return(sequences)



args = parserFunction()

debug   = args.debug
dicodons = args.dicodons
gtfFile = args.gtf[0]
fastaFile = args.fasta[0]

if dicodons:
    shift = 6
else:
    shift = 3

geneId, cdsStart, cdsLength = getCDS(gtfFile, debug)
sequences = getFASTA(fastaFile, debug)

for transcript_id in sorted( geneId ):      
    if transcript_id in cdsStart:
        if (cdsLength[transcript_id] % 3) == 0:
            start = cdsStart[transcript_id] - 1
            end = start + cdsLength[transcript_id]
            seqID = geneId[transcript_id] + '|' + transcript_id
            outID = transcript_id
            for p in range(start, end, 3):
                codon = sequences[seqID][p:p+shift]
                if len(codon) == shift:
                    print(outID, codon, p+1, '+', 1, sep="\t")
