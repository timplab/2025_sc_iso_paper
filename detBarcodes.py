#!/usr/bin/env python3
# Roger Volden
#Sheridan updates
'''
Demultiplex all possible reads by exact barcode matching to 10X whitelist; sort by #reads per barcode
'''

import sys
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', type=str)
parser.add_argument('-w', '--whitelist', type=str)
parser.add_argument('-o', '--output_path', type=str)

args = parser.parse_args()
input_file = args.input_file #barcodes from R2C2
whitelist = args.whitelist #10X whitelist of barcodes, most current vers
path = args.output_path

def fastaReader(inFile):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    readDict = {}
    for line in open(inFile):
        line = line.rstrip()
        if not line:
            continue
        if line.startswith('>'):
            readDict[line[1:]] = ''
            lastHead = line[1:]
        else:
            readDict[lastHead] += line
    return readDict

def readBarcodes(inFile):
    '''Reads in the canonical 10x barcodes'''
    barcodes = {}
    for line in open(inFile):
        line = line.rstrip()
        if not line:
            continue
        barcodes[line] = 0
    return barcodes

def tallyReads(reads, barcodes):
    '''Counts the number of occurrences of each barcode'''
    totalmatched=0
    total=0
    for read in reads:
        total += 1
        cellBC = read[:16]
        if cellBC in barcodes:
            barcodes[cellBC] += 1
            totalmatched +=1
    #print('%s out of %s total reads had exact matches to canonical barcodes'%(totalmatched, total))
    return barcodes,totalmatched

def main():
    name='all_bcs.fasta'
    out=open(path + name,'w+')

    readDict = fastaReader(input_file)
    reads = list(readDict.values())
    barcodes = readBarcodes(whitelist)
    barCounts, totalmatched = tallyReads(reads, barcodes)
    counts = []
    for bc, count in barCounts.items():
        counts.append((count, bc))
    counts = sorted(counts, reverse=True)

    n = 0
    toptotal=0
    for count in counts: #demux EVERYONE
        out.write('>barcode_' + str(n) + '_' + str(count[0])+'\n')
        out.write(count[1]+'\n')
        toptotal += count[0]
        n += 1
    out.close()

    #get rid of all the 0 counts
    dict={}
    linenum=1
    for line in open(path+name):
        #print(line)
        if line.startswith('>'):
            number=int(line.rstrip().split('_')[2])
            if number > 0:
                barcode=line.rstrip()
                dict[barcode]=[]
        else:
            dict[barcode]=line.rstrip()
        linenum+=1
    out=open(path + name,'w+')
    for barcode in dict:
        out.write('%s\n%s\n' %(barcode, dict[barcode]))

if __name__ == '__main__':
    main()
