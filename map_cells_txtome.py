#!/usr/bin/env python3
'''Script for mapping individual cell fasta files to a ref transcriptome and quantifying isoforms with salmon'''
'''dependencies:minimap2/salmon'''



import argparse
import editdistance
import numpy as np
import os
import sys
import shutil
import mappy as mm
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--input_path', type=str) #path UMI collapsed cell fastas (in the /demuxed_sorted/ directory)
parser.add_argument('-o', '--output_path', type=str) #where you want /cellsams_tx/ and /cell_salmon_quants/ dirs to go 
parser.add_argument('-c', '--config_file', type=str) #config file containing paths to executables if not in $PATH
parser.add_argument('-r', '--ref_tx', type=str) #reference TRANSCRIPTOME FASTA for minimap2 alignment


args = parser.parse_args()
input_path=args.input_path
path=args.output_path
config_file=args.config_file
ref_tx=args.ref_tx

'''Parse config'''
def configReader(configIn):
    progs = {}
    for line in open(configIn): 
        if line.startswith('#') or not line.rstrip().split():
            continue
        line = line.rstrip().split('\t') 
        progs[line[0]] = line[1] 
        
    possible = set(['poa', 'minimap2', 'gonk', 'consensus', 'racon', 'samtools','featureCounts', 'psl2pslx']) #don't think I actually need psl2pslx and emtrey...
    inConfig = set() 
    for key in progs.keys(): 
        inConfig.add(key)
        if key not in possible: 
            raise Exception('Check config file')

    for missing in possible-inConfig: 
        if missing == 'consensus': 
            path = 'consensus.py' 
        else:
            path = missing  
        progs[missing] = path 
        #sys.stderr.write('Using ' + str(missing)
                         #+ ' from your path, not the config file.\n')
    return progs

progs = configReader(config_file)
minimap2 = progs['minimap2']


'''Read each cell fasta in input path, create temp trimmed fasta for minimap2 alignment'''

def trim (inFile):
    readdict={}
    filename=inFile.split('/')
    filename=filename[len(filename)-1]
    cell=filename.split('_')[0]
    number=filename.split('_')[1]
    barcode=filename.split('_')[2]
    bc=barcode.split('.')[0]
    for line in open(inFile):
        if line.startswith('>'):
            name=line.split('_')[0]
            name=name+'_'+cell+'_'+number+'_'+bc
        else:
            sequence=line
            trim=len(sequence)-29
            sequence=sequence[0:trim]
            readdict[name]=sequence
    return readdict


def map_reads(name,inFile,ref_tx,input):
    #sys.stderr.write('Mapping reads to %s' %(ref_fasta))
    out= input + 'cellsams_tx/'+ name +'.sam'
    os.system('%s -ax map-ont -N 100 -p .99 %s %s > %s' % (minimap2, ref_tx, inFile, out))
    salmon_quant(out,ref_tx,input,name)

def salmon_quant(infile, ref_tx, input, name):
    os.system('salmon quant -t %s -l U -a %s --ont -o %s' % (ref_tx, infile, input))
    filename=input+'/quant.sf'
    shutil.move(filename,input+'/cell_salmon_quants/'+name+'.sf') 

def concat_fasta(input,ref_tx,path):
    #make directory for cell.sam files
    if not os.path.exists(path+'/cellsams_tx/'):
        os.makedirs(path+'/cellsams_tx/')
    else:
        shutil.rmtree(path+'/cellsams_tx/')
        os.makedirs(path+'/cellsams_tx/')

    #make directory for cell_quant.sf files
    if not os.path.exists(path+'/cell_salmon_quants/'):
        os.makedirs(path+'/cell_salmon_quants/')
    else:
        shutil.rmtree(path+'/cell_salmon_quants/')
        os.makedirs(path+'/cell_salmon_quants/')

    fileList=[]
    for file in os.listdir(input):
        if '.fasta' not in file:
            continue
        name=file.split('.')[0]
        final=open(path+'/temptrimmed.fasta','w')
        temp=path+'/temptrimmed.fasta'
        fasta_file = input+'/' + file
        readdict=trim(fasta_file)
        for i in readdict:
            sequence=readdict[i] 
            final.write('%s\n%s\n' % (i,sequence))
        #final.close()
        map_reads(name,temp,ref_tx,path)
    #return path+'/allfinalreads.fasta'

concat_fasta(input_path, ref_tx, path)


