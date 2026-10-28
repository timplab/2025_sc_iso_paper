#!/usr/bin/env python3
'''wrapper for splitting and processing big input fastqs with C3Poa'''

import os
import sys
import mappy as mm
from tqdm import tqdm
import multiprocessing as mp
import shutil
from glob import glob
import math
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', type=str) #path to big fastq file from catting fastq_pass
parser.add_argument('-o', '--output_path', type=str) #where you want temp output files and final combined R2C2_Consensus.fasta to go
parser.add_argument('-s', '--splint', type=str) #path to 10x_UMI_splint.fasta
parser.add_argument('-n', '--threads', type=int, default=20) #threads passed to each C3POa_sc.py run
parser.add_argument('-c', '--config', type=str, default='') #optional config file passed through to C3POa_sc.py

args = parser.parse_args()
input_file=args.input_file
path=args.output_path
splint=args.splint
threads=args.threads
config=args.config

#C3POa_sc.py lives next to this script
C3POA_SC = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/C3POa_sc.py'

#define max number of reads per file
groupsize=5000000
os.mkdir(path+'/R2C2_temp/')
temp_dir =path+'/R2C2_temp/'

#get number of reads in file and calculate number of iterations needed
def get_iters(input_file, groupsize):
    target=int(groupsize)*4 #number of fastq lines to target for each file (number of readsx4)
    rootname=(input_file.split('.')[0]).split('/')[-1]
    lines=0
    for line in open(input_file):
        lines+=1
    total_iters=math.ceil(lines/target)
    #print(lines)
    return total_iters, target, lines

def make_chunks(input_file,total_iters,target,lines,temp_dir):
    #get file rootname
    rootname=(input_file.split('.')[0]).split('/')[-1]

    read_dict={}
    linenum=0

    iters=1
    for line in open(input_file):
        linenum+=1
        if linenum%4 == 1:
            name=((line.rstrip()).split(' ')[0])[1:]
            read_dict[name]=[]
        if linenum%4 == 2:
            seq=line.rstrip()
            read_dict[name].append(seq)
        if linenum%4 == 3:
            strand=line.rstrip()
            read_dict[name].append(strand)
        if linenum%4 == 0:
            qual=line.rstrip()
            read_dict[name].append(qual)
        if linenum==iters*target or linenum==lines:
            iterate(iters,read_dict,rootname,temp_dir)
            read_dict={}
            iters+=1

def iterate(iteration, read_dict, rootname,temp_dir):
    #make file
    split_fq = temp_dir+'/'+rootname + '_' + str(iteration)+'.fastq'
    write_fastq = open(split_fq, 'w+')
    for entry in read_dict:
        #print(entry,read_dict[entry],iteration)
        write_fastq.write('@%s\n%s\n%s\n%s\n' %(entry,read_dict[entry][0],read_dict[entry][1],read_dict[entry][2]))
    write_fastq.close()

def c3poa(input_path,splint):
    list=os.listdir(input_path)
    for file in list:
        os.mkdir(input_path+file+'_/')
        reads=input_path+file
        output=input_path+file+'_/'
        cmd = '%s -n %d -r %s -o %s -s %s' %(C3POA_SC, threads, reads, output, splint)
        if config:
            cmd += ' -c %s' %(config)
        os.system(cmd)

total_iters,target,lines=get_iters(input_file, groupsize)
make_chunks(input_file,total_iters,target,lines,temp_dir)
c3poa(temp_dir,splint)
#cat files into one large consensus file
rootname=(input_file.split('.')[0]).split('/')[-1]
os.system('cat %s/*/R2C2_Consensus.fasta > %s/%s_R2C2_Consensus.fasta' %(temp_dir, path, rootname))
shutil.rmtree(temp_dir)


