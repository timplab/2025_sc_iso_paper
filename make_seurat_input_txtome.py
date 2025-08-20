#!/usr/bin/env python3

'''creates barcodes.tsv, features.tsv and matrix file from single cell salmon quant.sf files'''

import os
import sys
import argparse
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--input_path', type=str) #path to directory with individual cell salmon quant.sf files
parser.add_argument('-o', '--outpath', type=str) #where output files go

args = parser.parse_args()
input_path=args.input_path
path=args.outpath

def parsedir(input_path,path):
    '''Makes barcodes.tsv, barcode:position dictionary for all cell quant.sf files in input_dir and makes features.tsv from first file'''

    #initialize barcodes, matrix, and features files (Seurat needs to believe they are gzipped but they arent shhhhh)
    if not os.path.exists(path+'/Seurat_Input/'):
        os.makedirs(path+'/Seurat_Input/')
    else:
        shutil.rmtree(path+'/Seurat_Input/')
        os.makedirs(path+'/Seurat_Input/')

    input=input_path+'/'
    bcOut = open(path + '/Seurat_Input/barcodes.tsv.gz', 'w')
    mtxOut = open(path + '/Seurat_Input/matrix.mtx.gz', 'w')
    ftsOut= open(path + '/Seurat_Input/features.tsv.gz', 'w')

    feature_dict={}
    barcode_dict={}
    fileList=[]
    pos=0
    events=0 #will be used in matrix file

    #makes the features.tsv file from the first quant.sf file
    feature_list=[]
    for file in os.listdir(input):
        for line in open(input+file):
            if line.startswith('Name'):
               continue
            else:
                tx=(line.rstrip().split('/t')[0]).split('|')[0]
                feature_list.append(tx)
        break
    for item in feature_list:
        ftsOut.write(item+'\t'+item+'\tGene Expression\n')
    ftsOut.close()

    for file in os.listdir(input):
        pos += 1 #will be used in matrix file for barcode number
        barcode=(file.rstrip().split('_')[2]).split('.')[0]
        #write out barcode file which is just a single column of barcodes
        bcOut.write(barcode + '-1\n')
        barcode_dict[pos]=barcode

        line_num=-1
        for line in open(input+file):
            line_num += 1
            if line.startswith('Name'):
                continue
            numReads=float(line.rstrip().split('\t')[4]) #just take integer count values from salmon for DeSeq2 compatibility
            if float(numReads) > 0:
                events+=1
                feature_dict[events]=[]
                feature_dict[events].append(str(line_num))
                feature_dict[events].append(str(pos))
                feature_dict[events].append(str(int(numReads)))

    #write matrix file that encodes the counts for each feature in each cell
    mtxOut.write("%%MatrixMarket matrix coordinate integer general\n%\n")
    mtxOut.write( str(len(feature_list))+ ' ' + str(pos) + ' ' + str(events) + '\n')
    for entry in feature_dict:
        mtxOut.write(str(feature_dict[entry][0]) + ' ' + str(feature_dict[entry][1]) + ' ' + str(feature_dict[entry][2])+'\n')
    mtxOut.close()
    bcOut.close()


parsedir(input_path,path)

