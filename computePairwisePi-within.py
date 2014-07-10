#!/usr/bin/env python
#PBS -q cmb
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=2000mb,pmem=2000mb,vmem=2000mb
#PBS -t 0-41

import gzip
import sys
import glob
import os

arrayID = int(os.getenv('PBS_ARRAYID'))
print arrayID

pathIn = '/home/cmb-panasas2/jcrisci/geno-posts/minInd2-withCounts/'
pathOut = '/home/cmb-panasas2/jcrisci/new-pairwise-pi/minInd2-withCounts/'
listing = glob.glob(pathIn + '*.geno.gz')
print listing
#listing2 = [i for i in listing if not 'pvar' in i]
fileIn = listing[arrayID]
print fileIn

def computePi(aR, aH, aA, bR, bH, bA):
    pi = aR*bA + 0.5*(aR*bH + aH*bR + aH*bA + aA*bH + aH*bH) + aA*bR
    return pi

def generateHeader(nInd):
    header = []
    header.append('{0}   {1}'.format('chrom', 'position'))
    x = 2
    for i in range(1, nInd+1):
        for j in range(x, nInd+1):
            header.append('{0}-{1}'.format(i, j))
        x += 1
    return header

nInd = 12

fin = gzip.open(fileIn, 'r')
fileOut = fileIn.split('/')[-1][:-8] + '.pwp.gz'
print fileOut
fout = gzip.open(pathOut + fileOut, 'w')
print pathOut+fileOut
header = generateHeader(nInd)
fout.writelines('\t'.join(header) + '\n')
for line in fin:
    rowPis = []
    col = line.split('\t')
    chrom = col[0]
    pos = col[1]
    rowPis.append('{0}   {1}'.format(str(col[0]), str(col[1])))
    x = 1
    for i in range(0, nInd):
        for j in range(x, nInd):
            aR = float(col[i*3+2])
            aH = float(col[i*3+3])
            aA = float(col[i*3+4])
            bR = float(col[j*3+2])
            bH = float(col[j*3+3])
            bA = float(col[j*3+4])
            #print aR, aH, aA, bR, bH, bA
            rowPis.append(str(computePi(aR, aH, aA, bR, bH, bA)))
        x += 1
    #print rowPis
    fout.writelines('\t'.join(rowPis) + '\n')
fout.close()
fin.close()
