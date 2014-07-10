#!/home/rcf-40/jcrisci/bin/python
#PBS -q cmb
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=2000mb,pmem=2000mb,vmem=2000mb
#PBS -t 0-20

import gzip
import sys
import os
import glob

arrayID = int(os.getenv('PBS_ARRAYID'))
print arrayID

pathIn = '/home/cmb-panasas2/jcrisci/geno-posts/R-W_combined/'
pathOut = '/home/cmb-panasas2/jcrisci/new-pairwise-pi/minInd2-with-counts-R-W_combined/'

listing = glob.glob(pathIn + '*.geno.gz')
#listing2 = [i for i in listing if not 'pvar' in i]
fileIn = listing[arrayID]

def computePi(aR, aH, aA, bR, bH, bA):
    pi = aR*bA + 0.5*(aR*bH + aH*bR + aH*bA + aA*bH + aH*bH) + aA*bR
    return pi

nInd = 12
whitStart = nInd*3+2 # 38
whitEnd = whitStart + (nInd*3) # 73

def generateHeader(nInd):
    header = []
    header.append('{}   {}'.format('chrom', 'position'))
    for i in range(1, nInd+1):
        for j in range(1, nInd+1):
            header.append('1.{}-2.{}'.format(i, j))
    return header

#print generateHeader(nInd)
#print len(generateHeader(nInd))
            
with gzip.open(fileIn, 'r') as fin:
    fileOut = fileIn.split('/')[-1][:-8] + '.pwp.gz'
    fout = gzip.open(pathOut + fileOut, 'w')
    header = generateHeader(nInd)
    fout.writelines('\t'.join(header) + '\n')
    for line in fin:
        rowPis = []
        col = line.split('\t')
        #print col, len(col)
        #print col[38], col[39], col[40]
        chrom = col[0]
        pos = col[1]
        rowPis.append('{}   {}'.format(str(col[0]), str(col[1])))
        for i in range(0, nInd):
            #print i
            for j in range(whitStart, whitEnd-2, 3):
                #print j
                aR = float(col[i*3+2])
                aH = float(col[i*3+3])
                aA = float(col[i*3+4])
                bR = float(col[j+0])
                bH = float(col[j+1])
                bA = float(col[j+2])
                #print aR, aH, aA, bR, bH, bA
                rowPis.append(str(computePi(aR, aH, aA, bR, bH, bA)))
        #print rowPis
        fout.writelines('\t'.join(rowPis) + '\n')
fout.close()            
