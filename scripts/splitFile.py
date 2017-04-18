import sys
from itertools import islice

from Bio import  SeqIO
from Bio.Seq import Seq

fileName = sys.argv[1]
fileElements = open(fileName,"r")

splitAt = int(sys.argv[2])
index = 0
count = 0
outFile = open(fileName + "_" + str(index) + ".fasta", "w")
for i in SeqIO.parse(fileElements, "fasta"):
    if count > splitAt:
        count = 0
        index += 1
        outFile.close()
        outFile = open(fileName + "_" + str(index) + ".fasta", "w")
    else:
        count += len(i.seq)
    outFile.write(">" + str(i.id) +"\n")
    outFile.write(str(i.seq) + "\n")

