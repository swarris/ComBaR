import sys
from Bio import SeqIO

seqFile = sys.argv[1]
fileType = sys.argv[2]
color = sys.argv[3]

#chr - hs1 1 0 249250621 chr1

for s in SeqIO.parse(open(seqFile,"r"), fileType):
    print("chr - {id} {id} 0 {length} {color}".format(id=s.id, length=len(s.seq), color=color))
