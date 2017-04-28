import sys
from collections import defaultdict


chrs = sys.argv[1].split(",")
minimumWidth = int(sys.argv[2])
bundles = sys.argv[3:]


chromosomes = set()
chromosomesOrder = list(chrs)

for bFile in bundles:
    links = defaultdict(dict)
    for b in open(bFile, "r"):
        b = b.split(" ")
        if b[0] != b[3] : # different sequences
            chr1 = None
            chr2 = None
            width = 0
            if (b[0] in chrs and b[3] not in chrs) or (b[3] in chrs and b[0] not in chrs):
                if (b[0] in chrs and b[3] not in chrs):
                    chr1 = b[0]
                    chr2 = b[3]
                else:
                    chr1 = b[3]
                    chr2 = b[0]
                width = abs(int(b[5]) - int(b[4]))
            if chr1 != None and chr2 != None and width > minimumWidth:
                if chr2 not in links[chr1]:
                    links[chr1][chr2] = 0
                links[chr1][chr2] += width
  
    for l in links.keys():
        bundleWidths = []
        if l not in chromosomes:
            chromosomes.add(l)
        if l not in chromosomesOrder:
            chromosomesOrder.append(l)
        for t in links[l].keys():
            bundleWidths.append((links[l][t], t))
        bundleWidths.sort(reverse=True)
        for b in bundleWidths:
            chromosomes.add(b[1])
            #print(b)
            if b[1] not in chromosomesOrder:
                chromosomesOrder.append(b[1])
        
print("chromosomes = {}".format(";".join(list(chromosomes))))
print("chromosomes_order = {}".format(",".join(chromosomesOrder)))


 
                     
    
