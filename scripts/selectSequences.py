import sys
from collections import defaultdict

# command structure: selectSequences.py chr1,chr2,chr3:0-10000 1000 bundles*.csv

chrs = sys.argv[1].split(",")
minimumWidth = int(sys.argv[2])
bundles = sys.argv[3:]

# store command line regions and new regions based on links:
regions = defaultdict(list)
newRegions = defaultdict(list)

#check for regions: 'chr1:1000-1100'
for c in range(0,len(chrs)):
    if ":" in chrs[c]:
        region = chrs[c].split(":")[1].split("-")
        chrs[c] = chrs[c].split(":")[0]
        regions[chrs[c]] = [int(x) for x in region]
        
# list of chromosomes to show in circos and their order
chromosomes = set()
chromosomesOrder = list(chrs)

# parse bundle / link files
for bFile in bundles:
    # store links in dict. Structure: links[chr1][chr2] -> total width of bundles
    links = defaultdict(dict)
    for b in open(bFile, "r"):
        b = b.split(" ") # split on space (bundles)
        if len(b) == 1:
            b = b[0].split("\t") # split on tab (links)
        if b[0] != b[3] : # different sequences
            chr1 = None
            chr2 = None
            chr1xy = [0,0]
            chr2xy = [0,0]
            
            width = 0
            if (b[0] in chrs and b[3] not in chrs) or (b[3] in chrs and b[0] not in chrs):
                if (b[0] in chrs and b[3] not in chrs):
                    chr1 = b[0]
                    chr2 = b[3]
                    chr1xy = [int(b[1]),int(b[2])]
                    chr2xy = [int(b[4]),int(b[5])]

                else:
                    chr1 = b[3]
                    chr2 = b[0]
                    chr2xy = [int(b[1]),int(b[2])]
                    chr1xy = [int(b[4]),int(b[5])]
                chr1xy.sort()    
                chr2xy.sort()    
                width = abs(int(b[5]) - int(b[4]))
            if chr1 != None and chr2 != None and width >= minimumWidth:
                if ((chr1 in regions and regions[chr1][0] <= chr1xy[0] <= chr1xy[1] <= regions[chr1][1]) or chr1 not in regions) and \
                   ((chr2 in regions and regions[chr2][0] <= chr2xy[0] <= chr2xy[1] <= regions[chr2][1]) or chr2 not in regions):
                    if chr1 not in regions and len(regions.keys()) > 0:
                        if chr1 not in newRegions:
                            newRegions[chr1] = chr1xy
                        else:
                            if chr1xy[0] < newRegions[chr1][0]:
                                newRegions[chr1][0] = chr1xy[0]
                            if chr1xy[1] > newRegions[chr1][1]:
                                newRegions[chr1][1] = chr1xy[1]
                                
                    if chr2 not in regions and len(regions.keys()) > 0:
                        if chr2 not in newRegions:
                            newRegions[chr2] = chr2xy
                        else:
                            if chr2xy[0] < newRegions[chr2][0]:
                                newRegions[chr2][0] = chr2xy[0]
                            if chr2xy[1] > newRegions[chr2][1]:
                                newRegions[chr2][1] = chr2xy[1]
                                
                                
                    # add link
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

chromosomes = list(chromosomes)
for c in range(0, len(chromosomes)):
    chro = chromosomes[c]
    if chro in regions:
        chromosomes[c] = chro + ":" + "-".join([str(x) for x in regions[chro]])
    elif chro in newRegions:
        chromosomes[c] = chro + ":" + "-".join([str(x) for x in newRegions[chro]])
    
    
print("chromosomes = {}".format(";".join(list(chromosomes))))
print("chromosomes_order = {}".format(",".join(chromosomesOrder)))


 
                     
    
