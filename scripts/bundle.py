import sys
from collections import defaultdict

linkFile = open(sys.argv[1],"r")
bundleFile = open(sys.argv[2],"w")

minDistance = 500
maxDistance = 10**5
minLinks = 5

links = {}

for l in linkFile:
    l = l.strip().split("\t")
    chr1 = l[0]
    chr2 = l[3]
    chr1X = int(l[1])
    chr2X = int(l[4])
    distance = l[6]
    if chr1 < chr2:
        connect = (chr1,chr2)
        link = [chr1X, chr2X,0]
    else:
        connect = (chr2, chr1)
        link = [chr2X, chr1X,0]

    if connect not in links:
        links[connect] = defaultdict(list)
        
    links[connect][distance].append(link)

def compare(x,y):
    return x[2] - y[2]
    
for l in links:
    for d in links[l]:
        bundles = 0
        links[l][d].sort()
        #print(links[l][d])
        for i in range(0, len(links[l][d])):
            currentLink = links[l][d][i] 
            if currentLink[2] == 0:
                # not yet in a bundle:
                bundles += 1
                currentLink[2] = bundles
            if "-1" in d:
                findIndex = i+1
                while findIndex  < len(links[l][d]) and (links[l][d][findIndex][0] == currentLink[0] or links[l][d][findIndex][1] > currentLink[1]) and abs(links[l][d][findIndex][0]-currentLink[0]) < maxDistance and abs(links[l][d][findIndex][1]-currentLink[1]) < maxDistance:
                    
                    findIndex += 1
                if findIndex < len(links[l][d]) and abs(links[l][d][findIndex][0]-currentLink[0]) < maxDistance and abs(links[l][d][findIndex][1]-currentLink[1]) < maxDistance:
                    links[l][d][findIndex][2] = currentLink[2] # same bundle
            else:
                findIndex = i+1
                while findIndex < len(links[l][d]) and (links[l][d][findIndex][0] == currentLink[0] or links[l][d][findIndex][1] <= currentLink[1]) and abs(links[l][d][findIndex][0]-currentLink[0]) < maxDistance and abs(links[l][d][findIndex][1]-currentLink[1]) < maxDistance:
                    findIndex += 1
                if findIndex < len(links[l][d]) and abs(links[l][d][findIndex][0]-currentLink[0]) < maxDistance and abs(links[l][d][findIndex][1]-currentLink[1]) < maxDistance:
                    links[l][d][findIndex][2] = currentLink[2] # same bundle
#        print(links[l][d])
        # sort on bundles
        links[l][d].sort(compare)
        currentBundle = 1
        chr1 = [0,0]
        chr2 = [0,0]
        bundleSize = 1
        for b in links[l][d]:
            if b[2] == currentBundle:
                bundleSize += 1
                if chr1[0] == 0:
                    chr1[0] = b[0]
                    chr1[1] = b[0]
                    chr2[0] = b[1]
                    chr2[1] = b[1]
                    
                if b[0] > chr1[1]:
                    chr1[1] = b[0]
                if b[0] < chr1[0]:
                    chr1[0] = b[0]
                if b[1] > chr2[1]:
                    chr2[1] = b[1]
                if b[1] < chr2[0]:
                    chr2[0] = b[1]
#                if "-1" in distance:
#                    if b[1] < chr2[0]:
#                        chr2[0] = b[1]
#                else:
#                    if b[1] > chr2[1]:
#                        chr2[1] = b[1]
            if b[2] != currentBundle or b == links[l][d][-1]:
                if bundleSize >= minLinks and abs(chr1[0] - chr1[1]) > minDistance and abs(chr2[0] - chr2[1]) > minDistance:
                    if "-1" in d:
                        bundleFile.write("{}\t{}\t{}\t{}\t{}\t{}\t{},bundle={},size={}\n".format(l[0],chr1[0],chr1[1], l[1],chr2[1],chr2[0], d, currentBundle, bundleSize))
                    else:
                        bundleFile.write("{}\t{}\t{}\t{}\t{}\t{}\t{},bundle={},size={}\n".format(l[0],chr1[0],chr1[1], l[1],chr2[0],chr2[1], d, currentBundle, bundleSize))
                currentBundle = b[2]
                bundleSize = 1
                chr1 = [0,0]
                chr2 = [0,0]
