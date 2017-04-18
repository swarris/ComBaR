import sys

linksFile = open(sys.argv[1], "r")
bundlesFile = open(sys.argv[2], "w")
repeatsFile = open(sys.argv[3], "w")
links = []

maxGap = 10**6
minLinks = 10

maxGapRepeat = 10**4
minLinksRepeat = 10
repeatLinks = 10

numberOfBundles = 0
totalWidth = 0
numberOfRepeats = 0
totalWidthRepeats  = 0
# link structure: chr1 chr2 start1 start2 repeat?

for l in linksFile:
    l = l.strip().split("\t")
    links.append([l[0], l[3], int(l[1]), int(l[4]), l[6], False])
    links.append([l[3], l[0], int(l[4]), int(l[1]), l[6], False])

print("Number of links: {}".format(len(links)/2))

links.sort()

linksPresent = set()

def markRepeats(links, repeatLinks):
    global numberOfRepeats
    currentRepeat = []
    for l in links:
        if len(currentRepeat) > 0 and (l[0] != currentRepeat[-1][0] or l[2] != currentRepeat[-1][2]):
            if len(currentRepeat) >= repeatLinks:
                for cR in currentRepeat:
                    cR[5] = True
                numberOfRepeats += 1
            currentRepeat = []
        elif len(currentRepeat) > 0  and  l[0] == currentRepeat[-1][0] and l[2] == currentRepeat[-1][2]:
            currentRepeat.append(l)
        else:
            currentRepeat = [l]

def bundleLinks(links):
    global numberOfBundles
    global totalWidth
    global linksPresent
    global totalWidthRepeats
    global maxGapRepeat
    global minLinksRepeat
    
    remainingLinks = []
    bundle = []
    numberOfLinks = 1
    for l in links:
        if (l[0], l[2], l[1], l[3], l[4]) not in linksPresent:
            if bundle == []:
                linksPresent.add((l[0], l[2], l[1], l[3], l[4]))
                linksPresent.add((l[1], l[3], l[0], l[2], l[4]))
                bundle = [l[0], l[2], l[2], l[1], l[3], l[3], l[4],l[5]]
            else:
                if bundle[0] == l[0] and bundle[3] == l[1] and bundle[7] == l[5] and bundle[6] == l[4] and ((abs(bundle[2] - l[2]) < maxGap and abs(bundle[5] - l[3]) < maxGap and not bundle[7]) or (abs(bundle[2] - l[2]) < maxGapRepeat and abs(bundle[5] - l[3]) < maxGapRepeat and bundle[7])): # extend bundle
                    if bundle[2] != l[2] or bundle[5] != l[3]:
                        numberOfLinks += 1
                    bundle[2] = l[2]
                    bundle[5] = l[3]
                    linksPresent.add((l[0], l[2], l[1], l[3], l[4]))
                    linksPresent.add((l[1], l[3], l[0], l[2], l[4]))
                elif bundle[0] != l[0] or (bundle[0] == l[0] and bundle[3] == l[1] and bundle[6] == l[4] and (abs(bundle[2] - l[2]) > maxGap or abs(bundle[5] - l[3]) > maxGap)):
                    if not bundle[7] and numberOfLinks > minLinks:
                        bundlesFile.write(" ".join([str(x) for x in bundle[:-1]]) + "\n")
                        numberOfBundles += 1
                        totalWidth += bundle[2] - bundle[1] 
                    elif bundle[7] and numberOfLinks > minLinksRepeat:
                        repeatsFile.write(" ".join([str(x) for x in bundle[:-1]]) + "\n")
                        totalWidthRepeats += bundle[2] - bundle[1] 
                    numberOfLinks = 1
                    bundle = [l[0], l[2], l[2], l[1], l[3], l[3], l[4],l[5]]
                    linksPresent.add((l[0], l[2], l[1], l[3], l[4]))
                    linksPresent.add((l[1], l[3], l[0], l[2], l[4]))
                else:
                    remainingLinks.append(l)

    return remainingLinks

markRepeats(links, repeatLinks)
print("Number of repeated links: {}".format(numberOfRepeats))

while len(links)>0:
    links =bundleLinks(links)
print("Average width repeats: {}".format(totalWidthRepeats/numberOfRepeats))
print("Number of bundles: {}".format(numberOfBundles))
print("Average width: {}".format(totalWidth/numberOfBundles))
        
                 
