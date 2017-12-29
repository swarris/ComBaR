import sys

linksFile = open(sys.argv[1], "r")
noRepeatsFile = open(sys.argv[3], "w")
links = []

repeatLinks = int(sys.argv[2])

numberOfRepeats = 0

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

def writeNoRepeats(links):
    for l in links:
        bundle = [l[0], l[2], l[2], l[1], l[3], l[3], l[4],l[5]]
        if not bundle[7] :
           noRepeatsFile.write(" ".join([str(x) for x in bundle[:-1]]) + "\n")

markRepeats(links, repeatLinks)
print("Number of repeated links: {}".format(numberOfRepeats))

writeNoRepeats(links)
                 
