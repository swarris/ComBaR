import sys

combine = set()

for i in sys.argv[2:]:
    for j in sys.argv[2:]:
        if i != j and (i,j) not in combine and (j,i) not in combine:
            combine.add((i,j))

out = 0
for c in combine:
    print("python combar.py --program=linker -L /tmp/log.txt --loglevel=debug -o {dir}/links/{out}.csv {in1} {in2}".format(out=out, dir=sys.argv[1], in1=c[0],in2=c[1]))
    out += 1
                                    
