import os,sys
import fastq as fq
inp=open(sys.argv[1],"r")
oup=open(sys.argv[1]+".enzyme.fa","w")

for fo in fq.read(sys.argv[1]):
    #b=fo.body.count("TTCGAA")#BstBI
    b=fo.body.count("ATCGAT")#ClaI
    if b>=1:
        oup.write("%s\t%s\n"%(fo.head,b))
        #fq.write(fo,sys.argv[1]+".filter",mode="a")
inp.close()
oup.close()
