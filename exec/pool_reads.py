import fileinput
import sys

fileToProcess = open(sys.argv[1])
ids_1={}
ids_2={}

for line in fileToProcess:
    prts=line.split("\n")[0].split("\t")
    tmpid=prts[0]+"\t"+prts[1]+"\t"+prts[2]+"\t"+prts[3]
    if(ids_1.has_key(tmpid)):
        ids_1[tmpid]+=int(prts[4])
        ids_2[tmpid]+=int(prts[5])
    else:
        ids_1[tmpid]=int(prts[4])
        ids_2[tmpid]=int(prts[5])
fileToProcess.close()

for k in ids_1.keys():
    print k+"\t"+str(ids_1[k])+"\t"+str(ids_2[k])
