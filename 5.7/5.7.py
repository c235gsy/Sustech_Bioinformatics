# -*- coding: utf-8 -*-
from copy import deepcopy


def strToInt(list):
    for m in range (0,len(list)):
        list[m] = int(list[m])
    return list

def numTochr(num):
    chrome = "chr"+num
    if num == "MT":
        chrome = "chrM"
    return chrome

def getFeature_KB(gene):
    name = gene[0]
    Flength = 0
    for i in range (1,len(gene)):
        Flength += gene[i][1]-gene[i][0]
    return [name,Flength/1000.00]

def totalNum(DIC):
    l = 0
    for m in DIC.values():
        l += len(m)
    return l

def RPKM(chrome, gene):
    L = float(getFeature_KB(gene)[1])
    #print (L)
    N = float(mappable_reads)
    #print (N)
    C = 0.00
    result = 0.00
    gene = gene[1:]
    for read in obj[chrome]:
        #print (read)
        location = read
        for ex in gene:
            #print (ex)
            down = int(ex[0])
            up = int(ex[1])
            if location in range(down,up+1) or (location + read_length) in range(down,up+1):
                C += 1
                #print (C)
                break
    if C != 0:
        result = C / (L * N)
    return result



global pool
pool =["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"]



refFile = "./gene.txt"

refInput = open(refFile,"r")
genome = refInput.readlines()
del (genome[0])
#print (genome[0])

global num_lab
num_lab=len(genome)

for p in range (0,num_lab):
    genome[p] = genome[p].split("\t")


global lab
lab={}
for n in range(1,23):
    lab["chr%s"%n] = []
    lab["chrM"] = []
    lab["chrX"] = []
    lab["chrY"] = []

global chromes
chromes = list(lab.keys())

global obj
obj = deepcopy(lab)

for gene in genome:
    name = gene[1]
    chr = gene[2]
    #print (chr)
    if "_" in chr:
        continue
    #if chr == "chrM":
        #print (chr)
    exons = [name]
    starts = gene[9].split(",")
    del(starts[-1])
    starts = strToInt(starts)
    #print (starts)
    ends = gene[10].split (",")
    del(ends[-1])
    ends = strToInt(ends)
    #print (ends)
    for i in range (0,len(starts)):
        exons.append([starts[i],ends[i]])
    #print (exons)
    lab[chr].append(exons)

print ("get the lab")
#print (lab["chr21"][0])
#print (lab["chr22"][0])
#print (lab["chrX"][0])
#print (lab["chrM"][0])

#lab = {}
#lab[chr][]= [geneName,[start1,end1],[start2,end2],[start3,end3],.......]

samFile = "./SRR1793917_5MReads.sam"
samInput = open(samFile,"r")
RNA = samInput.readlines()

output = open("./5.7_output.txt","w")
output.write("NAME    RPKM    TPM\n")


global mappable_reads
mappable_reads = 0

global read_length

for i in range (0,len(RNA)):
    if RNA[i][0] != "@":
        RNA[i]=RNA[i].split("\t")
        if RNA[i][1] != "4" and RNA[i][2] in pool:
            mappable_reads += 1
            chr = numTochr(RNA[i][2])
            location = int(RNA[i][3])
            read_length =len(RNA[i][8])
            obj[chr].append(location)
print ("get the obj")


mappable_reads = mappable_reads/1000000.00


#print (obj["chrM"])

flag = 1.0

for j in range(0,len(chromes)):
    chrome = chromes[j]
    for p in range(0,len(lab[chrome])):
        gene = lab[chromes[j]][p]
        name = gene[0]
        rpkm= RPKM (chrome, gene)
        if rpkm == 0:
            continue
        TPM = rpkm * totalNum(obj) * read_length / (1000.00 * totalNum(lab))
        print(name + "\t" + str(rpkm) + "\t" + str(TPM) + "\n")
        output.write (name + "\t" + str(rpkm) + "\t" + str(TPM) + "\n")
        flag += 1.0
        print("flag: " + str(flag / totalNum(lab) * 1000) + "  /1000 ")



'''
RPKM1 = RPKM("chr9",lab["chr9"][14])
print(RPKM1)

RPKM2 = RPKM("chr6",lab["chr6"][14])
print(RPKM1)

RPKM3 = RPKM("chr5",lab["chr5"][14])
print(RPKM1)

RPKM4 = RPKM("chr11",lab["chr11"][14])
print(RPKM1)

RPKM5 = RPKM("chr19",lab["chr19"][14])
print(RPKM1)




flag = 1
for chrome in chromes:
        for gene in lab[chrome]:
            name = gene[0]
            print(chrome)
            print(gene)
            #print(RPKM(chrome,gene))
            RPKM = RPKM(chrome,gene)
            if RPKM == 0:
                continue
            TPM = RPKM * totalNum(obj) * read_length / (1000.00 * totalNum(lab))
            #print (RPKM)
            #print (TPM)
            output.write(name + "\t" + str(RPKM) + "\t" + str(TPM) + "\n")
            flag += 1
            print (flag)
'''











