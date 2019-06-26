#!/usr/bin/python
# -*- coding: utf-8 -*-

def avg(l):
    return  (sum (l)) / len (l);

def median(l):
    l = sorted (l);  # 先排序
    if len (l) % 2 == 1:
        return l[len (l) / 2.00];
    else:
        return (l[len (l) / 2 - 1] + l[len (l) / 2]) / 2.00;

def getCG(seq):
    CG = 0
    length = float(len(seq))
    for char in seq:
        str(char)
        if char == "C" or char == "G":
            CG = CG + 1
    return float(CG)/length

def getQofRead(seq):
    length =len(seq)
    sum_Q = 0.00
    for bp in seq:
        Q = float(ord(bp))-33
        sum_Q = sum_Q + Q
    return  sum_Q/length

def getTheNumberOfElementInRange(list,floor,top):
    num = 0
    for element in list:
        if element >= floor and element <= top:
            num = num + 1
    return num

def outputTheDic(dic):
    output = []
    keys = dic.keys()
    keys.sort()
    for key in keys:
        output.append(str(key)+"    "+str(dic[key]))
    return output

inputFile = './HomoTestRead.fq'## 选择输入文件
outputFile = './output_QandGC.txt'    ## 选择输出文件
#CG_outputFile = './output_CG.txt'    ## 选择输出文件

file = open(inputFile, "r") ##解压gzip文件得到内容
file = file.readlines()
listOfseqs = []##构造新的输入源文件
listOfreads = []

i = 1
for lines in file:##取输入文件中所有的行
    if i % 4 == 2:##如果是每四行中的第一第二行
        listOfseqs.append(str(lines))## 则加入新的文件中
    if i % 4 == 0: ##如果是每四行中的第一第二行
        listOfreads.append(str(lines)) ## 则加入新的文件中
    i=i+1
####################################################
####################################################
####################################################
listOfQMean = []
listOfQMedian = []

for ii in range (0,len(listOfseqs[0])):
    listOfQ_site= []
    for read in listOfreads:
        listOfQ_site.append(int(ord(read[ii]))-33.00)

    listOfQMean.append(str(avg(listOfQ_site)))
    listOfQMedian.append(str(median(listOfQ_site)))

print listOfQMean
print listOfQMedian

outFile = open (outputFile, "w")##打开输出文件

strOfQMean = "  ".join(listOfQMean)
strOfQMedian = "  ".join(listOfQMedian)

outFile.write("strOfQMean:\n"+strOfQMean+"\n")
outFile.write("strOfQMedian:\n"+strOfQMedian+"\n")
##################################################
##################################################
##################################################
listOfQ_seq = []

for seq in listOfseqs:
    listOfQ_seq.append(getQofRead(seq))

print listOfQ_seq
Distribution = {}

for key in range(0,50,5):
    #print key
    #print getTheNumberOfElementInRange(listOfCG_seq,key,key+5)
    Distribution[str(key)+"--"+str(key+5)] = getTheNumberOfElementInRange(listOfQ_seq,key,key+5)

outFile.write("\nDistribution of quality score: \n")

for line in outputTheDic(Distribution):
    outFile.write(line + "\n")

#outFile.write("\nDistribution of quality score: \n"+str(Distribution)+"\n")


##################################################
##################################################
##################################################
listOfCG_seq = []
for seq in listOfseqs:
    listOfCG_seq.append((getCG(seq)))

print listOfCG_seq

dicCG = {}
for percent in range(0,100,10):
    #print key
    #print getTheNumberOfElementInRange(listOfCG_seq,key,key+5)
    dicCG[str(percent+0.00)+"% -- "+str(percent+10.00)+"%"] = getTheNumberOfElementInRange(listOfCG_seq,percent/100.00,percent/100.00+0.10)

outFile.write("\nCG distribution of over all reads: \n")
for line in outputTheDic(dicCG):
    outFile.write(line + "\n")

outFile.close()
#outFile.write("\nCG distribution of over all reads: \n"+str(dicCG)+"\n")