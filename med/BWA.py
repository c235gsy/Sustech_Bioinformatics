# -*- coding: utf-8 -*-
import time


def get_number(char):
    if char == "A":
        num = num_A
    if char == "C":
        num = num_C
    if char == "G":
        num = num_G
    if char == "T":
        num = num_T
    return num
def get_c(char):
    if char == "A":
        c = 2
    if char == "C":
        c = 3
    if char == "G":
        c = 4
    if char == "T":
        c = 5
    return c
def location_to_no(char,lo):
    if char == "A":
        no = lo
    if char == "C":
        no = lo - num_A
    if char == "G":
        no = no - num_A - num_C
    if char == "T":
        no = no - num_A - num_C - num_G
    return no
def no_to_location(char,no):
    if char == "A":
        lo = no
    if char == "C":
        lo = no + num_A
    if char == "G":
        lo = no + num_A + num_C
    if char == "T":
        lo = no + num_A + num_C + num_G
    return lo
def seq_score(seq1,seq2):
    score = 0.00
    for i in range (0,len(seq1)):
        if seq1[i] == seq2[i]:
            score += 1
    return score/len(seq1)
def add_read_to_LIST(count,read,quality,table):
    #if genome[count-1:count+8] == read[0:9] and genome[count+len(read)-5:count+len(read)-1] == read[-5:-1]:
    for p in range(0,len(read)):
        table[count+p].append([read[p],quality[p]])
    return table

def getResult(location,data):
    vir = ""
    qua = ""
    for n in range(1,len(data)):
        qua = qua + data[n][1]
        if data[n][0] == data[0][1]:
            vir = vir + "."
        if data[n][0] != data[0][1]:
            vir = vir + str(data[n][0])
    result = "EcoliGenome "+str(location)+" "+data[0][1]+" "+str(len(data)-1)+" "+vir+" "+qua
    return result

def score(quality):
    length  = len(quality)
    Q_sum = 0.00
    for e in quality:
        Q_sum = Q_sum + float(ord(e))-33
    return float(Q_sum/length)

def judge(score):
    if score >= 30:
        judge = "PASS"
    else:
        judge = "Fail"
    return judge

def ALT(alt):
    ALT = ""
    for a in alt:
        if a != ".":
            ALT = ALT + "," +a
    return ALT

def GT_AD_DP(data):
    num = len(data)
    output = [""]*num
    for i in range(0,num) :
        if data[i] == ".":
            output[i] = output[i] + "0:"
        else:
            output[i] = output[i] + "1:"

# lo 表示在table中的位子， no 表示在同一碱基中的顺序，从table中得到的都是 no
def BWA(seq):
    no_up = get_number(seq[-1])
    no_down = 0
    lo_up = no_to_location(seq[-1],no_up)
    lo_down = no_to_location(seq[-1],no_down)
    #print([seq[-1],lo_up,lo_down,lo_up-lo_down])
    result = []

    i = -2
    while i >= -len(seq):
        element = seq[i]
        c = get_c(element)
        no_up = int(table[lo_up][c])
        no_down = int(table[lo_down][c])
        #print ([element,no_up,no_down,no_up-no_down])
        lo_up = no_to_location (element, no_up)
        lo_down = no_to_location (element, no_down)

        if no_up - no_down == 1:
            SA = int(table[lo_up][1])
            seq2 = genome[SA-i+1-len(seq):SA-i+1]
            #print ("SA",SA)
            #print (seq)
            #print (seq2)
            flag = True
            score = seq_score (seq, seq2)
            if score <= 0.5:
                flag = False
            result = [flag,SA - i + 1 - len(seq),len(seq),score,i]
            break

        if no_up - no_down == 0:
            SA = int (table[lo_up][1])
            seq2 = genome[SA - i + 1 - len(seq):SA - i + 1]
            result = [False, SA - i + 1 - len(seq),len(seq),seq_score(seq, seq2), i]
            break
        i = i - 1

    #print (table[lo_up][1])
    return result

def TransToSam(read):
    name = read[0].strip("@")
    seq = read[1]
    Qseq = read[2]
    location = str(read[3][1])
    score = str(int(read[3][3]*100))
    CIGAT = str(read[3][2]) + "M"
    output = name+"\t"+"0\t"+genome_info[0].strip(">")+"\t"+location+"\t"+score+"\t"+CIGAT+"\t*\t0\t0\t"+seq+"\t"+Qseq
    return output


################################################################
################################################################
################################################################
start = time.time()

prepare1 = open("./L_SA_count.txt","r")
#prepare2 = open("")
global table
table = prepare1.readlines()
for l in range(0,len(table)) :
    table[l] = table[l].strip("\n")
    table[l] = table[l].split("\t")
    #print(table[l])
    #if line[0] != "$":
        #print (line)

genome_file = open("./EcoliGenome.fa",'r')
global genome
global genome_info
genome = genome_file.readlines()
genome_info = genome[0]
genome.pop(0)
genome = "".join(genome)
genome = genome.replace("\n","")
genome = "$" + genome
genome_info = genome_info.split(" ")

Reads_file = open("./EcoTestRead2.fq","r")
Reads_data = Reads_file.readlines()
for data in range(0,len(Reads_data)) :
    Reads_data[data] = Reads_data[data].replace("\n","")
#print (Reads_data)
reads = []
for i in range(0,len(Reads_data)-4,4):
    #print (i)
    reads.append([Reads_data[i],Reads_data[i+1],Reads_data[i+3]])


global num_A
num_A = int(table[-1][2])
#print(num_A)
global num_C
num_C = int(table[-1][3])
#print(num_C)
global num_G
num_G = int(table[-1][4])
#print(num_G)
global num_T
num_T = int(table[-1][5])

#print(num_T)
#print(num_A+num_C+num_G+num_T)
#print (len(genome))
#print (genome.count("A"))
#print (genome.count("C"))
#print (genome.count("G"))
#print (genome.count("T"))
#result = []
#Seq = "GATAAAGCAGGAATTACTACTGCTTGTTTCGAATTAAATCGAAGTGGACTGCTGGCGGAAAATGAGAAA"
#re_Seq  = Seq[::-1]
#print (len(Seq))
#print (len(table))

for i in range (0,len(reads)):
    #print (reads[i])
    Seq = reads[i][1]
    test = BWA(Seq)
    if test[0]:
        reads[i].append(test)
    else :
        test = BWA(Seq[:test[4]])
        reads[i].append (test)
    #print (reads[i])

sam_write = open("./output/EcoTestRead2.bwt.sam","w")
reads = sorted(reads,reverse=True)
for read in reads:
    if read[3][0] == True:
        output = TransToSam(read)
        sam_write.write(output+"\n")
sam_write.close()

LIST = [""]
for c in range (1,len(genome)):
    LIST.append([["ref:",genome[c]]])

sam_file = open("./output/EcoTestRead2.bwt.sam","r")
sam = sam_file.readlines()
print (LIST)
print (sam)
Reads = []
NewReads = []
for i in range (0,len(sam)):
    if sam[i][0] != "@":
        Reads.append(sam[i].split("\t"))

print(Reads)

for Read in Reads:
    NewReads.append([Read[3],Read[9],Read[10]])

print (NewReads)

for NewRead in NewReads:
    print (NewRead)
    LIST = add_read_to_LIST(int(NewRead[0]),NewRead[1],NewRead[2],LIST)





result = []




for m in range(1,len(LIST)):
    if len(LIST[m]) >= 2:
        result.append(getResult(m,LIST[m]))

output1 = open("./output/EcoTestRead2.bwt.pileup","w")

for re in result:
    output1.write( re + "\n")
output1.close()

print (len(result))

for a in range(0,len(result)):
    result[a]=result[a].split(" ")
    print (result[a])

output2 = open("./output/EcoTestRead2.bwt.vcf","w")
output2.write("#genome\tposition\tID\tREF\tALT\tFILTER\tINFO\n")
VCF = []
for resu in result:
    if resu[3] != "1" and resu[4] != "." * int (resu[3]):
        gene = resu[0]
        pos = resu[1]
        ID = "."
        REF = resu[2]
        ALT = []
        for g in resu[4]:
            if g != "." :
                ALT.append(g)
        #ALT = ",".join(ALT)
        FILTER = "PASS"
        if len(ALT) > 0.5 * float(resu[3]):
            FILTER = "FALSE"
        INFO = "."
        ALT = ",".join(ALT)
        output2.write (gene+"\t"+pos+"\t"+ID+"\t"+REF+"\t"+ALT+"\t"+FILTER+"\t"+INFO+"\n")

output2.close()

print (time.time()-start)