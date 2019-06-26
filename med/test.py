
def add_read_to_LIST(count,read,quality,LIST):
    #if genome[count-1:count+8] == read[0:9] and genome[count+len(read)-5:count+len(read)-1] == read[-5:-1]:
    for p in range(0,len(read)):
        LIST[count+p].append([read[p],quality[p]])
    return LIST

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

sam_file = open("./Reads.bwt.sam")
sam = sam_file.readlines()

genome_file = open("./EcoliGenome.fa")
genome = genome_file.readlines()
genome.pop(0)
genome = "".join(genome)
genome = genome.replace("\n","")
print (genome)
LIST = [""]
for c in range (1,len(genome)):
    LIST.append([["ref:",genome[c]]])

print (LIST[1])

Reads = []
NewReads = []
for i in range (0,len(sam)):
    if sam[i][0] != "@":
        Reads.append(sam[i].split("\t"))

for read in Reads:
    NewReads.append([read[3],read[9],read[10]])

for newRead in NewReads:
    print (newRead)
    LIST = add_read_to_LIST(int(newRead[0]),newRead[1],newRead[2],LIST)

result = []

for m in range(1,len(LIST)):
    if len(LIST[m]) >= 2:
        result.append(getResult(m,LIST[m]))

output1 = open("./output.pileup","w")
for re in result:
    output1.write( re + "\n")
output1.close()

print (len(result))

for a in range(0,len(result)):
    result[a]=result[a].split(" ")
    print (result[a])

output2 = open("./output.vcf","w")
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

