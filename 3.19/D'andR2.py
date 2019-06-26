#!/usr/bin/python
import time
from collections import Counter
def get_D_point_and_R2(PA, PB, Pa, Pb, PAB) :
    D = PAB - PA * PB
    D_point = 0
    if D > 0:
        Dmax = min(PA * Pb, Pa * PB)
    if D < 0:
        Dmax = max(-PA * PB, -Pa * Pb)
    if Dmax == 0:
        print "the Dmax == 0, so the D' does not exist "
    D_point = D / Dmax

    R2 = D * D / (PA * PB * Pa * Pb)

    return [D_point,R2]

startTime = time.time()

with open("./lecture4_data_15.txt") as data:##读取SNPs的文件，准备抓取三条SNPs序列
    data = data.readlines()
    line = [None]*3

for i in range (0,3):##抓取任意的三条SNPs序列
    line[i] = data[i+9].split()

Haplotype = []#构造列表来储存取得的haplotypes

for i in range (1,len(line[0])):#通过寻找纯合子来获得haplotypes
    X = line[0][i]
    Y = line[1][i]
    Z = line[2][i]
    if X[0]==X[1] and Y[0]==Y[1] and Z[0]==Z[1]:
        Haplotype.append(X[0] + Y[0] + Z[0])
        Haplotype.append(X[0] + Y[0] + Z[0])

count = Counter(Haplotype) #统计从纯合子中获得的haplotyoes并且计数他们的出现次数
Haplotype = count.keys() #同字典中获得haplotypes

new_lines = [[],[],[]]#构造新的二维列表来储存phased的SNPs

for m in range(0,3):#将新的二维列表的第一个列写入各行SNPs的编号
    new_lines[m].append(line[m][0])

#print Haplotype
#print len(line[0])

for i in range (1,len(line[0])):
    threeSNPs = [list(line[0][i]),list(line[1][i]),list(line[2][i])]#每一次取二维列表中的一列三个SNPs
    for k in range (0,len(Haplotype)):#在取得的haplotypes里一个一个与当前三个SNPs比对
        element = Haplotype[k]
        new_SNPs = [[],[],[]]#新的SNPs如果生成，会被储存在这里
        for j in range (0, 3):##一对一对的比较SNPs和haplotypes
            X = threeSNPs[j]
            x = [None] * 2
            if X[0] == element[j]:
                #print "X[0] == element[i]:"
                x[0] = X[0]
                x[1] = X[1]
                #print x
                new_SNPs[j].append (x[0])
                new_SNPs[j].append (x[1])
                continue
            if X[1] == element[j]:
                #print "X[1] == element[i]:"
                x[0] = X[1]
                x[1] = X[0]
                #print x
                new_SNPs[j].append (x[0])
                new_SNPs[j].append (x[1])
                continue
            if X[0] != element[j] and X[1] != element[j]:
                break
        ##如果对应位置的SNPs和haplotypes对应位置对应，那么排序后的SNPs会被加入到列表中
        #print new_SNPs
        if len(new_SNPs[0]) == 2 and len(new_SNPs[1]) == 2 and len(new_SNPs[2]) == 2:
            for p in range(0,3):
                new_lines[p].append(new_SNPs[p])##如果新的序列每一个位置都是满的，那么就会被加入到对应的列表中
            break##当找到了符合条件的haplotypes，跳出循环找下一组SNPs
#for line in new_lines:
    #print line
    #print len(line)

listOfhaplotype = []  #构造list来储存三条序列组合后的每一个haplotype
listOfgenotpe = []  #构造list来储存三条序列组合后的每一个genotype

for b in range (1,41): ##在上述两个list中存入对应的信息
    genotype = [[],[]]
    listOfhaplotype.append (new_lines[0][b][0] + new_lines[1][b][0] + new_lines[2][b][0])
    listOfhaplotype.append (new_lines[0][b][1] + new_lines[1][b][1] + new_lines[2][b][1])
    genotype = [listOfhaplotype[2*b-2],listOfhaplotype[2*b-1]]
    #print genotype
    listOfgenotpe.append(genotype)

#print listOfhaplotype
#print len(listOfhaplotype)
#print listOfgenotpe
#print len(listOfgenotpe)

count_listOfhaplotype = Counter(listOfhaplotype)#统计各个haplotypes出现的次数

Tow_listOfhaplotype = []#取三条序列前两行的信息进行计算
Tow_haplotyoe = []
Tow_listOfgenotype = []

for ht in listOfhaplotype: #写入前两行所有的haplotypes
    ht = ht[:-1]
    Tow_listOfhaplotype.append(ht)

for gt in listOfgenotpe:#写入前两行的genotypes
    GT = []
    gt[0] = gt[0][:-1]
    gt[1] = gt[1][:-1]
    GT.append(gt[0])
    GT.append(gt[1])
    Tow_listOfgenotype.append(GT)

for HT in Haplotype:#写入前两行的haplotypes
    HT = HT[:-1]
    Tow_haplotyoe.append(HT)

Tow_haplotyoe = list(set(Tow_haplotyoe))#对这个list去重
count_listOfhaplotype = Counter(Tow_listOfhaplotype)#统计前两行所有的haplotypes和其出现个数

print Tow_listOfhaplotype
print Tow_haplotyoe
print count_listOfhaplotype

listOfP = {}#构造一个字典来储存各个haplotypes以及对应的概率
for m in range (0,4):#由上面获得的出现个数来统计频率
    mark = count_listOfhaplotype[Tow_haplotyoe[m]]
    #print mark
    listOfP[Tow_haplotyoe[m]]= mark/40.00

print listOfP
#print Tow_listOfgenotype

num_AB = 0
for geno in Tow_listOfgenotype:#随机选取基因型并计算出现的频数
    #print geno
    if geno[0] == Tow_haplotyoe[3] and geno[1] == Tow_haplotyoe[1]:
        num_AB = num_AB + 1
    if geno[1] == Tow_haplotyoe[3] and geno[0] == Tow_haplotyoe[1]:
        num_AB = num_AB + 1

PAB = num_AB/40.00#计算频率
print PAB

PA = listOfP[Tow_haplotyoe[3]]
PB = listOfP[Tow_haplotyoe[1]]
Pa = listOfP[Tow_haplotyoe[2]]
Pb = listOfP[Tow_haplotyoe[0]]

result = get_D_point_and_R2(PA , PB , Pa , Pb , PAB)#计算想要的数据结果

print "D': " , result[0]
print "R_power: " , result[1]

with open("./output.txt","w") as output:##写入对应的结果到文件

    output.write ("Unphased SNPs:\n")
    for old_line in line:
        output.write (old_line[0])
        for n in range (1, len (old_line)):
            output.write ("  ")
            output.write ("".join (old_line[n]))
        output.write ("\n")

    output.write("\nPhased SNPs:\n")
    for new_line in new_lines:
        output.write (new_line[0])
        for n in range (1,len(new_line)):
            output.write ("  ")
            output.write ("".join(new_line[n]))
        output.write("\n")

    output.write("\nHaplotypes:\n")
    for HT in Haplotype:
        output.write(HT+" ")

    output.write("\n")

    output.write("\n"+"Tow phased SNPs: \n"+"  ".join(Tow_listOfhaplotype))
    output.write ("\n")
    output.write("\n"+"Their haplotypes: \n"+"  ".join(Tow_haplotyoe))
    output.write ("\n\n")
    output.write("the number of these haploytpes: " + str(count_listOfhaplotype))
    output.write ("\n\n")
    output.write("the P of these haplotypes: " + str(listOfP))
    output.write ("\n\n")
    output.write("the number of %s and %s : %s" %(Tow_haplotyoe[3],Tow_haplotyoe[1],num_AB))
    output.write ("\n\n")
    output.write("the P of %s and %s : %s" %(Tow_haplotyoe[3],Tow_haplotyoe[1],PAB))
    output.write ("\n\n")
    output.write("D': " + str(result[0]) + "   " + "R_power: " + str(result[1]))

endTime = time.time()
print endTime - startTime