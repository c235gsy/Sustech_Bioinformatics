# -*- coding: UTF-8 -*-

fileName='./GA7411_TGACCA_L007_R1_001.fasta'
file = open(fileName,"r")
output = "./GA7411_TGACCA_L007_R1_001.fasta.output.txt"
outFile = open(output,"w")

name = []
read = []

i=1
for line in file:
    if i % 2 == 0 and i <= 20:
        read.append(line)
    if i % 2 == 1 and i <= 20:
        name.append(line)
    i = i + 1

file.close()

def getNewRow(row,length):
    newRow = [None] * length
    newRow[0] = row[length-1]
    for i in range (0,length-1):
        newRow[i+1]=row[i]
    return  newRow

def BWT(seq,NAME):
    seq = str(seq)
    seq = "$" + seq
    seq = seq.strip ('\n')
    length = len(seq)
    row = [None] * length
    result = ""
    for i in range (0,length):
        row[i]=seq[i]
    #print row

    matrix = [None] * length
    matrix[0] = row

    outFile.write ("\n序列编号：\n")
    outFile.write (NAME)
    outFile.write ("\n序列：\n")
    outFile.write (seq)

    for n in range(0,length-1):
        matrix[n+1]=getNewRow(matrix[n],length)

    outFile.write ("\n未排列的的矩阵：\n")

    for p in matrix:
        outFile.write ('\n')
        outFile.write (''.join (p).replace ("\r\n", " "), )
        
    matrix.sort()

    outFile.write ("\n排列后的的矩阵：\n")
    for rows in matrix:
        key = rows[-1]
        result = result + key
        outFile.write (''.join(rows).replace("\r\n"," "),)
        outFile.write ("\n")

    outFile.write ("\n结果：\n")
    outFile.write (result)


#BWT("BALALA","BALALA")

for i in range (0,len(read)):
    read[i]=str(read[i])
    name[i]=str(name[i])
    BWT(read[i],name[i])

outFile.close()


