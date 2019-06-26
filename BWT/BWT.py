# -*- coding: UTF-8 -*-

fileName='./GA7411_TGACCA_L007_R1_001.fasta'
file = open(fileName,"r")
output = "./GA7411_TGACCA_L007_R1_001.fasta.output.txt"
outFile = open(output,"w")
##读取文件
name = []
read = []
##构造空的列表来储存序列和序号
i=1
for line in file:
    if i % 2 == 0 and i <= 20:
        read.append(line)
    if i % 2 == 1 and i <= 20:
        name.append(line)
    i = i + 1
##选择10条序列和他们的序号分别加入之前构造的列表中
file.close()
##筛取序列和编号

def getNewRow(row,length): #输入序列和它的长度
    newRow = [None] * length #构造新的空的行
    newRow[0] = row[length-1] # 新行的第一个字符等于输入行的最后一个字符
    for i in range (0,length-1):
        newRow[i+1]=row[i]##新行的每一个字符相当于输入行的其他字符向后移动一位
    return  newRow  #返回新的行
##把所有的字符向后挪动一位

def BWT(seq,NAME):
    seq = str(seq)  ##将序列转化为字符串
    seq = "$" + seq ##在序列的前面加上 $ 符号
    seq = seq.strip ('\n') ##去掉序列自带的换行符
    length = len(seq)  ##得到序列的长度
    row = [None] * length  ##构造空的列表来储存序列里面的每一个字符
    result = "" ## 处理输入的序列并测量长度，定义变量result
    for i in range (0,length):
        row[i]=seq[i] ##拆分输入序列
    #print row

    matrix = [None] * length
    matrix[0] = row ##构造矩阵并定义第一行

    outFile.write ("\n序列编号：\n")
    outFile.write (NAME)
    outFile.write ("\n序列：\n")
    outFile.write (seq)  ##按格式写入取得的序列及其序号

    for n in range(0,length-1):
        matrix[n+1]=getNewRow(matrix[n],length)##用之前定义的仿佛来写剩下的行

    outFile.write ("\n未排列的的矩阵：\n")

    for p in matrix:
        outFile.write ('\n')
        outFile.write (''.join (p).replace ("\r\n", " "), )
        ##将矩阵中的每一行都合并成字符串并去换行符和空格，按原始顺序写入文件
    matrix.sort()
        ##将矩阵中的每一行都按第一个字符的字符表顺序排列，得到变换后的矩阵
    outFile.write ("\n排列后的的矩阵：\n")
    for rows in matrix:
        key = rows[-1] ##输出的结果需要取每一行的最后一个元素
        result = result + key ##迭代得到最终结果
        outFile.write (''.join(rows).replace("\r\n"," "),)
        outFile.write ("\n") ##将变换后的矩阵的每一行都合并成字符串并去换行符和空格，按顺序写入文件

    outFile.write ("\n结果：\n")
    outFile.write (result) ##将变换后的结果写入文件


#BWT("BALALA","BALALA")

for i in range (0,len(read)):##按顺序取文件中的序列
    read[i]=str(read[i]) ##取序列的序列
    name[i]=str(name[i]) ##取序列的序号
    BWT(read[i],name[i]) ##带入方法中

outFile.close()


