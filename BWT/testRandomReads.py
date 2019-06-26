import os
import os.path
import gzip
import random

strZipFile = './GA7411_TGACCA_L007_R1_001.fastq.gz'##选择输入文件
strDstFile = './Random_Guosiyuan.fastq'##选择输出文件
strOtherFile = './no_of_reads.txt'##选择输出所选随机序列的编号的文件

file = gzip.GzipFile (strZipFile, "r")##解压输入文件
NEWfile = []##构造文件用于输出到输出文件
FileList = []##构造列表来储存输入文件中的每一行
Dictionary = {}##构造字典来储存每个read的4行并给他们编号

for lines in file:##构造列表来储存输入文件中的每一行
    lines = str(lines)##格式化为字符串
    lines = lines.strip ('\n')##去掉换行符
    FileList.append(lines)##加入到列表中

#print FileList[5:500]

for n in range (3,len(FileList),4):#在列表中每隔四行选取一行
    fastq = FileList[n-3:n+1] ##选取与被选取行相关的四行
    Dictionary[(n+1)/4] = fastq ##将这四行加入字典中，key值为编号

#print Dictionary[500]

#print len(Dictionary) ## for test


noOfReads = open(strOtherFile,'w')##打开选择的随机序列的序列号的输出文件

i = 1
while i <= 500:##进行500次循环
    p = random.randint(1,len(Dictionary)+1)##在字典中所有的key值中随机选择
    noOfReads.write (str(p))##将改key值写入序号文件中
    noOfReads.write ("\n")
    #print "\n"
    #print Dictionary[p]
    NEWfile.append( Dictionary[p] )##把key值对应的四行作为一个列表加入到输出文件中
    i=i+1


outFile = open ( strDstFile, "w")##打开输出文件

for NEWreads in NEWfile:##取输出文件中的每一个小列表
    for NEWlines in NEWreads:##取每一个小列表中的每一行
        outFile.write ("\n"+NEWlines)##一行一行的写入文件中


outFile.close ()##关闭输出文件

