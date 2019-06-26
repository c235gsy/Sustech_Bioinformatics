import os
import os.path
import gzip


strZipFile = './GA7411_TGACCA_L007_R1_001.fastq.gz'## 选择输入文件
strDstFile = './GA7411_TGACCA_L007_R1_001.fasta'    ## 选择输出文件

file = gzip.GzipFile (strZipFile, "r") ##解压gzip文件得到内容
NEWfile = []##构造新的输入源文件

i = 1
for lines in file:##取输入文件中所有的行
    if i % 4 == 1 or i % 4 == 2:##如果是每四行中的第一第二行
        NEWfile.append(lines) ## 则加入新的文件中
        #print lines
    i=i+1


outFile = open (strDstFile, "w ") ##打开输出文件

for NEWlines in NEWfile:##取新文件中所有的行
    outFile.write (NEWlines)##按顺序写入文件中


outFile.close ()##关闭文件

print 'Finished!'


