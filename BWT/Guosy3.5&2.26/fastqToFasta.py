import os
import os.path
import gzip


strZipFile = './GA7411_TGACCA_L007_R1_001.fastq.gz'
strDstFile = './GA7411_TGACCA_L007_R1_001.fasta'

file = gzip.GzipFile (strZipFile, "r")
NEWfile = []

i = 1
for lines in file:
    if i % 4 == 1 or i % 4 == 2:
        NEWfile.append(lines)
        #print lines
    i=i+1


outFile = open (strDstFile, "w ")

for NEWlines in NEWfile:
    outFile.write (NEWlines)


outFile.close ()

print 'Finished!'


