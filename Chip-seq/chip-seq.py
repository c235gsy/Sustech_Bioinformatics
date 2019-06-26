import xlrd
import xlwt

gene =open('./gene.txt',"r")
pre_lab = gene.readlines()
pre_lab.pop(0)
for i in range(0,len(pre_lab)):
    pre_lab[i] = pre_lab[i].split("\t")
    chr = pre_lab[i][2]
    if "_" not in chr:
        pre_lab[i] = [chr,int(pre_lab[i][6]),int(pre_lab[i][7]),pre_lab[i][1],"ref"]
lab = {}

for gene in pre_lab:
    if gene[0] in lab:
        lab[gene[0]].append(gene[1:])
    else:
        lab[gene[0]] = []


ExcelFile = xlrd.open_workbook('/Users/guosiyuan/PycharmProjects/Chip-seq/GSM1246686_GA4961_MEF_H3K4me3_peaks.xlsx')
sheet=ExcelFile.sheets()[0]
ExcelFile_new = xlwt.Workbook()
sheet_new = ExcelFile_new.add_sheet('sheet 1')
sheet_new.write(0,0,"mapped gene")



chromesome = sheet.col_values(0)
start = sheet.col_values(1)
end = sheet.col_values(2)
output=[]

for i in range (1,len(start)):
    lab[chromesome[i]].append([int(start[i])-000,int(end[i])+000,i,"obj"])

for chrm in lab.keys():
    lab[chrm].sort()
    for m in range (0,len(lab[chrm])):
        if lab[chrm][m][-1] == "obj":
            #print (lab[chrm][m])
            a = m
            b = m

            while lab[chrm][a][-1] != "ref":
                if a == 0:
                    a = False
                    break
                a = a - 1

            while lab[chrm][b][-1] != "ref":
                if b == len(lab[chrm]) - 1:
                    b = False
                    break
                b = b + 1

            if a != False:
                if lab[chrm][m][0] <= lab[chrm][a][1]:
                    output.append([lab[chrm][m][-2],lab[chrm][a][-2]])

            if b != False:
                if lab[chrm][m][1] <= lab[chrm][b][1] and lab[chrm][m][1] >= lab[chrm][b][0]:
                    output.append([lab[chrm][m][-2],lab[chrm][b][-2]])


#print (output)
print (len(output))
#for m in output:
    #print(m)

p = 0
for i in range(0,len(output)-1):
    if output[i][0] == output[i+1][0]:
        p += 1
        print (output[i],output[i+1])
        output[i][0] = -1
print (p)

for result in output:
    if result[0] != -1:
        sheet_new.write(result[0], 0, result[1])


ExcelFile_new.save ("Mapped_GSM1246686_GA4961_MEF_H3K4me3_peaks.xls")
#sheet_new.write

#print(lab["chrX"])


#output= open("./output_DAVID.txt","w")

