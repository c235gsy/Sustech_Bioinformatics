import xlrd


ExcelFile = xlrd.open_workbook('./Human_CD4_stim_refGene_RPKM_Count_FDR.xlsx')
sheet=ExcelFile.sheets()[0]
names=sheet.col_values(0)
names=names[1:4841]
output= open("./output_DAVID.txt","w")


for i in range (0,len(names)):
    if "+" in str(names[i]):
        names[i] = ""
    else:
        names[i] = str(names[i])+"\n"

print(len(names))

output.writelines(names[0:2000])


