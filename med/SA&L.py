
#["A",num,[A,C,G,T]]

genome_file = open("./EcoliGenome.fa")
genome = genome_file.readlines()
head = genome[0]
genome.pop(0)
genome = "".join(genome)
genome = genome.replace("\n","")
#print (genome)

out = open("./re_L_SA_count.txt","w")

len_genome = len(genome)
genome = genome[::-1]

added_genome = "$" + genome + "$" + genome[0:48]

L_and_SA = []
for i in range(0,len_genome+1):
    L_and_SA.append([added_genome[-49-i:-1-i]+added_genome[-50-i],str(len_genome-i)])

#for line in L_and_SA:
    #print(line)
L_and_SA.sort()

count =[0,0,0,0]
for line in L_and_SA:
    line[0] = line[0][-1]

    if line[0] == "A":
        count[0] = count[0] + 1
    if line[0] == "C":
        count[1] = count[1] + 1
    if line[0] == "G":
        count[2] = count[2] + 1
    if line[0] == "T":
        count[3] = count[3] + 1

    line.append (str (count[0]))
    line.append (str (count[1]))
    line.append (str (count[2]))
    line.append (str (count[3]))

    print (str(line))
    out.write(str("\t".join(line)+"\n"))

print (len(L_and_SA))
print (len_genome)