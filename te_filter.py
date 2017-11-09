#a program that filters out all TEs of ref_id 1
F = open("TE_annotations_sorted.gff3","r")
W = open("TE_1_sorted.gff3","w")

while True:
    line = F.readline()
    if not line:
        break
    if line.split("\t")[0] == "1":
        W.write(line)
F.close()
W.close()
