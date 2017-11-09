#same as blast_sort, but for the te file by columns 3 and 4
import radix_sort
import sys

def buildlist(L,filename):
    F = open(filename,"r")
    for line in F:
        splitline = line.split("\t")
        index = min(int(splitline[3]),int(splitline[4]))
        L.append((index,line))
    F.close()
    return

def writefile(L,new_filename):
    F = open(new_filename,"w")
    for entry in L:
        F.write(entry[1])

def sort_file(filename, newfilename):
    L = []
    print("Initializing list")
    buildlist(L,filename)
    print("List Initialized. Length:",len(L))
    print("Sorting List")
    radix_sort.radix_sort(L)
    print("List sorted")
    print("Writing outfile")
    writefile(L,newfilename)
    print("Outfile written")
    print("Done!")
    return

if len(sys.argv) == 3:
    filename = sys.argv[1]
    newfilename = sys.argv[2]
    sort_file(filename,newfilename)
elif len(sys.argv) == 2:
    filename = sys.argv[1]
    newfilename = filename.split(".")[0]+"_sorted.out"
    sort_file (filename,newfilename)
else:
    print("te sort: INCORRECT NUMBER OF ARGUMENTS")
