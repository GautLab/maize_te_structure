#!/bin/bash
#module load samtools/1.3

TEGFF3="B73v4.TE.filtered.gff3"
REFASM="Zea_mays.AGPv4.dna.toplevel.fa"

#extract types while ignoring lines that start with #'s
awk '!/^#/{print}' ${TEGFF3} | awk '{print $3}' | sort | uniq > $(basename ${TEGFF3} .gff3).types

#generate seperate files based on types
for i in $(cat $(basename ${TEGFF3} .gff3).types)
do
    awk -v awki=${i} '$3 == awki { print $0 }' ${TEGFF3} > $(basename ${TEGFF3} .gff3)_${i}.gff3
    #the first statement ignores lines that start with #, 2nd statement formats a gff3 file into bed format.
    awk '!/^#/{print}' $(basename ${TEGFF3} .gff3)_${i}.gff3 | awk '{print $1"\t"$4"\t"$5"\t"$3"::"$1":"$4-1"-"$5"\t42\t"$7"\t"$5-$4+1"M"}' > $(basename ${TEGFF3} .gff3)_${i}.bed
done

#create fasta files for each feature type
#samtools faidx ${REFASM}
#echo "te file generated"

#for i in $(awk '{print $1":"$2"-"$3}' TE_posonly.tsv)
#  do samtools faidx ${REFASM} ${i} >> maizetes.fasta
#done
#echo "indexes created"

