#!/bin/bash
#processing order of operations
#trim fastqfiles then extract 21,22 & 24 lenghts

#convert 21,22,24 lenghts to ucounts with # of copies then generate fasta file
qsub qsub_convertucounts2fa.sh #calls covertucounts2fa.py
#align libs to reference
qsub qsub_bowtiebuild.sh
qsub qsub_bowtie2.sh
qsub qsub_samtools.sh
#get a list of unique feature types
cut -f3 Zea_mays.AGPv4.39.chr.gff3 | sort | uniq > Zea_mays.AGPv4.39.chr.types
cut -f3 B73v4.TE.filtered.noheaders.gff3 | sort | uniq > B73v4.TE.filtered.noheaders.types
bash bash_seperate_zmays4.39.sh
bash bash_seperate_zmaystes.sh
bash bash_extract_TEs.sh
#Look for duplicate features. If numbers dont match then examine and reduce.
#for i in $(ls *.bed); do echo $i; awk '{print $4}' $i | wc -l; awk '{print $4}' $i | sort | uniq | wc -l; done
#awk '{print $4}' B73v4.TE.filtered_SINE_element.bed | sort | uniq -d
