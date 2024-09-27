#!/bin/bash

source /networkshare/.mybashrc

GFF="Zea_mays.AGPv4.39.chr.gff3"

awk '$3 == "CDS" { print $0 }' ${GFF} > $(basename ${GFF} .gff3)_CDS.gff3
#awk '$3 == "chromosome" { print $0 }' ${GFF} > $(basename ${GFF} .gff3)_chromosome.gff3
awk '$3 == "exon" { print $0 }' ${GFF} > $(basename ${GFF} .gff3)_exon.gff3
awk '$3 == "five_prime_UTR" { print $0 }' ${GFF} > $(basename ${GFF} .gff3)_five_prime_UTR.gff3
awk '$3 == "gene" { print $0 }' ${GFF} > $(basename ${GFF} .gff3)_gene.gff3
awk '$3 == "lnc_RNA" { print $0 }' ${GFF} > $(basename ${GFF} .gff3)_lnc_RNA.gff3
awk '$3 == "miRNA" { print $0 }' ${GFF} > $(basename ${GFF} .gff3)_miRNA.gff3
awk '$3 == "mRNA" { print $0 }' ${GFF} > $(basename ${GFF} .gff3)_mRNA.gff3
awk '$3 == "ncRNA_gene" { print $0 }' ${GFF} > $(basename ${GFF} .gff3)_ncRNA_gene.gff3
awk '$3 == "pre_miRNA" { print $0 }' ${GFF} > $(basename ${GFF} .gff3)_pre_miRNA.gff3
awk '$3 == "rRNA" { print $0 }' ${GFF} > $(basename ${GFF} .gff3)_rRNA.gff3
awk '$3 == "snoRNA" { print $0 }' ${GFF} > $(basename ${GFF} .gff3)_snoRNA.gff3
awk '$3 == "snRNA" { print $0 }' ${GFF} > $(basename ${GFF} .gff3)_snRNA.gff3
awk '$3 == "SRP_RNA" { print $0 }' ${GFF} > $(basename ${GFF} .gff3)_SRP_RNA.gff3
awk '$3 == "three_prime_UTR" { print $0 }' ${GFF} > $(basename ${GFF} .gff3)_three_prime_UTR.gff3
awk '$3 == "tRNA" { print $0 }' ${GFF} > $(basename ${GFF} .gff3)_tRNA.gff3

bash bash_bedtools.sh Zea_mays.AGPv4.39.chr.types 4.39 .chr
