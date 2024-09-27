#!/bin/bash
head -n 5 Zea_mays.AGPv4.39.chr_snRNA.gff3 | awk '{print $1"\t"$4"\t"$5"\t"$3"::"$1":"$4-1"-"$5"\t42\t"$7"\t"$5-$4+1"M"}'
#head -n 5 Zea_mays.AGPv4.39.chr_snRNA.gff3 | awk '{print $1"\t"$4"\t"$5"\t"$3"::"$1":"$4-1"-"$5"\t42\t"$7"\t"$5-$4+1"M"}'
