# TESecondaryStructure
Scripts for determining secondary structure among feature types

### Secondary Structure Analysis and MFE Statistics
1. Seperate annotation by featuretypes, and generate unique feature id's, and generate bed files.

### siRNA Mappings
1. Download siRNA libraries
2. Trim using sbatch_trimmomatic.sh
3. Extract 21, 22, 24 length reads sbatch_extract_readsonlen.sh
4. Convert 21, 22, 24 fasta's to unique fasta files with counts added to header. sbatch_generate_ucountsfasta
5. Build reference fasta file index. sbatch_bowtiebuild.sh
6. Align exact matches to reference fasta file. sbatch_bowtie2.sh
7. Convert sam to bam files. sbatch_samtools.sh
8. Filter out mapped and unmapped alignments. sbatch_samtools_filter.sh
9. Sort and index mapped files. sbatch_samtools_sort.sh
9. Convert bam file to bed file and compress. sbatch_bam2bed.sh
10. Sort bed files. sbatch_sortbed.sh
11. Generate coverage files. sbatch_bedtoolscoverage.sh
12. Generate count files. sbatch_processcovfiles.sh
13. Generate data tables for each annotation feature type. sbatch_processcounts.sh
