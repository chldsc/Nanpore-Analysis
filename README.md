Analyzed Data
* WT-ITR (no PmlI)
dBC-ITR (no PmlI)
hybrid-ITR (no PmlI)
Software Used
Minimap2: 
Long-read DNA/mRNA sequence aligner for fast and accurate alignment of Pacbio and Nanopore reads
Inputs fasta/fastq files and outputs bam file
Integrative Genomics Viewer (IGV): 
Genomic data visualization tool 
 Displays coverage histogram 
Describes number of reads aligned to a specific location
Libraries Used
Samtools: 
Sequence data manipulation tools for command line post-processing/analysis
BAM and SAM sequence read alignment formats
Sorts and indexes alignment file for faster data analysis
Pysam:
Python module for writing scripts for analyses with sequence data
.AlignmentFile() formats bam file into processable format
.fetch() allows access to information such as read lengths, individual base pairs, etc.
Pipeline 
Align in Minimap2
Adjust and threshold with samtools on command line
Process through analysis pipeline in jupyter notebook with pysam
Visualize in IGV and jupyter notebook using matplotlib


Current Analyses: 
Looking at the sense and anti-sense strands
