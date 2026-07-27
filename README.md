# Analyzed Data
* WT-ITR (no PmlI)
* dBC-ITR (no PmlI)
* hybrid-ITR (no PmlI)
  
# Major Tools/Software Used
Minimap2: 
* CLI level Long-read DNA/mRNA sequence aligner for fast and accurate alignment of Pacbio and Nanopore reads
* Inputs fasta/fastq files and outputs bam file
  
Integrative Genomics Viewer (IGV): 
* Genomic data visualization tool 
* Displays coverage histogram 
* Describes number of reads aligned to a specific location

Libraries Used
Samtools: 
* Sequence data manipulation tools for command line post-processing/analysis
* BAM and SAM sequence read alignment formats
* Sorts and indexes alignment file for faster data analysis

Pysam:
* Python module for writing scripts for analyses with sequence data
* .AlignmentFile() formats bam file into processable format
* .fetch() allows access to information such as read lengths, individual base pairs, etc.

Pipeline 
1. Align in Minimap2
2. Adjust and threshold with samtools on command line
3. Process through analysis pipeline in jupyter notebook with pysam
4. Visualize in IGV and jupyter notebook using matplotlib

Environment Setup: 
- Install miniconda and install libraries through environment.yml file
- Install IGV following instructions on https://igv.org/doc/desktop/#QuickStart/

Current Analyses: 
Looking at the sense and anti-sense strands
