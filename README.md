# DNAMarkMaker

## Table of contents
 - [Introduction of modified QTL-seq](#Introduction-of-modified-QTL-seq)
 - [Installation](#Installation)
   + [Dependencies](#Dependencies)
   + [Installation using bioconda](#Installation-using-bioconda)
 - [Usage](#Usage)
   + [command : mqtlseq](#command-:-mqtlseq)
   + [command : mqtlseq_mpileup](#command-:-mqtlseq_mpileup)
   + [command : mqtlseq_snpindex](#command-:-mqtlseq_snpindex)
   + [command : mqtlseq_plot](#command-:-mqtlseq_plot)
 - [The input file format](#The-input-file-format)
 - [The example of execution](#The-example-of-format-execution)
 - [The output file format](#The-output-file-format)

## Introduction of modified QTL-seq
  
  modified QTL-seq is an improved version of QTL-seq that introduces a filter system that can be applied to plant species with heterozygous genomes. 

#### Citation
  - [Itoh et al. (2019). Next-generation sequencing-based bulked segregant analysis for QTL mapping in the heterozygous species Brassica rapa. Theoretical and Applied Genetics, 132, 2913-2925.](https://link.springer.com/article/10.1007/s00122-019-03396-z)

## Installation
### Dependencies
   - python
   - [bwa](https://github.com/lh3/bwa)
   - [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)
   - [samtools](https://github.com/samtools/samtools)
   - [bcftools](https://github.com/samtools/bcftools)

### Installation using bioconda
  You can install modified QTL-seq using bioconda.
  ```
  conda install -c bioconda modifiedqtlseq
  ```
  Alternatively, if you want to create modified QTL-seq specific environment.
  ```
  conda create -n modifiedqtlseq -c bioconda modifiedqtlseq
  conda activate modifiedqtlseq 
  ```
  
## Usage
  'modified QTL-seq' offers four commands. mqtlseq, mqtlseq_mpileup, mqtlseq_snpindex, mqtlseq_plot.

  test data can be downloaded from XXX.

### command : mqtlseq
  
  Description: This is the Standard commands that can be started from reference fasta and fq files.

  Required options
  ```
  -r FASTA                           Full path of reference sequence fasta.
  -rt TXT                            Full path of txt file discribed sample and fastq path.
  ```
  
  Additional options
  ```
  -p1 STR                            P1 name. (P1)
  -p2 STR                            P1 name. (P2)
  -a STR                             A-bulk name. (A_bulk)
  -b STR                             B-bulk name. (B_bulk)
  -o STR                             Output directry name. (QTL-seq)
  -t INT                             Number of parallel processing. (1)
  -bwa STR                           Alignment tools "bwa" or "bwa-mem2". (bwa)
  -bq INT                            Minimum base quality for calling SNP. (13)
  -mq INT                            Minimum mapping quality for calling SNP. (40)
  -mindepth INT                      Minimum depth of target SNP. (10)
  -maxdepth INT                      Maximum depth of target SNP. (99)
  -sim_replication INT               Simulation replication counts of SNP-index. (10000)
  -individual INT                    Number of bulked individuals. (20)
  -generation STR                    Bulked generation "F2", "RIL" or "BC1F1". (F2)
  -window_size INT                   Window size for sliding window analysis. (2000000)
  -step_size INT                     Step size for sliding window analysis. (50000)
  -min_plot INT                      The minimum number of SNPs in the window. (10)
  -filt_bias FLOAT                   If both bulks show set value<(SNP-index) and set value>1-(SNP-index), the position is exclueded. (0.3)
 ```

### command : mqtlseq_mpileup
  
  Description: This is the commands that can be started from P1 SNPs substituted fasta and bam data files aligned to substituted fasta.

  Required options
  ```
  -r FASTA                           Full path of P1 SNPs substituted fasta
  -bt TXT                            Full path of txt file discribed sample and bam path.
  ```
  
  Additional options
  ```
  -p1 STR                            P1 name (P1)
  -p2 STR                            P1 name (P2)
  -a STR                             A-bulk name (A_bulk)
  -b STR                             B-bulk name (B_bulk)
  -o STR                             Output directry name (QTL-seq)
  -t INT                             Number of parallel processing (1)
  -mindepth INT                      Minimum depth of target SNP. (10)
  -maxdepth INT                      Maximum depth of target SNP. (99)
  -sim_replication INT               Simulation replication counts of SNP-index (10000)
  -individual INT                    Number of bulked individuals (20)
  -generation STR                    Bulked generation "F2", "RIL" or "BC1F1" (F2)
  -window_size INT                   Window size for sliding window analysis (2000000)
  -step_size INT                     Step size for sliding window analysis (50000)
  -min_plot INT                      The minimum number of SNPs in the window (10)
  -filt_bias FLOAT                   If both bulks show set value<(SNP-index) and set value>1-(SNP-index), the position is exclueded. (0.3)
 ```

### command : mqtlseq_snpindex
  
  Description: This is the commands that can be started from vcf file called SNPs.

  Required options
  ```
  -vcf VCF                           Full path of vcf file called SNPs
  ```
  
  Additional options
  ```
  -p1 STR                            P1 name (P1)
  -p2 STR                            Specify the name if using P2 (False)
  -a STR                             A-bulk name (A_bulk)
  -b STR                             Specify the name if using B-bulk (False)
  -f1 STR                            Select True if you are using F1 (False)
  -o STR                             Output directry name (QTL-seq)
  -t INT                             Number of parallel processing (1)
  -mindepth INT                      Minimum depth of target SNP (10)
  -maxdepth INT                      Maximum depth of target SNP (99)
  -sim_replication INT               Simulation replication counts of SNP-index (10000)
  -individual INT                    Number of bulked individuals (20)
  -generation STR                    Bulked generation "F2", "RIL" or "BC1F1" (F2)
  -window_size INT                   Window size for sliding window analysis (2000000)
  -step_size INT                     Step size for sliding window analysis (50000)
  -min_plot INT                      The minimum number of SNPs in the window (10)
  -filt_bias FLOAT                   If both bulks show set value<(SNP-index) and set value>1-(SNP-index), the position is exclueded. (0.3)
 ```

### command : mqtlseq_plot
  
  Description: This is the commands that can be started from txt file calculated SNP-index and filtered with P1, P2 and F1.

  Required options
  ```
  -txt TXT                           Full path of txt file calculated SNP-index and filtered with P1, P2 and F1
  ```
  
  Additional options
  ```
  -p2 STR                            Specify the name if using P2 (False)
  -a STR                             A-bulk name (A_bulk)
  -b STR                             Specify the name if using B-bulk (False)
  -f1 STR                            Select True if you are using F1 (False)
  -o STR                             Output directry name (QTL-seq)
  -t INT                             Number of parallel processing (1)
  -mindepth INT                      Minimum depth of target SNP (10)
  -maxdepth INT                      Maximum depth of target SNP (99)
  -sim_replication INT               Simulation replication counts of SNP-index (10000)
  -individual INT                    Number of bulked individuals (20)
  -generation STR                    Bulked generation "F2", "RIL" or "BC1F1" (F2)
  -window_size INT                   Window size for sliding window analysis (2000000)
  -step_size INT                     Step size for sliding window analysis (50000)
  -min_plot INT                      The minimum number of SNPs in the window (10)
  -filt_bias FLOAT                   If both bulks show set value<(SNP-index) and set value>1-(SNP-index), the position is exclueded. (0.3)
 ```

## The input file format
### -rt TXT
Please enter paired-end reads separated by tabs.
'''
P1 reads
/Full/Path/P1.1.fq	/Full/Path/P1.2.fq
P2 reads
/Full/Path/P2.1.fq	/Full/Path/P2.2.fq
A bulk reads
/Full/Path/A_bulk.1.fq	/Full/Path/A_bulk.2.fq
B bulk reads
/Full/Path/B_bulk.1.fq	/Full/Path/B_bulk.2.fq
F1 reads
/Full/Path/F1.1.fq	/Full/Path/F1.2.fq
'''

If you do not use P2, B_bulk, or F1, please do not write fq. In that case, please do not delete P2 reads, B bulk reads, and F1 reads, but fill them in.
'''
P1 reads
/Full/Path/P1.1.fq	/Full/Path/P1.2.fq
P2 reads
A bulk reads
/Full/Path/A_bulk.1.fq	/Full/Path/A_bulk.2.fq
B bulk reads
F1 reads
'''

If the read is divided into multiple FQs, it can be written on a new line.
'''
P1 reads
/Full/Path/P1_1.1.fq	/Full/Path/P1_1.2.fq
/Full/Path/P1_2.1.fq	/Full/Path/P1_2.2.fq
/Full/Path/P1_3.1.fq	/Full/Path/P1_3.2.fq
P2 reads
A bulk reads
/Full/Path/A_bulk.1.fq	/Full/Path/A_bulk.2.fq
B bulk reads
F1 reads
'''

### -bt TXT
Please enter bam files.
'''
P1 bam
/Full/Path/P1.bam
P2 bam
/Full/Path/P2.bam
A bulk bam
/Full/Path/A_bulk.bam
B bulk bam
/Full/Path/B_bulk.bam
F1 bam
/Full/Path/F1.bam
'''

If you do not use P2, B_bulk, or F1, please do not write bam. In that case, please do not delete P2 bam, B bulk bam, and F1 bam, but fill them in.
'''
P1 bam
/Full/Path/P1.bam
P2 bam
A bulk bam
/Full/Path/A_bulk.bam
B bulk bam
F1 bam
'''

## The example of format execution
example1: mqtlseq
'''
mqtlseq -r /Full/Path/refarence.fa \
        -rt /Full/Path/reads_PATH.txt\
        -p1 High_cultivar \
        -p2 Low_cultivar \
        -a High_bulk\
        -b Low_nulk\
        -o High_vs_Low\
        -window_size 5000\
        -step_size 1000\
        -min_plot 1
'''

example2: mqtlseq_mpileup
'''
mqtlseq_mpileup -r /Full/Path/refarence.fa \
                -bt /Full/Path/bams_PATH.txt\
                -p1 High_cultivar \
                -p2 Low_cultivar \
                -a High_bulk\
                -b Low_nulk\
                -o High_vs_Low\
                -window_size 5000\
                -step_size 1000\
                -min_plot 1
'''

example3: mqtlseq_snpindex
'''
mqtlseq_snpindex -vcf /Full/Path/SNP_called.vcf \
                 -p1 High_cultivar \
                 -p2 Low_cultivar \
                 -a High_bulk\
                 -b Low_nulk\
                 -f1 True\
                 -o High_vs_Low\
                 -window_size 5000\
                 -step_size 1000\
                 -min_plot 1
'''

example4: mqtlseq_plot
'''
mqtlseq_plot -txt /Full/Path/SNP_called.vcf \
             -p1 High_cultivar \
             -p2 True \
             -a High_bulk\
             -b Low_nulk\
             -f1 True\
             -o High_vs_Low\
             -window_size 5000\
             -step_size 1000\
             -min_plot 1
'''

## The output file format
### output directory
'''
QTL-seq---log
        |-10_ref
        |-11_bam
        |-12_vcf
        |-13_snpindex
        |-14_filt
        |-20_ref
        |-21_bam
        |-22_vcf
        |-23_snpindex
        |-24_filt
        |-30_sim
        |-40_plot
'''

### 10_ref
The index files of reference are output.
'''
ref.fa              linked reference file
ref.fa.bwt          index files
ref.fa.amb          index files
ref.fa.ann          index files
ref.fa.pac          index files
ref.fa.sa           index files
ref.fa.fai          index files
'''

### 11_bam
The aligned results of P1 are output.
'''
mem.txt             command file 1
view.txt            command file 2
merge.txt           command file 3
sort.txt            command file 4
rmdup.txt           command file 5
index.txt           command file 6
coverage.txt        command file 7
rm.txt              command file 8
P1.bam              P1 bam file
P1.bam.bai          P1 bam index file
P1.coverage.txt     P1 coverage file
'''

### 12_vcf
The result of calling SNP with bcftools is output.
'''
mpileup.txt         command file 1
call.txt            command file 2
filter.txt          command file 3
rm.txt              command file 4
P1.vcf              vcf file that called SNP of P1
'''

### 13_snpindex
The result of calculating SNP-index is output.
'''
P1.snp_index.txt    Text file containing SNP-index of P1
'''
Chromosome name, position, reference base, SNP base, depth, and SNP-index are listed in tab-separated format.

### 14_filt
The result of selecting the SNP to be substituted with the reference will be output.
'''
P1_homo.snp_index.txt    A file in which the depth>=mindepth and SNP-index=1 is selected from P1.snp_index.txt
'''
The file format is the same as P1.snp_index.txt.

### 20_ref
The index files of reference are output.
'''
P1.fa              Fasta file with P1 SNPs substituted
P1.fa.bwt          index files
P1.fa.amb          index files
P1.fa.ann          index files
P1.fa.pac          index files
P1.fa.sa           index files
P1.fa.fai          index files
'''

### 21_bam
The aligned results of P1 are output.
'''
mem.txt                 command file 1
view.txt                command file 2
merge.txt               command file 3
sort.txt                command file 4
rmdup.txt               command file 5
index.txt               command file 6
coverage.txt            command file 7
rm.txt                  command file 8
P1.bam                  P1 bam file
P1.bam.bai              P1 bam index file
P1.coverage.txt         P1 coverage file
P2.bam                  P2 bam file
P2.bam.bai              P2 bam index file
P2.coverage.txt         P2 coverage file
F1.bam                  F1 bam file
F1.bam.bai              F1 bam index file
F1.coverage.txt         F1 coverage file
A_bulk.bam              A_bulk bam file
A_bulk.bam.bai          A_bulk bam index file
A_bulk.coverage.txt     A_bulk coverage file
B_bulk.bam              B_bulk bam file
B_bulk.bam.bai          B_bulk bam index file
B_bulk.coverage.txt     B_bulk coverage file
'''

### 22_vcf
The result of calling SNP with bcftools is output.
'''
mpileup.txt                   command file 1
call.txt                      command file 2
filter.txt                    command file 3
rm.txt                        command file 4
P1.P2.F1.A_bulk.B_bulk.vcf    vcf file that called SNP of P1, P2, F1, A_bulk and B_bulk
'''
P1, P2, F1, A_bulk and B_bulk are called in this order, and samples not in the input are omitted.

### 23_snpindex
The result of calculating SNP-index is output.
'''
P1.P2.F1.A_bulk.B_bulk.snp_index.txt    Text file containing SNP-index of P1, P2, F1, A_bulk and B_bulk
'''
Chromosome name, position, reference base, SNP base, depth, and SNP-index are listed in tab-separated format.
P1, P2, F1, A_bulk and B_bulk are called in this order, and samples not in the input are omitted.

### 24_filt
The selection results of the 3 filter steps are output.

Step1: P1 homo
'''
P1_homo.P2.F1.A_bulk.B_bulk.snp_index.txt    A file in which the depth>=mindepth and SNP-index=0 in P1 is selected
'''
The file format is the same as P1.P2.F1.A_bulk.B_bulk.snp_index.txt

Step2: P2 homo
'''
P1_homo.P2_homo.F1.A_bulk.B_bulk.snp_index.txt    A file in which the depth>=mindepth and SNP-index=1 in P2 is selected
'''
The file format is the same as P1.P2.F1.A_bulk.B_bulk.snp_index.txt

Step3: F1 hetero
'''
P1_homo.P2_homo.F1_hetero.A_bulk.B_bulk.snp_index.txt    The range of SNP-index value on heterozygous SNPs
'''
The file format is the same as P1.P2.F1.A_bulk.B_bulk.snp_index.txt

'''
F1_sim.txt        Range of SNP-index values for heterozygous SNPs(95% confidence interval) used in F1 hetero filter
'''
Depth, minimum value, and maximum value are listed separated by tabs.

### 30_sim
The calculation result of the confidence interval in bulk is output.
'''
F2_20_sim.txt     Range of SNP-index and delta SNP-index values for bulk sample
'''
Depth, 99% under interval of SNP-index, 99% top interval of SNP-index, 95% under interval of SNP-index, 95% top interval of SNP-index, 99% under interval of delta SNP-index, 99% Top interval of delta SNP-index, 95% under interval of delta SNP-index, 95% top interval of delta SNP-index are listed separated by tabs.

### 40_plot
Depth filter and graph results in bulk are output.
'''
A_bulk.snp_index.sim.txt    A file containing the SNP-index of the position showing A_bulk depth>=mindepth.
B_bulk.snp_index.sim.txt    A file containing the SNP-index of the position showing B_bulk depth>=mindepth.
'''
Chromosome, position, reference base, SNP base, depth of A_bulk, SNP-index of A_bulk, under and top of 99% confidence interval of SNP-index, under and top of 95% confidence interval of SNP-index, delta SNP The under and top of the 99% confidence interval for delta SNP-index and the under and top of the 95% confidence interval for delta SNP-index are tab-separated.

'''
A_bulk-B_bulk.snp_index.sim.txt    A file containing the delta SNP-index of the position.
'''
Chromosome, position, reference base, SNP base, higher depth of A_bulk and B_bulk, delta SNP-index(A_bulk-B_bulk), under and top of 99% confidence interval of SNP-index, under and top of 95% confidence interval of SNP-index, delta SNP The under and top of the 99% confidence interval for delta SNP-index and the under and top of the 95% confidence interval for delta SNP-index are tab-separated.

'''
A_bulk.2000000bp_window.50000bp_step.txt    A file containing the sliding window average of SNP-index in A-bulk.
B_bulk.2000000bp_window.50000bp_step.txt    A file containing the sliding window average of SNP-index in B-bulk.
'''
Chromosome, position, A_bulk SNP-index average values, under and top of 99% confidence interval of SNP-index, under and top of 95% confidence interval of SNP-index, delta SNP The under and top of the 99% confidence interval for delta SNP-index and the under and top of the 95% confidence interval for delta SNP-index are output in tab-separated format.

'''
A_bulk-B_bulk.2000000bp_window.50000bp_step.txt    A file containing the sliding window average of delta SNP-index.
'''
Chromosome, position, delta SNP-index average values, under and top of 99% confidence interval of SNP-index, under and top of 95% confidence interval of SNP-index, delta SNP The under and top of the 99% confidence interval for delta SNP-index and the under and top of the 95% confidence interval for delta SNP-index are output in tab-separated format.

'''
A_bulk.png
B_bulk.png
A_bulk-B_bulk.png
'''
Dots indicate SNP-index or delta SNP-index at each chromosomal position, and the red line shows the sliding window average value.
Orange and green lines indicate 99% and 95% confidence intervals, respectively.
