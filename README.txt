1.	Introduction
BS-SNPer is an ultrafast and memory-efficient package, a program for BS-Seq variation detection from alignments in standard BAM/SAM format using approximate Bayesian modeling.

2.	System requirement
BS-SNPer works on Unix (Linux, Ubuntu, Mac OS, etc) based systems. 

Hardware requirements
	One computing node equipped with at least 10 GB Memory
Software requirements
	GCC 4.6.0 or higher
	Perl 5.16.3 or higher
	zlib 1.2.8 or higher
	
3.	Getting started

Installing
Download BS-SNPer from https://github.com/hellbelly/BS-Snper by clicking the button “Download ZIP”. Run the commands below.
	1.	unzip BS-Snper-master.zip
	2.	cd BS-Snper-master
	3.	sh BS-Snper.sh
Make sure the executable file rrbsSnp is generated.

Usage
You can run BS-SNPer in Linux or MAC OS, using the command like:
perl BS-Snper.pl --fa <reference_file> --input <sorted_bam_file> --output <snp_result_file> --methcg <meth_cg_result_file> --methchg <meth_chg_result_file> --methchh <meth_chh_result_file> --minhetfreq 0.1 --minhomfreq 0.85 --minquali 15 --mincover 10 --maxcover 1000 --minread2 2 --errorate 0.02 --mapvalue 20 >SNP.out 2>ERR.log

Attention
Both of the input and output file arguments should be passed to BS-SNPer in the form of absolute paths. 

Options
	--fa: Reference genome file in fasta format
	--input: Input bam file
	--output: Temporary file storing SNP candidates
	--methcg: CpG methylation information
	--methchg: CHG methylation information
	--methchh: CHH methylation information
	--minhetfreq: Threshold of frequency for calling heterozygous SNP
	--minhomfreq: Threshold of frequency for calling homozygous SNP
	--minquali: Threshold of base quality
	--mincover: Threshold of minimum depth of covered reads
	--maxcover: Threshold of maximum depth of covered reads
	--minread2: Minimum mutation reads number
	--errorate: Minimum mutation rate
	--mapvalue: Minimum read mapping value
	SNP.out: Final SNP result file
	ERR.log: Log file

4.	Input file
Any alignments in standard sorted BAM/SAM format (see https://samtools.github.io/hts-specs/SAMv1.pdf for detailed information).
A bam file for evaluation is available at ftp://public.genomics.org.cn/BGI/BS-SNPer/example/

5.	Output files
The output files include an SNP output file and a methylation output file. 
The SNP output file has a standard VCF format.
The methylation output file has a tab-separated format same as MethylExtract (http://bioinfo2.ugr.es/MethylExtract/downloads/ManualMethylExtract.pdf):
	1. CHROM: Chromosome.
	2. POS: Sequence context most 5’ position on the Watson strand (1-based).
	3. CONTEXT: Sequence contexts with the SNVs annotated using the IUPAC nucleotide ambiguity code (referred to the Watson strand).
	4. Watson METH: The number of methyl-cytosines (referred to the Watson strand).
	5. Watson COVERAGE: The number of reads covering the cytosine in this sequence context (referred to the Watson strand).
	6. Watson QUAL: Average PHRED score for the reads covering the cytosine (referred to the Watson strand).
	7. Crick METH: The number of methyl-cytosines (referred to the Watson strand).
	8. Crick COVERAGE: The number of reads covering the guanine in this context (referred to the Watson strand).
	9. Crick QUAL: Average PHRED score for the reads covering the guanine (referred to the Watson strand).

6.	Contact information
If you have any problem please do not hesitate to contact:
gaoshengjie@genomics.org.cn
zoudan_001@163.com

If you use BS-Snper, please cite:
Gao Shengjie, Zou Dan, Mao Likai, et al. BS-SNPer: SNP calling in bisulfite-seq data[J]. Bioinformatics, 2015, 31(24): 4006-4008.



