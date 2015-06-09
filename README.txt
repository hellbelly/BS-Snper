1.	Introduction
The purpose of this program is to call SNP using bisulfite sequencing data.

2.	Install
Download BS-Snper from https://github.com/hellbelly/BS-Snper. After extracting the downloaded package, execute the command ./BS-Snper.sh. Make sure the executable file rrbsSnp is generated.

3.	Run
You can run BS-Snper in Linux or MAC OS. To run the program, use the command like:
perl BS-Snper.pl --interval hg19.len --fa hg19.fa --input merge.sort.bam --output SNP.candidate --methoutput Meth.out --minhetfreq 0.1 --minhomfreq 0.85 --minquali 15 --mincover 10 --maxcover 1000 --minread2 2 --errorate 0.02 --mapvalue 20 >SNP.out2 2>SNP.log

The meanings of the arguments are as follows.
	interval: A two-column file storing the length of reference sequence in the 2nd column (sequence name in the 1st column)
	fa: Fasta format of reference sequence
	input: Input bam file
	output: Preliminary SNP output
	methoutput: Methy output
	minhetfreq: The minimal frequency of heterozygous SNPs
	minhomfreq: The minimal frequency of homozygous SNPs
	minquali: The minimal sequencing quality of SNPs
	mincover: The minimal number of covered reads
	minread2: The minimal number of supporting reads
	errorate: The minimal frequency of SNPs
	mapvalue: The minimal threshold of mapping quality (Phred-scaled)

PS: 
All the chromosome names of the reference sequence file (eg. hg19.fa) must be specified in the interval file (eg. hg19.len) with correct length.

