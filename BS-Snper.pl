#!/usr/bin/perl -w

=head1 Name

        BS-Snper.pl  -- Detect the SNP according to the candidate positions


=head1 Version

        Author: Shengjie Gao, gaoshengjie@genomics.cn; Dan Zou, zoudan.nudt@gmail.com; Xuesong HU, huxs001@gmail.com.

	Version: 1.1,  Date: 2018-07-27 

=head1 Usage

        --variants-only         output variant sites only
        --regions-file <file>   restrict to regions listed in a file
        --help                  output help information to screen

=head1 Exmple
	perl BS-Snper.pl --fa hg19.fa --input BSMAP.sort.bam --output snp.candidate.out --methcg meth.cg --methchg meth.chg --methchh meth.chh --minhetfreq 0.1 --minhomfreq 0.85 --minquali 15 --mincover 10 --maxcover 1000 --minread2 2 --errorate 0.02 --mapvalue 20 >SNP.out 2>SNP.log
	#perl BS-Snper.pl --fa hg19.fa --input sort.bam --output result --methoutput meth.out --minhetfreq 0.1 --minhomfreq 0.85   --minquali 15 --mincover 10 --maxcover 1000 --minread2 2 --errorate 0.02 >SNP.out 2>SNP.log\n
=cut


use strict;
use Getopt::Long;
use POSIX;
use FindBin '$Bin';
use File::Path;  ## function " mkpath" and "rmtree" deal with directory
use File::Basename qw(basename dirname);
my ($Help,$fasta,$bam,$mapvalue,$minhetfreq,$minhomfreq,$minquali,$minread2,$mincover,$maxcover,$errorate,$pvalue,$interval,$output,$methcg,$methchg,$methchh);
my ($regFile,$varOnly);
GetOptions(
    "fa:s"=>\$fasta,
	"input:s"=>\$bam,
	"regions-file:s"=>\$regFile,
	"variants-only"=>\$varOnly,
	"output:s"=>\$output,
	"methcg:s"=>\$methcg,
	"methchg:s"=>\$methchg,
	"methchh:s"=>\$methchh,
	"minhetfreq:i"=>\$minhetfreq,
	"minhomfreq:i"=>\$minhomfreq,
	"minquali:i"=>\$minquali,
	"mincover:i"=>\$mincover,
	"maxcover:i"=>\$maxcover,
	"minread2:i"=>\$minread2,
	"errorate:i"=>\$errorate,
	"mapvalue:i"=>\$mapvalue,
    "help"=>\$Help
);
die `pod2text $0` if (@ARGV==0 || $Help);
$minhetfreq ||=0.1;
$minhomfreq ||=0.85;
$minquali ||=15;
$mincover ||=10;
$minread2 ||=2;
$maxcover ||=1000;
$errorate ||=0.02;
$mapvalue ||=20;
#$pvalue ||=0.01;

my $getFlagf = sub {return 1};
if (defined $varOnly) {
	$getFlagf = sub {return 0};
}
my %theRegion;
if (defined $regFile) {
	my ($cntL,$cntP) = (0,0);
	open R,'<',$regFile or die $!;
	while (<R>) {
		chomp;
		my ($chr,$begin,$end)=split /\s+/;
		next unless defined $end;
		++$cntL;
		for my $i ($begin .. $end) {
			$theRegion{$chr}{$i}=1;
			++$cntP;
		}
	}
	close R;
	warn "[!]Regions Load: $cntL lines of $cntP points.\n";
	$getFlagf = sub($$) {
		if (exists $theRegion{$_[0]} and exists $theRegion{$_[0]}{$_[1]}) {
			return 1;
		} else {
			return 0;
		}
	}
}
### Example ###
	my $chr='1'; my $pos=22345;
	my $flag = &{$getFlagf}($chr,$pos);
	warn "FLAG: $flag\n";
######

my $eee=2.7;
#$interval = $fasta . ".len";
#if(!(-e $interval)) {
#	system("$Bin/chrLenExtract $fasta");
#}

if(system("$Bin/rrbsSnp $fasta $bam $output $methcg $methchg $methchh $minquali $mincover $maxcover $minhetfreq $errorate $mapvalue") != 0) {
        die "Error!";
}


#if(system("$Bin/rrbsSnp $interval $fasta $bam $output $methcg $methchg $methchh $minquali $mincover $maxcover $minhetfreq $errorate $mapvalue") != 0) {
#	die "Error!";
#}
#
#
my $year_month_day=strftime("%Y%m%d",localtime());


print "##fileformat=VCFv4.3
##fileDate= $year_month_day
##bssnperVersion=1.1
##bssnperCommand=--fa $fasta   --input $bam --output $output --methcg $methcg --methchg $methchg --methchh $methchh --minhetfreq  $minhetfreq --minhomfreq  $minhomfreq --minquali $minquali --mincover $mincover --maxcover $maxcover --minread2 $minread2 --errorate $errorate --mapvalue $mapvalue
##reference=file://$fasta
##Bisulfite=directional>
";
open INTV,"$Bin/samtools-0.1.19/samtools view -H $bam|" or die $!;
#@SQ     SN:chr3 LN:198295559
while(<INTV>){
	chomp;
	my @a=split;
	my ($chr,$length);
	if(/SN\:(\S+)/){
		$chr=$1;		
	}
	if(/LN\:(\d+)/){
		$length=$1;
	}
 	print "##contig=<ID=$chr,length=$length>\n";
}
close INTV;

print "##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">
##INFO=<ID=IDV,Number=1,Type=Integer,Description=\"Maximum number of reads supporting an indel\">
##INFO=<ID=IMF,Number=1,Type=Float,Description=\"Maximum fraction of reads supporting an indel\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">
##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Total allelic depths\">
##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Total allelic depths on the forward strand\">
##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Total allelic depths on the reverse strand\">
##INFO=<ID=MQ0F,Number=1,Type=Float,Description=\"Fraction of MQ0 reads (smaller is better)\">
##FILTER=<ID=Low,Description=\"Low Quality\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes for each ALT allele, in the same order as listed\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases\">
##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Average mapping quality\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of high-quality bases\">
##FORMAT=<ID=SP,Number=1,Type=Integer,Description=\"Phred-scaled strand bias P-value\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">
##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths on the forward strand\">
##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths on the reverse strand\">
##FORMAT=<ID=BSD,Number=8,Type=Integer,Description=\"Depth, ATCG in watson strand and crick strand\">
##FORMAT=<ID=BSQ,Number=8,Type=Integer,Description=\"Avarage Base Quality, ATCG in watson strand and crick strand\">
##FORMAT=<ID=ALFR,Number=R,Type=Float,Description=\"Allele frequency\">
";

#print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGENOTYPE\tFREQUENCY\tNumber_of_watson[A,T,C,G]\tNumber_of_crick[A,T,C,G]\tMean_Quality_of_Watson[A,T,C,G]\tMean_Quality_of_Crick[A,T,C,G]\n";
#chr1    10583   G       70,0,0,24       0,2,0,243       33,0,0,32       0,35,0,33
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$bam\n";
open SNP,$output or die "no snp file\n";

while(<SNP>){
 next unless(/\w/);
 next if(/^#/);
 my @a=split;
 my $geno=&genotype(\@a);

}

#print "SNP finished\n";
#chr1    10583   G       70,0,0,24       0,2,0,243       33,0,0,32       0,35,0,33
#$minfreq,$minquali,$minread2,$mincover,$pvalue,$interval
#chr1    10006   .       C       T,A,G   1000    PASS    DP=2;ADF=0,1;ADR=0,1;AD=0,2;    GT:DP:ADF:ADR:AD:BSD:BSQ        0/1:6:0,0,0,0:54,10,1,1:54,10,1,1:0,0,0,0,1,54,10,1:0,0,0,0,37,36,38,22 0/2:6:2,0,10,0:0,0,0,0:2,0,10,0:0,0,0,0,1,54,10,1:0,0,0,0,37,36,38,22
#chr1    10105   .       A       C       27      Low     DP=2;ADF=0,1;ADR=0,1;AD=0,2;    GT:DP:ADF:ADR:AD:BSD:BSQ        0/1:6:2,0:2,2:4,2:0,0,0,0,69,9,0,0:0,0,0,0,36,32,0,0    0/1:6:2,0:2,2:4,2:0,0,0,0,1,54,10,1:0,0,0,0,37,36,38,22

#chr1    10006   .       C       T       1000    PASS    CT      0.152   0,0,0,0 1,10,54,1       0,0,0,0 37,38,36,22 

sub genotype
{
	my $line=shift;
	my @lines=@{$line};
	my $genoreturn=&Bayes($line);
	my @bayes=split /\t/,$genoreturn;
	my $genoqual=$bayes[1];
	my $genotypemaybe=$bayes[0];
	my @watson=split /\,/,$lines[3];
	my @crick=split /\,/,$lines[4];
	my @wsq=split /\,/,$lines[5];
	my @crq=split /\,/,$lines[6];
	#print "$lines[0]\t$lines[1]\t$lines[2]\t$genotypemaybe\n";
	my $totaldepth=$watson[0]+$watson[1]+$watson[2]+$watson[3]+$crick[0]+$crick[1]+$crick[2]+$crick[3];
	
	if($genotypemaybe eq "AA"){#genotypeis AA
		if($lines[2] =~/A/i){
			return "REF";	
		}
		if($lines[2] =~/T/i){#T>AA
			my $qvalue=($wsq[0]>$crq[0])?$wsq[0]:$crq[0];
			my $depth=$watson[0]+$crick[0]+$watson[1]+$crick[1];
			my $var=$watson[0]+$crick[0];
			my $adf="$watson[1]\,$watson[0]";
			my $adr="$crick[1]\,$crick[0]";
			my $ad0=$watson[1]+$crick[1];
			my $ad1=$watson[0]+$crick[0];
			my $ad="$ad0\,$ad1";
			my $ref=sprintf("%.3f",$ad0/$totaldepth);
			my $T2A=sprintf("%.3f",$var/$totaldepth);
			if($depth >= $mincover  && $qvalue >= $minquali && $var >=$minread2 ){
				#my $T2A=sprintf("%.3f",$var/$totaldepth);
				#my $T=sprintf("%.3f",$ad0/$totaldepth);
				if($T2A>=$minhomfreq){
					#print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tPASS\tAA\t$T2A\t".join("\t",@lines[3..6])."\n";
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tPASS\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2A\n";
				}else{
					#print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tAA\t$T2A\t".join("\t",@lines[3..6])."\n";
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2A\n";
				}				
			}else{
				my $T2A=sprintf("%.3f",$var/$totaldepth);				
				#print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tAA\t$T2A\t".join("\t",@lines[3..6])."\n";
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2A\n";
			}									
		}

		if($lines[2] =~/C/i){#C>AA
			my $qvalue=($wsq[0]>$crq[0])?$wsq[0]:$crq[0];
			my $depth=$watson[2]+$crick[2]+$watson[0]+$crick[0];
            my $var=$watson[0]+$crick[0];
			my $adf="$watson[2]\,$watson[0]";
            my $adr="$crick[2]\,$crick[0]";
            my $ad0=$watson[2]+$crick[2];
            my $ad1=$watson[0]+$crick[0];
            my $ad="$ad0\,$ad1";
			my $ref=sprintf("%.3f",$ad0/$totaldepth);
			if($depth >= $mincover  && $qvalue >= $minquali && $var >=$minread2 ){
                    my $C2A=sprintf("%.3f",$var/$totaldepth);
                    if($C2A>=$minhomfreq){
                                       #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tPASS\tAA\t$C2A\t".join("\t",@lines[3..6])."\n";
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tPASS\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2A\n";		
                                }else{
                                       # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tAA\t$C2A\t".join("\t",@lines[3..6])."\n";
                                	print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2A\n";
								}   
    
              }else{
				my $C2A=sprintf("%.3f",$var/$totaldepth);
                # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tAA\t$C2A\t".join("\t",@lines[3..6])."\n";
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2A\n";
                        } 					
		}

		if($lines[2] =~/G/i){##G>AA
			my $qvalue=$wsq[0];
            my $depth=$watson[3]+$watson[0]+$crick[0]+$crick[3];
            my $var=$watson[0]+$crick[0];
			my $G2A;
			my $adf="$watson[3]\,$watson[0]";
            my $adr="$crick[3]\,$crick[0]";
            my $ad0=$watson[3]+$crick[3];
            my $ad1=$watson[0]+$crick[0];
            my $ad="$ad0\,$ad1";	
			my $ref;
			if($depth>0){
				$G2A=sprintf("%.3f",$var/$depth);
				$ref=sprintf("%.3f",$ad0/$depth);
			}else{
				$G2A=0;
			}

            if($depth >= $mincover  && $qvalue >= $minquali && $var >=$minread2 ){
                 if($G2A>=$minhomfreq){
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tPASS\tAA\t$G2A\t".join("\t",@lines[3..6])."\n";
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tPASS\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2A\n";
                 }else{
                                       # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tAA\t".join("\t",@lines[3..6])."\n";
                    print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2A\n";
                 }

             }else{
                                #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tAA\t$G2A\t".join("\t",@lines[3..6])."\n";
                                print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2A\n";                        }
		}
	}
#genotype is AT

	if($genotypemaybe eq "AT"){ #genotype is AT
		if($lines[2] =~/A/i){#A>AT
			my $qvalue=($wsq[1]>$crq[1])?$wsq[1]:$crq[1];;
			my $depth=$watson[0]+$watson[1]+$crick[0]+$crick[1];
			my $var=$watson[1]+$crick[1];
	
			my $adf="$watson[0]\,$watson[1]";
            my $adr="$crick[0]\,$crick[1]";
            my $ad0=$watson[0]+$crick[0];
            my $ad1=$watson[1]+$crick[1];
            my $ad="$ad0\,$ad1";
			my $ref=sprintf("%.3f",$ad0/$totaldepth);
			if($depth >= $mincover  && $qvalue >= $minquali && $var >=$minread2 ){
               my $A2T=sprintf("%.3f",$var/$totaldepth);
               if($A2T>=$minhetfreq){
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tPASS\tAT\t$A2T\t".join("\t",@lines[3..6])."\n";
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tPASS\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2T\n";
               }else{
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2T\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tAT\t$A2T\t".join("\t",@lines[3..6])."\n";
               }
             }else{
				my $A2T=sprintf("%.3f",$var/$totaldepth);
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2T\n";
                               # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tAT\t$A2T\t".join("\t",@lines[3..6])."\n";
			} 
          }   

          if($lines[2] =~/T/i){#T>AT			
            		my $qvalue=($wsq[0]>$crq[0])?$wsq[0]:$crq[0];
                    my $depth=$watson[0]+$watson[1]+$crick[0]+$crick[1];
                    my $var=$crick[0]+$watson[0];
                    my $adf="$watson[1]\,$watson[0]";
                    my $adr="$crick[1]\,$crick[0]";
                    my $ad0=$watson[1]+$crick[1];
                    my $ad1=$watson[0]+$crick[0];
                    my $ad="$ad0\,$ad1"; 
					my $ref=sprintf("%.3f",$ad0/$totaldepth);
                        if($depth >= $mincover  && $qvalue >= $minquali && $var >=$minread2 ){
                                my $T2A=sprintf("%.3f",$var/$totaldepth);
                                if($T2A>=$minhetfreq){
									print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tPASS\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2A\n";
                                       # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tPASS\tAT\t$T2A\t".join("\t",@lines[3..6])."\n";
                                }else{
									print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2A\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tAT\t$T2A\t".join("\t",@lines[3..6])."\n";
                                }

                        }else{
							my $T2A=sprintf("%.3f",$var/$totaldepth);
                                #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tAT\t$T2A\t".join("\t",@lines[3..6])."\n";
							print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2A\n";
                        }
                }   
          if($lines[2] =~/C/i){#C>AT  #add the watson in the new version ##294
			my $qvalueA=($wsq[0]>$crq[0])?$wsq[0]:$crq[0];
			my $qvalueT=$crq[1];
			my $varA=$watson[0]+$crick[0];
			my $varT=$crick[1]+$watson[1];
			my $depth=$watson[0]+$crick[0]+$crick[1]+$watson[1];
			my $adf="$watson[2]\,$watson[0]\,$watson[1]";
            my $adr="$crick[2]\,$crick[0]\,$crick[1]";
            my $ad0=$watson[2]+$crick[2];
            my $ad1=$watson[0]+$crick[0];
			my $ad2=$watson[1]+$crick[1];
            my $ad="$ad0\,$ad1\,$ad2";
			my $ref=sprintf("%.3f",$ad0/$totaldepth);		
							
			if($depth >= $mincover  && $qvalueA >= $minquali && $qvalueT>=$minquali && $varA >=$minread2 && $varT>=$minread2 ){
                my $C2A=sprintf("%.3f",$varA/$totaldepth);
				my $C2T=sprintf("%.3f",$varT/$totaldepth);
                if($C2A>=$minhetfreq && $C2T>=$minhetfreq){
                      print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA,T\t$genoqual\tPASS\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2A\,$C2T\n";
                }else{
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA,T\t$genoqual\tLow\tAT\t$C2A\,$C2T\t".join("\t",@lines[3..6])."\n";
                                        print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA,T\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2A\,$C2T\n";
                }   

             }else{
                my $C2A=sprintf("%.3f",$varT/$totaldepth);
                my $C2T=sprintf("%.3f",$varT/$totaldepth);
				#print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tAT\t$genoqual\tLow\tAT\t$C2A\,$C2T\t".join("\t",@lines[3..6])."\n";
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA,T\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2A\,$C2T\n";
             }  				
          }   
          if($lines[2] =~/G/i){#G>AT  #add the crick in the new version ##326
			my $qvalueA=$wsq[0];
            my $qvalueT=($wsq[1]>$crq[1])?$wsq[1]:$crq[1];
            my $varA=$watson[0]+$crick[0];
            my $varT=$crick[1]+$watson[1];
            my $depth=$watson[0]+$crick[0]+$crick[1]+$watson[1];
			my $adf="$watson[3]\,$watson[0]\,$watson[1]";
            my $adr="$crick[3]\,$crick[0]\,$crick[1]";
            my $ad0=$watson[3]+$crick[3];
            my $ad1=$watson[0]+$crick[0];
            my $ad2=$watson[1]+$crick[1];
            my $ad="$ad0\,$ad1\,$ad2";
            my $ref=sprintf("%.3f",$ad0/$totaldepth);

            if($depth >= $mincover  &&  $qvalueA >= $minquali && $qvalueT>=$minquali && $varA >=$minread2 && $varT>=$minread2 ){
                   my $G2A=sprintf("%.3f",$varA/$totaldepth);
                   my $G2T=sprintf("%.3f",$varT/$totaldepth);
                     if($G2A>=$minhetfreq && $G2T>=$minhetfreq){
					 print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA,T\t$genoqual\tPASS\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2A\,$G2T\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tAT\t$genoqual\tPASS\tAT\t$G2A\,$G2T\t".join("\t",@lines[3..6])."\n";
                   }else{
					 print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA,T\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2A\,$G2T\n";
                                       # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tAT\t$genoqual\tLow\tAT\t$G2A\,$G2T\t".join("\t",@lines[3..6])."\n";
                   }

              }else{
                my $G2A=sprintf("%.3f",$varA/$totaldepth);
                my $G2T=sprintf("%.3f",$varT/$totaldepth);
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA,T\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2A\,$G2T\n";
                                #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tAT\t$genoqual\tLow\tAT\t$G2A\,$G2T\t".join("\t",@lines[3..6])."\n";
              }
            }   				
	}
#AC	#363
	if($genotypemaybe eq "AC"){ #genotype is AC # add watson T
      if($lines[2] =~/A/i){ #A>AC
            my $qvalueC=($wsq[2]>$crq[2])?$wsq[2]:$crq[2];
            my $varC=$watson[2]+$crick[2]+$watson[1];
            my $depth=$watson[0]+$crick[2]+$watson[2]+$crick[0];
			my $watsonC=$watson[1]+$watson[2];
			my $adf="$watson[0]\,$watsonC";
            my $adr="$crick[0]\,$crick[2]";
            my $ad0=$watson[0]+$crick[0];
            my $ad1=$watson[2]+$crick[2]+$watson[1];
            my $ad="$ad0\,$ad1";
			my $ref=sprintf("%.3f",$ad0/$totaldepth);							
			my $A2C=sprintf("%.3f",$varC/$totaldepth);
               if($depth >= $mincover  && $qvalueC >= $minquali  && $varC >=$minread2 ){
                   if($A2C>=$minhetfreq){
                                       # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tPASS\tAC\t$A2C\t".join("\t",@lines[3..6])."\n";
                     print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tPASS\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2C\n";
                   }else{
					 print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2C\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tAC\t$A2C\t".join("\t",@lines[3..6])."\n";
                   }   
				
                }else{
				  print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2C\n";
                                #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tAC\t$A2C\t".join("\t",@lines[3..6])."\n";
                }   			
         }
      if($lines[2] =~/T/i){#T>AC ############modify totaldepth 390
			my $qvalueA=($wsq[0]>$crq[0])?$wsq[0]:$crq[0];
            my $qvalueC=($wsq[2]>$crq[2])?$wsq[2]:$crq[2];
            my $varA=$watson[0]+$crick[0];
            my $varC=$crick[2]+$watson[2]+$watson[1];
            my $depth=$watson[0]+$crick[2]+$watson[2]+$crick[0];
			my $T2A=sprintf("%.3f",$varA/$totaldepth);
            my $T2C=sprintf("%.3f",$varC/$totaldepth);
			my $watsonC=$watson[1]+$watson[2];
			my $adf="$watson[1]\,$watson[0]\,$watsonC";
            my $adr="$crick[1]\,$crick[0]\,$crick[2]";
            my $ad0=$crick[1];
            my $ad1=$watson[0]+$crick[0];
            my $ad2=$watsonC+$crick[2];
            my $ad="$ad0\,$ad1\,$ad2";
            my $ref=sprintf("%.3f",$ad0/$totaldepth);
	
            if($depth >= $mincover  && $qvalueA >= $minquali && $qvalueC>=$minquali && $varA >=$minread2 && $varC>=$minread2 ){
                if($T2A>=$minhetfreq && $T2C>=$minhetfreq){
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA,C\t$genoqual\tPASS\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2A\,$T2C\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tAC\t$genoqual\tPASS\tAC\t$T2A\,$T2C\t".join("\t",@lines[3..6])."\n";
                                }else{
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA,C\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2A\,$T2C\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tAC\t$genoqual\tLow\tAC\t$T2A\,$T2C\t".join("\t",@lines[3..6])."\n";
                                }   

                }else{
				  print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA,C\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2A\,$T2C\n";
                               # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tAC\t$genoqual\tLow\tAC\t$T2A\,$T2C\t".join("\t",@lines[3..6])."\n";
                        }       
		
                }

        if($lines[2] =~/C/i){#C>AC 423 add wastonT
			  my $qvalueA=($wsq[0]>$crq[0])?$wsq[0]:$crq[0];
              my $varA=$watson[0]+$crick[0];
			  my $watsonC=$watson[1]+$watson[2]+$crick[2];
              my $depth=$watson[0]+$crick[2]+$watson[2]+$crick[0];
			
              my $adf="$watson[0]\,$watsonC";
              my $adr="$crick[0]\,$crick[2]";
              my $ad0=$watson[2]+$crick[2]+$watson[1];
              my $ad1=$watson[0]+$crick[0];
              my $ad="$ad0\,$ad1";
              my $ref=sprintf("%.3f",$ad0/$totaldepth);

              my $C2A=sprintf("%.3f",$varA/$totaldepth);
              if($depth >= $mincover  && $qvalueA >= $minquali  && $varA >=$minread2 ){
                    if($C2A>=$minhetfreq){
					   print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tPASS\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2A\n";	
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tPASS\tAC\t$C2A\t".join("\t",@lines[3..6])."\n";
                    }else{
					   print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2A\n";
                                       # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tAC\t$C2A\t".join("\t",@lines[3..6])."\n";
                    }

              }else{
				 print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2A\n";
                                #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tAC\t$C2A\t".join("\t",@lines[3..6])."\n";
              }			
          } 
    
          if($lines[2] =~/G/i){#G>AC 452
			my $qvalueA=$wsq[0];
                        my $qvalueC=($wsq[2]>$crq[2])?$wsq[2]:$crq[2];
                        my $varA=$watson[0]+$crick[0];
                        my $varC=$crick[2]+$watson[2]+$watson[1];
                        my $depth=$watson[0]+$crick[2]+$watson[2]+$crick[0];
                        my $G2A=sprintf("%.3f",$varA/$totaldepth);
                        my $G2C=sprintf("%.3f",$varC/$totaldepth);
			
						my $watsonC=$watson[1]+$watson[2];
                        my $adf="$watson[3]\,$watson[0]\,$watsonC";
                        my $adr="$crick[3]\,$crick[0]\,$crick[2]";
                        my $ad0=$crick[3]+$watson[3];
                        my $ad1=$watson[0]+$crick[0];
                        my $ad2=$watsonC+$crick[2];
                        my $ad="$ad0\,$ad1\,$ad2";
                        my $ref=sprintf("%.3f",$ad0/$totaldepth);

            if($depth >= $mincover  && $qvalueA >= $minquali && $qvalueC>=$minquali && $varA >=$minread2 && $varC>=$minread2 ){
                  if($G2A>=$minhetfreq && $G2C>=$minhetfreq){
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA,C\t$genoqual\tPASS\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2A\,$G2C\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tAC\t$genoqual\tPASS\tAC\t$G2A\,$G2C\t".join("\t",@lines[3..6])."\n";
                  }else{
					 print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA,C\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2A\,$G2C\n";
                                       # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tAC\t$genoqual\tLow\tAC\t$G2A\,$G2C\t".join("\t",@lines[3..6])."\n";
                                }

              }else{
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA,C\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2A\,$G2C\n";
                               # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tAC\t$genoqual\tLow\tAC\t$G2A\,$G2C\t".join("\t",@lines[3..6])."\n";
              }
           }
        }   
	
#AG genotype	
	if($genotypemaybe eq "AG"){ #only watson strand
          if($lines[2] =~/A/i){#A>AG modify depth  489##
			my $qvalueG=($wsq[3]>$crq[3])?$wsq[3]:$crq[3];
            my $varG=$watson[3];
            my $depth=$watson[0]+$watson[3];##
            my $A2G=sprintf("%.3f",$watson[3]/($watson[3]+$watson[0]));

			my $adf="$watson[0]\,$watson[3]";
            my $adr="0\,0";
            my $ad0=$watson[0];
            my $ad1=$watson[3];
            my $ad="$ad0\,$ad1";
            my $ref=sprintf("%.3f",$ad0/($ad0+$ad1));

            if($depth >= $mincover  && $qvalueG >= $minquali  && $varG >=$minread2 ){
                if($A2G>=$minhetfreq){
					 print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tPASS\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2G\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tPASS\tAG\t$A2G\t".join("\t",@lines[3..6])."\n";
                 }else{
					 print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2G\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tAG\t$A2G\t".join("\t",@lines[3..6])."\n";
                 }   

             }else{
				 print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2G\n";
                                #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tAC\t$A2G\t".join("\t",@lines[3..6])."\n";
             }    			
		}
    if($lines[2] =~/T/i){#T>AG 515 
			my $qvalueA=$wsq[0];
            my $qvalueG=($wsq[3]>$crq[3])?$wsq[3]:$crq[3];
            my $varA=$watson[0];
            my $varG=$watson[3];
            my $depth=$watson[0]+$watson[3];
            my $T2A=sprintf("%.3f",$varA/$depth);
            my $T2G=sprintf("%.3f",$varG/$depth);
			
			my $adf="$watson[1]\,$watson[0]\,$watson[3]";
            my $adr="0\,0\,0";
            my $ad0=$watson[1]+$crick[1];
            my $ad1=$watson[0];
            my $ad2=$watson[3];
            my $ad="$ad0\,$ad1\,$ad2";
            my $ref=sprintf("%.3f",$ad0/$depth);
						
            if($depth >= $mincover  && $qvalueA >= $minquali && $qvalueG>=$minquali && $varA >=$minread2 && $varG>=$minread2 ){
                if($T2A>=$minhetfreq && $T2G>=$minhetfreq){
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA,G\t$genoqual\tPASS\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2A\,$T2G\n";	
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tAG\t$genoqual\tPASS\tAG\t$T2A\,$T2G\t".join("\t",@lines[3..6])."\n";
                 }else{
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA,G\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2A\,$T2G\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tAG\t$genoqual\tLow\tAG\t$T2A\,$T2G\t".join("\t",@lines[3..6])."\n";
                 }

             }else{
				 print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA,G\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2A\,$T2G\n";
                               # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tAG\t$genoqual\tLow\tAG\t$T2A\,$T2G\t".join("\t",@lines[3..6])."\n";
             }
	
           }
     if($lines[2] =~/C/i){#C>AG  551
			  my $qvalueA=$wsq[0];
              my $qvalueG=($wsq[3]>$crq[3])?$wsq[3]:$crq[3];
              my $varA=$watson[0];
              my $varG=$watson[3];
              my $depth=$watson[0]+$watson[3];
              my $C2A=sprintf("%.3f",$varA/$depth);
              my $C2G=sprintf("%.3f",$varG/$depth);
			
			  my $adf="$watson[2]\,$watson[0]\,$watson[3]";
              my $adr="0\,0\,0";
              my $ad0=$watson[1];
              my $ad1=$watson[0];
              my $ad2=$watson[3];
              my $ad="$ad0\,$ad1\,$ad2";
              my $ref=sprintf("%.3f",$ad0/$depth);
              if($depth >= $mincover  && $qvalueA >= $minquali && $qvalueG>=$minquali && $varA >=$minread2 && $varG>=$minread2 ){
                 if($C2A>=$minhetfreq && $C2G>=$minhetfreq){ #568
                    print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA,G\t$genoqual\tPASS\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2A\,$C2G\n";
					# print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tAG\t$genoqual\tPASS\tAG\t$C2A\,$C2G\t".join("\t",@lines[3..6])."\n";
                  }else{
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA,G\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2A\,$C2G\n";
                    #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tAG\t$genoqual\tLow\tAG\t$C2A\,$C2G\t".join("\t",@lines[3..6])."\n";
                  }
               }else{
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA,G\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2A\,$C2G\n";
                                #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tAG\t$genoqual\tLow\tAG\t$C2A\,$C2G\t".join("\t",@lines[3..6])."\n";
               }
	
            }
     if($lines[2] =~/G/i){#G>AG
			my $qvalueA=$wsq[0];
            my $varA=$watson[0];
            my $depth=$watson[3]+$watson[0];
            my $G2A=sprintf("%.3f",$varA/$depth);
			
			my $adf="$watson[3]\,$watson[0]";
            my $adr="0\,0";
            my $ad0=$watson[3];
            my $ad1=$watson[0];
            my $ad="$ad0\,$ad1";
            my $ref=sprintf("%.3f",$ad0/($ad0+$ad1));
            if($depth >= $mincover  && $qvalueA >= $minquali  && $varA >=$minread2 ){
                  if($G2A>=$minhetfreq){
					 print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tPASS\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2A\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tPASS\tAG\t$G2A\t".join("\t",@lines[3..6])."\n";
                  }else{
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2A\n";
                                       # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tAG\t$G2A\t".join("\t",@lines[3..6])."\n";
                  }

             }else{
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2A\n";
                                #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tA\t$genoqual\tLow\tAG\t$G2A\t".join("\t",@lines[3..6])."\n";
             }
	
          }
        }   
#TT
	if($genotypemaybe eq "TT"){
        if($lines[2] =~/A/i){#A>TT 612 
			my $qvalueT=$crq[1];
            my $depth=$watson[0]+$crick[0]+$crick[1]+$watson[1];
            my $varT=$crick[1]+$watson[1];
			
			my $adf="$watson[0]\,$watson[1]";
            my $adr="$crick[0]\,$crick[1]";
                        my $ad0=$watson[0]+$crick[0];
                        my $ad1=$watson[1]+$crick[1];
                        my $ad="$ad0\,$ad1";
                        my $ref=sprintf("%.3f",$ad0/($ad0+$ad1));
			
            if($depth >= $mincover  && $qvalueT >= $minquali && $varT >=$minread2 ){
                                my $A2T=sprintf("%.3f",$varT/$depth);
                if($A2T>=$minhomfreq){ #627
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tPASS\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2T\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tPASS\tTT\t$A2T\t".join("\t",@lines[3..6])."\n";
                 }else{
					 print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2T\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tTT\t$A2T\t".join("\t",@lines[3..6])."\n";
                 }

             }else{
                                my $A2T=sprintf("%.3f",$varT/$depth);
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2T\n";
                               # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tTT\t$A2T\t".join("\t",@lines[3..6])."\n";
             }
          }
    if($lines[2] =~/T/i){
			return "REF"; 
     }
    if($lines[2] =~/C/i){#C>TT  
			my $qvalueT=$crq[1];
            my $depth=$crick[1]+$crick[2]+$watson[1]+$watson[2];
            my $varT=$crick[1]+$watson[1];
			my $C2T;
			if($depth>0){
				$C2T=sprintf("%.3f",$varT/$depth);

			}else{
				$C2T=0;				
			}

			my $adf="$watson[2]\,$watson[1]";
            my $adr="$crick[2]\,$crick[1]";
            my $ad0=$crick[2]+$watson[2];
            my $ad1=$crick[1]+$watson[1];
            my $ad="$ad0\,$ad1";
            my $ref=sprintf("%.3f",$ad0/($ad0+$ad1));
            if($depth >= $mincover  && $qvalueT >= $minquali && $varT >=$minread2 ){
                   if($C2T>=$minhomfreq){#665
					  print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tPASS\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2T\n";		
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tPASS\tTT\t$C2T\t".join("\t",@lines[3..6])."\n";
                   }else{
					  print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2T\n";
                                       # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tTT\t$C2T\t".join("\t",@lines[3..6])."\n";
                   }

             }else{
				  print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2T\n";
                               # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tTT\t$C2T\t".join("\t",@lines[3..6])."\n";
             }
	
	 }

     if($lines[2] =~/G/i){#G>TT
			my $qvalueT=$crq[1];
            my $depth=$watson[3]+$crick[3]+$crick[1]+$watson[1];
            my $varT=$crick[1]+$watson[1];
			my $adf="$watson[3]\,$watson[1]";
            my $adr="$crick[3]\,$crick[1]";
            my $ad0=$crick[3]+$watson[3];
            my $ad1=$crick[1]+$watson[1];
            my $ad="$ad0\,$ad1";
            my $ref=sprintf("%.3f",$ad0/($ad0+$ad1));

            if($depth >= $mincover  && $qvalueT >= $minquali && $varT >=$minread2 ){
                 my $G2T=sprintf("%.3f",$varT/$totaldepth);
                 if($G2T>=$minhomfreq){
					 print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tPASS\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2T\n";
                                       # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tPASS\tTT\t$G2T\t".join("\t",@lines[3..6])."\n";
                 }else{
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2T\n";
                                       # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tTT\t$G2T\t".join("\t",@lines[3..6])."\n";
                 }

            }else{
                my $G2T=sprintf("%.3f",$varT/$totaldepth);
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2T\n";
                                #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tTT\t$G2T\t".join("\t",@lines[3..6])."\n";
                        }	
                }
        }  

#CC     
        if($genotypemaybe eq "CC"){
            if($lines[2] =~/A/i){#A>CC
			   my $qvalueC=($crq[2]>$wsq[2])?$crq[2]:$wsq[2];
               my $depth=$watson[0]+$crick[0]+$crick[2]+$watson[2]+$watson[1];
               my $varC=$watson[2]+$crick[2]+$watson[1];
			   my $watsonC=$watson[1]+$watson[2];
			   my $adf="$watson[0]\,$watsonC";
               my $adr="$crick[0]\,$crick[2]";
               my $ad0=$crick[0]+$watson[0];
               my $ad1=$crick[2]+$watson[2]+$watson[1];
               my $ad="$ad0\,$ad1";
               my $ref=sprintf("%.3f",$ad0/($ad0+$ad1));

			   if($depth >= $mincover  && $qvalueC >= $minquali && $varC >=$minread2 ){
                  my $A2C=sprintf("%.3f",$varC/$depth);
                  if($A2C>=$minhomfreq){
					 print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tPASS\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2C\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tPASS\tCC\t$A2C\t".join("\t",@lines[3..6])."\n";
                  }else{
					 print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2C\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tCC\t$A2C\t".join("\t",@lines[3..6])."\n";
                  }

                }else{
                  my $A2C=sprintf("%.3f",$varC/$depth);
				  print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2C\n";
                                #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tCC\t$A2C\t".join("\t",@lines[3..6])."\n";
                }	
        }   
        if($lines[2] =~/T/i){#T>CC only Crick
			    my $qvalueC=($crq[2]>$wsq[2])?$crq[2]:$wsq[2];
                my $depth=$crick[1]+$crick[2]+$watson[1]+$watson[2];
                my $varC=$crick[2]+$watson[2]+$watson[1];
			    my $watsonC=$watson[2]+$watson[1];
			    my $adf="0\,$watsonC";
                my $adr="$crick[1]\,$crick[2]";
                my $ad0=$crick[1]+0;
                my $ad1=$crick[2]+$watsonC;
                my $ad="$ad0\,$ad1";
                my $ref=sprintf("%.3f",$ad0/($ad0+$ad1));
			
                if($depth >= $mincover  && $qvalueC >= $minquali && $varC >=$minread2 ){
                    my $T2C=sprintf("%.3f",$varC/$depth);
                    if($T2C>=$minhomfreq){
						print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tPASS\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2C\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tPASS\tCC\t$T2C\t".join("\t",@lines[3..6])."\n";
                    }else{
						print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2C\n";
                                       # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tCC\t$T2C\t".join("\t",@lines[3..6])."\n";
                    }

                }else{
					my $T2C=sprintf("%.3f",$varC/$depth);
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2C\n";
                               # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tCC\t$T2C\t".join("\t",@lines[3..6])."\n";
                }

        }   
        if($lines[2] =~/C/i){#C>CC  
				return "REF";
        }   
        if($lines[2] =~/G/i){#G>CC ######767
				my $qvalueC=($crq[2]>$wsq[2])?$crq[2]:$wsq[2];
                my $depth=$watson[3]+$crick[3]+$crick[2]+$watson[2]+$watson[1];
                my $varC=$watson[2]+$crick[2]+$watson[1];
				my $watsonC=$watson[2]+$watson[1];
				my $adf="$watson[3]\,$watsonC";
                my $adr="$crick[3]\,$crick[2]";
                my $ad0=$crick[3]+$watson[3];
                my $ad1=$crick[2]+$watsonC;
                my $ad="$ad0\,$ad1";
                my $ref=sprintf("%.3f",$ad0/($ad0+$ad1));
				
                if($depth >= $mincover  && $qvalueC >= $minquali && $varC >=$minread2 ){
					my $G2C=sprintf("%.3f",$varC/$depth);
                    if($G2C>=$minhomfreq){
						print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tPASS\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2C\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tPASS\tCC\t$G2C\t".join("\t",@lines[3..6])."\n";
                    }else{
						print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2C\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tCC\t$G2C\t".join("\t",@lines[3..6])."\n";
                    }

                 }else{
					my $G2C=sprintf("%.3f",$varC/$depth);
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2C\n";
                                #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tCC\t$G2C\t".join("\t",@lines[3..6])."\n";
                 }
             }   
        }   
#GG   804 hom add crick strand
        if($genotypemaybe eq "GG"){
			if($lines[2] =~/A/i){#A>GG hom add crick
				my $qvalueG=($crq[3]>$wsq[3])?$crq[3]:$wsq[3];
                my $depth=$watson[0]+$watson[3]+$crick[0]+$crick[3];#808
                my $varG=$watson[3]+$crick[3]+$crick[0];			
				my $adf="$watson[0]\,$watson[3]";
				my $crickG=$crick[3]+$crick[0];
                my $adr="0\,$crickG";
                my $ad0=$watson[0]+0;
                my $ad1=$watson[3]+$crickG;
                my $ad="$ad0\,$ad1";
                my $ref=sprintf("%.3f",$ad0/($ad0+$ad1));

                if($depth >= $mincover  && $qvalueG >= $minquali && $varG >=$minread2 ){
					my $A2G=sprintf("%.3f",$varG/$depth);
                    if($A2G>=$minhomfreq){
						print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tPASS\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2G\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tPASS\tGG\t$A2G\t".join("\t",@lines[3..6])."\n";
                    }else{
						print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2G\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tGG\t$A2G\t".join("\t",@lines[3..6])."\n";
                    }

				}else{
					my $A2G=sprintf("%.3f",$varG/$depth);
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2G\n";
                                #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tGG\t$A2G\t".join("\t",@lines[3..6])."\n";
                 }		
             }   
             if($lines[2] =~/T/i){#T>GG 	826	
				my $qvalueG=($crq[3]>$wsq[3])?$crq[3]:$wsq[3];
                my $depth=$watson[1]+$watson[3]+$crick[0]+$crick[3]+$crick[1];
                my $varG=$watson[3]+$crick[3]+$crick[0];			
				my $adf="$watson[1]\,$watson[3]";
				my $crickG=$crick[3]+$crick[0];
                my $adr="$crick[1]\,$crickG";
                my $ad0=$watson[1]+$crick[1];
                my $ad1=$watson[3]+$crickG;
                my $ad="$ad0\,$ad1";
                my $ref=sprintf("%.3f",$ad0/($ad0+$ad1));
			
                if($depth >= $mincover  && $qvalueG >= $minquali && $varG >=$minread2 ){
                    my $T2G=sprintf("%.3f",$varG/$depth);
                    if($T2G>=$minhomfreq){
						print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tPASS\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2G\n";
                                       # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tPASS\tGG\t$T2G\t".join("\t",@lines[3..6])."\n";
                    }else{
						print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2G\n";
                                       # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tGG\t$T2G\t".join("\t",@lines[3..6])."\n";
                                }

                 }else{
                    my $T2G=sprintf("%.3f",$varG/$depth);
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2G\n";
                    #            print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tGG\t$T2G\t".join("\t",@lines[3..6])."\n";
                 }	
			} 
            if($lines[2] =~/C/i){#C>GG 
				my $qvalueG=($crq[3]>$wsq[3])?$crq[3]:$wsq[3];
				my $depth=$watson[2]+$watson[3]+$crick[2]+$crick[3]+$crick[1]+$watson[1];
                my $varG=$watson[3]+$crick[3]+$crick[0];			
				my $adf="$watson[2]\,$watson[3]";
				my $crickG=$crick[3]+$crick[0];
                my $adr="$crick[2]\,$crickG";
                my $ad0=$watson[2]+$crick[2];
                my $ad1=$watson[3]+$crickG;
                my $ad="$ad0\,$ad1";
                my $ref=sprintf("%.3f",$ad0/($ad0+$ad1));

                if($depth >= $mincover  && $qvalueG >= $minquali && $varG >=$minread2 ){
                    my $C2G=sprintf("%.3f",$varG/$depth);
                    if($C2G>=$minhomfreq){
						print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tPASS\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2G\n";
                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tPASS\tGG\t$C2G\t".join("\t",@lines[3..6])."\n";
                    }else{
						print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2G\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tGG\t$C2G\t".join("\t",@lines[3..6])."\n";
                    }

                }else{
                    my $C2G=sprintf("%.3f",$varG/$depth);
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2G\n";
                                #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tGG\t$C2G\t".join("\t",@lines[3..6])."\n";
                }
            }  
            if($lines[2] =~/G/i){
				return "REF";
            }   
        }
#CT	884
	if($genotypemaybe eq "CT"){ #only use crick
       if($lines[2] =~/A/i){#A>CT
		my $qvalueC=($wsq[2]>$crq[2])?$wsq[2]:$crq[2];
        my $qvalueT=$crq[1];
        my $varC=$crick[2];
        my $varT=$crick[1];
		my $depth=$varC+$varT+$crick[0]+$watson[0];		 			
		my $adf="0\,0\,0";
        my $adr="$crick[0]\,$crick[1]\,$crick[2]";
        my $ad0=$crick[0];
        my $ad1=$crick[1];
        my $ad2=$crick[2];
        my $ad="$ad0\,$ad1\,$ad2";
        my $ref=sprintf("%.3f",$ad0/$depth);

        if($depth >= $mincover  && $qvalueC >= $minquali && $qvalueT>=$minquali && $varC >=$minread2 && $varT>=$minread2 ){
             my $A2C=sprintf("%.3f",$varC/$depth);
             my $A2T=sprintf("%.3f",$varT/$depth);
             if($A2C>=$minhetfreq && $A2T>=$minhetfreq){
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT,C\t$genoqual\tPASS\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2T\,$A2C\n";

                          #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tCT\t$genoqual\tPASS\tCT\t$A2C\,$A2T\t".join("\t",@lines[3..6])."\n";
             }else{
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT,C\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2T\,$A2C\n";
                          #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tCT\t$genoqual\tLow\tCT\t$A2C\,$A2T\t".join("\t",@lines[3..6])."\n";
             }
         }else{
                my $A2C= sprintf("%.3f",$varC/$depth);
                my $A2T= sprintf("%.3f",$varT/$depth);
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT,C\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2T\,$A2C\n";
                    #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tCT\t$genoqual\tLow\tCT\t$A2C\,$A2T\t".join("\t",@lines[3..6])."\n";
         }
    }
    if($lines[2] =~/T/i){ #T>CT only crick strand ##########
		my $qvalueC=($wsq[2]>$crq[2])?$wsq[2]:$crq[2];
        my $varC=$crick[2];
                        #my $depth=$crick[3]+$watson[3]+$watson[0]+$watson[1]+$crick[1]+$watson[2]+$crick[2];
		my $depth=$crick[1]+$crick[2];
        my $T2C= sprintf("%.3f",$varC/$depth);
		my $adf="0\,0";
        my $adr="$crick[1]\,$crick[2]";
        my $ad0=$crick[1];
        my $ad1=$crick[2];
        my $ad="$ad0\,$ad1";
        my $ref=sprintf("%.3f",$ad0/($ad0+$ad1));
		
        if($depth >= $mincover  && $qvalueC >= $minquali  && $varC >=$minread2 ){
			if($T2C>=$minhetfreq){
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tPASS\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2C\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tPASS\tCT\t$T2C\t".join("\t",@lines[3..6])."\n";
            }else{
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2C\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tCT\t$T2C\t".join("\t",@lines[3..6])."\n";
			}
		}else{
			print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2C\n";
                                #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tCT\t$T2C\t".join("\t",@lines[3..6])."\n";
        }
	}
    if($lines[2] =~/C/i){ #C>CT
		my $qvalueT=$crq[1];
        my $varT=$crick[1];
                        #my $depth=$crick[3]+$watson[3]+$watson[0]+$watson[1]+$crick[1]+$watson[2]+$crick[2];
        my $depth=$crick[1]+$crick[2];
		my $adf="0\,0";
        my $adr="$crick[2]\,$crick[1]";
        my $ad0=$crick[2];
        my $ad1=$crick[1];
        my $ad="$ad0\,$ad1";
        my $ref=sprintf("%.3f",$ad0/($ad0+$ad1));
        my $C2T=sprintf("%.3f",$varT/$depth);
        if($depth >= $mincover  && $qvalueT >= $minquali  && $varT >=$minread2 ){
			if($C2T>=$minhetfreq){
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tPASS\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2T\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tPASS\tCT\t$C2T\t".join("\t",@lines[3..6])."\n";
            }else{
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2T\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tCT\t$C2T\t".join("\t",@lines[3..6])."\n";
            }
		}else{
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tDP=$depth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2T\n";
 
                                #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tCT\t$C2T\t".join("\t",@lines[3..6])."\n";
		}
	}
	if($lines[2] =~/G/i){#G>CT
		my $qvalueC=($wsq[2]>$crq[2])?$wsq[2]:$crq[2];
		my $qvalueT=$crq[1];
        my $varC=$crick[2];
        my $varT=$crick[1];
		my $depth=$varC+$varT+$crick[3]+$watson[3];		 			
		my $adf="0\,0\,0";
        my $adr="$crick[3]\,$crick[1]\,$crick[2]";
        my $ad0=$crick[3];
        my $ad1=$crick[1];
        my $ad2=$crick[2];
        my $ad="$ad0\,$ad1\,$ad2";
        my $ref=sprintf("%.3f",$ad0/$depth);
		if($depth >= $mincover  && $qvalueC >= $minquali && $qvalueT>=$minquali && $varC >=$minread2 && $varT>=$minread2 ){
			my $G2C=sprintf("%.3f",$varC/$depth);
			my $G2T=sprintf("%.3f",$varT/$depth);
				if($G2C>=$minhetfreq && $G2T>=$minhetfreq){
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT,C\t$genoqual\tPASS\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2T\,$G2C\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tCT\t$genoqual\tPASS\tCT\t$G2C\,$G2T\t".join("\t",@lines[3..6])."\n";
				}else{
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT,C\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2T\,$G2C\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tCT\t$genoqual\tLow\tCT\t$G2C\,$G2T\t".join("\t",@lines[3..6])."\n";
                }
			}else{
				my $G2C= sprintf("%.3f",$varC/$depth);
				my $G2T= sprintf("%.3f",$varT/$depth);
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT,C\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2T\,$G2C\n";			
                                #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tCT\t$genoqual\tLow\tCT\t$G2C\,$G2T\t".join("\t",@lines[3..6])."\n";
            }			
         }
       }   
#GT
	if($genotypemaybe eq "GT"){ ###############1000 het add watsonA
		if($lines[2] =~/A/i){ #A>GT
			my $qvalueG=($wsq[3]>$crq[3])?$wsq[3]:$crq[3];
			my $qvalueT=$crq[1];
			my $varG=$watson[3]+$crick[3]+$watson[0];
            my $varT=$crick[1]+$watson[1];
            my $depth=$varG+$varT;

			my $crickG=$crick[0]+$crick[3];
            my $adf="$watson[0]\,$watson[1]\,$watson[3]";
            my $adr="0\,$crick[1]\,$crickG";
            my $ad0=$watson[0];
            my $ad1=$watson[1]+$crick[1];
            my $ad2=$watson[3]+$crickG;
            my $ad="$ad0\,$ad1\,$ad2";
            my $ref=sprintf("%.3f",$ad0/$totaldepth);

            if($depth >= $mincover  && $qvalueG >= $minquali && $qvalueT>=$minquali && $varG >=$minread2 && $varT>=$minread2 ){
				my $A2G= sprintf("%.3f",$varG/$depth);
				my $A2T= sprintf("%.3f",$varT/$depth);
				if($A2G>=$minhetfreq && $A2T>=$minhetfreq){
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT,G\t$genoqual\tPASS\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2T\,$A2G\n";
					#print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tGT\t$genoqual\tPASS\tGT\t$A2G\,$A2T\t".join("\t",@lines[3..6])."\n";
                }else{
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT,G\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2T\,$A2G\n";
					#print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tGT\t$genoqual\tLow\tGT\t$A2G\,$A2T\t".join("\t",@lines[3..6])."\n";
				}
             }else{
				my $A2G= sprintf("%.3f",$varG/$depth);
				my $A2T= sprintf("%.3f",$varT/$depth);
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT,G\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2T\,$A2G\n";
				#print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tGT\t$genoqual\tLow\tGT\t$A2G\,$A2T\t".join("\t",@lines[3..6])."\n";
             }
		}
        if($lines[2] =~/T/i){#T>TG
			my $qvalueG=($wsq[3]>$crq[3])?$wsq[3]:$crq[3];
			my $varG=$watson[3]+$crick[3]+$crick[0];
            my $depth=$watson[1]+$crick[1]+$varG; #1049
            my $adf="$watson[1]\,$watson[3]";
			my $crickG=$crick[0]+$crick[3];
            my $adr="$crick[1]\,$crickG";
            my $ad0=$watson[1]+$crick[1];
            my $ad1=$watson[3]+$crickG;
            my $ad="$ad0\,$ad1";
            my $ref=sprintf("%.3f",$ad0/$totaldepth);
			if($depth >= $mincover  && $qvalueG >= $minquali  && $varG >=$minread2 ){
				my $T2G= sprintf("%.3f",$varG/$depth);
				if($T2G>=$minhetfreq ){
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tPASS\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2G\n";	
										#print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tPASS\tGT\t$T2G\t".join("\t",@lines[3..6])."\n";
				}else{
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2G\n";	
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tGT\t$T2G\t".join("\t",@lines[3..6])."\n";
				}

			}else{
				my $T2G= sprintf("%.3f",$varG/$depth);
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2G\n";	
				#print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tGT\t$T2G\t".join("\t",@lines[3..6])."\n";
			}	
		}

		if($lines[2] =~/C/i){#C>GT ###########1058
			my $qvalueG=($wsq[3]>$crq[3])?$wsq[3]:$crq[3];
			my $qvalueT=$crq[1];
            my $varG=$watson[3]+$crick[3]+$watson[0];
            my $varT=$crick[1]+$watson[1];
            my $depth=$varG+$varT;
			my $crickG=$crick[0]+$crick[3];
            my $adf="$watson[2]\,$watson[1]\,$watson[3]";
            my $adr="$crick[2]\,$crick[1]\,$crickG";
            my $ad0=$watson[2]+$crick[2];
            my $ad1=$watson[1]+$crick[1];
            my $ad2=$watson[3]+$crickG;
            my $ad="$ad0\,$ad1\,$ad2";
            my $ref=sprintf("%.3f",$ad0/$totaldepth);
			
            if($depth >= $mincover  && $qvalueG >= $minquali && $qvalueT>=$minquali && $varG >=$minread2 && $varT>=$minread2 ){
				my $C2G= sprintf("%.3f",$varG/$depth);
                my $C2T= sprintf("%.3f",$varT/$depth);
                if($C2G>=$minhetfreq && $C2T>=$minhetfreq){#1076
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT,G\t$genoqual\tPASS\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2T\,$C2G\n";	
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tGT\t$genoqual\tPASS\tGT\t$C2G\,$C2T\t".join("\t",@lines[3..6])."\n";
                }else{
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT,G\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2T\,$C2G\n";						
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tGT\t$genoqual\tLow\tGT\t$C2G\,$C2T\t".join("\t",@lines[3..6])."\n";
                }
             }else{
				my $C2G= sprintf("%.3f",$varG/$depth);
				my $C2T= sprintf("%.3f",$varT/$depth);
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT,G\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2T\,$C2G\n";	
                #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tGT\t$genoqual\tLow\tGT\t$C2G\,$C2T\t".join("\t",@lines[3..6])."\n";
             }
        }
        if($lines[2] =~/G/i){#G>GT
			my $qvalueT=$crq[1];
            my $varT=$crick[1]+$watson[1];
            my $depth=$watson[1]+$crick[1]+$crick[3]+$watson[3]+$crick[0];
			my $qvalueG=($wsq[3]>$crq[3])?$wsq[3]:$crq[3];
			my $varG=$watson[3]+$crick[3]+$crick[0];
  
            my $adf="$watson[3]\,$watson[1]";
			my $crickG=$crick[0]+$crick[3];
            my $adr="$crickG\,$crick[1]";
            my $ad0=$watson[3]+$crickG;
            my $ad1=$watson[1]+$crick[1];
            my $ad="$ad0\,$ad1";
            my $ref=sprintf("%.3f",$ad0/$totaldepth);

            if($depth >= $mincover  && $qvalueT >= $minquali  && $varT >=$minread2 ){
				my $G2T= sprintf("%.3f",$varT/$depth);
                if($G2T>=$minhetfreq ){
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tPASS\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2T\n";						
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tPASS\tGT\t$G2T\t".join("\t",@lines[3..6])."\n";
                }else{
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2T\n";						

                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tGT\t$G2T\t".join("\t",@lines[3..6])."\n";
                }
              }else{
				my $G2T= sprintf("%.3f",$varT/$depth);
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2T\n";	
                #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tT\t$genoqual\tLow\tGT\t$G2T\t".join("\t",@lines[3..6])."\n";
              }
         }
     }   
#CG
	if($genotypemaybe eq "CG"){##a little difficulty
		if($lines[2] =~/A/i){#A>CG
			my $qvalueG=($wsq[3]>$crq[3])?$wsq[3]:$crq[3];
			my $qvalueC=($wsq[1]>$crq[1])?$wsq[1]:$crq[1];			
			my $varG=$watson[3]+$crick[3]+$crick[0];
			my $varC=$crick[2]+$watson[2]+$watson[1];
			my $depth=$varG+$varC;

			my $crickG=$crick[0]+$crick[3];
			my $watsonC=$watson[1]+$watson[2];
            my $adf="$watson[0]\,$watsonC\,$watson[3]";
            my $adr="0\,$crick[2]\,$crickG";

            my $ad0=$watson[0]+0;
            my $ad1=$watsonC+$crick[2];
            my $ad2=$watson[3]+$crickG;
            my $ad="$ad0\,$ad1\,$ad2";
            my $ref=sprintf("%.3f",$ad0/$totaldepth);

			if($depth >= $mincover  && $qvalueG >= $minquali && $qvalueC>=$minquali && $varG >=$minread2 && $varC>=$minread2 ){
				my $A2G= sprintf("%.3f",$varG/$depth);
				my $A2C= sprintf("%.3f",$varC/$depth);
				if($A2G>=$minhetfreq && $A2C>=$minhetfreq){
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC,G\t$genoqual\tPASS\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2C\,$A2G\n";	
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tCG\t$genoqual\tPASS\tCG\t$A2C\,$A2G\t".join("\t",@lines[3..6])."\n";
				}else{
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC,G\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2C\,$A2G\n";	

                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tCG\t$genoqual\tLow\tCG\t$A2C\,$A2G\t".join("\t",@lines[3..6])."\n";
				}

			}else{
				my $A2G= sprintf("%.3f",$varG/$depth);
				my $A2C= sprintf("%.3f",$varC/$depth);
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC,G\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$A2C\,$A2G\n";	
				#print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tCG\t$genoqual\tLow\tCG\t$A2G\,$A2C\t".join("\t",@lines[3..6])."\n";
			}
	
		}
        if($lines[2] =~/T/i){#T>CG #
			my $qvalueG=($wsq[3]>$crq[3])?$wsq[3]:$crq[3];
            my $qvalueC=($wsq[1]>$crq[1])?$wsq[1]:$crq[1];
			my $varG=$watson[3]+$crick[3]+$crick[0];
			my $varC=$crick[2]+$watson[2]+$watson[1];
			my $depth=$varG+$varC;

			my $crickG=$crick[0]+$crick[3];
			my $watsonC=$watson[1]+$watson[2];
            my $adf="0\,$watsonC\,$watson[3]";
            my $adr="$crick[1]\,$crick[2]\,$crickG";

            my $ad0=$crick[1]+0;
            my $ad1=$watsonC+$crick[2];
            my $ad2=$watson[3]+$crickG;
            my $ad="$ad0\,$ad1\,$ad2";
            my $ref=sprintf("%.3f",$ad0/$totaldepth);
            
            if($depth >= $mincover  && $qvalueG >= $minquali && $qvalueC>=$minquali && $varG >=$minread2 && $varC>=$minread2 ){
				my $T2G= sprintf("%.3f",$varG/$depth);
				my $T2C= sprintf("%.3f",$varC/$depth);
				if($T2G>=$minhetfreq && $T2C>=$minhetfreq){
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC,G\t$genoqual\tPASS\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2C\,$T2G\n";	
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tCG\t$genoqual\tPASS\tCG\t$T2C\,$T2G\t".join("\t",@lines[3..6])."\n";
                }else{
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC,G\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2C\,$T2G\n";
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tCG\t$genoqual\tLow\tCG\t$T2C\,$T2G\t".join("\t",@lines[3..6])."\n";
                }
             }else{
				my $T2G= sprintf("%.3f",$varG/$depth);
				my $T2C= sprintf("%.3f",$varC/$depth);
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC,G\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t1/2:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$T2C\,$T2G\n";
                                #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tCG\t$genoqual\tLow\tCG\t$T2G\,$T2C\t".join("\t",@lines[3..6])."\n";
             }
        }
        if($lines[2] =~/C/i){#C>CG
			my $qvalueG=$wsq[3];
  			my $varG=$watson[3]+$crick[3]+$crick[0];
			my $varC=$crick[2]+$watson[2]+$watson[1];
			my $depth=$varG+$varC;

			my $crickG=$crick[0]+$crick[3];
			my $watsonC=$watson[1]+$watson[2];
            my $adf="$watsonC\,$watson[3]";
            my $adr="$crick[2]\,$crickG";

            my $ad0=$varC;
            my $ad1=$varG;
            my $ad="$ad0\,$ad1";
            my $C2G= sprintf("%.3f",$varG/$depth);
			my $ref=sprintf("%.3f",$ad0/$depth);
			#########################1206

            if($depth >= $mincover  && $qvalueG >= $minquali  && $varG >=$minread2 ){
				if($C2G>=$minhetfreq){
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tPASS\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2G\n";	
                                       # print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tPASS\tCG\t$C2G\t".join("\t",@lines[3..6])."\n";
                }else{
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2G\n";	
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tCG\t$C2G\t".join("\t",@lines[3..6])."\n";
                }
             }else{
				print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$C2G\n";	
                                #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tG\t$genoqual\tLow\tCG\t$C2G\t".join("\t",@lines[3..6])."\n";
             }			
        }
        if($lines[2] =~/G/i){#G>CG
			my $qvalueC=($wsq[1]>$crq[1])?$wsq[1]:$crq[1];
			my $varG=$watson[3]+$crick[3]+$crick[0];
			my $varC=$crick[2]+$watson[2]+$watson[1];
			my $depth=$varG+$varC;

			my $crickG=$crick[0]+$crick[3];
			my $watsonC=$watson[1]+$watson[2];
            my $adf="$watson[3]\,$watsonC";
            my $adr="$crickG\,$crick[2]";
            my $ad0=$varG;
            my $ad1=$varC;
            my $ad="$ad0\,$ad1";
            my $G2C= sprintf("%.3f",$varC/$depth);
			my $ref=sprintf("%.3f",$ad0/$depth);
			if($depth >= $mincover  && $qvalueC >= $minquali  && $varC >=$minread2 ){
				if($G2C>=$minhetfreq){
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tPASS\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2C\n";	

                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tPASS\tCG\t$G2C\t".join("\t",@lines[3..6])."\n";
				}else{
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2C\n";	
                                        #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tCG\t$G2C\t".join("\t",@lines[3..6])."\n";
				}

				}else{
					print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tDP=$totaldepth\;ADF\=$adf\;ADR=$adr\;AD=$ad\;\tGT:DP:ADF:ADR:AD:BSD:BSQ:ALFR\t0/1:$depth:$adf:$adr:$ad:".join("\,",@lines[3,4]).":".join("\,",@lines[5,6]).":$ref\,$G2C\n";	
                                #print "$lines[0]\t$lines[1]\t\.\t$lines[2]\tC\t$genoqual\tLow\tCG\t$G2C\t".join("\t",@lines[3..6])."\n";
			}
      }
    }
	
	unless($lines[2]=~/[ACGT]/i){
		return 0;
		#print "$lines[0]\t$lines[1]\t\.\t$lines[2]\t0\tSuper\tNN\t\.\t".join("\t",@lines[3..6])."\n";
	}
}


###Bayes
sub Bayes
{
    my $line=shift;
    my @lines=@{$line};
    #print "$line[0]\n";
    my $refbase=$lines[2];
    my @watson=split /\,/,$lines[3];
    my @crick=split /\,/,$lines[4];
    my @wsq=split /\,/,$lines[5];
    my @crq=split /\,/,$lines[6];
    my $ptransition=0.00066;
    my $ptransversion=0.00033;
    #error of base
    
    my $baseqWA=sprintf("%.3e",0.1**($wsq[0]/10));
    my $baseqCA=sprintf("%.3e",0.1**($crq[0]/10));
    my $baseqWT=sprintf("%.3e",0.1**($wsq[1]/10));
    my $baseqCT=sprintf("%.3e",0.1**($crq[1]/10));
    my $baseqWC=sprintf("%.3e",0.1**($wsq[2]/10));
    my $baseqCC=sprintf("%.3e",0.1**($crq[2]/10));
    my $baseqWG=sprintf("%.3e",0.1**($wsq[3]/10));
    my $baseqCG=sprintf("%.3e",0.1**($crq[3]/10));
    my $gntpmaybe; my $qualerr;
    my @totalproduct = sort{$a<=>$b} ($baseqWA,$baseqCA,$baseqWT, $baseqCT, $baseqWC, $baseqCC, $baseqWG, $baseqCG);
    if($totalproduct[0]==0 ){
        $gntpmaybe="NN";
        $qualerr=0;
        return "$gntpmaybe\t$qualerr";
    }
    
    
    #P(Di|g=Gj)
    my $Anumber;
    if($refbase eq "A"){
        $Anumber=POSIX::ceil($watson[0]+$crick[0]*0.8);
    }else{
        $Anumber=$watson[0];
    }
    my $Tnumber;
    if($refbase eq "T"){
        $Tnumber=POSIX::ceil($watson[1]*0.8+$crick[1]);
    }else{
        $Tnumber=$crick[1];
    }
    
    my $Cnumber=int($watson[2]+$crick[2]);
    my $Gnumber=int($watson[3]+$crick[3]);
    
    my $nn=&Factorial($Anumber,$Tnumber,$Cnumber,$Gnumber);
    my ($aa,$ac,$at,$ag,$cc,$cg,$ct,$gg,$gt,$tt);
    $aa=$ac=$at=$ag=$cc=$cg=$ct=$gg=$gt=$tt=0;
    my $baseqA=($baseqCA>=$baseqWA)?$baseqWA:$baseqCA;
    my $baseqT=($baseqCT>=$baseqWT)?$baseqWT:$baseqCT;
    my $baseqC=($baseqCC>=$baseqWC)?$baseqWC:$baseqCC;
    my $baseqG=($baseqCG>=$baseqWG)?$baseqWG:$baseqCG;
    if($wsq[0]>0 || $crq[0]>0){#######A>0
        #my $baseqA=($baseqCA>=$baseqWA)?$baseqWA:$baseqCA;
        my $a_aa=1-$baseqA;
        my $other=$baseqA/3;
        #$aa = $nn+ $watson[0]*log($a_aa) + ($crick[1]+$watson[1]+$watson[2]+$crick[2]+$watson[3]+$crick[3])*log($other) ;
        $aa=$nn + $Anumber*log($a_aa) + ($Tnumber+$Cnumber+$Gnumber)*log($other);
        if($crick[1]>0 || $watson[1]>0){#AT
            
            my $a_at=(1-($baseqA+$baseqT)/2)/2;#provid genotype is AT, the probility to find A.
            $other=($baseqA+$baseqT)/4;#the probability of other 2 types opear if genotype is AT.
            #$at = $nn+ ($watson[0]+$crick[1])*log($a_at) + ($crick[2]+$watson[3]+$watson[2]+$crick[3])*log($other);
            $at = $nn + ($Anumber+$Tnumber)*log($a_at) + ($Cnumber+$Gnumber)*log($other);
        }
        if($crick[2]>0 || $watson[2]>0){#AC
            my $a_ac= (1-$baseqA)/2;
            $other=$baseqA/3;
            #$ac = $nn+ ($watson[0]+$watson[2]+$crick[2])*log($a_ac) +  ($crick[1]+$watson[3]+$crick[3])*log($other);
            $ac = $nn + ($Anumber+$Cnumber)*log($a_ac) + ($Tnumber+$Gnumber)*log($other);
            
        }
        if($crick[3]>0 || $watson[3]>0){#AG
            my $a_ag=(1-$baseqA)/2;
            $other=$baseqA/3;
            #$ag = $nn+ ($watson[0]+$watson[3]+$crick[3])*log($a_ag) + ($crick[2]+$crick[1]+$watson[2])*log($other);
            #$ag = $nn + ($watson[0]+$watson[3])*log($a_ag) + ($watson[1]+$watson[2])*log($other); ###GA only use watson
            $ag = $nn +($Anumber+$Gnumber)*log($a_ag) +($Cnumber+$Tnumber)*log($other);
        }
    }
    
    if($crq[1]>0 || $wsq[1]>0){###filter wsq[1]>0 but crq[1]==0
        my $t_tt=1-$baseqT;
        #my $t_gt=(1-($baseqCT+$baseqWG)/2)/2;
        my $other = $baseqT/3;
        #type TT
        #$tt= $nn + $crick[1]*log($t_tt)+ ($watson[0]+$crick[2]+$watson[3]+$crick[3]+$watson[2])*log($other);
        $tt=$nn+$Tnumber*log($t_tt)+($Anumber+$Cnumber+$Gnumber)*log($other);
        if($watson[2]>0 || $crick[2]>0){
            my $t_ct=(1-$baseqT)/2;
            $other=$baseqCT/3;
            ##typeCT
            #$ct= $nn + ($watson[2]+$crick[1]+$crick[2])*log($t_ct) + ($watson[0]+$watson[3]+$crick[3])*log($other);
            #$ct= $nn+($crick[1]+$crick[2])*log($t_ct) + ($crick[0]+$crick[3])*log($other);##CT type only use Crick!!
            $ct= $nn+($Cnumber+$Tnumber)*log($t_ct) + ($Anumber+$Gnumber)*log($other);##CT type only use Crick!!
        }
        if($watson[3]>0 || $crick[3]>0){ #GT
            my $t_gt=(1-$baseqT)/2;
            $other=$baseqT/3;
            #$gt= $nn + ($crick[1]+$watson[3]+$crick[3]) * log($t_gt) + ($watson[0]+$crick[2]+$watson[2])*log($other);
            $gt= $nn+($Tnumber+$Gnumber)*log($t_gt) + ($Anumber+$Cnumber)*log($other);
        }
    }
    if($crq[2]>0 || $wsq[2]>0){#CC CG
        #my $baseqC=($baseqCC>=$baseqWC)?$baseqWC:$baseqCC;
        my $c_cc=1-$baseqC;
        my $other = $baseqC/3;
        ###type CC
        #$cc = $nn + ($crick[2]+$watson[2])*log($c_cc) + ($watson[0]+$crick[0]+$crick[1]+$watson[3]+$crick[3])*log($other);
        $cc = $nn + $Cnumber*log($c_cc)+($Anumber+$Tnumber+$Gnumber)*log($other); ###CC only use crick
        if($watson[3]>0 || $crick[3]>0){#CG
            my $baseqG= ($baseqWG>=$baseqCG)?$baseqCG:$baseqWG;
            my $c_cg=(1-$baseqC)/2;
            $other=($baseqC)/2;
            #$cg = $nn + ($crick[2]+$watson[3]+$watson[2]+$crick[3])*log($c_cg) + ($watson[0]+$crick[1])*log($other);
            $cg=$nn + ($Cnumber+$Gnumber)*log($c_cg) + ($Anumber+$Tnumber)*log($other);
        }
    }
    if($wsq[3]>0 || $crq[3]>0 ){
        #my $baseqG = ($baseqWG>=$baseqCG)?$baseqCG:$baseqWG;
        my $g_gg=1-$baseqG;
        my $other = $baseqG/3;
        #$gg = $nn + ($watson[3]+$crick[3])*log($g_gg) + ($watson[0]+$crick[1]+$watson[1]+$crick[2]+$watson[2])*log($other);
        $gg=$nn+ $Gnumber*log($g_gg)+($Anumber+$Cnumber+$Tnumber)*log($other); ##GG only use Watson
    }
    
    #print "$aa\t$tt\t$cc\t$gg\t$at\t$ac\t$ag\t$ct\t$cg\t$gt\n";
    
    #P(D|g=Gj)
    #sum(P(Gj)*P(D|g=Gj))
    my $fenmu=0;
    my %hash=map{($_,eval('$'."$_"))}('aa','tt','cc','gg','at','ac','ag','ct','gt','cg');
    foreach my $type(keys %hash){
        if($hash{$type}==0){
            delete($hash{$type});
        }
    }
    
    if($refbase eq "A"){
        $aa+=log(0.985);
        $tt+=log(0.000083);
        $cc+=log(0.000083);
        $gg+=log(0.00033);
        $at+=log(0.00017);
        $ac+=log(0.00017);
        $ag+=log(0.000667);
        $ct+=(log(2.78) - 8*log(10));
        $gt+=(log(1.1) - 7*log(10));
        $cg+=(log(1.1) - 7*log(10));
    }
    if($refbase eq "T"){
        $aa+=log(0.000083);
        $tt+=log(0.985);
        $cc+=log(0.00033);
        $gg+=log(0.000083);
        $at+=log(0.00017);
        $ac+=(log(1.1) - 7*log(10));
        $ag+=(log(2.78) - 8*log(10));
        $ct+=log(0.000667);
        $gt+=log(0.00017);
        $cg+=(log(1.1) - 7*log(10));
    }
    if($refbase eq "C"){
        $aa+=log(0.000083);
        $tt+=log(0.00033);
        $cc+=log(0.985);
        $gg+=log(0.000083);
        $at+=(log(1.1) - 7*log(10));
        $ac+=log(0.00017);
        $ag+=(log(2.78) - 8*log(10));
        $ct+=log(0.000667);
        $gt+=(log(1.1) - 7*log(10));
        $cg+=log(0.00017);
    }
    if($refbase eq "G"){
        $aa+=log(0.00033);
        $tt+=log(0.000083);
        $cc+=log(0.000083);
        $gg+=log(0.9985);
        $at+=(log(1.1) - 7*log(10));
        $ac+=(log(1.1) - 7*log(10));
        $ag+=(log(6.67) - 4*log(10));
        $ct+=(log(2.78) - 8*log(10));
        $gt+=(log(1.67) - 4*log(10));
        $cg+=(log(1.67) - 4*log(10));
    }
    
    
    #print "$aa\t$tt\t$cc\t$gg\t$at\t$ac\t$ag\t$ct\t$cg\t$gt\n";
    my %hash2=map{($_,eval('$'."$_"))} (keys %hash);
    foreach my $type(keys %hash2){
        #$fenmu+=2.7**$hash2{$type};
        $fenmu+=2.7**$hash2{$type};
    }
    
    my @sort = sort {$hash2{$b}<=>$hash2{$a}} keys %hash2;
    my $genotypemaybe;my $qual;
    my $prob=0;
    if(@sort==0){
        $genotypemaybe="NN";
        $qual=0;
    }else{
        $genotypemaybe=uc($sort[0]);
        #my $first=2.7**$hash2{$sort[0]};
        my $first=2.7**$hash2{$sort[0]};
        if(@sort>1){
            if($fenmu==0){
                $qual=0;
            }else{
                $prob=1-$first/$fenmu;
                if($prob==0){
                    $qual=1000;
                }else{
                    $qual=-10*log($prob)/log(10);
                }
            }
        }elsif(@sort==1){
            #print STDERR $sort[0]."\n";
            if($sort[0] eq "aa"){
                my $hom=$watson[0]+$crick[0];
                $prob=1-1/(1+0.5**$hom);
            }
            elsif($sort[0] eq "tt"){
                my $hom=$crick[1]+$watson[1];
                $prob=1-1/(1+0.5**$hom);
            }
            elsif($sort[0] eq "cc"){
                my $hom=$watson[2]+$crick[2]+$watson[1];
                $prob=1-1/(1+0.5**$hom);
            }
            elsif($sort[0] eq 'gg'){
                my $hom=$watson[3]+$crick[3]+$crick[0];
                $prob=1-1/(1+0.5**$hom);
            }else{
                $prob=1;
            }
            
            if($prob==0){
                $qual=1000;
            }else{
                $qual=-10*log($prob)/log(10);
            }
        }
    }
    
    
    
    $qual=int($qual);
    #school CT GA
    if($genotypemaybe eq "CT"){
        if($crick[1]==0){
            $genotypemaybe = "CC";
        }elsif($crick[1]<3){
            $genotypemaybe = "NN";
        }
    }
    if($genotypemaybe eq "AG"){
        if($watson[0] == 0){
            $genotypemaybe = "GG";
        }elsif($watson[0] <3){
            $genotypemaybe = "NN";
        }
    }
    if($genotypemaybe eq "TT"){
        if($crick[2]>1){
            my $cfrq=$crick[2]/($crick[2]+$crick[1]);
            if($cfrq>0.02){
                $genotypemaybe = "NN";
            }
        }
    }
    if($genotypemaybe eq "AA"){
        if($watson[3]>1){
            my $gfrq=$watson[3]/($watson[0]+$watson[3]);
            if($gfrq>0.02){
                $genotypemaybe = "NN";
            }
        }
    }
    
    
    
    return "$genotypemaybe\t$qual";
    #return $genotypemaybe;
    #print "$sort[0]\t$hash{$sort[0]}\t$sort[1]\t$hash{$sort[1]}\n";
}



###Bayes


sub Bayesold
{
	my $line=shift;
    my @lines=@{$line};
        #print "$line[0]\n";
	my $refbase=$lines[2];
    my @watson=split /\,/,$lines[3];
    my @crick=split /\,/,$lines[4];
    my @wsq=split /\,/,$lines[5];
    my @crq=split /\,/,$lines[6];
	my $ptransition=0.00066;
	my $ptransversion=0.00033;
	#error of base
    
	my $baseqWA=sprintf("%.3e",0.1**($wsq[0]/10));
	my $baseqCA=sprintf("%.3e",0.1**($crq[0]/10));
	my $baseqWT=sprintf("%.3e",0.1**($wsq[1]/10));
	my $baseqCT=sprintf("%.3e",0.1**($crq[1]/10));
	my $baseqWC=sprintf("%.3e",0.1**($wsq[2]/10));
    my $baseqCC=sprintf("%.3e",0.1**($crq[2]/10));
	my $baseqWG=sprintf("%.3e",0.1**($wsq[3]/10));
    my $baseqCG=sprintf("%.3e",0.1**($crq[3]/10));
	my $gntpmaybe; my $qualerr;	
	my @totalproduct = sort{$a<=>$b} ($baseqWA,$baseqCA,$baseqWT, $baseqCT, $baseqWC, $baseqCC, $baseqWG, $baseqCG);
	if($totalproduct[0]==0 ){
		$gntpmaybe="NN";
		$qualerr=0;	
		return "$gntpmaybe\t$qualerr"; 
	}

	#P(Di|g=Gj)
	my $nn=&Factorial($watson[0],$crick[1],$crick[2],$watson[3]);
	my ($aa,$ac,$at,$ag,$cc,$cg,$ct,$gg,$gt,$tt);
	$aa=$ac=$at=$ag=$cc=$cg=$ct=$gg=$gt=$tt=0;
	if($wsq[0]>0 ){#######A>0
		my $a_aa=1-$baseqWA;
		my $other=$baseqWA/3;	
		my $nn=&Factorial($watson[0],$crick[1],$crick[2],$watson[3]);
		$aa = $nn+ $watson[0]*log($a_aa) + ($crick[1]+$watson[1]+$watson[2]+$crick[2]+$watson[3]+$crick[3])*log($other) ;
		if($crick[1]>0){#AT
			my $a_at=(1-($baseqWA+$baseqCT)/2)/2;#provid genotype is AT, the probility to find A.
			$other=($baseqWA+$baseqCT)/4;#the probability of other 2 types opear if genotype is AT. 
			$at = $nn+ ($watson[0]+$crick[1])*log($a_at) + ($crick[2]+$watson[3]+$watson[2]+$crick[3])*log($other);
		}
		if($crick[2]>0 || $watson[2]>0){#AC
			my $a_ac= (1-$baseqWA)/2;
			$other=$baseqWA/3;
			$ac = $nn+ ($watson[0]+$watson[2]+$crick[2])*log($a_ac) +  ($crick[1]+$watson[3]+$crick[3])*log($other);
		}  
		if($crick[3]>0 || $watson[3]>0){#AG
			my $a_ag=(1-$baseqWA)/2;
			$other=$baseqWA/3;
			$ag = $nn+ ($watson[0]+$watson[3]+$crick[3])*log($a_ag) + ($crick[2]+$crick[1]+$watson[2])*log($other);
		}
	}

	if($crq[1]>0 ){###filter wsq[1]>0 but crq[1]==0
		my $t_tt=1-$baseqCT;
		#my $t_gt=(1-($baseqCT+$baseqWG)/2)/2;
		my $other = $baseqCT/3;
		#type TT
		$tt= $nn + $crick[1]*log($t_tt)+ ($watson[0]+$crick[2]+$watson[3]+$crick[3]+$watson[2])*log($other);
		if($watson[2]>0 || $crick[2]>0){
			my $t_ct=(1-$baseqCT)/2;
			$other=$baseqCT/3; 
			##typeCT
			$ct= $nn + ($watson[2]+$crick[1]+$crick[2])*log($t_ct) + ($watson[0]+$watson[3]+$crick[3])*log($other);
		}
		if($watson[3]>0 || $crick[3]>0){
			my $t_gt=(1-$baseqCT)/2;
			$other=$baseqCT/3;
			$gt= $nn + ($crick[1]+$watson[3]+$crick[3]) * log($t_gt) + ($watson[0]+$crick[2]+$watson[2])*log($other);	 
		}
	}	
	if($crq[2]>0 || $wsq[2]>0){#CC CG
		my $baseqC=($baseqCC>=$baseqWC)?$baseqWC:$baseqCC;		
		my $c_cc=1-$baseqC;
		my $other = $baseqC/3;
		###type CC
		$cc = $nn + ($crick[2]+$watson[2])*log($c_cc) + ($watson[0]+$crick[0]+$crick[1]+$watson[3]+$crick[3])*log($other);
		if($watson[3]>0 || $crick[3]>0){#CG
			my $baseqG= ($baseqWG>=$baseqCG)?$baseqCG:$baseqWG;
			my $c_cg=(1-$baseqC)/2;
			$other=($baseqC)/2;
			$cg = $nn + ($crick[2]+$watson[3]+$watson[2]+$crick[3])*log($c_cg) + ($watson[0]+$crick[1])*log($other);
		}
	}
	if($wsq[3]>0 || $crq[3]>0 ){
		my $baseqG = ($baseqWG>=$baseqCG)?$baseqCG:$baseqWG;
		my $g_gg=1-$baseqG;
		my $other = $baseqG/3;			
		$gg = $nn + ($watson[3]+$crick[3])*log($g_gg) + ($watson[0]+$crick[1]+$watson[1]+$crick[2]+$watson[2])*log($other);
	}
	
	

	#P(D|g=Gj)					
	#sum(P(Gj)*P(D|g=Gj))
	my $fenmu=0;
	my %hash=map{($_,eval('$'."$_"))}('aa','tt','cc','gg','at','ac','ag','ct','gt','cg');
	foreach my $type(keys %hash){
		if($hash{$type}==0){
			delete($hash{$type});
		}
	}

	if($refbase eq "A"){
                $aa+=log(0.985);
                $tt+=log(0.000083);
                $cc+=log(0.000083);
                $gg+=log(0.00033);
                $at+=log(0.00017);
                $ac+=log(0.00017);
                $ag+=log(0.000667);
                $ct+=(log(2.78) - 8*log(10));
                $gt+=(log(1.1) - 7*log(10));
                $cg+=(log(1.1) - 7*log(10));    
        }    
    if($refbase eq "T"){
                $aa+=log(0.000083);
                $tt+=log(0.985);
                $cc+=log(0.00033);
                $gg+=log(0.000083);
                $at+=log(0.00017);
                $ac+=(log(1.1) - 7*log(10));
                $ag+=(log(2.78) - 8*log(10));
                $ct+=log(0.000667);
                $gt+=log(0.00017);
                $cg+=(log(1.1) - 7*log(10));
    }    
    if($refbase eq "C"){
                $aa+=log(0.000083);
                $tt+=log(0.00033);
                $cc+=log(0.985);
                $gg+=log(0.000083);
                $at+=(log(1.1) - 7*log(10));
                $ac+=log(0.00017);
                $ag+=(log(2.78) - 8*log(10));
                $ct+=log(0.000667);
                $gt+=(log(1.1) - 7*log(10));
                $cg+=log(0.00017);
    }    	
	if($refbase eq "G"){
                $aa+=log(0.00033);
                $tt+=log(0.000083);
                $cc+=log(0.000083);
                $gg+=log(0.9985);
                $at+=(log(1.1) - 7*log(10));
                $ac+=(log(1.1) - 7*log(10));
                $ag+=(log(6.67) - 4*log(10));
                $ct+=(log(2.78) - 8*log(10));
                $gt+=(log(1.67) - 4*log(10));
                $cg+=(log(1.67) - 4*log(10));
        }

	
    
	my %hash2=map{($_,eval('$'."$_"))} (keys %hash);	
	foreach my $type(keys %hash2){
               $fenmu+=2.7**$hash2{$type};
        }

    	my @sort = sort {$hash2{$b}<=>$hash2{$a}} keys %hash2;
	
	my $genotypemaybe;my $qual;
	my $prob=0;
	if(@sort==0){
		$genotypemaybe="NN";
		$qual=0;
	}else{
		$genotypemaybe=uc($sort[0]);
		my $first=2.7**$hash2{$sort[0]};
	#die "$first\n$fenmu\n";
		if(@sort>1){
			if($fenmu==0){
				$qual=1000;
			}else{
				$prob=1-$first/$fenmu;
				if($prob==0){

					$qual=1000;
				}else{
					$qual=-10*log($prob)/log(10);
				}
			}
		}elsif(@sort==1){
			#print STDERR $sort[0]."\n";
			if($sort[0] eq "aa"){
				my $hom=$watson[0];	
				$prob=1-1/(1+0.5**$hom);
			}
			elsif($sort[0] eq "tt"){
                                my $hom=$crick[1];     
                                $prob=1-1/(1+0.5**$hom);
                        }  	
			elsif($sort[0] eq "cc"){
				my $hom=$watson[2]+$crick[2];
				$prob=1-1/(1+0.5**$hom);
			}
			elsif($sort[0] eq 'gg'){
				my $hom=$watson[3]+$crick[3];
                                $prob=1-1/(1+0.5**$hom);
			}else{
				$prob=1;
			}

			if($prob==0){
				$qual=1000;	
			}else{
				$qual=-10*log($prob)/log(10);
			}
		}
	}
	
	
	
	$qual=int($qual);

        return "$genotypemaybe\t$qual";
	#return $genotypemaybe;
	#print "$sort[0]\t$hash{$sort[0]}\t$sort[1]\t$hash{$sort[1]}\n";
}

sub Factorial
{
	my $aa=shift;
	my $tt=shift;
	my $cc=shift;
	my $gg=shift;
	my $total=$aa+$tt+$cc+$gg;
	my ($naa,$ntt,$ncc,$ngg,$ntotal);
	if($aa<=1){
		$naa=0;	
	}else{
		foreach my $xx(1..$aa){
			$naa+=log($xx);
		}
	}
	if($tt<=1){
                $ntt=0; 
        }else{
                foreach my $xx(1..$tt){
                        $ntt+=log($xx);
                }   
        }   
	if($cc<=1){
                $ncc=0; 
        }else{
                foreach my $xx(1..$cc){
                        $ncc+=log($xx);
                }   
        }   
	if($gg<=1){
                $ngg=0; 
        }else{
                foreach my $xx(1..$gg){
                        $ngg+=log($xx);
                }   
        }   	
	if($total<=1){
                $ntotal=0; 
        }else{
                foreach my $xx(1..$total){
                        $ntotal+=log($xx);
                }   
        }   

	my $nn=$ntotal-$naa-$ncc-$ntt-$ngg;
	return $nn;	
}



