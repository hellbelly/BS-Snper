use strict;
die "perl $0 snpfile cgfile >cgfilterfile\n" unless(@ARGV==2);
my $snp=shift;
my $cpg=shift;

my %hash;

open SNP, $snp or die $!;
open CG,$cpg or die $!;

while(<SNP>){
	chomp;
	my @a=split;
	if($a[3] eq "C" && $a[4] eq "T"){
		$hash{$a[0]}{$a[1]}=1;
	}	
	if($a[3] eq "G" && $a[4] eq "A"){
		my $pos=$a[1]-1;
		$hash{$a[0]}{$pos}=1;

	}

}
close SNP;

while(<CG>){
	chomp;
	my @a=split;
	unless(exists($hash{$a[0]}{$a[1]})){

		print $_."\n";
	}
	
}
close CG;
