## Calculate genetic position of SNPs based on a genetic map
## Kaichi Huang 2020 Jun

use warnings;
use strict;

my $snp = $ARGV[0];
my $map = $ARGV[1];

my %map;
open MAP, $map;
my $first = 1;
while (<MAP>){
	if ($first) {
		$first = 0;
	}else{
		chomp;
		my @a = split(/\s/, $_);
		$map{$a[0]}={"bp" => [], "cm" => []} if !exists $map{$a[0]};
		push @{$map{$a[0]}{"bp"}}, $a[1];
		push @{$map{$a[0]}{"cm"}}, $a[2];
	}
}
close MAP;

print "position\trate\n";
open SNP, $snp;
while (<SNP>){
	chomp;
	my @b = split(/\t/,$_);
	my $chr = $b[0];
	my $bp = $b[1];
	my ($bp1,$bp2,$cm1,$cm2);
	for (my $i = 0; $i < $#{$map{$chr}{"bp"}}; $i = $i + 1){
		if ($map{$chr}{"bp"}[$i] <= $bp && $map{$chr}{"bp"}[$i+1] > $bp){
			$bp1=$map{$chr}{"bp"}[$i];
			$bp2=$map{$chr}{"bp"}[$i+1];
			$cm1=$map{$chr}{"cm"}[$i];
			$cm2=$map{$chr}{"cm"}[$i+1];
		}
	}
	my $cm = $cm1+(($cm2-$cm1)*(($bp-$bp1)/($bp2-$bp1)));
	print "$bp\t$cm\n";
}
close SNP;
