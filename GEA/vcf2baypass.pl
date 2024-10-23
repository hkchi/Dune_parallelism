## This script transforms SNP vcf to BayPass genotype input format
## Kaichi Huang 2021 Jun
## modified from Kaichi Huang 2019 Jul

use warnings;
use strict;
use List::MoreUtils qw(uniq);

my $out_prefix = $ARGV[1];
my %pop; # pop info
my @fields; # vcf columns, including individual names

# Population information
open POP, $ARGV[0];
while (<POP>){
	chomp;
	my @a = split(/\t/, $_);
	my $sample_name="NegPG_".$a[0]; # sample name conversion
	$pop{$sample_name} = $a[1];
}
close POP;
my @all_pop = uniq values %pop;
if ($ARGV[2] eq "ex") {
	@all_pop = grep {$_ ne "MON001" && $_ ne "MON300" && $_ ne "MON650"} @all_pop; # Optional argument to remove populations which have missing total.N data
}
open OUT_POP, ">$out_prefix.pop";
foreach my $pp (@all_pop) {
	print OUT_POP "$pp\n";
}
close OUT_POP;

# Traversal
open OUT, ">$out_prefix.txt";
open OUT_SNP, ">$out_prefix.snp";
while (<STDIN>) {
	if(/^##/){next;}
	chomp;
	if(/^#/){
		@fields = split(/\t/, $_);
		next;
	}
	my @a = split(/\t/, $_);
	my $chr = $a[0];
	my $pos = $a[1];
	my %pop_gt;
	foreach my $pp (@all_pop) {
		$pop_gt{$pp} = {"c0" => 0, "c1" => 0};
	}
	foreach my $i (9..$#a) {
		my @tmp = split(/:/,$a[$i]);
		my $gt = $tmp[0];
		if ($gt eq "0/0") {
			$pop_gt{$pop{$fields[$i]}}{"c0"} += 2;
		} elsif ($gt eq "1/1") {
			$pop_gt{$pop{$fields[$i]}}{"c1"} += 2;
		} elsif ($gt eq "0/1") {
			$pop_gt{$pop{$fields[$i]}}{"c0"} += 1;
			$pop_gt{$pop{$fields[$i]}}{"c1"} += 1;
		}
	}
	my @line;
	foreach my $pp (@all_pop) {
		push @line, "$pop_gt{$pp}{'c0'} $pop_gt{$pp}{'c1'}";
	}
	print OUT join(" ",@line),"\n";
	print OUT_SNP "$chr\t$pos\n";
}
close OUT;
close OUT_SNP;
