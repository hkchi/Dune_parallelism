## This script prepares environmental dataset for BayPass input
## Kaichi Huang 2021 Jun
## modified from Kaichi Huang 2019 Jul

use warnings;
use strict;

my $out_prefix = $ARGV[1];
my @all_pop;
my @all_cov;
my %all_info; # values of all covariates of all populations

# Read in population order from genetic data
open POP, $ARGV[0];
while (<POP>){
	chomp;
	push @all_pop, $_;
}
close POP;

# Store all info
my $first = 1;
while (<STDIN>) {
	chomp;
	my @a = split(/\t/, $_);
	if ($first) {
		foreach my $i (1..$#a) {
			push @all_cov, $a[$i];
			$all_info{$a[$i]} = {};
		}
		$first = 0;
		next;
	}
	my $pop = $a[0];
	foreach my $i (1..$#a) {
		$all_info{$all_cov[$i-1]}{$pop} = $a[$i];
	}
}

# Output
open OUT, ">$out_prefix.txt";
open OUT_COV, ">$out_prefix.cov";
foreach my $cov (@all_cov) {
	my @values = values $all_info{$cov};
	# print only suitable covariates
	if (($#all_pop == 14 && (grep /NA/, @values)) || ($#all_pop == 17 && !(grep /NA/, @values))) {
		print OUT_COV "$cov\n";
		my @line;
		foreach my $pp (@all_pop) {
			push @line, "$all_info{$cov}{$pp}";
		}
		print OUT join(" ",@line),"\n";
	}
}
close OUT;
close OUT_COV;
