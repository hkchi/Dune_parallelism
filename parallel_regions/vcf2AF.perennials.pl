## This script transforms SNP vcf to SweepFinder2 allele frequency input format (sites polarized)
## Kaichi Huang 2020 Nov

use warnings;
use strict;
use Getopt::Long;
use List::Util qw/sum/;
use List::MoreUtils qw/indexes/;

my $snp_pre = "";
GetOptions('snp=s' => \$snp_pre);
my @snp;

my @perennial_list = ("664647_GIG","DEC_1895","DIV_1956","GRO_2043");
my @perennial_index;

print "position\tx\tn\tfolded\n";
while (<STDIN>) {
	if(/^##/){next;}
	chomp;
	my @a = split(/\t/, $_);
	if (/^#/) {
		my %b = map{$_ => 1} @perennial_list;
		@perennial_index = indexes { $b{$_} } @a;
	} else {
		my $chr = $a[0];
		my $pos = $a[1];
		my (@gt, @perennial_gt);
		foreach my $i (9..$#a) {
			my @tmp = split(/:/,$a[$i]);
			my $gt = $tmp[0];
			my @push;
			if ($gt eq "0/0" || $gt eq "0|0") {
				@push = (0,0);
			} elsif ($gt eq "1/1" || $gt eq "1|1") {
				@push = (1,1);
			} elsif ($gt eq "0/1" || $gt eq "0|1") {
				@push = (0,1);
			}
			if (grep /^$i$/, @perennial_index) {
				push @perennial_gt, @push;
			} else {
				push @gt, @push;
			}
		}
		my $tally = @gt;
		# output
		if (@gt) {
			if ( @perennial_gt<8 ) {
				# not polarized (folded=1) if not all perennials called
				if ( sum(@gt)>0 && sum(@gt)<$tally ) {
					# polymorphic only
					my @out_array = ($pos, sum(@gt), $tally, 1);
					print join("\t",@out_array),"\n";
					#
					if ($snp_pre) {
						push @snp, "$chr\t$pos";
					}
				}
			} else {
				# all perennials called
				if ( (sum(@gt)==0 && sum(@perennial_gt)==8) || (sum(@gt)==$tally && sum(@perennial_gt)==0) ) {
					# a substitution
					my @out_array = ($pos, $tally, $tally, 0);
					print join("\t",@out_array),"\n";
					#
					if ($snp_pre) {
						push @snp, "$chr\t$pos";
					}
				} elsif ( sum(@gt)>0 && sum(@gt)<$tally ) {
					# polymorphic
					if (sum(@perennial_gt)==0) {
						my @out_array = ($pos, sum(@gt), $tally, 0);
						print join("\t",@out_array),"\n";
					} elsif (sum(@perennial_gt)==8) {
						my @out_array = ($pos, $tally-sum(@gt), $tally, 0);
						print join("\t",@out_array),"\n";
					} else {
						my @out_array = ($pos, sum(@gt), $tally, 1);
						print join("\t",@out_array),"\n";
					}
					#
					if ($snp_pre) {
						push @snp, "$chr\t$pos";
					}
				}
			}
		}
	}
}

if ($snp_pre) {
	open OUT, ">$snp_pre.snps";
	foreach my $i (@snp) {
		print OUT $i,"\n";
	}
	close OUT;
}
