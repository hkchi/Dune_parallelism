## Classify CLR outlier windows and merge outlier windows into regions
## Kaichi Huang 2021 Jan

use warnings;
use strict;
use List::Util qw/sum max min/;
use POSIX;

my $input = $ARGV[0];
my $out = $ARGV[1];

my $window_size = 20000;
my $step = 2000;
my $ids = 5;
my $num = 1;

my $f_switch = 1; # switch to output outlier classification for all windows

# Read in CLR results
open CLR, "<$input" or die $!;
my @scores = ();
my %rod;
while (<CLR>) {
	chomp;
	next if /LR/;
	my ($chr, $pos, $rod) = (split)[0, 1, 2];
	$rod{$chr}{$pos} = $rod;
	push @scores, $rod;
}
close CLR;
my @chr_array = keys %rod;

# Outlier classification
my @sorted = sort {$b <=> $a} @scores;
my $quantile = int(0.05 * @scores);
my $threshold = $sorted[$quantile];

my %outlier;
foreach my $chr (sort {$a <=> $b} @chr_array) {
	foreach my $pos (sort {$a <=> $b} keys %{$rod{$chr}}) {
		my $rod = $rod{$chr}{$pos};
		if ($rod > $threshold) {
			$outlier{$chr}{$pos} = $rod;
		}
	}
}

# Merge outlier windows into regions
my %merged;
my $buffer = $ids * $step;
foreach my $chr (sort {$a <=> $b} @chr_array) {
	my $group = 1;
	my @windows = sort {$a <=> $b} keys %{$outlier{$chr}};
	next if @windows == 0;
	my $query = shift @windows;
	push @{$merged{$chr}{$group}}, $query;
	my $query_buffer = $query + $buffer;
	foreach my $window (@windows) {
		if ($window <= $query_buffer) {
			push @{$merged{$chr}{$group}}, $window;
		}
		else {
			$group++;
			push @{$merged{$chr}{$group}}, $window;
		}
		$query = $window;
		$query_buffer = $query + $buffer;
	}
}

# Output the results
open REGION, ">$out.CLR.selected_region";
print REGION "#top5%: $threshold\n";
print REGION "#chr\tgroup\tstart\tend\tsize\tmax_score\taverage_score\n";
my %region_outlier;
foreach my $chr (sort {$a <=> $b} @chr_array) {
	foreach my $group (sort {$a <=> $b} keys %{$merged{$chr}}) {
		my @windows = sort {$a <=> $b} @{$merged{$chr}{$group}};
		next if @windows < $num;
		my $start = POSIX::floor($windows[0])-$window_size/2;
		my $end = POSIX::floor($windows[-1])+$window_size/2;
		my $size = int($end - $start + 1);
		my @scores;
		foreach my $window (@windows) {
			push @scores, $rod{$chr}{$window};
			$region_outlier{$chr}{$window} = 1;
		}
		my $max_score = max(@scores);
		my $average_score = sum(@scores)/@scores;
		print REGION "$chr\t$group\t$start\t$end\t$size\t$max_score\t$average_score\n" ;
	}
}
close REGION;

if ($f_switch) {
	open OUT, ">$out.CLR_outlier.txt";
	print OUT "chr\tlocation\tLR\toutlier\n";
	foreach my $chr (sort {$a <=> $b} @chr_array) {
		foreach my $pos (sort {$a <=> $b} keys %{$rod{$chr}}) {
			if (exists $region_outlier{$chr}{$pos}) {
				print OUT "$chr\t$pos\t$rod{$chr}{$pos}\tOutlier\n";
			} else {
				print OUT "$chr\t$pos\t$rod{$chr}{$pos}\tNon-outlier\n";
			}
		}
	}
	close OUT;
}
