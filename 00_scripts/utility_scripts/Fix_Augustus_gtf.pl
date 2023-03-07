#!/usr/bin/perl -w
#
#
# usage: cat INFILE.gtf ./Fix_Augustus_gtf.pl > OUTFILE.gtf
# https://github.com/Gaius-Augustus/BRAKER/issues/457
#
#
# what it does:
# Fixes lack of "gene_id" and "transcript_id" tags in column 9 of gtf file from BRAKER/Augustus
#
# CIKeeling 2022-02-24
#
#
use strict;

my @cols;
my $new_nine;
while (<>) {
	chomp;
	@cols = split /\t/, $_;
	$new_nine=$cols[8];
	if ($cols[2]eq"gene") {
		$new_nine="gene_id \"$cols[8]\"";
	} elsif ($cols[2]eq"transcript") {
		$new_nine="transcript_id \"$cols[8]\"";
	}
	print join "\t", $cols[0],$cols[1],$cols[2],$cols[3],$cols[4],$cols[5],$cols[6],$cols[7],$new_nine;
	print "\n";
}