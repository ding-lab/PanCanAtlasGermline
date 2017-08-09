#!/usr/bin/env perl
#
# (rjm), July 2017
# Syntax:   uncompressed vcf |  $thisScript


use strict;
use warnings;

my @myList=();
my @a;

while(<STDIN>) {
    chomp;
    if( /^#/ ) {
	# get list of input filenames given to merge
	if( /vcf-merge/ ) {
	    @a = split /\s+/;
	    for(my $i = 1; $i < scalar @a; $i++) {
		my %data = ('inputfile' => $a[ $i ], 'samplename' => "");
		push @myList, \%data;
	    }
	}
	if( /^#CHROM/ ) {
	    @a = split /\t/;
	    for(my $i = 9 ; $i < scalar @a; $i++) {
		$myList[ $i - 9 ]{'samplename'} = $a[ $i ];
	    }
	}
	
    } else {
	last;
    }
}

for(my $j=0 ; $j < scalar @myList; $j++) {
    print $j,"\t", $myList[$j]{'inputfile'},"\t", $myList[$j]{'samplename'}, "\n";
}
