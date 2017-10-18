#!/usr/bin/env perl
#
# Jay Mashl, July 2017
# Syntax:   uncompressed vcf |  $thisScript
# adopted by Kuan Oct 2017 for bcftools output and update directly to TCGA ID


use strict;
use warnings;

my @myList=();
my @a;
my $samples={};

my $fn = "/gscmnt/gc3020/dinglab/medseq/Germline/projects/PanCanAtlasGermline/TCGA_data/sampleQC/pca_table.20171017.tsv";
#my $fn = "../../../../TCGA_data/sampleQC/pca_table.20171017.tsv";
open ( IN , "<$fn" ) or die "Cannot open $fn: $!";
while ( <IN> )
{
	chomp;
	my @line = split "\t" , $_;
    $samples->{$line[0]}=$line[1];
}
close IN;

while(<STDIN>) {
    chomp;
    if( /^#/ ) {
	# get list of input filenames given to merge
	#if( /vcf-merge/ ) {
	if( /bcftools_mergeCommand/ ){
	    @a = split /\s+/;
	    #for(my $i = 1; $i < scalar @a; $i++) {
	    for(my $i = 5; $i < scalar @a; $i++) {
		my %data = ('inputfile' => $a[ $i ], 'samplename' => "");
		push @myList, \%data;
	    }
	}
	if( /^#CHROM/ ) {
	    @a = split /\t/;
	    for(my $i = 9 ; $i < scalar @a; $i++) {
#		$myList[ $i - 9 ]{'samplename'} = $a[ $i ];

                # in this application, extract unique identifier from first field
		my @b = split /\./, $myList[ $i - 9 ]{'inputfile'};
		my @c = split /\//, $b[0];
		my $TCGA = $samples->{$c[1]};
		$a[ $i ] = $TCGA;
	    }
	}

	#Print
	if( /^#CHROM/ ) {
	    print join("\t", @a),"\n";
	} else {
	    print $_,"\n";
	}
	
    } else  {
	last;
    }
}

#for(my $j=0 ; $j < scalar @myList; $j++) {
#    print $j,"\t", $myList[$j]{'inputfile'},"\t", $myList[$j]{'samplename'}, "\n";
#}
