#!/usr/bin/perl
#03 August 2017 - Adam D Scott - 
# adopted by Kuan @ WashU 11/1/2017 for PCA file

use strict;
use warnings;

use IO::File;
use FileHandle;

my $usage = 'perl prepFileForAllelicImbalance.pl <charger> <preAI>
';

die $usage , unless @ARGV == 2;
my ( $charger , $preAI ) = @ARGV;

# my $IN1 = FileHandle->new( "gunzip -c $charger |" , "r" );
# if ( not defined $IN1 ) { die "ADSERROR: Could not open/read $charger\n"; }
open(my $IN1, "gunzip -c $charger |") || die "can't open pipe to $charger";

my $OUTmissense = FileHandle->new( "ready.missense.".$preAI , "w" );
if ( not defined $OUTmissense ) { die "ADSERROR: Could not open/write $preAI\n"; }

my $OUTtruncation = FileHandle->new( "ready.truncation.".$preAI , "w" );
if ( not defined $OUTtruncation ) { die "ADSERROR: Could not open/write $preAI\n"; }

my $useTheseClasses = { 'missense' => 1 , 'frame_shift_del' => 1 , 'frame_shift_ins' => 1 , 'splice_site_del' => 1 , 'splice_site_ins' => 1 , 'nonsense' => 1 , 'splice_site' => 1 , 'nonstop' => 1 };
my $sites = {};
$OUTmissense->print( "HUGO_Symbol\tCHR\tSTART\tSTOP\tREF\tALT\tTYPE\tNormalRef\tNormalVar\tNormalVAF\tTumorRef\tTumorVar\tTumorVAF\tSample\n" );
$OUTtruncation->print( "HUGO_Symbol\tCHR\tSTART\tSTOP\tREF\tALT\tTYPE\tNormalRef\tNormalVar\tNormalVAF\tTumorRef\tTumorVar\tTumorVAF\tSample\n" );

while ( my $line = <$IN1> ) {
	next if ( $line =~ /HGVSg/ );
	chomp( $line );
	my @line = split( "\t" , $line );
	next if ( $line[9] =~ /synonymous_variant/ );
	next if ( $line[9] =~ /stop_retained_variant/ );
	next if ( $line[9] =~ /start_lost/ );
	next if ( $line[9] =~ /non_coding_transcript_exon_variant/ );
	next if ( $line[9] =~ /UTR/ );
	if ( $line[9] =~ /missense/ ) { $line[9] = "missense";}
	if ( $line[9] =~ /frameshift/ and $line[4] eq "-" ) { $line[9] = "frame_shift_ins"; }
	if ( $line[9] =~ /frameshift/ and $line[5] eq "-" ) { $line[9] = "frame_shift_del"; }
	if ( $line[9] =~ /splice/ and $line[7] eq "-" ) { $line[9] = "splice_site_ins"; }
	if ( $line[9] =~ /splice/ and $line[8] eq "-" ) { $line[9] = "splice_site_del"; }
	if ( $line[9] =~ /splice/ and $line[7] ne "-" and $line[8] ne "-" ) { $line[9] = "splice_site"; }
	if ( $line[9] =~ /stop_gained/ ) { $line[9] = "nonsense"; }
	if ( $line[9] =~ /stop_lost/ ) { $line[9] = "nonstop"; }
	
	# # expected fields: $gene, $chr, $strt, $stop, $ref, $var, $type,
    #       			$num_norm_refs, $num_norm_vars, $nvaf,
    #                   $num_tumr_refs, $num_tumr_vars, $tvaf,
    #                   $sample_name
	my $outString = join( "\t" , ( @line[3..9] , @line[99,100] , $line[104] , @line[102,103], $line[105], $line[1] ) );
	if ( $line[9] =~ /missense/) { $OUTmissense->print( $outString."\n" );} #gene site variant_class charger_pathogenicity
	else{$OUTtruncation->print( $outString."\n" );}
}
$IN1->close();

$OUTmissense->close();
$OUTtruncation->close();
