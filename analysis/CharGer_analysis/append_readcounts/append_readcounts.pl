#!/usr/bin/env perl

# Append/annotate readcounts info to CharGer output
# Authors: R. Jay Mashl, Kuan-Lin Huang @ WashU
# PanCanAtlas Germline
# v1 (10/2017, rjm): original version for single-base substitutions
# v2 (10/2017, kh): expanded to indels

#use strict;
use warnings;
use File::Basename;
use Data::Dumper;

my $myReadCountsSummary="";

if( scalar @ARGV != 2 ) {
    die "\nSyntax:  $0   charger_variants   readcounts_normals-tumors-master_table\n\n\n";
}

my ($chargerFile, $rcFile) = @ARGV;

# Lookup table
my %rcTable = ();
my %bestIndel = ();
my $cnt=0;
open(FH, "< $rcFile ") || die "ERROR: could not open $rcFile for reading\n\n";
while( my $line=<FH> ) {
    $cnt++;
    chomp $line;
    next if( $line =~ /^tag/);
    my ($tag, $chr, $pos, $refBase, $n_depth, $n_A_cnt, $n_C_cnt, $n_G_cnt, $n_T_cnt, $n_indel, $t_depth, $t_A_cnt, $t_C_cnt, $t_G_cnt, $t_T_cnt, $t_indel) = split /\t/, $line;
    my $loc = join("\t", $chr, $pos );

    $rcTable{$tag}{$loc}{'norm'}{'depth'}    = $n_depth;
    $rcTable{$tag}{$loc}{'norm'}{'A'}{'cnt'} = $n_A_cnt;
    $rcTable{$tag}{$loc}{'norm'}{'C'}{'cnt'} = $n_C_cnt;
    $rcTable{$tag}{$loc}{'norm'}{'G'}{'cnt'} = $n_G_cnt;
    $rcTable{$tag}{$loc}{'norm'}{'T'}{'cnt'} = $n_T_cnt;

    $rcTable{$tag}{$loc}{'tumor'}{'depth'}    = $t_depth;
    $rcTable{$tag}{$loc}{'tumor'}{'A'}{'cnt'} = $t_A_cnt;
    $rcTable{$tag}{$loc}{'tumor'}{'C'}{'cnt'} = $t_C_cnt;
    $rcTable{$tag}{$loc}{'tumor'}{'G'}{'cnt'} = $t_G_cnt;
    $rcTable{$tag}{$loc}{'tumor'}{'T'}{'cnt'} = $t_T_cnt;
    
    my $best_indel = "";
    my $best_indel_n_cnt = 0;
    my $best_indel_t_cnt = 0;
    if ( $n_indel ne "0" && $n_indel ne "-1" ){
	 my @n_indels  = split( /\|/ , $n_indel);
         foreach my $n_indel (@n_indels){
	      my @F = split( ",", $n_indel );
	      my $indel = $F[0]; my $indel_cnt = $F[1]; #print "Found ".$indel."; Count: ".$indel_cnt."\n";
              if ( $indel_cnt > $best_indel_n_cnt){
	          $best_indel = $indel; $best_indel_n_cnt = $indel_cnt;
              }
	 }
	 $bestIndel{$tag}{$loc}{'norm'} = $best_indel;
	 $rcTable{$tag}{$loc}{'norm'}{$best_indel}{'cnt'} = $best_indel_n_cnt;
    }
    if ( $t_indel ne "0" && $t_indel ne "-1" ){
         my @t_indels  = split( /\|/ , $t_indel);
         foreach my $t_indel (@t_indels){
              my ($indel, $indel_cnt) = split(",", $t_indel);
              if ( $indel eq $best_indel ){
                  $best_indel_t_cnt = $indel_cnt;
              }
         }
         $rcTable{$tag}{$loc}{'tumor'}{$best_indel}{'cnt'} = $best_indel_t_cnt;
    }
}
close(FH);
print "Done reading readcounts.\n";

my $out=*STDOUT;
my $header="";
my ($refBase, $altBase, $bestIndel, $normalDepth,  $normalRefCnt, $normalAltCnt,  $tumorDepth,  $tumorRefCnt, $tumorAltCnt);

open(FH, "< $chargerFile") || die "ERROR: could not open $chargerFile for reading\n\n";
while( my $line=<FH> ){
    chomp $line;
    if( $line =~ /^Sample/ ){
        $header = $line;
        print $out join("\t", $header, "refBase", "altBase", "normalDepth",  "normalRefCnt", "normalAltCnt", "tumorDepth", "tumorRefCnt", "tumorAltCnt"),"\n";
        next;
    }
    
# Sample	Genotype	HUGO_Symbol	Chromosome	Start	Stop	Reference	Alternate	Variant_Classification	HGVSg
    my ($sample, $chr, $start, $stop, $ref, $alt_charger, $varType, $HGVSg, $alt) = (split /\t/, $line)[ 0,3,4,5,6,7,8,9,87];
    my $donor = substr($sample, 0, 12);

    # Check for single substitution
    $ref = uc( $ref );
    $alt = uc( $alt );
    my $loc = join("\t", $chr, $start );
    # Default is case of missing entries
    ( $normalDepth, $normalRefCnt, $normalAltCnt, $tumorDepth, $tumorRefCnt, $tumorAltCnt ) = (-1) x 6;
    if( $start == $stop  && $ref =~ /^[ACGT]$/ &&  $alt =~ /^[ACGT]$/ ) {
        $refBase      = $ref;
        $altBase      = $alt;
        
        if( exists $rcTable{$donor}{$loc} ) {
            $normalDepth  = $rcTable{$donor}{$loc}{'norm'}{'depth'};
            $normalRefCnt = $rcTable{$donor}{$loc}{'norm'}{$refBase}{'cnt'};
            $normalAltCnt = $rcTable{$donor}{$loc}{'norm'}{$altBase}{'cnt'};
            
            $tumorDepth   = $rcTable{$donor}{$loc}{'tumor'}{'depth'};
            $tumorRefCnt  = $rcTable{$donor}{$loc}{'tumor'}{$refBase}{'cnt'};
            $tumorAltCnt  = $rcTable{$donor}{$loc}{'tumor'}{$altBase}{'cnt'};
        }
    # } elsif ($ref ne "-"  &&  $alt eq "-")c
    #        $altBase = $alt;
    #        $refBase = substr( $ref,0,1 ); #use first base for now
    #        $altBase = $alt;
    #        if( exists $rcTable{$donor}{$loc} ) {
    #             $normalDepth  = $rcTable{$donor}{$loc}{'norm'}{'depth'};
    #             $normalRefCnt = $rcTable{$donor}{$loc}{'norm'}{$refBase}{'cnt'};
    #             $normalAltCnt = $normalDepth - $normalRefCnt;

    #             $tumorDepth   = $rcTable{$donor}{$loc}{'tumor'}{'depth'};
    #             $tumorRefCnt  = $rcTable{$donor}{$loc}{'tumor'}{$refBase}{'cnt'};
    #             $tumorAltCnt = $tumorDepth - $tumorRefCnt;
    #        }
     } else{ # cases for all indels
            $refBase = substr( $ref,0,1 );
            if( exists $rcTable{$donor}{$loc} ) {
	         $normalDepth  = $rcTable{$donor}{$loc}{'norm'}{'depth'};
		 if ( exists $bestIndel{$donor}{$loc}{'norm'}){
		     $altBase = $bestIndel{$donor}{$loc}{'norm'};
                     $normalAltCnt = $rcTable{$donor}{$loc}{'norm'}{$altBase}{'cnt'};
		 }
	         if ( exists $rcTable{$donor}{$loc}{'norm'}{$refBase}{'cnt'}){
                     $normalRefCnt = $rcTable{$donor}{$loc}{'norm'}{$refBase}{'cnt'};
                 } elsif($normalDepth != -1 && $normalAltCnt != -1){
		     $normalRefCnt = $normalDepth - $normalAltCnt; #over-estimate ref count for insertions where we don't know the refBase
		 }

                 $tumorDepth   = $rcTable{$donor}{$loc}{'tumor'}{'depth'};
		 if ( exists $rcTable{$donor}{$loc}{'tumor'}{$altBase}{'cnt'} ) {
                     $tumorAltCnt = $rcTable{$donor}{$loc}{'tumor'}{$altBase}{'cnt'};
		 }
		 if ( exists $rcTable{$donor}{$loc}{'tumor'}{$refBase}{'cnt'}){
                     $tumorRefCnt = $rcTable{$donor}{$loc}{'tumor'}{$refBase}{'cnt'};
                 } elsif($tumorDepth != -1 && $tumorAltCnt != -1){
                     $tumorRefCnt = $tumorDepth - $tumorAltCnt; #over-estimate ref count for insertions where we don't know the refBase
                 }
	    }

     }
        
     print $out join("\t", $line, $refBase, $altBase, $normalDepth,  $normalRefCnt, $normalAltCnt,  $tumorDepth,  $tumorRefCnt, $tumorAltCnt),"\n";
}    
#            print join("\t", $line, $refBase, $altBase, $normalDepth,  $normalRefCnt, $normalAltCnt,  $tumorDepth,  $tumorRefCnt, $tumorAltCnt),"\n";
            
close(FH);
