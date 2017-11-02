#!/gsc/bin/perl
# Origin script and statistical tests by Mike Wendl @ WashU 
# further development by Kuan Huang @WashU 2016-17

#  CALCULATE P-VALUE ON A PER SITE BASIS FOR THE 11 SPECIFIC CANCER
#  GENES IDENTIFIED AS SIGNIFICANT BY THE BURDEN TEST
#
#  =================
#  PROGRAMMING NOTES
#  =================
#
#  (1) see hotspot_allelic_imbalance/pvalue_per_site_unified_background.pl
#      for the modified version of this code that is used for "hotspots
#      of allelic imbalance" analysis
#
#  (2) FISHER TEST
#
#      We are only testing whether tumor variant allele frequency (tvaf) is
#      "higher" than expected, not "higher or lower", so we need a 1-sided
#      test. We set-up the 2x2 table test as:
#
#                     | NORMAL   TUMOR |
#       --------------+----------------+---------------
#                     |                |
#       # REF ALLELES |   a        b   |   a + b
#                     |                |
#       # VAR ALLELES |   c        d   |   c + d
#                     |                |
#       --------------+----------------+---------------
#                     | a + c    b + d | a + b + c + d
#
#      whereby Statistics::FastFisher::one_tailed_non_symmetric_left becomes
#      the correct method to use, since it ramps "d" up and "b" down
#
#  (3) PERMUTATION TEST
#
#      We follow the method shown in script bootstrap_vaf_test.pl in directory
#      research/in_progress/ovarian_germline/vaf_differences/ for enumerating
#      the entire possible distribution of all tumor VAFs for the background
#
#      ACTUALLY: WE HAVE CORRECTED THE PERMUTATION TEST (12-MAY-14) AND THOSE
#      OLDER RESULTS MAY NOT BE STRICTLY CORRECT
#
#      variable: $bins_power  should be set at least to 3, i.e. 1000 bins, so
#                             that there is very little chance of any single
#                             variant event hitting the most extreme bin (and
#                             maybe then having a tail P-value of 0), for
#                             example, with 3 bins you'd need cases like
#                  992
#       TVAF =  -------- = 0.996  to actually realize this scenario
#                992 + 4
#
#  (4) FDR
#
#      We've used FDR many times and code is widely available, see for example
#
#       2011_pathscan/tsp_analysis/process_pathscan_file.pl
#       blood_clonal_expansion/statistics_blood_normal_test/complete_statistical_pipeline_matched.pl

############
#  SET UP  #
############

#__STANDARD PERL PACKAGES
   use strict;
   use warnings;
   use Time::Piece;

#__SPECIAL PACKAGES
   use Statistics::ConfidenceInterval 'binomial_wilson';
   use Statistics::FastFisher 'one_tailed_non_symmetric_left';
   use Statistics::CombinePvals;

###################
#  CONFIGURATION  #
###################

   my $date = localtime->strftime('%Y-%m-%d');

   if ($#ARGV != 2) {
      die "Usage: $#ARGV perl $0 <input file> <output file> <gene list file>\n\n";
   }

   my ($file, $outfile, $genefile) = @ARGV;
 
#__THRESHOLD FOR (MAXIMUM) CONFIDENCE INTERVAL FOR POINTS TO RETAIN IN THE
#  CALCULATION
#  (VALUES SUGGESTED BY COMPARING example_brca_vs_background_z_165.agr
#  TO example_brca_vs_background_z_165_filtered_0_1.agr IN DIRECTORY
#  choose_interv_size_4_exclusion)

   my $threshold_conf_interv = 0.12;

#__CONFIDENCE INTERVAL Z-SCORE (FILTER LOW-RELIABILITY POINTS BASED ON THIS)
   my $z = 1.65; # 90% INTERVAL
###my $z = 1.96; # 95% INTERVAL

#__PANCAN DATA FILE TO PROCESS
   # my $file = "files/pan8000_all_0.05truncation_AI_input.txt";
   # my $outfile = "files/20160502_pan8000_all_0.05truncation_AI_input_site_out.txt";
#__MUTATION FORM: SNPS OR INDELS (RELEVANT ONLY TO TRUNCATION)
   my $selected_mutation_form = "both";
#  my $selected_mutation_form = "snp";
#  my $selected_mutation_form = "indel";

#__RESOLUTION (NUMBER OF BINS) FOR EMPIRICAL NULL DISTRIBUTION OF TUMOR VAF
#
#  The number of bins is actually 10**$bins_power and the value is used directly
#  in 2 related ways: (1) as an "sprintf" formatting specification, it gives the
#  number of digits behind the decimal point that define the bins (and therefore
#  their total number), and (2) it is used to specify iteration over all bins
#  to make sure each has a value.

   my $bins_power = 3; #  1,000 BINS
###my $bins_power = 4; # 10,000 BINS

####################
#  PRE-PROCESSING  #
####################

#__LIST OF "GENES OF INTEREST" WHOSE EVENTS WILL BE TESTED AGAINST THE NULL
#
#  GENES SUGGESTED BY BURDEN TESTING TO BE FURTHER TESTED TO PRIORTIZE FIG 2
#
#  This is the final list from Mike McLellan and Charles Lu (12-NOV-14)
#
#  Date: Tue, 11 Nov 2014 16:49:44 -0600
#  From: Charles Lu <clu@genome.wustl.edu>
#  To: Mike Wendl <mwendl@genome.wustl.edu>
#  CC: Li Ding <lding@genome.wustl.edu>,
#      Mike McLellan <mmclella@watson.wustl.edu>,
#      Matthew Wyczalkowski <mwyczalk@genome.wustl.edu>,
#  Subject: Re: Validating significant genes using additional 1630 cases
#
#  Mike
#
#  See 34 genes below that are significant in burden test between pan12  vs
#  WHI1046 caucasians.  Mike M has removed all oncogenes.
   my $test_genes = {};
   #my $genefile = "files/04_30_2016_pan8000_MAF0.05_truncation_varOnly.vep.vcf_charger_pathogenic_gene_list.txt";
   open (F, $genefile) || die "cant open $genefile";
   while (<F>) {
      chomp;
      print "$_\n";
      $test_genes->{$_} = 1;
   }
   close(F);


#__ESTABLISH TYPE OF MUTATION DATA IN THE FILE (TRUNCATION OR MISSENSE)
   my $mutation_type;
   if ($file =~ /truncation/i) {
      $mutation_type = "truncation";
      warn "truncation mutations: you are analyzing mutations of " .
           "form '$selected_mutation_form'\n";
   } elsif ($file =~ /missense/i) {
      $mutation_type = "missense";
   } else {
      die "cannot infer mutation type from '$file'";
   }

#####################
#  MAIN PROCESSING  #
#####################

#__READ DATA AND PRE-FILTER INTO GENE-SPECIFIC AND BACKGROUND STRUCTS
   print "# ALLELIC IMBALANCE 'PER SITE' ANALYSIS\n#\n";
   print "# script: $0\n#\n";
   print "# file: $file\n#\n";
   open (F, $file) || die "cant open $file";
   open ( my $f_out, ">$outfile" ) or die "Cannot open $outfile: $!";
   my ($back_tumor_ref_counts, $back_tumor_var_counts, $rows) = ([], [], []);
   my ($gene_tumor_ref_counts, $gene_tumor_var_counts) = ([], []);
   while (<F>) {

   #__PRE-PROCESS AND PARSE: NVAF AND TVAF ARE IN PERCENT NOT FRACTIONS
      chomp;
      if ( /^HUGO/) {print $f_out "$_\tAI_event_type\tPERMUT_PVAL\tFISHER_PVAL\tCOMBINED_PVAL\tFDR\n";next;}
      
#     PARSER FOR CHARLES LU INPUT FILES
#
      my ($gene, $chr, $strt, $stop, $ref, $var, $type,
          $num_norm_refs, $num_norm_vars, $nvaf,
                      $num_tumr_refs, $num_tumr_vars, $tvaf,
                     $sample_name) = split /\t/;


   #__TRUNCATION MUTATIONS ARE FURTHER STRATIFIED BETWEEN SNPS AND INDELS
      if ($mutation_type eq "truncation" && $selected_mutation_form ne "both") {

      #__DETERMINE THE FORM OF THIS MUTATION (SNP OR INDEL)
         my $mutation_form;
         if ($type eq 'frame_shift_del' || $type eq 'frame_shift_ins' ||
             $type eq 'splice_site_del' || $type eq 'splice_site_ins') {
            $mutation_form = "indel";
         } elsif ($type eq 'nonsense' || $type eq 'splice_site' ||
                                         $type eq 'nonstop') {
            $mutation_form = "snp";
         } else {
            die "do not recognize form of mutation for line '$_'";
         }

      #__SKIP UNLESS THIS IS THE FORM WE WANT
         next unless $mutation_form eq $selected_mutation_form;
      }

   #__CHECK FOR "NA" IN THIS LINE
      next if $num_norm_refs eq "NA";
      next if $num_norm_vars eq "NA";
      next if $num_tumr_refs eq "NA";
      next if $num_tumr_vars eq "NA";

   #__DISCARD LOW RELIABILITY TUMOR EVENTS BASED ON DISTANCE TO INTERVAL EDGE
      my $nt = $num_tumr_refs + $num_tumr_vars;
      unless ($nt) {
         push @{$rows}, [$_, {"EVENT_DISCARDED" => "NO READS IN TUMOR"}];
         next;
      }
      my $pt = $num_tumr_vars / $nt;
      my $max_interval_tumor = binomial_wilson ($pt, $nt, $z);
      unless ($max_interval_tumor <= $threshold_conf_interv) {
         push @{$rows}, [
            $_,
            {"EVENT_DISCARDED" => "UNRELIABLE LEVEL OF TUMOR DATA: CONF INTERVAL $max_interval_tumor"}
         ];
         next;
      }

   #__DISCARD LOW RELIABILITY NORMAL EVENTS BASED ON DISTANCE TO INTERVAL EDGE
      my $nn = $num_norm_refs + $num_norm_vars;
      unless ($nn) {
         push @{$rows}, [$_, {"EVENT_DISCARDED" => "NO READS IN NORMAL"}];
         next;
      }
      my $pn = $num_norm_vars / $nn;
      my $max_interval_normal = binomial_wilson ($pn, $nn, $z);
      unless ($max_interval_normal <= $threshold_conf_interv) {
         push @{$rows}, [
            $_,
            {"EVENT_DISCARDED" => "UNRELIABLE LEVEL OF NORMAL DATA: CONF INTERVAL $max_interval_normal"}
         ];
         next;
      }

   #__PROCESS IF IT IS FROM A GENE OF INTEREST AND THEREFORE A "TEST EVENT"
      if (exists $test_genes->{$gene}) {
# print " TEST ", join (', ', ($gene, $num_norm_refs, $num_norm_vars, $num_tumr_refs, $num_tumr_vars)) , "\n";

      #__NOTE THAT WE DO NOT SKIP CASES WHERE TUMOR VAF IS SOMEWHAT LESS
      #  THAN NORMAL VAF BECAUSE WE ALSO HAVE TO TEST WHETHER CURRENT TUMOR
      #  VAF IS HIGHER THAN THE BACKGROUND TUMOR VAF
         #print join( "\t" , ($num_norm_refs, $num_tumr_refs,$num_norm_vars, $num_tumr_vars)."\n");
      #__FISHER 2X2 TABLE TEST OF WHETHER TVAF IS SIGNIFICANLY HIGHER THAN NVAF
         my $pval_fisher = one_tailed_non_symmetric_left (
            [
               [$num_norm_refs, $num_tumr_refs],
               [$num_norm_vars, $num_tumr_vars],
            ]
         );

      #__SAVE FISHER P-VALUE AND TVAF FOR PERMUTATION TEST
         push @{$rows}, [
            $_,
            {
               "FISHER_PVAL" => $pval_fisher,
               "TVAF" => $pt,
            }
         ];

      #__SAVE REF & VAR ALLELE COUNTS FOR MAKING EMPIRICAL DISTRIBUTION
         push @{$gene_tumor_ref_counts}, $num_tumr_refs;
         push @{$gene_tumor_var_counts}, $num_tumr_vars;

   #__OR IF THIS IS PART OF THE "BACKGROUND"
      } else {
# print " background ", join (', ', ($gene, $num_norm_refs, $num_norm_vars, $num_tumr_refs, $num_tumr_vars)) , "\n";

         push @{$rows}, [$_, {"BACKGROUND_EVENT" => ""}];

      #__SAVE REF & VAR ALLELE COUNTS FOR MAKING EMPIRICAL DISTRIBUTION
         push @{$back_tumor_ref_counts}, $num_tumr_refs;
         push @{$back_tumor_var_counts}, $num_tumr_vars;
# print "SAVING REFS $num_tumr_refs   VARS $num_tumr_vars\n";
      }
   }
   close (F);

#__CREATE RAW TUMOR VAF BACKGROUND (NULL) DISTRIBUTION BY FULL ENUMERATION
#
#  PERMUTATION TEST CREATES THE NULL DISTRIBUTION BY *POOLING* I.E. SATSIFIES
#  THE NULL HYPOTHESIS THAT THE TEST GENE HAS THE SAME DISTRIBUTION AS THE
#  BACKGROUND GENES. THE POOLING PROCESS FOLLOWS 2 STEPS
#
#     (1) EVERY OBSERVED TUMOR reference COUNT FROM THE TEST GENE IS COMBINED
#         WITH EVERY OBSERVED TUMOR variant COUNT FROM THE BACKGROUND GENES
#         TO GET ONE-HALF OF THE NULL DISTRIBUTION OF TUMOR VARIANT ALLELE FREQS
#     (2) EVERY OBSERVED TUMOR variant COUNT FROM THE TEST GENE IS COMBINED
#         WITH EVERY OBSERVED TUMOR reference COUNT FROM THE BACKGROUND GENES
#         TO GET OTHER HALF OF NULL DISTRIBUTION OF TUMOR VARIANT ALLELE FREQS
   print "BUILDING TUMOR VAF BACKGROUND DISTRIBUTION\n";
   my $format = join ('', '%.', $bins_power, 'f');
   my ($sample_space_size, $prob_mass) = (0, {});
#CHANGED
#  my $points_in_null = scalar @{$back_tumor_ref_counts} * scalar @{$back_tumor_var_counts};
   my $points_in_null = 2 * scalar @{$back_tumor_ref_counts} * scalar @{$gene_tumor_var_counts};
   print "   number of points in background: ", scalar @{$back_tumor_ref_counts}, "\n";
   print "   number of points in test gene set: ", scalar @{$gene_tumor_ref_counts}, "\n";
   print "   number of points in permutation set: $points_in_null\n";

#__BACKGROUND VARIANT TO TEST GENE REFERENCE
   foreach my $num_tumr_refs (@{$gene_tumor_ref_counts}) {
      foreach my $num_tumr_vars (@{$back_tumor_var_counts}) {

      #__CALCULATE TVAF AND BIN
      #
      #  NOTE: IT COULD HAPPEN BY CHANCE THAT A 0 REF COUNT FROM ONE EVENT
      #  IS COMBINED WITH A 0 VAR COUNT FROM ANOTHER: SKIP THESE
         my $nt = $num_tumr_refs + $num_tumr_vars;
         next unless $nt;
         my $pt = $num_tumr_vars / $nt;
         my $pt_bin = sprintf ($format, $pt);

      #__CATALOG
         $prob_mass->{$pt_bin}++;
         $sample_space_size++;
         unless ($sample_space_size % 100000) {
            my $pct = 100 * $sample_space_size / $points_in_null;
            warn "   ON EVENT $sample_space_size of $points_in_null ($pct %)\n"
         }
      }
   }

#__BACKGROUND REFERENCE TO TEST GENE VARIANT
   foreach my $num_tumr_refs (@{$back_tumor_ref_counts}) {
      foreach my $num_tumr_vars (@{$gene_tumor_var_counts}) {

      #__CALCULATE TVAF AND BIN
      #
      #  NOTE: IT COULD HAPPEN BY CHANCE THAT A 0 REF COUNT FROM ONE EVENT
      #  IS COMBINED WITH A 0 VAR COUNT FROM ANOTHER: SKIP THESE
         my $nt = $num_tumr_refs + $num_tumr_vars;
         next unless $nt;
         my $pt = $num_tumr_vars / $nt;
         my $pt_bin = sprintf ($format, $pt);

      #__CATALOG
         $prob_mass->{$pt_bin}++;
         $sample_space_size++;
         unless ($sample_space_size % 100000) {
            my $pct = 100 * $sample_space_size / $points_in_null;
            warn "   ON EVENT $sample_space_size of $points_in_null ($pct %)\n"
         }
      }
   }

#__NORMALIZE TO OBTAIN NULL PROBABILITY (MASS) DISTRIBUTION: FIND SMALLEST
   foreach my $pt_bin (sort _numerical_ keys %{$prob_mass}) {
      $prob_mass->{$pt_bin} /= $sample_space_size;
   }

#__MAKE SURE EACH BIN HAS AN ASSIGNED VALUE (ZERO IT IF NECESSARY)
#
#  FOLLOW EXACT FORMATTING AS ABOVE SO THAT YOU DO NOT CONFLATE SOMETHING
#  LIKE "0.01" WITH "0.010", THE FORMER ARISING FROM SKIPPING SPRINTF
   my $total_bins = 10**$bins_power;
   for (my $bin_number = 0; $bin_number <= $total_bins; $bin_number++) {
      my $pt = $bin_number / $total_bins;
      my $pt_bin = sprintf ($format, $pt);
      $prob_mass->{$pt_bin} = 0 unless exists $prob_mass->{$pt_bin};
   }

#__CALCULATE CDF FOR QUICK LOOKUP: ACCUMULATE FROM MOST EXTREME ON UP TO 1
#
#  AT THE VERY EXTREME END WE MAY HAVE CDF'S OF ZERO, I.E. IF THE MASSES IN
#  THE MOST EXTREME BINS ARE ZERO --- RESET THESE TO SOME NOMINALLY SMALL
#  VALUE --- BUT SUBTRACT THAT VALUE BACK OUT ONCE THE CDF STARTS ACCUMULATING
#  REAL VALUES

   my ($cdf, $ptail) = ({}, 0);
   foreach my $pt_bin (reverse sort _numerical_ keys %{$prob_mass}) {
      $cdf->{$pt_bin} = $prob_mass->{$pt_bin} + $ptail;
      $ptail = $cdf->{$pt_bin};
# print "$pt_bin   $prob_mass->{$pt_bin}   $cdf->{$pt_bin}\n";
   }

#__GO THROUGH ALL THE "TEST" DATA (THOSE WHICH ALREADY HAVE FISHER P-VALUE) AND
#  DO THE PERMUTATION TEST AND THEN ALSO COMBINE THE FISHER AND PERMUTATION
#  VALUES INTO AN OVERALL P-VALUE TAIL
   my ($p_test, $num_site_tests) = ({}, 0);
   foreach my $row_list (@{$rows}) {

   #__ACCESS INFO AND SKIP UNLESS THIS IS A TEST DATUM
      my $row           = $row_list->[0];
      my $row_info_hash = $row_list->[1];
      if (exists $row_info_hash->{'FISHER_PVAL'}) {

      #__PERMUTATION TEST P-VALUE
         my $tvaf = $row_info_hash->{'TVAF'};
         my $tvaf_bin = sprintf ($format, $tvaf);
         my $pval_permutation = $cdf->{$tvaf_bin};
         $row_info_hash->{'PERMUT_PVAL'} = $pval_permutation;

      #__DUMB REASSIGNMENT OF PVAL TO AVOID BUG WE STILL CANNOT FIND
         $row_info_hash->{'FISHER_PVAL'} = 0.999999999999 if $row_info_hash->{'FISHER_PVAL'} >= 1;
         $row_info_hash->{'PERMUT_PVAL'} = 0.999999999999 if $row_info_hash->{'PERMUT_PVAL'} >= 1;

      #__COMBINE FISHER AND PERMUTATION PVALS INTO SINGLE OVERALL P-VALUE
         my $obj = Statistics::CombinePvals->new ([
            $row_info_hash->{'FISHER_PVAL'},
            $row_info_hash->{'PERMUT_PVAL'}
         ]);
         my $pval_combined = $obj->fisher_chisq_transform;
         $row_info_hash->{'COMBINED_PVAL'} = $pval_combined;
         $num_site_tests++;
#print "$row\n";
#print "   computed tvaf $tvaf\n";
#print "   formatted tvaf $tvaf_bin\n";
#print "   permutation pval $pval_permutation\n";
#print "   fisher pval  $row_info_hash->{'FISHER_PVAL'}\n";
#print "   COMBINED PVAL  $row_info_hash->{'COMBINED_PVAL'}\n";
#print "\n";

      #__SAVE COMBINED P-VALUE WITH POINTER TO ITS ROW FOR DOWNSTREAM FDR
         push @{$p_test->{$pval_combined}}, $row;
      }
   }

#__CALCULATE FDR MAKING SURE P-VAL ORDER IS MOST EXTREME TO LEAST EXTREME
   print "FILE: $file\n";
   print "TEST RESULTS ONLY --- ORDERED BY FDR\n";
   print "\n";
   my ($results, $rank) = ({}, 1);

   my ($allgenes_sig_events, $allgenes_tot_events) = (0, 0);

   my ($brca1_sig_events, $brca2_sig_events) = (0, 0);
   my ($rad51d_sig_events, $palb2_sig_events) = (0, 0);
   my ($rad51c_sig_events, $bap1_sig_events) = (0, 0);
   my ($atm_sig_events, $brip1_sig_events) = (0, 0);
   my ($atr_sig_events, $fancm_sig_events) = (0, 0);
   my ($brca1_tot_events, $brca2_tot_events) = (0, 0);
   my ($rad51d_tot_events, $palb2_tot_events) = (0, 0);
   my ($rad51c_tot_events, $bap1_tot_events) = (0, 0);
   my ($atm_tot_events, $brip1_tot_events) = (0, 0);
   my ($atr_tot_events, $fancm_tot_events) = (0, 0);
   foreach my $pval (sort _numerical_ keys %{$p_test}) {
      foreach my $row (@{$p_test->{$pval}}) {
         my $fdr = $pval * $num_site_tests / $rank;
         $fdr = 1 if $fdr > 1;
#        print "$row: PVAL $pval X $num_site_tests / RANK $rank =  FDR $fdr\n";
##       print "$row	FDR $fdr\n";
#-------------- ADDED 12-NOV-14 TO JUST PRINT EVENTS FOR GENES SIGNIFICANT
         my @fields = split /\t/, $row;
         my $gene = $fields[13];
         if ($gene eq 'BRCA1') {
            $brca1_tot_events++;
            $brca1_sig_events++ if $fdr <= 0.15;
         }
         if ($gene eq 'BRCA2') {
            $brca2_tot_events++;
            $brca2_sig_events++ if $fdr <= 0.15;
         }
         if ($gene eq 'RAD51D') {
            $rad51d_tot_events++;
            $rad51d_sig_events++ if $fdr <= 0.15;
         }
         if ($gene eq 'PALB2') {
            $palb2_tot_events++;
            $palb2_sig_events++ if $fdr <= 0.15;
         }
         if ($gene eq 'RAD51C') {
            $rad51c_tot_events++;
            $rad51c_sig_events++ if $fdr <= 0.15;
         }
         if ($gene eq 'BAP1') {
            $bap1_tot_events++;
            $bap1_sig_events++ if $fdr <= 0.15;
         }
         if ($gene eq 'ATM') {
            $atm_tot_events++;
            $atm_sig_events++ if $fdr <= 0.15;
         }
         if ($gene eq 'BRIP1') {
            $brip1_tot_events++;
            $brip1_sig_events++ if $fdr <= 0.15;
         }
         if ($gene eq 'ATR') {
            $atr_tot_events++;
            $atr_sig_events++ if $fdr <= 0.15;
         }
         if ($gene eq 'FANCM') {
            $fancm_tot_events++;
            $fancm_sig_events++ if $fdr <= 0.15;
         }
         $allgenes_sig_events++ if $fdr <= 0.15;
         $allgenes_tot_events++;
#     my ($gene, $num_norm_refs, $num_norm_vars, $num_tumr_refs,
#         $num_tumr_vars) = @fields[13, $#fields - 3, $#fields - 2, $#fields - 7, $#fields - 6];

#######  print "GENE IS ----------->$gene<-------------\n";
#        if ($gene eq 'PALB2'){
#        if ($gene eq 'BRCA1' || $gene eq 'BRCA2' || $gene eq 'RAD51D' ||
#            $gene eq 'PALB2' || $gene eq 'BAP1'  || $gene eq 'RAD51C' ||
#            $gene eq 'ATM'   || $gene eq 'BRIP1' || $gene eq 'FANCM'  ||
#            $gene eq 'ATR') {
#-------------- FROM THE GENE-BASED ALLELIC IMBALANCE TEST
# print "\n";
         printf ("	%s	FDR %g\n", $row, $fdr);
# my $nvaf = 100 * $num_norm_vars / ($num_norm_refs + $num_norm_vars);
# my $tvaf = 100 * $num_tumr_vars / ($num_tumr_refs + $num_tumr_vars);
# print " ", join (', ', ($gene, $num_norm_refs, $num_norm_vars, $num_tumr_refs, $num_tumr_vars)) , "    NVAF = $nvaf   TVAF = $tvaf\n\n";
#--------------
#        }
#--------------
###      push @{$results->{$fdr}}, [$row, $pval]; # ORDER BY FDR
###      push @{$results->{$pval}}, [$row, $fdr]; # ORDER BY PVAL
         $results->{$row} = $fdr;
         $rank++;
      }
   }

# die "no full verbose output";

##################
#  POST PROCESS  #
##################

#__OUTPUT
   print "\nFULL VERBOSE OUTPUT FOR EVERY LINE --- ORIGINAL ORDER OF FILE\n";
   print "\n";
   foreach my $row_list (@{$rows}) {
      my $row           = $row_list->[0];
      my $row_info_hash = $row_list->[1];

   #__DEBUG OUTPUT
      print $f_out $row;
      if ($row =~ /HUGO_Symbol/){
         
         next;
      } 

      
      if (exists $row_info_hash->{'FISHER_PVAL'}) {
         print $f_out "\tTEST_EVENT";
         foreach my $key (reverse sort keys %{$row_info_hash}) { # handle discarded and non discarded events
            printf $f_out ("\t%g", $row_info_hash->{$key}) unless $key eq 'TVAF';  
         }
         printf $f_out ("\t%g", $results->{$row}) if exists $results->{$row};
      } else {
         foreach my $key (reverse sort keys %{$row_info_hash}) { 
            print $f_out "\t$key: $row_info_hash->{$key}\tNA\tNA\tNA\tNA";
         }
      }
      
      # foreach my $key (reverse sort keys %{$row_info_hash}) { # handle discarded and non discarded events
      #    print "$key\n";
      #    if (exists $row_info_hash->{'FISHER_PVAL'}) {
      #       #printf $f_out ("	%s	%g", $key, $row_info_hash->{$key}) unless $key eq 'TVAF';
      #       printf $f_out ("\t%g", $row_info_hash->{$key}) unless $key eq 'TVAF';
      #    } else {
      #       print $f_out "\t$key: $row_info_hash->{$key}";
	     #  }
      # }
      # print calculated FDR to the end of the row
#     print "	FDR	$results->{$row}" if exists $results->{$row};
      print $f_out "\n";
   }
##       my $pt_bin = sprintf ("%.3f", $pt);

   print "TRUNCATIONS: TOTAL AND NUMBER OF SIGNIFICANT (WITHIN 15% FDR)\n";
   print "  BRCA1:  $brca1_sig_events / $brca1_tot_events\n";
   print "  BRCA2:  $brca2_sig_events / $brca2_tot_events\n";
   print "  RAD51D: $rad51d_sig_events / $rad51d_tot_events\n";
   print "  PALB2:  $palb2_sig_events / $palb2_tot_events\n";
   print "  RAD51C: $rad51c_sig_events / $rad51c_tot_events\n";
   print "  BAP1:   $bap1_sig_events / $bap1_tot_events\n";
   print "  ATM:    $atm_sig_events / $atm_tot_events\n";
   print "  BRIP1:  $brip1_sig_events / $brip1_tot_events\n";
   print "  ATR:    $atr_sig_events / $atr_tot_events\n";
   print "  FANCM:  $fancm_sig_events / $fancm_tot_events\n";
   print "\n";
   print "  ALL GENES: $allgenes_sig_events / $allgenes_tot_events\n";

################################################################################
#                                                                              #
#                            S U B R O U T I N E S                             #
#                                                                              #
################################################################################

sub _numerical_ {$a <=> $b}
