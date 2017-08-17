# Kuan Huang @ WashU 2017 Aug1
# reference: # https://github.com/zhuochenbioinfo/VcfStat
# filter variants based on a specific frequency using the cohort freq
# update AN/AC fields

use strict;
use warnings;
use Getopt::Long;

my($vcf,$keepList,$all,$CHROM);

my $usage = "USAGE:\nperl $0 --vcf <vcf> --out <out>\n";
$usage .= "<vcf> is the input vcf fle. [Necessary].\n";
#$usage .= "<out> is the output file. [Necessary]\n";

GetOptions(
	"vcf=s" => \$vcf,
	#"out=s" => \$out,
	"chr=s" => \$CHROM,
) or die $usage;

die $usage unless(defined $vcf);
my $log = $vcf."cohortAF.filter.log";

# add in a module to open both vcf and vcf.gz
if ($vcf =~ /.gz$/) {
	open(IN, "gunzip -c $vcf |") || die "can’t open pipe to $vcf";
}
else {
	open(IN, $vcf) || die "can’t open $vcf";
}

# open(IN,"<$vcf") or die $!;
open(LOG,">$log");

# arbitrary cut-off for now
my $sample_size = 9401;
my $AF_threshold = 0.05;

my @samples = ();
my @keepRanks = ();
my $num_pass_alleles = 0;
my $filtered_alleles = 0;
my $filtered_vars = 0;
my $nonexisting_alleles = 0;
my $nonexisting_vars = 0;

while(<IN>){
	chomp;
	if($_ =~ /^##/){
		print $_."\n";
		next;
	}
	my($chr,$pos,$id,$ref,$alts_join,$qual,$filter,$info,$format,@datas) = split/\t/;
	my @alts = split/,/,$alts_join;
	my @alleles = ($ref,@alts);
	my @alleles_count = (0 x $#alleles);
	my @pass_alleles = (); 
	my @pass_alleles_index = (0); #keeping the ref index
	my @pass_alleles_count = (); 
	my @pass_alleles_freq = (); 

	if($_ =~ /^#CHROM/){
		print $_."\n";		
		#print "##samplenum=$num\n";
		#print OUT "CHROM\tPOS\tREF\tALT\tFILTER\tALLELEnum\tHETnum\tNAnum\tCOVfreq\tALLELEfreq\n";
		next;
	}

	foreach my $spot(@datas){
		if($spot =~ /^(\d+)\/(\d+)/){ # add count for each allele here
			$alleles_count[$1]++;
			$alleles_count[$2]++;
		}
	}
	
	# keep only alleles that are rare in the cohort
	for(my $i = 1; $i < @alleles; $i++){
		my $count = $alleles_count[$i];
		if (!defined($count)){
			$nonexisting_alleles++;
			next;
		}
		
		my $AF = $count/$sample_size;
		if ( $AF < $AF_threshold ){
			$num_pass_alleles++;
			push @pass_alleles, $alleles[$i];
			push @pass_alleles_index, $i;
			push @pass_alleles_count, $alleles_count[$i];
			push @pass_alleles_freq, $AF;
			#print $AF."(AF)\t"
		} else {
			$filtered_alleles++;
		}
		
	}

	# modify the geno field;
	# reorder the index to the new index; replace essentially
	for(my $i=1; $i < @pass_alleles_index; $i++){
		my $old_index = $pass_alleles_index[$i];
		foreach my $spot(@datas){
			$spot =~ s/^$old_index\//$i\//;
			$spot =~ s/\/$old_index/\/$i/;
		}
	}

	my $count_of_pass_alleles = scalar @pass_alleles;
	# only print the line if we have more than one passed allele
	if($count_of_pass_alleles > 0){
			# update alternative allele
			my $alts_join = join(",",@pass_alleles);
			# update AC, AN, AF in info field of vcf; AF=0.5,0.5;MLEAC=1,1;MLEAF=0.5,0.5;AN=326;AC=169,23
			my $AC = join(",",@pass_alleles_count);
			my $AF = join(",",@pass_alleles_freq);
			$info =~ s/"AF=.*;"/"AF=".$AF.";"/;
			$info =~ s/"AC=.*;"/"AC=".$AC.";"/;
			$info =~ s/"AN=.*;"/"AF=".$sample_size.";"/;
			print join("\t",($chr,$pos,$id,$ref,$alts_join,$qual,$filter,$info,$format,@datas))."\n";
	} else{
		$filtered_vars++;
		#print LOG "Filtered variant (including nonexisting variants): ".join(",",($chr,$pos,$id,$ref,$alts_join))."\n";
	}		
	
}

print LOG "Pass alleles: $num_pass_alleles\n\n";
print LOG "Nonexisting alleles: $nonexisting_alleles\n";
print LOG "Filtered alleles: $filtered_alleles\n";
print LOG "Total filtered variants: $filtered_vars\n";
my $total_filtered_alleles = $nonexisting_alleles + $filtered_alleles;
print LOG "Total filtered alleles: $total_filtered_alleles\n";
close LOG;
