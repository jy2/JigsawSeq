#!/usr/bin/perl

# Detect initial/terminal seeds.
# Requirement: bwa 0.7.5a-r405
# last modified: Apr-22-2015
# Developed by Jung-Ki Yoon

use strict;
use JigsawSeq;
use Benchmark ':hireswallclock';

our ($in_fname, $vector_seq, $k_mer_len, $step_size, $cut_seed, $min_seed, $num_thread, $out_fname,);
ReadArgument();

our $kmer_fname = "TMP_$in_fname\.kmer.fa";
my $len_ref_seed = int($k_mer_len*1.5);				# allow a half size of k-mer as insertion on seeds
our $ref_seed_fname = "TMP_k$k_mer_len\_" . $vector_seq;
our $seeds_sam_fname = "TMP_$in_fname\.seed.sam";
my $t_begin = new Benchmark;
print "[Report:detect_seeds] input: $in_fname vector: $vector_seq k-mer: $k_mer_len, step_size: $step_size, cut_seed: $cut_seed, min_seed: $min_seed, output: $out_fname\n";

Kmer2Fa();
Vector2Fa();

print "[Report:detect_seeds] Align nodes to vector using bwa\nBWA----------------------\n";
system("./bwa index $ref_seed_fname > $ref_seed_fname\.log");
system("./bwa mem -t $num_thread -v 1 -O2 -E1 $ref_seed_fname $kmer_fname > $seeds_sam_fname");
my $t_end = new Benchmark;
print "-----------------------BWA\n", 
      "[Report:detect_seeds] Alignment was completed; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n";

AnalyzingSeeds();
exit;

#----------------------------
sub Kmer2Fa{
	open(IN, "<$in_fname") or die "[Error:detect_seeds] Can't open $in_fname.\n";
	open(OUT, ">$kmer_fname");
	my $num_lines=0;
	while(<IN>){ 
		if (substr($_, 0, 1) eq "#"){ next;	}
		my ($str, $DP,) = split /\s+/, $_;
		print OUT ">$num_lines\_$DP\n";
		print OUT "$str\n";
		$num_lines++;
	}
	close(IN);
	close(OUT);
	print "[Report:detect_seeds] $num_lines nodes were ready for running dectect seeds.\n";
}

sub Vector2Fa{
	open(IN, "<$vector_seq") or die "[Error:detect_seeds] Can't open $vector_seq.\n";
	open(OUT, ">$ref_seed_fname");
	<IN>;
	my ($seq,)=split /\s+/, <IN>;
	print OUT ">initial\n", uc(substr($seq,-$len_ref_seed)), "\n";
	print OUT ">terminal\n", uc(substr($seq,0,$len_ref_seed)), "\n";
	close(OUT);
	close(IN);
}

sub AnalyzingSeeds {
	my %init_SEEDS;
	my %term_SEEDS;
	my $num_lines=my $num_align_init=my $num_init=my $final_init=my $num_align_term=my $num_term=my $final_term=0;
	my $max_DP_init=my $max_DP_term=-1;
	open(IN, "<$seeds_sam_fname") or die "[Error:detect_seeds] Can't open $seeds_sam_fname.\n";
	open(OUT, ">$out_fname");
	while(<IN>){ 
		next if (substr($_, 0, 1) eq "@");
		$num_lines++;
		if ($num_lines % 1000000 == 0) { print "[Process:detect_seeds] $num_lines lines were processed.\n";}
		my ($qname, $flag, $rname, $pos, $mapQ, $CIGAR, $rnext, $pnext, $tlen, $str, $qual, $NM, $MD, $AS, $XS,) = split /\s+/, $_;
		next if ($flag != 0);
		my @d = split /[MIDNSHPX=]/, $CIGAR;
		my @l = split /[0-9]+/, $CIGAR;
		shift @l;
		my $sum_d=0;
		if ($rname eq "initial"){
			$num_align_init++;
			($qname, my $DP) = split /_/, $qname;
			for(my $i=0; $i<=$#d; $i++){
				$sum_d+=$d[$i];
			}
			if (($pos + $sum_d) == ($len_ref_seed+1) ) {
		    	if ($l[$#l] eq "S"){
					next;
			    }
		    	$init_SEEDS{$str} = $DP;
		    	if ($max_DP_init<$DP) {$max_DP_init=$DP;}
			}
		}elsif ($rname eq "terminal"){
			$num_align_term++;
			($qname, my $DP) = split /_/, $qname;
			for(my $i=0; $i<=$#d; $i++){
				if (($l[$i] eq "M")||($l[$i] eq "I")){
					$sum_d+=$d[$i];
				}
			}
			if (($pos == 1)&&($sum_d==($k_mer_len-$step_size))){
				$term_SEEDS{$str} = $DP;
				if ($max_DP_term<$DP) {$max_DP_term=$DP;}
			}
		}
	}
	close(IN);

	foreach my $str (keys %init_SEEDS){
		my $DP = $init_SEEDS{$str};
		$num_init++;
		next if ( (($max_DP_init/$DP)>$cut_seed)||($DP<$min_seed) );
		print OUT ">initial $DP\n$str\n";
		$final_init++;
	}	
	foreach my $str (keys %term_SEEDS){
		my $DP = $term_SEEDS{$str};
		$num_term++;
		next if ( (($max_DP_term/$DP)>$cut_seed)||($DP<$min_seed) );
		print OUT ">terminal $DP\n$str\n";
		$final_term++;
	}	
	close(OUT);

	my $t_end = new Benchmark; 
	print "[Report:detect_seeds] Cutoff ratio of seeds: $cut_seed\tMax_inti_DP: $max_DP_init\tMax_term_DP: $max_DP_term\n";
	print "[Report:detect_seeds] $num_lines K-mers were processed.\n";
	print "[Report:detect_seeds] $num_align_init K-mers were alinged to initial seq. of backbone.\t$num_init were well matched.\t$final_init were passed the cutoff(=initial seeds).\n";
	print "[Report:detect_seeds] $num_align_term K-mers were aligned to terminal seq. of backbone.\t$num_term were well matched.\t$final_term were passed the cutoff(=terminal seeds).\n";
	print "                      Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n\n";
}

sub ReadArgument{
	my @command_line = @ARGV;               #command line argument

	#set default values;
	$k_mer_len = 120;
	$step_size = 3;
	$cut_seed = 100;
	$min_seed = 100;
	$num_thread = 3;

	#parse arguments
    while (defined(my $arg=shift(@ARGV))){
    	if (($arg eq "-I")||($arg eq "--input")){ $in_fname = shift @ARGV; next; }
    	if (($arg eq "-V")||($arg eq "--vector")){ $vector_seq = shift @ARGV; next; }
    	if (($arg eq "-O")||($arg eq "--output")){ $out_fname = shift @ARGV; next; }
    	if (($arg eq "-k")||($arg eq "--kmer")){ $k_mer_len = shift @ARGV; next; }
    	if (($arg eq "-s")||($arg eq "--step")){ $step_size = shift @ARGV; next; }
    	if ($arg eq "--cut_seed") { $cut_seed = shift @ARGV; next; }
    	if ($arg eq "--min_seed") { $min_seed = shift @ARGV; next; }
    	if (($arg eq "-t")||($arg eq "--thread")){ $num_thread = shift @ARGV; next; }
       	if (($arg eq "-h")||($arg eq "-?")||($arg eq "--help")){ printError(); next; }
    	PrintError("Can't identify the argument: $arg");
    }

    #check validity of arguments
    PrintError("Arguments -I, -O, -V are required.") unless ( (defined $in_fname)&&(defined $vector_seq)&&(defined $out_fname) );
	PrintError("Length of k-mer (--kmer) must be 40 <= k-mer <=150.") unless ((40<=$k_mer_len)&&($k_mer_len<=150));
	PrintError("Step size (--step) must be 1, 2, or 3.") unless ((1<=$step_size)&&($step_size<=3));
	PrintError("Length of k-mer (--kmer) must be divisible by step size (--step).") if ($k_mer_len % $step_size);
	PrintError("Input filename ($in_fname) does not exist.") unless (-e $in_fname);
	PrintError("Input filename ($vector_seq) does not exist.") unless (-e $vector_seq);
}

sub PrintError{
	print "\n[Error:detect_seeds] ", join("\n", @_), "\n\n";

	print "DESCRIPTION:\n    detect_seeds.pl is to detect initial and terminal seeds prior to explore graph.\n\n";
	print "USAGE: detect_seeds.pl -I [input:graph] -O [output:fasta] -V [input:vector]\n\n";
	print "OPTIONS:\n",
		"    -I, --input       FILE    Input file name of graph [Required]\n",
		"    -O, --output      FILE    Output filen name (fasta format) [Required]\n",
		"    -V, --vector      FILE    Input filename of vector sequence (fasta format) [Required]\n",
		"    -k, --kmer        INT     Length of k-mer [Default: $k_mer_len]\n",
		"    -s, --step        INT     Step size for exploring the graph [Default: $step_size]\n",
		"    --cut_seed        INT     Cutoff ratio for seeds [Default: $cut_seed]\n",
	        "    --min_seed        INT     Miniumn depth of seeds [Default: $min_seed]\n",
		"    -t, --thread      INT     Number of threads [Default: $num_thread]\n";
#		"    -h, -?, --help       This help message\n",
	print "\n";
	exit -1;
}
