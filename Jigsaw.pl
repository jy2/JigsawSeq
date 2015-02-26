#!/usr/bin/perl

use strict;
use JigsawSeq;
use Benchmark ':hireswallclock';

our $version = '$Version: r2.01$';
our $LAST_MODIFIED = 'LastModified: Feb-26-2014';
    	
our $verbose = 0;
our ($input_F, $input_R, $vector_seq, $output_prefix); 
our ($exp_contig_size, $k_mer_len, $step_size, $min_depth, $cut_edge, $cut_seed, $cut_CV, $read_length);
our ($bin_size, $num_thread);
my $t_begin = new Benchmark;

# Parse arguments
ReadArgument();

# Build de Bruijn graph
our $num_files = 0;
our $prefix = "TMP_input_$output_prefix";
BuildGraph();

# Detect initial/terminal seeds
system("./detect_seeds.pl $output_prefix\.graph.clean.sort $vector_seq $k_mer_len $step_size $cut_seed $num_thread $output_prefix\.seeds.fa");

# Search candidate contigs
system("./explore_graph.pl $output_prefix\.graph.clean.sort $output_prefix\.seeds.fa $k_mer_len $step_size $exp_contig_size $output_prefix\.candidates.fa");

# Map raw reads to candidates and find contigs
system("./mapping_contigs.pl $output_prefix\.candidates.fa $input_F $input_R $k_mer_len $num_thread $output_prefix\.candidates");
system("./select_contigs.pl $output_prefix\.candidates.fa $output_prefix\.candidates.DP $k_mer_len $step_size $cut_CV $read_length $output_prefix\.contigs");
DelTemp();
exit;

my $t_end=new Benchmark;
print "[Report:main] Whole process of JigsawSeq analysis were successfully completed; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n\n";
exit;

#----------------------------
sub ReadArgument{
	my @command_line = @ARGV;               #command line argument

	#set default values;
	$verbose = 0;
	$k_mer_len = 120;
	$step_size = 3;
	$min_depth = 2;
	$cut_seed = 100;
	$cut_edge = 150;
	$cut_CV = 0.2163;
	$num_thread = 3;
	$read_length = 150;
	$bin_size = 15000000;

	#parse arguments
    while (defined(my $arg=shift(@ARGV))){
    	if (($arg eq "-v")||($arg eq "--verbose")){ $verbose = 1; next; } # verbose mode on
       	if (($arg eq "-h")||($arg eq "-?")||($arg eq "--help")){ printError(); next; }
    	if (($arg eq "-F")||($arg eq "--forward")){ $input_F = shift @ARGV; next; }
    	if (($arg eq "-R")||($arg eq "--reverse")){	$input_R = shift @ARGV; next; }
    	if (($arg eq "-V")||($arg eq "--vector")){ $vector_seq = shift @ARGV; next; }
    	if (($arg eq "-O")||($arg eq "--output")){ $output_prefix = shift @ARGV; next; }
    	if (($arg eq "-L")||($arg eq "--length")){ $exp_contig_size = shift @ARGV; next; }
    	if (($arg eq "-k")||($arg eq "--kmer")){ $k_mer_len = shift @ARGV; next; }
    	if (($arg eq "-s")||($arg eq "--step")){ $step_size = shift @ARGV; next; }
    	if (($arg eq "-m")||($arg eq "--min_depth")){ $min_depth = shift @ARGV; next; }
    	if (($arg eq "-C")||($arg eq "--cut_CV")){ $cut_CV = shift @ARGV; next; }
    	if ($arg eq "--cut_edge") { $cut_edge = shift @ARGV; next; }
    	if ($arg eq "--cut_seed") { $cut_seed = shift @ARGV; next; }
    	if ($arg eq "--read_length") { $read_length = shift @ARGV; next; }
    	if (($arg eq "-t")||($arg eq "--thread")){ $num_thread = shift @ARGV; next; }
    	if (($arg eq "-b")||($arg eq "--bin_size")){ $bin_size = shift @ARGV; next; }
    	PrintError("Can't identify the argument: $arg");
    }

    #check validity of arguments
    PrintError("Arguments -F, -R, -V, -O, -L are required.") unless ( (defined $input_F)&&(defined $input_R)&&(defined $vector_seq)&&(defined $exp_contig_size)&&(defined $output_prefix) );
	PrintError("Length of k-mer (--kmer) must be 40 <= k-mer <=150.") unless ((40<=$k_mer_len)&&($k_mer_len<=150));
	PrintError("Step size (--step) must be 1, 2, or 3.") unless ((1<=$step_size)&&($step_size<=3));
	PrintError("Length of k-mer (--kmer) must be divisible by step size (--step).") if ($k_mer_len % $step_size);
	PrintError("Number of reads in a bin (--bin_size) must be divisible by 4.") if ($bin_size % 4);
	PrintError("Input filename ($input_F) does not exist.") unless (-e $input_F);
	PrintError("Input filename ($input_R) does not exist.") unless (-e $input_R);
	PrintError("Input filename ($vector_seq) does not exist.") unless (-e $vector_seq);
}

sub PrintError{
	print "\n[Error:main] ", join("\n", @_), "\n\n";

	print "DESCRIPTION:\n    Jigsaw.pl is the main module to analyze JigsawSeq data\n";
	print "    Developed by Jung-Ki Yoon (dr.jkyoon\@gmail.com)\n    $LAST_MODIFIED\n\n";
	print "USAGE: Jigsaw.pl -F [input:forward] -R [input:reverse] -V [input:vector] -L [input: expected length] -O [output:prefix]\n\n";
	print "OPTIONS:\n",
		"    -F, --forward     FILE    Input filename of forward reads (fastq format) [Required]\n",
		"    -R, --reverse     FILE    Input filename of reverse reads (fastq format) [Required]\n",
		"    -V, --vector      FILE    Input filename of vector sequence (fasta format) [Required]\n",
		"    -O, --output      STR     Prefix of output file [Required]\n",
		"    -L, --length      INT     Expected length of contig [Required]\n",
		"    -k, --kmer        INT     Length of k-mer [Default: 120]\n",
		"    -s, --step        INT     Step size for exploring the graph [Default: 3]\n",
		"    -m, --min_depth   INT     Minimum depth of nodes and edges [Defulat: 2]\n",
		"    -C, --cut_CV      FLOAT   Cutoff for coefficient of variation [Default: 0.2163]\n",
		"    --cut_edge        INT     Cutoff ratio for edges [Default: 100]\n",
		"    --cut_seed        INT     Cutoff ratio for seeds [Default: 150]\n",
		"	 --read_length     INT     Read length of raw reads [Default: 150]\n",
		"    -t, --thread      INT     Number of threads [Default: 3]\n",
		"    -b, --bin_size    INT     Number of reads in a bin [Default: 15000000]\n";
#		"    -h, -?, --help       This help message\n",
	print "\n";
	exit -1;
}


sub BuildGraph{

	#Chop Fastq
	ChopFastq();

	#Construct graph for each chopped fastq
	for (my $i=1; $i<=$num_files; $i++){
		system("./construct_graph.pl $prefix.$i\.fastq $k_mer_len $step_size $prefix.$i\.graph");
	}

	#Merge graphs
	MergeGraph();

	#Clean up and sort graph
	system("./cleanup_graph.pl $output_prefix\.graph $min_depth $cut_edge $output_prefix\.graph.clean");
}


sub ChopFastq{
	print "[Report:main] Start cutting input files (bin size = $bin_size).\n";
	my $num_lines=0;

	open(IN, "<$input_F") or die "[Error] Can't open $input_F";
	while(<IN>){
		if (($num_lines % $bin_size) == 0){
			close(OUT) if ($num_files > 0);
			$num_files++;
			open(OUT, ">$prefix\.$num_files\.fastq");
		}
		print OUT $_;
		$num_lines++;
	}
	close(IN);
	open(IN, "$input_R") or die "[Error] Can't open $input_R";
	while(<IN>){
		if (($num_lines % $bin_size) == 0){
			close(OUT);
			$num_files++;
			open(OUT, ">$prefix\.$num_files\.fastq");
		}	
		print OUT $_;
		$num_lines++;
	}
	close(IN);
	close(OUT);
	print "[Report:main] $num_files temporary files are ready.\n\n";
}

sub MergeGraph{
	if ($num_files > 1){
		system("./merge_graph.pl $prefix\.1.graph $prefix\.2.graph $prefix\.1.merged");
		for (my $i=2; $i<$num_files; $i++){
			system("./merge_graph.pl $prefix\." . ($i-1) . ".merged $prefix\." . ($i+1) . ".graph $prefix\." . ($i) . ".merged");
		}
		print "[Report:main] $prefix\." . ($num_files-1) . ".merged -> $output_prefix\.graph\n";
		system("mv $prefix\." . ($num_files-1) . "\.merged $output_prefix\.graph");
	} else{
		print "[Report:main] $prefix\.1.merged -> $output_prefix\.graph\n";
		system("mv $prefix\.1.graph $output_prefix\.graph");
	}
}

sub DelTemp{
	if ($verbose==0){
		`rm $prefix.\*`;
		`rm $output_prefix\.graph`;
		`rm $output_prefix\.graph.clean`;
		`rm TMP_$output_prefix\.graph.clean.sort.kmer.fa`;
		`rm TMP_$output_prefix\.graph.clean.sort.seed.sam`;
		`rm TMP_k$k_mer_len\_$vector_seq\*`;
		`rm $output_prefix.candidates.fa.\*`;
		`rm $output_prefix.candidates.DP\*`;
		`rm TMP_$output_prefix.candidates.sam`;
		`rm TMP_$output_prefix.candidates.exact.\*`;
		`rm TMP_$output_prefix.candidates.DP.stat`;
		print "[Report:main] Temporary files were deleted.\n\n";
	}
}
