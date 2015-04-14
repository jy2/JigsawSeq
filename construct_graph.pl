#!/usr/bin/perl

# Construct de Bruijn graph from Jigsaw-Seq fastq data. 
# last modified: Apr-13-2015
# Developed by Jung-Ki Yoon

use strict;
use JigsawSeq;
use Benchmark ':hireswallclock';

our $t_begin = new Benchmark;
our $t_end;
our ($in_fname, $k_mer_len, $step_size, $out_fname,);

ReadArgument();

my %node; 
my $a_node = \%node;
my $num_node=my $num_read=my $num_filter_node=0;
print "[Report:construct_graph] input: $in_fname, k-mer: $k_mer_len, step_size: $step_size, output: $out_fname\n";

# Build graph
build_graph($in_fname);

$t_end = new Benchmark;
print "[Report:construct_graph] $num_read reads were processed; $num_node nodes were detected; Mem. Used = ", JigsawSeq::memcheck(), " Gb; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n";

# Recording
open(OUT, ">$out_fname");
for my $c_node (sort(keys %node)){
	$num_filter_node++;	
	unless ($num_filter_node % 5000000){
		$t_end = new Benchmark;
		print "[Process:construct_graph] $num_filter_node nodes were recorded; Mem. Used = ", JigsawSeq::memcheck(), " Gb; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n";
	}
	my $str_back_case="";
	my $str_back_case_depth="";
	my $num_edge=0;
	for my $c_edge (keys %{$a_node->{$c_node}}){
		next if (($c_edge eq "dep_node") || ($c_edge eq "num_edge"));
		$num_edge++; 
		$str_back_case .= ($c_edge . ",");
		$str_back_case_depth .= ($node{$c_node}{$c_edge} . ",");
	}
	print OUT join("\t", $c_node, $node{$c_node}{'dep_node'}, $num_edge, $str_back_case, $str_back_case_depth), "\n";
	$str_back_case="";
	$str_back_case_depth="";
}
close(OUT);
$t_end = new Benchmark;

print "[Report:construct_graph] $num_filter_node nodes were remained and recorded; Mem. Used = ", JigsawSeq::memcheck(), " Gb; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n\n";
exit;

sub build_graph(){
	my $fname = shift;
	open(IN, "<$fname") or die "[Error:construct_graph] Can't open $fname.\n";
#	print "[Process:construct_graph] $fname was opened.\n";
	while(<IN>){          # 1st line: header
		$num_read++;
		unless (($num_read % 100000)) {
			$t_end = new Benchmark;
			print "[Process:construct_graph] $num_read reads were processed; Mem. Used = ", JigsawSeq::memcheck(), " Gb; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n";
		}
		die "[Error:construct_graph] header was not detected. Probably input file was truncated or not fastq file.\n$_" unless (substr($_, 0, 1) eq "@");
		my $seq = <IN>;   # 2nd line; sequence
		chop($seq);
		<IN>;             # 3rd line; +
		<IN>;             # 4th line; score

		my $rev_seq = JigsawSeq::rev_comp($seq);
		next if ($rev_seq eq "Error"); # if the reverse complementary of this read could not be retrieved, skip the read. 
		# the reads having only A,T,C,G will be processed.

		# Index K-mer and build de Bruijn graph
		for (my $i=0; $i<=(length($seq)-$k_mer_len); $i++){
			my $cur_kmer = substr($seq, $i, $k_mer_len);
			my $back_node = substr($cur_kmer, $step_size, $k_mer_len-$step_size);
			my $back_case = substr($cur_kmer, -$step_size);
			my $front_node = substr($cur_kmer, 0, $k_mer_len-$step_size);
			$num_node++ if ($node{$front_node}{'dep_node'} eq "");
			$node{$front_node}{'dep_node'}++;		
			$node{$front_node}{$back_case}++;

			# for reverse complementary sequences
			$cur_kmer = substr($rev_seq, $i, $k_mer_len);	
			$back_node = substr($cur_kmer, $step_size, $k_mer_len-$step_size);
			$back_case = substr($cur_kmer, -$step_size);
			$front_node = substr($cur_kmer, 0, $k_mer_len-$step_size);
			$num_node++ if ($node{$front_node}{'dep_node'} eq "");
			$node{$front_node}{'dep_node'}++;		
			$node{$front_node}{$back_case}++;
		}
	}
	close(IN);
}


sub ReadArgument{
	my @command_line = @ARGV;               #command line argument

	#set default values;
	$k_mer_len = 120;
	$step_size = 3;

	#parse arguments
    while (defined(my $arg=shift(@ARGV))){
    	if (($arg eq "-I")||($arg eq "--input")){ $in_fname = shift @ARGV; next; }
    	if (($arg eq "-O")||($arg eq "--output")){ $out_fname = shift @ARGV; next; }
    	if (($arg eq "-k")||($arg eq "--kmer")){ $k_mer_len = shift @ARGV; next; }
    	if (($arg eq "-s")||($arg eq "--step")){ $step_size = shift @ARGV; next; }
    	PrintError("Can't identify the argument: $arg");
    }

    #check validity of arguments
    PrintError("Arguments -I, -O are required.") unless ( (defined $in_fname)&&(defined $out_fname));
	PrintError("Length of k-mer (--kmer) must be 40 <= k-mer <=150.") unless ((40<=$k_mer_len)&&($k_mer_len<=150));
	PrintError("Step size (--step) must be 1, 2, or 3.") unless ((1<=$step_size)&&($step_size<=3));
	PrintError("Length of k-mer (--kmer) must be divisible by step size (--step).") if ($k_mer_len % $step_size);
	PrintError("Input filename ($in_fname) does not exist.") unless (-e $in_fname);
}

sub PrintError{
	print "\n[Error:construct_graph] ", join("\n", @_), "\n\n";

	print "DESCRIPTION:\n    This script is to generate de Bruijn graph from short reads.\n\n";
	print "USAGE: construct_graph.pl -I [input:fastq format] -O [output:prefix]\n\n";
	print "OPTIONS:\n",
		"    -I, --input       FILE    Input file name of reads (fastq format) [Required]\n",
		"    -O, --output      FILE    Output file name (graph) [Required]\n",
		"    -k, --kmer        INT     Length of k-mer [Default: $k_mer_len]\n",
		"    -s, --step        INT     Step size for exploring the graph [Default: $step_size]\n";
#		"    -v, --verbose             Verbose-mode On\n"; 
#		"    -h, -?, --help       This help message\n",
	print "\n";
	exit -1;
}
