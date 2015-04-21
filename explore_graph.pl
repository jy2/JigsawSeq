#!/usr/bin/perl

# Search candidates of contigs from the de bruijn graph with seeds.
# last modified: Apr-13-2015
# Developed by Jung-Ki Yoon

use strict;
use JigsawSeq;
use Benchmark ':hireswallclock';

our ($in_fname, $in_seeds, $k_mer_len, $step_size, $exp_contig_size, $out_fname,);
ReadArgument();

our $t_begin = new Benchmark;
our $t_end;
our %node; 
our %e_DP;
our $a_node = \%node; 
our $num_node=0;
our @initial_nodes; 
our @terminal_nodes;
our $max_stroll = int(($exp_contig_size + 50 + $k_mer_len)/$step_size);  # allow a half of expected contig length as insertion size.

print "[Report:explore_graph] input: $in_fname seeds: $in_seeds k-mer_len: $k_mer_len step_size: $step_size expected_contig_size: $exp_contig_size max_num_stroll: $max_stroll output: $out_fname\n";

LoadSeeds();
LoadGraph();

# Explore graph with each initial seed.
my %contigs;
for(my $i=0; $i<=$#initial_nodes; $i++){
	print "[Process:explore_graph] Searching with ",$i+1,"-th initial nodes: $initial_nodes[$i]\n";
	explore_contig($i);
	print "[Report:explore_graph] ", scalar(keys %contigs), " contigs were detected.\n";
}

# Print out (fasta + result)
my $fa_fname = $out_fname . ".fa";
my $result_fname = $out_fname . ".result";
my $num_contigs=0;
open(OUT, ">$result_fname");
print OUT "#number_of_initial\t", $#initial_nodes+1, "\n";
for(my $i=0; $i<=$#initial_nodes; $i++){ print OUT "#initial_node\_$i\t$initial_nodes[$i]\n"; }
print OUT "#number_of_terminal\t", $#terminal_nodes+1, "\n";
for(my $i=0; $i<=$#terminal_nodes; $i++){ print OUT "#terminal_node\_$i\t$terminal_nodes[$i]\n"; }
print OUT join("\t", "#contig_name", "length", "depth", "initial_node", "terminal_node", "sequence"), "\n";
open(FA, ">$fa_fname");
for my $key (keys %contigs){
	$num_contigs++;
	print FA ">$out_fname\_$num_contigs\n", $initial_nodes[$contigs{$key}{'n_init'}], $key, $terminal_nodes[$contigs{$key}{'n_term'}], "\n";
	print OUT join("\t", "$out_fname\_$num_contigs", length($key), $contigs{$key}{'num'}, $contigs{$key}{'n_init'}, $contigs{$key}{'n_term'}, $key), "\n";
}
close(FA);
close(OUT);

$t_end = new Benchmark;
print "[Report:explore_graph] $fa_fname and $result_fname were recorded\n", 
      "                       Mem. Used = ", JigsawSeq::memcheck(), " Gb; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n\n";
exit;

#------------------
sub explore_contig(){
    my $c_n_init = shift;
	my $c_n_term = my $n_stroll=0;	
	my @jobs;		   # job-queue
	@{$jobs[0]} = split(/\,/, $node{$initial_nodes[$c_n_init]});
	while($#{$jobs[0]}>-1){ # while any job remains
		my $pointer;
		my $c_status = "M"; # M: over the maximum stroll trials, T: arrived at terminal node, E: arrived at dead-end node
		my $c_contig = "";
		my $c_node = $initial_nodes[$c_n_init];

		for($pointer=0; $pointer<$max_stroll; $pointer++){

  			# move to the next node and add the next edges to job-queue;
			$c_node .= $jobs[$pointer][0];

			for(my $i=0; $i<=$#terminal_nodes; $i++){
				if ($c_node =~ /$terminal_nodes[$i]/){
					$c_n_term = $i;
					$c_status="T";
					last;
				}
			}

			substr($c_node, 0, $step_size) = "";
			if ($c_status eq "T"){last;}

			if ($#{$jobs[$pointer+1]}==-1){
				if (!exists($node{$c_node})){   
					$c_status = "E";					# Dead-end node
					last;
				}
				push @{$jobs[$pointer+1]}, split(/\,/, $node{$c_node});
			}	
		}

		# Construct contigs and track the minimun depth of edges
		$c_node = $initial_nodes[$c_n_init];
		my $min_depth = 9999999;
		my $idx_edge = -1;
		for(my $i=0; $i<=$pointer; $i++){
			my @temp_edges = split(/\,/, $node{$c_node});
			my @temp_DP = split(/\,/, $e_DP{$c_node});
			for(my $j=0; $j<=$#temp_edges; $j++){
				if ($temp_edges[$j] eq $jobs[$i][0]){
					$idx_edge = $j;
					last;
				}
			}
			if ($idx_edge == -1){die "Error\n";}
			my $cur_edge_DP = $temp_DP[$idx_edge];
			if ($min_depth > $cur_edge_DP){$min_depth = $cur_edge_DP;}
			substr($c_node, 0, 3) = "";
			$c_node .= $jobs[$i][0];

			$c_contig .= $jobs[$i][0];
		}

		# Check terminal nodes
		if ($c_status eq "T") {
			($c_contig, ) = split /$terminal_nodes[$c_n_term]/, $c_contig;
		}

		if ($c_status eq "T") {
			if (!exists($contigs{$c_contig}{'num'})){
				$contigs{$c_contig}{'num'} = $min_depth;
				$contigs{$c_contig}{'status'} = $c_status;
				$contigs{$c_contig}{'n_init'} = $c_n_init if ($contigs{$c_contig}{'n_init'} eq "");
				$contigs{$c_contig}{'n_term'} = $c_n_term if (($contigs{$c_contig}{'n_term'} eq "")||($contigs{$c_contig}{'n_term'}>$c_n_term));
			}
		}

		# Erase used nodes from job-queue;
		for (my $e=$pointer; $e>=0; $e--){
			if ($#{$jobs[$e+1]}==-1) {
				shift @{$jobs[$e]};
			} else{	
				last;
			}
		}
	}
	return;
}

sub LoadSeeds{
	open(IN, "<$in_seeds") or die "[Error:explore_graph] Can't open $in_seeds.\n";
	while(my $head=<IN>){
		die "[Error:explore_graph] Seed file must be fasta format file.\n" if (substr($head, 0, 1) ne ">");
		my $str = <IN>;	
		chop($str);
		if ($head =~ "initial"){
			push @initial_nodes, $str;
		}elsif ($head =~ "terminal"){
			push @terminal_nodes, $str;
		}else{
			die "[Error:explore_graph] header ($str) must be initial or terminal.\n";
		}
	}
	close(IN);
	print "[Report:explore_graph] ", ($#initial_nodes+1), " initial nodes and ", ($#terminal_nodes+1), " terminal_nodes were loaded.\n";
}

sub LoadGraph{
	open(IN, "<$in_fname") or die "[Error:explore_graph] Can't open $in_fname.\n";
	open(OUT, ">$out_fname");
	while(<IN>){ 
		if (substr($_, 0, 1) eq "#"){
			print OUT $_;
			next;
		}
		my ($cur_node, $nodeDP, $num_edge, $cur_edge, $edgeDP,) = split /\s+/, $_;
		$node{$cur_node}=$cur_edge;
		$e_DP{$cur_node}=$edgeDP;
		$num_node++;
		if ($num_node % 5000000 == 0) {
			$t_end = new Benchmark;
	    	print "[Process:explore_graph] $num_node nodes were loaded; Mem. Used = ", JigsawSeq::memcheck(), " Gb; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n";
		}
	}
	close(IN);
	$t_end = new Benchmark;
	print "[Report:explore_graph] $num_node nodes were loaded; Mem. Used = ", JigsawSeq::memcheck(), " Gb; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n";
}

sub ReadArgument{
	my @command_line = @ARGV;               #command line argument

	#set default values;
	$k_mer_len = 120;
	$step_size = 3;

	#parse arguments
    while (defined(my $arg=shift(@ARGV))){
    	if (($arg eq "-I")||($arg eq "--input")){ $in_fname = shift @ARGV; next; }
    	if (($arg eq "-S")||($arg eq "--seeds")){ $in_seeds = shift @ARGV; next; }
    	if (($arg eq "-O")||($arg eq "--output")){ $out_fname = shift @ARGV; next; }
    	if (($arg eq "-L")||($arg eq "--length")){ $exp_contig_size = shift @ARGV; next; }
    	if (($arg eq "-k")||($arg eq "--kmer")){ $k_mer_len = shift @ARGV; next; }
    	if (($arg eq "-s")||($arg eq "--step")){ $step_size = shift @ARGV; next; }
       	if (($arg eq "-h")||($arg eq "-?")||($arg eq "--help")){ printError(); next; }
    	PrintError("Can't identify the argument: $arg");
    }

    #check validity of arguments
    PrintError("Arguments -I, -S, -O, -L are required.") unless ( (defined $in_fname)&&(defined $in_seeds)&&(defined $exp_contig_size)&&(defined $out_fname) );
	PrintError("Length of k-mer (--kmer) must be 40 <= k-mer <=150.") unless ((40<=$k_mer_len)&&($k_mer_len<=150));
	PrintError("Step size (--step) must be 1, 2, or 3.") unless ((1<=$step_size)&&($step_size<=3));
	PrintError("Length of k-mer (--kmer) must be divisible by step size (--step).") if ($k_mer_len % $step_size);
	PrintError("Input filename ($in_fname) does not exist.") unless (-e $in_fname);
	PrintError("Input filename ($in_seeds) does not exist.") unless (-e $in_seeds);
}

sub PrintError{
	print "\n[Error:explore_graph] ", join("\n", @_), "\n\n";

	print "DESCRIPTION:\n    explore_graph.pl is to search contigs from the de bruijn graph with seeds.\n\n";
	print "USAGE: explore_graph.pl -I [input:graph] -S [input:seeds] -L [input:expected length] -O [output:prefix]\n\n";
	print "OPTIONS:\n",
		"    -I, --input       FILE    Input file name of graph [Required]\n",
		"    -S, --seeds       FILE    Seeds file name (fasta format) [Required]\n",
		"    -O, --output      STR     Prefix of output file [Required]\n",
		"    -L, --length      INT     Expected length of contig [Required]\n",
		"    -k, --kmer        INT     Length of k-mer [Default: 120]\n",
		"    -s, --step        INT     Step size for exploring the graph [Default: 3]\n";
#		"    -h, -?, --help       This help message\n",
	print "\n";
	exit -1;
}
