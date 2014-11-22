#!/usr/bin/perl

# Search candidates of contigs from the de bruijn graph with seeds.
# last modified: Nov-22-2014
# Developed by Jung-Ki Yoon

use strict;
use JigsawSeq;
use Benchmark ':hireswallclock';

my $usage = "[Usage] ./explore_graph.pl [input: graph] [input: seeds fasta] [k: k-mer length] [s: step size] [L: expected length of contigs] [output: candidate fasta]\n";
our ($in_fname, $in_seeds, $in_kmer, $step_size, $exp_length, $out_fname,) = @ARGV;
die $usage unless ($#ARGV == 5);

our $t_begin = new Benchmark;
our $t_end;
our %node; 
our $a_node = \%node; 
our $num_node=0;
our @initial_nodes; 
our @terminal_nodes;
our $max_stroll = int(($exp_length + 50 + $in_kmer)/$step_size);  # allow the half size of expected contig length as insertion size.

print "[Report:explore_graph] input: $in_fname seeds: $in_seeds k-mer_len: $in_kmer step_size: $step_size expected_contig_size: $exp_length max_num_stroll: $max_stroll output: $out_fname\n";

LoadSeeds();
LoadGraph();

# Explore graph with each initial seed.
my %contigs;
for(my $i=0; $i<=$#initial_nodes; $i++){
	print "[Process:explore_graph] Searching with ",$i+1,"-th initial nodes: $initial_nodes[$i]\n";
	explore_contig($i);
}
close(OUT);

# Print out
my $num_contigs=0; 
open(OUT, ">$out_fname");
for my $key (keys %contigs){
	$num_contigs++;
	print OUT ">$in_fname\_$num_contigs\n", $initial_nodes[$contigs{$key}{'n_init'}], $key, $terminal_nodes[$contigs{$key}{'n_term'}], "\n";
}
close(OUT);

$t_end = new Benchmark;
print "[Report:explore_graph] $num_contigs contigs (T+M) were detected.\n", 
	  "                       Mem. Used = ", JigsawSeq::memcheck(), " Gb; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n\n";
exit;


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

		# Construct contigs
		for(my $i=0; $i<=$pointer; $i++){
			$c_contig .= $jobs[$i][0];
		}

		# Check terminal nodes
		if ($c_status eq "T") {
			($c_contig, ) = split /$terminal_nodes[$c_n_term]/, $c_contig;
		}

		if (($c_status eq "T")) {
#		if (($c_status eq "T")||($c_status eq "M")) {
			$contigs{$c_contig}{'num'}++;
			$contigs{$c_contig}{'status'} = $c_status;
			$contigs{$c_contig}{'n_init'} = $c_n_init if ($contigs{$c_contig}{'n_init'} eq "");
			$contigs{$c_contig}{'n_term'} = $c_n_term if (($contigs{$c_contig}{'n_term'} eq "")||($contigs{$c_contig}{'n_term'}>$c_n_term));
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