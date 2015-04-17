#!/usr/bin/perl

# This script is to merge de Bruijn graph files with index files
# developed by Jung-Ki Yoon
# last modified: Apr-12-2015

use strict;
use JigsawSeq;
use Benchmark ':hireswallclock';

my $usage = "[Usage] ./merge_graph.pl [input: graph] [input:graph] [output: graph]\n";
die $usage unless ($#ARGV == 2);
my ($in_fname1, $in_fname2, $out_fname,) = @ARGV;
my $index_fname1 = $in_fname1 . ".index";
my $index_fname2 = $in_fname2 . ".index";

# If there is no indexed files, run graph_indexing.pl first
if (!(-e $index_fname1)){
	die "[Error:merge_graph] Can't find $in_fname1.\n" unless (-e $in_fname1);
	system("./graph_indexing.pl -I $in_fname1 -S");
}
if (!(-e $index_fname2)){
	die "[Error:merge_graph] Can't find $in_fname2.\n" unless (-e $in_fname1);
	system("./graph_indexing.pl -I $in_fname2 -S");
}

# Load index files
my %IDX;
open(IN, "<$index_fname1");
my ($temp, $word_len, ) = split /\s+/, <IN>;
if ($temp ne "#WORD_LEN") {die "[Error:merge_graph] $index_fname1 is not index file.\n";}
while(<IN>){
	my ($word, $line,) = split/\s+/, $_;
	$IDX{$word}{'1'} = $line;
}
close(IN);
open(IN, "<$index_fname2");
my ($temp, $word_len2, ) = split /\s+/, <IN>;
if ($temp ne "#WORD_LEN") {die "[Error:merge_graph] $index_fname2 is not index file.\n";}
if ($word_len != $word_len2) {die "[Error:merge_graph] The sizes of index letters must be same. ($word_len <> $word_len2)\n";}
while(<IN>){
	my ($word, $line,) = split/\s+/, $_;
	$IDX{$word}{'2'} = $line;
}
close(IN);

# Merge
print "[Report:merge_graph] Merge $in_fname1 and $in_fname2 into $out_fname\n";
our $t_begin = new Benchmark;
our $t_end;
my $num_kmer=my $num_merged=0;
open(IN1, "<$in_fname1") or die "[Error:merge_graph] Can't open $in_fname1.\n";
open(IN2, "<$in_fname2") or die "[Error:merge_graph] Can't open $in_fname2.\n";
open(OUT, ">$out_fname");
my $p_1=my $p_2=0;
foreach my $idx_str (sort (keys %IDX)){
	my $load_nodes_1=my $load_nodes_2=0;
	my %data;
	for(;$p_1<=$IDX{$idx_str}{'1'};$p_1++){
		my @foo = split /\s+/, <IN1>;
	    $data{$foo[0]} = join("\t", @foo[1..4]);
	    $load_nodes_1++;
	    $num_kmer++;
	}
#	print "[Report:merge_graph] From $in_fname1, $load_nodes_1 nodes ($idx_str) were loaded.\n";
	for(;$p_2<=$IDX{$idx_str}{'2'};$p_2++){
    	my ($cur_kmer, @info2) = split /\s+/, <IN2>;
    	$load_nodes_2++;
		if (exists $data{$cur_kmer}) { # merge
			my %BackNode;
			my @info1 = split /\t/, $data{$cur_kmer};
			my @oriBackL = split /\,/, $info1[2];
			my @oriBackD = split /\,/, $info1[3];
			for(my $i=0; $i<=$#oriBackL; $i++){ $BackNode{$oriBackL[$i]} = $oriBackD[$i]; }

			my @curBackL = split /\,/, $info2[2];
			my @curBackD = split /\,/, $info2[3];
			for(my $i=0; $i<=$#curBackL; $i++){ $BackNode{$curBackL[$i]} += $curBackD[$i]; }

			$info1[2] = $info1[3] = "";
			$info1[1] = 0;
			foreach my $k (keys %BackNode){
		    	$info1[1]++;
		    	$info1[2] .= ($k . ",");
	    		$info1[3] .= ($BackNode{$k} . ",");
			}
			$info1[0] += $info2[0];

			$data{$cur_kmer} = join("\t", @info1[0..3]);
			$num_merged++;
    	}else{
			$data{$cur_kmer} = join("\t", @info2[0..3]);
			$num_kmer++;
	    }
	}
#	print "[Report:merge_graph] From $in_fname2, $load_nodes_2 nodes ($idx_str) were loaded.\n";

	foreach my $k (keys %data){
    	print OUT join("\t", $k, $data{$k}), "\n";
	}
#	$t_end = new Benchmark; 
#	print "[Report:merge_graph] Total $num_kmer nodes has beed processed, and $num_merged nodes has been merged.\n", 
#	      "                     Mem. Used = ", JigsawSeq::memcheck(), " Gb; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n";
}
close(IN1);
close(IN2);
close(OUT);

$t_end = new Benchmark; 
print "[Report:merge_graph] In $out_fname, $num_kmer nodes were recorded.\n", 
      "                     Mem. Used = ", JigsawSeq::memcheck(), " Gb; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n\n";
exit;
