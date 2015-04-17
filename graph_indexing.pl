#!/usr/bin/perl

# Build index file for graph
# last modified: Apr-13-2015
# Developed by Jung-Ki Yoon

use strict;
use JigsawSeq;
use Benchmark ':hireswallclock';

our ($in_fname, $out_fname, );
our $word_len;
our $sort;

my %INDEX;
my $num_lines=0;

# Parse arguments
ReadArgument();

my $t_begin = new Benchmark;
print "[Report:graph_indexing] input: $in_fname output: $out_fname\n";

# Sorting
if ($sort == 0){
    my $sort_fname = $in_fname . ".tmp";
    print "[Report:graph_indexing] Sorting\n";
    `sort $in_fname > $sort_fname`;
    `mv $sort_fname $in_fname`;
}

# Indexing
print "[Report:graph_indexing] Indexing\n";
open(IN, "<$in_fname") or die "Can't open $in_fname\n";
while(<IN>){
    my ($cur_node, $nodeDP, $num_edge, $cur_edge, $edgeDP,) = split /\s+/, $_;
    my $cur_idx = substr($cur_node, 0, $word_len);
    $INDEX{$cur_idx} = $num_lines;
    $num_lines++;
}
close(IN);

# Recording
open(OUT, ">$out_fname");
print OUT "#WORD_LEN $word_len\n";
foreach my $key (sort(keys %INDEX)){
    print OUT join("\t", $key, $INDEX{$key}), "\n";
}
close(OUT);

my $t_end = new Benchmark;
print "[Report:graph_indexing] Indexing was completed; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n";
exit;

#----------------------------
sub ReadArgument{
    my @command_line = @ARGV;               #command line argument

    #set default values;
    $sort = 0;
    $word_len = 4;

    #parse arguments
    while (defined(my $arg=shift(@ARGV))){
        if (($arg eq "-I")||($arg eq "--input")){ $in_fname = shift @ARGV; next; }
        if (($arg eq "-O")||($arg eq "--output")){ $out_fname = shift @ARGV; next; }
        if (($arg eq "-W")||($arg eq "--word")){ $word_len = shift @ARGV; next; }
        if (($arg eq "-S")||($arg eq "--sorted")){ $sort = 1; next; }
    }

    PrintError("Arguments -I is required.") unless (defined $in_fname);
    if (!(defined $out_fname)){ $out_fname = $in_fname . ".index"; }
}

sub PrintError{
    print "\n[Error:graph_indexing] ", join("\n", @_), "\n\n";
    print "USAGE: graph_indexing.pl \n";
    print "OPTIONS:\n",
        "    -I, --input     FILE    Input file name (graph file) [Required]\n",
        "    -O, --output    FILE    Prefix of index filename [Default: input file name]\n",
        "    -W, --word      INT     Length of index letters [Default: $word_len]\n",
        "    -S, --sorted            Input file is already alphabetically sorted. Skip the sorting step.\n"; 
    print "\n";
    exit -1;
}
