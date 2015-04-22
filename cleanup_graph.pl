#!/usr/bin/perl

# Clean up the de Bruijn graph to speed up the searching step and reduce false postive.
# last modified: Apr-22-2015
# Developed by Jung-Ki Yoon

use strict;
use JigsawSeq;
use Benchmark ':hireswallclock';

our ($in_fname, $exc_fname, $min_depth, $cut_edge, $out_fname,);
#my $polyA = "AAAAAAAA";   # due to illumina seq with short insertion size sequences


ReadArgument();
my $minEdgeDepth = my $minNodeDepth = $min_depth;

my $t_begin = new Benchmark;
print "[Report:cleanup_graph] input: $in_fname min_depth: $minNodeDepth cut_edge: $cut_edge output:$out_fname\n";

# Load exclusion.fa for filtering problematic reads
my @stranges;
open(IN, "<$exc_fname") or die "[Error:cleanup_graph] Can't open $exc_fname.\n";
while(<IN>){
    if (substr($_, 0, 1) ne ">"){ die "[Error:cleanup_graph] $exc_fname is not a fasta file.\n $_";}
    my ($cur_str, ) = split/\s+/, <IN>;
    push(@stranges, $cur_str);
}
close(IN);
print "[Report:cleanup_graph] ", ($#stranges+1), " sequences were loaded as filtering seq.\n";

# Clean-up process
open(IN, "<$in_fname") or die "[Error:cleanup_graph] Can't open $in_fname.\n";
open(OUT, ">$out_fname");
my @edges;my @DPs;
my $num_skipped_edge=my $num_DP_skipped_node=my $num_polyA_skipped_node=my $num_Gibbs_skipped_node=0;
my $num_record_nodes=my $num_record_edges=0;
while(<IN>){
    if (substr($_, 0, 1) eq "#") { print OUT $_; next; }
    chop($_);
    (my $seq, my $NodeDepth, my $num_edge, my $edge_str, my $DP_str,) = split(/\t/, $_);
    if ($NodeDepth<$minNodeDepth){
        $num_DP_skipped_node++;
        next;
    }

#    if (($seq =~ $polyA)){
#        $num_polyA_skipped_node++;
#        next;
#    }

    my $idx=0;
    for(my $i=0; $i<=$#stranges; $i++){
        if ($seq =~ $stranges[$i]){
            $num_Gibbs_skipped_node++;
            $idx=1;
            last;
        }
    }
    if ($idx==1){next;}

    @edges = split /\,/, $edge_str;
    @DPs = split /\,/, $DP_str;

    my $max=-1;
    for (my $i=0; $i<=$#DPs; $i++){
        $max = $DPs[$i] if ($max < $DPs[$i]);
    }

    my @new_DP_edges;
    my $new_edge_str = "";
    my $new_DP_str = "";

    for (my $i=0; $i<=$#DPs; $i++){
        if ( (($max/$DPs[$i])<=$cut_edge) && ($DPs[$i]>=$minEdgeDepth) ){
            push @new_DP_edges, ($DPs[$i] . " " . $edges[$i]);
        }else{
            $num_skipped_edge++;
        }
    }

    @new_DP_edges = sort { &csort($b, $a) } @new_DP_edges;
    $NodeDepth=0;
    for (my $i=0; $i<=$#new_DP_edges; $i++){
        my ($c_DP, $c_edge) = split /\s/, $new_DP_edges[$i];
        $new_edge_str .= ($c_edge . ",");
        $new_DP_str .= ($c_DP . ",");
        $NodeDepth += $c_DP;
    }

    if ($#new_DP_edges >= 0){
        print OUT join("\t", $seq, $NodeDepth, ($#new_DP_edges+1), $new_edge_str, $new_DP_str), "\n";
        $num_record_nodes++;
        $num_record_edges += ($#new_DP_edges+1);
    }
}
close(OUT);
close(IN);

my $num_skipped_node = $num_DP_skipped_node + $num_polyA_skipped_node + $num_Gibbs_skipped_node;
my $t_end = new Benchmark;
print "[Report:cleanup_graph] $num_DP_skipped_node nodes were neglected since the depth was less than $minNodeDepth.\n";
#print "[Report:cleanup_graph] Poly-A node ($polyA): $num_polyA_skipped_node.\n";
print "[Report:cleanup_graph] $num_Gibbs_skipped_node nodes were neglected since the vector node had at least one of ", ($#stranges+1), " problematic sequences.\n";
print "[Report:cleanup_graph] $num_skipped_edge edges were cleaned up (maxEdgeDepth/EdgeDepth > $cut_edge or EdgeDepth is less than $minEdgeDepth.)\n";
print "[Report:cleanup_graph] Remaining $num_record_nodes nodes and $num_record_edges edges were recorded.\n",
      "                       Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n\n";

exit;

#----------------------------
sub ReadArgument{
    my @command_line = @ARGV;               #command line argument

    #set default values;
    $min_depth = 2;
    $cut_edge = 150;
    $exc_fname = "exclusion.fa";

    #parse arguments
    while (defined(my $arg=shift(@ARGV))){
        if (($arg eq "-I")||($arg eq "--input")){ $in_fname = shift @ARGV; next; }
        if (($arg eq "-O")||($arg eq "--output")){ $out_fname = shift @ARGV; next; }
        if (($arg eq "-m")||($arg eq "--min_depth")){ $min_depth = shift @ARGV; next; }
        if (($arg eq "-e")||($arg eq "--exclude")){ $exc_fname = shift @ARGV; next; }
        if ($arg eq "--cut_edge") { $cut_edge = shift @ARGV; next; }
        if (($arg eq "-h")||($arg eq "-?")||($arg eq "--help")){ printError(); next; }
        PrintError("Can't identify the argument: $arg");
    }

    #check validity of arguments
    PrintError("Arguments -I, -O are required.") unless ( (defined $in_fname)&&(defined $out_fname) );
    PrintError("Input filename ($in_fname) does not exist.") unless (-e $in_fname);
}

sub PrintError{
    print "\n[Error:cleanup_graph] ", join("\n", @_), "\n\n";

    print "DESCRIPTION:\n    cleanup_graph.pl is to clean up the graph for reducing false positive of contigs.\n\n";
    print "USAGE: cleanup_graph.pl -I [input:graph] -O [output:graph]\n\n";
    print "OPTIONS:\n",
        "    -I, --input       FILE    Input file name of graph [Required]\n",
        "    -O, --output      FILE    Output file name [Required]\n",
        "    -m, --min_depth   INT     Minimum depth of nodes and edges [Default: $min_depth]\n",
        "    -e, --exclude     FILE    List of exclusion sequences while constructing graph [Default: $exc_fname]\n",
        "    --cut_edge        INT     Cutoff ratio for edges [Default: $cut_edge]\n";
#       "    -h, -?, --help       This help message\n",
    print "\n";
    exit -1;
}

sub csort{
    my $ta = shift;
    my $tb = shift;
    (my $ca, ) = split /\s/, $ta;
    (my $cb, ) = split /\s/, $tb;
    $ca <=> $cb;
}

