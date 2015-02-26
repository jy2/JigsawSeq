#!/usr/bin/perl

# Clean up the de Bruijn graph to speed up the searching step.
# last modified: Nov-22-2014
# Developed by Jung-Ki Yoon

use strict;
use JigsawSeq;
use Benchmark ':hireswallclock';

my $usage = "[Usage] ./cleanup_graph.pl [input: graph] [m: min.depth of node/edge] [r: cutoff ratio of edge] [output: graph]\n";
die $usage unless ($#ARGV == 3);
my ($in_fname, $minNodeDepth, $cutoff_ratio, $out_fname,) = @ARGV;
my $minEdgeDepth = $minNodeDepth;
#my $polyA = "AAAAAAAA";   # due to illumina seq with short insertion size sequences
my $strange1 = "CAGTAATACAAGGGGTGTTGCGGAGTGTATACT";    
my $strange2 =   "GTAATACAAGGGGTGTTACATAAACAGTAATACAA";    # due to Gibbs cloning breakpoint

print "[Report:cleanup_graph] input: $in_fname min_depth: $minNodeDepth cut_edge: $cutoff_ratio output:$out_fname\n";

my $t_begin = new Benchmark;
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
    if (($seq =~ $strange1)||($seq =~ $strange2)){
        $num_Gibbs_skipped_node++;
        next;
    }

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
        if ( (($max/$DPs[$i])<=$cutoff_ratio) && ($DPs[$i]>=$minEdgeDepth) ){
            push @new_DP_edges, ($DPs[$i] . " " . $edges[$i]);
        }else{
            $num_skipped_edge++;
        }
    }
#    print "Cleaned\n";
#    print join("|", @new_DP_edges), " $#new_DP_edges\n";
    @new_DP_edges = sort { &csort($b, $a) } @new_DP_edges;
    $NodeDepth=0;
    for (my $i=0; $i<=$#new_DP_edges; $i++){
        my ($c_DP, $c_edge) = split /\s/, $new_DP_edges[$i];
        $new_edge_str .= ($c_edge . ",");
        $new_DP_str .= ($c_DP . ",");
        $NodeDepth += $c_DP;
    }
#    print "Sorted\n";
#    print join("|", @new_DP_edges), " $#new_DP_edges\n";
    if ($#new_DP_edges >= 0){
#        print join("\t", $seq, $NodeDepth, ($#new_DP_edges+1), $new_edge_str, $new_DP_str), "\n";
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
print "[Report:cleanup_graph] $num_Gibbs_skipped_node nodes were neglected since the vector node was problematic sequences ($strange1, $strange2).\n";
print "[Report:cleanup_graph] $num_skipped_edge edges were cleaned up (maxEdgeDepth/EdgeDepth > $cutoff_ratio or EdgeDepth is less than $minEdgeDepth.)\n";
print "[Report:cleanup_graph] Remaining $num_record_nodes nodes and $num_record_edges edges were recorded.\n",
      "                       Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n\n";

$t_begin = new Benchmark;
print "[Report:sorting] Sorting initiated.\n";
`sort -r -n -k 2 $out_fname > $out_fname\.sort`;
$t_end = new Benchmark;
print "[Report:sorting] Sorting was completed.\n",
      "                 Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n\n";
exit;

sub csort{
    my $ta = shift;
    my $tb = shift;
    (my $ca, ) = split /\s/, $ta;
    (my $cb, ) = split /\s/, $tb;
    $ca <=> $cb;
}
