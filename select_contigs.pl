#!/usr/bin/perl
use strict;

# Search contigs from candidates using depth distribution
# requirement: perl version 5.10.1 or higher
# last modified: Nov-22-2014
# Developed by Jung-Ki Yoon

my $usage = "./select_contigs.pl [input: candidate.fa] [input: samtools depth] [input: k-mer length] [input: step-size] [input: cutoff of CV] [input: read length] [output_prefix: contigs]\n";
die $usage if ($#ARGV != 6);
our ($candi_fa, $in_DP, $in_kmer, $step_size, $criteria, $read_length, $out_fname,) = @ARGV;
my $clipped = 25;

$read_length -= $clipped;
my $candi_stat = "TMP_$in_DP\.stat";
my $contig_stat = "$out_fname\.stat";
my $contig_fa = "$out_fname\.fa";
our %seq;
our %len;

LoadCandidates();
CalculateStat();

# Select contigs from candidates
open(IN, "<$candi_stat") or die "[Error:select_contigs] Can't open $candi_stat.\n";
open(OUT, ">$contig_stat");
open(OUTFA, ">$contig_fa");
my $num_process=my $num_pass=0;
my $sum_ave=0;
while(<IN>){
    if (substr($_, 0, 1) eq "#"){ print OUT $_; next; }
    my ($name, $len, $ave, $std, $cv,) = split /\s+/, $_;
    if ($cv <= $criteria){
    	print OUT join("\t", $name, $len-($in_kmer-$step_size)*2, $ave, $std, $cv), "\n";
        print OUTFA ">$name\t$ave\t$cv\n$seq{$name}\n";
    	$sum_ave += $ave;
        $num_pass++;
    }
    $num_process++;
}
close(IN);
close(OUT);
close(OUTFA);

print "[Report:select_contigs] $num_process contigs were processed and $num_pass contigs were passed (cutoff: CV <= $criteria); Average depth of passed contigs = ", ($sum_ave / $num_pass), "\n\n";
exit;

#-----------------------
sub LoadCandidates{
    open(IN, "<$candi_fa") or die "[Error:select_contigs] Can't open $candi_fa.\n";
    while(my $header = <IN>){
        die "[Error:select_contigs] candidate.fa must be fasta format.\n$_\n" unless (substr($header, 0, 1) eq ">");
        substr($header, 0, 1)  = "";
        chop($header);
        (my $cur_seq,) = split /\s+/, <IN>;

        $len{$header} = length($cur_seq);
        substr($cur_seq, 0, ($in_kmer - $step_size))=""; #trim seeds
        substr($cur_seq, -($in_kmer-$step_size)) = "";
        $seq{$header} = $cur_seq;
    }
    close(IN);
}

sub CalculateStat{
    open(IN, "<$in_DP") or die "[Error:select_contigs] Can't open $in_DP.\n";
    open(OUT, ">$candi_stat");
    print OUT "#Regions: $read_length ~ (length of candidates)-$read_length\n#candidate_name\tcandidate_length\taveDP\tstdDP\tCV\n";
    my $prev_candi="";
    my @curDPs;

    while(<IN>){
        my ($cur_candi, $pos, $depth) = split /\s+/, $_;
        if ($cur_candi ne $prev_candi){
            if ($prev_candi ne ""){
                for (my $i=1; $i<=$len{$prev_candi}; $i++){ $curDPs[$i] = 0 unless (defined $curDPs[$i]); }
                print OUT join("\t", $prev_candi, $len{$prev_candi}, calc_stat($read_length, $len{$prev_candi}-$read_length, @curDPs)), "\n";
                @curDPs = qw();
            }
            $prev_candi = $cur_candi;
        }

        $curDPs[$pos] = $depth;
    }
    close(IN);

    for (my $i=1; $i<=$len{$prev_candi}; $i++){ $curDPs[$i] = 0 unless (defined $curDPs[$i]); }
    print OUT join("\t", $prev_candi, $len{$prev_candi}, calc_stat($read_length, $len{$prev_candi}-$read_length, @curDPs)), "\n";
    close(OUT);
}

sub calc_stat(){
    my $init = shift @_;
    my $term = shift @_;
    my $len_ = $term - $init + 1;
    die "[Error:select_contigs] Failed to calc_stat(); init ($init) >= term ($term)\n" if ($len_ <=0);
    my @D = @_;

    my $c_sum=0;
    for (my $i=$init; $i<=$term; $i++){
        $c_sum += $D[$i];
    }
    my $c_ave = $c_sum / $len_;

    my $c_std=0;
    for (my $i=$init; $i<=$term; $i++){
        $c_std += ($D[$i]-$c_ave)*($D[$i]-$c_ave);
    }
    $c_std = sqrt($c_std/$len_);

    my $CV = $c_std/$c_ave;

    return ($c_ave, $c_std, $CV);
}
