#!/usr/bin/perl

# Mapping pair-end read (Jigsaw-Seq data) into candidate contigs. 
# Requirement: bwa 0.7.5a-r405     samtools 0.1.19-44428cd
# last modified: Nov-22-2014
# Developed by Jung-Ki Yoon

use strict;
use JigsawSeq;
use Benchmark ':hireswallclock';

my $usage = "[Usage] ./mapping_contigs.pl [input: candidates fasta] [input: fastq format (forward)] [input: fastq format (reverse)] [k: k-mer length] [t: num threads] [output: prefix of sam/bam]\n";
die $usage unless ($#ARGV == 5);
our ($candi_fa, $in_fnameF, $in_fnameR, $in_kmer, $threads, $out_fname,) = @ARGV;
our $t_begin = new Benchmark;
our $t_end;

print "[Report:mapping_contigs] Align raw reads to candidate contigs using bwa\nBWA----------------------\n";
system("./bwa index $candi_fa");
system("./bwa mem -t $threads -B100 -O100 -E50 $candi_fa $in_fnameF $in_fnameR > TMP_$out_fname\.sam");
$t_end = new Benchmark;
print "-----------------------BWA\n",
      "[Report:mapping_contigs] Alignment was completed; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n";


my $num_EM = ExtractExactMatch();
print "[Report:mapping_contigs] $num_EM reads were aligned to contigs exactly.\n";


print "[Report:mapping_contigs] Calculate depths for each candidates using samtools\nSAMTOOLS----------------------\n";
system("./samtools view -bS TMP_$out_fname\.exact.sam -o TMP_$out_fname\.exact.unsorted.bam");
system("./samtools sort TMP_$out_fname\.exact.unsorted.bam TMP_$out_fname\.exact.sorted");
system("./samtools index TMP_$out_fname\.exact.sorted.bam");
system("./samtools depth TMP_$out_fname\.exact.sorted.bam > $out_fname\.DP");
$t_end = new Benchmark;
print "-----------------------SAMTOOLS\n", 
	  "[Process:mapping_contigs] Depth distributions were calculated; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n\n";


exit;

#---------------------------

sub ExtractExactMatch{
	open(IN, "<TMP_$out_fname\.sam") or die "[Error] Can't open TMP_$out_fname\.sam.\n";
	open(OUT, ">TMP_$out_fname\.exact.sam");
	my $num_lines=my $num_exact=0;
	while(<IN>){ 
		if (substr($_, 0, 1) eq "@"){ print OUT $_; next; }
		$num_lines++;
		unless ($num_lines % 5000000){
			$t_end = new Benchmark;
			print "[Process:mapping_contigs] $num_lines lines were processed; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n";
		}

		my ($qname, $flag, $rname, $pos, $mapQ, $CIGAR, $rnext, $pnext, $tlen, $str, $qual, $NM, $MD, $AS, $XS,) = split /\s+/, $_;
		my @d = split /[MIDNSHPX=]/, $CIGAR; # number in CIGAR
		my @l = split /[0-9]+/, $CIGAR;      # letter in CIGAR
		shift @l;
		substr($MD, 0, 5) = ""; 			 #MD:Z:150 -> 150

		# ignore the fusion reads
		next unless (($rnext eq "=")||($rnext eq "*"));

		# use exact match in contigs.
		if (($#l==0)&&($l[0] eq "M")&&($MD eq $d[0])) {
			print OUT $_;
			$num_exact++;
		}

		# non-unique reads (XS = AS)
	}
	close(IN);
	close(OUT);
	return $num_exact;
}
