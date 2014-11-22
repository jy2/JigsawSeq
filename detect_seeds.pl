#!/usr/bin/perl

# Detect initial/terminal seeds.
# Requirement: bwa 0.7.5a-r405
# last modified: Nov-22-2014
# Developed by Jung-Ki Yoon

use strict;
use JigsawSeq;
use Benchmark ':hireswallclock';

#my $usage = "[Usage] ./detect_seeds.pl [input: Kmer fasta] [input: vector fasta] [k: k-mer length] [s: step size] [r: cutoff ratio] [output: seeds]\n";
my $usage = "[Usage] ./detect_seeds.pl [input: graph] [input: vector fasta] [k: k-mer length] [s: step size] [r: cutoff ratio of seeds] [t: number of threads] [output: seeds]\n";
die $usage unless ($#ARGV == 6);
our ($in_fname, $vector_fname, $in_kmer, $step_size, $cutoff_ratio, $threads, $out_fname,) = @ARGV;
our $kmer_fname = "TMP_$in_fname\.kmer.fa";
my $len_ref_seed = int($in_kmer*1.5);				# allow a half size of k-mer as insertion on seeds
our $ref_seed_fname = "TMP_k$in_kmer\_" . $vector_fname;
our $seeds_sam_fname = "TMP_$in_fname\.seed.sam";
my $t_begin = new Benchmark;
print "[Report:detect_seeds] input: $in_fname vector: $vector_fname k-mer: $in_kmer, step_size: $step_size, cut_seeds: $cutoff_ratio, output: $out_fname\n";

Kmer2Fa();
Vector2Fa();

print "[Report:detect_seeds] Align nodes to vector using bwa\nBWA----------------------\n";
system("./bwa index $ref_seed_fname > $ref_seed_fname\.log");
system("./bwa mem -t $threads -O2 -E1 $ref_seed_fname $kmer_fname > $seeds_sam_fname");
my $t_end = new Benchmark;
print "-----------------------BWA\n[Report:detect_seeds] Alignment was completed; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n";

AnalyzingSeeds();
exit;

#----------------------------
sub Kmer2Fa{
	open(IN, "<$in_fname") or die "[Error:detect_seeds] Can't open $in_fname.\n";
	open(OUT, ">$kmer_fname");
	my $num_lines=0;
	while(<IN>){ 
		if (substr($_, 0, 1) eq "#"){ next;	}
		my ($str, $DP,) = split /\s+/, $_;
		print OUT ">$num_lines\_$DP\n";
		print OUT "$str\n";
		$num_lines++;
	}
	close(IN);
	close(OUT);
	print "[Report:detect_seeds] $num_lines nodes were ready for running dectect seeds\n";
}

sub Vector2Fa{
	open(IN, "<$vector_fname") or die "[Error:detect_seeds] Can't open $vector_fname.\n";
	open(OUT, ">$ref_seed_fname");
	<IN>;
	my ($seq,)=split /\s+/, <IN>;
	print OUT ">initial\n", uc(substr($seq,-$len_ref_seed)), "\n";
	print OUT ">terminal\n", uc(substr($seq,0,$len_ref_seed)), "\n";
	close(OUT);
	close(IN);
}

sub AnalyzingSeeds {
	my $num_lines=my $num_align_init=my $num_init=my $num_align_term=my $num_term=0;
	my $max_DP_init=my $max_DP_term=-1;
	open(IN, "<$seeds_sam_fname") or die "[Error:detect_seeds] Can't open $seeds_sam_fname.\n";
	open(OUT, ">$out_fname");
	while(<IN>){ 
		next if (substr($_, 0, 1) eq "@");
		$num_lines++;
		if ($num_lines % 1000000 == 0) { print "[Process:detect_seeds] $num_lines lines were processed.\n";}
		my ($qname, $flag, $rname, $pos, $mapQ, $CIGAR, $rnext, $pnext, $tlen, $str, $qual, $NM, $MD, $AS, $XS,) = split /\s+/, $_;
		next if ($flag != 0);
		my @d = split /[MIDNSHPX=]/, $CIGAR;
		my @l = split /[0-9]+/, $CIGAR;
		shift @l;
		my $sum_d=0;
		if ($rname eq "initial"){
			$num_align_init++;
			($qname, my $DP) = split /_/, $qname;
			if ($max_DP_init == -1){$max_DP_init=$DP;}
			next if (($max_DP_init / $DP > $cutoff_ratio)||($DP < $cutoff_ratio));
			for(my $i=0; $i<=$#d; $i++){
				$sum_d+=$d[$i];
			}
			if (($pos + $sum_d) == ($len_ref_seed+1) ) {
		    	if ($l[$#l] eq "S"){
				next;
			    }
		    	print OUT ">initial $DP\n$str\n";
			    $num_init++;
			}
		}elsif ($rname eq "terminal"){
			$num_align_term++;
			($qname, my $DP) = split /_/, $qname;
			if ($max_DP_term == -1){$max_DP_term=$DP;}
			next if (($max_DP_term / $DP > $cutoff_ratio)||($DP < $cutoff_ratio));
			for(my $i=0; $i<=$#d; $i++){
				if (($l[$i] eq "M")||($l[$i] eq "I")){
					$sum_d+=$d[$i];
				}
			}
			if (($pos == 1)&&($sum_d==($in_kmer-$step_size))){
				print OUT ">terminal $DP\n$str\n";
				$num_term++;
			}
		}
	}
	close(IN);
	close(OUT);

	my $t_end = new Benchmark; 
	print "[Report:detect_seeds] Cutoff ratio of seeds: $cutoff_ratio\tMax_inti_DP: $max_DP_init\tMax_term_DP: $max_DP_term\n";
	print "[Report:detect_seeds] $num_lines K-mers were processed.\n";
	print "[Report:detect_seeds] $num_align_init K-mers were alinged to initial node.\t$num_init were considered as initial nodes.\n";
	print "[Report:detect_seeds] $num_align_term K-mers were aligned to terminal node.\t$num_term were considered as terminal nodes.\n";
	print "                      Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n\n";
}