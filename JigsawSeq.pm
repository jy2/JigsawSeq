#!/usr/bin/perl

package JigsawSeq;

sub rev_comp($){
	my @in_str = split //, uc(shift);
	my %complemBase = ('A'=>'T', 'T'=>'A', 'C'=>'G', 'G'=>'C');
	my $re_str;
	for(my $i=0; $i<=$#in_str; $i++){
		return "Error" unless (exists $complemBase{$in_str[$i]});
		$re_str = $complemBase{$in_str[$i]} . $re_str;
	}
	return $re_str;
}

sub codon2aa($){
	my $codon = uc(shift);
	my %genetic_codes = ('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');;
	return "Error"	unless (exists $genetic_codes{$codon});
	return $genetic_codes{$codon};
}

sub memcheck{
	my $memory_used = `ps h -o size $$`;			#my $memory_used = `ps h -o sz $$`;
	chop($memory_used);
	$memory_used /= (1024*1024);
	$memory_used = substr($memory_used, 0, 4);
	return $memory_used;
}

1;