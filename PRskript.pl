#!usr/bin/perl

open(FH,"piRanhaoutput") or die print "can't open file: $!";
while(<FH>)
	{
	if($_=~/tRF-Glu-CTC/) # or tRF-Gly-GCC
		{
		@a=split(":", $_);
		$a[9]=~s/ kJ.+$//;
		
		@b=split(" ", $_);
		$hash{$b[1]}+=$a[9];	#hier wird dem key(Element1) sein spezifischer value(Element3) zugeordnet
		}
	}
close FH;

open(OUT,">tRNAscore_PR_Glu"); # or >tRNAscore_PR_Gly
foreach(sort{$hash{$a}<=>$hash{$b}}keys%hash)	#spezielle Variablen: $a ist die, mit dem kleinsten Wert, $b mit dem größten.
	{
	print OUT"$_\t$hash{$_}\n";
	}
close OUT;

