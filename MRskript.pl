#!usr/bin/perl

open(FH,"miRandaoutput_Glu") or die print "can't open file: $!"; # or miRandaoutput_Gly
while(<FH>)
	{
	if($_=~/^>>/)
		{
		@a=split("\t", $_);
		$hash{$a[1]}=$a[3];	#hier wird dem key(Element1) sein spezifischer value(Element3) zugeordnet
		}
	}
close FH;

open(OUT,">tRNAscores_Glu"); # or >tRNAscores_Gly
foreach(sort{$hash{$a}<=>$hash{$b}}keys%hash)	#spezielle Variablen: $a ist die, mit dem kleinsten Wert, $b mit dem größten.
	{
	print OUT"$_\t$hash{$_}\n";
	}
close OUT;












 


