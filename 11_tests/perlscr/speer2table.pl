#! /usr/bin/perl -w

@ARGV ||
    die "Usage:  $0  <file name> \n";

( -e "capra.table") || die "need capra.table as a ref.\n";


$filename = $ARGV[0];
open (IF, "<$filename" ) 
    || die "Cno $filename: $!.\n";

while ( <IF> ) {
    last if ( /SPEER Results/);
}
while ( <IF> ) {
   next if ( ! /\S/ );
   next if (  /^#/ );
    chomp;
    @aux = split;
    #($aux[8] == -2) && ($aux[8] = 2);
    #$score{$aux[0]} = -$aux[2];
    $score{$aux[0]} = $aux[4];
}

close IF;

@lines = split "\n", `cat capra.table`; 



foreach $line (@lines) {
    ($pm, $type) = split " ", $line;
    $pos = $pm+1;
    if ( defined $score{$pos} ) {
	printf  " %3d  %s  %8.3f\n", $pm, $type, $score{$pos};
    } else {
	printf " %3d  %s  %8.3f\n", $pm, $type, 2.0;
    }
	
}
