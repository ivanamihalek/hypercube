#! /usr/bin/perl -w

@ARGV ||
    die "Usage:  $0  <file name> \n";

$filename = $ARGV[0];
open (IF, "<$filename" ) 
    || die "Cno $filename: $!.\n";

while ( <IF> ) {
    next if ( ! /\S/ );
    chomp;
    @aux = split;
    #($aux[8] == -2) && ($aux[8] = 2);
    printf "%4d   %s   %8.3f\n", $aux[0]-1, $aux[1], $aux[8];
}

close IF;
