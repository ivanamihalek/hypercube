#! /usr/bin/perl -w

@ARGV ||
    die "Usage:  $0  <file name> \n";

$filename = $ARGV[0];
open (IF, "<$filename" ) 
    || die "Cno $filename: $!.\n";

while ( <IF> ) {
    next if ( ! /\S/ );
    next if (  /^#/ );
    chomp;
    ($id, $score, $type) = split;
    ( $score eq "None" )  && ($score = -1);
    printf "%4d   %s   %8.5f\n", $id, (substr $type, 0, 1), -$score;
}

close IF;
