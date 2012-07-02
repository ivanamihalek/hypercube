#! /usr/bin/perl -w

@ARGV ||
    die "Usage:  $0  <file name> \n";

$filename = $ARGV[0];
open (IF, "<$filename" ) 
    || die "Cno $filename: $!.\n";
$tag = "";

while ( <IF> ) {
    next if ( ! /\S/);
    chomp;
    @aux = split;

    if ($aux[0] eq "name") {
	$tag = $aux[1];
    } else {
	print "$aux[0]  $aux[0]\_$tag\n";  
    }

}

close IF;
