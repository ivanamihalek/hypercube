#! /usr/bin/perl -w
use IO::Handle;         #autoflush
# FH -> autoflush(1);


(@ARGV < 2)  &&
    die "Usage: $0  <fasta >  <list of new names.\n"; 

($fasta, $new_names) = @ARGV;

@names = ();
open ( FASTA, "<$fasta") ||
    die "Cno $fasta: $!\n";

while ( <FASTA> ) {
    chomp;
    if (/^>(.+)/ ) {
	$name = $1;
	$name =~ s/\s//g;
	push @names,$name;
	$sequence{$name} = "";
    } else  {
	$sequence{$name} .= $_;
    } 
}
close FASTA;

open (NAMES, "$new_names")  ||
    die "Cno $new_names: $!\n";

while (<NAMES>) {
    next if ( !/\S/);
    chomp;
    ($old, $new) = split;
    $new_name{$old} = $new;
}
close NAMES;
	

foreach $seq_name ( @names ) {
	
    @seq = split ('', $sequence{$seq_name});
    print  ">$new_name{$seq_name} \n";
    $ctr = 0;
    for $i ( 0 .. $#seq ) {
	print  $seq[$i];
	$ctr++;
	if ( ! ($ctr % 100) ) {
	    print "\n";
	}

    }
    print "\n";
}

 

