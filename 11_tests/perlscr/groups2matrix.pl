#! /usr/bin/perl -w

# matrix is the input fromat used by Xdet

@ARGV ||
    die "Usage:  $0  <file name> \n";

$filename = $ARGV[0];
open (IF, "<$filename" ) 
    || die "Cno $filename: $!.\n";

@names = ();
@seqs  = ();
while ( <IF> ) {
    next if ( ! /\S/ );
    chomp;
    @aux = split;
    if ( $aux[0] eq "name") {
	$name = $aux[1];
	push @names, $name;
    } else {
	$belongs_to{$aux[0]}= $name;
	push @seqs, $aux[0];
    }
}

close IF;

foreach $seq1 (@seqs) {
    foreach $seq2 (@seqs) {
	print "$seq1\t$seq2\t";
	if ( $belongs_to{$seq1} eq  $belongs_to{$seq2}) {
	    print "1.0";
	} else {
	    print "0.0";
	}
	print "\n";
    }
}
