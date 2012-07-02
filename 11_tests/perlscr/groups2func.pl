#! /usr/bin/perl -w

# func is the input fromat used by s3det

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

print "\t";
print join "\t", @names;
print "\n\n";
foreach $seq (@seqs) {
    print "$seq";
    foreach $name (@names){
	if ( $belongs_to{$seq} eq $name ) {
	    print "\t1";
	} else {
	    print "\t0";
	}
    }
    print "\n";
}
