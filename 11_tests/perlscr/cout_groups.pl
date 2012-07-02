#! /usr/bin/perl -w

# func is the input fromat used by s3det

@ARGV ||
    die "Usage:  $0  <file name> \n";

$filename = $ARGV[0];
open (IF, "<$filename" ) 
    || die "Cno $filename: $!.\n";

@names = ();

while ( <IF> ) {
    next if ( ! /\S/ );
    chomp;
    @aux = split;
    if ( $aux[0] eq "name") {
	$name = $aux[1];
	push @names, $name;
	$count{$name} = 0;
    } else {
	$count{$name} ++;
    }
}

close IF;


foreach $name (@names){

    #print "$name  $count{$name}\n";
    print "  $count{$name} ";
}
print "\n";
