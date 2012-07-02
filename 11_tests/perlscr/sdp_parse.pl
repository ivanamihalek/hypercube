#! /usr/bin/perl -w

( -e "capra.table") || die "need capra.table as a ref.\n";

while (<>) {
    last if (/Project/);
}
while (<>) {
    last if (/Alignment position /);
}

$method = "sdp";

while (<>) {
    next if (!/\d/);
    @aux = split;
    $aux[2] =~ /\d(\D)/;
    $type = $1;
    $pos = $aux[1];
    printf  "%4d  %s   %8.3f\n", $pos, $type,  $aux[4];
    $score{$method}{$pos} = $aux[4];
    $type{$method}{$pos}  = $type;

}


@lines = split "\n", `cat capra.table`; 

foreach $method ("sdp") {
    $file = "$method.table";
    print "$file\n";
    open (OF, ">$file") || die "Cno $file: $!\n";
    foreach $line (@lines) {
	($pos, $type) = split " ", $line;
	if ( defined $score{$method}{$pos} ) {
	    ($type eq $type{$method}{$pos}) || die "type mmatch:$line\n";
	    printf OF " %3d  %s  %8.3f\n", $pos, $type, $score{$method}{$pos};
	} else {
	    printf OF " %3d  %s  %8.3f\n", $pos, $type, 2.0;
	}
	
    }
    close OF;
} 
