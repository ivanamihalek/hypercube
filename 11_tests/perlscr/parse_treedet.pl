#! /usr/bin/perl -w

(@ARGV > 1) ||
    die "Usage:  $0  <treedet file>  <rep seqs>\n";

($filename, $repseq) = @ARGV;

@lines = split "\n", `cat $filename`;
 
( -e "capra.table") || die "need capra.table as a ref.\n";

$ctr = 0;
foreach $line (@lines) {
    $ctr++;
    last if ($line =~ /<methods>/);
}


$reading = 0;
$prot_ctr = 0;
foreach $line (@lines[$ctr .. $#lines]) {

    if ( $line  =~ /method name/ ) {
	$line   =~ /method name=\"(.+?)\"/;
	$method =  $1;
	print "\n$method\n";
	$prot_ctr = 0;
    } elsif ( $line =~ /protein name/) {
	$prot_ctr ++;
	if ( $line =~ $repseq) {
	    $reading = 1;
	    print "$line\n";
	}
    } elsif ( $line =~ /<\/protein>/) {
	$reading = 0;

    } elsif ( $line =~ /\/protein_list/ ) {

    } elsif ( $reading && $line =~ /residue pos=/ ) {
	$line =~ /residue pos=\"(\d+?)\" type=\"(\w)\"/;
	$pos = $1; $type = $2;
	printf "\t %3d  %s ",  $pos, $type;
	if ($method eq "S3DET") {
	    print "\n";
	    $score{$method}{$pos} = -1;
	}
	$type{$method}{$pos} = $type;

    } elsif ( $reading && $line =~ /score type=/ ) {
	$line =~/value=\"(.+?)\"/;
	$evalue = $1;
	printf "  %s  \n", $evalue;
	$score{$method}{$pos} = $evalue;
	
    }
}


@lines = split "\n", `cat capra.table`; 

foreach $method (keys %score) {
    $file = "$method.table";
    print "$file\n";
    open (OF, ">$file") || die "Cno $file: $!\n";
    foreach $line (@lines) {
	($pos, $type) = split " ", $line;
	$pos++;

	if ( defined $score{$method}{$pos} ) {
	    ($type eq $type{$method}{$pos}) || die "type mmatch:$line\n";
	    printf OF " %3d  %s  %8.3f\n", $pos-1, $type, $score{$method}{$pos};
	} else {
	    printf  OF " %3d  %s  %8.3f\n", $pos-1, $type, 2.0;
	}
	
    }
    close OF;
} 
