#! /usr/bin/perl -w
sub print_hdr();

@ARGV || die "Usage: $0 <tp_list> <cube score (filename root)> [-c]  [method list]\n";

$tp_list    = shift @ARGV;
$score_file = shift @ARGV;

$check_type = 0;
if ( @ARGV  && $ARGV[0] eq "-c") {
    $check_type = 1;
    shift @ARGV;
}

if ( @ARGV ) {
    $method_subset = join " ", @ARGV;
} else {
    $method_subset = "";
}


@target_residues = ();
open (IF, "<$tp_list") || die "Cno $tp_list: $!\n";
while ( <IF> ) {
    next if ( !/\S/);
    chomp;
    @aux = split;
    push @target_residues, $aux[1];
    $type{$aux[1]} = $aux[0];
}
close IF;


@header_cols = split " ", `head -n1 $score_file.comp_score`;
shift @header_cols;  # get rid of the comment sign

@methods = @header_cols;
$offset = 1;
while ( (shift @methods) ne "gaps") {$offset++};


foreach $i(0 .. $#methods) {
    $col{$methods[$i]} = $i + $offset;

}
foreach $method (@methods) {
    $avg{$method} = 0;
    $std{$method} = 0;
}


#print "$offset    @methods\n"; exit;
print_hdr();

$tot = 0;
foreach $tgt (@target_residues) {

    @aux = split " ", `awk \'\$1==$tgt\' $score_file.comp_score`;
    next if ($aux[3] =~ "-");

    if ( $check_type ) {
	( $type{$tgt} eq $aux[1] ) ||
	    die "type mismatch for $tgt $type{$tgt} in $tp_list\n";
    }

    printf " %4d %s  ", $tgt, $type{$tgt};
    foreach $method (@methods) {
	next if ( $method_subset && $method_subset !~ $method);
	$cvg = $aux[$col{$method}];
	($cvg eq "-") && ($cvg = 1.0);
	printf "%7.2f", $cvg;
	$avg{$method} += $cvg;
	$std{$method} += $cvg*$cvg
    }


    $tot ++;
    printf "\n";
    
}

printf "%% ----------------------------------\n";
print_hdr();

printf "%% %4s   ", "avg";
foreach $method (@methods) {
    next if ( $method_subset && $method_subset !~ $method);
     $avg{$method} /=  $tot;
     printf "%7.2f", $avg{$method};
}
printf "\n";
printf "%% %4s   ", "std";
foreach $method (@methods) {
    next if ( $method_subset && $method_subset !~ $method);
     $std{$method} /= $tot;
     printf "%7.2f", 
     sqrt($std{$method} - $avg{$method}*$avg{$method});
}
printf "\n";



############################
sub print_hdr() {
    printf "%% %4s    ", "";
    foreach $method (@methods, "rvet") {
	next if ( $method_subset && $method_subset !~ $method);
	printf "%7s", $method; 
    }
    printf "\n";
}
