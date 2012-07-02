#!/usr/bin/perl -w
@ARGV || die "Usage: $0 <tp_list> <tn_list> <capra out>  \n";

$tplist    = shift @ARGV;
$tnlist     = shift @ARGV;
$score_file = shift @ARGV;


#################################
$file = $tplist;
open (IF, "<$file")
    || die "cno $file: $!\n";
@tp_ids = ();
%pos_aa_type = ();
while ( <IF>) {
    next if ( !/\S/);
    chomp;
    ($type, $id) = split;
    $id --;
    push @tp_ids, $id;
    $pos_aa_type{$id} = $type;
}
close IF;

#################################
if ( $tnlist ne "na" ) { # otherwise, we'll take that what
                       # is not true positive is true negative

    $file = $tnlist;
    open (IF, "<$file")
	|| die "cno $file: $!\n";
    @tn_ids = ();
    %neg_aa_type = ();
    while ( <IF>) {
	next if ( !/\S/);
	chomp;
	($type, $id) = split;
	$id --;
	push @tn_ids, $id;
	$neg_aa_type{$id} = $type;
    }
    close IF;
}



#################################

open (IF, "<$score_file") || die "Cno $score_file: $!\n";
@positions = ();
while (<IF> ) {
    next if (!/\S/);
    next if (/^#/);
    ($pos, $score) = split;
    if ($score eq "None") {
	$score = -1;
    }
    $capra[$pos] = $score;
    push @positions, $pos;
}

@sorted = sort { $capra[$a] <=> $capra[$b] } @positions;


foreach $pos (@sorted) {
    if ( $tnlist ne "na" ) {
	if ( defined   $pos_aa_type{$pos} ) {

	    $positive[$pos] = 1;
	} else {
	    $positive[$pos] = 0;
	}

	if ( defined  $neg_aa_type{$pos} ) {

	    $negative[$pos] = 1;
	} else {
	    $negative[$pos] = 0;
	}
    } else {
	if ( defined   $pos_aa_type{$pos} ) {

	    $positive[$pos] = 1;
	    $negative[$pos] = 0;
	} else {
	    $positive[$pos] = 0;
	    $negative[$pos] = 1;
	}
    }
    
}

#################################
$tp[0] = $positive[0];
$fp[0] = $negative[0];
for $ctr ( 1 .. $#sorted) {

    $tp[$ctr] = $positive[$ctr] + $tp[$ctr-1]; 
    $fp[$ctr] = $negative[$ctr] + $fp[$ctr-1]; 

}


#################################
$file = "capra.roc";
$file =~ s/\.score//;

open (OF, ">$file") || 
    die "Cno $file:$! \n";

 
$prev_fp_fraction = 0.0;
$roc_area = 0.0;
for $ctr (0 .. $#sorted) {
    $tp_fraction = $tp[$ctr]/$tp[$#sorted];
    $fp_fraction = $fp[$ctr]/$fp[$#sorted];
    $roc_area += $tp_fraction*($fp_fraction-$prev_fp_fraction); 
    $prev_fp_fraction = $fp_fraction;
    printf OF " %8.3lf  %8.3lf \n", $fp_fraction, $tp_fraction;
    
}
close OF;
