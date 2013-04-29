#!/usr/bin/perl -w

(@ARGV >= 3 ) ||
    die "Usage: $0 <true_pos_list>  <true_neg_list>  <cube_out>  <outname.pdf> [-l] [methods to plot] \n";



$tplist = shift @ARGV;
$tnlist = shift @ARGV;
$cube   = shift @ARGV;
$pdf    = shift @ARGV;
$legend = 0;
if ( $ARGV[0] eq "-l") {
     shift @ARGV;
     $legend = 1;
}



foreach ( $tplist, $cube) {
    (-e $_) || die "$_ not found\n";
}

@methods_to_plot = ();
if ( @ARGV ) {
    @methods_to_plot = @ARGV;
}

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
	push @tn_ids, $id;
	$neg_aa_type{$id} = $type;
    }
    close IF;
}

################################
################################

$title = $cube;
$title =~ s/\.score//;
 
(-e "rocs") || `mkdir rocs`;
(-e "sorted") || `mkdir sorted`;

#################################
@header_cols = split " ", `head -n1 $cube`;
shift @header_cols;  # get rid of the comment sign

@methods = @header_cols;
$offset = 1;
do {
    $offset++;
} while ( (shift @methods) ne "gaps");

foreach $i(0 .. $#methods) {
    $col{$methods[$i]} = $i + $offset;
    #print "$i  $offset  $methods[$i] \n";
}


if (scalar @methods_to_plot){

    foreach $method (@methods_to_plot) {
	(defined $col{$method}) ||
	    die "$method not found in  $cube.\n";
	#print "$method $col{$method} \n";
    }

    @methods = @methods_to_plot;


}


#################################
#################################
foreach $method (@methods ) {

    #################################
    $file   = "$title.$method.sorted";

    $column = $col{$method};


    # 2 indicates that this was a position that was gap in the majority of considered families
    
    $cmd = "grep -v gaps  $cube | awk \'\$$column != \"-\" && \$$column < 2.0 \' |  sort -gk $column ".
	" | awk \'{printf \"%3d\", NR}{print}\' > $file";

    (system $cmd ) && die "Error running $cmd\n";
 
    $number_of_points = 0;

    open (IF, "<$file")
    || die "cno $file: $!\n";
    $ctr = 0;
    ($rank, $id, $aa) = (0);
    while ( <IF>) {
	next if ( !/\S/);
	chomp;
	@aux = split;
	($rank, $id, $aa) = @aux;

	if ( $tnlist ne "na" ) {
	    if ( defined   $pos_aa_type{$id} ) {

		$positive[$ctr] = 1;
	    } else {
		$positive[$ctr] = 0;
	    }

	    if ( defined  $neg_aa_type{$id} ) {
		
		$negative[$ctr] = 1;
	    } else {
		$negative[$ctr] = 0;
	    }
	} else {
	    if ( defined   $pos_aa_type{$id} ) {

		$positive[$ctr] = 1;
		$negative[$ctr] = 0;
	    } else {
		$positive[$ctr] = 0;
		$negative[$ctr] = 1;
	    }
	}


	$ctr++;
    }
    $number_of_points = $ctr;
    $last_index = $number_of_points -1;
    close IF;
    `mv $file sorted`;


    #################################
    $tp[0] = $positive[0];
    $fp[0] = $negative[0];
    for $ctr ( 1 .. $last_index) {

	$tp[$ctr] = $positive[$ctr] + $tp[$ctr-1]; 
	$fp[$ctr] = $negative[$ctr] + $fp[$ctr-1]; 

    }


    #################################
    $file = "$title.$method.roc";
    $file =~ s/\.score//;

    open (OF, ">$file") || 
	die "Cno $file:$! \n";

 
    $prev_fp_fraction = 0.0;
    $roc_area = 0.0;
    for $ctr (0 ..$last_index) {
	$tp_fraction = $tp[$ctr]/$tp[$last_index];
	$fp_fraction = $fp[$ctr]/$fp[$last_index];
	$roc_area += $tp_fraction*($fp_fraction-$prev_fp_fraction); 
	$prev_fp_fraction = $fp_fraction;
	printf OF " %8.3lf  %8.3lf \n", $fp_fraction, $tp_fraction;
    
    }
    close OF;
    `mv $file rocs`;


    ##############################
    $auc{$method} =  $roc_area;
}

############################
# make graph using gnuplot


$file = "$title.gscr";   
open (OF, ">$file") || 
    die "Cno $file:$! \n";

print OF "set size square \n";
if ( $legend) {
    print OF "set key bottom right \n";
} else {
    print OF "unset key\n";
}
print OF "set style data lines \n";
print OF "set xlabel \"FP rate\" \n";
print OF "set ylabel \"TP rate\" \n";
#print OF "set title \"$title\" \n";
print OF "set term post color \n";
print OF "set output \'tmp.ps\' \n";
print OF "plot ";

################################
$first = 1;
foreach $method (@methods ) {
    $file = "$title.$method.roc";
   
    ( -e "rocs/$file" )|| next;
    ($first) || print OF ", ";
    print OF " \'rocs/$file\' title \'$method\'";
    $first = 0;
}

print OF "\n";
print OF "quit\n";
close OF;


`gnuplot $title.gscr`;
`ps2pdf tmp.ps $pdf`;
`rm tmp.ps`;


################################
# output the list of methods sorted by auc
@meth_sorted_by_auc = sort { $auc{$b} <=> $auc{$a} } @methods;
    
$file = "$cube.auc";
open (OF, ">$file") || 
    die "Cno $file:$! \n";
foreach $method (@meth_sorted_by_auc) {
    printf OF " %20s   %8.3f\n", $method, $auc{$method};
}


close OF;

