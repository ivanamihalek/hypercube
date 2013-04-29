#! /usr/bin/perl -w
 

use Switch;
sub set_palettes();


(@ARGV >= 3) ||
    die "Usage:  $0   <specs score file>  ".
    "<pdb_file_full_path>  <output name> [-c <chain> ] [-g <group>]   \n".
    "-g (optional) group for which to find conservation and determinants\n"; 

$method = "rvet";

$ranks_file = shift @ARGV;
$pdb_file   = shift @ARGV;
$outname    = shift @ARGV;

$chain = "";
$group = "";

if ( @ARGV ) {

    $switch = shift @ARGV;

    if ( $switch eq "-g" ) {
 	$group =  shift @ARGV;
    } elsif ( $switch eq "-c" ) {
	$chain =  shift @ARGV;
    }
}

($group) && ( $group = lc $group);


$COLOR_RANGE = 20;

set_palettes();


##################################################
# input/ouput
##################################################

# open the output file
$filename = "tmp";

open (FPTR, ">$filename") || die "cno $filename\n";
#print FPTR "load $pdb_file, conservation\n";
print FPTR "open $pdb_file\n";
#print FPTR "color white, conservation\n";
print FPTR "color white\n";


print FPTR "\n";

#print FPTR "bg_color white\n";
print FPTR "set bg_color white\n";

print  FPTR "\n";


for $ctr ( 0 .. $#color ) {
    #print  FPTR "set_color c$ctr = $color[$ctr]\n";
    my $color_scheme = $color[$ctr];
    $color_scheme =~ s/(\[|\]|\,)/ /g;
    print FPTR "colordef c$ctr $color_scheme\n";
}
print  FPTR "\n";


for ( $ctr= 0; $ctr <= $COLOR_RANGE/2; $ctr++ ) {
    #print  FPTR "set_color orange$ctr = $orange_range[$ctr]\n";
    my $color_scheme = $orange_range[$ctr];
    $color_scheme =~ s/\[|\]|\,/ /g;
    print FPTR "colordef orange$ctr $color_scheme\n";
}
print  FPTR "\n";


for ( $ctr= 0; $ctr <= $COLOR_RANGE/2; $ctr++ ) {
    #print  FPTR "set_color berry$ctr  = $berry_range[$ctr]\n";
    my $color_scheme = $berry_range[$ctr];
    $color_scheme =~ s/\[|\]|\,/ /g;
    print FPTR "colordef berry$ctr $color_scheme\n";
}
print  FPTR "\n";

for ( $ctr= 0; $ctr <= $COLOR_RANGE/2; $ctr++ ) {
    #print  FPTR "set_color blue$ctr   = $blue_range[$ctr]\n";
    my $color_scheme = $blue_range[$ctr];
    $color_scheme =~ s/\[|\]|\,/ /g;
    print FPTR "colordef blue$ctr $color_scheme\n";
}
print  FPTR "\n";

close FPTR;


$cons_file = "$outname\_cons.com";
$spec_file = "$outname\_spec.com";

`cp tmp $cons_file`;
`cp tmp $spec_file`;

open (CONS_FPTR, ">>$cons_file") || die "cno $cons_file: $!\n";
open (SPEC_FPTR, ">>$spec_file") || die "cno $spec_file: $!\n";


##########################################################
# read in the score file and output the colors accordingly
open (RANKS_FILE, "<$ranks_file") || 
    die "cno $ranks_file\n";
    

$line = 1;

$number_of_groups = 0;
%column = ();

while ( <RANKS_FILE> ) {
    next if ( !/\S/ );
    if ( /\%/ ){
	@aux = split;
	shift @aux;
	
	for  $col (0 .. $#aux) {

	    $title = lc $aux[$col];

	    foreach ("discr",  "$method", "pdb_id", "pdb_aa" ) {
		( $title eq $_) || next;
		$column{$_} = $col;
		last;
	    }
	    if ( $title =~ $group) {
		if ( $title eq lc $group) {
		    $column{"cons_$group"} = $col;
		} elsif ( $title =~ "dets"){
		    $column{"dets_$group"} = $col;
		}
	    }
	}
	
	if ( $group ) {
	( defined  $column{"dets_$group"} &&  defined  $column{"cons_$group"} ) ||
	    die "$group not found in $ranks_file.\n";
	}

	(defined $column{"pdb_id"}) || die "No pdb_id in the output (?!).\n";
	(defined $column{"pdb_aa"}) || die "No pdb_aa in the output (?!).\n";
	(defined $column{$method})  || die "No $method in the output (?!).\n";
	(defined $column{"discr"})  || die "No discr in the output (?!).\n";

	$number_of_groups  =  scalar ( grep {/dets_/} @aux);

 	$line = 2;

   } else {

	chomp;
	@aux = split;

 	# the rest of the input file
	$i      = $column{"pdb_id"};
	$pdb_id = $aux[$i];
	next if ($pdb_id eq ".");

  	

	
	###################################
	# conservation
	if ( $group ) {    
	    # conservation in the  group
	    $i   = $column{"cons_$group"};
	} else {
	    # overall conservation:
	    $i   = $column{$method};
	}
	$cvg = $aux[$i];
	$color_index = int ( $cvg*$COLOR_RANGE );
	$cn  = "c$color_index";
	#print FPTR "color $cn, cns_$group and resid $pdb_id ";
	#( $chain ) && print  FPTR  "and chain $chain";
	#print  FPTR "\n";
	print CONS_FPTR "color $cn : $pdb_id";
	($chain) && print CONS_FPTR ".$chain";
	print CONS_FPTR "\n";

	###################################
	# specialization
	if ( $group ) {    
	    # determinants
	    $i   = $column{"dets_$group"};
	} else {
	    # discriminants:
	    $i   = $column{"discr"};
	}
	$cvg = $aux[$i];
	if ( $cvg <= 0.5) {
	    $color_index = int ( (0.5-$cvg)*$COLOR_RANGE );
	    $cn = "orange$color_index";
	} else {
	    $color_index = int ( ($cvg-0.5)*$COLOR_RANGE );
	    $cn = "blue$color_index";
	}
	print SPEC_FPTR "color $cn : $pdb_id";
	($chain) && print SPEC_FPTR ".$chain";
	print SPEC_FPTR "\n";


   }

    $line++;

}


close RANKS_FILE;


printf CONS_FPTR "save $outname\_cons.com\n";
printf SPEC_FPTR "save $outname\_spec.com\n";

close CONS_FPTR; 
close SPEC_FPTR; 




##################################################
##################################################
##################################################

##################################################
##################################################
##################################################



sub set_palettes () {

    ##################################################
    #set the pallette:
    $green = $blue = $red = 0;


    $N = 5;
    $C1 = $COLOR_RANGE-1;

    $red   = 1.00;
    $green = 0.87;
    $blue  = 0.0;
    $color[0]    = "[$red, $green, $blue]"; 
 

    $bin_size = $C1/$N;
    for ( $ctr=1; $ctr <= int ($COLOR_RANGE/$N); $ctr++ ) {

	$ratio =  ( int ( 100*($bin_size- $ctr+1)/$bin_size) ) /100;
	$red   = $ratio;
	$green = $blue = 0;
		 
	$color[$ctr] = "[$red, $green, $blue]"; 
    }

    for ( $ctr= int ($COLOR_RANGE/$N)+1 ; $ctr <= $COLOR_RANGE; $ctr++ ) {

	$ratio =  ( $ctr -  $COLOR_RANGE/$N)/ ($COLOR_RANGE*($N-1)/$N);
	$red = $ratio;
	$green = $blue = $red;


  
	$color[$ctr] = "[$red, $green, $blue]"; 
    }

    $var_color_space_size = $ctr;

    ######## specificity colors
    # this is still  not general - 
    # for now number_of_groups = 4;


    $color_entry = $var_color_space_size+8;

    $color_entry ++;
    $orange_range[0] = $blue_range[0] = $berry_range[0] = "[1.0, 1.0, 1.0]";
	

    for ( $ctr= 1; $ctr <= $COLOR_RANGE/2; $ctr++ ) {

	$ratio = $ctr/($COLOR_RANGE/2) ;

	# orange
	$red   = 255;
	$green = 255 - (255-153)*$ratio;
	$blue  = 255 - (255- 51)*$ratio ;

	$color_entry ++;
	$orange_range[$ctr] = sprintf "[%6.3f, %6.3f, %6.3f]", $red/255, $green/255, $blue/255;
	

	# blue
	$red   = 255 - (255 -   0)*$ratio;
	$green = 255 - (255 -   0)*$ratio;
	$blue  = 255 - (255 - 128)*$ratio ;

	$color_entry ++;
	$blue_range[$ctr] = sprintf "[%6.3f, %6.3f, %6.3f]", $red/255, $green/255, $blue/255;	
  
	# berry
	$red   = 255 - (255 - 199)*$ratio;
	$green = 255 - (255 -  21)*$ratio;
	$blue  = 255 - (255 - 133)*$ratio ;

	$color_entry ++;
	$berry_range[$ctr] = 	sprintf "[%6.3f, %6.3f, %6.3f]", $red/255, $green/255, $blue/255;
  
    }

}

