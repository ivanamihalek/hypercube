#! /usr/bin/perl -w
 

use Switch;
sub set_palettes();


(@ARGV >= 3) ||
    die "Usage:  $0   <specs score file>  ".
    "<pdb_file_full_path>  <output name> [-c <chain> ] [ -g <group1> <group2> ...]   \n".
    "-g (optional) groups for which to find conservation and determinants\n"; 

$method = "entropy";

$ranks_file = shift @ARGV;
$pdb_file   = shift @ARGV;
$outname    = shift @ARGV;

$chain = "";
@groups = ();

if ( @ARGV ) {

    $switch = shift @ARGV;

    if ( $switch eq "-c" ) {
	$chain = shift @ARGV;
	if ( @ARGV ) {
	    $switch = shift @ARGV;
	     if ( $switch eq "-g" ) {
		 @groups = map {lc $_} @ARGV; # whatever is left
	     }
	}

    } elsif ( $switch eq "-g" ) {
	while ( @ARGV ) {
	    $tmp = shift @ARGV;
	    if ( $tmp  eq "-c" ) {
		$chain =  shift @ARGV;
		last;
	    } else {
		push @groups, lc $tmp;
	    }
	}
    }
}




$COLOR_RANGE = 20;

set_palettes();


##################################################
# input/ouput
##################################################

# open the output file
$filename = $outname;

open (FPTR, ">$filename") || die "cno $filename\n";

if ($chain) {
    print FPTR "load $pdb_file, the_whole_thing\n";
    print FPTR "create conservation, the_whole_thing and chain $chain\n";
} else {
    print FPTR "load $pdb_file, conservation\n";
}


print FPTR "color white, conservation\n";
print FPTR "hide lines,  conservation\n";
print FPTR "create specificity, conservation\n";
foreach $group (@groups) {
    print FPTR "create cns_$group, conservation\n";
    print FPTR "create dts_$group, conservation\n";
}

print FPTR "\n";
=pod
foreach (@group_name) {
    print FPTR  "create $_, conservation\n";
    print FPTR "color white, $_\n";
    print FPTR "hide lines, $_\n";
}
=cut
print FPTR "\n";
print FPTR "hide everything\n";
print FPTR "bg_color white\n";

print  FPTR "\n";

for $ctr ( 0 .. $#color ) {
    print  FPTR "set_color c$ctr = $color[$ctr]\n";
}
print  FPTR "\n";


for ( $ctr= 0; $ctr <= $COLOR_RANGE/2; $ctr++ ) {
    print  FPTR "set_color orange$ctr = $orange_range[$ctr]\n";
}
print  FPTR "\n";

for ( $ctr= 0; $ctr <= $COLOR_RANGE/2; $ctr++ ) {
    print  FPTR "set_color berry$ctr  = $berry_range[$ctr]\n";
}
print  FPTR "\n";

for ( $ctr= 0; $ctr <= $COLOR_RANGE/2; $ctr++ ) {
    print  FPTR "set_color blue$ctr   = $blue_range[$ctr]\n";
}
print  FPTR "\n";

print FPTR "# color by conservation \n";


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
	    foreach $group (@groups) {
		( $title =~ $group) || next;
		if ( $title eq lc $group) {
		    $column{"cons_$group"} = $col;
		} elsif ( $title =~ "dets"){
		    $column{"dets_$group"} = $col;
		}
	    }
	}

	foreach $group (@groups) {
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
	next if ( $pdb_id eq ".");

        # conservation:
	$i   = $column{$method};
	$cvg = $aux[$i];
	$color_index = int ( $cvg*$COLOR_RANGE );
	$cn  = "c$color_index";
	print FPTR "color $cn, conservation and resid $pdb_id ";
	( $chain ) && print  FPTR  "and chain $chain";
	print  FPTR "\n";
	

	# conservation in each group
	foreach $group (@groups) {
	    $i   = $column{"cons_$group"};
	    $cvg = $aux[$i];
	    $color_index = int ( $cvg*$COLOR_RANGE );
	    $cn  = "c$color_index";
	    print FPTR "color $cn, cns_$group and resid $pdb_id ";
	    ( $chain ) && print  FPTR  "and chain $chain";
	    print  FPTR "\n";
	}


	# discriminants:
	$i   = $column{"discr"};
	$cvg = $aux[$i];
	if ( $cvg <= 0.5) {
	    $color_index = int ( (0.5-$cvg)*$COLOR_RANGE );
	    $cn = "orange$color_index";
	} else {
	    $color_index = int ( ($cvg-0.5)*$COLOR_RANGE );
	    $cn = "blue$color_index";
	}
	print FPTR "color $cn, specificity and resid $pdb_id ";
	( $chain ) && print  FPTR  "and chain $chain";
	print  FPTR "\n";

	# determinants
	foreach $group (@groups) {
	    $i   = $column{"dets_$group"};
	    $cvg = $aux[$i];
	    if ( $cvg <= 0.5) {
		$color_index = int ( (0.5-$cvg)*$COLOR_RANGE );
		$cn = "orange$color_index";
	    } else {
		$color_index = int ( ($cvg-0.5)*$COLOR_RANGE );
		$cn = "blue$color_index";
	    }
	    print FPTR "color $cn, dts_$group and resid $pdb_id ";
	    ( $chain ) && print  FPTR  "and chain $chain";
	    print  FPTR "\n";
	}

   }

    $line++;

}

close RANKS_FILE;


print  FPTR "\n";
print  FPTR  "deselect\n";
print  FPTR  "zoom complete=1\n";
foreach $group (@groups) {
    print  FPTR  "show spheres, cns_$group\n";
    print  FPTR  "show spheres, dts_$group\n";
}
print  FPTR  "show spheres, conservation\n";
print  FPTR  "show spheres, specificity\n";



if ( $chain ) {
 
    print  FPTR "select heteroatoms,  hetatm and not solvent \n";
    print  FPTR "select other_chains, not chain $chain \n";
    print  FPTR "select struct_water, solvent and chain $chain \n";
    print  FPTR "select metals, symbol  mg+ca+fe+zn+na+k+mn+cu+ni+cd+i\n";

    print  FPTR "cartoon putty \n";
    print  FPTR "show  cartoon,  other_chains \n";
    print  FPTR "hide  spheres,   heteroatoms \n";
    print  FPTR "show  sticks,   heteroatoms \n";
    print  FPTR "show  spheres,  struct_water \n";
    print  FPTR "show  spheres,  metals \n";
    print  FPTR "color palecyan, struct_water \n";
    print  FPTR "color lightteal, other_chains or heteroatoms \n";
    print  FPTR "color magenta, metals\n";
    print  FPTR "zoom  chain $chain\n";
}

print  FPTR "deselect \n";

$session_file = $outname;
if ($outname =~ /\.pml$/) {
    $session_file =~ s/\.pml$/\.pse/;
} else {
    $session_file .= ".pse";
    
}
print  FPTR "save $session_file \n";



close FPTR; 




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

