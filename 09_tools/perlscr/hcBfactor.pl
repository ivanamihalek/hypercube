#! /usr/bin/perl -w
 

use Switch;
sub set_palettes();


(@ARGV >= 3) ||
    die "Usage:  $0   <specs score file>  ".
    "<pdb_file_full_path>  <output name> [-c <chain> ] [ -g <group1> <group2> ...]   \n".
    "-g (optional) groups for which to find conservation and determinants\n"; 

$method = "rvet";

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
		 @groups = grep {/lc $_/} @ARGV; # whatever is left
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



##################################################
# input/ouput pdb
##################################################


# read in the score file and output the colors accordingly
open (RANKS_FILE, "<$ranks_file") || 
    die "cno $ranks_file\n";
    

$line = 1;

$number_of_groups = 0;
%column = ();

%cons_cvg = ();

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
	$cons_cvg{$pdb_id} = $aux[$i];
	

=pod
	# conservation in each group
	foreach $group (@groups) {
	    $i   = $column{"cons_$group"};
	    $cvg = $aux[$i];
	}

	# discriminants:
	$i   = $column{"discr"};
	$cvg = $aux[$i];
	if ( $cvg <= 0.5) {
	    $color_index = int ( (0.5-$cvg)*$COLOR_RANGE );
	    #orange
	} else {
	    $color_index = int ( ($cvg-0.5)*$COLOR_RANGE );
	    # blue
	}


	# determinants
	foreach $group (@groups) {
	    $i   = $column{"dets_$group"};
	    $cvg = $aux[$i];
	    if ( $cvg <= 0.5) {
		# orange
	    } else {
		# blue
	    }
	}
=cut
   }

    $line++;

}

close RANKS_FILE;


# open the output file
$filename = $pdb_file;
open (IN_PDB, "<$filename") || die "cno $filename\n";


$filename = $outname;
open (OUT_PDB, ">$filename") || die "cno $filename\n";

while(<IN_PDB>){
    if(!/^ATOM/){
        print  OUT_PDB "$_";
    }
    else{
        $chain_pdb = substr($_, 21,1);
        if($chain ne ""){
            if($chain eq $chain_pdb){
                $pos = substr($_, 22,4);
                $pos =~ s/\s*//g;
		if (defined $cons_cvg{$pos} ) {
		    $score = $cons_cvg{$pos}*100;
		} else {
		    $score = 100.0;
		}
		$newBfactor = sprintf(" %6.2f", $score);
                $newline = $_;
                substr($newline, 59,7, $newBfactor);
                print OUT_PDB $newline; 
           } else{
               print OUT_PDB "$_";
           } 
        } else{
            $pos = substr($_, 22,4);
            $pos =~ s/\s*//g;
	    if (defined $cons_cvg{$pos} ) {
		$score = $cons_cvg{$pos}*100;
	    } else {
		$score = 100.0;
	    }
            $newBfactor = sprintf(" %6.2f", $score);
            $newline = $_;
            substr($newline, 59,7,$newBfactor);
            print OUT_PDB $newline;
        }  
    }
}


close IN_PDB; 
close OUT_PDB; 




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

