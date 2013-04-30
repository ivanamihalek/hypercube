#! /usr/bin/perl -w
 
# perl -MCPAN -e "install 'Spreadsheet::WriteExcel'"
# alternatively (in Ubuntu)
# sudo apt-get install libSpreadsheet-WriteExcel-perl
# doc: 
# http://search.cpan.org/dist/Spreadsheet-WriteExcel/lib/Spreadsheet/WriteExcel.pm

use Switch;

use Spreadsheet::WriteExcel;                             # Step 0

(defined $ARGV[1]) ||
    die "Usage:  hc2xls.pl  <score file>  <output name>  \n"; 


($ranks_file, $outname) = @ARGV;

my $workbook = Spreadsheet::WriteExcel->new("$outname"); # Step 1
$worksheet   = $workbook->add_worksheet();                   # Step 2
$worksheet->keep_leading_zeros();


##################################################
#set the pallette:
$COLOR_RANGE = 20;
$green = $blue = $red = 0;


$N = 5;
$C1 = $COLOR_RANGE-1;

$red   = 1.00;
$green = 0.87;
$blue  = 0.0;
$color[0] = "[$red, $green, $blue]"; 
$xl_color[0] = $workbook->set_custom_color(8, $red*255, $green*255, $blue*255 );	


$bin_size = $C1/$N;
for ( $ctr=1; $ctr <= int ($COLOR_RANGE/$N); $ctr++ ) {

    $ratio =  ( int ( 100*($bin_size- $ctr+1)/$bin_size) ) /100;
    $red   = $ratio;
    $green = $blue = 0;
		 
    $color[$ctr] = "[$red, $green, $blue]"; 
    $color_name[$ctr] = "c$ctr";

    $xl_color[$ctr] = $workbook->set_custom_color($ctr+8, $red*255, $green*255, $blue*255 );	
    #printf "  %d   %d    %d , $ctr, $xl_color[$ctr] \n",$red*255, $green*255, $blue*255;
}

for ( $ctr= int ($COLOR_RANGE/$N)+1 ; $ctr <= $COLOR_RANGE; $ctr++ ) {

    $ratio =  ( $ctr -  $COLOR_RANGE/$N)/ ($COLOR_RANGE*($N-1)/$N);
    $red = $ratio;
    $green = $blue = $red;

    $xl_color[$ctr] = $workbook->set_custom_color($ctr+8, $red*255, $green*255, $blue*255 );	
    #printf "  %d   %d    %d , $ctr, $xl_color[$ctr] \n",$red*255, $green*255, $blue*255;

  
    $color[$ctr] = "[$red, $green, $blue]"; 
    $color_name[$ctr] = "c$ctr";
}

$var_color_space_size = $ctr;

######## specificity colors
# this is still  not general - 
# for now number_of_groups = 4;


$color_entry = $var_color_space_size+8;

$color_entry ++;
$orange_range[0] = $blue_range[0] = $berry_range[0] 
    = $workbook->set_custom_color($color_entry, 255, 255, 255);

for ( $ctr= 1; $ctr <= $COLOR_RANGE/2; $ctr++ ) {

    $ratio = $ctr/($COLOR_RANGE/2) ;

    # orange
    $red   = 255;
    $green = 255 - (255-153)*$ratio;
    $blue  = 255 - (255- 51)*$ratio ;

    $color_entry ++;
    $orange_range[$ctr] = $workbook->set_custom_color($color_entry, $red, $green, $blue);
    #printf "  %d   %d    %d , %8.3lf  $ctr\n",$red, $green, $blue, $ratio;
	

    # blue
    $red   = 255 - (255 -   0)*$ratio;
    $green = 255 - (255 -   0)*$ratio;
    $blue  = 255 - (255 - 128)*$ratio ;

    $color_entry ++;
    $blue_range[$ctr] = $workbook->set_custom_color($color_entry, $red, $green, $blue);	
    #printf "  %d   %d    %d , %8.3lf  $ctr\n",$red, $green, $blue, $ratio;
  
    # berry
    $red   = 255 - (255 - 199)*$ratio;
    $green = 255 - (255 -  21)*$ratio;
    $blue  = 255 - (255 - 133)*$ratio ;

    $color_entry ++;
    $berry_range[$ctr] = $workbook->set_custom_color($color_entry, $red, $green, $blue);	
    #printf "  %d   %d    %d , %8.3lf  $ctr\n",$red, $green, $blue, $ratio;
  
}


##################################################
# input/ouput

open (RANKS_FILE, "<$ranks_file") || 
    die "cno $ranks_file\n";
    

$format_centered = $workbook->add_format();
$format_centered->set_align('center');
#
$format_right = $workbook->add_format();
$format_right->set_align('right');
#
$format_left = $workbook->add_format();
$format_left->set_align('left');
#
$format_float = $workbook->add_format();
$format_float->set_num_format('0.00');
#
$format_hdr = $workbook->add_format(); # Add a format
$format_hdr->set_bold();
$format_hdr->set_border();
$format_hdr->set_align('center');
#
$format_hdr_merged = $workbook->add_format(); # Add a format
$format_hdr_merged->set_bold();
$format_hdr_merged->set_border();
$format_hdr_merged->set_align('center');
$format_hdr_merged->set_align('vcenter');


$line = 1;
$last_column  = 0;
$rvet_column  = 6;
$discr_column = 7;
$dets_column  = 8;
$pdbid_column  = 9;
$pdbaa_column  = 10;
$number_of_groups = 0;

while ( <RANKS_FILE> ) {
    next if ( !/\S/ );
    if ( /\%/ ){
	@aux = split;
	shift @aux;
	
	$rvet_column = 0;
	while ($rvet_column < @aux  && 
	       $aux[$rvet_column] ne "rvet"  && 
	       $aux[$rvet_column] ne "entropy" ) {$rvet_column++};
	($rvet_column == @aux) && die "No rvet in the output (?!).\n";
	$aux[$rvet_column] = "overall"; # the word "conservation" will appear above
	                                # it in the table
	
	$discr_column = 0;
	while ($discr_column < @aux  && $aux[$discr_column] !~ /discr/ ) {$discr_column++};
	($discr_column == @aux) && die "No discr in the output (?!).\n";

	$dets_column = 0;
	while ($dets_column < @aux  && $aux[$dets_column] !~ /dets_/ ) {$dets_column++};
	($dets_column == @aux) && die "No dets in the output (?!).\n";


	$pdbid_column = 0;
	while ($pdbid_column < @aux  && $aux[$pdbid_column] !~ /pdb_id/ ) {$pdbid_column++};
	($pdbid_column == $#aux+1 )  && ($pdbid_column = -1);

	$pdbaa_column = 0;
	while ($pdbaa_column < @aux  && $aux[$pdbaa_column] !~ /pdb_aa/ ) {$pdbaa_column++};
	($pdbaa_column == $#aux+1 )  && ($pdbaa_column = -1);

	$surf_column = 0;
	while ($surf_column < @aux  && $aux[$surf_column] !~ /surf/ ) {$surf_column++};
	($surf_column == $#aux+1 )  && ($surf_column = -1);

	$annot_column = 0;
	while ($annot_column < @aux  && $aux[$annot_column] !~ /annot/ ) {$annot_column++};
	($annot_column == $#aux+1 )  && ($annot_column = -1);


	# how many groups
	$number_of_groups =  scalar ( grep {/dets_/} @aux) ;
	


 	# unshift @aux, "cons";
        # header
	$last_column = @aux;

	$col_num = $last_column-1;
	
	$offset = int($col_num/26);
	$leftover = int($col_num%26);

	if(!$offset && $leftover){
	    $column = chr(ord('A') + $leftover);
	}
	else{
	   
	    $first_letter = chr(ord('A') + $offset -1);
	    #print $offset;
	    #print "++\n";
	    $second_letter = chr(ord('A') + $leftover);
	    $column = "$first_letter$second_letter";
	}
	#$column =  chr (ord('A') + $last_column-1);

        $worksheet->set_column("$column:$column", 20);
        # merge_range($first_row, $first_col, $last_row, $last_col, $token, $format, $utf_16_be)
	$first_col = 2;
	$last_col  = $first_col + $number_of_groups;
	$worksheet->merge_range(0, $first_col, 0, $last_col, "CONSERVATION", $format_hdr_merged);
	#
	$first_col = $last_col  + 1; # the previous last column
	$last_col  = $first_col + $number_of_groups;
	$worksheet->merge_range(0, $first_col, 0, $last_col, "SPECIFICITY", $format_hdr_merged);
	#
	$first_col = $last_col  + 1; # the previous last column
	$last_col  = $first_col + 2*$number_of_groups - 1;
	$worksheet->merge_range(0, $first_col, 0, $last_col, "REPRESENTATIVE SEQUENCES", $format_hdr_merged);
	#
	@sub_aux =  @aux[2..$#aux]; 
	$worksheet->write_row(1, 2, \@sub_aux, $format_hdr);
        #
	$first_row = 0;
	$last_row  = 1;
	$worksheet->merge_range($first_row, 0, $last_row, 0, $aux[0], $format_hdr_merged);
	$worksheet->merge_range($first_row, 1, $last_row, 1, $aux[1], $format_hdr_merged);
	#
	($pdbid_column > 0) &&
	    $worksheet->merge_range($first_row, $pdbid_column, $last_row, $pdbid_column, 
				$aux[$pdbid_column], $format_hdr_merged);
	#
	($pdbaa_column > 0) && $worksheet->merge_range($first_row, $pdbaa_column, $last_row, $pdbaa_column, 
				$aux[$pdbaa_column], $format_hdr_merged);

	($surf_column > 0)  && $worksheet->merge_range($first_row, $surf_column, $last_row, $surf_column, 
				"surface", $format_hdr_merged);

	($annot_column > 0) && $worksheet->merge_range($first_row, $annot_column, $last_row, $annot_column, 
				$aux[$annot_column], $format_hdr_merged);


	
	$line = 2;

   } else {

	chomp;
	@aux = split;

 	# the rest of the input file
	foreach $i (0 .. $#aux) {
	    $col_num = $i;
	   
	    
	    $offset = int($col_num/26);
	    $leftover = int($col_num%26);

	    if(!$offset && $leftover){
		$column = chr(ord('A') + $leftover);
	    }
	    else{

		
	    
		$first_letter = chr(ord('A') + $offset -1);
		$second_letter = chr(ord('A') + $leftover);
		$column = "$first_letter$second_letter";
	    }
	    #$column =  chr (ord('A') + $i);
	    if ($i == 0) {
		$worksheet->write("$column$line", $aux[$i], $format_right);
		if($aux[$i] =~ /\d+-\d+/){
		    
		    $worksheet->merge_range($line-1, 1, $line-1, $#aux, 
					    "insert in non-human sequences", $format_hdr_merged);
		    $last_used_column =  $#aux;
		    last;
		}
		    
	    } elsif( $i == 1) {
		$worksheet->write("$column$line", $aux[$i], $format_centered);  
  
	    } elsif( $i == $pdbid_column ) {
		$worksheet->write("$column$line", $aux[$i], $format_right);  

	    } elsif( $i == $pdbaa_column ) {
		$worksheet->write("$column$line", $aux[$i], $format_centered);  

	    } elsif( $i == $surf_column ) {
		if ($aux[$i] == 1) {
		    $worksheet->write("$column$line", "surface", $format_centered);  
		}

	    } elsif( $i == $annot_column ) {
		if ($aux[$i] eq  "none")  {
		    $aux[$i] = "  ";
		} else {
		    $aux[$i]  =~ s/_/ /g;
		}
		$worksheet->write("$column$line", $aux[$i], $format_left);  

	    } elsif( $i >= $discr_column &&  $i <= $discr_column + $number_of_groups ) {
		$cvg =  $aux[$i];
		if ( $cvg <= 0.5) {
		    $color_index = int ( (0.5-$cvg)*$COLOR_RANGE );
		    $format = $workbook->add_format(
			bg_color => $orange_range[$color_index],
			pattern  => 1,
			border   => 1
			);
		} else {
		    $color_index = int ( ($cvg-0.5)*$COLOR_RANGE );
		    $format = $workbook->add_format(
			bg_color => $blue_range[$color_index],
			pattern  => 1,
			border   => 1
			);
		    
		}
		$format ->set_num_format('0.00');

		$offset = int($col_num/26);
		$leftover = int($col_num%26);

		if(!$offset && $leftover){
		    $column = chr(ord('A') + $leftover);
		}
		else{
	    
		    $first_letter = chr(ord('A') + $offset -1);
		    $second_letter = chr(ord('A') + $leftover);
		    $column = "$first_letter$second_letter";
		}

		$worksheet->write("$column$line", $aux[$i], $format);
		

	    } elsif( $i >= $rvet_column  &&  $i <= $rvet_column + $number_of_groups ) {
		# color the  columns representing conservation
		if ( $aux[$i] =~ /[A-Z]/) {
		    print " $i   $aux[$i] \n";
		    exit;
		}
		$cvg =  $aux[$i];
		$color_index = int ($cvg*$COLOR_RANGE );
		$format = $workbook->add_format(
		    bg_color => $xl_color[$color_index],
		    pattern  => 1,
		    border   => 1
		    );
		$format ->set_num_format('0.00');

		$offset = int($col_num/26);
		$leftover = int($col_num%26);

		if(!$offset && $leftover){
		    $column = chr(ord('A') + $leftover);
		}
		else{
	    
		    $first_letter = chr(ord('A') + $offset -1);
		    $second_letter = chr(ord('A') + $leftover);
		    $column = "$first_letter$second_letter";
		}
		
		#$column = chr (ord('A') + $i);

		$worksheet->write("$column$line", $aux[$i], $format);

	    } else {
		$worksheet->write("$column$line", $aux[$i], $format_centered);    
	    }
	}
	$last_used_column =  $#aux;


   }

    $line++;

}


##################################################
##################################################
##################################################
##################################################
# colorbar

$column =  $last_used_column + 2;
$row = 2;
$worksheet->write($row, $column, "conserved", $format_right);   
$row = 2+$COLOR_RANGE;
$worksheet->write($row, $column, "variable", $format_right);   



$column =  $last_used_column + 3;
$row = 0;

$worksheet->write($row, $column, "cons colorbar", $format_hdr);   

for ($ctr = 0; $ctr <= $COLOR_RANGE; $ctr++) {
    $format = $workbook->add_format(
	bg_color => $xl_color[$ctr],
	pattern  => 1,
	border   => 1
	);
    $row = $ctr+2;
    $worksheet->write($row, $column, '',  $format);
}
##################################################

$column =  $last_used_column + 7;
$row = 2;
$worksheet->write($row, $column, "specific", $format_left);   
$row = 2+$COLOR_RANGE;
$worksheet->write($row, $column, "non-specific", $format_left);   

$column =  $last_used_column + 6;
$row = 0;
$worksheet->write( $row, $column,"spec colorbar", $format_hdr);   

for ($ctr = 0; $ctr <= $COLOR_RANGE/2; $ctr++) {
    $format = $workbook->add_format(
	bg_color => $orange_range[$ctr],
	pattern  => 1,
	border   => 1
	);
    $row = 2+$COLOR_RANGE/2-$ctr;
    $worksheet->write( $row, $column, '', $format);
}

$row = 1+$COLOR_RANGE/2;


for ($ctr = 0; $ctr <= $COLOR_RANGE/2; $ctr++) {
    $format = $workbook->add_format(
	bg_color => $blue_range[$ctr],
	pattern  => 1,
	border   => 1
	);
    $row ++;
    $worksheet->write($row, $column, '', $format);
}



close RANKS_FILE;
