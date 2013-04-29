#!/usr/bin/perl 
@ARGV || die "Usage: $0  <score_file>  <alm_msf> \n";

my ($score_file,  $alm_msf) = @ARGV;
my $score_annot_file = $score_file . ".annotated";
my $tmp_file = $score_file . ".tmp";

my %seq_annotation = ();
$seq_annotation{"MAC_EUG_LALBA"}[56] = "disord  region";
$seq_annotation{"MAC_EUG_LALBA"}[68] = "disord  region";

my @annotated_seqs = keys %seq_annotation; 



############################################
# read in the alignment
open(MSF, "<$alm_msf")||html_die("Cno:$alm_msf");
while(<MSF>){
    last if(/\/\//);
}
my %sequence = ();
my @names = ();
do {
   if ( /\w/ ) {
       my @aux = split;
       my $name = $aux[0];
       my $aux_str = join ('', @aux[1 .. $#aux] );

       if ( defined $sequence{$name} ) {
	   $sequence{$name} .= $aux_str;
       } else {
	   push @names, $name;
	   $sequence{$name}  = $aux_str;
       }
   } 
} while ( <MSF>);


############################################
# turn the msf into a table (first index= sequence, 2nd index= position
my $seq = 0;
my (%array,%counter,%seqno);
   


# recalculate the position numbers for each sequence
foreach my $name ( @annotated_seqs) {
    
    my @aux = split '', $sequence{$name};

    foreach my $pos ( 0 .. $#aux ) {
	$array{$name}[$pos] = $aux[$pos];
	if ( $aux[$pos] !~ /[\.\-]/ ) {
	    (defined  $counter{$name}) || ($counter{$name} = 0);
	    # take care of the fact that people count from 1
	    # by incrementing first
	    $counter{$name} ++; 
	    $seqno{$name}[$pos] = $counter{$name};
	} else {
	    $seqno{$name}[$pos] = "-";
	}
    }
    $seq++;
}

    ##############################################
    #  output
    open(OFH,">$tmp_file")||return("Cno:$tmp_file");
    print OFH " annotation \n"; # header

    # [ZH: ] here random chose $annot_name[0],since all seq are the same length with gappes
    my $alignment_length = length $sequence{$annotated_seqs[0]};


    for my $pos  (0 .. $alignment_length-1) { 
	    
	my $annot_str = "";
	foreach  my $name (keys %seq_annotation){
	    my $pos_map  = $seqno{$name}[$pos];
	    if (defined $seq_annotation{$name}[$pos_map]) {
		$annot_str && ($annot_str.="; ");
		$annot_str .= $seq_annotation{$name}[$pos_map];
	    }
	}
	$annot_str || ($annot_str="none");
	$annot_str  =~ s/\s/_/g;
	$annot_str .= "\n";

	print OFH $annot_str;
    }
    close (OFH);


    `pr -m -t -s $score_file $tmp_file > $score_annot_file`;
    `rm $tmp_file`;

 
    print  $score_annot_file."\n";
