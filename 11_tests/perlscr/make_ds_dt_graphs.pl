#! /usr/bin/perl -w

$file = "all_methods.auc";

(-e $file) || die "$file not found.\n";

@lines = split "\n", `grep _dt $file`;
$name2 = "";

(-e "ds_dt.table") && `rm ds_dt.table`;

$out = "dt_ds.table";
(open OF, ">$out") || die "Cno $out: $!.\n";
foreach $line (@lines) {
    ($name, $val) = split " ", $line;
    $name =~ s/_dt//;
    $ret = `grep $name\_ds $file`;
    ($name2, $val2) = split " ", $ret;
    
    print  OF "$name   $val  $val2\n";
}
$ret = `grep CAP  $file`;
chomp $ret;
($name, $capval) = split " ", $ret;

$ret = `grep MI  $file`;
chomp $ret;
($name, $mival) = split " ", $ret;

$ret = `grep SPR  $file`;
chomp $ret;
($name, $sprval) = split " ", $ret;


close OF;


################################################################
$infile = $out;
$out = "ds_dt.gscr";
(open OF, ">$out") || die "Cno $out: $!.\n";
print OF 
"
unset key
set style data lines 

set xlabel \"Area under roc\"
set ylabel \"Method\"

set term  post color
set output \"tmp.ps\"

set size ratio 0.3
set xtics nomirror rotate by -45
mi(x) = $mival
cap(x) = $capval
spr(x) = $sprval
plot \'$infile\' using 2,  \'$infile\' using 3:xticlabels(1), mi(x), cap(x), spr(x)


";


`gnuplot ds_dt.gscr`;
`ps2pdf tmp.ps auc.pdf && rm tmp.ps `;

################################################################


$file = "dt_ds.table";

($name, $ds, $dt) = ();

@lines = split "\n", `cat $file`;
@cs = ();
foreach $line (@lines) {
    ($name, $dt, $ds) = split " ", $line;
    $cons_spec = substr $name, 0, 2;
    $dist      = substr $name, 2, 1;
    if ( ! defined $seen{$cons_spec} ) {
	push @cs, $cons_spec;
	$seen{$cons_spec} = 1;
    } 
    $score{$cons_spec}{$dist} = $dt;
}

$out = "dist_comp.table";
(open OF, ">$out") || die "Cno $out: $!.\n";

foreach $cons_spec (@cs) {
    print OF " $cons_spec   ", $score{$cons_spec}{"l"}, "  ", $score{$cons_spec}{"e"}, "\n";
}
close OF;


################################################################
$infile = $out;

$out = "ds_dt.gscr";
(open OF, ">$out") || die "Cno $out: $!.\n";
print OF 
"
unset key
set style data lines 

set xlabel \"Area under roc\"
set ylabel \"Method\"

set term  post color
set output \"tmp.ps\"

set size ratio 0.3
set xtics nomirror rotate by -45

plot \'$infile\' using 2,  \'$infile\' using 3:xticlabels(1)


";


`gnuplot ds_dt.gscr`;
`ps2pdf tmp.ps auc.dist.pdf && rm tmp.ps `;
