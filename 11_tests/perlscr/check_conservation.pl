#! /usr/bin/perl -w

@ARGV ||
    die "Usage:  $0  <tp list>  <score>\n";

foreach (@ARGV) {
    (-e $_) || die "$_ not found.\n";
}
($tp, $score) = @ARGV;

@lines = split "\n", `cat $tp`;

foreach  (@lines) {
    next if ( ! /\S/);
    chomp;
    @aux = split;
    $ret =  `awk \' \$1 == $aux[1] \' $score`;
    print $ret;
}
