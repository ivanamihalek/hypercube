#!/usr/bin/perl -w

while (<>) {
    if ( /^>/ ) {
	/>(\S+)/;
	print ">$1\_THRB\n";
    } else {
	print;
    }
}
