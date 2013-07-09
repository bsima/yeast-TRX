#! /usr/bin/perl
#

use strict;
use warnings;
use diagnostics;


open(FILE, "src.txt") || die "Cannot open $!\n";

my @text = <FILE>;

my $scer = qr/Scer\s{12}([atgcATGC-]{1,})/;
 
foreach my $line (@text) {

	if ( $line =~ $scer) {
		print "Match" . "\n";
	} else {
		print "No match.\n";
	}
}


close FILE;
