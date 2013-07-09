#! /usr/bin/perl
#

use strict;
use warnings;
use diagnostics;


open(FILE, "src.txt") || die "Cannot open $!\n";

my @text = <FILE>;

my $scer = qr/Scer\s{12}([atgcATGC-]{1,})/;


foreach my $line (@text) {
	match($scer,$line);
}


close FILE;


sub match {
	
	my $re = $_[0];
	my $text = $_[1];

	if ( $text =~ $re ) {
		print "Match\n";
	} else {
		print "No match.\n";
	}
}


