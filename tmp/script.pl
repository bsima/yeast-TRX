#! /usr/bin/perl
#

use strict;
use warnings;
use diagnostics;

use File::Copy;

my @dirFiles = glob "./src/*";

# my $findScer = /^(Scer)/;

foreach my $fileName (@dirFiles) {
		open(FILE, $fileName) || die "Cannot open file: $!";
		
		while ( my $line = <FILE> ) {
				
				if ( $line =~ m/^Scer\s{12}([atgcATGC]{1,})/ ) {
						print $1 . "\n";
				};
		};
		close FILE;
};
