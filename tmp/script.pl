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
			
				# If it's a Crick file, reverse it
				if ( $fileName =~ /.{1,}(C.aln)/ ) {
						$line = reverseCompliment($line);
						print "Line reversed\n";
						print $line . "\n";
				};

				# Print the matched line
				if ( $line =~ /^Scer\s{12}([atgcATGC]{1,})/ ) {
						print $1 . "\n";
				};
		};
		close FILE;
};


sub reverseCompliment {

		# Get DNA to work on
		my ($dna) = @_;

		# First reverse the DNA strand
		$dna = reverse $dna;

		# Now translate the DNA
		$dna =~ tr/ACGTacgt/TGCAtgca/;

		# Return the output
		return $dna;
}
