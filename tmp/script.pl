#! /usr/bin/perl


use strict;
use warnings;
use diagnostics;

use File::Copy;

my @dirFiles = glob "./src2/*";

foreach my $fileName (@dirFiles) {
		
    # First, open each file.	
	open(FILE, $fileName) || die "Cannot open file: $!";
	
	# For each line of the file, extract the genetic code for Scer	
	while ( my $line = <FILE> ) {	
	
		# Extract the genetic code using my super awesome match sub	
		my $scer = qr/Scer\s{12}([atgcATGC-]{1,})/;
		$line = match($scer,$line);
		print "Match\n";
		print $line . "\n";

		# If it's a Crick file, we have to reverseCompliment it
		if ( $fileName =~ /.{1,}(C.aln)/ ) {
			$line = reverseCompliment($line);
			print $fileName . " is a Crick File! Line reversed\n";
			#print $line . "\n";
		};
		
	};
		
	close FILE;
};

# @name match
# @description Matches a given text to a regular expression. First pass the regular expression, then the text as params. Returns the matched text.
sub match {

	my $re = $_[0];
	my $text = $_[1];
	
	# Let's do some error checking
	# if ( defined($re) ) { print "Regular expression is " . $re . "\n"; } else { print "Regular expression not defined.\n"; }
	# if ( defined($text) ) { print "Text is " . $text . "\n"; } else { print "Text is not defined.\n"; }	

	if ( $text =~ $re ) {
		return $1;
	} else {
		return $1;
		print "No match.\n";
	}
}

# @name reverseCompliment
# @description Given a string of text of coded DNA, this outputs the reverse compliment of the strand.
# @TODO Right now this can only handle one string at a time. I wonder if it would be worth it to make it handle various types of input...?
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
