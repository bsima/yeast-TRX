#! /usr/bin/perl


use strict;
use warnings;
use diagnostics;

#use File::Copy;

# Var matey!
my @dirFiles = glob "./src2/*";
my $scer = qr/Scer\s{12}([atgcATGC-]{1,})/;
my $spar = qr/Spar\s{12}([atgcATGC-]{1,})/;
my $smik = qr/Smik\s{12}([atgcATGC-]{1,})/;
my $sbay = qr/Sbay\s{12}([atgcATGC-]{1,})/;
my $crickRegex = qr/.{1,}C.aln/;

my %Saccharomyces = (
	'cerevisiae' => qr/Scer\s{12}([atgcATGC-]{1,})/,
	'paradoxus'  => qr/Spar\s{12}([atgcATGC-]{1,})/,
    'martinae'   => qr/Smik\s{12}([atgcATGC-]{1,})/,
	'bayanus'    => qr/Sbay\s{12}([atgcATGC-]{1,})/
);	

# Start building the R script
open(RSCRIPT, ">>rscript.r");

# First create a NULL vector for each species to store the data
for (keys %Saccharomyces) {
	print RSCRIPT $_ . "<-c()\n";
}

foreach my $fileName (@dirFiles) {

    # First, open each file.	
	open(FILE, $fileName) || die "Cannot open file: $!";
	
	my @text = <FILE>;

	# For each line of the file, extract the genetic code (using my awesome match sub)
	foreach my $line (@text) {

		# Lets cycle through each species
		for (keys %Saccharomyces) {
			my $regex = $Saccharomyces{$_};

			if ( match($regex,$line) ) {
				# If we've found a match, extract the data and assign it to the $line variable
				$line = match($regex,$line);
			
				# Don't forget to reverseCompliment the Crick files!
				if ( $fileName =~ $crickRegex ) {
					$line = reverseCompliment($line);
				}

				# Now let's do something with the raw data...
				# --------------------------------------------
				# We need to separate each species. Then we can input the raw data from
				# each into an R data frame and run the TRX calculation separately on
				# each species.
				#
				# I think the way to do this is to generate a temporary R script for each
				# respective species that passess the $species variable to the script as
				# the key for the data. The genetic data should all be stored together.
				#
				# First I will just pass it to R and calculate the average, then figure
				# out the TRX calculation later.
	
				# First create a NULL vector for each species to store the data
				# Now append the data to each specie's respective vectors 
				print RSCRIPT $_ . "<-append(" . $_ . ",\"" . $line . "\")\n";			
			} 
		}
	}
	close FILE;
};

close RSCRIPT;


# @name match
# @description Matches a given text to a regular expression. First pass the regular expression, then the text as params. Returns the matched text.
sub match {

	my $re = $_[0];
	my $text = $_[1];

	if ( $text =~ $re ) {
		return $1;
	} else {
		return $1;
	}
}

# @name reverseCompliment
# @description Given a string of text of coded DNA, this outputs the reverse compliment of the strand.
sub reverseCompliment {

	# Get DNA to work on
	my ($dna) = @_;

	# First reverse the DNA strand
	$dna = reverse $dna;

	# Now translate the DNA
	$dna =~ tr/ACGTacgt-/TGCAtgca-/;

	# Return the output
	return $dna;
}
