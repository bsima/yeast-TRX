#! /usr/bin/perl
# 
# @TODO Modify $trxData to create the proper data frame that we need
# @TODO Put $trxData into a tmp folder for deletion after it's run in R
# @TODO Write the trx.r script

use strict;
use warnings;
use diagnostics;

#use File::Copy;

# Var matey!
my @yeastGenome = glob "./YeastGenome/*";
my $isItCrick   = qr/.{1,}C.aln/;
my %Saccharomyces = (
	'cerevisiae' => qr/Scer\s{12}([atgcATGC-]{1,})/,
	'paradoxus'  => qr/Spar\s{12}([atgcATGC-]{1,})/,
    'martinae'   => qr/Smik\s{12}([atgcATGC-]{1,})/,
	'bayanus'    => qr/Sbay\s{12}([atgcATGC-]{1,})/
);	
my %trxSequences = (
	'CpGCpG' => qr/(CG.{1,}CG)/,
	'CpATpG' => qr/(CA.{1,}TG)/,
	'GpGCpC' => qr/(GG.{1,}CC)/,
	'GpCGpC' => qr/(GC.{1,}GC)/,
	'GpATpC' => qr/(GA.{1,}TC)/,
	'TpATpA' => qr/(TA.{1,}TA)/,
	'ApGCpT' => qr/(AG.{1,}CT)/,
	'ApATpT' => qr/(AA.{1,}TT)/,
	'ApCGpT' => qr/(AC.{1,}GT)/,
	'ApTApT' => qr/(AT.{1,}AT)/
);
my $trxData   = "trxData.r";

print "\nInitializing...\n\n";

# Start building the R script to make the TRX data frame
open(TRXDATA, ">>$trxData") || die "Cannot open file: $!";

# First create a NULL vector for each species so we
# can append the data later
foreach (keys %Saccharomyces) {
	print TRXDATA $_ . "<-c()\n";
}

print "Currently writing data to $trxData\n";
print "------------------------------------\n\n";

# We're gonna have to run the operation for each
# file in the YeastGenome
foreach my $fileName (@yeastGenome) {

	# First, open each file.	
	open(FILE, $fileName) || die "Cannot open file: $!";

	# Then, load the text of the file into an array	
	my @text = <FILE>;

	# For each line of the file, extract the genetic code (using my awesome match sub)
	foreach my $line (@text) {

		# Lets cycle through each species
		foreach my $species (keys %Saccharomyces) {

			# Load the regex for the species from the hash into $regex
			my $regex   = $Saccharomyces{$species};

			if ( match($regex,$line) ) {
				# If the line from the data file matches the species we 
				# are currently looking for, then extract the data and 
				# assign it to the $line variable
				$line = match($regex,$line);
			
				# Don't forget to reverseCompliment the Crick files!
				if ( $fileName =~ $isItCrick ) {
					$line = reverseCompliment($line);
				}

				# Not done yet... We have to run the TRX matches now
				# We cycle through each $sequence in the %trxSequences hash.
				# If a match is found, we print it to TRXDATA in a format ready
				# to run as an R script that creates a data frame.
				foreach my $sequence (keys %trxSequences) {
					if ( match($trxSequences{$sequence},$line) ) {
						my $trxLine = match($trxSequences{$sequence},$line);
						print TRXDATA $species . "<-append(" . $sequence . ",\"" . $trxLine . "\")\n";
					}
				}
			} 
		}
	}
	close FILE;
};

# Finally, combine all of the species vectors into one data frame
# I don't know how to do this... 
close TRXDATA;
print "\nData has been written to $trxData.\n"; 

# ##################################################################################################
# Now we have a file, $trxData, that, when run as an R script via source() function, will (sort of) 
# create a data frame of all the matched TRX sequences. With some editing (i.e. better R scripting),
# we could easily pair the matched TRX sequence with the species. Then, an R script could calculate
# the TRX values for each sequence in each species. Said R script will be trx.r
# ##################################################################################################

# @name match
# @description Matches a given text to a regular expression. First pass the regular expression, then the text as params. Returns the matched text.
sub match {

	my ( $re, $text ) = @_;

	if ( $text =~ $re ) {
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
