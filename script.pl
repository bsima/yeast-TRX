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
my @yeastGenome = glob "./YeastGenome-tmp/*";
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
my $geneRe  = qr/\/([\d|\w]+)\.[\d|\w]+/;



print "\nInitializing...\n\n";

# This creates a separate CSV file for each species
# and prints a header line with the species name
foreach (keys %Saccharomyces) {
	open(FILE, ">", "$_.csv");
	print FILE "Saccharomyces $_\n";
	print FILE "gene,code\n";
	close FILE;
}

print "Currently writing data\n";
print "----------------------\n\n";

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
				
				# Get the specific gene name from the $fileName
				my $geneName = match($geneRe,$fileName);

				# Report the data to the respective CSV file in the
				# following format:
				# 		gene,code	
				open(FILE, ">>$species.csv");
				print FILE $geneName . "," . $line . "\n";
				close FILE;

				# Not done yet... We have to run the TRX matches now
				# We cycle through each $sequence in the %trxSequences hash.
				# If a match is found, we print it to TRXDATA in a format ready
				# to run as an R script that creates a data frame.
				#foreach my $sequence (keys %trxSequences) {
				#	if ( match($trxSequences{$sequence},$line) ) {
				#		my $trxLine = match($trxSequences{$sequence},$line);
				#		print TRXDATA $species . "<-append(" . $sequence . ",\"" . $trxLine . "\")\n";
				#	}
				#}
			} 
		}
	}
	close FILE;
};

print "\nData has been written to CSV files.\n"; 

foreach (keys %Saccharomyces) {
	
	# First, open the respective file
	# Second, note the position on the genome -> $position
	# 		Probably use a while loop and an increasing $i var
	# 		Just make sure to record this into the resultant file with the TRX value
	# Third, calculate the TRX value and write to a new file
	# Fourth, move on to the next position
	# Finally, close the file
	open(FILE,"<$_.csv");

	my @text = <FILE>;

	foreach my $line (@text) {
		
	}	

};

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

# @name trxScore
# @description Calculates the TRX value of a given phosphate linkage
sub trxScore {

	my %trxScores = (
		'CpG' => 43,
		'CpA' => 42,
		'TpG' => 42,
		'GpG' => 42,
		'CpC' => 42,
		'GpC' => 25,
		'GpA' => 22,
		'TpC' => 22,
		'TpA' => 14,
		'ApG' =>  9,
		'CpT' =>  9,
		'ApA' =>  5,
		'TpT' =>  5,
		'ApC' =>  4,
		'GpT' =>  4,
		'ApT' =>  0
	);

	


}
