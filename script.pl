#! /usr/bin/perl
# 
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
my $trxData    = "trxData.r";
my $geneRe     = qr/([\d|\w]+)[.|,][\d|\w]+/;
my $lineNumber = 0;

print "\nInitializing...\n\n";

# This creates a separate CSV file for each species
# and prints a header line with the species name
foreach (keys %Saccharomyces) {
	open(FILE, ">", "$_.csv");
	print FILE "gene,code\n";
	close FILE;

    open(TRXSCORE,">$_-TRXscore.csv");
    print TRXSCORE "gene,dinucleotide,position,trx score\n";
    close TRXSCORE;
}

print "Currently writing data\n";
print "----------------------\n\n";

# We're gonna have to run the operation for each
# file in the YeastGenome
foreach my $fileName (@yeastGenome) {

	# First, open each file.	
	open(GENEFILE, $fileName) || die "Cannot open file: $!";

	# Then, load the text of the file into an array	
	my @text = <GENEFILE>;
    my $geneName = match($geneRe,$fileName);

    # Print to each of the species' CSV files the gene name we are currently analyzing
    printToSpecies("$geneName,");

	# For each line of the gene file, extract the genetic code (using my awesome match sub)
	foreach my $line (@text) {

        # Lets cycle through each species
		foreach my $species (keys %Saccharomyces) {

			# Load the regex for the species from the hash into $regex
			my $regex = $Saccharomyces{$species};

			if ( match($regex,$line) ) {
				# If the line from the data file matches the species we 
				# are currently looking for, then extract the data and 
				# assign it to the $line variable
				$line = match($regex,$line);
			
				# Don't forget to reverseCompliment the Crick files!
				if ( $fileName =~ $isItCrick ) {
					$line = reverseCompliment($line);
				}
        
				# Append the data to the species CSV file 
                open(SPECIES, ">>$species.csv") || die "Cannot open file: $!\n";
                print SPECIES $line; 
                close SPECIES;
		    } 
		}
	}
  
    # End the line of genetic code for this gene with a
    # line break character
    printToSpecies("\n");

	close GENEFILE;
};

print "\nData has been written to CSV files.\nNow calculating TRX Score.";

# Now I have to loop through each species CSV file and read them straight away
# I will have to track my position, but not line number.
#

foreach my $species (keys %Saccharomyces) {

    open(SPECIES,"<$species.csv") || die "Cannot open file: $!";
    my @text = <SPECIES>;
   
    foreach my $line (@text) {
        
        if ( $line =~ qr/.+,[acgtACGT]+/ ) {
       
            my $geneName = match($geneRe,$line);
       
            # Print the data to the respective CSV file in the
            # following format:
            #       gene,dinucleotide,position,trx score
            open(TRXSCORE, ">>$species-TRXscore.csv") || die "Cannot open file: $!";
        
            for ( my $position = 0; $position < length($line); $position++ ) {
            
                my $dinucleotide = substr($line,$position,2); 
            
                if ( $dinucleotide =~ qr/[actgACTG]{2}/ ) {     
                    my $trxValue = trxScore($dinucleotide);
                    print TRXSCORE $geneName . "," . $dinucleotide . "," . $position . "," . $trxValue . "\n";
                }
            }
            close TRXSCORE;
        }
    }
    close SPECIES;
}


# @name match
# @description Matches a given text to a regular expression. First pass the regular expression, then the text as params. Returns the matched text.
# @param $re {string} Regular expression
# @param $text {string} Text to be run against the regular expression
# @return {string} The matched text
sub match {

	my ( $re, $text ) = @_;

	if ( $text =~ $re ) {
		return $1;
	}
}

# @name reverseCompliment
# @description Given a string of text of coded DNA, this outputs the reverse compliment of the strand.
# @param $dna {string} DNA string to be reversed
# @return {string} Reversed DNA
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
# @description Calculates and returns the TRX value of a given phosphate linkage
# @param $dinucleotide {string} The nucleotide to be checked
# @return {integer} The TRX value
sub trxScore {
   
    my ( $dinucleotide ) = @_;

    # @TODO Find the paper where I got these scores from... Can't remember 
	my %trxScores = (
		qr/(CG)/ => 43,
		qr/(CA)/ => 42,
		qr/(TG)/ => 42,
		qr/(GG)/ => 42,
		qr/(CC)/ => 42,
		qr/(GC)/ => 25,
		qr/(GA)/ => 22,
		qr/(TC)/ => 22,
		qr/(TA)/ => 14,
		qr/(AG)/ =>  9,
		qr/(CT)/ =>  9,
		qr/(AA)/ =>  5,
		qr/(TT)/ =>  5,
		qr/(AC)/ =>  4,
		qr/(GT)/ =>  4,
		qr/(AT)/ =>  0
	);

    foreach my $re (keys %trxScores) {
        if ( match($re,$dinucleotide) ) {
            return $trxScores{$re};
        } 
    }
    return 0;
}

# @name printToSpecies
# @description Loops through the $species.csv files and writes something to them
# @param $str {string} Text to be written
sub printToSpecies {

    my ( $str ) = @_;

	foreach my $species (keys %Saccharomyces) {
        # Get the specific gene name from the $fileName
        open(SPECIES, ">>$species.csv") || die "Cannot open file: $!\n";
        print SPECIES $str; 
        close SPECIES;
    }
    return 1;
}
# @name position
# @description A loop that updates and returns the current position on the genome
sub position {
    
}
