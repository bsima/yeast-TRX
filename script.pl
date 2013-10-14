#! /usr/bin/perl
# 
use strict;
use warnings;
use diagnostics;
use Statistics::Descriptive;


# Var matey!
my @yeastGenome = glob "./data/YeastGenome-tmp/*";
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
my %deltaESequences = (
    'rA/rA' => qr/(AA)/,
	'rA/rC' => qr/(AC)/,
	'rA/rG' => qr/(AG)/,
	'rA/rU' => qr/(AU)/,
	'rC/rA' => qr/(CA)/,
    'rC/rC' => qr/(CC)/,
	'rC/rG' => qr/(CG)/,
	'rC/rU' => qr/(CU)/,
	'rG/rA' => qr/(GA)/,
	'rG/rC' => qr/(GC)/,
    'rG/rG' => qr/(GG)/,
    'rG/rU' => qr/(GU)/,
    'rU/rA' => qr/(UA)/,
    'rU/rC' => qr/(UC)/,
    'rU/rG' => qr/(UG)/,
    'rU/rU' => qr/(UU)/
);
my $geneNameRe  = qr/([\d\w]{5,})[.,][\d|\w]+/;
my $geneRe      = qr/,([atgcATGC-]+)/;
my $fileTracker = 0;

print "\nInitializing...\n\n";

# Setup the directories for each species and
# the data files
foreach my $species (keys %Saccharomyces) {
	
    mkdir("./data/$species");
    
    open(FILE, ">", "./data/$species/genome.csv");
	print FILE "gene,code\n";
	close FILE; 
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
    my $geneName = match($geneNameRe,$fileName);

    # Print to each of the species' CSV files the gene name we are currently analyzing
    printToSpecies("$geneName,");

    # Setup the data files for writing the TRX and delta-E values later.
    # We have to do this now while we have the `$geneName` varialbe available.
    foreach my $species (keys %Saccharomyces) {
        
        mkdir("./data/$species/$geneName");
        
        open(my $raw,">./data/$species/$geneName/raw.csv");
            print $raw "gene,dinucleotide,position,trx.score,energy.score\n";
        close $raw; 
        open(my $smooth,">./data/$species/$geneName/smooth.csv");
            print $smooth "gene,position,trx.score,energy.score\n";
        close $smooth;
    }

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
                open(SPECIES, ">>./data/$species/genome.csv") || die "Cannot open file: $!\n";
                print SPECIES $line; 
                close SPECIES;
		    } 
		}
	}
  
    # End the line of genetic code for this gene with a
    # line break character
    printToSpecies("\n");

	close GENEFILE;

    $fileTracker++;
    print "File number: $fileTracker\n";
};

print "\nGenomes have been written to CSV files.\nNow calculating TRX and delta-E Scores.\n";

# Calculate the TRX and delta-E scores and then input them to the 
# respective files in the `$species/$geneName` directory.
#
#
# Let's loop through each `$species/genome.csv` and read them straigt away.
foreach my $species (keys %Saccharomyces) {

    open(SPECIES,"<./data/$species/genome.csv") || die "Cannot open file: $!";
    my @text = <SPECIES>;    

    print "Working on $species right now.\n";

    # Each line of the `$species/genome.csv` contains the entire code of one gene.
    # We will loop through each line (gene) and calculate two things: the raw TRX
    # and delta-E scores, and the smoothed TRX and delta-E.
    foreach my $line (@text) {
    
        # Only run the calculation is we have genetic code.
        # Also, we strip out the prepended gene name in the file.
        if ( defined $line && $line =~ qr/.+,[acgtACGT-]+/ ) {
       
            my $geneName = match($geneNameRe,$line);
            my $gene = match($geneRe,$line);

            # Raw data
            # 
            # This is the original formulation that calculates the trxValue and
            # energyScore on a per-dinucleotide basis. It then prints the data 
            # to the respective CSV file in the following format:
            #       gene,dinucleotide,position,trx.score,energy.score 
            open(RAW, ">>./data/$species/$geneName/raw.csv") || die "Cannot open file: $!";
            for ( my $position = 0; $position < length($gene); $position++ ) {
                my $dinucleotide = substr($gene,$position,2); 
            
                if ( $dinucleotide =~ qr/[actgACTG]{2}/ ) {     
                    my $trxValue = trxScore($dinucleotide);
                    my $energyScore = energyScore($dinucleotide);
                   
                    print RAW $geneName . "," . $dinucleotide . "," . $position . "," . $trxValue . "," . $energyScore . "\n";
                }
            }
            close RAW;
           
            # Smoothing function
            #
            # This moves through the gene data and counts the position until it arrives
            # at the end of the smoothing window (e.g. 200). Then it calculates the average
            # of the selected data set and outputs it into the respecive `smooth.csv` file 
            # in the following format, to be read later by R's graphing functions:
            #       gene,position,trx.mean,energy.mean
            open(SMOOTH, ">>./data/$species/$geneName/smooth.csv") || die "Cannot open file $!";
            my $start = 0;
            my $smoothingWindow = 200; 
            for ( my $position = $start; $position < length($gene); $position++ ) {
                
                for ( my $smoothing = 0; $smoothing < $smoothingWindow; $smoothing++ ) {
                
                    my @trxValues;
                    my @energyScores;

                    my $dinucleotide = substr($gene,$position,2); 
            
                    if ( $dinucleotide =~ qr/[actgACTG]{2}/ ) {     
                        my $trxValue = trxScore($dinucleotide);
                        my $energyScore = energyScore($dinucleotide);
                       
                        # print $trxValue; 
                        push @trxValues,$trxValue;
                        push @energyScores,$energyScore;
                    }
                   
                    #print "trxValues = @trxValues \n";

                    my $trxStat = Statistics::Descriptive::Full->new();
                    $trxStat->add_data($trxValues[$smoothing]);
                    #print "trxStat = $trxStat \n";
                    my $trxMean = $trxStat->mean();

                    my $energyStat = Statistics::Descriptive::Full->new();
                    $energyStat->add_data($energyScores[$smoothing]);
                    my $energyMean = $energyStat->mean();
                    
                    print SMOOTH $geneName . "," . $position . "," . $trxMean . "," . $energyMean . "\n";

                }
            }
            close SMOOTH;
        }
    } 
}


# Initiate R scripts
#
print "Initating R scripts";
open(R, "| R --vanilla") || die "Could not start R command line\n";
print R 'source("plot.R")';
close R;

print "\n\nAll Done!\nNow enjoy a cup of tea.\n\n";

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

# @name energyScore
# @description Calculates and returns the delta-E value of a given phosphate linkage
# @param $dinucleotide {string} The nucleotide to be checked
# @return {integer} The delta-E value
sub energyScore {
   
    my ( $dinucleotide ) = @_;

	my %energyScores = (
		qr/(AA)/ => -18.5,
		qr/(AC)/ => -19.0,
		qr/(AG)/ => -23.6,
		qr/(AU)/ => -15.7,
        qr/(CA)/ => -20.0,
		qr/(CC)/ => -21.4,
		qr/(CG)/ => -26.9,
		qr/(CU)/ => -17.2,
        qr/(GA)/ => -23.7,
        qr/(GC)/ => -22.9,
        qr/(GG)/ => -24.3,
		qr/(GU)/ => -18.9,
		qr/(UA)/ => -19.6,
		qr/(UC)/ => -28.2,
		qr/(UG)/ => -23.3,
		qr/(UU)/ => -15.8
	);

    foreach my $re (keys %energyScores) {
        if ( match($re,$dinucleotide) ) {
            return $energyScores{$re};
        } 
    }
    return 0;
}

# @name printToSpecies
# @description Loops through the $species/genome.csv files and writes something to them
# @param $str {string} Text to be written
sub printToSpecies {

    my ( $str ) = @_;

	foreach my $species (keys %Saccharomyces) {
        # Get the specific gene name from the $fileName
        open(SPECIES, ">>./data/$species/genome.csv") || die "Cannot open file: $!\n";
        print SPECIES $str; 
        close SPECIES;
    }
    return 1;
}
