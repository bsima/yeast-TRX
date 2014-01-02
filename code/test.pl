#!/usr/bin/perl
use v5.14;
use strict;
use warnings;
use diagnostics;
use feature "switch";

my %Saccharomyces = (
	'cerevisiae' => qr/Scer\s{12}([atgcATGC-]{1,})/,
	'paradoxus'  => qr/Spar\s{12}([atgcATGC-]{1,})/,
    'martinae'   => qr/Smik\s{12}([atgcATGC-]{1,})/,
	'bayanus'    => qr/Sbay\s{12}([atgcATGC-]{1,})/
);
my $geneNameRe  = qr/(eY\w{5}[CW])/;
my $geneRe      = qr/,([atgcATGC-]+)/;
# Create an array for each dinucleotide in the %sequences
# that I'm interested in.
my %sequences = {
    'AA' => qr/(AA)/;
    'AC' => qr/(AC)/;
    'AG' => qr/(AG)/;
    'AU' => qr/(A[UT])/;
    'CA' => qr/(CA)/;
    'CC' => qr/(CC)/;
    'CG' => qr/(CG)/;
    'CU' => qr/(C[UT])/;
    'GA' => qr/(GA)/;
    'GC' => qr/(GC)/;
    'GG' => qr/(GG)/;
    'GU' => qr/(G[UT])/;
    'UA' => qr/([UT]A)/;
    'UC' => qr/([UT]C)/;
    'UG' => qr/([UT]G)/;
    'UU' => qr/(UU|TT)/
}
for ( my $keys in %sequences ) {
    my @{$keys};
}

foreach my $species (keys %Saccharomyces) {

    open(SPECIES,"<../data/$species/genome.csv") || die "Cannot open file: $!";
    my @text = <SPECIES>;    

    print "Working on $species right now.\n";

    # Each line of the `$species/genome.csv` contains the entire code of one gene.
    # We will loop through each line (gene) and calculate two things: the raw TRX
    # and delta-E scores, and the smoothed TRX and delta-E.
    foreach my $line (@text) {

        chomp($line);
        # Only run the calculation is we have genetic code.
        if ( defined $line && $line =~ m/e.{4,},[acgtACGT-]+/ ) {
       
            # Break up the CSV file, get useful info
            my $geneName = match($geneNameRe,$line);
            my $gene     = match($geneRe,$line);

            # Raw data
            # 
            # Calculates the trxValue and energyScore on a per-dinucleotide
            # basis. Then print the data to the respective CSV file in the
            # following format:
            #     gene,dinucleotide,position,trx.score,energy.score 
            # open(my $raw, ">>../data/$species/$geneName/raw.csv") || die "Cannot open file: $!";
            for ( my $position = 0; $position < length($gene); $position++ ) {
                
                my $dinucleotide = substr($gene,$position,2); 
            
                if ( $dinucleotide =~ m/[actgACTG-]{2}/ ) {     
                    
                    my $trxValue    = trxScore($dinucleotide);
                    my $energyScore = energyScore($dinucleotide);
                   
                    print $geneName . "," . $dinucleotide . "," . $position . "," . $trxValue . "," . $energyScore . "\n";

                    # Now that we've printed the dinucleotide, let's begin calculating the weight
                    for ( $dinucleotide ) {
                        for ( my \$keys in %Sequences ) {
                            @{$keys}++ when /${keys}/;
                            default { return }
                        }
                    }
                } 
            }
            # close $raw;           
        }
    }
}

sub match {

	my ( $re, $text ) = @_;

	if ( $text =~ $re ) {
		return $1;
	}
}
