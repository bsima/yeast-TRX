#!/usr/bin/perl -w
#use strict;

use threads;
use File::Copy; 
##########################################################################
my $threadcount = 1; # enter number of computer cores from 1 to 8
open (THR, ">threadcount.txt")|| die "could not open\n";
print THR "threads\t"."$threadcount\n";
close THR; 
##########################################################################
#goto SKIP;
print "\n\n\n";
print "\n\n\n";
print "----------------------------------------------------------------------------\n";
print "----------------------------------------------------------------------------\n";
print "  WC.pl is a Perl script invoking R, TRX scale and delta IE \n";
print "----------------------------------------------------------------------------\n";
print "----------------------------------------------------------------------------\n";
print "\n";
print " input alignments (homologous sequences) are stored in GeneDirectory folder\n";
print "----------------------------------------------------------------------------\n";
print "----------------------------------------------------------------------------\n";
print "\n\n\n";
print "\n\n\n";  
sleep (5);
# create ancestral sequences
#system "perl ReverseCompASR.pl\n";
# calculate and plot weighted correlations
system "weighted_correlation.pl\n";

print "END WC.pl PROGRAM\n";
exit;


