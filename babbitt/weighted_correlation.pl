#!/usr/bin/perl -w
#use strict;
use File::Copy;
use Descriptive();

$windowsize = 200;    # size of window to shuffle within
$windowstep = 1;       # step size for sliding window 

mkdir "WCoutput" or die "\nWCoutput folder already exists. Delete it and try again.\n";
mkdir "WCoutputGraphic" or die "\nWCoutputGraphic folder already exists. Delete it and try again.\n";

@extantlist = ("Scer", "Spar", "Smik", "Sbay");
@ancestorlist = ("asr7", "asr7", "asr6", "asr5");

$filecount = 0;
$dir = "./GeneDirectoryRC1/";
opendir(DIR,$dir)or die "can't open directory $genes:$!";
print"\n";
print "Alignment files in $dir are:\n";

while ($filename = readdir DIR){ # loop through alignment files

$ORFname = substr($filename, 0, 6);
print "ORFname = "."$ORFname\n";

### obtain sequences to work with in this implementation ###
$filelocation = "./GeneDirectoryRC1/"."$filename";
#$filename = "/Users/gbabbitt/Desktop/RandTestTRX/R-2.10.1/bin/GeneDirectory2/"."eYBL029W.aln";
#$filename = "/Users/gbabbitt/Desktop/RandTestTRX/R-2.10.1/bin/GeneDirectory2/"."eYAL018C.aln";

if (length $ORFname == 6){
open(FILE, $filelocation) or die "Cannot open file";
}
else {next;}


while ($filename = readdir DIR){ # loop through alignment files

$ORFname = substr($filename, 1, 7);
print "ORFname = "."$ORFname\n";

$filelocation = "./GeneDirectoryRC1/"."$filename";

if (length $ORFname == 7){
open(INFILE, $filelocation) or die "Cannot open file";
}
else {next;}
   


# find ATG start codon

# calculate index position relative to start codon (do not count gaps)

################################################  
  
  while(<INFILE>){
    chomp;
    @a = split(/\s+/, $_);
    #print "$_\n";
    
    for (my $l = 0; $l <= scalar @extantlist; $l++){
        $species = $extantlist[$l];
        $speciesread = $a[0];
        if ($species eq $speciesread) {
            open (OUTFILE, ">"."./WCoutput/"."$ORFname"."_"."$species".".txt") || die " could not open output file\n";
            open (TEMPFILE, ">"."graph_"."$species".".txt") || die " could not open output file\n";
            $sequence = $a[1];
            # sequence without gap
            $sequenceNG = $sequence;
            $sequenceNG =~ s/-//g; 
            #########################
            # find ATG position
            #########################
            $countATG = 1;
            for (my $m = 0; $m <= length $sequenceNG; $m++){
            $triplet = substr($sequenceNG, $m, 3);
            if ($triplet eq "ATG"){$ATGposition = $countATG; goto SKIP;}
            else {$countATG = $countATG+1;}
            }
            SKIP:
            #print "ATG position at "."$ATGposition\n";
            # print column headers
            print OUTFILE "position\t"."weightedCORR\t"."weightedR2\t"."GCcontent\n";
            print TEMPFILE "position\t"."weightedCORR\t"."weightedR2\t"."GCcontent\n";
            # to screen
            print "analyzing "."$species\t"."$ORFname\n";
            #print "$species\n"."$sequenceNG\n\n";
            

           # sliding window
           for (my $ll = 0; $ll <= length $sequenceNG; $ll = $ll + $windowstep){
              $windowSEQUENCE = substr($sequenceNG, $ll, $windowsize);
              if (length $windowSEQUENCE == $windowsize) {
                #print "window\t"."$windowSEQUENCE\n";
                # compute dinucleotide frequency
                   weights();
                   #print "\n";
                   #print @weights;
                   #print "\n";
                # compute GC content
                   basecomp();                  
                # compute weighted correlation
                   weightedR();
                   #print "weighted mean X = "."$weightedmeanX\n";
                   #print "weighted mean Y = "."$weightedmeanY\n";
                   #print "weighted covariance XY = "."$weightedCOV\n";
                   #print "weighted covariance XX = "."$weightedCOVxx\n";
                   #print "weighted covariance YY = "."$weightedCOVyy\n";
                   #print "weighted correlation = "."$weightedCORR\n";
                   $position = $ll - $ATGposition;
                   print OUTFILE "$position\t"."$weightedCORR\t"."$weightedR2\t"."$GCcontent\n";
                   print TEMPFILE "$position\t"."$weightedCORR\t"."$weightedR2\t"."$GCcontent\n";
              }
            }
           
         close OUTFILE;
         close TEMPFILE;
        }
        
    }
  
  }
################################################


close INFILE;


# RUN R graphics AND COPY FILE TO OUTPUT FOLDER
print "making graphic\n";
graphics();  # make graphics
$oldfilename = "Rplots.pdf";
$newfilename = "WCoutputGraphic/"."$ORFname"."graph.pdf";
copy($oldfilename, $newfilename);



}

}



print "END weighted_correlation.pl PROGRAM\n";
exit;


sub basecomp{
  $total = 0;
  $GorC = 0;
  for (my $b = 0; $b < length $windowSEQUENCE; $b++){
      $total = $total+1;
      $BASE = substr(uc $windowSEQUENCE, $b, 1);
      if ($BASE eq "C" || $BASE eq "G") {$GorC = $GorC+1;}
           
    }
    $GCcontent = $GorC/($total+0.0001) 
}


sub weights{

@weights = ();
@trx = ();
@deltaIE = ();

$CGcount = 0;    
$CAcount = 0;
$TGcount = 0;
$GGcount = 0;
$CCcount = 0;
$GCcount = 0;
$GAcount = 0;
$TCcount = 0;
$TAcount = 0;
$AGcount = 0;
$CTcount = 0;
$AAcount = 0;
$TTcount = 0;
$ACcount = 0;
$GTcount = 0;
$ATcount = 0;
$TOTALcount = 0;

for (my $t = 0; $t < length $windowSEQUENCE; $t++){
      $DIMER = substr(uc $windowSEQUENCE, $t, 2);
      $TOTALcount = $TOTALcount + 1;
      if ($DIMER eq "CG"){$CGcount = $CGcount+1};
      if ($DIMER eq "CA"){$CAcount = $CAcount+1};
      if ($DIMER eq "TG"){$TGcount = $TGcount+1};
      if ($DIMER eq "GG"){$GGcount = $GGcount+1};
      if ($DIMER eq "CC"){$CCcount = $CCcount+1};
      if ($DIMER eq "GC"){$GCcount = $GCcount+1};
      if ($DIMER eq "GA"){$GAcount = $GAcount+1};
      if ($DIMER eq "TC"){$TCcount = $TCcount+1};
      if ($DIMER eq "TA"){$TAcount = $TAcount+1};
      if ($DIMER eq "AG"){$AGcount = $AGcount+1};
      if ($DIMER eq "CT"){$CTcount = $CTcount+1};
      if ($DIMER eq "AA"){$AAcount = $AAcount+1};
      if ($DIMER eq "TT"){$TTcount = $TTcount+1};
      if ($DIMER eq "AC"){$ACcount = $ACcount+1};
      if ($DIMER eq "GT"){$GTcount = $GTcount+1};
      if ($DIMER eq "AT"){$ATcount = $ATcount+1};
      
    }
    # local dinucleotide frequency 
    $CGfreq = $CGcount/($TOTALcount+0.0001);
    $CAfreq = $CAcount/($TOTALcount+0.0001);
    $TGfreq = $TGcount/($TOTALcount+0.0001);
    $GGfreq = $GGcount/($TOTALcount+0.0001);
    $CCfreq = $CCcount/($TOTALcount+0.0001);
    $GCfreq = $GCcount/($TOTALcount+0.0001);
    $GAfreq = $GAcount/($TOTALcount+0.0001);
    $TCfreq = $TCcount/($TOTALcount+0.0001);
    $TAfreq = $TAcount/($TOTALcount+0.0001);
    $AGfreq = $AGcount/($TOTALcount+0.0001);
    $CTfreq = $CTcount/($TOTALcount+0.0001);
    $AAfreq = $AAcount/($TOTALcount+0.0001);
    $TTfreq = $TTcount/($TOTALcount+0.0001);
    $ACfreq = $ACcount/($TOTALcount+0.0001);
    $GTfreq = $GTcount/($TOTALcount+0.0001);
    $ATfreq = $ATcount/($TOTALcount+0.0001);
    
    # weights array
    push (@weights, $CGfreq);
    push (@weights, $CAfreq);
    push (@weights, $TGfreq);
    push (@weights, $GGfreq);
    push (@weights, $CCfreq);
    push (@weights, $GCfreq);
    push (@weights, $GAfreq);
    push (@weights, $TCfreq);
    push (@weights, $TAfreq);
    push (@weights, $AGfreq);
    push (@weights, $CTfreq);
    push (@weights, $AAfreq);
    push (@weights, $TTfreq);
    push (@weights, $ACfreq);
    push (@weights, $GTfreq);
    push (@weights, $ATfreq);
    
    # trx scores array
    push (@trx, 43);
    push (@trx, 42);
    push (@trx, 42);
    push (@trx, 42);
    push (@trx, 42);
    push (@trx, 25);
    push (@trx, 22);
    push (@trx, 22);
    push (@trx, 14);
    push (@trx, 9);
    push (@trx, 9);
    push (@trx, 5);
    push (@trx, 5);
    push (@trx, 4);
    push (@trx, 4);
    push (@trx, 0);
    
    # delta IE scores
    push (@deltaIE, -26.9);
    push (@deltaIE, -20.0);
    push (@deltaIE, -23.3);
    push (@deltaIE, -24.3);
    push (@deltaIE, -21.4);
    push (@deltaIE, -22.9);
    push (@deltaIE, -23.7);
    push (@deltaIE, -28.2);
    push (@deltaIE, -19.6);
    push (@deltaIE, -23.6);
    push (@deltaIE, -17.2);
    push (@deltaIE, -18.5);
    push (@deltaIE, -15.8);
    push (@deltaIE, -19.0);
    push (@deltaIE, -18.9);
    push (@deltaIE, -15.7);
}


sub weightedR {

# calculate weighted mean X
  $weightedmeanX = 0;
  for (my $c = 0; $c < scalar @weights; $c++){
   $weight = @weights[$c];
   $trx = @trx[$c];
   $weightedmeanX = $weightedmeanX + ($weight*$trx);
   }
# calculate weighted mean Y
  $weightedmeanY = 0;
  for (my $cc = 0; $cc < scalar @weights; $cc++){
   $weight = @weights[$cc];
   $deltaIE = @deltaIE[$cc];
   $weightedmeanY = $weightedmeanY + ($weight*$deltaIE);
   }
# calculate weighted covariance X to Y
  $weightedCOV = 0;
  for (my $ccc = 0; $ccc < scalar @weights; $ccc++){
   $weight = @weights[$ccc];
   $trx = @trx[$ccc];
   $deltaIE = @deltaIE[$ccc];
   $weightedCOV = $weightedCOV + ($weight*($trx-$weightedmeanX)*($deltaIE-$weightedmeanY));
   }
# calculate weighted correlation
  $weightedCOVxy = $weightedCOV;
  $weightedCOVxx = 0;
  $weightedCOVyy = 0;
  $weightedCORR = 0;
  for (my $cccc = 0; $cccc < scalar @weights; $cccc++){
   $weight = @weights[$cccc];
   $trx = @trx[$cccc];
   $deltaIE = @deltaIE[$cccc];
   $weightedCOVxx = $weightedCOVxx + ($weight*($trx-$weightedmeanX)*($trx-$weightedmeanX));
   $weightedCOVyy = $weightedCOVyy + ($weight*($deltaIE-$weightedmeanY)*($deltaIE-$weightedmeanY));
   }
   $weightedCORR = $weightedCOVxy/(sqrt($weightedCOVxx*$weightedCOVyy)+0.0001);
   $weightedR2 = $weightedCORR*$weightedCORR;
}


sub graphics{
# open a pipe
open (Rinput, "| R --vanilla")||die "could not start R command line\n";
# read data into R
#print Rinput "library(bbmle)\n";

#print Rinput "data\n"; # print data to screen option
print Rinput "data1 = read.table('graph_Scer.txt', header = TRUE)\n"; 
$position_Scer = "data1\$position"; # window position
$weightedCORR_Scer = "data1\$weightedCORR"; # r
$weightedR2_Scer = "data1\$weightedR2"; # r square
$GCfreq_Scer = "data1\$GCcontent"; # GC content

#print Rinput "data\n"; # print data to screen option
print Rinput "data2 = read.table('graph_Spar.txt', header = TRUE)\n"; 
$position_Spar = "data2\$position"; # window position
$weightedCORR_Spar = "data2\$weightedCORR"; # r
$weightedR2_Spar = "data2\$weightedR2"; # r square
$GCfreq_Spar = "data2\$GCcontent"; # GC content

#print Rinput "data\n"; # print data to screen option
print Rinput "data3 = read.table('graph_Smik.txt', header = TRUE)\n"; 
$position_Smik = "data3\$position"; # window position
$weightedCORR_Smik = "data3\$weightedCORR"; # r
$weightedR2_Smik = "data3\$weightedR2"; # r square
$GCfreq_Smik = "data3\$GCcontent"; # GC content


#print Rinput "data\n"; # print data to screen option
print Rinput "data4 = read.table('graph_Sbay.txt', header = TRUE)\n"; 
$position_Sbay = "data4\$position"; # window position
$weightedCORR_Sbay = "data4\$weightedCORR"; # r
$weightedR2_Sbay = "data4\$weightedR2"; # r square
$GCfreq_Sbay = "data4\$GCcontent"; # GC content


print Rinput "par(mfrow=c(2,1))\n";
print Rinput "plot($position_Scer, $weightedCORR_Scer, mfg = c(1,1), type = 'l', col = 2, ylim = c(-0.8, -0.4), main = 'weighted correlation - $ORFname',  sub = 'red=Scer, green=Spar, blue=Smik, mag=Sbay', xlab = 'RELATIVE POSITION (ATG)', ylab = 'R value', new = TRUE)\n";
print Rinput "lines($position_Spar, $weightedCORR_Spar, col = 3, type = 'l')\n";
print Rinput "lines($position_Smik, $weightedCORR_Smik, col = 5, type = 'l')\n";
print Rinput "lines($position_Sbay, $weightedCORR_Sbay, col = 6, type = 'l')\n";


print Rinput "plot($position_Scer, $weightedR2_Scer, mfg = c(2,1), type = 'l', col = 2, ylim = c(0.16, 0.64), main = 'weighted R square - $ORFname',  sub = 'red=Scer, green=Spar, blue=Smik, mag=Sbay, black=GCcontent', xlab = 'RELATIVE POSITION (ATG)', ylab = 'R square' , new = TRUE)\n";
print Rinput "lines($position_Spar, $weightedR2_Spar, col = 3, type = 'l')\n";
print Rinput "lines($position_Smik, $weightedR2_Smik, col = 5, type = 'l')\n";
print Rinput "lines($position_Sbay, $weightedR2_Sbay, col = 6, type = 'l')\n";
print Rinput "lines($position_Scer, $GCfreq_Scer, col = 1, type = 'l')\n";

# write to output file and quit R
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";
}









