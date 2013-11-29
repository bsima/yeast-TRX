#!/usr/bin/perl

open (THR, "<threadcount.txt")|| die "could not open\n";
my @threadnumber = <THR>;
#print @threadnumber;
my $threadnumber_row = join ('',@threadnumber);
#print "$threadnumber_row\n";
my @threadnumber_row = split (/\s+/, $threadnumber_row);
my $threadnumber = @threadnumber_row[1];
print "threads = "."$threadnumber\n";
sleep(2);
$filecount = 0;
$file_iteration = 0;

for (my $th = 1; $th < $threadnumber + 1; $th++){
mkdir "GeneDirectoryRC"."$th" or die "\nGeneDirectoryRC folder already exists. Delete it and try again.\n";
$filecount = 0;
$file_iteration = 0;

$sp1_heading = "Scer";
$sp2_heading = "Spar";
$sp3_heading = "Smik";
$sp4_heading = "Sbay";

$dir = "GeneDirectory";
opendir(DIR,$dir)or die "can't open directory $genes:$!";
print"\n";
print "Alignment files in $dir are:\n";

while ($filename = readdir DIR){ # loop through alignment files

$filecount = $filecount + 1;
$file_iteration = $file_iteration + 1;
if ($file_iteration > $threadnumber){$file_iteration = 1;}

### take slice of files ...e.g. every nth file
print "$th\t"."$filecount\t"."$file_iteration\n";
if ($th != $file_iteration){next;}

$ORFname = substr($filename, 1, 7);
print "ORFname = "."$ORFname\n";

##reverse complement if Crick gene
$WatsonCrick = substr ($filename, 7, 1);
#print "$WatsonCrick\n";

$filelocation = "./GeneDirectory/"."$filename";

if (length $ORFname == 7){
open(INFILE, $filelocation) or die "Cannot open file";
}
else {next;}
  
  open (OUTFILE, ">"."./GeneDirectoryRC$th/"."$filename") || die " could not open output file\n";
  
  $flat1 = "";
  $flat2 = "";
  $flat3 = "";
  $flat4 = "";
  
  while(<INFILE>){
    chomp;
    @a = split(/\s+/, $_);
    #print "$_\n";
    if($_ =~ m/^$sp1_heading/){
      $flat1 = $flat1.$a[1];
      }
    elsif($_ =~ m/^$sp2_heading/){
      $flat2 = $flat2.$a[1];
      }
    elsif($_ =~ m/^$sp3_heading/){
      $flat3 = $flat3.$a[1];
      }
    elsif($_ =~ m/^$sp4_heading/){
      $flat4 = $flat4.$a[1];
      }
  }
  # sequence 1
  if ($WatsonCrick eq "C"){my $revcom = reverse $flat1; $revcom =~ tr/ACGTacgt/TGCAtgca/; $sequenceA = $revcom;}
  if ($WatsonCrick eq "W"){$sequenceA = $flat1;}
  print OUTFILE "$sp1_heading\t"."$sequenceA\n";
  # sequence 2
  if ($WatsonCrick eq "C"){my $revcom = reverse $flat2; $revcom =~ tr/ACGTacgt/TGCAtgca/; $sequenceB = $revcom;}
  if ($WatsonCrick eq "W"){$sequenceB = $flat2;}
  print OUTFILE "$sp2_heading\t"."$sequenceB\n";
  # sequence 3
  if ($WatsonCrick eq "C"){my $revcom = reverse $flat3; $revcom =~ tr/ACGTacgt/TGCAtgca/; $sequenceC = $revcom;}
  if ($WatsonCrick eq "W"){$sequenceC = $flat3;}
  print OUTFILE "$sp3_heading\t"."$sequenceC\n";
  # sequence 4
  if ($WatsonCrick eq "C"){my $revcom = reverse $flat4; $revcom =~ tr/ACGTacgt/TGCAtgca/; $sequenceD = $revcom;}
  if ($WatsonCrick eq "W"){$sequenceD = $flat4;}
  print OUTFILE "$sp4_heading\t"."$sequenceD\n";
  
  $sequencelengthA = length $sequenceA;
  
  # CREATE ANCESTRAL SEQUENCE RECONSTRUCTIONS
  print "GENERATING ANCESTRAL SEQUENCE RECONSTRUCTIONS\n";
  $sequencelengthPAML = $sequencelengthA;
  open (PHYLIP, ">phylip.txt")||die "could not open phylip.txt\n";
  print PHYLIP "      "."4   "."$sequencelengthPAML\n";
  print PHYLIP "sequence 1\n";
  print PHYLIP "$sequenceA\n";
  print PHYLIP "sequence 2\n";
  print PHYLIP "$sequenceB\n";
  print PHYLIP "sequence 3\n";
  print PHYLIP "$sequenceC\n";
  print PHYLIP "sequence 4\n";
  print PHYLIP "$sequenceD\n";
  close PHYLIP;
  system baseml;
  open (RST, "<rst")|| die "could not open rst file\n";
  @rst = <RST>;
  for (my $r = 0; $r < scalar @rst; $r++){
  $rst_row = @rst[$r];
  chomp $rst_row;
  #print "$rst_row\n";
  @rst_row = split (/\s+/,$rst_row);
  if (@rst_row[0] eq "node" && @rst_row[1] eq "#5"){
    #print "$rst_row\n\n";
    $asr_join = join('',@rst_row);
    #print "$asr_join\n\n";
    $asr5 = substr (lc $asr_join, 6, length ($asr_join) - 6);
    #print "$asr5\n\n";
    }
  if (@rst_row[0] eq "node" && @rst_row[1] eq "#6"){
    #print "$rst_row\n\n";
    $asr_join = join('',@rst_row);
    #print "$asr_join\n\n";
    $asr6 = substr (lc $asr_join, 6, length ($asr_join) - 6);
    #print "$asr6\n\n";
    }
  if (@rst_row[0] eq "node" && @rst_row[1] eq "#7"){
    #print "$rst_row\n\n";
    $asr_join = join('',@rst_row);
    #print "$asr_join\n\n";
    $asr7 = substr (lc $asr_join, 6, length ($asr_join) - 6);
    #print "$asr7\n\n";
    }  
}

# create caps for coding regions in asr5
$asr5new = '';
for (my $k = 0; $k <= length $asr5; $k++){
  $testbase = substr($sequenceD, $k, 1);
  $ASRbase = substr($asr5, $k, 1);
  if ($testbase eq "A" || $testbase eq "C" || $testbase eq "G" || $testbase eq "T"){$ASRbase = uc $ASRbase;}
  else {$ASRbase = $ASRbase;}
  $asr5new = $asr5new.$ASRbase;
}

# create caps for coding regions in asr6
$asr6new = '';
for (my $k = 0; $k <= length $asr6; $k++){
  $testbase = substr($sequenceC, $k, 1);
  $ASRbase = substr($asr6, $k, 1);
  if ($testbase eq "A" || $testbase eq "C" || $testbase eq "G" || $testbase eq "T"){$ASRbase = uc $ASRbase;}
  else {$ASRbase = $ASRbase;}
  $asr6new = $asr6new.$ASRbase;
}

# create caps for coding regions in asr7
$asr7new = '';
for (my $k = 0; $k <= length $asr7; $k++){
  $testbase = substr($sequenceB, $k, 1);
  $ASRbase = substr($asr7, $k, 1);
  if ($testbase eq "A" || $testbase eq "C" || $testbase eq "G" || $testbase eq "T"){$ASRbase = uc $ASRbase;}
  else {$ASRbase = $ASRbase;}
  $asr7new = $asr7new.$ASRbase;
}


# fix indels - Scer Spar
$asr7newfix = '';
for (my $k = 0; $k <= length $asr7new; $k++){
  $testbaseA = substr($sequenceA, $k, 1);
  $testbaseB = substr($sequenceB, $k, 1);
  $ASRbase = substr($asr7new, $k, 1);
  if ($testbaseA eq "-" && $testbaseB ne "-"){$ASRbase = $ASRbase;}
  if ($testbaseA ne "-" && $testbaseB eq "-"){$ASRbase = "-";}
  if ($testbaseA eq "-" && $testbaseB eq "-"){$ASRbase = "-";}
  if ($testbaseA ne "-" && $testbaseB ne "-"){$ASRbase = $ASRbase;}
  $asr7newfix = $asr7newfix.$ASRbase;
}


# fix indels - Spar Smik
$asr6newfix = '';
for (my $k = 0; $k <= length $asr6new; $k++){
  $testbaseB = substr($sequenceB, $k, 1);
  $testbaseC = substr($sequenceC, $k, 1);
  $ASRbase = substr($asr6new, $k, 1);
  if ($testbaseB eq "-" && $testbaseC ne "-"){$ASRbase = $ASRbase;}
  if ($testbaseB ne "-" && $testbaseC eq "-"){$ASRbase = "-";}
  if ($testbaseB eq "-" && $testbaseC eq "-"){$ASRbase = "-";}
  if ($testbaseB ne "-" && $testbaseC ne "-"){$ASRbase = $ASRbase;}
  $asr6newfix = $asr6newfix.$ASRbase;
}

# fix indels - Smik Sbay
$asr5newfix = '';
for (my $k = 0; $k <= length $asr5new; $k++){
  $testbaseC = substr($sequenceC, $k, 1);
  $testbaseD = substr($sequenceD, $k, 1);
  $ASRbase = substr($asr5new, $k, 1);
  if ($testbaseC eq "-" && $testbaseD ne "-"){$ASRbase = $ASRbase;}
  if ($testbaseC ne "-" && $testbaseD eq "-"){$ASRbase = "-";}
  if ($testbaseC eq "-" && $testbaseD eq "-"){$ASRbase = "-";}
  if ($testbaseC ne "-" && $testbaseD ne "-"){$ASRbase = $ASRbase;}
  $asr5newfix = $asr5newfix.$ASRbase;
}

print "\nlocal alignments\n";
print "sequence A\n";
#print "$sequenceA\n";
print "sequence B\n";
#print "$sequenceB\n";
print "sequence C\n";
#print "$sequenceC\n";
print "sequence D\n";
#print "$sequenceD\n";
print "ancestral sequence - node 5\n";
#print "$asr5\n";
print "ancestral sequence - node 6\n";
#print "$asr6\n";
print "ancestral sequence - node 7\n";
#print "$asr7\n";
print "\n";
$asr5length = length ($asr5);
$asr6length = length ($asr6);
$asr7length = length ($asr7);
print "sequence lengths\n";
print "$sequencelengthA\t"."$sequencelengthB\t"."$sequencelengthC\t"."$sequencelengthD\t"."$asr5length\t"."$asr6length\t"."$asr7length\n";
print "\n\n";
# add ASR to output file
print OUTFILE "asr7\t"."$asr7newfix\n";  
print OUTFILE "asr6\t"."$asr6newfix\n";  
print OUTFILE "asr5\t"."$asr5newfix\n";  
close OUTFILE;

} # end while loop
} # end threads loop
print "end program\n";
exit;