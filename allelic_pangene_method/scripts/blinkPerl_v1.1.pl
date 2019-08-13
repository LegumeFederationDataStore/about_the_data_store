#!/usr/bin/env perl
#use strict;
use warnings;
use Getopt::Long;

my $usage_brief = <<USAGE_BRIEF;
  Synopsis: blinkPerl_v1.1.pl -i infile [-o outfile -log logFile -s summaryFile]
 
  Required:
  -in        Input file name. A file with identifiers in the first two columns.
  
  Options:
  -sum       Summary file name. 
  -out       Output file name. If not provided, output is to stdout
  -log       Filename for log file
  -help      This message and additional information. 
USAGE_BRIEF

my $usage_rest = <<USAGE_REST;
  script to construct single-linkage clusters

  input:  file with two columns, each entry one sequence id, each line indicating a "hit",
   Format of the sequence ids is arbitrary, as long as each entry is followed by white space.
   This script can be used for any input where each line is an edge and the pair of words are nodes.
 
  output: file with two columns, first column cluster number, in "clXX" format, 
   second column sequence id. In the same format as in the input file.

  strategy:  build data structures "on the fly" (i.e., as the data is read into memory),
   linking clusters as necessary, making sure to merge smaller into larger to save time.
 
  Reference:  This is based on "blinkPerl.pl by Michelle McMahon, which in turn is based on 
   the algorithm developed in "blink.c" by MJ Sanderson, which see for further documentation.
   All three "blink" programs take the same input format and produce the same basic output, but have
   slight differences in the way the data is summarized, e.g., in the summary file here.
USAGE_REST

my ($infile, $outfile, $sumfile, $logfile, $help);

GetOptions (
  "in=s"   =>  \$infile,
  "sum:s"  =>  \$sumfile,
  "out:s"  =>  \$outfile,
  "log:s"  =>  \$logfile,
  "help"   =>  \$help,
);

die "\n$usage_brief\n$usage_rest\n" if ( $help );

die "\n$usage_brief\n" if (!($infile));

my %siClH = ();    # keys are seqids pointing to cluster numbers
my %clSiHoA = ();  # keys are cluster numbers pointing to arrays of seqids
my $nextCl = 1;

my ($FH_IN, $FH_OUT, $FH_LOG, $FH_SUM); # filehandles
open ($FH_IN, '<', $infile) or die "can't open in $infile: $!\n";
if ($logfile) {open ($FH_LOG, '>', $logfile) or die "can't open out $logfile: $!\n";}
while (<$FH_IN>) {
  chomp;
  next if (/^#|^\s*$/); # skip blank or commented lines
  my ($si1, $si2) = /\S+/g; 
  # if the two are equal, start a new cluster if necessary
  if ($si1 eq $si2)  { 
      if (! (exists $siClH{$si1})) {
          $siClH{$si1}=$nextCl; 
          $clSiHoA{$nextCl} = [ $si1 ]; 
          ++$nextCl;
          next;
      }
  }
  # if both seqids have been encountered, the two clusters merge, smallest merging into largest
  if (exists $siClH{$si1} && exists $siClH{$si2}) {
      ($cl1, $cl2) = ($siClH{$si1}, $siClH{$si2});
      next if ($cl1 == $cl2);
      $size1 = scalar @{$clSiHoA{$cl1}};
      $size2 = scalar @{$clSiHoA{$cl2}};
      if ($size1 < $size2) {$oldCl = $cl2; $newCl = $cl1;}
      else {$oldCl = $cl1; $newCl = $cl2;}
      if ($logfile) {print $FH_LOG "merging cluster $oldCl into $newCl because of $si1 or $si2\n";}
      for $si (@{$clSiHoA{$oldCl}}) 
        {$siClH{$si} = $newCl;}
      push @{$clSiHoA{$newCl}}, @{$clSiHoA{$oldCl}};
      delete $clSiHoA{$oldCl};
  }
  # if only one has been encountered, add the other seqid to the cluster 
  elsif (exists $siClH{$si1}) {
      $cl = $siClH{$si1};
      $siClH{$si2} = $cl;
      push @{$clSiHoA{$cl}}, $si2;
  } 
  elsif (exists $siClH{$si2}) {
      $cl = $siClH{$si2};
      $siClH{$si1} = $cl;
      push @{$clSiHoA{$cl}}, $si1;
  }
  # if neither is listed, start a new cluster.
  else {
    if ($logfile) {print $FH_LOG "making new cluster $nextCl of $si1 and $si2\n";}
    $siClH{$si1} = $nextCl;  
    $siClH{$si2} = $nextCl;
    $clSiHoA{$nextCl} = [ $si1, $si2 ];
    ++$nextCl;     
  }
}

# remove dups from cluster lists ------------------

%globalSiH = ();  # this is to make sure that no seqid is in 2 dif clusters
$globalSiHref = \%globalSiH;

for $cl (keys %clSiHoA) {
  removeDups ($clSiHoA{$cl}, $globalSiHref);
  $length2 = scalar @{$clSiHoA{$cl}};
  die ("cl $cl has no seqids\n") if ($length2 == 0);
}

# rename clusters with sequential numbers ---------
#  (creating new clHoA, old structures remain the same)

$newCl = 1;
for $oldCl (keys %clSiHoA) {
  $newClSiHoA{$newCl} = $clSiHoA{$oldCl};
  ++$newCl;
}

# print output ------------------------------------
@allSis = (keys %siClH);
@allCls = sort {$a <=> $b} (keys %newClSiHoA);
$totalCls = scalar @allCls;
$totalSis = scalar @allSis;
$checkTot = 0;

if ($sumfile) { 
  open ($FH_SUM, '>', $sumfile) or die "can't open out $sumfile: $!\n";
  print $FH_SUM "Total number of sequences: $totalSis.\nTotal number of clusters: $totalCls.\n";
  print $FH_SUM "\nCluster\tsequences\n";
}

if ($outfile) {
  open ($FH_OUT, '>', $outfile) or die "can't open out $outfile: $!\n";
}

for $cl (@allCls) {
  for $si (@{$newClSiHoA{$cl}}) {
    if ($outfile) {
      print $FH_OUT "cl$cl\t$si\n"; 
    }
    else {
      print "cl$cl\t$si\n"; 
    }
  }
  $numSis = scalar @{$newClSiHoA{$cl}};
  die ("NO seqids in cl $cl") if ($numSis == 0);
  
  if ($sumfile) {
    print $FH_SUM "$cl\t$numSis\n";
  }
  $checkTot = $checkTot + $numSis;
}
if ($sumfile) { 
  print $FH_SUM "Number of seqids assigned to clusters is $checkTot and total seqids in was $totalSis\n";
}

# -------------------------------------------------
# sub removeDups ----------------------------------
#  takes ref to an array; modifies original array
#  also takes ref to global si hash and adds to it... returns nothing

sub removeDups {
  my ($aref, $href) = @_;
  my @oldA = @{$aref};
  my %tempH = ();
  my @newA = ();

  for $x (@oldA) {
    if (! (exists $tempH{$x})) {
      if (! (exists ${$href}{$x}))
        {${$href}{$x}=1;}
      else {die ("$x is in mulitple clusters\n");} 
      push @newA, $x; 
      $tempH{$x} = 1;
      } 
   }
   @{$aref} = @newA;
   return;
}

__END__
2005-03-14 Michelle McMahon  blinkPerl.pl v. 1.0 
2017-02-02 Steven Cannon  some code restructuring: add GetOpts, etc. Make output files 
              optional (sending main output to stdout if no outfile is provided)
            Call it blinkPerl_v1.1.pl 

