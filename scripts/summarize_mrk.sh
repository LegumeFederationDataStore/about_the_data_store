#!/bin/sh

# For the gwas files, derive html-format summary, from marker sets

set -o errexit
set -o nounset 

BASEDIR="/usr/local/www/data/public" # base-level working directory
SPECIESDIR=$1   # e.g. Glycine_max
COLECTIONDIR=$2 # e.g. Wm82.gnm2.mrk.XJP2m

WORKDIR=$BASEDIR/$SPECIESDIR/$COLECTIONDIR
OUTFILE="_h5ai.header.html"

cat /dev/null > $WORKDIR/$OUTFILE
echo "<p><b><big>Overview of data in this directory</big></b></p>" > $WORKDIR/$OUTFILE

for path in $WORKDIR/*mrk*.gff3.gz; do
  file=`basename $path .gff3.gz`
  export platform=`echo $file | perl -pe 's/.+\.(\w+)$/$1/'`
  gzcat $path | perl -lne '
    $PLATFORM=$ENV{"platform"};
    $line=$_; 
    $line =~ s/^# r*(\w+): /$1\t/; 
    ($key, $value) = split(/\t/, $line); 
    if($key=~/Synopsis/){$synopsis=$value}; 
    END{
      if (defined($synopsis)){ print "<b>$PLATFORM</b>: $synopsis <br>" }
      else { print "<b>$PLATFORM</b>: [description to be added] <br>" }
    }' >> $WORKDIR/$OUTFILE
done

echo "<hr>" >> $WORKDIR/$OUTFILE

echo "  Output written to "
echo "    $WORKDIR/$OUTFILE"
echo ""

