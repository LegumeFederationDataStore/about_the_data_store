#!/bin/sh

# For the gwas files, derive html-format summary, from gwas.tsv

set -o errexit
set -o nounset 

BASEDIR="/usr/local/www/data/public" # base-level working directory
SPECIESDIR=$1   # e.g. Glycine_max
COLECTIONDIR=$2 # e.g. mixed.gwas.1W14

WORKDIR=$BASEDIR/$SPECIESDIR/$COLECTIONDIR
OUTFILE="TEST_h5ai.header.html"

cat /dev/null > $WORKDIR/$OUTFILE
echo "<p><b><big>Overview of data in this directory</big></b></p>" > $WORKDIR/$OUTFILE

for path in $WORKDIR/*gwas.tsv.gz; do
  gzcat $path |
		awk 'BEGIN{FS="\t"; ORS=""}
				 $1~/^Identifier/ {print "<b>" $2 "</b>: "}
				 $1~/^Name/ {print $2 "; "}
				 $1~/^Synopsis/ {print $2 "<br>\n"}' >> $WORKDIR/$OUTFILE
done

echo "<hr>" >> $WORKDIR/$OUTFILE

echo "  Output written to "
echo "    $WORKDIR/$OUTFILE"
echo ""

