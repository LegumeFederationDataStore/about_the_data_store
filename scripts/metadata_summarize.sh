#!/bin/sh

# For each collection, derive html-format summary, from yml files

set -o errexit
set -o nounset 

BASEDIR="/usr/local/www/data/public" # base-level working directory
OUTERDIR=$1  # e.g. Glycine_max or Gene_families

echo "<p><b><big>Overview of data in this directory</big></b></p>" \
  > $BASEDIR/$OUTERDIR/_h5ai.header.html

echo "  Making header"
find $BASEDIR/$OUTERDIR -name "README*.yml" -print0 | sort -z | xargs -0 | perl -pe 's/ /\n/g' | grep -v "\/\." |
  xargs -I{} grep -H "subject:" {} | 
  perl -pe 's{(.+[^/]+)/([^/]+)/(README.\w+.yml:subject: +)(.+)}{$2 $4}' |
  perl -pe 's/(^\S+) (.+)/  <b>$1:<\/b> $2<br>/' \
  >> $BASEDIR/$OUTERDIR/_h5ai.header.html

echo "<hr>" >> $BASEDIR/$OUTERDIR/_h5ai.header.html

echo ""

