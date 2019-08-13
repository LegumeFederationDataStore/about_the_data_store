#!/usr/bin/env ksh93
set -o errexit
set -o nounset

FASTADIR1=$1
FASTADIR2=$2
JOBMAX=$3

for qry_path in $FASTADIR1/*fna; do
  qry_base=`basename $qry_path .fna`
  for sbj_path in $FASTADIR2/*fna; do
    sbj_base=`basename $sbj_path .fna`
   if [[ "$qry_base" > "$sbj_base" ]]; then
     echo $qry_base.x.$sbj_base
      blastn -query $qry_path -db blastdb/$sbj_base -perc_identity 95 -evalue 1e-10 -outfmt 6 \
        -num_threads 6 -out blastout/$qry_base.x.$sbj_base.bln &
   fi
  done
done
wait


