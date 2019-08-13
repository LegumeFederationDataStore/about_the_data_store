# Objective: test methods for pan-gene-set construction
# Test using gene-based synteny + homology, using BLAST and DAGChainer
# Started: 2019-08-08; S. Cannon

# The method below was tested with six assemblies+annotations for soybean: two from Glycine soja and four from Glycine max.
# All data sets are public, and came from the Legume Federation / LegumeInfo Data Store: https://legumeinfo.org/data/public/

########## 
# 0. Method outline:
#
# 1.  Assemble the data:
#       CDS sequences for each assembly, as the cds_primaryTranscript.fna form.
#       Gene models for each assembly, as the gene_models_main.gff3 form.
# 
# 2.  Derive simple BED-format files from the gene model (GFF) files, and add the gene positional
#       information in the BED files to the CDS IDs.
#
# 3.  Run BLAST comparisons of each CDS file against each other (skipping self-comparisons).
#       To speed the computation, the comparisons can be made just for the upper-triangular N*(N+1)/2 pairs ...
#       ... and then derive the corresponding remaining pairs by swapping query and subject and sorting.
#
# 4.  Filter.
#     Filter to top hits in each direction (i.e. top hit for each query sequence), and apply identity threshold.
#     OPTIONALLY: Filter out local, repetitive matches.
#
# 5. Calculate synteny.
#     Run DAGChainer, to identify collections of matches in syntenic regions.
#
# 6.  Derive a file of all gene pairs per the above criteria (synteny + top match per query + minimum identity).
#
# 7.  Generate single-linkage clusters from the gene pairs.
#
# 8. Calculate summary stats.
#


##########
# 1.  Assemble the data:
#       CDS sequences for each assembly, as the cds_primaryTranscript.fna form.
#       Gene models for each assembly, as the gene_models_main.gff3 form.

# Get CDS sets (_primaryTranscript) from five assemblies, from LegFed Data Store
  ls 01_CDS
    glyma.Lee.gnm1.ann1.6NZV.cds_primaryTranscript.fna
    glyma.Wm82.gnm2.ann1.RVB6.cds_primaryTranscript.fna
    glyma.Wm82.gnm4.ann1.T8TQ.cds_primaryTranscript.fna
    glyma.Zh13.gnm1.ann1.8VV3.cds_primaryTranscript.fna
    glyso.PI483463.gnm1.ann1.3Q3Q.cds_primaryTranscript.fna
    glyso.W05.gnm1.ann1.T47J.cds_primaryTranscript.fna

  # glyma.Zh13 has splice-form IDs like ".m1". For consistency with the other annotations, delete the "m"
    perl -pi -e 's/\.m(\d+)$/.$1/' 01_CDS/glyma.Zh13.gnm1.ann1.8VV3.cds_primaryTranscript.fna

# Get GFFs
  mkdir 01_GFFs
  cp /data/Glycine_max/Wm82.gnm2.ann1.RVB6/glyma.Wm82.gnm2.ann1.RVB6.gene_models_main.gff3.gz 01_GFFs/
  cp /data/Glycine_max/Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.gff3.gz 01_GFFs/
  cp /data/Glycine_max/Lee.gnm1.ann1.6NZV/glyma.Lee.gnm1.ann1.6NZV.gene_models_main.gff3.gz 01_GFFs/
  cp /data/Glycine_max/Zh13.gnm1.ann1.8VV3/glyma.Zh13.gnm1.ann1.8VV3.gene_models_main.gff3.gz 01_GFFs/
  cp /data/Glycine_soja/PI483463.gnm1.ann1.3Q3Q/glyso.PI483463.gnm1.ann1.3Q3Q.gene_models_main.gff3.gz 01_GFFs/
  cp /data/Glycine_soja/W05.gnm1.ann1.T47J/glyso.W05.gnm1.ann1.T47J.gene_models_main.gff3.gz 01_GFFs/

  for file in 01_GFFs/*gz; do gunzip $file &
  done

##########
# 2.  Derive simple BED-format files from the gene model (GFF) files, and add the gene positional
#       information in the BED files to the CDS IDs.

# Make BED files from GFF
# NOTE: for later matching of gene file to the corresponding genome file, STRIP e.g. ".ann1.KEY"
  mkdir 02_gene_bed
  GENE_BED_DIR=02_gene_bed

  for path in 01_GFFs/*; do 
    base=`basename $path .gene_models_main.gff3`
    base_no_ann=`echo $base | perl -pe 's/\.ann\d+\.\w+//'`
    cat $path | awk -v OFS="\t" '$3=="gene" {print $1, $4, $5, $9}' | 
      perl -pe 's/ID=([^;]+);.+/$1/' > $GENE_BED_DIR/$base_no_ann.bed
  done

# Add positional information to gene names in CDS files, to allow easily 
# carrying this information through the BLAST search and into DAGChainer.
# The defline format needs to have this structure:
#    >CHROMOSOME__GENENAME__START__STOP
  mkdir 03_gene_pos_hash 04_CDS_posn

  for path in 02_gene_bed/*; do  
    base=`basename $path .bed`; 
    cat $path | awk '{print $4 "\t" $1 "__" $4 "__" $2 "__" $3 }' > 03_gene_pos_hash/$base.hsh
    hash_into_fasta_id.pl -fasta 01_CDS/$base.*fna -hash 03_gene_pos_hash/$base.hsh -suff_regex "\.\d+$" > 04_CDS_posn/$base.fna
  done

  # Move splice variant number back onto gene ID. 
  # Example: >glyma.Lee.gnm1.Gm01__glyma.Lee.gnm1.ann1.GlymaLee.01G057700.1__11256607__11259706
  perl -pi -e 's/(>.+)__(.+)__(\d+)__(\d+)\.(\d+)/$1__$2.$5__$3__$4/' 04_CDS_posn/*

##########
# 3.  Run BLAST comparisons of each CDS file against each other (skipping self-comparisons).
#       To speed the computation, the comparisons can be made just for the upper-triangular N*(N+1)/2 pairs ...
#       ... and then derive the corresponding remaining pairs by swapping query and subject and sorting.

# Prep and run the BLAST comparisons
  mkdir blastdb blastout scripts

  for path in 04_CDS_posn/*fna; do
    base=`basename $path .fna`
    echo $base
    makeblastdb -in $path -dbtype nucl -hash_index -parse_seqids -title $base -out blastdb/$base &
  done

  FASTADIR1=04_CDS_posn
  FASTADIR2=04_CDS_posn
  JOBMAX=16
  nohup scripts/do_blast_all.ksh $FASTADIR1 $FASTADIR2 $JOBMAX &


# Double the BLAST results (given qry-subj, print subj-qry)
  for path in blastout/*; do
    file=`basename $path .bln`
    part1=`echo $file | perl -pe 's/(.+)\.x\..+/$1/'`
    part2=`echo $file | perl -pe 's/.+\.x\.(.+)/$1/'`
    # echo "$part1 $part2"
    cat $path | awk -v OFS="\t" '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' |
      sort -k1,1 -k12nr,12nr \
      > blastout/$part2.x.$part1.bln
  done

##########
# 4.  Filter.
#     Filter to top hits in each direction (i.e. top hit for each query sequence), and apply identity threshold.
#     OPTIONALLY: Filter out local, repetitive matches. In tests with five soybean annotation sets, this step 
#     reduces average cluster sizes by about 1%, and excludes about 0.1% of all genes.

# Filter to top hit, and prepare for DAGChainer
  mkdir 05_dag_in
  MATCHLIST_DIR=05_dag_in

  for path in blastout/*; do 
    file=`basename $path .bln`
    cat $path | 
      awk -v OFS="\t" '$3>=95 {print $1, $2, $11}' | 
      top_line.awk | perl -pe 's/__/\t/g' > $MATCHLIST_DIR/$file.matchList &
  done

# NOTE/Checkpoint: The format at this point is like this: CHROMOSOME  GENENAME  START  STOP  CHROMOSOME  GENENAME  START  STOP  
#  glyma.Lee.gnm1.Gm01 glyma.Lee.gnm1.ann1.GlymaLee.01G000200.1  47692 49374 glyma.Wm82.gnm2.Gm01  glyma.Wm82.gnm2.ann1.Glyma.01G000100.1  27355 28320 6.04e-166
#  glyma.Lee.gnm1.Gm01 glyma.Lee.gnm1.ann1.GlymaLee.01G000300.1  72220 77056 glyma.Wm82.gnm2.Gm01  glyma.Wm82.gnm2.ann1.Glyma.01G000200.1  58975 67527 4.43e-127
#  glyma.Lee.gnm1.Gm01 glyma.Lee.gnm1.ann1.GlymaLee.01G000400.1  79900 81484 glyma.Wm82.gnm2.Gm15  glyma.Wm82.gnm2.ann1.Glyma.15G276800.1  51600225  51651653  1.50e-88


## OPTIONAL. Reduce repetitive matches. See note above regarding effect. This step is not recommended unless the intent
## is to produce clusters with few local paralogs.
#  mkdir 06_dag_in_nonrep
#  MATCHLIST_DIR=05_dag_in
#  MATCHLIST_NR_DIR=06_dag_in_nonrep
#
#  for path in $MATCHLIST_DIR/*; do 
#    file=`basename $path`
#    ~/bin/DAGCHAINER/accessory_scripts/filter_repetitive_matches.pl 50000 < $path > $MATCHLIST_NR_DIR/$file & 
#  done

##########
# 5. Calculate synteny.
#     Run DAGChainer, to identify collections of matches in syntenic regions.

# Run DAGChainer. NOTE: if running on the "filter_repetitive" results, set MATCHLIST_DIR=06_dag_in_nonrep
  MATCHLIST_DIR=05_dag_in

  for path in $MATCHLIST_DIR/*; do 
    file=`basename $path`
    run_DAG_chainer.pl -i $path 
  done

##########
# 6.  Derive a file of all gene pairs per the above criteria (synteny + top match per query + minimum identity).

# Generate file of all pairs
  cat 05_dag_in/*.aligncoords | awk '$1!~/^#/ {print $2 "\t" $6}' > tmp.all_synt_pairs.tsv

##########
# 7.  Generate single-linkage clusters from the gene pairs.

# Cluster
  blinkPerl_v1.1.pl -in tmp.all_synt_pairs.tsv -sum soy_pangene.sumry -out soy_pangene.clst.tsv

# Also make a version in which all genes in a set are on a single line:
  cat soy_pangene.clst.tsv | 
    awk -v ORS="" '$1==prev {print $2 "\t"} 
                   $1!=prev && NR==1 {print $2 "\t"; prev=$1} 
                   $1!=prev && NR>1 {print "\n" $2 "\t"; prev=$1} 
                   END{print "\n"}' |
      perl -pe 's/\t$//' > soy_pangene.sets.tsv


##########
# 8. Calculate summary stats.

# Stats
  head -2 soy_pangene.sumry
    Total number of sequences: 292941.
    Total number of clusters: 53504. 
  # Average cluster size:
    perl -le 'print 292941/53504'
      5.4751

  cut -f1 soy_pangene.clst.tsv | uniq -c | awk '{print $1}' | sort -n | histogram -n | head -20
    2 4934
    3 3435
    4 2593
    5 3410
    6 36391
    7 864
    8 369
    9 291
    10  271
    11  437
    12  285
    13  62
    14  33
    15  25
    16  29
    17  9
    18  11
    19  4
    20  8
    21  4

    # Max cluster size: 344

  cut -f1 soy_pangene.clst.tsv | uniq -c | 
    awk '$1==6 {eq++} $1<6 {lt++} $1>6 {gt++} 
         END{print "lt: " lt "\t" lt/(lt+eq+gt); print "eq: " eq"\t" eq/(lt+eq+gt); print "gt: " gt "\t" gt/(lt+eq+gt)}'
    lt: 14372 0.268615
    eq: 36391 0.680155
    gt: 2741  0.051230

##########
# Compare results with and without reducing repetitive matches with the filter_repetitive_matches.pl script

  wc -l tmp.all_synt_pairs_nonrep.tsv tmp.all_synt_pairs.tsv 
    1288186 tmp.all_synt_pairs_nonrep.tsv
    1290098 tmp.all_synt_pairs.tsv
      # 1912 more pairs in the full (nonreduced) result

  wc -l soy_pangene_dag_nonrep.clst.tsv soy_pangene.clst.tsv
    292136 soy_pangene_dag_nonrep.clst.tsv
    292941 soy_pangene.clst.tsv
      # 805 more genes in the full (nonreduced) result

# Numbers less than, equal to, or greater than the modal cluster size (6):
   # nonrep               # full monty
    lt: 14906 0.276298     14372 0.268615
    eq: 36950 0.684906     36391 0.680155
    gt: 2093  0.038796     2741  0.051230

# Conclude: the unfiltered results are probably a bit "better", in that they include more genes.
# The average cluster sizes are ~1% larger in the unfiltered results: 5.4751 vs. 5.4150

  
####################################################################################################
######### shell scripts / command files ###########

cat scripts/do_blast_all.ksh 
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

##########
cat /home/scannon/bin/top_line.awk
#!/usr/bin/awk -f
#
# NAME
#   top_line.awk - Filter tabular blast output to top hit per query.
# 
# SYNOPSIS
#   ./top_line.awk [INPUT_FILE(S)]
#     
# OPERANDS
#     INPUT_FILE
#         A file containing tabular blast output (or similar); the important
#         field is 0 (query ID in blast output).
#
# EQUIVALENT TO
#   awk '$1==prev && ct<2 {print; ct++} $1!=prev {print; ct=1; prev=$1}' FILE(s)

BEGIN { MAX = 1 } 

$1 == prev && count < MAX { print; count++ }
$1 != prev { print; count = 1; prev = $1 }


