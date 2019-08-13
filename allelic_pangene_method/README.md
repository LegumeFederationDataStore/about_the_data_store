## allelic_pangene_method
Method for calculating a pan-gene set from a collection of annotations in a species (or related species).
This method uses a combination of homology information among CDS sequences from the included accessions, and
genomic synteny data. The objective is to identify sets of genes that are homologous and probably allelic 
(i.e. come from similar genomic neighborhoods). The method employs BLAST, DAGChainer, single-linkage clustering, 
and some special processing of coordinate data and BLAST results to make it all work. See details in 
the "methods" notes : [notes/allelic_pangene_methods.sh] (notes/allelic_pangene_methods.sh)
and the scripts (which are extremely simple, except for the clustering script by Michelle McMahon and
Mike Sanderson [modified by S. Cannon]).

Note that this repository doesn't contain the raw or intermediate data. The raw/initial files come
from the Data Store (https://legumeinfo.org/data/public/Glycine_max/ and https://legumeinfo.org/data/public/Glycine_soja/).



