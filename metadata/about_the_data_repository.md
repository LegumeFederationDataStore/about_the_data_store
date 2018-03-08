## DESCRIPTION OF DATA ORGANIZATION AND METADATA IN THE LEGUME FEDERATION PRIMARY REPOSITORY AND FEDERATION SITES

Return to the [metadata folder](.).


### OBJECTIVES AND PRINCIPLES 
1. A clear system (for curators and users)
  - Consistent, prescribed prefixes, abbreviations, etc.
2. Easy navigation for users
  - Limited directory depth (generally, only two directories deep)
  - Organization by species, then by genotype and data/analysis type
3. Borrow naming schemes from NCBI and Phytozome (but adapt in some places for clarity or 
    for other file types). Examples of similarities and differences:
  - NCBI doesn't maintain _primaryTranscript.fna files. Adopt this convention from Phytozome.
  - NCBI uses a _genomic.fna suffix for genome assembly (as well as a fairly complex directory 
    hierarchy for assembly components and tiling path, etc.). We use genome_main.fna.gz for
    genome assemblies, gene_models_main.gff3.gz for annotations.
  - Phytozome annotation names are versioned only per the assembly name, which doesn't 
    allow for the possibility of multiple annotation versions on an assembly.
4. In cases where there are multiple analyses (assembly, then annotation), be explicit 
   about naming the type of analysis, and the analysis version. This could be done using 
   "assembly" and "annotation", but assembly could refer to transcriptome or genomic 
   assemblies, and "assembly" and "annotation" are also lengthy. So, shorten generally to
   a three-letter abbreviation, e.g. "gnm" and "ann". 
5. Provide a unique, 4 character key for each dataset (a file or directory). This key will
   link the dataset unambiguously to its metadata.
6. Provenance is important. Retain the original README for data coming from other 
   repositories, and provide a structured README for this project - with same key name as 
   for the described files. Note the source repository clearly in the keyed project README. 
   Also record changed filenames (original and new), in the new README.


### PROCEDURE FOR MAKING A DATA REPOSITORY AT A FEDERATION SITE
In this case, the federation site is legumeinfo.org. Similar procedures may be used for other federation sites, but the paths used in the protocol described below will be different. 
Some naming conventions may also differ - but as long as the conventions are consistent 
and transparent at each site, mappings or translations should be possible.

The basic, simple concepts of federated repositories are:
  - There are repositories (filesystems, with ftp and/or http access) at member LegFed sites
  - A set of standards or practices for directory and file naming facilitates sharing of 
    file resources among projects and repositories, where appropriate. For example, 
    Cool Season Food Legume and LegumeInfo may both benefit from holding the reference 
    chickpea genome assemblies and annotations - particularly if each projects makes changes
    to any files for chickpea (e.g. GFFs for display in browsers) - and wants to provide those
    files to users, for the purpose of reproducibility, transparency, and convenience.
  - Any modifications of files vs. original forms in a primary repository should be noted in
    a way that can be reproduced. 
  - File or resource provenance, restrictions on use, etc., should be clearly noted.
  - Metadata should be stored in a way that makes it easy for people to keep the metadata with
    the files that the metadata describes.
  - Given the practices above (metadata records provenance, and stays with described files or 
    projects, and modifications are noted, and modified files get unique names), files or 
    directories can be re-hosted when appropriate - for example, at a data repository at CyVerse 
    AND at LegumeInfo or CoolSeasonFoodLegume.
  - The primary level of granularity of metadata is a directory, containing files that are 
    related in some meaningful biological way: a directory with a genome assembly and the 
    associated metadata; a set of annotations for a genome assembly and the associated 
    metadata, etc. 
  - Each such directory is keyed with a unique four-letter key name (as part of the 
    directory name), and the metadata files contained within the directory contain the 
    same key as the directory.
  - The keyed directories are typically organized within a human-readable containing directory,
    which is named by project convention. The top-level directory is typically named by
    genus and species, e.g. `Medicago_truncatula` or `Glycine_max`
  - The second-level directory/collection names have either three or four elements (depending 
    on content type), with elements dot-separated. 
      - The first element is genotype (e.g. Lee or Wm82 for the Lee or cultivars in soybean, or 
        V14167 for that accession identifier in Arachis duranensis). 
      - The second element indicates the genome and version: gnm1, gnm2, etc. 
      - For collections that depend on other collections (e.g. an annotation on an assembly, 
        or a genotyping dataset with variants called with respect to an assembly), the third element 
        indicates this secondary data type and version (annotation: ann1; markers: mrk1). 
      - The final element is the randomly-generated unique key (e.g. BXNC).
  - Directories containing dataset files should have md5 checksum files using the standard file
    naming prefixes associated with the directory and using an md5sum suffix. The contents should
    be in the format generated by the command: md5 -r FILE(s) > CHECKSUM.KEY#.md5
      <checksum> <filename>
    allowing validation as (e.g. on a Linux system): md5sum -c CHECKSUM.nnV9.md5


### DIRECTORY STRUCTURE AND NAMING
```
  Genus_species/GENOTYPE.ANALYSIS#.KEY/   or 
  Genus_species/GENOTYPE.ANALYSIS#.ANALYSIS#.KEY/ (for e.g. genome assembly, annotation)
    GENSP is a five-letter abbreviation of the GENus and SPecies, lower case,
      e.g. glyma for Glycine max or medtr for Medicago truncatula
    GENOTYPE is genotype name
    ANALYSIS# is e.g. "gnm2" for a genome assembly 2.0
    KEY is a 4-character unique identifier
    where analysis_type is abbreviated (generally by removing vowels and truncating):
      ann => annotation
      gnm => genome assembly
      tcp => transcriptome
      syn => synteny
      div => diversity
      gws => GWAS
      map => genetic map 
    In the version (v#.#), v1 implies v1.0. The decimal is used for subversions, e.g. v1.1
```
#### Example of directory names for a species (see examples of other data types below):
```
    Lupinus_angustifolius/
      P27255.tcp1.ycvS/
      Tanjil.gnm1.Qq0N/
      Tanjil.gnm1.ann1.nnV9/
      Tanjil.gnm1.syn1.TwnC/
      Tanjil.tcp1.p27w/
      Unicrop.tcp1.YVT4/
```
  
### METADATA: README, MANIFEST, CHECKSUM:
Retain the original README file(s) from the original repository, if any. Prefix
such files with "original_readme.", followed by the repository name and ".txt"

For the remaining prescribed metadata files, there should be a README file, two MANIFEST files, and a CHECKSUM file:

```
There should be four metadata files in each directory:
  - README.KEY#.md
  - MANIFEST.KEY#.correspondence
  - MANIFEST.KEY#.descriptions
  - CHECKSUM.KEY#.md5

The templated README file should have these 26 fields/sections, each prefixed with "#### ":
  - Identifier 
  - Provenance 
  - Source
  - Subject
  - Related To
  - Scientific Name
  - Taxid
  - BioProject
  - Scientific Name Abbrev
  - Genotype
  - Description
  - Dataset DOI 
  - GenBank Accession
  - Original File Creation Date
  - Local File Creation Date
  - Publication DOI 
  - Dataset Release Date
  - Publication Title
  - Contributors
  - Data Curators
  - Public Access Level
  - License
  - Keywords 
  - Citations 
  - File Transformation
  - Changes
```

### FILE NAMING SCHEME
  
Principles in the file naming are:
- Include brief, prescribed, descriptive prefixes:
  - gensp, e.g. glyma
  - analysis type and version, e.g. gnm1 or tcp1
  - four-character key# from the registry, e.g. RVB6
- Use consistent suffixes:
  - faa - protein FASTA
  - fna - nucleotide FASTA
  - gff3 - GFF
  - md - README files
  - tsv - ascii table file
  - txt - ascii files other than READMEs
  - md5 - md5 checksum file for contents of folders

Within the GENOTYPE.ANALYSIS#/ directories, files have this format:
gensp.GENOTYPE.ANALYSIS#.KEY#.file_type.ext.gz
where file_type is a prescribed, abbreviated description of the contents, 
and where KEY# is a four-letter character string, uniquely identifying a file
or files referred to by the README file containing the same key.
Examples of file_type abbreviations (these generally follow Phytozome and GenBank patterns,
but with some exceptions - for example, GenBank doesn't have _primaryTranscript forms, and
neither subdivides GFF files into subtypes as needed in legumeinfo and peanutbase for 
database loading and browsers):

  Annotation files:
```
    gensp.GENOTYPE.ANALYSIS#.ANALYSIS#.KEY#.gene_model_main.gff3.gz 

    gensp.GENOTYPE.ANALYSIS#.ANALYSIS#.KEY#.annot.ahrd.txt.gz
    gensp.GENOTYPE.ANALYSIS#.ANALYSIS#.KEY#.cds.fna.gz
    gensp.GENOTYPE.ANALYSIS#.ANALYSIS#.KEY#.cds_primaryTranscript.fna.gz
    gensp.GENOTYPE.ANALYSIS#.ANALYSIS#.KEY#.gene_for_gbrowse.gff3.gz  -- e.g. where this specialized file has been generated
    gensp.GENOTYPE.ANALYSIS#.ANALYSIS#.KEY#.protein.faa.gz
    gensp.GENOTYPE.ANALYSIS#.ANALYSIS#.KEY#.protein_primaryTranscript.faa.gz
```
  Genome assembly files:
```
    gensp.GENOTYPE.ANALYSIS#.KEY#.genome_main.fna.gz   -- the primary assembly, with pseudomolecules & remaining scaffolds
    gensp.GENOTYPE.ANALYSIS#.KEY#.agp.txt.gz  -- Tiling path (order and orientation of scaffolds), in AGP format
    gensp.GENOTYPE.ANALYSIS#.KEY#.scaffolds.fna.gz  -- All scaffolds, or as otherwise specified in the MANIFEST
```
  Transcriptome assembly files:
```
    gensp.GENOTYPE.ANALYSIS#.KEY#.METHOD.fna.gz -- e.g. .ycvS.Trinity.fna
```
  README files
```
    README.KEY#.md  -- template-based README file for this repository, in markdown format
```

### EXAMPLES OF DIRECTORIES AND FILES:
  
#### Genome assembly
```
  Lupinus_angustifolius/
    Tanjil.gnm1.Qq0N/
      CHECKSUM.Qq0N.md5
      MANIFEST.Qq0N.correspondence
      MANIFEST.Qq0N.descriptions
      README.Qq0N.md
      lupan.Tanjil.gnm1.Qq0N.genome_main.fna.gz
      lupan.Tanjil.gnm1.Qq0N.scaffolds.fna.gz
      lupan.Tanjil.gnm1.Qq0N.scaffolds_unassigned.fna.gz
      original_readme.lupinexpress.txt

  Trifolium_pratense/
    MilvusB.gnm2.gNmT/
      CHECKSUM.gNmT.md5
      MANIFEST.gNmT.correspondence
      MANIFEST.gNmT.descriptions
      README.gNmT.md
      tripr.MilvusB.gnm2.gNmT.genome_main.fna.gz

  Vigna_angularis/
    Gyeongwon.gnm3.JyYC/
      CHECKSUM.JyYC.md5
      MANIFEST.JyYC.correspondence
      MANIFEST.JyYC.descriptions
      README.JyYC.md
      vigan.Gyeongwon.gnm3.JyYC.genome_main.fna.gz
```

#### Genome annotation
```
  Lupinus_angustifolius/
    Tanjil.gnm1.ann1.nnV9/
      CHECKSUM.nnV9.md5
      MANIFEST.nnV9.correspondence
      MANIFEST.nnV9.descriptions
      README.nnV9.md
      lupan.Tanjil.gnm1.ann1.nnV9.cds_all.fna.gz
      lupan.Tanjil.gnm1.ann1.nnV9.gene_models.gff3.gz
      lupan.Tanjil.gnm1.ann1.nnV9.gene_models_main.gff3.gz
      lupan.Tanjil.gnm1.ann1.nnV9.gene_models_miRNA.gff3.gz
      lupan.Tanjil.gnm1.ann1.nnV9.gene_models_pchr_plus_scaff.gff3.gz
      lupan.Tanjil.gnm1.ann1.nnV9.gene_models_snRNA.gff3.gz
      lupan.Tanjil.gnm1.ann1.nnV9.gene_models_tRNA.gff3.gz
      lupan.Tanjil.gnm1.ann1.nnV9.protein.faa.gz
      lupan.Tanjil.gnm1.ann1.nnV9.rgene_models_RNA.gff3.gz
      lupan.Tanjil.gnm1.ann1.nnV9.scaffolds_in_pchr.gff3.gz
      original_readme.lupinexpress.txt
  
  Trifolium_pratense/
    MilvusB.gnm2.ann1.DFgp/
      CHECKSUM.DFgp.md5
      MANIFEST.DFgp.correspondence
      MANIFEST.DFgp.descriptions
      README.DFgp.md
      tripr.MilvusB.gnm2.ann1.DFgp.cds.fna.gz
      tripr.MilvusB.gnm2.ann1.DFgp.cds_primaryTranscript.fna.gz
      tripr.MilvusB.gnm2.ann1.DFgp.gene_models_gbrowse.gff3.gz
      tripr.MilvusB.gnm2.ann1.DFgp.gene_models_main.gff3.gz
      tripr.MilvusB.gnm2.ann1.DFgp.info_annot_AHRD.txt.gz
      tripr.MilvusB.gnm2.ann1.DFgp.protein.faa.gz
      tripr.MilvusB.gnm2.ann1.DFgp.protein_primaryTranscript.faa.gz

  Vigna_angularis/
    Gyeongwon.gnm3.ann1.3Nz5/
      CHECKSUM.3Nz5.md5
      MANIFEST.3Nz5.correspondence
      MANIFEST.3Nz5.descriptions
      README.3Nz5.md
      vigan.Gyeongwon.gnm3.ann1.3Nz5.cds.fna.gz
      vigan.Gyeongwon.gnm3.ann1.3Nz5.cds_primaryTranscript.fna.gz
      vigan.Gyeongwon.gnm3.ann1.3Nz5.gene_models_main.gff3.gz
      vigan.Gyeongwon.gnm3.ann1.3Nz5.info_annot_AHRD.txt.gz
      vigan.Gyeongwon.gnm3.ann1.3Nz5.protein.faa.gz
      vigan.Gyeongwon.gnm3.ann1.3Nz5.protein_primaryTranscript.faa.gz
```

#### Transcriptome assembly
```
  Lupinus_angustifolius/
    P27255.tcp1.ycvS/
    MANIFEST.ycvS.correspondence
    MANIFEST.ycvS.descriptions
    README.ycvS.md
    lupan.P27255.tcp1.ycvS.Trinity.fna.gz
    original_readme.lupinexpress.txt
```

### CURATION PROTOCOLS
Notes below describe procedures specifically at the data filesystem associated with legumeinfo.
Presumably, similar practices would be used for other "federation member" repositories.

Data is initially populated and managed at legumeinfo - though might similarly 
initially be managed at another repository.

Data are organized by species (or major cross-cutting data type such as "gene_families", 
in directories at /usr/local/www/data/, in either public/ or private/. 
These directories are accessed at lis-stage.agron.iastate.edu, and are  
mounted in production at legumeinfo.org at the same directory address, and are accessible
via browser at http://lis-stage.agron.iastate.edu/data/ or http://legumeinfo.org/data/ 
or http://peanutbase.org/data/  (and soon at http://legumefederation.org/data/).

Data for this data directory will typically come from another repository, sometimes via 
a local data repository. For LIS development at Ames, this repository is maintained at the
machines "cicer" and "cicer" -- where the /data directory is mounted from the  NAS "mora".

The primary purpose for the repository is to expose data sets that are used at LIS or the other
respective genomic data portals (GDPs). To that end, not every piece of genomic data for
a species needs to be housed in the repository. (There may be some instances, though, where
extensive, multifaceted data sets are held for some species).

### MISCELLANEOUS DATA STANDARDS
This is a set of conventions and standards that we have agreed to adopt with respect to certain types 
of data.

- Peptides SHALL be given the same names and identifiers as the transcripts from which they were derived.
This means, for example, that in the case of providers who have suffixed the names of transcripts with a
".p" to discriminate the peptides derived from them, we will have slightly different naming within the 
federation. 
