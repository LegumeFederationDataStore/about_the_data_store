## TEMPLATES FOR GWAS, QTL, AND ASSOCIATED DATA

If you have genetic data that you would like to submit (maps, markers, association data, QTL results), 
we have templates to help you get your data into formats that we can incorporate into the Data Store
and associated applications. See these templates for GWAS and QTL data (note that these are somewhat 
simplified relative to the forms that are maintained in the Data Store).

  * [GWAS template](http://bit.ly/template__LIS_gwas)
  * [QTL template](http://bit.ly/template__LIS_qtl) - including map and marker worksheets 

## FULL DATA STORE SPECIFICATIONS FOR GENETIC DATA FILES

See the following specifications for the Data Store files related to genetic, association, and phenotype data.
These are the formats that determine file naming and data structure for these files in the Data Store.

[Genetic map files](https://docs.google.com/document/d/1d7O9K_5CTcRt_z8dlIXu1_Alh9RHpe3f5jy6xpKEerw/edit)
  * gensp.mixed.map.KEY4.identifier.expt.tsv
    - meta: genetic map experiment
    - data: linkage groups

  * gensp.mixed.map.KEY4.identifier.mrk.tsv
    - meta: genetic map name
    - data:	markers and linkage group positions

  * gensp.mixed.map.KEY4.identifier.qtl.tsv
    - meta: genetic map name, interval description
    - data: QTLs and linkage group start/end

  * gensp.mixed.map.KEY4.identifier.phen.tsv


[Genetic marker GFFs](https://docs.google.com/document/d/1hUwpu0ImMmKYb4tAftL250qLUF0rk2XSOUjNxbvc-yE/edit#)
  * gensp.strain.gnm#.mrk.KEY4.source.gff3
    - marker positions on a genome with additional attributes


[GWAS file](https://docs.google.com/document/d/1R10QpE-OIl0FEMamugcnL7-Z9difD6UwajbE4_vZtGg/edit#)
  * gensp.mixed.gwas.KEY4.dataset_name.gwas.tsv
    - meta: GWAS experiment details
    - data: QTL, marker, p-value
  * gensp.mixed.gwas.KEY4.dataset_name.phen.tsv

[QTL files](https://docs.google.com/document/d/1FWhhoL9iW7lgnNqI0VvEWtgNiUIe_1P4-Yqed4TT9wg/edit) (reported from other than map or GWAS expt)
  * gensp.mixed.qtl.KEY4.experimentidentifier.expt.tsv
    - meta: experiment details
    - data: QTLs
  * gensp.mixed.qtl.KEY4.experimentidentifier.marker.tsv
    - data: marker and QTL
  * gensp.mixed.qtl.KEY4.experimentidentifier.phen.tsv

[Phenotypes files](https://docs.google.com/document/d/10pX8nbDmrTp9FEcapsEe1tf711POiGh_ZzSzBIoFIc8/edit) (placed in appropriate directories)
  * gensp.mixed.xxx.KEY4.xxx.phen.tsv
    - data: phenotypes, ontology identifiers


