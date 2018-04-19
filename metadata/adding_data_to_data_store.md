## PROCEDURE FOR ADDING A NEW DATA SET TO THE LEGUME FEDERATION DATA STORE

NOTE: The instructions below are for curators working on any instance of
LegFed Data Store - at e.g. soybase.org, peanutbase.org, legumeinfo.org etc. 
If you are a researcher or user of and you have a data set that you would like
to contribute, PLEASE <a href="https://legumeinfo.org/contact">CONTACT US!</a>. 
We would love to work with you. You are welcome to use the templates in this 
directory and begin preparing your data for submission, but the final uploading
will need to be done by curators with the affiliated database projects.

### Some resources
First, look at the documentation and templates at the following locations: 
  - https://legumeinfo.org/data/metadata/
  - Registry: http://bit.ly/LegFed_registry
  - Data Store protocols: https://legumeinfo.org/data/metadata/about_the_Data_Store.md


### Uploading the data to the local Data Store file system
The data store is accessible via command line from several servers.
As of winter, 2018, any of these servers can be used:
  - lis-stage.usda.iastate.edu 
  - soybase-stage.usda.iastate.edu 
  - legumefederation.usda.iastate.edu 
  - peanutbase-stage.usda.iastate.edu

Upload (scp) data to the private directory (and appropriate subdirectory) here:
  /usr/local/www/data/private/
  e.g.
  /usr/local/www/data/private/Glycine_max

### Naming the directories and files
Apply directory names, following the patterns described in 
  https://legumeinfo.org/data/metadata/about_the_Data_Store.md
... and including a key from http://bit.ly/LegFed_registry (and make appropriate changes at the registry)
For example, for methylation data initially at methylation/kim_et_al_2015/ 
```
  mv Gmax2.0 Wm82.gnm2.met1.K8RT
```

Apply file names, following the patterns described in
  https://legumeinfo.org/data/metadata/about_the_Data_Store.md
  
```
  cd  Wm82.gnm2.met1.K8RT
  mv W82_L1_Gm275.ctable.gz   glyma.Wm82.gnm2.met1.K8RT.methylation_L1.txt.gz 
  mv W82_L2_Gm275.ctable.gz   glyma.Wm82.gnm2.met1.K8RT.methylation_L2.txt.gz
  mv W82_RH1_Gm275.ctable.gz  glyma.Wm82.gnm2.met1.K8RT.methylation_RH1.txt.gz
  mv W82_RH2_Gm275.ctable.gz  glyma.Wm82.gnm2.met1.K8RT.methylation_RH2.txt.gz
  mv W82_SR1_Gm275.ctable.gz  glyma.Wm82.gnm2.met1.K8RT.methylation_SR1.txt.gz
  mv W82_SR2_Gm275.ctable.gz  glyma.Wm82.gnm2.met1.K8RT.methylation_SR2.txt.gz
```

### Filling out the README and MANIFEST files
Fill out the README file. The empty template is at https://legumeinfo.org/data/metadata/
but it is often easiest to copy a README from another data collection for the 
same species and then change the fields that need to be changed
```
  cp /usr/local/www/data/public/Glycine_max/Wm82.gnm2.ann1.RVB6/README.RVB6.md .
  mv README.RVB6.md README.K8RT.md
  perl -pi -e 's/RVB6/K8RT/' README.K8RT.md
```

Fill out the MANIFEST files.
```
  cp /usr/local/www/data/metadata/template__M* .
  rename 's/template__//' tem*
  rename 's/KEY/K8RT/' *KEY*
```

### Calculate the CHECKSUMs
Note the -r flag for the md5 command.
```
  KEY=K8RT
  rm CHECKSUM*
  md5 -r * > CHECKSUM.$KEY.md5
```

### Move the collection from public to private
Move the directory from from public to private, e.g.
```
  DIR=MY_NEW_DIRECTORY
  mv /usr/local/www/data/private/Glycine_max/$DIR /usr/local/www/data/public/Glycine_max/$DIR
```
... and note the change in the status file in private/Species_dir/, e.g. 
```
  echo $'\n'"Moved Wm82.gnm2.met1.K8RT to public on 2018-04-19 by YOUR NAME"$'\n' \
    >> private/Glycine_max/status.glyma.txt
```

### Update the CyVerse Data Store
Change to the species directory where your new collection sits and then 
also use the irods icommands to cd to the corresponding location at CyVerse:

```
  cd /usr/local/www/data/public/Glycine_max/
  ipwd     # see your path at CyVerse
  icd /iplant/home/shared/Legume_Federation/Glycine_max    # cd into the LegFed directory at CyVerse
  ils
```
If the directory exists at CyVerse, then icd into the directory and just push the updated files to it.
If the directories don't exist, then recursively push the directory.
If the directories are wrong, then BLOW AWAY the CyVerse directories (CAREFULLY) and 
recursively push the correct ones.

```
  DIR=MY_NEW_DIRECTORY
  ipwd
  ils
  iput -rf $DIR  # copy directory & files - assuming we are at the correct location locally and at CyVerse
  ils
  ils $DIR
```
  
