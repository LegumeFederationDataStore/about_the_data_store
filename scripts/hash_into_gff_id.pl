#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $usage = <<EOS;
  Synopsis: hash_into_gff_id.pl [options] -gff FILE -hash FILE 
  
  Read a key-value hash file and a GFF file. 
  Swap the feature IDs with the values from the hash file.
  
  Required:
    -gff         (string) GFF filename
    -hash        (string) Key-value hash filename, where first column has IDs from GFF file
  Options:
    -suppress    (string) File with list of components (mRNAs, CDSs, exons) or genes to exclude.
                   To exclude a component, use a splice suffix (GENE_A.1). To also
                   exclude a gene record, use the bare gene name (GENE_A)
    -out_file    (string) write to this file; otherwise, to stdout.
    -help        (boolean) This message.
EOS

my ($gff_file, $hash_file, $suppress, $out_file, $help);

GetOptions (
  "gff_file=s" =>    \$gff_file,    # required
  "hash_file=s" =>   \$hash_file,   # required
  "out_file:s" =>    \$out_file,   
  "suppress:s" =>    \$suppress,   
  "help" =>          \$help,
);

die "$usage" unless (defined($gff_file) and defined($hash_file));
die "$usage" if ($help);

# Read in hash of names
open(my $HSH, '<', $hash_file) or die "can't open in input_hash, $hash_file: $!";
my %hash;
while (<$HSH>) {
  chomp;
  /(\S+)\s+(.+)/;
  my ($id, $hash_val) = ($1,$2);
  $hash{$id} = $hash_val;   
  # print "$id, $hash{$id}\n";
}

# Read in list of genes and splice variants to suppress, if provided
my %suppress_hsh;
if (defined($suppress)) {
  open(my $SUP, '<', $suppress) or die "can't open in suppress, $suppress: $!";
  while (<$SUP>) {
    chomp;
    next if /^$/;
    $suppress_hsh{$_}++;
    #print "$_\n";
  }
}

my $OUT;
if ($out_file){
  open ($OUT, '>', $out_file) or die "can't open out out_file, $out_file: $!";
}

# Read in the GFF;
my $comment_string = "";
my $printed_comment_flag=0;
my ($gene_name, $new_id);
open(my $GFF, '<', $gff_file) or die "can't open in gff_file, $gff_file: $!";
while (<$GFF>) {
  chomp;
  
  # print comment line 
  if (/(^#.+)/) {
    $comment_string = "$1\n";
    &printstr("$comment_string");
  }
  else { # body of the GFF
    my @fields = split(/\t/, $_);
    my $col9 = $fields[8];
    my @col9_k_v = split(/;/, $col9);
    my $col3 = $fields[2];
    my $attr_ct = 0;
    #print $fields[4]-$fields[3], "\t", $fields[2], "\t";
    if ( $col3 =~ /gene/i ) {
      $gene_name = $col9;
      $gene_name =~ s/ID=([^;]+);.+/$1/;
      $new_id = $hash{$gene_name};
      #print "GENE:$gene_name\t";
      #print "[$col9]\n";
      if ( $suppress_hsh{$gene_name} ) { 
        # skip; do nothing
      }
      else {
        &printstr(join("\t", @fields[0..7]));
        &printstr("\t");
        $col9 =~ s/$gene_name/$new_id/g;
        &printstr("$col9;");
        &printstr("\n");
      }
    }
    else { # one of the gene components: CDS, mRNA etc.
      my $part_name = $col9;
      $part_name =~ s/ID=([^;]+);.+/$1/;
      $part_name =~ s/([^:]+):.+/$1/;
      #print "PART:$part_name\t";
      #print "{$col9}\n";

      if ($suppress_hsh{$part_name}) {
        # skip; do nothing
      }
      else {
        &printstr(join("\t", @fields[0..7]));
        &printstr("\t");

        # Split parents into an array
        my $parents = $col9;
        $parents =~ s/.+;Parent=(.+)/$1/i;
        my @parents_ary = split(/,/, $parents);
        #print join("=====", @parents_ary), "\n";

        foreach my $parent (@parents_ary) {
          if ( $suppress_hsh{$parent} ) { # strip this parent from list of parents
            $col9 =~ s/$parent//g;
            $col9 =~ s/,,/,/g;
            if (length($col9)<10) {
              warn "WARNING: Check elements related to $part_name. The 9th field is suspiciously short. " .
                "All parents may have been removed.\n"
              }
          }
          else {
            # do nothing (don't remove parent from string) because it's OK.
          }
        }
        $col9 =~ s/$gene_name/$new_id/g;
        &printstr("$col9;");
        &printstr("\n");
      }
    }
  }
}

#####################
sub printstr {
  my $str_to_print = shift;
  if ($out_file) {
    print $OUT $str_to_print;
  }
  else {
    print $str_to_print;
  }
}

__END__

Steven Cannon
Versions
v01 2014-05-22 New script. Appears to work.
v02 2017-01-10 Handle commented lines in GFF (the header)
v03 2018-02-09 Add more flexibility in printing commented lines (e.g. interspersed comments)
v04 2018-09-17 Handle features with multiple parents (an exon or CDS can be a part of multiple mRNAs)
     Also remove option for taking in regex for stripping suffix such as "-mRNA-"
       because this can now be replaced afterwards with e.g. perl -pi -e 's/-mRNA-/./g'
     Also take in a list of features to suppress (genes and/or splice variants)

Tests:
test2.gff3
  We want to suppress maker-Arahy.17-snap-gene-43.4-mRNA-1 :
    grep Arahy.17-snap-gene-43.4 lis.short_mRNAs_and_genes_to_remove 
  This should leave six rows in the output:
    grep -v "maker-Arahy.17-snap-gene-43.4-mRNA-1" test2.gff3 | cut -f3,9
      gene	ID=maker-Arahy.17-snap-gene-43.4;Name=maker-Arahy.17-snap-gene-43.4
      mRNA	ID=maker-Arahy.17-snap-gene-43.4-mRNA-2;Parent=maker-Arahy.17-snap-gene-43.4;Name=maker-Arahy.17-snap-gene-43.4-mRNA-2;_AED=0.35;_eAED=0.35;_QI=429|0|1|1|0|1|2|139|32
      three_prime_UTR	ID=maker-Arahy.17-snap-gene-43.4-mRNA-2:three_prime_utr;Parent=maker-Arahy.17-snap-gene-43.4-mRNA-2
      exon	ID=maker-Arahy.17-snap-gene-43.4-mRNA-2:exon:10142;Parent=maker-Arahy.17-snap-gene-43.4-mRNA-2
      CDS	ID=maker-Arahy.17-snap-gene-43.4-mRNA-2:cds;Parent=maker-Arahy.17-snap-gene-43.4-mRNA-2
      five_prime_UTR	ID=maker-Arahy.17-snap-gene-43.4-mRNA-2:five_prime_utr;Parent=maker-Arahy.17-snap-gene-43.4-mRNA-2
      five_prime_UTR	ID=maker-Arahy.17-snap-gene-43.4-mRNA-2:five_prime_utr;Parent=maker-Arahy.17-snap-gene-43.4-mRNA-2


