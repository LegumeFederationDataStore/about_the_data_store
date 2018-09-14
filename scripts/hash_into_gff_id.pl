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
    -out_file    (string) write to this file; otherwise, to stdout.
    -help        (boolean) This message.
EOS

my ($gff_file, $hash_file, $out_file, $help);

GetOptions (
  "gff_file=s" =>    \$gff_file,    # required
  "hash_file=s" =>   \$hash_file,   # required
  "out_file:s" =>    \$out_file,   
  "help" =>          \$help,
);

die "$usage" unless (defined($gff_file) and defined($hash_file));
die "$usage" if ($help);

# read hash in
open(my $HSH, '<', $hash_file) or die "can't open in input_hash, $hash_file: $!";
my %hash;
while (<$HSH>) {
  chomp;
  /(\S+)\s+(.+)/;
  my ($id, $hash_val) = ($1,$2);
  $hash{$id} = $hash_val;   # print "$id, $hash{$id}\n";
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
  else {
    my @fields = split(/\t/, $_);
    my $col9 = $fields[8];
    my @col9_k_v = split(/;/, $col9);
    my $col3 = $fields[2];
    my $attr_ct = 0;
    &printstr(join("\t", @fields[0..7]));
    &printstr("\t");
    if ( $col3 =~ /gene/i ) {
      $gene_name = $col9;
      $gene_name =~ s/ID=([^;]+)\;.+/$1/;
      $new_id = $hash{$gene_name};
      #print "GENE:[$gene_name] --> [$new_id]";
    }
    else {
      $col9 =~ s/$gene_name/$new_id/g;
      &printstr("$col9;");
    }
  }
  &printstr("\n");
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
v04 2018-09-13 Handle features with multiple parents (an exon or CDS can be a part of multiple mRNAs)
     Also remove option for taking in regex for stripping suffix such as "-mRNA-"
     because this can now be replaced afterwards with e.g. perl -pi -e 's/-mRNA-/./g'


