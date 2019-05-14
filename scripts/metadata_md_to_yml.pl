#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

my $usage = <<EOS;
  Synopsis: metadata_md_to_yml.pl [options] FILENAME.md
  
  For markdown-format metadata used in the Legume Federation data store, circa 2018, 
  parse and generate a new yml-format file
  
  OPTIONS:
    -help         This message. 
EOS

my $help;

GetOptions (
  "help" =>           \$help,
);

die "\n$usage\n" if ( $help or not $ARGV[0] );

my $file_in = $ARGV[0];

my (%dataset_H, %file_HoH);
my (@dataset_A, @fileset_A);
my ($attr, $key, $val);
my $ct_dataset_attr=0;
my $ct_fileset_attr=0;

open (my $IN_FH, "<", $file_in) or die "can't open in $file_in: $!\n";

# open (my $OUT_YML_FH, ">", "$outbasefile.yml") or die "can't open out $outbasefile.yml: $!\n";

# Read input and store attributes and values in hashes.
# Also store attributes in an array, to preserve attribute order.
# Printing will happen at the end.
my $count_tics = 0;
my $count_this_multiline = 0;
my $count_attrs = 0;

print "%YAML 1.2\n---\n";

while (<$IN_FH>){
  chomp; 
  my $line = $_;

  # Attribute line
  if ($line =~ /^#+ *(\w.+\w) *$/){
  
    # convert e.g. "Publication DOI" to publication_doi
    $attr = lc($1); $attr =~ s/ /_/g; 
    if ( $count_attrs == 0 ) { print "$attr: " }
    else { print "\n$attr: " }
    $count_this_multiline = 0;
    $count_attrs++;
  } 
  elsif ($line =~ /<!---/ || $line =~ /^ *#/ ){
    next
  }
  elsif ($line =~ /^ *([^`].+)/ ){
    $val = $1;
    if ($count_tics % 2 == 0 ) { # outside a multi-line
      $count_this_multiline = 0;
      print "$val\n";
    }
    elsif ( $count_tics % 2 == 1 ) { # inside a multi-line
      if ( $count_this_multiline == 0 ) {print "\n"}
      print "  - $val\n";
      $count_this_multiline++;
    }
    else {die "WARNING: something unexpected with multi-lines; check the modulo\n" } 
  }
  elsif ($line =~ /```/){
    $count_tics++;
  }
  
}

__END__

Versions
2018 ...
v01 sc 07-11 started, based on metadata_md_to_csv.pl
v02 sc 07-25 print version directive at top (%YAML 1.2\n---\n)

