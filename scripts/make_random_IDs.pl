#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my ($help, $mixed_case, $vowels, $out_file, $in_file);
my $width=4;
my $count=10;
my $too_many_tries = 10000;

my $usage = <<EOS;
  Synopsis: make_random_IDs.pl [options]

  Generates a list of random short names.
  
  Options:
    -count #    Positive integer. Return this many IDs. Default 10.
    -width #    width of generated random string used in the new IDs (>=2). Default 4.
    -in_file    File with IDs to avoid (one per line)
    -mixed_case Boolean. Include upper and lower case (discouraged if you expect these 
                  might be used naively with "sort" or as URLs). 
    -vowels     Boolean. Include vowels (discouraged if you want to avoid words). 
    -too_many_tries    Positive integer. Try this many times to find new unique IDs 
                  before bailing. Default 10000
    -out_file   Name of output file. Default stdout
    -help       Boolean. This message. 
EOS

GetOptions (
  "count:i" =>           \$count,  
  "width:i" =>           \$width,  
  "in_file:s" =>         \$in_file,
  "mixed_case" =>        \$mixed_case,
  "vowels" =>            \$vowels,
  "out_file:s" =>        \$out_file, 
  "too_many_tries:i" =>  \$too_many_tries, 
  "help" =>              \$help,
);

die "\n$usage\n" if ( $help );
die "\n-width must be >= 2\n" if ( $width < 2 );
die "\n-count must be >= 1\n" if ( $count < 1 );

my @s; # hash with characters to use in random ID strings

my @upper_cons =   qw(B C D F G H J K L M N P Q R S T V W X Y Z);
my @upper_vowels = qw(A E I O U);
my @lower_cons =   qw(b c d f g h j k l m n p q r s t v w x y z);
my @lower_vowels = qw(a e i o u);
my @nums =         qw(0 1 2 3 4 5 6 7 8 9);

if ( $mixed_case ){
  push(@s, @upper_cons, @lower_cons, @nums);
  if ($vowels){ push(@s, @upper_vowels, @lower_vowels) }
} else { # all upper case
  push(@s, @upper_cons, @nums);
  if ($vowels){ push(@s, @upper_vowels) }
}

my $characters = scalar(@s); # number of characters in alphabet to use for IDs
#print "$characters: ", join(" ", @s), "\n";

my %id_hsh;
my $new_ID;
my @IDs;
my $number_of_tries = 0;
my $OUT;
if ($out_file) { 
  open $OUT, '>', $out_file or die "can't open $out_file for writing: $!\n"; 
}

# If there is an $in_file, read it into the hash of IDs
if ($in_file){
  open my $IN, '<', $in_file or die "can't open $in_file for reading: $!\n";
  while (<$IN>){
    chomp;
    my $line = $_;
    next if $line =~ /^$/;
    $id_hsh{$line}++;
  }
}

my $count_OK_IDs = 0;
do {
  $new_ID="";
  for (1..$width) { $new_ID .= $s[int(rand()*$characters)] }
  if ($id_hsh{$new_ID}) {
    #print "$count_OK_IDs, $count, $number_of_tries REPEAT: $new_ID\n";
    $number_of_tries++;
    if ($number_of_tries > $too_many_tries) {
      die "Bailing out after $number_of_tries tries. Number of IDs generated: $count_OK_IDs\n"; 
    }
  } 
  else {
    $id_hsh{$new_ID}++;
    $count_OK_IDs++;
    if ( $out_file ) {
      print $OUT "$new_ID\n";
    } 
    else {
      print "$new_ID\n";
      #print "$count_OK_IDs, $count, $number_of_tries, $new_ID\n";
    }
  }
} until $count_OK_IDs >= $count;

__END__

VERSIONS

v01 2017-07-07 S. Cannon. Initial version. There is a problem: the "bailing out" warning
      isn't reported in some cases, despite hitting "too_many_tries"
v02 2017-07-09 S.C. Add option to take in a list of IDs to avoid. Fix counting problem.


