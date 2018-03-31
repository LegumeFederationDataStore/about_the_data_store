#!/usr/bin/env perl
use strict;
use Getopt::Long;
my $gnm;
my $ann;
my $name_prefix;
my $use_ID_as_base = 0;
GetOptions(
    "gnm=s" => \$gnm,
    "ann=s" => \$ann,
    #since we use less yucky prefixing for names than IDs; best to have supplied explicitly, as the extent to which it is less yucky depends on how much disambiguation is needed.
    "name=s" => \$name_prefix,
    #in rare cases, what we are given as ID makes more sense to use than what is given for Name
    "use_ID_as_base" => \$use_ID_as_base,
);

my $USAGE = "$0 --gnm <genome_prefix> --ann <annotation_prefix> [--name <name_prefix> --use_ID_as_base] <gff3_file> \n";

if (! defined $gnm) {
    print STDERR $USAGE;
    die "must supply --gnm prefix; if you are prefixing a gene annotation file, you should also supply an --ann prefix\n";
}
my %id_lookup;
my $do_fasta = 0;
while (<>) {
    if (/^#/) {
        print;
    }
    else {
        chomp;
        #even if a gene annotation file, the src_feature will be gnm
        s/^/$gnm/;
        #the rest will be context dependent
        my $prefix = (defined $ann ? $ann : $gnm);
        #we want names and IDs, sometimes we only get IDs, so will add a name based on the original ID but will get whatever is designated as the name prefix
        #TODO: maybe restrict this to types of features for which we care (genes, mRNA, chromosomes, supercontigs)
        if (! /Name=/) {
            s/ID=([^;]+)/ID=\1;Name=\1/;
        }
        my ($name);
        if ($use_ID_as_base) {
            ($name) = /ID=([^;]+)/;
        }
        else {
            ($name) = /Name=([^;]+)/;
        }
        #some files like phytozome use different conventions for IDs/Names; we will rewrite IDs to use our prefixing based on the Name as the core component of both
        my ($orig_id) = /ID=([^;]+)/;
        my $new_id = $prefix.$name;
        $id_lookup{$orig_id} = $new_id;
        s/ID=$orig_id/ID=$new_id/;
        if (/Parent=/) {
            my ($parent_old_id) = /Parent=([^;]+)/;
            my $parent_new_id = $id_lookup{$parent_old_id};
            if (!defined $parent_new_id) {
                die "could not find new id for Parent=$parent_old_id\n";
            }
            s/Parent=$parent_old_id/Parent=$parent_new_id/;
        }
        if ($use_ID_as_base) {
                s/Name=[^;]+/Name=$name/;
        }
        if (defined $name_prefix) {
                s/Name=/Name=$name_prefix/;
        }
        print $_,"\n";
    }
    if (/^##FASTA/) {
        $do_fasta = 1;
        last;
    }
}
if ($do_fasta) {
    while (<>) {
        if (/^>/) {
            s/^>/>$gnm/;
        }
        print;
    }
}
