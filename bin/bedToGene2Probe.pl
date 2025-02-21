#!/usr/bin/perl

use strict;

use Getopt::Long;

my ($help, $gtfFeatureType, $bedFile, $geneGtfTag, $outputFile);

&GetOptions('help|h' => \$help,
            'gtfFeatureType=s' => \$gtfFeatureType,
            'bed=s' => \$bedFile,
            'geneGtfTag=s' => \$geneGtfTag,
            'outputFile=s' => \$outputFile
    );

if($help) {
    &usage();
}

sub usage {
    die "bedToProbe2Gene.pl --gtfFeatureType exon/CDS --bed probes.bed --geneGtfTag gene_id --outputFile probe2gene";
}

open(BED, $bedFile) or die "Cannot open $bedFile for reading: $!";
open(OUT, ">$outputFile") or die "Cannot open $outputFile for writing: $!";


my %dictionary;

my @a;
while(<BED>) {
    chomp;
    @a = split(/\t/, $_);

    my $gtfFeatureTypeCol = $a[14];
    my $probeName = $a[3];
    my $tags = $a[20];

    next unless($gtfFeatureTypeCol eq $gtfFeatureType);

    my ($gene) = $tags =~ /$geneGtfTag\s\"([^\"]+)/;

    $dictionary{$gene}->{$probeName} = 1;
}

foreach my $gene (keys %dictionary) {
    my @probes = keys %{$dictionary{$gene}};

    print OUT $gene . "\t" . join(",", @probes) . "\n";
}

close BED;
close OUT;
