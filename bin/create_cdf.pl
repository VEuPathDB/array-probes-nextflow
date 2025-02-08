#!/usr/bin/perl

# Written by Gregory R Grant
# University of Pennsylvania, 2010

if(@ARGV < 2) {
    die "
This is for match-only arrays.

Usage: create_cdf.pl <cdf> <gene2probes> <tbase-pbase>

Where: 
        <cdf> is the name of the cdf outfile which already exists and
        has the following header, with details changed as necessary:
--------------------------------------------------------
[CDF]
Version=GC3.0

[Chip]
Name=Tgondii Match Only Array
Rows=476
Cols=476
NumberOfUnits=5220
MaxUnit=5220
NumQCUnits=0
ChipReference=
--------------------------------------------------------

         <gene2probes> is the file mapping genes to probes, with gene id followed by
         a tab delimited list of all probe ids mapping to that gene.

         <tbase-pbase> is the file output from get_pbase-tbase.pl

         <minProbes> The minimum number of probes that are required for a gene to be
         included in the cdf file (default is 3)

";
}

# see http://biosun1.harvard.edu/complab/dchip/HG-U133_Plus_2_dchip_example.txt
# and http://www.stat.lsa.umich.edu/~kshedden/Courses/Stat545/Notes/AffxFileFormats/cdf.html


open(INFILE, $ARGV[2]);
while($line = <INFILE>) {
    chomp($line);
    @a = split(/\t/,$line);
    $pbase{$a[0]} = $a[1];
    $tbase{$a[0]} = $a[2];
}
close(INFILE);

my $minProbes;

if ($ARGV[3]) {
   $minProbes = $ARGV[3];
   }
else { 
  $minProbes = 3;
     }
open(OUTFILE, ">>$ARGV[0]");

# ARGV[0] is the name of the file with the header above

open(INFILE, "$ARGV[1]");

# the file mapping genes -> probes

$unit = 1;
while($line = <INFILE>) {
    chomp($line);
    @a = split(/\t/,$line,2);
    @probes = split(/\t/,$a[1]);
    $numprobes = @probes;
    if($numprobes >= $minProbes) {
	print OUTFILE "[Unit$unit]\n";
	print OUTFILE "Name=NONE\n";
	print OUTFILE "Direction=1\n";
# change this to 2 if they are anti-sense
	print OUTFILE "NumAtoms=$numprobes\n";
	print OUTFILE "NumCells=$numprobes\n";
	print OUTFILE "UnitNumber=$unit\n";
	print OUTFILE "UnitType=3\n";
	print OUTFILE "NumberBlocks=1\n";
	print OUTFILE "\n";
	
	print OUTFILE "[Unit$unit";
	print OUTFILE "_Block1]\n";
	print OUTFILE "Name=$a[0]\n";
	print OUTFILE "BlockNumber=1\n";
	print OUTFILE "NumAtoms=$numprobes\n";
	print OUTFILE "NumCells=$numprobes\n";
	print OUTFILE "StartPosition=0\n";
	$num = $numprobes-1;
	print OUTFILE "StopPosition=$num\n";
	print OUTFILE "CellHeader=X\tY\tPROBE\tFEAT\tQUAL\tEXPOS\tPOS\tCBASE\tPBASE\tTBASE\tATOM\tINDEX\tCODONIND\tCODON\tREGIONTYPE\tREGION\n";
	for($i=0; $i<$numprobes; $i++) {

            my ($x, $y) = $probes[$i] =~ /(\d+)[\|\.\-\_\:](\d+)$/;
#	    @c = split(/[\.\|]/, $probes[$i]);
#            print STDERR "probes[$i] = $probes[$i]\n";

	    $p = $pbase{$probes[$i]};
	    $t = $tbase{$probes[$i]};
	    $ii = $i+1;
	    print OUTFILE "Cell$ii=$x\t$y\tN\tx\t$a[0]\t$ii\t13\tx\t$p\t$t\tx\tx\t-1\t-1\t99\t \n";
	}
	print OUTFILE "\n";
	$unit++;
    }
}
