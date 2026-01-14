#!/usr/bin/perl -w

use strict;
use File::Copy;
use FindBin;

my $home = "$FindBin::Bin";
#my $home=`pwd`;
#chomp $home;

if ($#ARGV != 1) {
  print "runNMA.pl <PDB_file> <chain ids>\n";
  exit;
}

my $pdb = $ARGV[0];
my $pdbCode = $ARGV[0];
my $chains = $ARGV[1];

my $dirname = "$pdbCode.$chains";

mkdir $dirname or print "cannot create $dirname\n";
chdir $dirname or die "cannot change to $dirname\n";

if (!-e "../$pdb") {
	die "cannot find file $pdb\n";
} 

copy("../$pdb","$pdb");
`$home/getChain.Linux $chains $pdb > pdb`;
`$home/part1 10 18; $home/useblz; $home/part2`;
rename("gnm1anmvector", "$pdbCode.1vector");
rename("gnm2anmvector", "$pdbCode.2vector");
rename("hinges", "$pdbCode.hinge");
rename("newcoordinat.mds", "$pdbCode.new");
`$home/processHinges $pdbCode.new $pdbCode.hinge $pdbCode.1vector $pdbCode.2vector 15 14.0`;

rename("$pdbCode.new.moved1.pdb", "$pdbCode.mode1.pdb");
rename("$pdbCode.new.moved2.pdb", "$pdbCode.mode2.pdb");

my @filelist = ("coordinates","upperhessian", "sortedeigen", "sloweigenvectors", "anm_length", "alpha.cor", "mapping.out", "pdb", "$pdbCode.hinge", "$pdbCode.new");
unlink @filelist;
unlink <*.mds12>; unlink <cross*>; unlink <*coor>; unlink <*cross>; unlink <*anm.pdb>; unlink <mod*>; unlink <newcoor*>; unlink <slow*>; unlink <*vector>; unlink <*.new>; unlink <.*>; unlink <fort*>;
