#!/usr/bin/env perl

# Quick controller script to do a grid of cluster distribution runs at
# different densities and tumble rates in 2d

# Run as ./rungrid.pl <jobid> <action>

use v5.10;
use warnings;

$thisjob = shift || 0;
$mode = shift || ($thisjob ? 'dryrun' : 'count');

$datadir = 'data';
mkdir $datadir unless -d $datadir;

@Ls = ( [100, 100], [120, 120], [140, 140], [160, 160], [180, 180],
        [200, 200], [220, 220], [240, 240], [260, 260], [280, 280],
        [300, 300], [320, 320], [340, 340], [360, 360], [380, 380],
        [400, 400] );

@rho = ( 0.01, 0.02, 0.04, 0.08, 0.16, 0.32 );
@tumble = (10, 1, 0.1, 0.01, 0.001);

$burnin = 50000;
$samples = 250;
$every = 1000;

# Derived parameters
$until = $every * $samples;

$jobid = 0;

for my $Lvec (@Ls) {
  my @L = @{$Lvec};
  $V = 1;
  for my $L (@L) {
    $V *= $L;
  }
  $Lcmd = join(" ", @L);
  $Lfn = join("x", @L);


  for my $rho (@rho) {
    my $Nlo = int($rho*$V);
    my $Nhi = $V-$Nlo;
    for my $tumble (@tumble) {
      for my $N ($Nlo, $Nhi) {
        ++$jobid;

        if($jobid == $thisjob) {
          my $ofn = "$datadir/L$Lfn-N$N-t$tumble-clusters.dat";
          my $cmd = "./crumble -L$Lcmd -N$N -t$tumble -b$burnin -e$every -u$until clusters";
          if($mode eq 'clobber' || ! -e $ofn) {
            say localtime . ": $cmd > $ofn";
            system "$cmd > $ofn" if ($mode eq 'go' || $mode eq 'clobber');
          }
        }
      }
    }
  }
}

if($mode eq 'count') {
  say "Job IDs run from 1-$jobid";
} else {
  say localtime . ": finished";
}
