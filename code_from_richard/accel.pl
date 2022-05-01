#!/usr/bin/env perl

# Accelerate a timeseries by a factor of n by echoing only comments and every nth non-comment line

use warnings;

$n = (shift or 1);

$m = 0;
foreach (<>) {
  if(/^# interval = (\w+)/) {
    print "# interval = ".($1*$n)."\n";
  }
  elsif(/^#/ || ($m++ % $n) == 0) {
    print;
  }
}
