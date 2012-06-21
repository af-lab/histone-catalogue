#!/usr/bin/perl
## Copyright (C) 2011 Carnë Draug <carandraug+dev@gmail.com>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/>.

use 5.010;                      # Use Perl 5.10
use strict;                     # Enforce some good programming rules
use warnings;                   # Replacement for the -w flag, but lexically scoped
use Lingua::EN::Numbers;        # Turn a number into its 'english' form
use FindBin;                    # Locate directory of original perl script

use lib $FindBin::Bin;          # Add script directory to @INC to find 'package'
use MyVar;                      # Load variables
use MyLib;                      # Load functions

## this returns an array of references to an hash whose keys are the column names
my @data = MyLib::load_canonical;
my %canon;
foreach my $gene (@data) {
  my $symbol = $$gene{'gene symbol'};

  my $cluster = $1 if $symbol =~ m/^(HIST\d+)/;
  die "unable to identify cluster of $symbol" unless defined $cluster;

  ## find the start and end of each cluster (smallest and larger coordinate values
  ## of each gene on the cluster)
  foreach ( $$gene{'chromosome start coordinates'}, $$gene{'chromosome stop coordinates'} ) {
    ## start and end values may still be uninitialized so we check them first
    $canon{$cluster}{'start'} = $_ if !exists $canon{$cluster}{'start'} || $_ < $canon{$cluster}{'start'};
    $canon{$cluster}{'end'}   = $_ if !exists $canon{$cluster}{'end'}   || $_ > $canon{$cluster}{'end'};
  }

  ## count histones (by cluster, type and coding/pseudo)
  my $codness;
  if ($$gene{'pseudo'}) {
    $codness = "pseudo";
  } else {
    $codness = "coding";
  }
  my $histone = $1 if $symbol =~ m/^${cluster}($MyVar::histone_regexp)/;
  die "unable to identify histone type of $symbol" unless defined $histone;

  $canon{$cluster}{"total"}++;
  $canon{$cluster}{$histone}{"total"}++;
  $canon{$cluster}{$histone}{$codness}++;
  ## to create a hash specific for each gene
#  $canon{$cluster}{$symbol}{'start'} = $$gene{'chromosome start coordinates'};
}

for (keys %canon) {
  my $length  = abs ($canon{$_}{'start'} - $canon{$_}{'end'});
  $canon{$_}{'length'} = MyLib::pretty_length($length);
  say "$_ $canon{$_}{'length'}";
  say "Pseudo genes on $_ are $canon{$_}{'pseudo'}";
  say "Coding genes on $_ are $canon{$_}{'coding'}";
}

## write to results file
open($file, ">", $MyVar::results_clust) or die "Could not open $MyVar::results_clust for reading: $!";

close($file);
say "$_ is $canon{$_}{'length'}" for (keys %canon);
#for (keys %canon) {

#  say $canon{$_}{'start'};
#  my $remainder = length($length) % 3;
#  my $etc;
#  given ($remainder) {
#    when (0) {  $etc = sprintf("%1.3g", $length); }

#    default {  $etc = sprintf("%2.2g", $length); }
#  }
##  my $etc = sprintf("%1.${MyVar::size_precision}g", $length);
#  say "$_ has pot $pot from $length";
#say length( "23.4535433" );

#printf("%3.1f", value / pow(10, pot) )
#$pot = 0;
#while ( 10**($pot +3) ) {
#  $pot += 3;
#}
#printf("%3.1f", value / pow(10, pot) )
#use Data::Dumper;
#print Dumper %canon;
#my $number = Lingua::EN::Numbers::num2en($cluster =~ m/(\d*)$/);
#$nn =~ s/^(\w)/\U$1/g;
#say $nn;

