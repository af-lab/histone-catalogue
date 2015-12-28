#!/usr/bin/perl
use utf8;

## Copyright (C) 2011-2015 Carnë Draug <carandraug+dev@gmail.com>
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

## SYNOPSIS
##
##   protein_stats.pl path/for/dbfile
##
## DESCRIPTION
##
## It will read the store file of an HistoneSequencesDB object, and
## print to stdout Tex newcommands with protein stats (arg by lys ratios
## only at the moment).

use 5.010;
use strict;
use warnings;

use HistoneCatalogue;
use HistoneSequencesDB;

=func arg_lys_ratio
Compute arginine by lysine ratio.

Count all arginines and lysines in all of the gene proteins and return
their ratio to the other.

Args:
  $db (HistoneSequencesDB)
  @genes ([Gene])

Returns:
  float with arg/lys ratio
=cut
sub arg_lys_ratio
{
  my $db = shift;
  my @genes = @_;

  my $arg = 0; # count of arginine residues
  my $lys = 0; # count of lysine residues
  foreach my $products (map { $_->coding_products } @genes)
    {
      foreach my $acc (values %$products)
        {
          my $seq = $db->get_protein($acc)->seq;
          ## we know that some genes will encode proteins with the same sequence. We
          ## are counting those again on purpose. This gives us the Arg/Lys ratio for
          ## the proteins being expressed. Of course, we do not know the expression
          ## levels of each gene so we have to assume they are equal
          $arg++ while $seq =~ m/R/ig;
          $lys++ while $seq =~ m/K/ig;
        }
    }
  if ($arg == 0 || $lys == 0)
    { die "Found 0 arginine or lysine while calculating arg/lys ratio"; }
  return ($arg/$lys);
}

sub main
{
  if (@_ != 1)
    {
      print "Usage error -- no input arguments.\n";
      print "Correct usage is:\n";
      print "\n";
      print "  \$ protein_stats.pl path/for/dbfile\n";
      exit (1);
    }

  my $db = HistoneSequencesDB::read_db ($_[0]);

  my @core = $db->canonical_core();
  say HistoneCatalogue::latex_newcommand(
    "CoreArgLysRatio",
    arg_lys_ratio($db, @core),
    "Ratio of Total Number of Arginine and Lysines in all of the core histones"
  );

  foreach my $histone (@HistoneCatalogue::histones)
    {
      my @this_histones = grep {$_->histone_type eq $histone} @core;
      say HistoneCatalogue::latex_newcommand(
        "${histone}ArgLysRatio",
        arg_lys_ratio($db, @this_histones),
        "Ratio of Total Number of Arginine and Lysines in all of the ${histone} histones",
      );
    }

  my @linkers = grep {$_->isa ('CanonicalHistoneGene')} $db->linkers();
  say HistoneCatalogue::latex_newcommand(
    "LinkerArgLysRatio",
    arg_lys_ratio($db, @linkers),
    "Ratio of Total Number of Arginine and Lysines in all of the H1 histones"
  );

  return 0;
}

main (@ARGV);
