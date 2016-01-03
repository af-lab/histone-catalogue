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
## print to stdout Tex newcommands with the following protein stats:
##
##  * arg by lys ratio (per histone type, and core vs linker)
##  * number of unique protein sequences per histone type

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

  my $arg = 0;
  my $lys = 0;
  my $it = $db->protein_iterator (@genes);
  while (my $p = $it->())
    {
      my $s = $p->seq;
      $arg++ while $s =~ m/R/ig;
      $lys++ while $s =~ m/K/ig;
    };
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

      my %seqs;     # sequence as keys to count unique proteins
      my $arg = 0;  # count of arginine residues
      my $lys = 0;  # count of lysine residues

      my $it = $db->protein_iterator (@this_histones);
      while (my $p = $it->())
        {
          my $s = $p->seq;
          $seqs{$s} = 1;

          ## we know that some genes will encode proteins with the same sequence. We
          ## are counting those again on purpose. This gives us the Arg/Lys ratio for
          ## the proteins being expressed. Of course, we do not know the expression
          ## levels of each gene so we have to assume they are equal
          $arg++ while $s =~ m/R/ig;
          $lys++ while $s =~ m/K/ig;
        }
      if ($arg == 0 || $lys == 0)
        { die "Found 0 arginine or lysine while calculating arg/lys ratio"; }

      say HistoneCatalogue::latex_newcommand(
        "${histone}ArgLysRatio",
        ($arg / $lys),
        "Ratio of Total Number of Arginine and Lysines in all of the ${histone} histones",
      );

      say HistoneCatalogue::latex_newcommand(
        "${histone}UniqueProteins",
        scalar (keys (%seqs)),
        "Number of unique proteins encoded by all histone $histone genes",
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
