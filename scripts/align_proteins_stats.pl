#!/usr/bin/perl
use utf8;

## Copyright (C) 2011-2016 Carnë Draug <carandraug+dev@gmail.com>
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
##   align_proteins_stats.pl path/for/protein/align_1 path/for/protein/align_2 path/for/protein/align_X
##
## DESCRIPTION
##
## It read several alignment of proteins and prints to stdout stats about
## them.  At the moment, that's only the percentage identity for each.

use 5.010;
use strict;
use warnings;

use File::Spec;

use Bio::AlignIO;

use HistoneCatalogue;

sub main
{
  if (@_ < 1)
    {
      print "Usage error -- no input arguments.\n";
      print "Correct usage is:\n";
      print "\n";
      print "  \$ align_proteins_stats.pl path/for/protein/align_1 path/for/protein/align_2 path/for/protein/align_X \n";
      exit (1);
    }

  for my $aln_file (@_)
    {
      if ((File::Spec->splitpath($aln_file))[2]
           !~ /aligned_($HistoneCatalogue::histone_regexp)_proteins.fasta/)
        { die "unable to identify histone type from alignment filename '$aln_file'"; }
      my $histone = $1;

      my $align = Bio::AlignIO->new(-file => $aln_file)->next_aln();

      ##
      ## (PID) Percentage IDentity per histone type
      ##

      say HistoneCatalogue::latex_newcommand(
        $histone."PID",
        $align->overall_percentage_identity,
        "Overall percentage identity between all histone $histone proteins"
      );
    }
  return 0;
}

main (@ARGV);
