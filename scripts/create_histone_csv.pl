#!/usr/bin/perl
use utf8;

## Copyright (C) 2015 Carnë Draug <carandraug+dev@gmail.com>
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
##   create_histone_csv.pl path/for/dbfile dir/to/save/newfiles
##
## DESCRIPTION
##
## It will read the store file of an HistoneSequencesDB object, and
## will create a csv files for canonical core, variants, and H1
## histones in the specified directory.  The directory must exist.
##
## This is completely useless and is not used anywhere for the manuscript.
## It is even dangerous because csv is a really poor format for genes.
## We only have this because Andrew wants it for his other projects.
## See https://github.com/af-lab/histone-catalogue/issues/3

use strict;
use warnings;

use File::Spec;
use Text::CSV;

use HistoneSequencesDB;

## I wish there was a file format for gene information but there is not.
## Because histones are simple, we can get away with a CSV file.
sub gene2csv
{
  my $fpath = shift;
  my @genes = @_;

  my $csv = Text::CSV->new ({ binary => 1, eol => $/})
    or die "Cannot use Text::CSV: ". Text::CSV->error_diag ();
  open (my $fh, ">:encoding(utf8)", $fpath)
    or die "Could not open $fpath for writing: $!";

  $csv->print ($fh, ["Histone type", "Symbol", "NCBI UID",
                     "chr_acc", "chr_start", "chr_end",
                     "Transcript accession", "Protein accession"]);

  foreach my $gene (@genes)
    {
      my @line = ($gene->histone_type, $gene->symbol, $gene->uid,
                  $gene->chr_acc, $gene->chr_start, $gene->chr_end, "", "");

      my $products = $gene->products;
      if (! keys %{$products})
        { $csv->print ($fh, \@line); }
      else
        {
          foreach my $transcript (keys %{$products})
            {
              @line[-2,-1] = ($transcript, $products->{$transcript});
              $csv->print ($fh, \@line);
            }
        }
    }
  close $fh
    or warn "Could not close $fpath after writing: $!";
}

if (@ARGV != 2)
  {
    print <<'END';
Usage error -- no input arguments.
Correct usage is:

  $ create_histone_csv.pl path/for/dbfile dir/to/save/newfiles
END
    exit (1);
  }

my $db_file = $ARGV[0];
my $db = HistoneSequencesDB::read_db($db_file);

my $out_dir = $ARGV[1];
if (! -d $out_dir)
  {
    print "directory '$out_dir' to save csv files does not exist";
    exit (2);
  }

my $out_path = sub { File::Spec->catdir ($out_dir, $_[0]); };

gene2csv (&$out_path ("canonical_core_histones.csv"), $db->canonical_core);
gene2csv (&$out_path ("variant_core_histones.csv"), $db->variants_core);

## this is both H1 and H5 (when they exist), and variant H1
gene2csv (&$out_path ("linker_histones.csv"), $db->linkers);
