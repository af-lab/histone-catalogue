#!/usr/bin/perl
## Copyright (C) 2013-2014 Carnë Draug <carandraug+dev@gmail.com>
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

use 5.010;
use strict;
use warnings;
use File::Spec;

use MyLib; # FIXME: we should stop using this.

## This script compares the new data against some reference. In our case, the reference
## is Marzluff, W.F., Gongidi, P., Woods, K.R., Jin, J., Maltais, l.J. (2002) The human
## and mouse replication-dependent histone genes. Genomics (80) 5:487--498
## doi:10.1006/geno.2002.6850
##
## It will create the following files:
##    * table-reference_comparison.tex (LaTeX table with the list of differences
##      between the current analysis and a reference.
##    * variable-reference_comparison.tex
##
## Usage is:
##
## reference_comparison.pl --sequences path/for/sequences --results path/for/results --reference path/to/reference/

my %path = MyLib::parse_argv ("sequences", "results", "reference");

my %changes = (
  removed  => [], # were removed
  added    => [], # were added
  pseudo   => [], # changed to pseudo genes
  coding   => [], # changed to coding genes
  sequence => [], # sequence has changed
);

my @new_data  = MyLib::load_canonical ($path{sequences});
my @reference = load_marzluff ($path{reference});

foreach my $gene (@new_data) {
  my $symbol = $gene->{symbol};
  ## If there's a match, there will be only one, otherwise this will be empty
  my ($previous) = grep ($_->{symbol} eq $symbol, @reference);

  push (@{$changes{current}}, $symbol);
  ## not present on reference, therefore it's a new gene
  if (! $previous) {
    push (@{$changes{added}}, $symbol);
  ## used to be coding gene, now is a pseudo gene
  } elsif (! $previous->{pseudo} && $gene->{pseudo}) {
    push (@{$changes{pseudo}}, $symbol);
  ## used to be pseudo gene, now is a coding gene
  } elsif ($previous->{pseudo} && ! $gene->{pseudo}) {
    push (@{$changes{coding}}, $symbol);
  ## last difference we check is if the protein sequence has changed
  } elsif (! $previous->{pseudo} && ! $gene->{pseudo}) {
    my $old_seq = MyLib::load_seq ("protein", (keys %{ $previous->{proteins} })[0], $path{reference});
    my $new_seq = MyLib::load_seq ("protein", (keys %{ $gene->{proteins} })[0], $path{sequences});
    if ($old_seq->seq ne $new_seq->seq) {
      push (@{$changes{sequence}}, $symbol);
    }
  }
}

## Loop over the reference data to find any gene that may have been removed
foreach my $previous (@reference) {
  my $symbol = $previous->{symbol};
  if (! grep ($_->{symbol} eq $symbol, @new_data)) {
    push (@{$changes{removed}}, $symbol);
  }
}

## Write down results
my $tab_path = File::Spec->catdir($path{results}, "table-reference_comparison.tex");
open (my $tab_fh, ">", $tab_path) or die "Could not open $tab_path for writing: $!";

say {$tab_fh} "\\begin{tabular}{>{\\raggedright\\arraybackslash}p{\\dimexpr\\textwidth-2\\tabcolsep\\relax}}";
say {$tab_fh} "  \\toprule";

my %pairs = (
             "Changed sequences",              => \@{$changes{sequence}},
             "New genes",                      => \@{$changes{added}},
             "Removed genes",                  => \@{$changes{removed}},
             "Now identified as pseudo genes", => \@{$changes{pseudo}},
             "Now identified as coding genes", => \@{$changes{coding}},
             );
my @blocks;
foreach (keys %pairs) {
  push (@blocks, "  $_: \\\\\n  " . join (", ", sort (@{$pairs{$_}})) . "\\\\")
    if (@{$pairs{$_}});
}
say {$tab_fh} join ("\n  \\addlinespace\n", @blocks);

say {$tab_fh} "  \\bottomrule";
say {$tab_fh} "\\end{tabular}";

close ($tab_fh) or die "Couldn't close $tab_path after writing: $!";

my $var_path = File::Spec->catdir($path{results}, "variables-reference_comparison.tex");
open (my $var_file, ">", $var_path)
  or die "Could not open $var_path for writing: $!";

say {$var_file} HistoneCatalogue::latex_newcommand (
  "RemovedSinceReference",
  scalar @{$changes{removed}},
  "Number of removed genes since the reference"
);
say {$var_file} HistoneCatalogue::latex_newcommand (
  "AddedSinceReference",
  scalar @{$changes{added}},
  "Number of added genes since the reference"
);
say {$var_file} HistoneCatalogue::latex_newcommand (
  "PseudoSinceReference",
  scalar @{$changes{pseudo}},
  "Number of genes that have changed from coding to pseudo since the reference"
);
say {$var_file} HistoneCatalogue::latex_newcommand (
  "CodingSinceReference",
  scalar @{$changes{coding}},
  "Number of genes that have changed from pseudo to coding since the reference"
);
say {$var_file} HistoneCatalogue::latex_newcommand (
  "SequenceChangeSinceReference",
  scalar @{$changes{sequence}},
  "Number of genes whose sequence has changed since the reference"
);

close ($var_file)
  or die "Couldn't close $var_path after writing: $!";


## We can't use the normal load_canonical() because we didn't had enough
## data from his paper to rebuild the whole data. We only have the gene
## symbols and whether it's a coding or pseudo gene. But we want to use this
## script later with data downloaded by ourselves as the reference (to
## create comparisons between two specific dates) so we write load_marzluff
## an an exceptional case. The following hacks are required when using the
## Marzluff data we generated.
##    * use the gene symbol for the protein accession numbers. The filename
##      for the protein sequences are usually their accession numbers but they
##      are not mentioned in the paper. This could be a problem but on his
##      data, no gene encodes two proteins so we could use the gene symbol
##      to name the protein files.
##    * we can simplify a lot since we know there's no gene encoding multiple
##      proteins and there are only two fields on the CSV file;
sub load_marzluff {
  my $data_path = File::Spec->catdir($_[0], 'data.csv');
  my $csv = Text::CSV->new ({
    binary => 1,
    eol    => $/,
  }) or die "Cannot use Text::CSV: ". Text::CSV->error_diag ();
  open (my $file, "<", $data_path) or die "Could not open $data_path for reading: $!";

  $csv->column_names ($csv->getline ($file)); # read first line and sets it as the column name
  my $data = $csv->getline_hr_all ($file);    # reads all lines of file into an array of hashes (returns ref to array)
  close $file;

  my @reference = map {{
    'symbol'   =>  uc ($$_{'gene symbol'}),
    'pseudo'   =>  $$_{'pseudo'},
    'proteins' => {$$_{'gene symbol'} => ''},
  }} @$data;
  return @reference;
}

