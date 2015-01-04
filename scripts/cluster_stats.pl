#!/usr/bin/perl
## Copyright (C) 2011-2014 Carnë Draug <carandraug+dev@gmail.com>
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
use File::Spec;                 # Perform operation on file names
use List::Util;                 # Includes min and max
use Number::Format;             # pretty format of long numbers with SI prefix and precision

use FindBin;                    # Locate directory of original perl script
use lib $FindBin::Bin;          # Add script directory to @INC to find 'package'
use MyVar;                      # Load variables
use MyLib;                      # Load functions

## This script will calculate stats for each of the histone clusters. It will
## create the following files:
##    * table-histone_catalogue.tex (a very long LaTeX table with all canonical
##      histones, their UIDs, and transcript and protein accession numbers)
##    * variables-cluster_stats.tex (LaTeX variables with the number of histones
##      in each cluster, how many are pseudo and coding, the length of each
##      each cluster, and the location in the genome of each cluster)
##
## Usage is:
##
## cluster_stats.pl --sequences path/for/sequences --results path/for/results

## Check input options
my %path = MyLib::parse_argv("sequences", "results");

## Read the data and create a data structure for the analysis
my %canon;   # organized by cluster, with counts of histones and other info
my %types;   # histone types as keys for arrays of histones of that type
my @data = MyLib::load_canonical ($path{sequences});
foreach my $gene (@data) {
  my $symbol  = $$gene{'symbol'};
  my $cluster = "HIST" . $$gene{'cluster'};
  my $histone = $$gene{'histone'};
  push (@{$types{$histone}}, $gene);

  ## to find the start and end of each cluster, we just list all the start and
  ## end coordinates. At the end, we get the min and max of them.
  push (@{$canon{$cluster}{"coordinates"}}, $$gene{'start'}, $$gene{'end'});

  ## count histones (by cluster, type and coding/pseudo)
  $canon{$cluster}{"total"}++;
  $canon{$cluster}{$histone}{"total"}++;
  if ($$gene{'pseudo'}) {
    $canon{$cluster}{"pseudo"}++;
    $canon{$cluster}{$histone}{"pseudo"}++;
  } else {
    $canon{$cluster}{"coding"}++;
    $canon{$cluster}{$histone}{"coding"}++;
  }

  ## Get the locus.
  ## It is not possible to calculate it from the genomic coordinates (but
  ## should be possible to make bp_genbank_ref_extractor do it). However, we
  ## can find the value in the features of the transcripts files. Because of
  ## that, this only works on coding genes.
  my @transcripts = keys %{$$gene{"transcripts"}};
  if (@transcripts) {
    my  $seq      = MyLib::load_seq("transcript", $transcripts[0], $path{sequences});
    my ($feature) = $seq->get_SeqFeatures("source");
    my ($locus)   = $feature->get_tag_values("map");
    if ($locus =~ m/[\d]+[qp][\d\.]+/) {
      push (@{$canon{$cluster}{'locus'}}, $locus);
    } else {
      warn ("Could not find locus for $symbol in $transcripts[0]");
    }
  }
}

## Sort the histone by types, then sort them by their symbol, to
## make the table with all of the canonical histones, their gene symbols,
## UIDs, and protein and transcript accession numbers
my @sorted;
foreach my $histone (sort keys %types) {
  push (@sorted, sort {$$a{'symbol'} cmp $$b{'symbol'}} @{$types{$histone}});
}
my $tex_table_path = File::Spec->catdir($path{results}, "table-histone_catalogue.tex");
my $csv_table_path = File::Spec->catdir($path{results}, "table-histone_catalogue.csv");
MyLib::make_tex_catalogue ($tex_table_path, @sorted);
MyLib::make_csv_catalogue ($csv_table_path, @sorted);


my $stats_path = File::Spec->catdir($path{results}, "variables-cluster_stats.tex");
open (my $stats, ">", $stats_path) or die "Could not open $stats_path for writing: $!";

## Get the counts and stats for each of the clusters
my $formatter = Number::Format->new(-decimal_digits => $MyVar::size_precision);
my %totals;
foreach my $cluster_k (keys %canon) {
  my $cluster = $canon{$cluster_k};
  ## Calculate the length (in bp) of each cluster
  my $coord_start  = List::Util::min (@{$$cluster{'coordinates'}});
  my $coord_end    = List::Util::max (@{$$cluster{'coordinates'}});

  my $coord_length = $formatter->format_bytes (abs ($coord_start - $coord_end),
    base => 1000
  ) . "bp";
  say {$stats} MyLib::latex_newcommand ($cluster_k."Span", $coord_length);

  ## Get a nice LaTeX string showing the range of locus for each cluster,
  ## e.g, 6p21.3--6p22.2. The problem is that some locus have less precision
  ## than others. For example, 1q21 does not mean 1q21.0, it only means
  ## somewhere in 1q21. While it is tempting to just ignore it, it is not
  ## correct, we should use the one with lowest precision. Luckily, sorting
  ## seems to already do that for us. We are left with the problems of
  ## having some genes in the short and long arm of a chromosome but a cluster
  ## should not be spanning the two arms.
  my @locus = @{$$cluster{'locus'}};
  if (@locus == 0) {
    warn ("No locus information for cluster $cluster_k.");
  } else {
    my $locus_start = List::Util::minstr (@locus);
    my $locus_end   = List::Util::maxstr (@locus);
    my $locus = $locus_start eq $locus_end ?
                $locus_start : "$locus_start--$locus_end";
    say {$stats} MyLib::latex_newcommand ($cluster_k."Locus", $locus);
  }

  ## Some clusters may not have any coding or pseudo gene in which case
  ## these values were never initialized.
  $$cluster{'coding'} //= 0;
  $$cluster{'pseudo'} //= 0;

  say {$stats} MyLib::latex_newcommand ("CodingGenesIn$cluster_k", $$cluster{'coding'});
  say {$stats} MyLib::latex_newcommand ("PseudoGenesIn$cluster_k", $$cluster{'pseudo'});
  say {$stats} MyLib::latex_newcommand ("TotalGenesIn$cluster_k",  $$cluster{'total'});

  foreach my $histone (@MyVar::histones) {
    $$cluster{$histone}{"pseudo"} //= 0;
    $$cluster{$histone}{"coding"} //= 0;
    $$cluster{$histone}{"total"}  //= 0;
    say {$stats} MyLib::latex_newcommand ($histone."CodingIn$cluster_k", $$cluster{$histone}{'coding'});
    say {$stats} MyLib::latex_newcommand ($histone."PseudoIn$cluster_k", $$cluster{$histone}{'pseudo'});
    say {$stats} MyLib::latex_newcommand ($histone."TotalIn$cluster_k",  $$cluster{$histone}{'total'});
    $totals{pseudo} += $$cluster{$histone}{pseudo};
    $totals{coding} += $$cluster{$histone}{coding};
    $totals{total}  += $$cluster{$histone}{total};
  }
}

say {$stats} MyLib::latex_newcommand ("TotalGenes",       $totals{total});
say {$stats} MyLib::latex_newcommand ("TotalCodingGenes", $totals{coding});
say {$stats} MyLib::latex_newcommand ("TotalPseudoGenes", $totals{pseudo});

## Get the counts and stats for each of the histone types
foreach my $histone (@MyVar::histones) {
  my $coding = 0;
  my $pseudo = 0;
  foreach my $cluster(values %canon) {
    $coding += $$cluster{$histone}{"coding"};
    $pseudo += $$cluster{$histone}{"pseudo"};
  }
  say {$stats} MyLib::latex_newcommand ($histone."CodingGenes", $coding);
  say {$stats} MyLib::latex_newcommand ($histone."PseudoGenes", $pseudo);
  say {$stats} MyLib::latex_newcommand ($histone."TotalGenes", $coding + $pseudo);
}

close ($stats) or die "Couldn't close $stats_path after writing: $!";
