package MyLib;
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

use 5.010;                                  # use Perl 5.10
use strict;                                 # enforce some good programming rules
use warnings;                               # replacement for the -w flag, but lexically scoped
use Carp;                                   # alternative warn and die for modules
use Text::CSV 1.21;                         # Comma-separated values manipulator (require 1.21 for getline_hr_all)
use POSIX;                                  # Perl interface to IEEE Std 1003.1
use FindBin;                                # Locate directory of original perl script

use lib $FindBin::Bin;                      # Add script directory to @INC to find 'package'
use MyVar;                                  # Load variables

## load the gene information from all genes found
sub load_csv {
  ## To cover the widest range of parsing options, you will always want to set binary
  my $csv = Text::CSV->new ({
                              binary => 1,
                              eol    => $/,
                              }) or die "Cannot use Text::CSV: ". Text::CSV->error_diag ();
  open (my $file, "<", $MyVar::data_path) or die "Could not open $MyVar::data_path for reading: $!";

  $csv->column_names ($csv->getline ($file));   # read first line and sets it as the column name

  ## note that get_line_hr_all was only implemented on 1.21. If using 1.18, would
  ## need a while loop and use get_line_hr
  my $data_ref = $csv->getline_hr_all ($file);  # reads all lines of file into an array of hashes (returns ref to array)
  close $file;                                  # close file
  return @$data_ref;                            # dereference the array and return it
}

## rather than load information from all genes found and extracted, get only
## the canonical histones
sub load_canonical {
  my @data = load_csv;
  my @canon;
  foreach my $gene (@data) {
      my $symbol = $$gene{'gene symbol'};

      ## skip genes that don't look canonical and get cluster number
      next unless $symbol =~ m/^HIST(\d+)/;

      ## warn if a gene is found whose nomeclature mentions an unknown cluster
      if ($1 > $MyVar::cluster_number) {
        warn ("Update/Check the code, found possible NEW histone cluster $1 with gene '$symbol'");
      }

      ## skip genes without genomic information
      if ( !$$gene{'chromosome accession'}) {
        warn ("Gene '$symbol' has no genomic information. Skipping it!");
        next;
      }

    push (@canon, $gene);
  }
  return @canon;
}

## turn a distance in base pairs (a positive integer) into a string appropriate
## for text (in LaTeX format), e.g. 1000000 -> 1\,Mbp; 1540 -> 1.54\,kbp (the
## precision is defined in MyVar.pm
sub pretty_length {
  my $length   = $_[0];
  ## (ceil() +1 ) because when /3, if it's 3, floor will give 1 when we want 0
  my $power    = ( ceil(length($length) / 3) -1) * 3;
  my $dec_case = 0;
  my $number   = sprintf("%1.${dec_case}f", $length / (10 ** $power) );
  if ( length($length) < $MyVar::size_precision) {
    ## nothing, no decimal cases at all in these cases
  } elsif (length($number) < $MyVar::size_precision) {
    my $dec_case = $MyVar::size_precision - length ($number);
    $number      = sprintf("%1.${dec_case}f", $length / (10 ** $power) );
  }
  my $prefix;
  given ($power) {
    when  (0) { $prefix = ''  }
    when  (3) { $prefix = 'k' }
    when  (6) { $prefix = 'M' }
    when  (9) { $prefix = 'G' }
    when (12) { $prefix = 'T' }
    when (15) { $prefix = 'P' }
    when (18) { $prefix = 'E' }
    when (21) { $prefix = 'Z' }
    when (24) { $prefix = 'Y' }
    default   { $power -= 24; $prefix = "e+${power}Y" }
  }
  return "$number\\,${prefix}bp";
}

## make LaTeX tables (we probably should use one of the LaTeX modules in CPAN
## but our needs are so simple). First argument is the type of row, followed
## by the values for each column. Possible values are:
##    * start  - start table environment (second argument string with column alignments)
##    * header - values for header
##    * row    - values for each row
##    * end    - only closes table enviroment (use `row' to enter values for last row)
sub latex_table {
  my $tab;
  my $cols = @_ - 1;
  given ($_[0]) {
    when ("row") {
      $tab .= "  ";
      $tab .= "$_ & " for @_[1 .. $#_ -1];
      $tab .= "$_[-1] \\\\";
    }
    when ("start") {
      $tab .= "\\begin{tabular}{$_[1]}";
    }
    when ("header") {
      $tab .= latex_table ("row", @_[1 .. $#_]);
      $tab .= "\n";
      $tab .= '  \hline';
    }
    when ("end") {
      $tab .= '\end{tabular}';
      ## the end should be ONLY to close the environment. Let's be nice and warn
      ## if someone tries to give more arguments
      warn "MyLib::latex_table: ignoring some arguments. Do not use end for the last row" if $cols > 0;
    }
    default         { die "Unrecognized value $_[0] for MyLib::latex_table" }
  }
  return $tab;
}

1; # a package must return true
