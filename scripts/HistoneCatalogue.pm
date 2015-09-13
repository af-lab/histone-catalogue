package HistoneCatalogue;
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

use 5.010;
use strict;
use warnings;
use Carp;

=var cluster_number
An integer value with the current number of known clusters.
=cut
our $cluster_number = 4;

## histones that we care about (in case one day we start caring about H1)
our @histones       = qw(H2A H2B H3 H4);
our $histone_regexp = join ("|", @histones);
## how to call bp_genbank_ref_extractor
our $seq_extractor  = 'bp_genbank_ref_extractor';
## max number of significant figures (digits) for sizes (cluster length)
our $size_precision = 2;
## LaTeX code to mark genes as pseudo on tables
our $pseudo_mark    = '($\psi$)';
## max distance since end of CDS and start of stem-loop
our $stlp_dist      = 70;

=var stlp_seq
A regexp to identify stem-loop.  Currently, according to PMID:17531405.
=cut
our $stlp_seq = 'GG[CT][CT]CTT[CT]T[CTA]AG[GA]GCC';

=var stlp_length
Known length of the stem-loop
=cut
our $stlp_length    = 16;

=var tex_macro_name
The name of the LaTeX command that we should use for our automatically
generated values.  The idea is to control the appearance of these values
automatically from the LaTeX sources.
=cut
our $tex_macro_name = "ScriptValue";


=func get_sequences_date

Returns the date when sequences where retrieved.  This is parsed from
the log file of bp_genbank_ref_extractor.  Note that this date may differ
from the date when the sequences are analysed.

Params:
  fpath - path to the extractor.log by bp_genbank_ref_extractor

Returns:
  String with date.
=cut

sub get_sequences_date
{
  my $fpath = shift;
  open (my $data_log, "<", $fpath)
    or croak "Could not open $fpath for reading: $!";
  my $data_header = <$data_log>; # read the first line only
  close $data_log;

  $data_header =~ m/(?<=\[)([\d\-: ]+)(?=\])/;
  return $1;
}


=func mk_latex_string

Convert a string into a LaTeX usable string.  Basically, this means
escaping characters for LaTeX.  We have really basic needs for now.
Hopefully it will remain that way.

Params:
  string - A string to be converted for LaTeX.

Returns:
  The input string with escaped characters for use in LaTeX sources.
=cut

sub mk_latex_string
{
  my $string = shift;
  my %replace = (
    '&'   => '\&',
    '%'   => '\%',
    '$'   => '\$',
    '#'   => '\#',
    '_'   => '\_',
    '{'   => '\{',
    '}'   => '\}',
    '~'   => '\\textasciitilde{}',
    '^'  => '\\textasciicircum{}',
    '\\'  => '\\textbackslash{}',
  );
  $string =~ s/(&|%|\$|#|_|\{|\}|~|\^|\\)/$replace{$1}/g;
  return $string;
}


=func mk_latex_newcommand

Creates the LaTeX code that generates a new LaTeX command with that
value as argument to $HistoneCatalogue::tex_macro_name. Because of
TeX limitations on valid identifiers, commands must be alphabetic only.
Numeric characters are replaced by their capitalized English name.  Any
other character will cause it to fail.

Params:
  name - name of the command.  Numbers get replaced, e.g., the
    name "HIST1" will be converted to "HISTOne"
  value - value corresponding to the new command.

Returns:
  string - with \newcommand macro

Exception:
  dies if name has any non-alphanumeric character.
=cut

sub mk_latex_newcommand
{
  my $name  = shift;
  my $value = mk_latex_string (shift);

  if ($name =~ m/[^a-z0-9]/i)
    { croak "Unable to form LaTeX command '$name': must be alphanumeric"; }

  ## Replace numeric characters with words.  This is not meant complete,
  ## there is Lingua::EN::Nums2Words, Lingua::EN::Numbers or Number::Spell
  ## for that.  Simply replace 0-9 is enough for us.
  my %trans = (
    0 => "Zero",
    1 => "One",
    2 => "Two",
    3 => "Three",
    4 => "Four",
    5 => "Five",
    6 => "Six",
    7 => "Seven",
    8 => "Eight",
    9 => "Nine",
  );
  $name =~ s/([0-9])/$trans{$1}/g;

  ## we need to set the height of the color box manually otherwise \colorbox
  ## will change the height of the line
  return "\\newcommand{\\$name}{\\${tex_macro_name}{$value}}";
}


=func say_latex_newcommand

Writes to a file a new LaTeX command with documentation as a LaTeX
comment.  Command names must be alphanumeric (see mk_latex_newcommand
for details).

Params:
  file - file open to write.
  name - a string.  See mk_latex_newcommand
  value -  a string.  See mk_latex_newcommand
  docs - description of the command.  It will be added as a comment
    on the previous line for documentation of the generated sources.
    Defaults to "Not documented"

Returns:
  void
=cut

sub say_latex_newcommand
{
  my $file  = shift;
  my $name  = shift;
  my $value = shift;
  my $docs  = shift || "Not documented";

  say {$file} "%% $_" for split ("\n", $docs);
  say {$file} mk_latex_newcommand ($name, $value);
}

1;
