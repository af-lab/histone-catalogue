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

use strict;
use warnings;

use Test::More;
use File::Temp;

use HistoneCatalogue;

## To test functions whose first argument is an open file, and we want
## to test what gets written to it.
sub test_with_perlio
{
  my $foo = shift;
  my $file_contents;
  open (my $fh, '>', \$file_contents)
    or die "Can't open variable to write: $!";
  &{$foo}($fh, @_);
  close $fh;
  return $file_contents;
}

## Create temporary file with the contents of a variable.  The first
## argument to the function must be the expected file.
sub test_with_tmp_file
{
  my $foo = shift;
  my $content = shift;
  my $tmp = File::Temp->new(UNLINK => 1);
  open (my $fh, ">", $tmp->filename)
    or die "Can't open variable to write: $!";
  print {$fh} $content;
  close $fh;
  return &{$foo}($tmp->filename, @_);
}


ok (test_with_tmp_file (\&HistoneCatalogue::get_sequences_date, <<END)
This is bp_genbank_ref_extractor on Bioperl 1.006924 on [2015-09-13 14:10:39]
Entrez gene: searching with '"homo sapiens"[organism] (H2A*[gene name] OR H2B*[gene name])'
Entrez gene: query "homo sapiens"[organism] (H2A*[gene name] OR H2B*[gene name])"
END
    eq '2015-09-13');


is (HistoneCatalogue::mk_latex_string ('foobar'), 'foobar',
   "Escape (not) simple string for LaTeX");
is (HistoneCatalogue::mk_latex_string ('_foo_bar_'), '\_foo\_bar\_',
   "Escape underscores for LaTeX");
is (HistoneCatalogue::mk_latex_string ('$foobar'), '\\$foobar',
   "Escape dollar sign for LaTeX");
is (HistoneCatalogue::mk_latex_string ('^$foobar^'),
    '\\textasciicircum{}\\$foobar\\textasciicircum{}',
    'Escape ^ and $ for LaTeX');
is (HistoneCatalogue::mk_latex_string ('&%$#_{}~^\\' x2),
    '\\&\\%\\$\\#\\_\\{\\}\\textasciitilde{}\\textasciicircum{}\\textbackslash{}' x2,
    "Escape all escapable characters in string for LaTeX");
is (HistoneCatalogue::mk_latex_string ("\nfoo\nbar\n"),
    "\nfoo\nbar\n",
    "Handle newlines when escaping (not) characters for LaTeX");
is (HistoneCatalogue::mk_latex_string ("\nfoo^\nbar\n"),
    "\nfoo\\textasciicircum{}\nbar\n",
    "Handle newlines while escaping characters for LaTeX");


is ($HistoneCatalogue::tex_macro_name, 'ScriptValue',
    'Macro name for results is ScriptValue');

is (HistoneCatalogue::latex_newcommand ("foo", "bar"),
    "%% Not documented\n\\newcommand{\\foo}{\\ScriptValue{bar}}");
is (HistoneCatalogue::latex_newcommand ('HIST1H2A', '67%_p'),
    "%% Not documented\n\\newcommand{\\HISTOneHTwoA}{\\ScriptValue{67\\%\\_p}}");

is (HistoneCatalogue::latex_newcommand ("foo", "bar", "This is some serious documentation."),
    "%% This is some serious documentation.\n"
    . "\\newcommand{\\foo}{\\ScriptValue{bar}}");
is (HistoneCatalogue::latex_newcommand ("f0", "67%", "This is some serious\nMultiline documentation."),
    "%% This is some serious\n%% Multiline documentation.\n"
    . "\\newcommand{\\fZero}{\\ScriptValue{67\\%}}");

is (HistoneCatalogue::mk_latex_list_name_isoforms ("H2A", ("HIST2H2AF", "HIST1H2AB", "HIST1H2AC", "HIST1H2AP", "HIST4H2A")),
    "HIST1H2A --B, --C, --P; HIST2H2AF; HIST4H2A");
is (HistoneCatalogue::mk_latex_list_name_isoforms ("H4", ("HIST4H4")),
    "HIST4H4");
is (HistoneCatalogue::mk_latex_list_name_isoforms ("H2A", ("HIST2H2AA3", "HIST1H2AB", "HIST1H2AC", "HIST1H2AP", "HIST4H2A")),
    "HIST1H2A --B, --C, --P; HIST2H2AA3; HIST4H2A");
is (HistoneCatalogue::mk_latex_list_name_isoforms ("H4", ("HIST4H4")),
    "HIST4H4");

done_testing;
