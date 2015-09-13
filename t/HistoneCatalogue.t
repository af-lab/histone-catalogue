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

use Test::More tests => 7;
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


ok (HistoneCatalogue::mk_latex_string ('foobar')
    eq 'foobar');
ok (HistoneCatalogue::mk_latex_string ('_foo_bar_')
    eq '\_foo\_bar\_');
ok (HistoneCatalogue::mk_latex_string ('$foobar')
    eq '\\$foobar');
ok (HistoneCatalogue::mk_latex_string ('^$foobar^')
    eq '\\textasciicircum{}\\$foobar\\textasciicircum{}');
ok (HistoneCatalogue::mk_latex_string ('&%$#_{}~^\\' x2)
    eq '\\&\\%\\$\\#\\_\\{\\}\\textasciitilde{}\\textasciicircum{}\\textbackslash{}' x2);
ok (HistoneCatalogue::mk_latex_string ("\nfoo\nbar\n")
    eq "\nfoo\nbar\n");
ok (HistoneCatalogue::mk_latex_string ("\nfoo^\nbar\n")
    eq "\nfoo\\textasciicircum{}\nbar\n");


ok ($HistoneCatalogue::tex_macro_name eq 'ScriptValue');


ok (HistoneCatalogue::mk_latex_newcommand ('foo', 'bar')
    eq '\\newcommand{\\foo}{\\ScriptValue{bar}}');
ok (HistoneCatalogue::mk_latex_newcommand ('HIST1H2A', '67%_p')
    eq '\\newcommand{\\HISTOneHTwoA}{\\ScriptValue{67\%\_p}}');


ok (test_with_perlio (\&HistoneCatalogue::say_latex_newcommand, "foo", "bar")
    eq "%% Not documented\n\\newcommand{\\foo}{\\ScriptValue{bar}}\n");
ok (test_with_perlio (\&HistoneCatalogue::say_latex_newcommand, "foo", "bar",
                      "This is some serious documentation.")
    eq "%% This is some serious documentation.\n"
       . "\\newcommand{\\foo}{\\ScriptValue{bar}}\n");
ok (test_with_perlio (\&HistoneCatalogue::say_latex_newcommand, "f0", "67%",
                      "This is some serious\nMultiline documentation.")
    eq "%% This is some serious\n%% Multiline documentation.\n"
       . "\\newcommand{\\fZero}{\\ScriptValue{67\\%}}\n");

