package WebLogo;
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

use strict;
use warnings;
use Carp;

use File::Copy;
use File::Temp;
use File::Which;

use Moose;
use Moose::Util::TypeConstraints qw(subtype as where);

use namespace::autoclean;

# ABSTRACT: Simple interface to WebLogo http://weblogo.threeplusone.com/

## If this module is used, then the weblogo paper should be cited:
##
## Crooks, GE, Hon, G, Chandonia, JM, Brenner SE (2004). "WebLogo: a sequence
## logo generator", Genome Research, 14:1188-1190
##
##@article{weblogo2004,
##  title={{WebLogo}: a sequence logo generator},
##  author={Crooks, G.E. and Hon, G. and Chandonia, J.M. and Brenner, S.E.},
##  journal={Genome research},
##  volume={14},
##  number={6},
##  pages={1188--1190},
##  year={2004},
##}


subtype 'Bin',
  as 'Str',
  where { -f $_ && -x _ };

has ['bin_path'] =>
(
  is => 'ro',
  isa => 'Bin',
  lazy => 1,
  default => sub { File::Which::which("weblogo"); },
);

has ['default_args'] =>
(
  is => 'ro',
  isa => 'HashRef',
  lazy => 1,
  default => sub { {"--format" => "eps",
                    "--fineprint" => "",
                    "--color-scheme" => "monochrome",
                    }
                 },
);

has ['version'] =>
(
  is => 'ro',
  isa => 'Str',
  init_arg => undef,
  lazy => 1,
  builder => '_get_version',
);

=method call
Call weblogo with specific options

Args:
  fin (string) - input file with alignment file
  fout (string) - path for output file
  extra_args (HashRef) - extra options (possibly update from the object
    default_args)

Exception:
  dies on failing system call

Returns:
  void
=cut
sub call
{
  my $self = shift;
  my $fin = shift;
  my $fout = shift;
  my $extra_args = shift;

  my %args = (%{$self->default_args}, %$extra_args);
  system ($self->bin_path, %args, "--fin", $fin, "--fout", $fout) == 0
    or croak "Call to weblogo failed: $?";
}

=method _get_version
Get version of this weblogo.

Returns:
  string with WebLogo version identifier.

=cut
sub _get_version
{
  my $self = shift;
  my $bin = $self->bin_path;
  ## We expect something like "WebLogo 3.4 (2014-06-02)"
  `$bin --version` =~ m/^WebLogo ([0-9\.]+) \(/;
  if (! defined $1)
    { croak "Unable to identify version of weblogo"; }
  return $1;
}


=func seqlogo_2_fancy_align
Modify sequence logo (eps file) into a fancy looking alignment output.

A sequence logo displays a lot of interesting stuff.  However, we
really only want to see an alignment.  So we modify the sequence logo
(which is already a simple: y axis is probability and no error bars)
so that the sequence is on light gray except when all it's the same on
all sequences.  The end result is a mostly light gray sequence, with
only the locations where there's a variant sequence on black.

Args:
  fin - filepath for input seqlogo eps file
  fout - filepath for output seqlogo eps file

Exception:
  dies if fails to open or write the eps files, or if the eps file
    does not match what is expected (maybe a different version of
    weblogo.
=cut
sub seqlogo_2_fancy_align
{
  my $fin = shift;
  my $fout = shift;

  ## It seems to me that the content of postscript files generated by
  ## weblogo are (in this order):
  ##
  ##    1 - options set as variable at the top;
  ##    2 - followed by a piece of code common to all usage;
  ##    3 - the stacks where all characters and height values are.
  ##
  ## Technically, we could create the file ourselves, no need for weblogo.
  ## We would just copy the start of one, and use it as template.  The
  ## stacks should be trivial to create.  Or we could even use weblogo to
  ## create the stacks and replace everything that comes before before.
  ## But that's just duplicating code which should be avoided.  It means
  ## that changes on their side can break ours but we should be able to run
  ## a test at run time to get an alternative "fix".
  ##
  ## So we process the postscript file and change the color used for the
  ## letter depending on its height.  This is done in the operator DrawChar
  ## (called from ShowSymbol).  We change the color right before the call
  ## to show() to minimize problems if something changes in the future.
  ## There is only line with "tc show" in the whole file.  It would be easy
  ## to make sure we are in the DrawChar operator but if the code changes
  ## so much that there's two "tc show" in the file, then it's time for us
  ## to update this function.

  ## FIXME: in the case of amino acids that are only present in some of the
  ##        proteins, their height will be the height of the stack so not in
  ##        black.  We probably should highlight those as well.

  my $fh = File::Temp->new(UNLINK => 0);

  open (my $read, "<", $fin)
    or croak "Couldn't open $fin for reading: $!\n";

  my $fixed = 0;
  while (my $line = <$read>)
    {
      if ($line =~ m/tc show/)
        {
          ## There should be only one 'tc show' ocurrence on each file.
          ## If this is not true, something's wrong (maybe a different
          ## version of weblogo).
          if ($fixed)
            { croak "Trying to fix previously fixed sequence Logo. Probably new version of weblogo and we need to update our script."; }
          print {$fh} "        ysize stack_height stack_margin sub ge {\n" .
                      "            0.7 0.7 0.7 setrgbcolor\n" .
                      "        } {\n" .
                      "            0.0 0.0 0.0 setrgbcolor\n" .
                      "        } ifelse\n";
          $fixed = 1;
        }
      print {$fh} $line;
    }
  close($read)
    or croak "Couldn't close $fin after reading: $!";

  close($fh)
    or croak "Couldn't close $fh after reading: $!";
  File::Copy::move("$fh", $fout);
}


__PACKAGE__->meta->make_immutable;

1;