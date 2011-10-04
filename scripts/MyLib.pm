package Local::MODULE_NAME;
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

## setting up Exporter
our (@ISA, @EXPORT, @EXPORT_OK, $VERSION);  # must be package global variables for Exporter
use Exporter;                               # handle module's external interface
$VERSION    = 0.9;                          # set module version number
@ISA        = qw(Exporter);                 # inherit import method from Exporter
@EXPORT     = qw();                         # export nothing automatically
@EXPORT_OK  = qw();                         # export only these and by request


1;    # last line. A module must return a true value so that 'require' and 'use succeeds
