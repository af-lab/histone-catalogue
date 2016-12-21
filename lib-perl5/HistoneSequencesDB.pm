package HistoneSequencesDB;
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

=head1 DESCRIPTION

Handles the directory where bp_genbank_ref_extractor saved sequences,
and log files.  This is the path specified via the "--save" option.

We drop any sequence that we do not recognize as an histone.  Maybe
we shouldn't?

The genes property is at the core of this object, and it's an array
of CanonicalHistoneGene and HistoneGene objects.  This may also not
be the best of designs.

=cut

use strict;
use warnings;
use Carp;

use Storable;
use File::Spec;
use List::Util;

use Text::CSV 1.21; # require 1.21 for getline_hr_all

use Moose;
use Moose::Util::TypeConstraints qw(subtype as where);

use Bio::SeqIO;

use CanonicalHistoneGene;
use HistoneGene;
use Gene;

use namespace::autoclean;

subtype 'Dir',
  as 'Str',
  where { -d $_ };

has 'dir' =>
(
  is => 'ro',
  isa => 'Dir',
  required => 1
);

has 'genes' => (
  is => 'ro',
  isa => 'ArrayRef[Gene]',
  builder => '_build_genes_from_csv',
  init_arg => undef,
);

has '_canonical_core_idx' => (
  is        => 'ro',
  isa       => 'ArrayRef',
  lazy      => 1,
  builder   => '_build_canonical_core_idx',
  init_arg  => undef,
);
has '_linkers_idx' => (
  is        => 'ro',
  isa       => 'ArrayRef',
  lazy      => 1,
  builder   => '_build_linkers_idx',
  init_arg  => undef,
);
has '_variants_core_idx' => (
  is        => 'ro',
  isa       => 'ArrayRef',
  lazy      => 1,
  builder   => '_build_variants_core_idx',
  init_arg  => undef,
);

has 'log_path' =>
(
  is => 'ro',
  lazy => 1,
  default => "extractor.log",
);

has 'csv_path' =>
(
  is => 'ro',
  lazy => 1,
  default => "data.csv",
);

has 'genes_dir' =>
(
  is => 'ro',
  lazy => 1,
  default => "genes",
);

has 'proteins_dir' =>
(
  is => 'ro',
  lazy => 1,
  default => "proteins",
);

has 'transcripts_dir' =>
(
  is => 'ro',
  lazy => 1,
  default => "transcripts",
);

around BUILDARGS => sub
{
  my $orig = shift;
  my $class = shift;

  if (@_ == 1 && ! ref $_[0])
    { return $class->$orig(dir => $_[0]); }
  else
    { return $class->$orig(@_); }
};


=func _csv_line_2_Gene
Creates a new Gene object from a Text::CSV line (an hash ref with the
csv headers as keys).

Args:
  entry (HashRef)

Returns:
  Gene

Exception:
  dies if failing to create a Gene object (e.g., a line for a coding
    gene without products).
=cut
sub _csv_line_2_Gene
{
  my $entry = shift;

  my $uid = $entry->{'gene UID'};
  my $symbol = $entry->{'gene symbol'};

  my $gene_ctor;
  if ($symbol =~ m/^HIST(\d+)/i) # a canonical histone (belongs to a cluster)
    { $gene_ctor = sub { CanonicalHistoneGene->new(@_); }; }
  else # must be an histone variant
    {
      my $histone_type = HistoneGene::symbol2type($symbol);
      if (! $histone_type)
        { croak "Unable to find histone variant type on '$symbol'"; }
      $gene_ctor = sub { HistoneGene->new('histone_type' => $histone_type, @_); };
    }

  my $type;
  my %products = ();
  if ($entry->{'pseudo'})
    { $type = 'pseudo'; }
  else
    {
      $type = 'coding';
      my $nm_acc = $entry->{'transcript accession'};
      if ($nm_acc)
        { $products{$nm_acc} = $entry->{'protein accession'}; }
      else
        { croak ("Coding gene $symbol entry without a transcript"); }
    }

  my $gene = &$gene_ctor (
    'uid' => $uid,
    'symbol' => $symbol,
    'type' => $type,
    'species' => $entry->{"species"},
    'description' => $entry->{'gene name'},
    'chr_acc' => $entry->{'chromosome accession'},
    'chr_start' => $entry->{'chromosome start coordinates'},
    'chr_end' => $entry->{'chromosome stop coordinates'},
    'products' => \%products,
  );
  return $gene;
}

=func _updated_Gene_with_csv_line

Args:
  entry (HashRef)

Returns:
  Gene

Exception:
  dies if failing to create a Gene object (e.g., a line for a coding
    gene without products).
=cut
sub _updated_Gene_with_csv_line
{
  my $old_gene = shift;
  my $entry = shift;

  my $products = $old_gene->products;
  my $nm_acc = $entry->{'transcript accession'};

  if ($nm_acc)
    { $products->{$nm_acc} = $entry->{'protein accession'}; }
  else
    ## Hopefully we will never get here
    { croak 'Found two entries for ' . $old_gene->symbol . ' but no transcript-protein pair'; }

  my @ctor_args;

  ## If we were dealing a non-canonical histone gene, we must
  ## specify the histone type on the constructor.
  if (! $old_gene->isa("CanonicalHistoneGene"))
    { push (@ctor_args, 'histone_type' => $old_gene->histone_type); }

  my $gene = $old_gene->new
    (
      'uid' => $old_gene->uid,
      'symbol' => $old_gene->symbol,
      'type' => $old_gene->type,
      'species' => $old_gene->species,
      'description' => $old_gene->description,
      'chr_acc' => $old_gene->chr_acc,
      'chr_start' => $old_gene->chr_start,
      'chr_end' => $old_gene->chr_end,
      'products' => $products,
      @ctor_args,
    );
  return $gene;
}


=method _build_genes_from_csv
Read genes from the csv file create by bp_genbank_ref_extractor to build
this object "genes" property.

Returns:
  Array of Gene objects.
=cut
sub _build_genes_from_csv
{
  my $self = shift;
  my $csv_path = File::Spec->catfile($self->dir, $self->csv_path);

  ## To cover the widest range of parsing options, you will always
  ## want to set binary
  my $csv = Text::CSV->new ({binary => 1, eol => $/})
    or croak "Cannot use Text::CSV: ". Text::CSV->error_diag ();

  open (my $file, "<", $csv_path)
    or croak "Could not open '$csv_path' for reading: $!";

  ## read all into array of hashes with headers as keys
  $csv->column_names ($csv->getline ($file));
  my $data = $csv->getline_hr_all ($file);
  close $file;

  ## We fill the genes as we read the csv file, one line at a time.
  ## A gene may take more than one csv entry.  So we update gene
  ## entries when.  However, some entries on their own do not form
  ## a valid gene.  Our current case is a coding gene with one protein
  ## and another non-coding transcript.  If the first entry of the
  ## gene is the non-coding transcript, the Gene constructor will fail
  ## (a coding gene must have a protein product).  For those cases, we
  ## add the csv entries to this pending hash.
  my %pending; # keys are gene UID, values are {'entries' => \@csv_lines, 'err' => '$@'}

  my %genes;
  foreach my $entry (@$data)
    {
      my $symbol = $entry->{'gene symbol'};
      my $histone_type = HistoneGene::symbol2type($symbol);
      next unless $histone_type;

      my $uid = $entry->{'gene UID'};

      if (! $entry->{'chromosome accession'})
        {
          carp ("Gene with UID '$uid' has no genomic information. Skipping it!");
          next;
        }

      if (grep { $_ == $uid} keys %pending)
        {
          my $gene;
          eval { $gene = _csv_line_2_Gene($entry); };
          if ($@)
            { push (@{ $pending{$uid}->{entries} }, $entry); }
          else
            {
              my @failed_entries = @{ $pending{$uid}->{entries} };
              for (@failed_entries)
                { $gene = _updated_Gene_with_csv_line ($gene, $_); }
              $genes{$uid} = $gene;
              delete $pending{$uid};
            }
        }
      elsif (exists $genes{$uid})
        {
          my $old_gene = $genes{$uid};
          my $new_gene = _updated_Gene_with_csv_line($old_gene, $entry);
          $genes{$uid} = $new_gene;
        }
      else
        {
          my $gene;
          eval { $gene = _csv_line_2_Gene($entry); };
          if ($@)
            {
              $pending{$uid}->{err} = $@;
              $pending{$uid}->{entries} = [$entry];
            }
          else
            { $genes{$uid} = $gene; }
        }
    }
  if (%pending)
    { carp 'Unable to form Gene for UIDs '. join (" ", keys %pending); }

  return [values %genes];
}


sub _build_canonical_core_idx
{
  my $self = shift;
  my $genes = $self->genes;
  my @idx = grep {$genes->[$_]->isa ('CanonicalHistoneGene')
                  && $genes->[$_]->is_core_histone} 0 .. $#{$genes};
  return \@idx;
}
sub _build_linkers_idx
{
  my $self = shift;
  my $genes = $self->genes;
  my @idx = grep {$genes->[$_]->is_linker_histone} 0 .. $#{$genes};
  return \@idx;
}
sub _build_variants_core_idx
{
  my $self = shift;
  my $genes = $self->genes;
  my @idx = grep {! $genes->[$_]->isa ('CanonicalHistoneGene')
                  && $genes->[$_]->is_core_histone} 0 .. $#{$genes};
  return \@idx;
}


=method write_db

Save a HistoneSequencesDB object into file for later retrieval via
HistoneSequencesDB::read_db()

Args:
  fpath (string) - path to file to be saved.

Returns:
  void

=cut
sub write_db
{
  my $self = shift;
  my $store_fpath = shift;
  Storable::store ($self, $store_fpath)
    or croak "Unable to store HistoneSequencesDB in '$store_fpath'";
  return;
}

=func read_db

Reads a HistoneSequencesDB object from a file and returns it.  This
file is expected to be a HistoneSequencesDB object saved via
HistoneSequencesDB::write_db().

Args:
  fpath (string) - file to read.

Returns:
  HistoneSequencesDB object

=cut
sub read_db
{
  my $store_fpath = shift;
  my $db = Storable::retrieve($store_fpath);
  if (! defined $db)
    { croak "Unable to retrieve from HistoneSequencesDB from '$store_fpath'"; }
  if (! $db->isa("HistoneSequencesDB"))
    { croak "Retrieve object from '$store_fpath' is not an HistoneSequencesDB"; }
  return $db;
}


=method canonical_core

Return an array with all canonical and core histone genes.
=cut
sub canonical_core
{
  my $self = shift;
  return @{$self->genes}[@{$self->_canonical_core_idx}];
}

=method linkers

Return an array with all linker histone genes.  This includes
both canonical and variant linkers.
=cut
sub linkers
{
  my $self = shift;
  return @{$self->genes}[@{$self->_linkers_idx}];
}

=method variants_core

Return an array with all variant core histone genes.
=cut
sub variants_core
{
  my $self = shift;
  return @{$self->genes}[@{$self->_variants_core_idx}];
}


=method _get_seq

Args:
  sub_dir
  fname

Returns
  Bio::Seq
=cut
sub _get_seq
{
  my $self = shift;
  my $sub_dir = shift;
  my $fname = shift;
  my $path = File::Spec->catdir($self->dir, $sub_dir, "$fname.gb");
  ## We know there's only one sequence in those genbank files.
  my $seq = Bio::SeqIO->new(-file => $path)->next_seq;
  return $seq;
}

=method get_genomic
Args:
  uid - a string with the gene uid.

Returns:
  Bio::Seq object for the gene.
=cut
sub get_genomic
{
  my $self = shift;
  my $uid = shift;
  return $self->_get_seq($self->genes_dir, $uid);
}

=method get_transcript
Args:
  acc - a string with the transcript accession number.

Returns:
  Bio::Seq object for the transcript.
=cut
sub get_transcript
{
  my $self = shift;
  my $acc = shift;
  return $self->_get_seq($self->transcripts_dir, $acc);
}

=method get_transcript_cds
Args:
  acc - a string with the transcript accession number.

Returns:
  Bio::Seq object for the transcript CDS.
=cut
sub get_transcript_cds
{
  my $self = shift;
  my $acc = shift;
  my $seq = $self->_get_seq($self->transcripts_dir, $acc);
  my $cds = ($seq->get_SeqFeatures("CDS"))[0]->seq();

  ## Remove the start codon, for the same reasons we remove the methionine
  ## in get_protein().  We want the sequences to match.
  ##
  ## In addition, also remove the stop codon (so it matches with the
  ## protein sequence).  I'm not 100% sure this is correct but matches
  ## our current needs.
  $cds = $cds->trunc(4, $cds->length -3);

  return $cds;
}

=method get_protein
Args:
  acc - a string with the protein accession number.

Returns:
  Bio::Seq object for the protein.
=cut
sub get_protein
{
  my $self = shift;
  my $acc = shift;
  my $seq = $self->_get_seq($self->proteins_dir, $acc);

  ## XXX  it is standard to remove the initial methionine when dealing
  ##      with histone proteins (it does not even count when giving position
  ##      in the protein). This is really really not recommended by the HGVS
  ##      (one would think is also common sense), but the number of the
  ##      amino acids is so ingrained in the field that we can't start
  ##      now to give them other numbers (everyone knows H3 K4Me3, we can't
  ##      start to correct them and call it H3 K5Me3 because we are not
  ##      important enough to break the convention).  Comment and Uncomment
  ##      the following line, to keep or remove the N-terminal methionine

  ## Remove the first amino acid since in histones it's cleaved off
  $seq = $seq->trunc(2, $seq->length);

  return $seq;
}


=func sort_histones

Sorts an array of histone genes in the most expected order, i.e., by
histone type and then by gene symbol.

Args:
  Array of HistoneGene

Returns
  Array of HistoneGene
=cut
sub sort_histones
{
  return sort {$a->histone_type cmp $b->histone_type
               || $a->symbol cmp $b->symbol} @_;
}


=method protein_iterator

Given a list of genes, creates an iterator for all its proteins.  Note
that some genes may have more than one, or even zero proteins.  So the
iterator may return more or less proteins than the number of genes.

Args:
  @genes ([Gene])

Returns:
  iterator for Bio::Seq protein
  undef when iterator is exhausted
=cut
sub protein_iterator
{
  my $db = shift;
  my @genes = @_;

  my @acc;
  for my $g (@genes)
    {
      my $products = $g->coding_products();
      push (@acc, values %{$products});
    }

  return sub
    {
      while (my $p_acc = shift @acc)
        { return $db->get_protein($p_acc); }
      return undef;
    };
}


__PACKAGE__->meta->make_immutable;

1;
