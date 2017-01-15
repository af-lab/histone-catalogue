#!/usr/bin/perl
## Copyright (C) 2013 Carnë Draug <carandraug+dev@gmail.com>
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
use Text::CSV;                  # Comma-separated values manipulator
use File::Spec;                 # Perform operation on file names
use Bio::Seq;                   # Sequence object, with features
use Bio::SeqIO;                 # Handler for SeqIO formats
use FindBin;                    # Locate directory of original perl script
use Bio::DB::EUtilities;        # Retrieve entries from Entrez
use List::Util;                 # General-utility list subroutines

## This script builds the directory structure and CSV file with the data from
## our reference (Marzluff paper from 2002), as if it was downloaded with
## our tool. By doing this, our code to compare the latest data with this
## reference, can be used to compare against our own data later, or to see
## how much it changed over time.

my $email = ''; # used to download sequences from Genbank

## Because this is meant to be only used once anyway, we don't remove the
## directory if it's already there. Actually, exiting with an error means
## we are doing something funny. Remove it manually if that's really what
## you want to do.
my $proteins_dir = File::Spec->catdir($FindBin::Bin, "proteins");
mkdir ($proteins_dir) or die "Could not mkdir $proteins_dir: $!";

my %histone; # will hold all the histones

my %table1 = (
  ## from Table 1 (page 488 and 489)
  coding => {
    H2A => {
      Hist1h2aa  => "AY158921",
      Hist1h2ab  => "AY158920",
      Hist1h2ac  => "AY158919",
      Hist1h2ad  => "AY158918",
      Hist1h2ae  => "AY158917",
      Hist1h2af  => "AY158916",
      Hist1h2ag  => "AY158915",
      Hist1h2ah  => "AY158914",
      Hist1h2ai  => "AY158909",
      Hist1h2ak  => "AY158911",
      Hist1h2an  => "AY158912",
      Hist1h2ao  => "AY158913",
      Hist2h2aa1 => "AY158925",
      Hist2h2aa2 => "AY158924",
      Hist2h2ab  => "AY158922",
      Hist2h2ac  => "AY158923",
      Hist3h2a   => "AY158909",
    },
    H2B => {
      Hist1h2ba => "AY158939",
      Hist1h2bb => "AY158938",
      Hist1h2bc => "AY158937",
      Hist1h2be => "AY158936",
      Hist1h2bf => "AY158935",
      Hist1h2bg => "AY158934",
      Hist1h2bh => "AY158933",
      Hist1h2bj => "AY158932",
      Hist1h2bk => "AY158931",
      Hist1h2bl => "AY158927",
      Hist1h2bm => "AY158928",
      Hist1h2bn => "AY158929",
      Hist1h2bp => "AY158930",
      Hist2h2bb => "AY158940",
      Hist2h2be => "AY158941",
      Hist3h2ba => "AY158942",
      Hist3h2bb => "AY158943",
    },
    H3 => {
      Hist1h3a   => "AY158952",
      Hist1h3b   => "AY158951",
      Hist1h3c   => "AY158950",
      Hist1h3d   => "AY158949",
      Hist1h3e   => "AY158948",
      Hist1h3f   => "AY158947",
      Hist1h3g   => "AY158946",
      Hist1h3h   => "AY158944",
      Hist1h3i   => "AY158945",
      Hist2h3b   => "AY158955",
      Hist2h3ca1 => "AY158954",
      Hist2h3ca2 => "AY158953",

    },
    H4 => {
      Hist1h4a  => "AY158965",
      Hist1h4b  => "AY158964",
      Hist1h4c  => "AY158963",
      Hist1h4d  => "AY158962",
      Hist1h4f  => "AY158961",
      Hist1h4h  => "AY158960",
      Hist1h4i  => "AY158959",
      Hist1h4j  => "AY158956",
      Hist1h4k  => "AY158957",
      Hist1h4m  => "AY158958",
      Hist2h4   => "AY158966",
      Hist4h4   => "AY158967",
    },
  },

  ## from Table 1 (page 488 and 489)
  pseudo => {
    H2A => {
      Hist1h2aj  => "AY158910",
    },
    H2B => {},
    H3 => {},
    H4 => {},
  },
);

## H2A consensus sequence from Table 2: Histone H2A protein variants (page 494)
my $H2A_consensus = "SGRGKQGGKARAKAKTRSSRAGLQFPVGRVHRLLRKGNYSERVGAGAPVYLAAVLEYLTAEILE".
                    "LAGNAARDNKKTRIIPRHLQLAIRNDEELNKLLGRVTIAQGGVLPNIQAVLLPKKTESHHKAKGK";

my @table2 = qw(Hist1h2aa  Hist1h2ab  Hist1h2ac Hist1h2ad Hist1h2ae Hist1h2af
                Hist1h2ag  Hist1h2ah  Hist1h2ai Hist1h2ak Hist1h2an Hist1h2ao
                Hist2h2aa1 Hist2h2aa2 Hist2h2ab Hist2h2ac Hist3h2a);

## Fill all the H2A sequences with the consensus sequence, and then
## apply the differences to each of them.
$histone{$_} = $H2A_consensus foreach @table2;

for ("Hist1h2aa") {
  substr ($histone{$_},   2, 1) = "P";
  substr ($histone{$_},   3, 1) = "T";
  substr ($histone{$_},   5, 1) = "R";
  substr ($histone{$_},  13, 1) = "V";
  substr ($histone{$_},  15, 1) = "S";
  substr ($histone{$_},  35, 1) = "Q";
  substr ($histone{$_},  39, 1) = "A";
  substr ($histone{$_},  40, 1) = "Q";
  substr ($histone{$_},  42, 1) = "I";
  substr ($histone{$_},  61, 1) = "V";
  substr ($histone{$_},  78, 1) = "T";
  substr ($histone{$_}, 123)    = "KSQTK";
}
substr ($histone{Hist1h2af}, 125, 1) = "P";
substr ($histone{Hist1h2ah}, 124)    = "KAK";
for ("Hist2h2aa1", "Hist2h2aa2") {
  substr ($histone{$_}, 15, 1) = "S";
  substr ($histone{$_}, 39, 1) = "A";
  substr ($histone{$_}, 50, 1) = "M";
  substr ($histone{$_}, 98, 1) = "K";
}
for ("Hist2h2ab") {
  substr ($histone{$_}, 15, 1) = "S";
  substr ($histone{$_}, 39, 1) = "A";
  substr ($histone{$_}, 50, 1) = "M";
  substr ($histone{$_}, 98, 1) = "K";
  substr ($histone{$_}, 123)   = "KPGKNK";
}
for ("Hist2h2ac") {
  substr ($histone{$_}, 15, 1) = "S";
  substr ($histone{$_}, 39, 1) = "A";
  substr ($histone{$_}, 50, 1) = "M";
  substr ($histone{$_}, 98, 1) = "K";
  substr ($histone{$_}, 123)   = "KAKSK";
}
substr ($histone{Hist3h2a}, 15, 1) = "S";

## H2B consensus sequence from Table 3: Histone H2B protein variants (page 495)
my $H2B_consensus = "PEPAKSAPAPKKGSKKAVTKAQKKDGKKRKRSRKESYSVYVYKVLKQVHPDTGISSKAMGIM".
                    "NSFVNDIFERIASEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVTKYTSSK";

my @table3 = qw(Hist1h2ba Hist1h2bb Hist1h2bc Hist1h2be Hist1h2bf Hist1h2bg
                Hist1h2bh Hist1h2bj Hist1h2bk Hist1h2bm Hist1h2bn Hist1h2bp
                Hist2h2bb Hist2h2be Hist3h2ba Hist3h2bb);

# Fill all the H2B sequences with the consensus sequence, and then
# apply the differences to each of them.
$histone{$_} = $H2B_consensus foreach @table3;

for ("Hist1h2ba") {
  substr ($histone{$_},   5, 1) = "G";
  substr ($histone{$_},   7, 3) = "TIS";
  substr ($histone{$_},  13, 1) = "F";
  substr ($histone{$_},  20, 1) = "T";
  substr ($histone{$_},  24, 1) = "E";
  substr ($histone{$_},  26, 1) = "R";
  substr ($histone{$_},  31, 1) = "C";
  substr ($histone{$_},  38, 1) = "I";
  substr ($histone{$_},  40, 1) = "I";
  substr ($histone{$_},  59, 1) = "S";
  substr ($histone{$_},  66, 1) = "T";
  substr ($histone{$_},   2, 2) = "VAV"; # at the end because of the addition
}
for ("Hist1h2bb") {
  substr ($histone{$_},  3, 1) = "S";
  substr ($histone{$_}, 17, 1) = "I";
  substr ($histone{$_}, 18, 1) = "S";
}
for ("Hist1h2bc", "Hist1h2be", "Hist1h2bg") {
  substr ($histone{$_}, 74, 1) = "G";
}
for ("Hist1h2bh") {
  substr ($histone{$_}, 17, 1) = "L";
  substr ($histone{$_}, 74, 1) = "G";
}
substr ($histone{Hist1h2bk}, 123, 1) = "A";
for ("Hist1h2bm") {
  substr ($histone{$_},  3, 1) = "T";
  substr ($histone{$_}, 74, 1) = "G";
}
for ("Hist1h2bp") {
  substr ($histone{$_},  3, 1) = "V";
  substr ($histone{$_},  6, 1) = "V";
}
for ("Hist2h2bb") {
  substr ($histone{$_},  1, 1) = "D";
  substr ($histone{$_}, 20, 1) = "V";
  substr ($histone{$_}, 74, 1) = "G";
}
for ("Hist2h2be") {
  substr ($histone{$_},   2, 1) = "L";
  substr ($histone{$_},  38, 1) = "I";
  substr ($histone{$_},  74, 1) = "N";
  substr ($histone{$_},  96, 1) = "S";
  substr ($histone{$_}, 123, 1) = "A";
}
for ("Hist3h2ba") {
  substr ($histone{$_},   3, 1) = "S";
  substr ($histone{$_},   4, 1) = "R";
  substr ($histone{$_},   6, 1) = "T";
  substr ($histone{$_},  17, 1) = "I";
  substr ($histone{$_},  31, 1) = "G";
  substr ($histone{$_},  38, 1) = "I";
  substr ($histone{$_},  93, 1) = "V";
}
for ("Hist3h2bb") {
  substr ($histone{$_},   1, 1) = "D";
  substr ($histone{$_},   3, 1) = "S";
  substr ($histone{$_},  31, 1) = "G";
  substr ($histone{$_},  38, 1) = "I";
  substr ($histone{$_},  93, 1) = "V";
}

## There is no protein sequence for H3. In human it was easier because there's
## a description for the variants, and each one is limited to one of the
## clusters. However, in mice this does not happen. Instead, we have:
##
##    ... rodents have much more H3.2 protein, and four of the nine histone
##    H3 genes in the mouse Hist1 cluster encode the H3.2 protein [13]. In
##    addition there are multiple genes encoding the H3.2 protein in the mouse
##    Hist2 cluster [14].
##
## There is no reference on which of the listed genes encode which protein, only
## how many. The references point to papers before the nomenclature so care
## must be taken when interpreting results from the two together.
##
## Figuring out H3.2 genes in Hist1
##
##    [13]  Wang, Z.-F., et al. (1996). Characterization of the mouse histone
##          gene cluster on chromosome 13: 45 histone genes in three patches
##          spread over one megabase. Genome Res. 6: 688–701.
##
## This paper has a doi:10.1101/gr.6.8.688 and says the following about the
## the H3 genes in Hist1
##
##    Of the five previously characterized histone H3 genes from mouse
##    chromosome 13, two encode the histone H3.1 protein and three encode the
##    histone H3.2 protein (Taylor et al. 1986; Gruber et al. 1990; Brown
##    et al. 1996). The remaining four histone H3 genes from these YACs were
##    isolated and sequenced. Two of these genes encoded the H3.1 protein
##    (H3-D and H3-I), and two encoded the H3.2 protein (H3-B and H3-F).
##
## Using Marzluff's Table 1, we can map H3-B and H3-F to Hist1h3d and Hist1h3e
## (note that the end letters do not match at all) as H3.2, and H3-D to
## Hist1h3a as H3.1. We do not know which one was H3-I in the new nomenclature.
## Looking at its Figure 6, we can also see that H3-H, mapped to Hist1h3i
## encodes an H3.1.
##
## Following the references in [13], both Gruber et al. 1990 and Brown et al.
## 1996 were behind a paywall. In Taylor, 1986, we find that H3-291, mapped to
## Hist1h3h, encodes H3.1:
##
##    Of the two newly sequenced genes, H3.291 codes for an H3.1 protein
##
## and that H3.1-221, mapped to Hist1h3g also encodes H3.1. We can't see the
## reports to map all of the H3.2, but we already mapped all four H3.1
## (Hist1h3a, Hist1h3i, Hist1h3h, Hist1h3g), so the missing two, Hist1h3b,
## and Hist1h3c must encode H3.2 (the Hist1h3f aso encodes H3.2 from Taylors
## paper)
##
## Figuring out H3.2 genes in Hist2
##
##    [14]  Wang, Z.-F., et al. (1996). Characterization of the 55 kilobase
##          mouse histone gene cluster on chromosome 3. Genome Res. 6: 702–714.
##
## Figure 1 from [14] matches Figure 2 of Marzluff paper, and the H3 protein
## sequences appear on Figure 4. Using Marzluff's Table 1, we can map Hist2h3b,
## Hist2h3ca1, and Hist2h3ca2 to H3-616, H3-615, and H3-614. While Marzluff
## says that are "multiple genes enconding the H3.2 protein in the mouse Hist2
## cluster", the Wang paper says, about this specific cluster, "All three
## histone H3 genes encode the histone H3.2 protein". The sequence on its
## Figure 4 is equivalent to what we deduced for human H3.
##
##
## HUMAN DEDUCTION OF H3 PROTEIN SEQUENCE
##
## There was no actual H3 protein sequence on the paper, we had to deduce and
## track it down on the references.
##
## It is written on page 489:
##
##    ... the histone H3 genes encode the two variants previously described [15] ...
##
## The reference [15] is "Zweidler, A. (1984). Core histone variants of the mouse: primary
## structure and expression. In Histone Genes: Structure, Organization and Regulation (G. Stein,
## W.Stein W. F. Marzluff, Eds), pp.373-395. John Wiley and Sons, New York" (the correct page
## numbers are actually 339-371).
##
## Later, on page 494:
##
##    There are two variants encoded by the histone H3 genes, H3.1 and H3.2. These
##    differ in a single amino acid change at position 96, a cysteine in H3.1 genes
##    and a serine in H3.2 genes. The functional histone H3 gene (HIST2H3C) in the
##    HIST2 cluster on chromosome 1 encodes the H3.2 protein and all 11 histone H3
##    genes in the HIST1 cluster encode the H3.1 protein.
##
## However, there are only 10 histone H3 genes for the HIST1 cluster on Table 1. Also, it says
## nothing about the protein encoded by HIST3H3. That is mentioned on page 492:
##
##    HIST3H3, described previously as a solitary gene [19], encodes a novel histone
##    H3 protein, with the serine at position 96 and four additional changes (A24V,
##    V71M, A98S, A111V)
##
## This leads to a total of 3 unique H3 protein sequences: all the H3 on HIST1 encoding the same
## thing; a single coding H3 gene on HIST2 with just one amino acid difference; and the lone
## H3 gene on HIST3 with the five changes.

foreach (keys $table1{'coding'}{'H3'}) {
  ## sequence from Zweidler, A (1984) chapter on the histone genes book
  $histone{$_} = "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKST".
                 "ELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEACEAYLVGLFEDTNLCAIHAKR".
                 "VTIMPKDIQLARRIRGERA";
}
for ("Hist1h3b", "Hist1h3c", "Hist1h3d", "Hist1h3e", "Hist1h3f",
     "Hist2h3b", "Hist2h3ca1", "Hist2h3ca2") {
  substr ($histone{$_}, 95, 1) = "S";
}

## There is no protein sequence for H4. But it is said that it is the same as
## the human sequence on page 494:
##
##    Even more impressive are the 26 nonallelic human and mouse histone H4 genes
##    that encode exactly the same protein, ...
##
## The protein sequence for the human version was also not given. Instead,
## on the caption of Table 1 it is said
##
##    The human H4 gene described by Akasaka and coworkers [56] is HISTH4I[sic]
##
## Then, on page 489, it is written
##
##    all the histone H4 genes encode the same protein
##
## thus logically stating that all have the same sequence from reference 56 which targets "Akasaka, T.,
## et al. (1997). A recurring translocation, t(3;6)(q27;p21), in non-Hodgkin's lymphoma results in
## replacement of the 5' regulatory region of BCL6 with a novel H4 histone gene. Cancer Res. 57: 7--12."
foreach (keys $table1{'coding'}{'H4'}) {
  ## sequence from Akasaka, T et al (1997)
  $histone{$_} = "SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVL".
                 "KVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG";
}

## Fill the proteins directory with sequence files
while (my ($symbol, $seq) = each %histone) {
  ## Histones sequences are usually given without the first methionine.
  ## Put it back before saving it as a file.
  $seq = "M" . $seq;
  ## Remove newlines from the sequences because it confuses Bio::Seq
  ## We name the file based on the symbol because we don't have UIDs or
  ## accession numbers for all of them.
  Bio::SeqIO->new(
    -format => 'genbank',
    -file   => ">" . File::Spec->catdir ($proteins_dir, $symbol.".gb"),
  )->write_seq(Bio::Seq->new(
    -seq        => $seq,
    -display_id => $symbol,
  ));
}

## Create a data.csv file as if generated by out script. But this time, the
## only data we have is the gene symbol and whether they are coding or pseudo
## genes.
my $data_path = File::Spec->catdir($FindBin::Bin, "data.csv");
open (my $data, ">", $data_path) or die "Could not open $data_path for writing: $!";
my $csv = Text::CSV->new ({
  binary => 1,
  eol    => $/,
}) or die "Cannot use Text::CSV: ". Text::CSV->error_diag ();

$csv->print ($data, ['gene symbol', 'pseudo']);
$csv->print ($data, [$_, 0]) foreach keys %histone;
foreach my $type (keys $table1{'pseudo'}) {
  $csv->print ($data, [$_, 1]) foreach keys $table1{'pseudo'}{$type};
}
close ($data) or die "Couldn't close $data_path after writing: $!";

## Just for curiosity, we compared the sequences mentioned throughout the
## paper against the sequences obtained from their accession numbers.
## Does the sequences from the accession number on table 1, match what was
## actually written on the paper? The answer is no, there are 16 mistakes.
my $mistakes = 0;
foreach my $type (sort keys $table1{'coding'}) {
  my $fetcher = Bio::DB::EUtilities->new(
    -eutil   => "efetch",
    -db      => "nucleotide",
    -id      => [sort values $table1{'coding'}{$type}],
    -email   => $email,
    -retmode => "text",
    -rettype => "gb",
  );
  my $response = $fetcher->get_Response->content;
  open (my $seq_fh, "<", \$response)
    or die "Could not open sequences string for reading: $!";
  my $seqs = Bio::SeqIO->new(
    -fh     => $seq_fh,
    -format => "genbank",
  );
  while (my $seq = $seqs->next_seq()) {
    my $accession = $seq->display_id();
    my $symbol = List::Util::first {
      $table1{'coding'}{$type}{$_} eq $accession
    } keys $table1{'coding'}{$type};
    my $protein = (($seq->get_SeqFeatures("CDS"))[0]->get_tag_values("translation"))[0];
    if (! exists $histone{$symbol}) {
      say "coding gene $symbol with accession $accession has no protein sequence";
      $mistakes++;
    ## substr - because when dealing with histones, we never count the first
    ##residue and so we need to "unfix" the annotation
    } elsif (substr ($protein, 1) ne $histone{$symbol}) {
      say "different protein sequence for $symbol";
      $mistakes++;
    }
  }
}

## reverse check. Do all the proteins mentioned on this table, appear on the
## big table 1 with all genes and teir accession numbers? We already know this
## is not true. Table 1 says that HIST1H2AF is not  present in humans but
## Table 2 says it has the same sequence as the consensus. And HIST1H2AG
## appears as coding protein on Table 1 but is absent from Table 2 where its
## sequence should be.
sub find_on_table1 {
  my $type = shift;
  foreach my $gene (@_) {
    if (! grep {$gene eq $_} keys $table1{'coding'}{$type}) {
      say "gene $gene has protein sequence but no accession number";
      $mistakes++;
    }
  }
}
find_on_table1 ("H2A", @table2);
find_on_table1 ("H2B", @table3);

say "Found a total of $mistakes mistakes";
