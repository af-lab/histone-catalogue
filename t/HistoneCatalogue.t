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

use Test::Output;
use Test::More;
use File::Temp;

use Test::Exception;

use Bio::Seq;
use Bio::SimpleAlign;
use Bio::LocatableSeq;

use HistoneCatalogue;
use HistoneSequencesDB;

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

## Creates a temporary directory to mimic one for HistoneSequencesDB.
## Takes the contents of csv file to be created.
sub create_test_db
{
  ## A sample from actual data for testing.  It has canonical histones,
  ## both pseudo and coding, and with multiple products, the CENPA variant
  ## with its atypical gene symbol, a non histone gene, and a coding gene
  ## with non-coding transcripts.
  my $csv_contents = <<END;
"gene symbol",species,"gene UID","EnsEMBL ID","gene name",pseudo,"transcript accession","protein accession",locus,"chromosome accession","chromosome start coordinates","chromosome stop coordinates",assembly
HIST1H4K,"Homo sapiens",8362,ENSG00000273542,"histone cluster 1, H4k",0,NM_003541,NP_003532,6p22.1,NC_000006,27830674,27832027,"Reference GRCh38.p2 Primary Assembly"
HIST1H2APS4,"Homo sapiens",8333,,"histone cluster 1, H2a, pseudogene 4",1,,,6p21.3,NC_000006,26271693,26273040,"Reference GRCh38.p2 Primary Assembly"
HIST1H4L,"Homo sapiens",8368,ENSG00000275126,"histone cluster 1, H4l",0,NM_003546,NP_003537,6p22.1,NC_000006,27872648,27874011,"Reference GRCh38.p2 Primary Assembly"
HIST1H3F,"Homo sapiens",8968,ENSG00000277775,"histone cluster 1, H3f",0,NM_021018,NP_066298,6p22.2,NC_000006,26249642,26251107,"Reference GRCh38.p2 Primary Assembly"
HIST1H2BD,"Homo sapiens",3017,ENSG00000158373,"histone cluster 1, H2bd",0,NM_138720,NP_619790,6p21.3,NC_000006,26157343,26171849,"Reference GRCh38.p2 Primary Assembly"
HIST1H2BD,"Homo sapiens",3017,ENSG00000158373,"histone cluster 1, H2bd",0,NM_021063,NP_066407,6p21.3,NC_000006,26157343,26171849,"Reference GRCh38.p2 Primary Assembly"
CENPA,"Homo sapiens",1058,ENSG00000115163,"centromere protein A",0,NM_001042426,NP_001035891,2p23.3,NC_000002,26785514,26795089,"Reference GRCh38.p2 Primary Assembly"
CENPA,"Homo sapiens",1058,ENSG00000115163,"centromere protein A",0,NM_001809,NP_001800,2p23.3,NC_000002,26785514,26795089,"Reference GRCh38.p2 Primary Assembly"
RNASE9,"Homo sapiens",390443,ENSG00000188655,"ribonuclease, RNase A family, 9 (non-active)",0,NM_001001673,NP_001001673,14q11.2,NC_000014,20555593,20561431,"Reference GRCh38.p2 Primary Assembly"
RNASE9,"Homo sapiens",390443,ENSG00000188655,"ribonuclease, RNase A family, 9 (non-active)",0,NM_001110357,NP_001103827,14q11.2,NC_000014,20555593,20561431,"Reference GRCh38.p2 Primary Assembly"
RNASE9,"Homo sapiens",390443,ENSG00000188655,"ribonuclease, RNase A family, 9 (non-active)",0,NM_001110361,NP_001103831,14q11.2,NC_000014,20555593,20561431,"Reference GRCh38.p2 Primary Assembly"
HIST1H4I,"Homo sapiens",8294,09109,"histone cluster 1, H4i",0,NM_003495,NP_003486,6p21.33,NC_000006,27138809,27140178,"Reference GRCh38.p2 Primary Assembly"
H2AFJ,"Homo sapiens",55766,ENSG00000246705,"H2A histone family, member J",0,NM_177925,NP_808760,12p12.3,NC_000012,14773836,14778502,"Reference GRCh38.p2 Primary Assembly"
H2AFJ,"Homo sapiens",55766,ENSG00000246705,"H2A histone family, member J",0,NR_027716,,12p12.3,NC_000012,14773836,14778502,"Reference GRCh38.p2 Primary Assembly"
H2AFZP4,"Homo sapiens",100462795,,"H2A histone family, member Z pseudogene 4",1,,,11,NC_000011,70278415,70279797,"Reference GRCh38.p2 Primary Assembly"
H2AFZ,"Homo sapiens",3015,ENSG00000164032,"H2A histone family, member Z",0,NM_002106,NP_002097,4q24,NC_000004,99947587,99950855,"Reference GRCh38.p2 Primary Assembly"
H1F0,"Homo sapiens",3005,ENSG00000189060,"H1 histone family, member 0",0,NM_005318,NP_005309,22q13.1,NC_000022,37804607,37807936,"Reference GRCh38.p2 Primary Assembly"
H1FX,"Homo sapiens",8971,ENSG00000184897,"H1 histone family, member X",0,NM_006026,NP_006017,3q21.3,NC_000003,129314271,129316777,"Reference GRCh38.p2 Primary Assembly"
HIST1H1T,"Homo sapiens",3010,ENSG00000187475,"histone cluster 1, H1t",0,NM_005323,NP_005314,6p21.3,NC_000006,26106912,26108636,"Reference GRCh38.p2 Primary Assembly"
HIST1H1D,"Homo sapiens",3007,ENSG00000124575,"histone cluster 1, H1d",0,NM_005320,NP_005311,6p21.3,NC_000006,26233712,26235488,"Reference GRCh38.p2 Primary Assembly"
END

  my $dir = File::Temp->newdir();
  my $csv_path = File::Spec->catfile($dir->dirname, "data.csv");
  open (my $fh, ">", $csv_path)
    or die "Can't open '$csv_path' to write: $!";
  print {$fh} $csv_contents;
  close $fh;
  my $db = HistoneSequencesDB->new($dir->dirname);
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

is (HistoneCatalogue::mk_latex_list_name_isoforms ("HIST2H2AF", "HIST1H2AB", "HIST1H2AC", "HIST1H2AP", "HIST4H2A"),
    "HIST1H2A --B, --C, --P; HIST2H2AF; HIST4H2A",
    "list of isoforms gene symbols");
is (HistoneCatalogue::mk_latex_list_name_isoforms ("HIST2H2AA3", "HIST1H2AB", "HIST1H2AC", "HIST1H2AP", "HIST4H2A"),
    "HIST1H2A --B, --C, --P; HIST2H2AA3; HIST4H2A",
    "list of isoforms gene symbols with a weird 'AA3' descriptor");
is (HistoneCatalogue::mk_latex_list_name_isoforms ("HIST1H4A", "HIST1H4B", "HIST1H4P", "HIST4H4"),
    "HIST1H4 --A, --B, --P; HIST4H4",
    "list of isoforms gene symbols with a symbol without descriptor");
is (HistoneCatalogue::mk_latex_list_name_isoforms ("HIST2H2AB"),
    "HIST2H2AB",
    "list of gene symbols with a single element");
is (HistoneCatalogue::mk_latex_list_name_isoforms ("HIST4H4"),
    "HIST4H4",
    "list of gene symbols with a single element and without descriptor");

my $db = create_test_db ();

stdout_is (sub {HistoneCatalogue::say_histone_catalogue($db->canonical_core) },
    <<'END', 'Test canonical core latex table');
\begin{ctabular}{l l l l l}
  \toprule
  Type & Gene name & Gene UID & Transcript accession & Protein accession \\
  \midrule
  H2A & HIST1H2APS4 & 8333 & pseudogene & pseudogene\\
  H2B & HIST1H2BD & 3017 & NM\_021063 & NP\_066407\\
  H2B & HIST1H2BD & 3017 & NM\_138720 & NP\_619790\\
  H3 & HIST1H3F & 8968 & NM\_021018 & NP\_066298\\
  H4 & HIST1H4I & 8294 & NM\_003495 & NP\_003486\\
  H4 & HIST1H4K & 8362 & NM\_003541 & NP\_003532\\
  H4 & HIST1H4L & 8368 & NM\_003546 & NP\_003537\\
  \bottomrule
\end{ctabular}
END

stdout_is (sub {HistoneCatalogue::say_histone_catalogue($db->variants_core) },
    <<'END', 'Test variant histones with non-coding transcripts latex table');
\begin{ctabular}{l l l l l}
  \toprule
  Type & Gene name & Gene UID & Transcript accession & Protein accession \\
  \midrule
  H2A & H2AFJ & 55766 & NM\_177925 & NP\_808760\\
  H2A & H2AFJ & 55766 & NR\_027716 & non-coding\\
  H2A & H2AFZ & 3015 & NM\_002106 & NP\_002097\\
  H2A & H2AFZP4 & 100462795 & pseudogene & pseudogene\\
  H3 & CENPA & 1058 & NM\_001042426 & NP\_001035891\\
  H3 & CENPA & 1058 & NM\_001809 & NP\_001800\\
  \bottomrule
\end{ctabular}
END

stdout_is (sub { HistoneCatalogue::say_histone_counts($db) },
    <<'END', 'Test print of histone counts');
%% Total number of canonical core histone genes in the genome
\newcommand{\TotalCoreGenes}{\ScriptValue{6}}
%% Total number of canonical core histone protein coding genes in the genome
\newcommand{\TotalCoreCodingGenes}{\ScriptValue{5}}
%% Total number of canonical core histone protein pseudogenes in the genome
\newcommand{\TotalCorePseudoGenes}{\ScriptValue{1}}
%% Total number of histone H2A genes
\newcommand{\HTwoATotalGenes}{\ScriptValue{1}}
%% Number of histone H2A coding genes
\newcommand{\HTwoACodingGenes}{\ScriptValue{0}}
%% Number of histone H2A pseudogenes
\newcommand{\HTwoAPseudoGenes}{\ScriptValue{1}}
%% Total number of histone H2B genes
\newcommand{\HTwoBTotalGenes}{\ScriptValue{1}}
%% Number of histone H2B coding genes
\newcommand{\HTwoBCodingGenes}{\ScriptValue{1}}
%% Number of histone H2B pseudogenes
\newcommand{\HTwoBPseudoGenes}{\ScriptValue{0}}
%% Total number of histone H3 genes
\newcommand{\HThreeTotalGenes}{\ScriptValue{1}}
%% Number of histone H3 coding genes
\newcommand{\HThreeCodingGenes}{\ScriptValue{1}}
%% Number of histone H3 pseudogenes
\newcommand{\HThreePseudoGenes}{\ScriptValue{0}}
%% Total number of histone H4 genes
\newcommand{\HFourTotalGenes}{\ScriptValue{3}}
%% Number of histone H4 coding genes
\newcommand{\HFourCodingGenes}{\ScriptValue{3}}
%% Number of histone H4 pseudogenes
\newcommand{\HFourPseudoGenes}{\ScriptValue{0}}
%% Number of core histone coding genes in the histone cluster 1
\newcommand{\CoreCodingGenesInHISTOne}{\ScriptValue{5}}
%% Number of core histone pseudogenes in the histone cluster 1
\newcommand{\CorePseudoGenesInHISTOne}{\ScriptValue{1}}
%% Total number of core histone genes in the histone cluster 1
\newcommand{\TotalCoreGenesInHISTOne}{\ScriptValue{6}}
%% Number of H2A coding genes in the histone cluster 1
\newcommand{\HTwoACodingInHISTOne}{\ScriptValue{0}}
%% Number of H2A pseudogenes in the histone cluster 1
\newcommand{\HTwoAPseudoInHISTOne}{\ScriptValue{1}}
%% Total Number of H2A genes in the histone cluster 1
\newcommand{\HTwoATotalInHISTOne}{\ScriptValue{1}}
%% Number of H2B coding genes in the histone cluster 1
\newcommand{\HTwoBCodingInHISTOne}{\ScriptValue{1}}
%% Number of H2B pseudogenes in the histone cluster 1
\newcommand{\HTwoBPseudoInHISTOne}{\ScriptValue{0}}
%% Total Number of H2B genes in the histone cluster 1
\newcommand{\HTwoBTotalInHISTOne}{\ScriptValue{1}}
%% Number of H3 coding genes in the histone cluster 1
\newcommand{\HThreeCodingInHISTOne}{\ScriptValue{1}}
%% Number of H3 pseudogenes in the histone cluster 1
\newcommand{\HThreePseudoInHISTOne}{\ScriptValue{0}}
%% Total Number of H3 genes in the histone cluster 1
\newcommand{\HThreeTotalInHISTOne}{\ScriptValue{1}}
%% Number of H4 coding genes in the histone cluster 1
\newcommand{\HFourCodingInHISTOne}{\ScriptValue{3}}
%% Number of H4 pseudogenes in the histone cluster 1
\newcommand{\HFourPseudoInHISTOne}{\ScriptValue{0}}
%% Total Number of H4 genes in the histone cluster 1
\newcommand{\HFourTotalInHISTOne}{\ScriptValue{3}}
%% Total number of core histone variants genes
\newcommand{\TotalCoreVariantGenes}{\ScriptValue{4}}
END

sub compare_variant_description
{
  my $o_seq = shift;
  my $v_seq = shift;
  my $expected  = shift;
  my $test_name = shift;

  my $o = Bio::Seq->new(-seq => $o_seq);
  my $v = Bio::Seq->new(-seq => $v_seq);
  is(HistoneCatalogue::describe_protein_variant($o, $v), $expected, $test_name);
}

compare_variant_description(
  "TESHHKAK",
  "TESHHKTK",
  "A7T",
  "describe single aa substitution");

compare_variant_description(
  "TESHHKAK",
  "TESHEATK",
  "H5_A7EAT",
  "describe chunk of aa substitution");

compare_variant_description(
  "TESHHKAK",
  "TASHEATK",
  "E2A H5_A7EAT",
  "describe multiple aa substitution");

compare_variant_description(
  "MKMGHQQQCC",
  "MKMGHQQ-CC",
  "Q8del",
  "describe single aa deletion");

compare_variant_description(
  "MK---MGHQQQCC",
  "MKQSKMGHQQQCC",
  "K2_M3insQSK",
  "describe insertion of several aa");

## p.Met1ext-5
## a variant in the 5' UTR activates a new upstream translation initiation
## site starting with amino acid Met-5 (Methionine -5)
compare_variant_description(
  "-----MTAS",
  "MGTASMTAS",
  "M1ext-5",
  "describe extension on N terminus because of new upstream start codon");

## p.Met1Valext-12
## amino acid Met1 is changed to Val activating an upstream translation
## initiation site at position -12 (Methionine -12)
compare_variant_description(
  "------------MTAS",
  "MGTASHRDNKKTVTAS",
  "M1Vext-12",
  "describe extension on N terminus because of change on start codon");

## p.*110Glnext*17 (alternatively p.Ter110GlnextTer17 or p.*110Qext*17)
## describes a variant in the stop codon (Ter/*) at position 110, changing it
## to a codon for Glutamine (Gln, Q) and adding a tail of new amino acids to
## the protein's C-terminus ending at a new stop codon (Ter17/*17)
compare_variant_description(
  "GTASHRDN*-----",
  "GTASHRDNQAEIL*",
  "*9Qext*5",
  "describe extension on C terminus because of change on stop codon");

## This case is tricky because it's a deletion of the C-terminus
## but the original sequence also has a deletion inside (note that
## sequence come from a multi sequence alignment so this is normal)
compare_variant_description(
  "TESHHKAKG-K",
  "TESHHKTK---",
  "A7T G9_K10del",
  "describe protein variant with deletions at the end in both sequences");

compare_variant_description(
  "GTASHRDN*------",
  "GTASHRDNQAE-IL*",
  "*9Qext*5",
  "describe extension on C terminus with gap on the variant");

compare_variant_description(
  "SGR--GKGGKGLGKGGAKR",
  "SMRLYGK-GK--GKLLAKR",
  "G2M R3_G4insLY G6del G9_L10del G13_G14LL",
  "describe several types of protein changes");

compare_variant_description(
  "SGRGKG",
  "SGRGKG",
  "",
  "describe nothing when sequences are equal");

## we check this because our sequences may not be complete
compare_variant_description(
  "MGRGKG",
  "SGRGKG",
  "M1S",
  "do not confuse change at the start with extension of N terminus");

compare_variant_description(
  "TESHH-K--A-K",
  "TESHE-A--T-K",
  "H5_A7EAT",
  "do not get distracted by equal gaps on the sequences");

throws_ok
  { HistoneCatalogue::describe_protein_variant(Bio::Seq->new(-seq => 'SMRLY'),
                                               Bio::Seq->new(-seq => 'SMRLY-'))}
  qr/length of common and variant sequences is different/,
  'die with sequences of different length';

{
  my $aln = Bio::SimpleAlign->new();
  $aln->add_seq(Bio::LocatableSeq->new(-display_id => 'f1', -seq => 'ART-IRG-RA'));
  $aln->add_seq(Bio::LocatableSeq->new(-display_id => 'f2', -seq => 'ARTKIRGERA'));
  $aln->add_seq(Bio::LocatableSeq->new(-display_id => 'f3', -seq => 'QRTKIRGARV'));
  $aln->add_seq(Bio::LocatableSeq->new(-display_id => 'f4', -seq => 'ARTKIRGERA'));
  is(HistoneCatalogue::most_common_seq_in_alignment($aln)->seq,
     'ARTKIRGERA',
     'find most common sequence in alignment');
}
{
  my $aln = Bio::SimpleAlign->new();
  $aln->add_seq(Bio::LocatableSeq->new(-display_id => 'f1', -seq => 'ART-IRG-RA'));
  $aln->add_seq(Bio::LocatableSeq->new(-display_id => 'f2', -seq => 'QRTKIRGARV'));
  $aln->add_seq(Bio::LocatableSeq->new(-display_id => 'f3', -seq => 'ARTKIRGERA'));
  $aln->add_seq(Bio::LocatableSeq->new(-display_id => 'f4', -seq => 'QRTKIRGARV'));
  $aln->add_seq(Bio::LocatableSeq->new(-display_id => 'f5', -seq => 'ARTKIRGERA'));
  is(HistoneCatalogue::most_common_seq_in_alignment($aln)->seq,
     'ARTKIRGERA',
     'return "lowest" sequence for multiple most common sequences');
}
{
  my $aln = Bio::SimpleAlign->new();
  throws_ok { HistoneCatalogue::most_common_seq_in_alignment($aln) }
    qr/Unable to find most common sequence \(maybe empty alignment\)/,
    'die properly with empty alignment';
}


done_testing;
