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

use FindBin;
use File::Temp;
use File::Spec;

use Test::More;
use Test::Exception;

use HistoneSequencesDB;

## Creates a temporary directory to mimic one for HistoneSequencesDB.
## Takes the contents of csv file to be created.
sub create_seq_dir
{
  my $csv_contents = shift;

  my $dir = File::Temp->newdir();
  my $csv_path = File::Spec->catfile($dir->dirname, "data.csv");
  open (my $fh, ">", $csv_path)
    or die "Can't open '$csv_path' to write: $!";
  print {$fh} $csv_contents;
  close $fh;
  return $dir;
}

sub test_db
{
  my $db = shift;

  my @genes = @{$db->genes};

  is (scalar @genes, 16, "Read right number of genes");

  is_deeply ([sort (map {$_->symbol} @genes)],
             ['CENPA', 'H1F0', 'H1FX', 'H2AFJ', 'H2AFZ', 'H2AFZP4', 'H2BFFD',
              'HIST1H1D', 'HIST1H1T', 'HIST1H2APS4', 'HIST1H2BD',
              'HIST1H2BTTF', 'HIST1H3F', 'HIST1H4I', 'HIST1H4K', 'HIST1H4L'],
             "Got the right gene symbols");


  is_deeply ([sort {$a <=> $b} (map {$_->uid} @genes)],
             [1058, 3005, 3007, 3010, 3015, 3017, 8294, 8333, 8362,
              8368, 8968, 8971, 55766, 99999998, 99999999, 100462795],
             "Got the right gene UID");

  is_deeply ([map {$_->species} @genes], [('Homo sapiens') x 16],
             "Read the right species name");

  {
    my ($hist1h2bd) = grep {$_->symbol eq 'HIST1H2BD'} @genes;
    my %products = %{$hist1h2bd->products};
    is (scalar keys %products, 2,
        "Read gene canonical histone with multiple products correctly");
    is_deeply([sort (keys %products)], ["NM_021063", "NM_138720"],
              "Read transcripts acc from multiple products");
    is_deeply([sort (values %products)], ["NP_066407", "NP_619790"],
              "Read protein acc from multiple products");
    is ("NC_000006", $hist1h2bd->chr_acc, "Read chromosome accession");
    is ("26157343", $hist1h2bd->chr_start, "Read chromosome start coordinates");
    is ("26171849", $hist1h2bd->chr_end, "Read chromosome end coordinates");
    is ("histone cluster 1, H2bd", $hist1h2bd->description, "Read description");
    is ($hist1h2bd->histone_type, 'H2B', "Read histone type of coding multiple products");
  }

  my ($h2afzp4) = grep {$_->symbol eq 'H2AFZP4'} @genes;
  is ($h2afzp4->type, 'pseudo', "Read pseudo variant as pseudo");
  is ($h2afzp4->histone_type, 'H2A', "Read histone type of pseudo variant");

  my ($cenpa) = grep {$_->symbol eq 'CENPA'} @genes;
  is (scalar keys %{$cenpa->products}, 2,
      "Read gene CENPA variant with multiple products correctly");
  is ($cenpa->histone_type, 'H3', "Read histone type of CENPA");

  my ($hist1h2aps4) = grep {$_->symbol eq 'HIST1H2APS4'} @genes;
  is ($hist1h2aps4->type, 'pseudo', "Read pseudo gene correctly");
  is ($hist1h2aps4->histone_type, 'H2A', "Read histone type of pseudo canonical");

  my ($hist1h4ai) = grep {$_->symbol eq 'HIST1H4I'} @genes;
  is ($hist1h4ai->type, 'coding', "Read coding gene correctly");
  is ($hist1h4ai->histone_type, 'H4', "Read histone type of coding canonical");

  my ($h2afj) = grep {$_->symbol eq 'H2AFJ'} @genes;
  is ($h2afj->type, 'coding', "Read coding gene with non-coding transcript");
  is_deeply ($h2afj->products, {'NM_177925' => 'NP_808760', 'NR_027716' => ''},
    'read products with non-coding transcripts');

  my @canonical_core = $db->canonical_core;
  is_deeply ([sort map {$_->symbol} @canonical_core],
             ['HIST1H2APS4', 'HIST1H2BD', 'HIST1H2BTTF', 'HIST1H3F', 'HIST1H4I',
              'HIST1H4K', 'HIST1H4L'],
             "Check grep of canonical core histones");

  my @linkers = $db->linkers;
  is_deeply ([sort map {$_->symbol} @linkers],
             ['H1F0', 'H1FX', 'HIST1H1D', 'HIST1H1T'],
             "Check grep of linker histones");

  my @variants = $db->variants_core;
  is_deeply ([sort map {$_->symbol} @variants],
             ['CENPA', 'H2AFJ', 'H2AFZ', 'H2AFZP4', 'H2BFFD'],
             "Check grep of variant histones");

  is_deeply ([map {$_->symbol} HistoneSequencesDB::sort_histones(@canonical_core)],
             ['HIST1H2APS4', 'HIST1H2BD', 'HIST1H2BTTF', 'HIST1H3F', 'HIST1H4I',
              'HIST1H4K', 'HIST1H4L'],
             "Check sorting of array of histones");

  is_deeply ([map {$_->symbol} HistoneSequencesDB::sort_histones(@variants)],
             ['H2AFJ', 'H2AFZ', 'H2AFZP4', 'H2BFFD', 'CENPA'],
             "Check sorting of array of variant histones with CENPA");

  my ($h2bffd) = grep { $_->symbol eq 'H2BFFD' } @genes;
  is ($h2bffd->description, "H2B fake for test",
      'Confirm we got the right fake gene for testing');
  is_deeply($h2bffd->products, {'NR_027999' => '', 'NR_027998' => '',
                                'NM_177999' => 'NP_808999'},
            "Find products of gene with non-coding transcript defined first");

  my ($hist1h2bttf) = grep { $_->symbol eq 'HIST1H2BTTF' } @genes;
  is ($hist1h2bttf->description, "H2B fake for test 2",
      'Confirm we got the right other fake gene for testing');
  is_deeply($hist1h2bttf->products, {'NR_027997' => '', 'NR_027996' => '',
                                     'NM_177998' => 'NP_808998'},
            "Find products of gene with non-coding transcript defined first");

}

## A sample from actual data for testing.  It has canonical histones,
## both pseudo and coding, and with multiple products, the CENPA variant
## with its atypical gene symbol, a non histone gene, and a coding gene
## with non-coding transcripts.
## Also has the fake H2BFD and HIST1H2BTTF, which are coding genes with
## some non-coding transcripts.  They differ on H2AFJ for testing in that
## their first csv lines are the non-coding transcript.
my $dir = create_seq_dir(<<END);
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
H2BFD,"Homo sapiens",99999999,ENSG99999999999,"H2B fake for test",0,NR_027999,,12p12.3,NC_000012,14773836,14778502,"Reference GRCh38.p2 Primary Assembly"
H2BFFD,"Homo sapiens",99999999,ENSG99999999999,"H2B fake for test",0,NR_027998,,12p12.3,NC_000012,14773836,14778502,"Reference GRCh38.p2 Primary Assembly"
H2BFFD,"Homo sapiens",99999999,ENSG99999999999,"H2B fake for test",0,NM_177999,NP_808999,12p12.3,NC_000012,14773836,14778502,"Reference GRCh38.p2 Primary Assembly"
HIST1H2BTTF,"Homo sapiens",99999998,ENSG99999999998,"H2B fake for test 2",0,NR_027997,,12p12.3,NC_000012,14773836,14778502,"Reference GRCh38.p2 Primary Assembly"
HIST1H2BTTF,"Homo sapiens",99999998,ENSG99999999998,"H2B fake for test 2",0,NM_177998,NP_808998,12p12.3,NC_000012,14773836,14778502,"Reference GRCh38.p2 Primary Assembly"
HIST1H2BTTF,"Homo sapiens",99999998,ENSG99999999998,"H2B fake for test 2",0,NR_027996,,12p12.3,NC_000012,14773836,14778502,"Reference GRCh38.p2 Primary Assembly"
H2AFZP4,"Homo sapiens",100462795,,"H2A histone family, member Z pseudogene 4",1,,,11,NC_000011,70278415,70279797,"Reference GRCh38.p2 Primary Assembly"
H2AFZ,"Homo sapiens",3015,ENSG00000164032,"H2A histone family, member Z",0,NM_002106,NP_002097,4q24,NC_000004,99947587,99950855,"Reference GRCh38.p2 Primary Assembly"
H1F0,"Homo sapiens",3005,ENSG00000189060,"H1 histone family, member 0",0,NM_005318,NP_005309,22q13.1,NC_000022,37804607,37807936,"Reference GRCh38.p2 Primary Assembly"
H1FX,"Homo sapiens",8971,ENSG00000184897,"H1 histone family, member X",0,NM_006026,NP_006017,3q21.3,NC_000003,129314271,129316777,"Reference GRCh38.p2 Primary Assembly"
HIST1H1T,"Homo sapiens",3010,ENSG00000187475,"histone cluster 1, H1t",0,NM_005323,NP_005314,6p21.3,NC_000006,26106912,26108636,"Reference GRCh38.p2 Primary Assembly"
HIST1H1D,"Homo sapiens",3007,ENSG00000124575,"histone cluster 1, H1d",0,NM_005320,NP_005311,6p21.3,NC_000006,26233712,26235488,"Reference GRCh38.p2 Primary Assembly"
END

my $db = HistoneSequencesDB->new($dir->dirname);
test_db($db);

my $db_store = File::Temp->new();
$db->write_db($db_store->filename);
my $read_db = HistoneSequencesDB::read_db($db_store->filename);
isa_ok($read_db, 'HistoneSequencesDB', "Read back HistoneSequencesDB from store");

subtest "Test everything again with the read HistoneSequencesDB"
  => sub { test_db($read_db); };

## Tests that require actual sequence files
{
  my $db_dir = File::Spec->catfile($FindBin::Bin, "test-data", "test-sequences");
  my $db2 = HistoneSequencesDB->new($db_dir);

  ## we pick two pseudo genes and a gene with two transcripts, so we can
  ## check the number of iterations.
  my @genes = grep { $_->symbol =~ /HIST1H2APS4|HIST1H2BD|HIST1H3F|HIST1H4I|H2AFZP4/ } @{$db2->genes};

  my @lengths;
  my $c = 0;

  my $it = $db2->protein_iterator(HistoneSequencesDB::sort_histones (@genes));
  while (my $p = $it->())
    {
      $c++;
      push @lengths, $p->length;
    }
  is ($c, 4, "iterate over proteins right number of times");
  is_deeply (\@lengths, [125, 125, 135, 102],
             "iterate over the right proteins in right order");
  ok (! defined $it->(), "exhausted iterator stays undef");
}


## One of the problems of using a csv file with one product per line when
## what we want is genes, is that a gene gets split over lines and then we
## have to join them.  We had a weird bug that was triggered only when one
## such gene was the first row in the csv file and the second row could not
## make a valid gene (H2AFJ has a non-protein coding transcript).  This test
## case triggers it.
sub test_split_gene_first_2_rows
{
  my $db_dir = shift;
  my $db2 = HistoneSequencesDB->new($db_dir->dirname);
  isa_ok($db2, 'HistoneSequencesDB', "Created HistoneSequencesDB");

  my @genes = @{$db2->genes};
  is (scalar @genes, 2, "Read right number of genes");

  is_deeply ([sort (map {$_->symbol} @genes)], ['H2AFJ', 'HIST1H4K',],
             "Got the right gene symbols");

  is_deeply ([sort {$a <=> $b} (map {$_->uid} @genes)], [8362, 55766],
             "Got the right gene UID");

  my ($h2afj) = grep {$_->symbol eq 'H2AFJ'} @genes;
  is ($h2afj->type, 'coding', "Read coding gene with non-coding transcript");
  is_deeply ($h2afj->products, {'NM_177925' => 'NP_808760', 'NR_027716' => ''},
    'read products with non-coding transcripts');
}

#my $dir2 = ;
subtest "Special case with 'bad' gene split on first two rows"
  => sub { test_split_gene_first_2_rows(create_seq_dir(<<END)); };
"gene symbol",species,"gene UID","EnsEMBL ID","gene name",pseudo,"transcript accession","protein accession",locus,"chromosome accession","chromosome start coordinates","chromosome stop coordinates",assembly
H2AFJ,"Homo sapiens",55766,ENSG00000246705,"H2A histone family, member J",0,NM_177925,NP_808760,12p12.3,NC_000012,14773836,14778502,"Reference GRCh38.p2 Primary Assembly"
H2AFJ,"Homo sapiens",55766,ENSG00000246705,"H2A histone family, member J",0,NR_027716,,12p12.3,NC_000012,14773836,14778502,"Reference GRCh38.p2 Primary Assembly"
HIST1H4K,"Homo sapiens",8362,ENSG00000273542,"histone cluster 1, H4k",0,NM_003541,NP_003532,6p22.1,NC_000006,27830674,27832027,"Reference GRCh38.p2 Primary Assembly"
END

done_testing;
