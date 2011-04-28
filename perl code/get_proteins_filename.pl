#!/usr/bin/perl -w
use 5.010;                      # Use Perl 5.10
use strict;                     # Enforce some good programming rules
use Bio::SeqIO;                 # Handler for SeqIO Formats
use File::Spec;                 # Portably perform operations on file names

my $dirname = $ARGV[0];

opendir(DIR, $dirname) or die "Can't opendir $dirname: $!\n";
while (defined(my $filename = readdir(DIR))) {

  if ($filename =~ /protein/) {
    say ">$filename";
    my $in    = Bio::SeqIO->new(-file   => "$dirname/$filename", -format=>"fasta");
    while (my $seq = $in->next_seq) {

#      $seq->display_id($filename);
#      $seq->desc("");
#      my $out = Bio::SeqIO->newFh(-format => "fasta");
#      $out    ->write_seq($seq);
      say $seq->seq;
    }
  }
}
