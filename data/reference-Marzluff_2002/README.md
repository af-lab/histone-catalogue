Marzluff 2002 as reference for histone catalogue
================================================

The sequences found here were built automatically from the data in our
reference, the published "Marzluff, W.F., Gongidi, P., Woods, K.R., Jin, J.,
Maltais, l.J. (2002) The human and mouse replication-dependent histone genes.
Genomics (80) 5:487--498 doi:10.1006/geno.2002.6850"

There were some issues while doing this since the text, tables, and references
for online data (accession numbers on Table 1), were not always in agreement.
In addition, the sequences for the H3 and H4 proteins were not given, they
were deduced from references describing sequences and notes scattered throught
the text.

The companion script `build_marzluff_data_2002.pl` has all the data deduced
from the paper and is able to build the `data.csv` file as well as the
individual sequence files for each protein. Simply run it without any extra
arguments.

Aside building the reference data, it also performs check of the data against
the online database. Table 1 of the paper refers accession number for all the
genes. We compared those sequences against the protein sequences used in the
text for H2A (Table 2), H2B (Table 3), H3 (reference 15 and text), and H4
(reference 56). Comments on the source code have more detail. This comparison
found the following discrepancies within the paper and its data:

* different protein sequence for HIST2H2AB
* different protein sequence for HIST2H2AC
* different protein sequence for HIST3H2A
* different protein sequence for HIST1H2AA
* different protein sequence for HIST1H2AB
* different protein sequence for HIST1H2AD
* different protein sequence for HIST1H2AE
* coding gene HIST1H2AG with accession AY131987 has no protein sequence
* different protein sequence for HIST1H2AH
* different protein sequence for HIST1H2AJ
* different protein sequence for HIST1H2BA
* different protein sequence for HIST1H2BB
* different protein sequence for HIST1H2BJ
* different protein sequence for HIST3H3
* different protein sequence for HIST1H4G
* gene HIST1H2AF has protein sequence but no accession number

Found a total of 16 mistakes.

This check is interesting, but is not part of the publication. For use in the
comparison, it was used the data from the tables and text of the paper. The
reason is that such table of differences serves two purposes:

1. offer a view that online databases are not static, that things do change
with time, and it's the reference paper that should update itself. Which one
is used is irrelevant for this.
2. for anyone that used the previous publication as basis for work, a quick
list of changes. Anyone that has done this is more likely to have used the
data given in the text than the online sequences, and even if the online
sequences were used, the differences would have already been noticed. For this
purpose a comparison against the text is more helpful.

