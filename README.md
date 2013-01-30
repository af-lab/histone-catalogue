Human Histone Catalogue
=======================

A catalogue of all the human histone genes and proteins, with discussion
on the definition of histone variables, histone isoforms, histone clusters,
and what it means to its nomenclature.

This paper is generated automatically. The concept is that when the databases
are updated, this paper can be updated automatically. There is no actual
data on this repository, only code to download the actual date and build a
paper with the current information. Only the text can be found on the LaTeX
source, all the tables and figures are created during the build of the pdf.

Building instructions
---------------------


Dependencies
------------

The following needs to be installed for a successful build:

* [scons](www.scons.org) - the build system
* a TeX distribution with pdfTeX such as [TeX Live](http://www.tug.org/tex-live/)
* perl - at least version 5.010 and the following modules
  * Lingua::EN::Numbers
  * Text::CSV
  * Bio::EUtilities
  * Bio::DB::SeqIO
* [weblogo](weblogo.threeplusone.com) - generate the sequence logos
