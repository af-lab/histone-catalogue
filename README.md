Human Histone Catalogue
=======================

A catalogue of all the human histone genes and proteins, with discussion
on the definition of histone variables, histone isoforms, histone clusters,
and what it means to its nomenclature.

This paper is generated automatically. The concept is that when the databases
are updated, this paper updates itself automatically. There is no actual data on
this repository (data used for the paper version is on a separate branch. See
below for details), only code that automates the downloading and analysis of
current data, the generation of figures and tables, and building of a PDF.

Directory structure
-------------------

Note that most of the directories *should* be empty before the build. Even some
of the directories that start with some files will fill up after running scons.


* figs/ - figures
* results/ - data after processing (all generated automatically so should be empty)
  * results/sequences/ - automatically downloaded sequences
* scripts/ - directory for scripts that will analyse the data
* tables/ - LaTeX tables which will be generated automatically
* histone_catalog.tex - LaTeX source for the paper
* README.md - this file, in Markdown syntax
* references.bib - BibTeX database for the references in the paper
* SConstruct - build instructions for SCons

Building instructions
---------------------

Running `scons` on the root directory with a valid e-mail address should take
care of everything:

    scons --email example@domain.com

NCBI requires a valid e-mail address to download sequences in this automated
fashion. This is to prevent overloading of their servers. See documentation of
`bp_genbank_ref_extractor` script for details. If you see the following
message, then you already know what's wrong.

    --------------------- WARNING ---------------------
    MSG: The -email parameter is now required, per NCBI E-utilities policy
    ---------------------------------------------------

It is possible to pass commands to SCons, thus running only certain parts of the
build. Not passing a command is actually the same as running them in the
following order

    scons data --email example@domain.com
    scons analyse
    scons publish

Dependencies
------------

The following needs to be installed for a successful build:

* [SCons](www.scons.org) - the build system
* a TeX distribution with pdfTeX such as [TeX Live](http://www.tug.org/tex-live/)
* perl - at least version 5.010 and the following modules
  * Lingua::EN::Numbers
  * Text::CSV
  * Bio::EUtilities
  * Bio::DB::SeqIO
  * Bio::Tools::Run::Alignment::TCoffee
* [weblogo](weblogo.threeplusone.com) - generate the sequence logos

The plan is to add a tag for the source version that got accepted for
publication, and a branch with the data used for it. A checkout of that branch
should only require pdfTeX and BibTeX.
