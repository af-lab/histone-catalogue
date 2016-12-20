#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path
import subprocess

env = Environment()
env.Append(PERL5LIB=['lib-perl5'])
env.Tool('perl5')

## Add a NoShellCommand builder to be used like Command()
##
## This has the advantage that there's no shell involved, saving us
## from having to escape quotes, spaces, wildcards, and whatsnot.
##
## See the whole mailing list thread at
##  https://pairlist4.pair.net/pipermail/scons-users/2015-October/004150.html
## with the solution found much later with
##  https://pairlist4.pair.net/pipermail/scons-users/2015-November/004193.html

def no_shell_command(target, source, env):
  return subprocess.call(env['action'])
def no_shell_command_strfunc(target, source, env):
  args = env['action']
  return "$ %s " % (args[0]) + " ".join(["'%s'" % (arg) for arg in args[1:]])
no_shell_command_action = Action(no_shell_command, strfunction=no_shell_command_strfunc)
env.Append(BUILDERS={'NoShellCommand' : Builder(action=no_shell_command_action)})


## FIXME very temporary while we move all of MyLib to our modules
env.Append(PERL5LIB=['scripts'])


env.Help("""
By default, SCons will download new data from the Entrez Gene database,
analyse it, and prepare a ready to publish PDF of the paper. Each of these
is a different target and SCons can be used to process only part of them.
The different target names are:

    data        - download new data
    analysis    - analyse data
    manuscript  - prepare PDF of the published paper

Not specifying any target, calling `scons` on its own, is equivalent to
selecting all of them. That means, if all dependencies are installed, the
following command will take care of everything:

    scons --email=example@domain.com

Note than in order to prevent overload of the NCBI servers, a valid email is
required to connect to the Entrez Gene database using E-utilities. Other
limitations may be in place. See section *"Usage guidelines and requirements"*,
on [A General Introduction to the E-utilities](http://www.ncbi.nlm.nih.gov/books/NBK25497/).

> In order not to overload the E-utility servers, NCBI recommends that users
> post no more than three URL requests per second and limit large jobs to either
> weekends or between 9:00 PM and 5:00 AM Eastern time during weekdays. Failure
> to comply with this policy may result in an IP address being blocked from
> accessing NCBI.
>
> ...
>
> The value of email will be used only to contact developers if NCBI observes
> requests that violate our policies, and we will attempt such contact prior to
> blocking access.

Each target may require specific tools to be installed in your system, or
certain options to be specified.

OPTIONS

  --email=ADDRESS
      Set email to be used when connecting to the NCBI servers. This can
      be anything that conforms to RCF822. The following are valid:

          scons --email="Your Name <your.name@domain.here>"
          scons --email="<your.name@domain.here>"

  --organism=NAME
      Organism species name to use when searching RefSeq for histone
      sequences.  Defaults to Homo Sapiens.


  --verbose
      LaTeX and BibTeX compilers are silenced by default using the
      batchmode and terse options. Set this option to revert it.

""")

AddOption(
  "--email",
  dest    = "email",
  action  = "store",
  type    = "string",
  help    = "E-mail provided to NCBI when connecting to Entrez."
)
AddOption(
  "--verbose",
  dest    = "verbose",
  action  = "store_true",
  default = False,
  help    = "Print LaTeX and BibTeX output."
)
AddOption(
  "--organism",
  dest    = "organism",
  action  = "store",
  type    = "string",
  default = "Homo sapiens",
  help    = "Organism to search for histones."
)


if not env.GetOption("verbose"):
  env.AppendUnique(PDFLATEXFLAGS  = "-interaction=batchmode")
  env.AppendUnique(PDFTEXFLAGS    = "-interaction=batchmode")
  env.AppendUnique(TEXFLAGS       = "-interaction=batchmode")
  env.AppendUnique(LATEXFLAGS     = "-interaction=batchmode")
  env.AppendUnique(BIBTEXFLAGS    = "--terse")  # some ports of BibTeX may use --quiet instead

env.Help("""
TARGETS

  data
    Get raw data required for the analysis from the Entrez system.  This
    data is then used for analysis.

  update
    Remove previously downloaded raw data and download fresh one.

  analysis
    Run all the scripts to analyse the data such as: sequence alignments,
    search for anomalies on the sequence annotations, LaTeX tables listing
    all genes and proteins, sequence differences between isoforms. The
    analysis results are required for the publication.

  manuscript
    Build PDF for publication.

""")

scripts_dir   = "scripts"
results_dir   = "results"
figures_dir   = "figs"
seq_dir       = os.path.join (results_dir, "sequences")
reference_dir = os.path.join ("data", "reference-Marzluff_2002")

def path4lib(name):
  return os.path.join("lib-perl5", name)
def path4script (name):
  return os.path.join (scripts_dir, name)
def path4result (name):
  return os.path.join (results_dir, name)
def path4figure (name):
  return os.path.join (figures_dir, name)
def path4seq (name):
  return os.path.join (seq_dir, name)


## TARGET data
##
## SCons does not like it when the target is a directory and will always
## consider it up to date, even if the source changes.  Because of this,
## we set the data.csv and data.asn1 files as target

def create_extract_sequences_args():
  ## Gene names to use on the entrez query.  Note that:
  ##    "Right side truncation with wild card does work for gene symbol"
  ##                    ~ NCBI helpdesk via email (September 2011)
  gene_names = [
    "H1*[gene name]",
    "H2A*[gene name]",
    "H2B*[gene name]",
    "H3*[gene name]",
    "H4*[gene name]",
    "HIST1*[gene name]",
    "HIST2*[gene name]",
    "HIST3*[gene name]",
    "HIST4*[gene name]",
    "CENPA[gene name]"
  ]
  entrez_query = '"%s"[organism] AND (%s)' % (env.GetOption("organism"),
                                              " OR ".join (gene_names))

  bp_genbank_ref_extractor_call = [
    "bp_genbank_ref_extractor",
    "--assembly",     "Reference GRC",
    "--genes",        "uid",
    "--pseudo",
    "--non-coding",
    "--upstream",     "500",
    "--downstream",   "500",
    "--transcripts",  "accession",
    "--proteins",     "accession",
    "--limit",        "300",
    "--format",       "genbank",
    "--save",         seq_dir,
    "--save-data",    "csv",
    ## Should we check the email is valid?  It is important that the
    ## email is correct since this script allows one to abuse the NCBI
    ## servers who may block access. With an email address they will
    ## contact the user first.
    "--email",        env.GetOption("email"),
    entrez_query
  ]
  return bp_genbank_ref_extractor_call

## Ideally we would set target to the sequences directory and it would be
## automatically removed when it gets rebuild.  However, a Dir as targets
## means it's always out of date so it would force a rebuild every time.
## That is why we set the csv and log files as targets.
raw_data = env.NoShellCommand(
  source = None,
  target = [path4seq("data.csv"), path4seq("extractor.log")],
  action = create_extract_sequences_args())

## AddPreAction() is required so that the directory is removed when rebuilding.
## Clean() is required so that it's removed when calling "scons -c".
env.AddPreAction(raw_data, 'rm -r ' + seq_dir)
env.Clean(raw_data, seq_dir)

db_store = path4result("histones_db.store")
data_store = env.PerlSub(
  target = db_store,
  source = path4lib("HistoneSequencesDB.pm"),
  action = 'HistoneSequencesDB->new("%s")->write_db("%s")' %(seq_dir, db_store)
)
env.Depends(data_store, raw_data)

## old Storable files which we are replacing by HistoneSequencesDB
seq_store = env.PerlScript(
  target = [path4seq("canonical.store"), path4seq("variant.store"),
            path4seq("h1.store")],
  source = path4script("extract_sequences.pl"),
  action = [seq_dir]
)

env.Alias("data", [raw_data, data_store, seq_store])


## TARGET csv_data
##
## This is completely useless and is not required by any other target.
## It is even dangerous because csv is a really poor format for genes.
## We only have this because Andrew wants it for his other projects.
## See https://github.com/af-lab/histone-catalogue/issues/3

def path4csv (name=""):
  return os.path.join(results_dir, "csv", name)

csv_data = env.PerlScript(
  target = [path4csv ("canonical_core_histones.csv"),
            path4csv ("variant_histones.csv"),
            path4csv ("linker_histones.csv")],
  source = path4script ("create_histone_csv.pl"),
  action = [db_store, path4csv()]
)
env.Depends(csv_data, [data_store])
env.Alias("csv", csv_data)


## TARGET update
##
## Remove the previously downloaded data forcing a rebuild.
if "update" in COMMAND_LINE_TARGETS:
  env.AlwaysBuild(raw_data)
  env.Alias("update", raw_data)

## TARGET analysis
##
## For analysis, each script is its own target. We then set an alias that
## groups all of them. Each of these scripts generate a large number of
## files, the targets, we need to make lists of them all

perl_db_var = "HistoneSequencesDB::read_db('%s')" % db_store

clust_targets = list ()
refer_targets = list ()
utr_targets   = list ()

refer_targets += [path4result ("table-reference_comparison.tex")]

utr_targets += [
  path4result ("variables-utr.tex"),
  path4result ("aligned_stem_loops.fasta"),
  path4figure ("seqlogo_stem_loops.eps"),
  path4result ("aligned_HDEs.fasta"),
  path4figure ("seqlogo_HDEs.eps"),
]


analysis = [
  env.PerlOutput(
    target = path4result("table-histone_catalogue.tex"),
    source = path4lib("HistoneCatalogue.pm"),
    M      = ["HistoneCatalogue", "HistoneSequencesDB"],
    eval   = ("HistoneCatalogue::say_histone_catalogue(%s->canonical_core)"
              % perl_db_var),
  ),
  env.PerlOutput(
    target = path4result("table-variant_catalogue.tex"),
    source = path4lib("HistoneCatalogue.pm"),
    M      = ["HistoneCatalogue", "HistoneSequencesDB"],
    eval   = ("HistoneCatalogue::say_histone_catalogue(%s->variants_core)"
              % perl_db_var),
  ),
  env.PerlOutput(
    target = path4result("variables-histone_counts.tex"),
    source = path4lib("HistoneCatalogue.pm"),
    M      = ["HistoneCatalogue", "HistoneSequencesDB"],
    eval   = ("HistoneCatalogue::say_histone_counts(%s)" % perl_db_var),
  ),
  env.PerlOutput(
    target = path4result("variables-cluster_stats.tex"),
    source = path4script("cluster_stats.pl"),
    args   = [db_store],
  ),
  env.PerlOutput(
    target = path4result("variables-protein_stats.tex"),
    source = path4script("protein_stats.pl"),
    args   = [db_store],
  ),
  env.PerlScript(
    target = refer_targets,
    source = path4script ("reference_comparison.pl"),
    action = ["--sequences", seq_dir, "--results", results_dir,
              "--reference", reference_dir],
  ),
  env.PerlOutput(
    target = path4result("histone_insanities.tex"),
    source = path4script("histone_sanity_checks.pl"),
    args   = [db_store],
  ),
  env.PerlScript(
    target = utr_targets,
    source = path4script ("utr_analysis.pl"),
    action = ["--sequences", seq_dir, "--figures", figures_dir,
              "--results", results_dir],
  ),
  env.PerlOutput(
    target = path4result("variables-configuration.tex"),
    source = path4lib("HistoneCatalogue.pm"),
    M      = ["HistoneCatalogue"],
    eval   = "HistoneCatalogue::write_config_variables('%s')" % path4seq("extractor.log")
  ),
  env.PerlOutput(
    target = path4result("table-codon_usage.tex"),
    source = path4script("codon_usage.pl"),
    args   = [db_store],
  ),
]

## The whole mess of alignment targets and their dependencies:
##
##  protein aligns ----> table describing isoforms
##       |           |--> protein sequence logo
##       |           |--> protein align stats -- >alignment percentage identity
##      \_/
##  transcript aligns ----> cds sequence logo
##                     |--> transcript align stats --> dn/ds

cds_aligns = []
protein_aligns = []

for histone in ["H2A", "H2B", "H3", "H4"]:
  protein_align_f = path4result("aligned_%s_proteins.fasta" % histone)
  protein_align = env.PerlScript(
    target = protein_align_f,
    source = path4script("align_proteins.pl"),
    action = [db_store, histone, protein_align_f],
  )
  env.Depends(protein_align, [data_store])
  protein_aligns += protein_align

  cds_align_f = path4result("aligned_%s_cds.fasta" % histone)
  cds_align = env.PerlScript(
    target = cds_align_f,
    source = path4script("align_transcripts.pl"),
    action = [db_store] + protein_align + [cds_align_f],
  )
  env.Depends(cds_align, [data_store] + protein_align)
  cds_aligns += cds_align

  isoforms_desc = env.PerlOutput(
    target = path4result("table-%s-proteins-align.tex" % histone),
    source = path4script("describe_isoforms.pl"),
    args   = [db_store] + protein_align,
  )
  env.Depends(isoforms_desc, protein_align)
  analysis += [isoforms_desc]

  protein_logo_f = path4figure("seqlogo_%s_proteins.eps" % histone)
  cds_logo_f = path4figure("seqlogo_%s_cds.eps" % histone)
  for aln, logo_f in zip ([protein_align, cds_align], [protein_logo_f, cds_logo_f]):
    logo = env.PerlScript(
      target = logo_f,
      source = path4script("mk_histone_seqlogo.pl"),
      action = aln + [logo_f],
    )
    Depends(logo, [aln])
    analysis += [logo]

cds_align_stats = env.PerlOutput(
  target = path4result("variables-align_transcripts_stats.tex"),
  source = path4script("align_transcripts_stats.pl"),
  args   = cds_aligns,
)
Depends(cds_align_stats, cds_aligns)
analysis += [cds_align_stats]

protein_align_stats = env.PerlOutput(
  target = path4result("variables-align_proteins_stats.tex"),
  source = path4script("align_proteins_stats.pl"),
  args   = protein_aligns,
)
Depends(protein_align_stats, protein_aligns)
analysis += [protein_align_stats]


env.Alias ("analysis", analysis)
env.Depends (
  analysis,
  [data_store, seq_store, path4script ("MyLib.pm")]
)

## Our figures, converted to pdf as required for pdflatex
figures = env.PDF(source = Glob(os.path.join(figures_dir, "*.eps")))

## TARGET catalogue
##
## A simpler LaTeX document with the most important tables and figures.
## This works for all organisms, since the actual manuscript is only
## available for the select organisms that we bothered to write.

catalogue = env.PDF(
  target = "catalogue.pdf",
  source = "catalogue.tex"
)
env.Alias("catalogue", catalogue)
Depends(catalogue, [figures, analysis])
env.Default(catalogue)


## TARGET manuscript

manuscript = env.PDF (
  target = "manuscript.pdf",
  source = "manuscript.tex"
)
env.Alias ("manuscript", manuscript)
Depends (manuscript, [figures, analysis])

## Because the manuscript is not built by default, then it's also not
## removed by default by doing `scons -c`.  So we add it to the default
## targets when cleaning
## See http://dcreager.net/2010/01/08/default-scons-clean-targets/
if env.GetOption("clean"):
  env.Default(manuscript)


## TARGET check
##
## Only runs if specified from command line.

if "check" in COMMAND_LINE_TARGETS:
  test_suite = []
  for test_file in env.Glob("t/*.t"):
    unit = env.PerlScript(source=test_file, target=None, action=[])
    test_suite.append(unit)
  check = Alias ("check", [test_suite])


## Build configuration (check if dependencies are all installed)
##
## The really really really right way to do the checks would be to set up a
## scanner that finds the required LaTeX packages and perl modules. But that's
## something that should be done upstream in SCons (the scan for LaTeX source
## is already being worked on so this may not be necessary in the future)

def CheckLaTeXPackage(context, package):
  context.Message("Checking for LaTeX package %s..." % package)
  is_ok = 0 == subprocess.call(["kpsewhich", package + ".sty"],
                               stdout = open(os.devnull, "wb"))
  context.Result(is_ok)
  return is_ok

def CheckLaTeXClass(context, doc_class):
  context.Message("Checking for LaTeX document class %s..." % doc_class)
  is_ok = 0 == subprocess.call(["kpsewhich", doc_class + ".cls"],
                               stdout = open(os.devnull, "wb"))
  context.Result(is_ok)
  return is_ok

def CheckBibTeXStyle(context, style):
  context.Message("Checking for BibTeX style %s..." % style)
  is_ok = 0 == subprocess.call(["kpsewhich", style + ".bst"],
                               stdout = open(os.devnull, "wb"))
  context.Result(is_ok)
  return is_ok

def CheckEmail(context, email):
  context.Message("Checking e-mail address...")
  ## Don't check email validity, the program that uses it will do it
  is_ok = email
  context.Result(is_ok)
  return is_ok

def CheckProg(context, app_name):
  context.Message("Checking for %s..." % app_name)
  is_ok = context.env.WhereIs(app_name)
  context.Result(is_ok)
  return is_ok

def CheckCommand(context, command, message):
  context.Message("Checking %s..." % message)
  is_ok = context.TryAction(command)[0]
  context.Result(is_ok)
  return is_ok

def CheckPerlModule(context, module_name):
  context.Message("Checking for perl module %s..." % module_name)
  is_ok = 0 == subprocess.call(["perl", "-M" + module_name, "-e 1"],
                               stdout = open(os.devnull, "wb"))
  context.Result(is_ok)
  return is_ok

conf = Configure(
  env,
  custom_tests = {
    "CheckLaTeXClass"   : CheckLaTeXClass,
    "CheckLaTeXPackage" : CheckLaTeXPackage,
    "CheckBibTeXStyle"  : CheckBibTeXStyle,
    "CheckPerlModule"   : CheckPerlModule,
    "CheckEmail"        : CheckEmail,
    "CheckProg"         : CheckProg,
    "CheckCommand"      : CheckCommand,
  }
)

## grep -rh '^use ' t/ lib-perl5/ scripts/| sort | uniq
## and then remove the core modules and pragmas
perl_module_dependencies = [
  "Bio::AlignIO",
  "Bio::Align::Utilities",
  "Bio::CodonUsage::Table",
  "Bio::DB::EUtilities",
  "Bio::LocatableSeq",
  "Bio::Root::Version",
  "Bio::Seq",
  "Bio::SeqIO",
  "Bio::SeqUtils",
  "Bio::SimpleAlign",
  "Bio::Tools::CodonTable",
  "Bio::Tools::EUtilities",
  "Bio::Tools::Run::Alignment::Clustalw",
  "Bio::Tools::Run::Alignment::TCoffee",
  "Bio::Tools::Run::Phylo::PAML::Codeml",
  "Bio::Tools::SeqStats",
  "File::Which",
  "Moose",
  "Moose::Util::TypeConstraints",
  "MooseX::StrictConstructor",
  "namespace::autoclean",
  "Statistics::Basic",
  "Test::Exception",
  "Test::Output",
  "Text::CSV",
]

latex_package_dependencies = [
  "fontenc",
  "inputenc",
  "graphicx",
  "url",
  "todonotes",
  "natbib",
  "color",
  "kpfonts",
  "seqsplit",
  "eqparbox",
  "capt-of",
  "hyperref",
  "fp",
  "afterpage",
  "isodate",
  "etoolbox",
  "stringstrings",
  "intcalc",
  "siunitx",
]

env.Help("""
DEPENDENCIES

  Programs:
    * bp_genbank_ref_extractor - Distributed with the perl module
      Bio-EUtilities version 1.74 or later.
    * weblogo - Available at http://weblogo.threeplusone.com/

  Perl modules:
""")
for module in perl_module_dependencies:
  env.Help("    * %s\n" % module)

env.Help("""
  LaTeX document class
    * memoir

  LaTeX packages
""")
for package in latex_package_dependencies:
  env.Help("    * %s\n" % package)

env.Help("""
  BibTeX style
    * agu
""")

## Seriously, this should be the default.  Otherwise, users won't even get
## to see the help text  unless they pass the configure tests.
if not (env.GetOption('help') or env.GetOption('clean')):
  for prog in ["bp_genbank_ref_extractor", "weblogo"]:
    if not conf.CheckProg(prog):
      print ("Unable to find `%s' installed" % prog)
      Exit(1)

  ## We need this option in weblogo to remove the numbering from the X axis.
  ## See issue #33.  This option was added to weblogo version 3.5.0.
  if not conf.CheckCommand("printf '>1\\nAC\\n>2\\nAC\\n'"
                           + " | weblogo --number-interval 50"
                           + " > %s" % os.devnull,
                           "if weblogo supports --number-interval"):
    print "weblogo has no --number-interval option (added in weblogo 3.5.0)"
    Exit(1)

  for module in perl_module_dependencies:
    if not conf.CheckPerlModule(module):
      print "Unable to find perl module %s." % module
      Exit(1)

  for package in latex_package_dependencies:
    if not conf.CheckLaTeXPackage(package):
      print "Unable to find required LaTeX package %s." % package
      Exit(1)

  if not conf.CheckLaTeXClass("memoir"):
    print "Unable to find the LaTeX document class memoir."
    Exit(1)

  if not conf.CheckBibTeXStyle("agu"):
    print "Unable to find the BibTeX style agu."
    Exit(1)

  if not conf.CheckEmail(env.GetOption("email")):
    print ("Per NCBI policy, an email is required when using EUtilities to retrieve data\n"
           "from the Entrez system. Run `scons -h' for details.")
    Exit(1)

env = conf.Finish()
