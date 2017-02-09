#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path
import subprocess

env = Environment()
env.Append(PERL5LIB=[env.Dir('lib-perl5')])
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
  ## Some args may be File or Dir which need to be converted to str.
  args = [str(x) for x in env['action']]
  return subprocess.call(args)
def no_shell_command_strfunc(target, source, env):
  args = env['action']
  return "$ %s " % (args[0]) + " ".join(["'%s'" % (arg) for arg in args[1:]])
no_shell_command_action = Action(no_shell_command, strfunction=no_shell_command_strfunc)
env.Append(BUILDERS={'NoShellCommand' : Builder(action=no_shell_command_action)})


## FIXME very temporary while we move all of MyLib to our modules
env.Append(PERL5LIB=[env.Dir('scripts')])

env.Help("""
TARGETS

  catalogue (default)
    Build pdf for catalogue --- the pdf with tables of all genes
    and their products, counts of genes per cluster, sequence
    alignment, and list of anomalies.

  data
    Get raw data required for the analysis from the Entrez system
    This will not download new data if there is only data already
    downloaded.  See target 'update' for that.

  update
    Remove previously downloaded raw data and download fresh one.
    It will not rebuild a catalogue or manuscript unless those targets
    are also specified.

  csa_data
    Prepare a csv file with the histone gene information downloaded
    from the Entrez gene database.

  analysis
    Run all the scripts to analyse the data such as: sequence alignments,
    search for anomalies on the sequence annotations, LaTeX tables listing
    all genes and proteins, sequence differences between isoforms. The
    analysis results are required for the publication.

  manuscript
    Build manuscript.pdf.


OPTIONS

  --email=ADDRESS
      Set email to be used when connecting to the NCBI servers.:

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
  default = "",
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

def CheckVariable(context, variable_name):
  context.Message("Checking for variable %s..." % variable_name)
  ## Because env is not really a dict, we can't use 'is in env'
  ## and need to resort to a try catch.
  try:
    context.env[variable_name]
    is_ok = True
  except KeyError:
    is_ok = False
  context.Result(is_ok)
  return is_ok

## Many modules bioperl-run are mainly a wrapper to executable.  Checking
## for the presence of the module is not enough, we need to check if they
## are working (well, we just check if bioperl finds the programs).
def CheckBioperlRunExecutable(context, module):
  context.Message("Checking for working %s..." % module)
  command = ("perl -M%s -e '%s->new()->executable()'" % (module, module))
  is_ok = context.TryAction(command)[0]
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
    "CheckVariable"     : CheckVariable,
    "CheckBioperlRunExecutable" : CheckBioperlRunExecutable,
  }
)


## this is needed by the scons perl tool
perl_dependencies = [
  "Module::ScanDeps",
]

## grep -rh '^use ' lib-perl5/ scripts/| sort | uniq
## and then remove the core modules and pragmas
perl_analysis_dependencies = [
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
  "Bio::Tools::SeqStats",
  "File::Which",
  "Moose",
  "Moose::Util::TypeConstraints",
  "MooseX::StrictConstructor",
  "namespace::autoclean",
  "Statistics::Basic",
  "Text::CSV",
]

## grep -rh '^use ' t/ | sort | uniq
## and then remove the core modules and pragmas
perl_test_dependencies = [
  "Bio::LocatableSeq",
  "Bio::Seq",
  "Bio::SimpleAlign",
  "Test::Exception",
  "Test::Output",
]

## This is a dict where the key is perl module and value the likely program
## (or suite of programs) that is likely to be missing to make it work.
bioperl_run_dependencies = {
  "Bio::Tools::Run::Alignment::Clustalw" : "clustalw",
  "Bio::Tools::Run::Alignment::TCoffee" : "t-coffee",
  "Bio::Tools::Run::Phylo::PAML::Codeml" : "PAML",
}

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
for module in (perl_dependencies + perl_analysis_dependencies
               + bioperl_run_dependencies.keys()):
  env.Help("    * %s\n" % module)

env.Help("""
  Perl modules for test suite:
""")
for module in perl_test_dependencies:
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

## Seriously, this should be the default.  Otherwise, users won't even
## get to see the help text  unless they pass the configure tests.
## And Configure(..., clean=False,help=False) does not really work,
## it just makes all configure tests fail.
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

  for module in set (perl_dependencies + perl_analysis_dependencies
                     + bioperl_run_dependencies.keys()):
    if not conf.CheckPerlModule(module):
      print "Unable to find perl module %s." % module
      Exit(1)

  for module, program in bioperl_run_dependencies.iteritems():
    if not conf.CheckBioperlRunExecutable(module):
      print "bioperl's %s is not working (did you install %s?)" % (module, program)
      Exit(1)

  if "check" in COMMAND_LINE_TARGETS:
    for module in set (perl_test_dependencies):
      if not conf.CheckPerlModule(module):
        print "Unable to find perl module %s." % module
        Exit(1)

  if not conf.CheckVariable('EPSTOPDF'):
    print "SCons EPSTOPDF not configured.  Do you have epstopdf installed (part of texlive)"
    Exit(1)

  ## We need this so we can then use CheckLatex* or the user gets a
  ## pretty cryptic error message.
  if not conf.CheckProg("kpsewhich"):
    print "Unable to find `kpsewhich' (part of texlive) installed"
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

  ## If the users does not want to set an email, let him.  We warn
  ## here and Bio-EUtilities warns again, but don't force it.
  if not conf.CheckEmail(env.GetOption("email")):
    print ("WARNING: Per NCBI policy, an email is required when using EUtilities\n"
           "         to retrieve data from the Entrez system.  Consider using\n"
           "         '--email' and see README for details why.")

env = conf.Finish()

##
## Actual TARGETS from this point on
##

scripts_dir   = env.Dir ("scripts")
results_dir   = env.Dir ("results")
figures_dir   = env.Dir ("figs")
seq_dir       = env.Dir (os.path.join (str (results_dir), "sequences"))

reference_dirname = "reference-" + env.GetOption("organism").replace(" ", "-").lower()
reference_dir = env.Dir (os.path.join ("data", reference_dirname))

def path4lib(name):
  return os.path.join("lib-perl5", name)
def path4script (name):
  return os.path.join (str (scripts_dir), name)
def path4result (name):
  return os.path.join (str (results_dir), name)
def path4figure (name):
  return os.path.join (str (figures_dir), name)
def path4seq (name):
  return os.path.join (str (seq_dir), name)


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
env.AddPreAction(raw_data, Delete (seq_dir))
env.Clean(raw_data, seq_dir)

db_store = File (path4result("histones_db.store"))
data_store = env.PerlSub(
  target = db_store,
  source = path4lib("HistoneSequencesDB.pm"),
  action = ('HistoneSequencesDB->new("%s")->write_db("%s")'
            %(seq_dir.path, db_store.path))
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
  return os.path.join(str (results_dir), "csv", name)

csv_data = env.PerlScript(
  target = [path4csv ("canonical_core_histones.csv"),
            path4csv ("variant_core_histones.csv"),
            path4csv ("linker_histones.csv")],
  source = path4script ("create_histone_csv.pl"),
  action = [db_store, Dir (path4csv()).path]
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

perl_db_var = "HistoneSequencesDB::read_db('%s')" % db_store.path

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
    eval   = ("HistoneCatalogue::write_config_variables('%s')"
              % File (path4seq("extractor.log")).path)
  ),
  env.PerlOutput(
    target = path4result("table-codon_usage.tex"),
    source = path4script("codon_usage.pl"),
    args   = [db_store],
  ),
]

if os.path.isdir(str (reference_dir)):
  analysis.append (env.PerlScript(
    target = refer_targets,
    source = path4script ("reference_comparison.pl"),
    action = ["--sequences", seq_dir, "--results", results_dir,
              "--reference", reference_dir],
  ))
else:
  ## This is only needed for manuscript.pdf anyway.  If we ever get to
  ## support a manuscript for multiple organisms, it may be possible
  ## that they will have no reference, so we will have to figure out a
  ## way to handle this better then.
  print ("WARNING: no reference data found for %s.\n"
         "         Skipping comparison against reference."
         % env.GetOption("organism"))

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
  protein_align_f = File (path4result("aligned_%s_proteins.fasta" % histone))
  protein_align = env.PerlScript(
    target = protein_align_f,
    source = path4script("align_proteins.pl"),
    action = [db_store, histone, protein_align_f],
  )
  env.Depends(protein_align, [data_store])
  protein_aligns += protein_align

  cds_align_f = File (path4result("aligned_%s_cds.fasta" % histone))
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

  protein_logo_f = File (path4figure("seqlogo_%s_proteins.eps" % histone))
  cds_logo_f = File (path4figure("seqlogo_%s_cds.eps" % histone))
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
  [data_store, seq_store, File (path4script ("MyLib.pm"))]
)

## Our figures, converted to pdf as required for pdflatex
figures = env.PDF(source = Glob(os.path.join(str(figures_dir), "*.eps")))

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
