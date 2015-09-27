#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path
import subprocess

env = Environment()
env.Append(PERL5LIB=['lib-perl5'])
env.Tool('perl5')

## FIXME very temporary while we move all of MyLib to our modules
env.Append(PERL5LIB=['scripts'])

## FIXME temporary while we create our perl5 builders, and clean how the
##        scripts create output.
perl_inc = ["-I" + inc for inc in env['PERL5LIB']]
perl_args = ["perl"] + perl_inc
perl_command = " ".join(perl_args)

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

  --email=address
      Set email to be used when connecting to the NCBI servers. This can
      be anything that conforms to RCF822. The following are valid:

          scons --email="Your Name <your.name@domain.here>"
          scons --email="<your.name@domain.here>"

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

if not env.GetOption("verbose"):
  env.AppendUnique(PDFLATEXFLAGS  = "-interaction=batchmode")
  env.AppendUnique(PDFTEXFLAGS    = "-interaction=batchmode")
  env.AppendUnique(TEXFLAGS       = "-interaction=batchmode")
  env.AppendUnique(LATEXFLAGS     = "-interaction=batchmode")
  env.AppendUnique(BIBTEXFLAGS    = "--terse")  # some ports of BibTeX may use --quiet instead

env.Help("""
TARGETS

  data
    Connect to the Entrez system to download new sequences. This data is
    required for the analysis. To prevent conflicts during analysis,
    previously downloaded data will first be removed.

  analysis
    Run all the scripts to analyse the data such as: sequence alignments,
    search for anomalies on the sequence annotations, LaTeX tables listing
    all genes and proteins, sequence differences between isoforms. The
    analysis results are required for the publication.

  publication
    Build PDF for publication.

""")

scripts_dir   = "scripts"
results_dir   = "results"
figures_dir   = "figs"
seq_dir       = os.path.join (results_dir, "sequences")
reference_dir = os.path.join ("data", "reference-Marzluff_2002")

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
data = env.Command (
  target = [path4seq ("data.csv"),
            path4seq ("canonical.csv"), path4seq ("canonical.store"),
            path4seq ("variant.csv"), path4seq ("variant.store"),
            path4seq ("h1.csv"), path4seq ("h1.store")],
  source = path4script ("extract_sequences.pl"),
  action = "%s $SOURCE --email %s %s" % (perl_command, env.GetOption ("email"), seq_dir)
)
env.Alias("data", data)
env.Clean(data, seq_dir)


## TARGET analysis
##
## For analysis, each script is its own target. We then set an alias that
## groups all of them. Each of these scripts generate a large number of
## files, the targets, we need to make lists of them all

align_targets = list ()
clust_targets = list ()
prot_targets  = list ()
refer_targets = list ()
check_targets = list ()
utr_targets   = list ()
var_targets   = list ()
check_targets = list ()

for histone in ["H2A", "H2B", "H3", "H4"]:
  align_targets += [
    path4result ("aligned_%s_proteins.fasta"   % histone),
    path4figure ("seqlogo_%s_proteins.eps"     % histone),
    path4result ("aligned_%s_cds.fasta"        % histone),
    path4figure ("seqlogo_%s_cds.eps"          % histone),
    path4result ("table-%s-proteins-align.tex" % histone),
  ]
align_targets += [path4result ("variables-align_results.tex")]

clust_targets += [
  path4result ("table-histone_catalogue.tex"),
  path4result ("table-histone_catalogue.csv"),
  path4result ("variables-cluster_stats.tex"),
]

prot_targets += [path4result ("variables-protein_stats.tex")]

refer_targets += [path4result ("table-reference_comparison.tex")]

check_targets += [path4result ("histone_insanities.tex")]

utr_targets += [
  path4result ("variables-utr.tex"),
  path4result ("aligned_stem_loops.fasta"),
  path4figure ("seqlogo_stem_loops.eps"),
  path4result ("aligned_HDEs.fasta"),
  path4figure ("seqlogo_HDEs.eps"),
]

var_targets += [path4result ("table-variant_catalogue.tex")]

check_targets += [path4result ("histone_insanities.tex")]

analysis = [
  env.Command (
    target = align_targets,
    source = path4script ("align_sequences.pl"),
    action = "%s $SOURCE --sequences %s --figures %s --results %s" % (perl_command, seq_dir, figures_dir, results_dir)
  ),
  env.Command (
    target = clust_targets,
    source = path4script ("cluster_stats.pl"),
    action = "%s $SOURCE --sequences %s -results %s" % (perl_command, seq_dir, results_dir)
  ),
  env.Command (
    target = prot_targets,
    source = path4script ("protein_stats.pl"),
    action = "%s $SOURCE --sequences %s --results %s" % (perl_command, seq_dir, results_dir)
  ),
  env.Command (
    target = refer_targets,
    source = path4script ("reference_comparison.pl"),
    action = "%s $SOURCE --sequences %s --results %s --reference %s" % (perl_command, seq_dir, results_dir, reference_dir)
  ),
  env.Command (
    target = check_targets,
    source = path4script ("histone_sanity_checks.pl"),
    action = "%s $SOURCE --sequences %s --results %s" % (perl_command, seq_dir, results_dir)
  ),
  env.Command (
    target = utr_targets,
    source = path4script ("utr_analysis.pl"),
    action = "%s $SOURCE --sequences %s --figures %s --results %s" % (perl_command, seq_dir, figures_dir, results_dir)
  ),
  env.Command (
    target = var_targets,
    source = path4script ("variants.pl"),
    action = "%s $SOURCE --sequences %s --results %s" % (perl_command, seq_dir, results_dir)
  ),
  env.Command(
    target = path4result ("variables-configuration.tex"),
    source = path4seq ("extractor.log"),
    action = ("%s -MHistoneCatalogue -e " % (perl_command)
              + "\"HistoneCatalogue::write_config_variables "
              + "(\'$TARGET\', \'$SOURCE\')\"")
  ),
]

env.Alias ("analysis", analysis)
env.Depends (
  analysis,
  [data, path4script ("MyLib.pm")]
)


## TARGET manuscript
##
## Dependent on the figures being converted into PDF.
figures = env.PDF (source = Glob (os.path.join (figures_dir, "*.eps")))

manuscript = env.PDF (
  target = "histone_catalog.pdf",
  source = "histone_catalog.tex"
)
env.Alias ("manuscript", manuscript)
Depends (manuscript, [figures, analysis])

env.Default(manuscript)


## TARGET check
##
## Only runs if specified from command line.

if "check" in COMMAND_LINE_TARGETS:
  def build_test_unit (test_file):
    test_args = perl_args + [str(test_file)]
    ## target must be "None" so it always builds (an empty string is always
    ## up to date so that's no good).  And we need source, otherwise it
    ## would also be always up to date.
    test_unit = env.Command (
      source = test_file,
      target = None,
      action = (lambda target="fff", source="fff", env="fff"
        : subprocess.call (test_args, env = env['ENV']))
    )
    return test_unit

  test_suite = []
  for tf in env.Glob("t/*.t"):
    test_suite.append(build_test_unit(tf))
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
    "CheckProg"         : CheckProg
  }
)

perl_module_dependencies = [
  "Email::Valid",
  "Bio::Root::Version",
  "Bio::Tools::EUtilities",
  "Bio::SeqIO",
  "Bio::Tools::Run::Alignment::TCoffee",
  "Text::CSV",
  "Statistics::Basic",
  "Moose",
  "Moose::Util::TypeConstraints",
  "Test::More",
  "Test::Exception",
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


## Seriously, this should be the default.  Otherwise, users won't even get
## to see the help text  unless they pass the configure tests.
if not env.GetOption('help'):
  for prog in ["bp_genbank_ref_extractor", "weblogo"]:
    if not conf.CheckProg(prog):
      print ("Unable to find `%s' installed" % prog)
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

  if not conf.CheckEmail(env.GetOption("email")):
    print ("Per NCBI policy, an email is required when using EUtilities to retrieve data\n"
           "from the Entrez system. Run `scons -h' for details.")
    Exit(1)

env = conf.Finish()
