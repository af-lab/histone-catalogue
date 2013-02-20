## coding: utf-8
import os
import re
import subprocess

## reuse the building instructions already written on the README file
def readme_help():
    readme  = open("README.md", "read").read()
    ini_str = ("Building instructions\n"
               "---------------------\n")
    end_str = "Run `scons -h` for more details on the building targets and options\."
    match   = re.search("(?<=^%s).*(?=%s)" % (ini_str, end_str),
                        readme,
                        re.M | re.S | re.I)
    if match:
      return match.group()
    else:
      return ""
Help(readme_help())
Help("""
Each target may require specific tools to be installed in your system, or
certain options to be specified.

OPTIONS

    --email=address
            Set email to be used when connecting to the NCBI servers. This can
            be anything that conforms to RCF822. The following is valid:
            
                scons --email="Your Name <your.name@domain.here>"
    
    --verbose
            LaTeX and BibTeX compilers are silenced by default using the
            batchmode and terse options. Set this option to revert it.

""")

AddOption("--email",
          type = "string",
          dest = "email",
          help = "E-mail provided to NCBI when connecting to Entrez."
          )
AddOption("--verbose",
          action  = "store_true",
          dest    = "verbose",
          default = False,
          help    = "Print LaTeX and BibTeX output."
          )

## Set environment
env = Environment()
if not GetOption("verbose"):
    env.AppendUnique(PDFLATEXFLAGS  = "-interaction=batchmode")
    env.AppendUnique(PDFTEXFLAGS    = "-interaction=batchmode")
    env.AppendUnique(TEXFLAGS       = "-interaction=batchmode")
    env.AppendUnique(LATEXFLAGS     = "-interaction=batchmode")
    env.AppendUnique(BIBTEXFLAGS    = "--terse")  # some ports of BibTeX may use --quiet instead

Help("""
TARGETS

    data
    
        Connect to the Entrez system to download new sequences. This data is
        required for the analysis. To prevent conflicts during analysis,
        previously downloaded data will first be removed.
    
    analysis
    
        Run all the scripts to analyse the data such as: sequence alignments,
        search for anomalies on the sequence annotations, LaTeX tables listing
        all genes and proteins, sequence differences between isoforms. The
        analysis results are required for the report and publication.
    
    report
    
        Build short PDF with the tables and figures generated from analysis.
    
    publication
    
        Build PDF for publication.

""")

## scons does not like it when target is a directory. It will always consider it
## up to date, even the source changes. So we use the data.csv file as target
datapath    = os.path.join("results", "sequences")
datatarget  = os.path.join(datapath, "data.csv")
env.Alias("data", datatarget)
data = env.Command(target = datatarget,
                   source = "scripts/extract_sequences.pl",
                   action = "$SOURCES --email %s" % GetOption('email'))
env.Clean(data, datapath)

## we will probably need to create our own builder for this
env.Alias("analysis",  "")
analysis = env.Command(target = [],
                       source = [],
                       action = [])

env.Alias("report", "report.pdf")
report = env.PDF(target = "report.pdf",
                 source = "report.tex")

env.Alias("publication", "histone_catalog.pdf")
publication = env.PDF(target = "histone_catalog.pdf",
                      source = "histone_catalog.tex")


## The really really really right way to do the checks would be to set up a
## scanner that finds the required LaTeX packages and perl modules. And that's
## something that should be done upstream in SCons.

def CheckLaTeXPackage(context, package):
    test_text = ("\documentclass{article}\n"
                 "\usepackage{%s}\n"
                 "\\begin{document}\n"
                 "dummy text\n"
                 "\end{document}") % package
    context.Message("Checking for LaTeX package %s..." % package)
    is_ok = context.TryBuild(context.env.PDF, test_text, ".tex")
    context.Result(is_ok)
    return is_ok

def CheckLaTeXClass(context, doc_class):
    test_text = ("\documentclass{%s}\n"
                 "\\begin{document}\n"
                 "dummy text\n"
                 "\end{document}") % doc_class
    context.Message("Checking for LaTeX document class %s..." % doc_class)
    is_ok = context.TryBuild (context.env.PDF, test_text, ".tex")
    context.Result(is_ok)
    return is_ok

def CheckEmail(context, email):
    context.Message("Checking e-mail address...")
    is_ok = False
    if email:
      ## leave the actual email validation to the script. email.utils.parseaddr
      ## does not actually make sure the address is valid
      is_ok = email
    context.Result(is_ok)
    return is_ok

def CheckApp(context, app_name):
    context.Message("Checking for %s..." % app_name)
    is_ok = context.env.WhereIs(app_name)
    if not is_ok:
        ## because https://bitbucket.org/scons/scons/pull-request/67
        is_ok = False
    context.Result(is_ok)
    return is_ok

def CheckPerlModule(context, module_name):
    context.Message("Checking for perl module %s..." % module_name)
    is_ok = True
    if (subprocess.call(["perl", "-M%s" % module_name, "-e 1"],
                        stderr = open(os.devnull, "wb"))):
      is_ok = False
    context.Result(is_ok)
    return is_ok

conf = Configure(env,
                 custom_tests = {
                                 "CheckLaTeXClass"   : CheckLaTeXClass,
                                 "CheckLaTeXPackage" : CheckLaTeXPackage,
                                 "CheckPerlModule"   : CheckPerlModule,
                                 "CheckEmail"        : CheckEmail,
                                 "CheckApp"          : CheckApp,
                                 })

Help("""
DEPENDENCIES

    for DATA
    
        * bp_genbank_ref_extractor - Distributed with bioperl's Bio-EUtilities
          version 1.73 or later.
        * email address set with --email
        * perl module Valid::Email
    
    for ANALYSIS
        
        * weblogo - Available at http://weblogo.threeplusone.com/
        * the following perl modules:
              * Bio::SeqIO
              * Bio::Tools::Run::Alignment::TCoffee
              * Lingua::EN::Numbers
              * Text::CSV
    
    for REPORT and PUBLICATION
        * LaTeX document class memoir
        * the following LaTeX packages:
              * url
              * todonotes
              * natbib
              * palatino

""")

if "data" in map (str, BUILD_TARGETS):
    if not conf.CheckApp("bp_genbank_ref_extractor"):
        print ("Unable to find `bp_genbank_ref_extractor' installed which is required to\n"
               "download new sequences. It is distributed with bioperl's Bio-EUtilities.\n")
        Exit(1)
    
    if not conf.CheckEmail(GetOption("email")):
        print ("Per NCBI policy, an email is required when using EUtilities to retrieve data\n"
               "from the Entrez system. Run `scons -h' for details.")
        Exit(1)
    
    required_module = "Email::Valid"
    if not conf.CheckPerlModule(required_module):
        print "Unable to find perl module %s." % required_module
        Exit(1)

if "analysis" in map (str, BUILD_TARGETS):
    if not conf.CheckApp("weblogo"):
        print "Unable to find weblogo installed."
        Exit(1)
    
    required_modules = [
                        "Bio::SeqIO",
                        "Bio::Tools::Run::Alignment::TCoffee",
                        "Lingua::EN::Numbers",
                        "Text::CSV",
                        ]
    for module in required_modules:
        if not conf.CheckPerlModule(module):
            print "Unable to find perl module %s." % module
            Exit(1)

if ("publication" or "report") in map (str, BUILD_TARGETS):
    if not conf.CheckLaTeXClass("memoir"):
        print "Unable to find the LaTeX document class memoir."
        Exit(1)
    
    required_packages = [
                         "fontenc",
                         "graphicx",
                         "url",
                         "todonotes",
                         "natbib",
                         "palatino",
                          ]
    for package in required_packages:
        if not conf.CheckLaTeXPackage(package):
            print "Unable to find required LaTeX package %s." % package
            Exit(1)

env = conf.Finish()
