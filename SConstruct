paper = PDF(target = 'histone_catalog.pdf', source = 'histone_catalog.tex')

Depends(paper, ['library.bib',
                'figs/nomenclature-schematic.pdf',
                'figs/H2A_weblogo.pdf', 'figs/H2B_weblogo.pdf',
                'results/tables.tex', 'results/variables.tex'])
Decider('content')  # same as MD5

## we'll use inkscape to generate pdfs of image from their svgs
SVGS = """
figs/nomenclature-schematic
figs/H2A_weblogo
figs/H2B_weblogo
figs/H3_weblogo
figs/H4_weblogo
""".split()

for SVG in SVGS:
  ## should skip if svg has not changed since last time
  x    = Command(SVG + '.pdf', SVG +'.svg', 'inkscape `basename ' + SVG + '.svg` --export-pdf=`basename ' + SVG + '.pdf`', chdir=1)

## should create option to get data again

## should create option to recreate tables and figs from data

## specific option to reconvert svg and add inkscape as dependency

## should create option to do everything

## should mark Latex dependencies:
##      * bibtex style
##      * url package
##      * todonotes package
##      * graphicx

## dependencies for perl modules, only if running the specific sections so that building paper from data does not need bioperl

