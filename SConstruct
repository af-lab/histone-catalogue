env = Environment()
paper = env.PDF(target = 'histone_catalog.pdf', source = 'histone_catalog.tex')

Depends(paper, ['references.bib',
                ])
Decider('content')

def BuildPackageTest(package):
    package_test_text = "\documentclass{article}\n\usepackage{%s}\n\\begin{document}\n\end{document}" % package
    def PackageTest(context):
        context.Message("Checking for LaTeX package %s..." % package)
        try:
            b = context.env.PDF
        except AttributeError:
            return False
        is_ok = context.TryBuild(b, package_test_text, '.tex')
        context.Result(is_ok)
        return is_ok
    return PackageTest

required_packages = ["graphicx", "multirow", "url", "todonotes", "natbib", "hyperref", "ifdraft", "palatino"]
package_tests = dict()
for package in required_packages:
    package_tests[package] = BuildPackageTest(package)

conf = Configure(env, package_tests)

for package in required_packages:
    if not getattr(conf, package)():
        print "Unable to find required LaTeX package %s" % package
        Exit(1)

env = conf.Finish()
