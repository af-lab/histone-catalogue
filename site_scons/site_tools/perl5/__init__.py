#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Copyright (C) 2011-2015 Carnë Draug <carandraug+dev@gmail.com>
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, see
## <http:##www.gnu.org/licenses/>.

"""SCons tool for perl5 scripts and modules
"""

import os
import os.path
import subprocess
import distutils.spawn

import SCons.Script
import SCons.Errors

def get_perl_I(env):
  return ["-I" + str(inc) for inc in env['PERL5LIB']]


##
## PerlScript
##

def create_args_perl5_script(source, env):
  return (['perl'] + get_perl_I(env)
          + [str(source[0])]
          + [str(a) for a in env['action']])

def perl5_script_strfunc(target, source, env):
  args = create_args_perl5_script(source, env)
  cmd = "$ %s " % (args[0])
  cmd += " ".join(["'%s'" % (arg) for arg in args[1:]])
  return cmd

def perl5_script(target, source, env):
  """The actual PerlScript Builder action.
  """
  args = create_args_perl5_script(source, env)
  return subprocess.call(args)


##
## PerlSub
##

def create_args_perl5_sub(source, env):
  f_name = os.path.basename(str(source[0]))
  m_name = os.path.splitext(f_name)[0]
  return ['perl'] + get_perl_I(env) + ["-M"+m_name, "-e", env.get('action')]

def perl5_sub_strfunc(target, source, env):
  args = create_args_perl5_sub(source, env)
  cmd = "$ %s " % (args[0])
  cmd += " ".join(["'%s'" % (arg) for arg in args[1:]])
  return cmd

def perl5_sub(target, source, env):
  """The actual PerlSub Builder action.
  """
  args = create_args_perl5_sub(source, env)
  return subprocess.call(args)


##
## PerlOutput
##

def create_args_perl5_output(env, source):
  """Create list of arguments for the perl command.

  Kinda silly that we need to do this twice, once to generate the
  arguments for subprocess, and another for the displayed string.
  """
  args = ['perl'] + get_perl_I(env)
  if env.has_key('M'):
    args += ["-M" + m for m in env['M']]
  if env.has_key('eval'):
    args += ["-e", env.get('eval')]
  else:
    args += [str(source[0])]
    if env.has_key('args'):
      args += [str(a) for a in env.get('args')]
  return args

def perl5_output_strfunc(target, source, env):
  """ Create string to be displayed for the builder.
  """
  args = create_args_perl5_output(env, source)
  cmd = "$ %s " % (args[0])
  cmd += " ".join(["'%s'" % (arg) for arg in args[1:]])
  cmd += " > " + str(target[0])
  return cmd

def perl5_output(target, source, env):
  """The actual PerlOutput Builder action.
  """
  args = create_args_perl5_output(env, source)
  with open(str(target[0]), "w") as target_file:
    rv = subprocess.call(args, stdout=target_file)
  return rv



##
## Scanner
##

def perl5_scanner(node, env, path):
  perl_args = ["perl", "-MModule::ScanDeps"] + get_perl_I(env)

  ## We handle the recursion ourselves instead of using the recursion
  ## option from Module::ScanDeps, because that will descend deep into
  ## the whole system libraries and our dependencies (from the project),
  ## are likely to be much smaller.
  our_deps = set()
  cwd = os.path.realpath(os.getcwd())
  def scan_deps(fpaths):
    my_deps = set()
    for fpath in fpaths:
      ## note that perl sorts the keys because the order matters (it will
      ## rebuild the target even if only the order changed), and the order
      ## of keys in a hash is random.
      code = """
        $deps = scan_deps (files => ['%s'],
                           recurse => 0);
        foreach (sort keys %%{$deps})
          { print $deps->{$_}->{file} . "\\n"; }
      """ % str(fpath)
      deps = subprocess.check_output(perl_args + ["-e", code], env=env['ENV'])
      for dep in deps.splitlines():
        dep = os.path.realpath(dep)
        if dep.startswith(cwd) and dep not in our_deps:
          my_deps.add(dep)

    old_len = len(our_deps)
    my_deps.difference_update(our_deps)
    our_deps.update(my_deps)
    if len(our_deps) != old_len:
      scan_deps(my_deps)

  scan_deps([node])
  return env.File(list(our_deps))


##
## Configure checks
##

## This is currently unused because there's no way to share configure tests.
def CheckPerlModule(context, module_name):
  context.Message("Checking for perl module %s..." % module_name)
  is_ok = 0 == subprocess.call(["perl", "-M" + module_name, "-e 1"],
                               stdout = open(os.devnull, "wb"))
  context.Result(is_ok)
  return is_ok



def generate(env):
  vars = SCons.Script.Variables()
  vars.Add('PERL5LIB', """
List of directories in which to look for Perl library files before
looking in the standard library and the current directory.  It is
akin to the environment variable of the same name.
""")
  vars.Update(env)
  SCons.Script.Help(vars.GenerateHelpText(env))

  scanner = SCons.Script.Scanner(perl5_scanner, skeys=['.pl', '.pm'])
  env.Append(SCANNERS=scanner)

  perl5_output_action = SCons.Script.Action(perl5_output,
                                            strfunction=perl5_output_strfunc)
  bld_perl5_output = SCons.Script.Builder(action=perl5_output_action)
  env.Append(BUILDERS={'PerlOutput' : bld_perl5_output})

  perl5_script_action = SCons.Script.Action(perl5_script,
                                            strfunction=perl5_script_strfunc)
  bld_perl5_script = SCons.Script.Builder(action=perl5_script_action)
  env.Append(BUILDERS={'PerlScript' : bld_perl5_script})

  perl5_sub_action = SCons.Script.Action(perl5_sub,
                                         strfunction=perl5_sub_strfunc)
  bld_perl5_sub = SCons.Script.Builder(action=perl5_sub_action)
  env.Append(BUILDERS={'PerlSub' : bld_perl5_sub})

def exists(env):
  if not distutils.spawn.find_executable("perl"):
    return 0

  rc = subprocess.call(["perl", "-MModule::ScanDeps", "-e", "1;"])
  return rc == 0
