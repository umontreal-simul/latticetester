#! /usr/bin/env python
# encoding: utf-8

top = '.'
out = 'build'

from waflib import Utils

import imp
def waftool(name):
    return imp.load_module('waf_' + name, *imp.find_module(name, ['./waftools']))

NUMTYPES = {'LLDD': 1, 'ZZDD': 2, 'ZZRR': 3}

version = waftool('version')
compiler = waftool('compiler')
deps = waftool('deps')

def options(ctx):
    # ctx.recurse('latticetester')
    ctx.load('compiler_c compiler_cxx gnu_dirs waf_unit_test')
    ctx.add_option('--link-static', action='store_true', help='statically link with Boost and NTL')
    ctx.add_option('--build-docs', action='store_true', default=False, help='build documentation')
    ctx.add_option('--boost', action='store', help='prefix under which Boost is installed')
    ctx.add_option('--ntl', action='store', help='prefix under which NTL is installed')
    ctx.add_option('--ntltypes', action='store', help='set NTL integer/real types: {}'.format(', '.join(NUMTYPES)))
    ctx.add_option('--gmp', action='store', help='prefix under which GMP is installed')
    ctx.add_option('--testu01', action='store', help='prefix under which testu01 is installed')


def configure(ctx):
    build_platform = Utils.unversioned_sys_platform()
    ctx.msg("Build platform", build_platform)

    # detect NTL types
    if ctx.options.ntltypes in NUMTYPES:
        ctx.define('WITH_NTL', 1)
        ctx.define('NTL_TYPES_CODE', NUMTYPES[ctx.options.ntltypes])
        ctx.env.LATTICETESTER_SUFFIX = str(ctx.options.ntltypes)
    elif ctx.options.ntltypes is None:
        ctx.env.LATTICETESTER_SUFFIX = ''
    else:
        ctx.fatal("invalid NTL number types {}; must be one of: {}".format(ctx.options.ntltypes, ', '.join(NUMTYPES)))

    ctx.msg("Use NTL", ctx.options.ntltypes and "yes" or "no")
    if ctx.options.ntltypes:
        ctx.msg("NTL integer/real types", ctx.options.ntltypes)

    ctx.load('compiler_c compiler_cxx gnu_dirs waf_unit_test')
    ctx.check(features='cxx', cxxflags='-std=c++11')
    ctx.env.append_unique('CXXFLAGS', ['-std=c++11', '-Wall'])
    ctx.check(features='c', cflags='-std=c99')
    ctx.env.append_unique('CFLAGS', ['-std=c99', '-Wall'])

    # suppress Boost ublas warnings
    compiler.add_cxx_flag_if_supported(ctx, '-Wno-unused-local-typedefs')
    compiler.add_cxx_flag_if_supported(ctx, '-Wno-unused-function')
    compiler.add_cxx_flag_if_supported(ctx, '-Wnon-virtual-dtor')
    compiler.add_cxx_flag_if_supported(ctx, '-Wshorten-64-to-32') # clang

    if ctx.options.link_static:
        #flags = ['-static', '-static-libgcc', '-static-libstdc++']
        flags = ['-static-libgcc', '-static-libstdc++']
        if ctx.check(features='cxx cxxprogram',
                linkflags=flags,
                mandatory=False):
            ctx.env.append_unique('LINKFLAGS', flags)

    ctx.version_file()

    # options
    if ctx.options.boost:
        deps.add_deps_path(ctx, 'Boost', ctx.options.boost)
    if ctx.options.ntl:
        deps.add_deps_path(ctx, 'NTL', ctx.options.ntl)
    if ctx.options.testu01:
        deps.add_deps_path(ctx, 'TestU01', ctx.options.testu01)

    ctx_check = deps.shared_or_static(ctx, ctx.check)

    if ctx.options.ntltypes:
        # NTL
        ctx_check(features='cxx cxxprogram', header_name='NTL/vector.h')
        ctx_check(features='cxx cxxprogram', lib='ntl', uselib_store='NTL')

    boost_version = (1,57,0)
    boost_version_str = '.'.join(str(x) for x in boost_version)
    ctx_check(features='cxx',
            msg='Checking for Boost version >= %s' % boost_version_str,
            fragment=
            '#include <boost/version.hpp>\n'
            '#if BOOST_VERSION < %d\n'
            '#error "Boost >= %s is required"\n'
            '#endif' %
            (boost_version[0] * 100000 + boost_version[1] * 100 + boost_version[2], boost_version_str))


    ctx.env['VERSION'] = 1.1
    #ctx.version_file()
    ctx.define('LATTICETESTER_VERSION', ctx.set_version())
    ctx.msg("Setting LatticeTester version", version.VERSION)


    if not ctx.env.LATTICETESTER_SUFFIX:
        ctx.fatal('Please specify the NTL number types with --ntltypes (try ./waf --help)')

    # GMP
    ctx_check(features='cxx cxxprogram', header_name='gmp.h')
    ctx_check(features='cxx cxxprogram', lib='gmp', uselib_store='GMP')

    # TestU01
    ctx_check(features='cxx cxxprogram', header_name='ulcg.h')
    ctx_check(features='cxx cxxprogram', lib='testu01', uselib_store='TESTU01', use='GMP')

    # mylib (part of TestU01)
    ctx_check(features='cxx cxxprogram', header_name='num.h')
    ctx_check(features='cxx cxxprogram', lib='mylib', uselib_store='MYLIB')


    # build variants
    env = ctx.env.derive()
    env.detach()

    # release (default)
    ctx.env.append_unique('CXXFLAGS', ['-O2'])
    ctx.define('NDEBUG', 1)

    ctx.setenv('debug', env)
    ctx.env.append_unique('CXXFLAGS', ['-g'])
    ctx.define('DEBUG', 1)


def distclean(ctx):
    ctx.recurse('latticetester')

    verfile = ctx.path.find_node('VERSION')
    if verfile:
        verfile.delete()
    from waflib import Scripting
    Scripting.distclean(ctx)


def build(ctx):

    if ctx.variant:
        print("Building variant `%s'" % (ctx.variant,))
    # if ctx.variant:
    #     print("Building variant `%s'" % (ctx.variant,))

    ctx.recurse('src')

    # ctx.recurse('analysis')
    ctx.recurse('progs')    

    # ctx.recurse('test')
    if ctx.env.BUILD_DOCS:
         ctx.recurse('doc')


# build variants

from waflib.Build import BuildContext
class debug(BuildContext):
    cmd = 'debug'
    variant = 'debug'
