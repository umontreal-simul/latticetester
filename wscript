#! /usr/bin/env python
# encoding: utf-8

top = '.'
out = 'build'

from waflib import Utils

import imp
def waftool(name):
    return imp.load_module('waf_' + name, *imp.find_module(name, ['./waftools', './latticetester/waftools']))

version = waftool('version')
compiler = waftool('compiler')
deps = waftool('deps')

def options(ctx):
    ctx.load('compiler_c compiler_cxx gnu_dirs waf_unit_test')
    ctx.add_option('--link-static', action='store_true', help='statically link with dependencies')
    ctx.add_option('--build-docs', action='store_true', default=False, help='build documentation')
    ctx.add_option('--boost', action='store', help='prefix under which Boost is installed')
    ctx.add_option('--ntl', action='store', help='prefix under which NTL is installed')
    ctx.add_option('--gmp', action='store', help='prefix under which GMP is installed')


def configure(ctx):
    build_platform = Utils.unversioned_sys_platform()
    ctx.msg("Build platform", build_platform)

    ctx.load('compiler_c compiler_cxx gnu_dirs waf_unit_test')
    ctx.check(features='cxx', cxxflags='-std=c++14')
    ctx.env.append_unique('CXXFLAGS', ['-std=c++14', '-Wall'])
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

    # options
    if ctx.options.boost:
        deps.add_deps_path(ctx, 'boost', ctx.options.boost)
    if ctx.options.ntl:
        deps.add_deps_path(ctx, 'NTL', ctx.options.ntl)
    if ctx.options.gmp:
        deps.add_deps_path(ctx, 'GMP', ctx.options.gmp)

    ctx_check = deps.shared_or_static(ctx, ctx.check)


    # Boost
    ctx_check(features='cxx cxxprogram',
            header_name='boost/config.hpp')

    # NTL
    # if ctx.options.ntltypes:  # ntlttypes are now mandatory   
    ctx_check(features='cxx cxxprogram', header_name='NTL/vector.h')
    ctx_check(features='cxx cxxprogram', lib='ntl', uselib_store='NTL')

    # GMP
    ctx_check(features='cxx cxxprogram', header_name='gmp.h')
    ctx_check(features='cxx cxxprogram', lib='gmp', uselib_store='GMP')

    # Doxygen
    if ctx.options.build_docs:
        ctx.env.BUILD_DOCS = True
        if not ctx.find_program('doxygen', var='DOXYGEN', mandatory=False):
            ctx.fatal('Doxygen is required for building documentation.\n' +
                      'Get it from http://www.stack.nl/~dimitri/doxygen/')


    ctx.version_file('latticetester')
    version_tag = ctx.set_version('latticetester')
    ctx.define('LATTICETESTER_VERSION', version_tag)
    ctx.msg("Setting LatticeTester version", version_tag)

    # definitions of DEBUG and NDEBUG flags: not used in this version
    # if not hasattr(ctx.options, 'nested') or not ctx.options.nested:
    #     # build variants
    #     env = ctx.env.derive()
    #     env.detach()

    #     # release (default)
    #     ctx.env.append_unique('CXXFLAGS', ['-O2'])
    #     ctx.define('NDEBUG', 1)

    #     ctx.setenv('debug', env)
    #     ctx.env.append_unique('CXXFLAGS', ['-g'])
    #     ctx.define('DEBUG', 1)


def distclean(ctx):

    verfile = ctx.path.find_node('VERSION')
    if verfile:
        verfile.delete()
    from waflib import Scripting
    Scripting.distclean(ctx)


def build(ctx):

    if ctx.variant:
        print("Building variant `%s'" % (ctx.variant,))

    ctx.recurse('src')

    if not hasattr(ctx.options, 'nested') or not ctx.options.nested:
        ctx.recurse('progs')
        if ctx.env.BUILD_DOCS:
            ctx.recurse('doc') 

    ctx.recurse('data')    

    # ctx.recurse('analysis') # probably broken, not used in the current version


# # build variants not used in this version

# from waflib.Build import BuildContext
# class debug(BuildContext):
#     cmd = 'debug'
#     variant = 'debug'
