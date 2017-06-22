#!/usr/bin/env python
# coding: utf-8

VERSION = '1.0'

from waflib.Configure import conf

@conf
def version_file(ctx):
    """Update the VERSION file by reading Git tags"""
    try:
        version = ctx.cmd_and_log("git describe --tags --match 'v[0-9]*'", cwd=ctx.path.abspath()).strip()[1:]

        # extract revision
        pos = version.rfind('-g')
        if pos > 0:
            version = version[:pos]
        pos = version.rfind('-')
        if pos > 0:
            r = version[pos:]
            version = version[:pos]
        else:
            r = ''

        # branch
        branch = ctx.cmd_and_log("git rev-parse --abbrev-ref HEAD", cwd=ctx.path.abspath()).strip()
        if branch == 'HEAD':
            version += '-unknown'
        elif branch != 'master':
            version += '-' + branch

        version += r

        ctx.path.make_node('VERSION').write(version)
    except:
        pass

@conf
def set_version(ctx):
    """Returns VERSION from environment or by reading the VERSION file"""
    if 'VERSION' in ctx.env:
        version = ctx.env['VERSION']
    else:
        verfile = ctx.path.find_node('VERSION')
        if verfile:
            version = verfile.read().strip()
        else:
            version = 'unknown'
    global VERSION
    VERSION = version
    return version

