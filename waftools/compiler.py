#!/usr/bin/env python
# coding: utf-8

from waflib import Utils, Context, Errors

def check_compiler_flag(ctx, flag, error_markers=['error:', 'warning:'], mandatory=True):
    ctx.start_msg("checking for compiler flag %s" % flag)
    cmd = ctx.env.CXX + ['-xc++', '-E', flag, '-']
    result = True
    try:
        out, err = ctx.cmd_and_log(cmd, output=Context.BOTH, stdin=Utils.subprocess.PIPE)
        if any([m in err for m in error_markers]):
            result = False
    except Errors.WafError:
        result = False
    ctx.end_msg(result and "yes" or "no")
    if not result and mandatory:
        ctx.fatal("%s doesn't support flag %s" % (ctx.env.get_flat('CXX'), str(flag)))
    return result

def add_cxx_flag_if_supported(ctx, flag):
    #if ctx.check(features='cxx', cxxflags=flag, mandatory=False):
    if check_compiler_flag(ctx, flag, mandatory=False):
        ctx.env.append_value('CXXFLAGS', flag)

