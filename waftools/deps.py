#!/usr/bin/env python
# coding: utf-8

def add_deps_path(ctx, what, where):
    try:
        prefix = ctx.root.find_node(where)
        ctx.env.append_value('INCLUDES', prefix.find_dir('include').abspath())
        ctx.env.append_value('LIBPATH',  prefix.find_dir('lib').abspath())
    except:
        ctx.fatal("cannot locate %s installation under `%s'" % (what, where))

def shared_or_static(ctx, f):
    if ctx.options.link_static:
        def f_static(*a, **kw):
            if 'lib' in kw:
                kw['stlib'] = kw['lib']
                del kw['lib']
            if 'shlib' in kw:
                kw['lib'] = kw['shlib']
                del kw['shlib']
            return f(*a, **kw)
        return f_static
    else:
        def f_shared(*a, **kw):
            if 'shlib' in kw:
                kw['lib'] = kw.get('lib', []) + kw['shlib']
                del kw['shlib']
            return f(*a, **kw)
        return f_shared


