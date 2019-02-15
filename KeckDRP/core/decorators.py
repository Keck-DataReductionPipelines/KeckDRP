from .. import conf


def repeat_processing(f):
    def g(*args, **kwargs):
        p = args[0]
        frame = args[1]
        p.set_frame(frame)
        p.read_proctab()
        if p.in_proctab() and conf.OVERWRITE is False:
            p.log.warning("Already processed")
            return None
        else:
            return f(*args, **kwargs)
    g.__name__ = f.__name__
    g.__doc__ = f.__doc__
    return g


