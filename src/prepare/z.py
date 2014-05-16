from __future__ import division

import contextlib
import cPickle as pickle
import logging
import marshal
import os
import sys
import string
import traceback

def format_exception(e):
    exception_list = traceback.format_stack()
    exception_list = exception_list[:-2]
    exception_list.extend(traceback.format_tb(sys.exc_info()[2]))
    exception_list.extend(traceback.format_exception_only(sys.exc_info()[0], sys.exc_info()[1]))

    exception_str = "Traceback (most recent call last):\n"
    exception_str += "".join(exception_list)
    # Removing the last \n
    exception_str = exception_str[:-1]
    return exception_str

def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

@contextlib.contextmanager
def info(doc):
    logging.info("Started " + doc)
    yield
    logging.info("Finished " + doc)

def error(info):
    try:
        raise RuntimeError(info)
    except RuntimeError as e:
        logging.exception(format_exception(e))

    ## maybe exit?
def info(msg):
    logging.info(msg)

def warning(msg):
    logging.warning(msg)

def assertion(condition, msg):
    if not condition:
        return None
    error(msg)

def set_similarity(s, t):
    lst = len(s.union(t))
    if lst == 0:
        return 0
    return 1 - len(s.symmetric_difference(t)) / lst

def save(version):
    def dec(fn):
        def ret(*args):
            foldername = 'zsave'
            filename = fn.__name__
            extra = ':'.join([str(a) for a in args])
            trantab = string.maketrans('/', '|')
            extra = extra.translate(trantab)
            fullfile = foldername + '/' + filename + ':' + str(version) + ':' + extra
            lastfile = foldername + '/' + filename + ':' + str(version - 1) + ':' + extra
            if not os.path.exists(foldername):
                os.makedirs(foldername)

            if os.path.exists(lastfile):
                os.remove(lastfile)

            try:
                if os.path.exists(fullfile):
                    fd = open(fullfile, "rb")
                    return marshal.load(fd)
            except Exception as e:
                logging.exception(format_exception(e))

            v = fn(*args)
            fd = open(fullfile, "wb")
            marshal.dump(v, fd)
            return v
        return ret
    return dec

def pickledump(obj, filename):
    out = file(filename, "w+")
    pickle.dump(obj, out, pickle.HIGHEST_PROTOCOL)

def pickleload(filename):
    fd = file(filename, "r")
    return pickle.load(fd)
