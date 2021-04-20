import re
import os
import os.path

def get_bsmdir():
    pythonpath = os.environ['PYTHONPATH']
    bsmdir = re.sub('(^[^:]+/bsm)/src:.*$', '\\1', pythonpath)
    return(bsmdir)
