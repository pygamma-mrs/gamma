# This file is renamed to __init__.py by setup.py before it is copied to
# the pygamma package. The comments below apply to this file's role as
# an __init__.py.

# SWIG creates a Python wrapper with the same name as our package (pygamma)
# and as a result the correct way to import it would be like so:
#    import pygamma.pygamma

# This isn't difficult but the duplication is kinda pointless so here we
# move the entire pygamma namespace up into this level which allows users
# to simply `import pygamma`

from __future__ import division

import sys

if sys.version_info < (3,0):
    from pygamma import *
else:
    from pygamma.pygamma import *

# I also need to get the version number from pygamma_version. 
# 'import *' doesn't import anything beginning with a double underscore, 
# so we have to import the important __version__ string explicitly.
if sys.version_info < (3,0):
    from pygamma_version import __version__
else:
    from .pygamma_version import __version__

if sys.version_info < (3,0):
    del pygamma

