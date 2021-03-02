from __future__ import division

import shutil
import sys
import os
import distutils.util

"""A utility script that handles the post-build events for the PyGAMMA
library. It performs several tasks --

1) Copies the PyGAMMA binary (_pygamma.pyd) and the SWIG shell (pygamma.py)
to gamma/pygamma/dist_staging for use by setup.py.
2) Writes pygamma_version.py to dist_staging
3) Cleans up (deletes) the pythonxx.lib copy created by clean_python_lib.py

MSVC runs this script as a post-build step.
"""

# Ensure destination dir exists.
path = "../../../pygamma/dist_staging"
if not os.path.exists(path):
    os.mkdir(path)

for filename in ('_pygamma.pyd', 'pygamma.py'):
    source = "../../../i686-pc-msvc/" + filename
    destination = "../../../pygamma/dist_staging/" + filename
    print( "Copying %s to %s..." % (source, destination))
    shutil.copyfile(source, destination)

version = open("../../../VERSION").read().strip()

version = "__version__ = '%s'" % version

open("../../../pygamma/dist_staging/pygamma_version.py", "w").write(version)

# Clean up the pythonxx.lib copy created by clean_python_lib.py.
major, minor = sys.version_info[:2]

lib_filename = "python%d%d.lib" % (major, minor)

lib_filename = os.path.join(os.getcwd(), "temp_python_lib", lib_filename)

if os.path.exists(lib_filename):
    os.remove(lib_filename)
