# Python modules
from __future__ import division
import os
import shutil
import sys

"""A utility script that copies the Python import library (.lib) file into
a predictable location so that MSVC can find it. 

MSVC runs this script as a pre-build step. It's useful because it means the
user doesn't have to manually enter the the Python lib path as a build option.

Copying the file is a somewhat crude solution but the best we can currently
come up with. Although it's not too hard to pass programmatically-generated
options to the compile step of a Visual Studio build, it's AFAICT impossible
to do so with the linker.

See temp_python_lib/readme.txt for more info
"""

# Python doesn't provide a standard method for finding its import library.
# Furthermore, the name changes based on the Python version. However, both
# the path and name are predictable. The python.org version of Python,
# EPD's Python and all of the docs I can find place the .lib file in a
# subdirectory called libs below the directory that contains the Python 
# executable.
# The .lib file name is pythonNN.lib where NN are the major and minor 
# Python version numbers, e.g. python25.lib


# sys.executable is the path to the Python executable. It will be something 
# like "c:\Python25\python.exe"
path = sys.executable

# I trim the latter part
path, _ = os.path.split(path)

path = os.path.join(path, "libs")

major, minor = sys.version_info[:2]

lib_filename = "python%d%d.lib" % (major, minor)

# Copy the .lib file to ./temp_python_lib
source_filename = os.path.join(path, lib_filename)

destination_filename = os.path.join(os.getcwd(), "temp_python_lib", lib_filename)

shutil.copyfile(source_filename, destination_filename)
