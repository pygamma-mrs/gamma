from __future__ import division

import distutils.sysconfig

"""A utility script that writes Python's include directory into a compiler
/I switch. That is then written to a file and that file is passed by MSVC
to cl.exe via the @ option.

MSVC runs this script as a pre-build step. It's useful because it means the
user doesn't have to manually enter the the Python include path as a build 
option.
"""

path = distutils.sysconfig.get_python_inc()

path = '/I "%s"' % path 

open("python_include_path.rsp", "w").write(path)
