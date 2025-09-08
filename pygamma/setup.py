# Python modules
from __future__ import division
from __future__ import print_function

import sys
import os
import shutil

# 3rd party imports
# setuptools is distributed with Python.org Python, but it's still not part of the standard
# library so it might not be present (especially on older Pythons). We can install PyGamma
# without it. It's only required if you want to build a wheel.
try:
    import setuptools as distutools
except ImportError:
    import distutils.core as distutools


# Check Python version ASAP
major, minor = sys.version_info[:2]

if not (major,minor) in [(2,7), (3,9), (3,8), (3,7), (3,6)]:
    msg = "Sorry, PyGAMMA requires Python 2.7 or 3.6-3.9 and this is Python %d.%d." % (major, minor)
    print( msg )
    sys.exit(-1)

# The VERSION file will be in the current dir if this is being run from an
# expanded tarball/ZIP file (i.e. a normal end user distro) or in the parent
# dir if this is being run from an SVN tree. The latter is not just a
# developer-only thing. End users will need to use it if they're compiling
# their own binaries.
version_path = "VERSION"
if not os.path.exists(version_path):
    version_path = os.path.join("..", version_path)

version = open(version_path, "rb").read().strip()
version = version.decode('utf-8')

NAME = "pygamma"
DESCRIPTION = "A Python wrapper for the GAMMA C++ Library."
LONG_DESCRIPTION = """
GAMMA is a C++ library for the simulation of magnetic resonance (NMR, MRS) experiments. It
provides a simple and intuitive means to construct simulation programs to suit researchers'
individual needs. GAMMA is an acronym for a General Approach to Magnetic resonance
Mathematical Analysis.

PyGamma is a Python wrapper around GAMMA that makes almost all of GAMMA's API available via
Python.

Both GAMMA and PyGamma work on OS X, Linux, and Windows.

As of version 4.3.4 only 64bit libraries are provided. This will also be the last release 
that includes libraries for Python 2.
"""
# Note that Python's distutils writes a PKG-INFO file that replaces the author metadata with
# the maintainer metadata. As a result, it's impossible (AFAICT) to get correct author metadata
# to appear on PyPI.
# https://bugs.python.org/issue16108
AUTHOR = "Dr. Scott A. Smith and Dr. Tilo Levante"
AUTHOR_EMAIL = "no_known_address@example.com"
MAINTAINER = "Brian Soher"
MAINTAINER_EMAIL = "bsoher@briansoher.com"
URL = "https://pygamma-mrs.github.io/gamma.io/"
# http://pypi.python.org/pypi?:action=list_classifiers
CLASSIFIERS = ['Development Status :: 5 - Production/Stable',
               'Intended Audience :: Science/Research',
               "License :: OSI Approved :: BSD License",
               "Operating System :: MacOS :: MacOS X",
               "Operating System :: POSIX :: Linux",
               "Operating System :: POSIX :: SunOS/Solaris",
               "Operating System :: Microsoft :: Windows",
               "Operating System :: Unix",
               ]
LICENSE = "http://creativecommons.org/licenses/BSD/"
KEYWORDS = "GAMMA pygamma MRI MRS MR magnetic resonance spectroscopy"
PLATFORMS = 'Linux, OS X, Windows, Solaris, POSIX'


def list_non_py_files(path):
    """Given a path, returns a list of the non-Python files found there. The files
    returned are just names; they're not fully qualified.

    Directories, files that start with '.', .py and .pyc files are excluded.

    This function doesn't travel into subdirectories.
    """
    filenames = []

    for filename in os.listdir(path):
        if os.path.isfile(os.path.join(path, filename)) and \
           not filename.startswith(".") and \
           not filename.endswith(".pyc") and \
           not filename.endswith(".py"):
            filenames.append(filename)

    return filenames


# dist_staging is where files reside in preparation for creating an installable package and/or a
# distributable wheel file. The build process (which is Makefile execution under *nix, the
# post_build.py step under Windows) places files here, as does some code below.
package_dir = {"pygamma": 'dist_staging'}

# package_data tells setup about files other than .py that I want to be in the package. I grab
# everything in the staging directory that's not a .py file and stuff it into package_data.
package_data = {"pygamma": list_non_py_files('dist_staging')}

# The package needs an __init__.py file to make pygamma a proper package.
# Many __init__.py files are blank; this one is not. In order to manage the
# source code somewhat sanely, I put the code in a file called init.py and
# copy it to dist_staging/__init__.py during this setup. That way the file is
# under source code management but no one will think that it's meant to make
# the setup.py directory a package.
shutil.copy("init.py", os.path.join('dist_staging', "__init__.py"))

# The current package for building wheels (v0.26.0) has some shortcomings when it comes to
# building wheels out of binaries that aren't compiled by the setup process itself. That's
# documented in wheel ticket 128:
# https://bitbucket.org/pypa/wheel/issues/128/allow-forcing-pure-platform-specific-wheel
# If setup is being invoked to build a wheel, we subclass the bdist_wheel command to get
# the behavior we want.
try:
    from wheel.bdist_wheel import bdist_wheel

    class my_bdist_wheel(bdist_wheel):
        def finalize_options(self):
            """Subclass finalize_options() to fix the value of root_is_pure"""
            bdist_wheel.finalize_options(self)
            # This *is not* a pure Python package
            self.root_is_pure = False

        def get_tag(self):
            """Subclass get_tag() to fix the platform tag"""
            tag = bdist_wheel.get_tag(self)

            if sys.platform == 'darwin':
                # Under OS X, bdist_wheel.get_tag() returns a platform tag that contains the min
                # OS X version with which this Python is compatible (e.g. 'macosx_10_5_x86_64' for
                # Anaconda Python). That's appropriate for C/C++ extensions built via the Python
                # build process, but not for PyGamma which is built via a custom build process.
                # PyGamma's min OS X version is determined by a value in the Makefile, so we give
                # the make process responsiblity for writing a file that contains an appropriate
                # platform tag (e.g. 'macosx_10_9_x86_64').
                platform_tag = open('dist_staging/wheel_platform_tag.txt', 'rb').read().strip()
                if sys.version_info[0] > 2:
                    # python 3 reads text file as array of bytes, but we need str, so we decode
                    platform_tag = platform_tag.decode()
                tag = (tag[0], tag[1], platform_tag)
            # else:
                # Do nothing. On Linux and Windows, bdist_wheel.get_tag() generates an
                # appropriate platform tag.

            return tag

    cmdclass = {'bdist_wheel': my_bdist_wheel}
except ImportError:
    # This is an error we can ignore -- only people who want to create wheels need the
    # wheel package installed. Others can still use this setup.py to install homegrown
    # versions of PyGamma.
    cmdclass = {}


distutools.setup(name=NAME,
                 version=version,
                 package_dir=package_dir,
                 packages=["pygamma"],
                 package_data=package_data,
                 cmdclass=cmdclass,
                 url=URL,
                 author=AUTHOR,
                 author_email=AUTHOR_EMAIL,
                 maintainer=MAINTAINER,
                 maintainer_email=MAINTAINER_EMAIL,
                 classifiers=CLASSIFIERS,
                 license=LICENSE,
                 keywords=KEYWORDS,
                 description=DESCRIPTION,
                 long_description=LONG_DESCRIPTION,
                 platforms=PLATFORMS,
                 )
