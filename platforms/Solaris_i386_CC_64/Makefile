#########################################################################
#
#          Solaris Makefile for the GAMMA and PyGAMMA Libraries
#
# This builds GAMMA and PyGAMMA. It assumes that one is using the Solaris
# Studio C++ compiler. To build PyGAMMA and to run the tests one must have
# Python available.
#
# ---------------------------------------------------------------------
#
#                             Quick Info
#
#  - Type "make" or "make all" to build GAMMA static and shared libraries.
#  - Type "make install" to install the libraries, gamma script, and
#    the necessary include files. May require sudo or root privileges.
#  - Type "make pysg" to compile the PyGAMMA shared library.
#  - Type "make pysgdist" to compile the PyGAMMA shared library and copy
#    the output files to the pygamma/dist_staging directory to prepare for
#    installation or building a distributable wheel via setup.py.
#  - Type "make test" to compile and run the GAMMA test program(s).
#  - Type "make clean" to remove all object files.
#  - Type "make distclean" to remove all object files libraries, exes.
#  - The sources are assumed to be in SRCDIR
#  - The object files are produced in BINDIR
#  - The lib files are produced in BINDIR
#  - Any executable files are produced in BINDIR.
#  - The important static library file is libgamma.a
#  - The important dynamic library file is libgamma.so
#  - The important PyGAMMA shared library file is _pygamma.so
#  - Produces the following: a.) libgamma.a - The GAMMA static library
#                            b.) libgamma.so - The GAMMA shared library
#                            c.) _pygamma.so - The PyGAMMA shared library
#
# ---------------------------------------------------------------------
#
#                                 Notes
#
# This is NOT a standalone makefile. It includes an additional file
# ../../make/MakeMods which is part of standard GAMMA. This latter file
# defines which files are part of the GAMMA source codes. Changes in
# MakeMods will evoke changes in the build herein.
#
# make all: These both run the default target all (first target). It
#           will compile all source codes in GAMMA as defined in file
#           MakeMods into object files. It will then take the object
#           files and convert them into a static library named
#           libgamma.a. Similarly, other object files object files will
#           be compiled and then be combined into a shared library of
#           name libgamma.so. All libraries and object files will be
#           made in the directory specified by BINDIR. The shared libs
#           will be relocated with added soft links when installed.
#
# make install: This will convert the file gamma.in into the working
#               script gamma and move it into a system file so that it
#               is accessible by all users. Additionally it will move
#               the shared library for GAMMA into a system directory
#		place appropriate soft links so that it is generally
#               workable. Lastly, if PyGAMMA has been built, its
#               shared library _pygamma.so will be moved into a Python
#               directory.
#
# make clean: This will delete all object files from any earlier build.
#             It is important to run make clean since the object files
#             can use lots of disk space as well as interfere with
#             other GAMMA builds (not part of this Makefile)
#
# make test: This will build and run a set of programs that test the
#            the modules provided with GAMMA.
#
#########################################################################

#########################################################################
#
#                Semi-Standard Gnuish Makefile Variables
#
# CXX         = C++ compiler (used in default compiler rules)
# CXXFLAGS    = C++ compilation flags (used in default C++ builds)
#
# The default rule in make that uses the variables CXX and
# CXXFLAGS to generate object files from .cc files is as follows:
#
# %.o : %.cc
#        $(CXX) -c $(CXXFLAGS) $< -o $@
#
#########################################################################

SHELL         = /bin/bash
PYTHON        = python
SWIG          = swig
# If you have SWIG >= 2.0.4, you can add -builtin to the SWIG_OPTIONS for
# a (very) modest performance gain. See here for details:
# http://scion.duhs.duke.edu/vespa/gamma/wiki/SwigBuiltinTypesBenefit
SWIG_OPTIONS  = -c++ -python

GAMMA_ARCH    = i386-pc-solaris2.10-CC64
INSTALLDIR    = /nmrnet/gamma
BININSTALLDIR = ${INSTALLDIR}/bin
LIBINSTALLDIRV= ${INSTALLDIR}/lib/gamma-${MAJORVERSION}.${MINORVERSION}.${RELEASEVERSION}
LIBINSTALLDIR = ${INSTALLDIR}/lib/gamma-${MAJORVERSION}.${MINORVERSION}.${RELEASEVERSION}/${GAMMA_ARCH}
HDRINSTALLDIR = ${INSTALLDIR}/include/gamma-${MAJORVERSION}.${MINORVERSION}.${RELEASEVERSION}

ROOTDIR = ../..
SRCDIR  = ../../src

# comment out these next two lines if not using blas/sunperflib
# and uncomment both of the following two lines
BLASINCL_FLAG   = -I/opt/local/include -D_USING_BLAS_ -D_USING_SUNPERFLIB_
LAPACKINCL_FLAG =
#BLASINCL_FLAG   =
#LAPACKINCL_FLAG =

SRCDIRINCL = -I$(SRCDIR)

BINDIR = ./../../${GAMMA_ARCH}

DEBUG    = -g
OPTIMIZE = -fast -fsimple=1 -xtarget=pentium4 -m64
WARNINGS =
PLATFORM_FLAGS =

CC       =    /opt/SUNWspro/bin/cc
CFLAGS   =    $(WARNINGS) $(PLATFORM_FLAGS) $(SRCDIRINCL) $(DEBUG) \
              $(BLASINCL_FLAG) $(LAPACKINCL_FLAG) $(OPTIMIZE)

CXX      =    /opt/SUNWspro/bin/CC
CXXFLAGS =    $(WARNINGS) $(PLATFORM_FLAGS) $(SRCDIRINCL) $(DEBUG) \
              $(BLASINCL_FLAG) $(LAPACKINCL_FLAG) $(OPTIMIZE)

LIBCXXFLAGS = $(WARNINGS) $(PLATFORM_FLAGS) $(SRCDIRINCL) $(DEBUG) \
              $(BLASINCL_FLAG) $(LAPACKINCL_FLAG) $(OPTIMIZE)

SOCXXFLAGS =  -KPIC \
              $(WARNINGS) $(PLATFORM_FLAGS) $(SRCDIRINCL) $(DEBUG) \
              $(BLASINCL_FLAG) $(LAPACKINCL_FLAG) $(OPTIMIZE)

PYSGCXXFLAGS = -KPIC \
              $(WARNINGS) $(PLATFORM_FLAGS) $(SRCDIRINCL) $(DEBUG) \
              $(BLASINCL_FLAG) $(LAPACKINCL_FLAG) $(OPTIMIZE)


LDFLAGS  = -m64

AR       = ar
ARFLAGS  = r
RANLIB   = ranlib

LDCONFIG = /sbin/ldconfig

# comment out the next three lines if not using blas/sunperflib
# and uncomment the last line.
LAPACK_LIBS = -L/opt/solstudio12.2/lib/amd64 -R/opt/solstudio12.2/lib/amd64 -lsunperf -lsunmath
BLAS_LIBS  = -L/opt/local/lib/amd64 -lcblas
LAB_LIBS = $(BLAS_LIBS) $(LAPACK_LIBS)
#LAB_LIBS =

#########################################################################
#
#	           Names Of Various Targets To Build
#
# These names are associated with the library we plan to build.
#
# LIB_NAME  = Name of library to be generated
# SO_NAME   = Name of shared lib to generate
# PYS_NAME  = Name of SWIG shared lib to generate
#
# Note: The shared library name should be libgamma.#.#.# where the 3
#       numbers are MAJORVERSION, MINORVERSION, and RELEASEVERSION
#       respectively. For GAMMA we
#       are currently on major version 4 and minor version 2 so that it
#       will have an soname of libgamma.so.4 whereas the library name
#       produced in the build will be libgamma.so.4.2. When the library
#       libgamma.so.4.2 is put into a standard directory (e.g. into
#       /usr/lib .... done with make install) then the program ldconfig
#       will make the symbolic link between this file and the soname.
#       That is, libgamma.so.4 -> libgamma.so.4.2. A second link is
#       set manually, libgamma.so -> libgamma.so.4 so that ld will find
#       the right GAMMA library (with soname) to use at link time. Then
#       if GAMMA is upgraded but not to the point that any previous
#       programs are affected, the library may be rebuilt with the same
#       soname (that which all previously made executable try to link
#       to) but change the MINORVERSION for the build. If the GAMMA changes
#       will break older programs, then adjust the major version and
#       the soname (e.g. libgamma.so.5.0 and libgamma.so.5) for the
#       rebuild. Then (after the install to /usr/local/lib) then switch the
#       symbolic link to point to the new soname library:
#       libgamma.so -> libgamma.so.5.
#
#########################################################################

LIB_NAME    = libgamma.a

SO_NAME    = libgamma.so
SO_NAMEX   = ${SO_NAME}.${MAJORVERSION}
SO_NAMEXX  = ${SO_NAMEX}.${MINORVERSION}
SO_NAMEXXX = ${SO_NAMEXX}.${RELEASEVERSION}

PYS_NAME    = _pygamma.so
PYS_NAMEX   = ${PYS_NAME}.${MAJORVERSION}
PYS_NAMEXX  = ${PYS_NAMEX}.${MINORVERSION}
PYS_NAMEXXX = ${PYS_NAMEXX}.${RELEASEVERSION}

#########################################################################
#
#                       Gamma Related Variables
#
#
# These variables are used to specify the GAMMA sources. Exactly what
# sources are associated with any one GAMMA module is contained in the
# external makefile MakeMods. From these we can construct the names of
# module source files and object files. These are directory specific.
#
# MakeMods    = File listing all current GAMMA source files.
# SRCDIR      = Base directory of GAMMA source codes.
# BINDIR      = Base directory of binaries.
# LIB_SRCS    = All code files to be put into the GAMMA static library.
# LIB_OBJS    = All object files to be put in the GAMMA static library.
# SO_SRCS     = All code files to be put into the GAMMA shared library.
# SO_OBJS     = All object files to be put in the GAMMA shared library.
#
# The variable ALLCCFILES should contain all GAMMA source .cc files.
# The variable ALLOFILES should contain all GAMMA object .o files.
# The source file names will have the source sub-directory as a prefix.
#     but not the full path relative to this directory.
#
# The included makefile supplies NO targets or rules of any kind.
#
#########################################################################


include ../../make/MakeMods


#########################################################################
#
#                              Main Target
#
# This target produces the dynamic and static library.
#
# LIB_NAME    = Name of library to be generated/targeted.  (libgamma.a)
# SO_NAME     = Name of library to be generated/targeted.  (libgamma.so)
#
#########################################################################


all: so lib
lib: $(BINDIR) ${BINDIR}/${LIB_NAME}
so:  $(BINDIR) $(BINDIR)/$(SO_NAME)
pysg: $(BINDIR) $(BINDIR)/$(PYS_NAME)

${LIB_NAME}: ${BINDIR}/${LIB_NAME}
${SO_NAME}:  ${BINDIR}/${SO_NAME}
${PYS_NAME}: ${BINDIR}/${PYS_NAME}

big: so lib pysg test


#########################################################################

LIB_SRCS = $(ALLCCFILES)
LIB_OBJS = $(addsuffix .o, $(addprefix $(BINDIR)/, $(ALLFILES)))

SO_SRCS = $(ALLCCFILES)
SO_OBJS = $(addsuffix .soo, $(addprefix $(BINDIR)/, $(ALLFILES)))

# Python Swig
PYSG_SRCS = $(ALLCCFILES)
PYSG_SRCS =
PYSG_OBJS = $(addsuffix .pso, $(addprefix $(BINDIR)/, $(ALLFILES)))

GAM_HDRS = $(ALLHFILES)
INS_HDRS = $(addprefix $(HDRINSTALLDIR)/, $(ALLHFILES))

PYSG_WRAP_SRC = $(addprefix $(BINDIR)/, pygamma_wrap.cxx)
PYSG_WRAP_SRC1 = $(addprefix $(BINDIR)/, pygamma.py)

PYSG_WRAP_SRCS = $(PYSG_WRAP_SRC) $(PYSG_WRAP_SRC1)
PYSG_WRAP_OBJS = $(addprefix $(BINDIR)/, pygamma_wrap.o)

# If add a new .i file, it needs to be added to gamma/trunk/make/Makemods
# You'll see examples there of how to do this.

OTHERDOTIFILES = $(addprefix $(SRCDIR)/, $(ALLIFILES))

$(PYSG_WRAP_SRC) : ${SRCDIR}/pygamma.i $(OTHERDOTIFILES)
	@echo
	${SWIG} ${SWIG_OPTIONS} -outdir $(BINDIR) -o $@ \
         -I../../src $<

#########################################################################
#
#                Build GAMMA Static Library Pattern Rules
#
# These rules specify how we build object file (*.o) from C++ code
# files (*.cc). Gnu make has a built in implicit pattern rule to do
# this task:
#
# %.o : %.cc
#        $(CXX) -c $(CXXFLAGS) $< -o $@
#
# But that does not suffice for two reasons. First, GAMMA source files
# reside in module directories off of the source directory SRCDIR.
# Second, we want all the object files to be built in the binaries
# directory, BINDIR. Thus we must build ${BINDIR}/*.o from a code file
# ${SRCDIR}/${MODDIR}/*.cc. The main GAMMA  makefile handles this by
# recursing through each module directory and letting it decide how to
# build itself. Here we will define a rule for each module directory.
#
#########################################################################

BASICS    = Basics
MATRIX    = Matrix
HSPACE    = HSLib
LEVEL1    = Level1
LSPACE    = LSLib
LEVEL2    = Level2
BLOCH     = Bloch
GAMIO     = GamIO
BWR       = BWRRelax
PULSES    = Pulses
INTRANK2  = IntRank2
FLOQUET   = Floquet
GRADIENTS = Gradients
MULTISYS  = MultiSys
ESRLIB    = ESRLib
DEPREC    = Deprecated
TESTING   = Testing
MODULEDIRS = ${ALLMODDIRS}

#########################################################################
#
#                Build GAMMA Static Library Pattern Rules
#
# These rules specify how we build object file (*.o) from C++ code
# files (*.cc). Gnu make has a built in implicit pattern rule to do
# this task:
#
# %.o : %.cc
#        $(CXX) -c $(CXXFLAGS) $< -o $@
#
# But that does not suffice for two reasons. First, GAMMA source files
# reside in module directories off of the source directory SRCDIR.
# Second, we want all the object files to be built in the binaries
# directory, BINDIR. Thus we must build ${BINDIR}/*.o from a code file
# ${SRCDIR}/${MODDIR}/*.cc. The main GAMMA  makefile handles this by
# recursing through each module directory and letting it decide how to
# build itself. Here we will define a rule for each module directory.
#
#########################################################################

${BINDIR}/%.o : ${SRCDIR}/${BASICS}/%.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${MATRIX}/%.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${HSPACE}/%.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${LEVEL1}/%.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${LSPACE}/%.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${LEVEL2}/%.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${BLOCH}/%.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${GAMIO}/%.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${BWR}/%.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${PULSES}/%.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${INTRANK2}/%.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${FLOQUET}/%.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${GRADIENTS}/%.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${MULTISYS}/%.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${ESRLIB}/%.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${DEPREC}/%.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${TESTING}/%.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

#########################################################################
#
#                Build GAMMA Shared Library Pattern Rules
#
# These rules specify how we build object file (*.soo) from C++ code
# files (*.cc). Gnu make has a built in implicit pattern rule to do
# this task:
#
# %.o : %.cc
#        $(CXX) -c $(CXXFLAGS) $< -o $@
#
# But that does not suffice for a few reasons. First, GAMMA sources
# reside in module directories off of the source directory SRCDIR.
# Second, we want all the object files to be built in the binaries
# directory, BINDIR. Thus we must build ${BINDIR}/*.o from a code file
# ${SRCDIR}/${MODDIR}/*.cc. Third, we want to clearly distinguish the
# object files for the shared library from those that are static. For
# that we will use the suffix .soo rather than .o. Here we will define
# define a rule for each module directory.
#
#########################################################################

${BINDIR}/%.soo : ${SRCDIR}/${BASICS}/%.cc
	$(CXX) -c $(SOCXXFLAGS) $< -o $@

${BINDIR}/%.soo : ${SRCDIR}/${MATRIX}/%.cc
	$(CXX) -c $(SOCXXFLAGS) $< -o $@

${BINDIR}/%.soo : ${SRCDIR}/${HSPACE}/%.cc
	$(CXX) -c $(SOCXXFLAGS) $< -o $@

${BINDIR}/%.soo : ${SRCDIR}/${LEVEL1}/%.cc
	$(CXX) -c $(SOCXXFLAGS) $< -o $@

${BINDIR}/%.soo : ${SRCDIR}/${LSPACE}/%.cc
	$(CXX) -c $(SOCXXFLAGS) $< -o $@

${BINDIR}/%.soo : ${SRCDIR}/${LEVEL2}/%.cc
	$(CXX) -c $(SOCXXFLAGS) $< -o $@

${BINDIR}/%.soo : ${SRCDIR}/${BLOCH}/%.cc
	$(CXX) -c $(SOCXXFLAGS) $< -o $@

${BINDIR}/%.soo : ${SRCDIR}/${GAMIO}/%.cc
	$(CXX) -c $(SOCXXFLAGS) $< -o $@

${BINDIR}/%.soo : ${SRCDIR}/${BWR}/%.cc
	$(CXX) -c $(SOCXXFLAGS) $< -o $@

${BINDIR}/%.soo : ${SRCDIR}/${PULSES}/%.cc
	$(CXX) -c $(SOCXXFLAGS) $< -o $@

${BINDIR}/%.soo : ${SRCDIR}/${INTRANK2}/%.cc
	$(CXX) -c $(SOCXXFLAGS) $< -o $@

${BINDIR}/%.soo : ${SRCDIR}/${FLOQUET}/%.cc
	$(CXX) -c $(SOCXXFLAGS) $< -o $@

${BINDIR}/%.soo : ${SRCDIR}/${GRADIENTS}/%.cc
	$(CXX) -c $(SOCXXFLAGS) $< -o $@

${BINDIR}/%.soo : ${SRCDIR}/${MULTISYS}/%.cc
	$(CXX) -c $(SOCXXFLAGS) $< -o $@

${BINDIR}/%.soo : ${SRCDIR}/${ESRLIB}/%.cc
	$(CXX) -c $(SOCXXFLAGS) $< -o $@

${BINDIR}/%.soo : ${SRCDIR}/${DEPREC}/%.cc
	$(CXX) -c $(SOCXXFLAGS) $< -o $@

${BINDIR}/%.soo : ${SRCDIR}/${TESTING}/%.cc
	$(CXX) -c $(SOCXXFLAGS) $< -o $@


#########################################################################
#
#         Build _PGSGAMMA (SWIG) Shared Library Pattern Rules
#
########################################################################

${BINDIR}/%.pso : ${SRCDIR}/${BASICS}/%.cc
	$(CXX) -c $(PYSGCXXFLAGS) $< -o $@

${BINDIR}/%.pso : ${SRCDIR}/${MATRIX}/%.cc
	$(CXX) -c $(PYSGCXXFLAGS) $< -o $@

${BINDIR}/%.pso : ${SRCDIR}/${HSPACE}/%.cc
	$(CXX) -c $(PYSGCXXFLAGS) $< -o $@

${BINDIR}/%.pso : ${SRCDIR}/${LEVEL1}/%.cc
	$(CXX) -c $(PYSGCXXFLAGS) $< -o $@

${BINDIR}/%.pso : ${SRCDIR}/${LSPACE}/%.cc
	$(CXX) -c $(PYSGCXXFLAGS) $< -o $@

${BINDIR}/%.pso : ${SRCDIR}/${LEVEL2}/%.cc
	$(CXX) -c $(PYSGCXXFLAGS) $< -o $@

${BINDIR}/%.pso : ${SRCDIR}/${BLOCH}/%.cc
	$(CXX) -c $(PYSGCXXFLAGS) $< -o $@

${BINDIR}/%.pso : ${SRCDIR}/${GAMIO}/%.cc
	$(CXX) -c $(PYSGCXXFLAGS) $< -o $@

${BINDIR}/%.pso : ${SRCDIR}/${BWR}/%.cc
	$(CXX) -c $(PYSGCXXFLAGS) $< -o $@

${BINDIR}/%.pso : ${SRCDIR}/${PULSES}/%.cc
	$(CXX) -c $(PYSGCXXFLAGS) $< -o $@

${BINDIR}/%.pso : ${SRCDIR}/${INTRANK2}/%.cc
	$(CXX) -c $(PYSGCXXFLAGS) $< -o $@

${BINDIR}/%.pso : ${SRCDIR}/${FLOQUET}/%.cc
	$(CXX) -c $(PYSGCXXFLAGS) $< -o $@

${BINDIR}/%.pso : ${SRCDIR}/${GRADIENTS}/%.cc
	$(CXX) -c $(PYSGCXXFLAGS) $< -o $@

${BINDIR}/%.pso : ${SRCDIR}/${MULTISYS}/%.cc
	$(CXX) -c $(PYSGCXXFLAGS) $< -o $@

${BINDIR}/%.pso : ${SRCDIR}/${ESRLIB}/%.cc
	$(CXX) -c $(PYSGCXXFLAGS) $< -o $@

${BINDIR}/%.pso : ${SRCDIR}/${DEPREC}/%.cc
	$(CXX) -c $(PYSGCXXFLAGS) $< -o $@

${BINDIR}/%.pso : ${SRCDIR}/${TESTING}/%.cc
	$(CXX) -c $(PYSGCXXFLAGS) $< -o $@

$(PYSG_WRAP_OBJS) : $(PYSG_WRAP_SRC)
	$(CXX) -c $(PYSG_WRAP_SRC)  \
        $(PYSGCXXFLAGS) $(PYSGINCFLAGS) -o $(PYSG_WRAP_OBJS)

#########################################################################
#
#                      Copying GAMMA Headers Rules
#
# These rules specify how to update the GAMMA header files in the
# installation directory from those same header files in a source
# directory. If the headers are to be kept in the source directory then
# these targets should never be used.
#
#########################################################################

${HDRINSTALLDIR}/%.h : ${SRCDIR}/%.h
	cp -f $< $@

#########################################################################
#
#                       Library Related Variables
#
# These variables are used exclusively when building the static lib.
#
# LIB_SRCS    = All code files to be put into the library.
# LIB_OBJS    = All object files to be put into the library.
# LIB_NAME    = Name of library to be generated/targeted.  (libgamma.a)
# LIB_OBJS    = Name of object files used in library.
#
# These variables are used exclusively when building the shared library.
#
# SO_LDFLAGS = Flags used in .so build.                      (-shared)
# SO_LDLIBS  = Any libraries the shared lib. to be built depends upon.
# SO_SRCS    = All code files to be put into the shared library.
# SO_OBJS    = All object files to be put into the shared lib. Note
#              they are NOT the same as LIB_OBJ since we use special
#              compiler flags to generate them over the static library
#
# These variables are used generally.
#
# BINDIR      = Where we will ultimately place or binaries. This is
#               set to replace a CygWin build.
#
#########################################################################

SO_LDFLAGS = -G $(LDFLAGS)
SO_LDLIBS  = $(LAB_LIBS)

#########################################################################
#
#                      Build The GAMMA Static Library
#
# LIB_SRCS    = All code files to be put into the library.
# LIB_OBJS    = All object files to be put into the library.
# LIB_NAME    = Name of library to be generated/targeted.  (libgamma.a)
# LIB_OBJS    = Name of object files used in library.
#
# The static library build depends on the object files in LIB_OBJS. If
# or when these are up to date then the target uses ar to form the
# library from the object files and then ranlib to adjust the library.
#
# The LIB_OBJS variable is set from LIB_SRCS earlier in this file.
# These are apparently implicitly dependent upon the source codes in
# LIB_SRCS (i.e. the .cc file with the same base name as the .o file)
# & have an implicit pattern rule as to how the .cc files are compiled
# into .o files. This pattern rule, as specified in the documentation
# for GNU make, is
#
# %.o : %.cc
#        $(CXX) -c $(CXXFLAGS) $< -o $@
#
# This implicit rule can be replaced by defining a new pattern rule
# with the same target and prerequisites, but different commands.
# If there are no commands following the rule, the rule is cancelled.
# When the new rule is defined, the built-in one is replaced. The new
# rule's position in the sequence of implicit rules is determined by
# where you write the new rule and affects how they are chained.
#
#########################################################################

${BINDIR}/${LIB_NAME}: ${BINDIR} ${LIB_OBJS}
	@echo
	@echo "Making Static Library $(LIB_NAME) From Object Files"
	@echo
	${AR} ${ARFLAGS} ${BINDIR}/${LIB_NAME} $(LIB_OBJS)
	@echo
	@echo "Adjusting Static Library $(LIB_NAME) For Random Access"
	@echo
	${RANLIB} ${BINDIR}/${LIB_NAME}

#########################################################################
#
#                       Build The GAMMA Shared Library
#
# SO_NAME    = Name of shared library to be generated.   (libgamma.so)
# SO_OBJS    = Name of object files used in shared lib.
# SO_LDFLAGS = Flags used in shared library linking.         (-shared)
# SO_LDLIBS  = Libraries GAMMA shared library depends upon.
#
# The shared library build depends on the object files in SO_OBJS. If
# or when # these are up to date then the target uses g++ to form the
# library from the object files using appropriate linking commands.
#
# The SO_OBJS variable is set from SO_SRCS earlier in this file.
# These are NOT the same as LIB_OBJS used for the static library build
# because different flags are used when the source files are compiled.
# As with LIB_OBJS, the objects files are implicitly dependent upon the
# source codes in SO_SRCS (i.e. the .cc file with the same base name as
# the .soo file) & have an implicit pattern rule as to how the .cc
# files are compiled into .soo files. This pattern rule, as specified
# in the documentation for GNU make, is
#
# %.o : %.cc
#        $(CXX) -c $(CXXFLAGS) $< -o $@
#
# Since we switched from .o to .soo and use particular sub-directories
# we have made our own rules to handle the compilation of sources.
#
#########################################################################

${SO_NAME}: ${BINDIR}/${SO_NAME}

${BINDIR}/${SO_NAME}: ${BINDIR}/${SO_NAMEXXX}
	@echo "Making Links For Shared Library"
	@rm -f ${BINDIR}/${SO_NAME}
	@rm -f ${BINDIR}/${SO_NAMEX}
	@rm -f ${BINDIR}/${SO_NAMEXX}
	@ln -s ${SO_NAMEXXX} ${BINDIR}/${SO_NAMEXX}
	@ln -s ${SO_NAMEXXX} ${BINDIR}/${SO_NAMEX}
	@ln -s ${SO_NAMEXXX} ${BINDIR}/${SO_NAME}

$(BINDIR)/$(SO_NAMEXXX): $(SO_OBJS)
	@echo
	@echo "Making Shared Library $(SO_NAME) From Object Files"
	@echo
	@$(CXX)  ${SO_LDFLAGS} -h ${SO_NAMEX} -o ${BINDIR}/${SO_NAMEXXX} $(SO_OBJS) $(SO_LDLIBS)

########################################################################
#
#  Python Swig build shared library details
#
########################################################################

# Using python to find the default library and include path for "python"
# executable on the local machine.

PYTHONINCDIRS := ${shell ${PYTHON} -c "import distutils.sysconfig; print distutils.sysconfig.get_python_inc()"}
PYTHONLIBDIRS := ${shell ${PYTHON} -c "import distutils.sysconfig; print distutils.sysconfig.get_python_lib()"}

PYSGINCFLAGS = -I$(PYTHONINCDIRS)
PYSG_LDFLAGS = -L$(PYTHONLIBDIRS) -R$(PYTHONLIBDIRS)

PYSGAMMA_CHK = $(addprefix $(BINDIR)/, python_lib_incl.chk)

$(PYSGAMMA_CHK):
	@echo
	@echo Include directory: $(PYTHONINCDIRS)
	@echo Library directory: $(PYTHONLIBDIRS)
	@echo
	@echo Checking for python include and lib directories
	@echo
	@if (test -d ${PYTHONINCDIRS}); \
	then echo Good. The Python Includes Directory Exists. ${PYTHONINCDIRS}; \
	else echo Cannot Find The Python Include Directory! ${PYTHONINCDIRS}; \
	     echo Python Must Be Installed On The System For PyGAMMA To Work.; \
	     echo Sorry, Cannot Build PyGAMMA Until This Is Taken Care Of.; \
	     echo ; \
	     false; \
	fi
	@if (test -d ${PYTHONLIBDIRS}); \
	then echo Good. The Python Library Directory Exists. ${PYTHONLIBDIRS}; \
	else echo Cannot Find The Python Library Directory! ${PYTHONLIBDIRS}; \
		echo Python Must Be Installed On The System For PyGAMMA To Work.; \
		echo Sorry, Cannot Build PyGAMMA Until This Is Taken Care Of.; \
		echo ; \
		false; \
	fi
	@echo
PYS_LDFLAGS = -G $(LDFLAGS)
PYS_LDLIBS  = $(LAB_LIBS) -lCstd -lCrun -lm


${BINDIR}/$(PYS_NAME): $(PYSGAMMA_CHK) $(PYSG_OBJS) $(PYSG_WRAP_OBJS)
	@echo
	@echo "Making PYSGAMMA Shared $(PYS_NAME) From Pythonized and Swigged Object Files"
	@echo
	$(CXX) ${PYS_LDFLAGS} ${PYSG_LDFLAGS} -o ${BINDIR}/${PYS_NAME} $(PYSG_OBJS) $(PYSG_WRAP_OBJS) $(PYS_LDLIBS)

#########################################################################
#
#               Ensure Needed Binary Directory Exists
#
#########################################################################

${BINDIR}: force
	@if (test ! -d ./$(BINDIR)) \
        then \
          echo "First Constructing Binary Directory For GAMMA"; \
          echo "Making Directory ./$(BINDIR)"; \
          mkdir $(BINDIR); \
          chmod ug+w $(BINDIR); \
	fi

force:

#########################################################################
#
#                     Build And Run The Tests
#
# The purpose of this target is to run the GAMMA test program(s). The
# test programs painstakingly run through each GAMMA module and do some
# rather extensive checking for validity.
#
# The program will be named TestGAMMA and should be found in the GAMMA
# test directory. This target should always run the test program. It
# will insure that the GAMMA library is up to date, that the
# GAMMA library exists in the PATH, and that the test executable is up
# to date.
#
# $(exeext)   = Executable extension on system.
# CXX         = C++ compiler name
# CXXFLAGS    = C++ compiler flags
# LDFLAGS     = Linker flags
#
# -o $@       : Sets output executable name to match the name of
#               the target
# -L          : Look in local directory for libraries
# -lgamma    : Link to library lgamma
#
#########################################################################

#########################################################################
#
# New Test routines
#

TEST_PGM = $(BINDIR)/runtests $(BINDIR)/hstat \
       $(BINDIR)/sl_exchange7_L $(BINDIR)/par_xixA $(BINDIR)/xix_dec

test:  ${BINDIR}/${LIB_NAME} $(TEST_PGM) run


TESTS       = Tests
TESTSBASE   = testsuite runtests
TESTSHDR    = testsuite.h

TESTCCFILES   += $(addsuffix .cc, $(addprefix $(TESTS)/, $(TESTSBASE)))
TESTOBJFILES  += $(addsuffix .too,  $(addprefix $(BINDIR)/, $(TESTSBASE)))

ZTESTSBASE   = homo_stat sl_exchange7_L par_xixA xix_dec
ZTESTSHDR    =

ZTESTCCFILES += $(addsuffix .cc, $(addprefix $(TESTS)/, $(ZTESTSBASE)))
ZTESTOBJFILES  += $(addsuffix .too,  $(addprefix $(BINDIR)/, $(ZTESTSBASE)))

ALLTESTOBJFILES = $(TESTOBJFILES) $(ZTESTOBJFILES)

CXX_TESTFLAGS = $(SRCDIRINCL) $(OPTIMIZE) $(BLASINCL_FLAG) $(LAPACKINCL_FLAG)

LD_TESTFLAGS = $(LDFLAGS)


${BINDIR}/%.too : ${SRCDIR}/${TESTS}/%.cc
	$(CXX) -w -c $(CXX_TESTFLAGS) $< -o $@

$(BINDIR)/runtests : $(TESTOBJFILES)
	$(CXX) -w -Bstatic -o $@ $(LD_TESTFLAGS) $(TESTOBJFILES) -L${BINDIR}/ -lgamma -Bdynamic $(LAB_LIBS)

$(BINDIR)/hstat : $(BINDIR)/homo_stat.too
	$(CXX) -w -Bstatic -o $@ $(LD_TESTFLAGS) $(BINDIR)/homo_stat.too -L${BINDIR}/ -lgamma -Bdynamic $(LAB_LIBS)

$(BINDIR)/sl_exchange7_L : $(BINDIR)/sl_exchange7_L.too
	$(CXX) -w -Bstatic -o $@ $(LD_TESTFLAGS) $(BINDIR)/sl_exchange7_L.too -L${BINDIR}/ -lgamma -Bdynamic $(LAB_LIBS)

$(BINDIR)/par_xixA : $(BINDIR)/par_xixA.too
	$(CXX) -w -Bstatic -o $@ $(LD_TESTFLAGS) $(BINDIR)/par_xixA.too -L${BINDIR}/ -lgamma -Bdynamic $(LAB_LIBS)

$(BINDIR)/xix_dec : $(BINDIR)/xix_dec.too
	$(CXX) -w -Bstatic -o $@ $(LD_TESTFLAGS) $(BINDIR)/xix_dec.too -L${BINDIR}/ -lgamma -Bdynamic $(LAB_LIBS)

RUNDIR = ./../../src/Tests

run: force
	@echo
	@echo Running tests:
	@echo
	@cd $(RUNDIR); $(PYTHON) run_tests.py -v -p $(BINDIR)/
	@echo


# Python test suite.

PY_RUNDIR = ./../../src/pyTests

pytest: force
	@echo
	@echo Running Python tests:
	@echo
	@cd $(PY_RUNDIR); $(PYTHON) ../Tests/run_tests.py -v
	@echo


#########################################################################
#
#                           Installation
#
# This is the typical Gnuish install target with some added GAMMA
# features. There are 3 items of concern here. 1.) The GAMMA script,
# 2.) The GAMMA shared library and 3.) The GAMMA headers. Here are
# details regarding each.
#
# 1.) The GAMMA script (gamma)
# ============================
# The GAMMA script (gamma) is that which is run when programs linking
# to GAMMA are compiled. This script is tailored for a particular
# version and installation, that is, it knows which GAMMA version it
# is associated with & where the library + headers for that version are
# located. For example, the version may be 4.2.0 & compiled for Linux.
# The library headers are either placed in a Linux system directory
# (e.g. /usr/include/gamma-4.2.0) or they are left in a GAMMA directory
# (e.g. /E/gamma/gamma-4.2.0/src). Similarly the library itself is put
# in either the Linux system directory (e.g. /usr/lib) or left in a
# GAMMA directory (e.g. /E/gamma/gamma-4.2.0/i686-pc-linux-gnu). Since
# shared obeject libraries are preferable on Linux systems, the system
# directories are best because then their location is standard. If, for
# reasons unknown, a static library is preferred then the either the
# system or GAMMA directory setup can be used.  The only reason to have
# the GAMMA directory settings is if the user shares disks between
# different systems and wishes them all to have similar locations.
#
# When using the Gnu configuration process all of this information is
# set into the GAMMA script automatically. Since we do not run
# ./configure for this Makefile build, the GAMMA script must be
# configured directly herein. All that entails is setting the variable
# GAMMA_PATH with the script. This is done using sed on the precursor
# file gamma.in. When "make install" is run the target will convert
# gamma.in into gamma and set the proper directory information therein.
#
# In addition, the GAMMA script must be placed into a directory that
# is commonly known to most users. So, gamma is copied to ${INSTALLDIR}
# which is typically /usr/local/bin.
#
# 2.) The GAMMA shared library (libgamma.so)
# ==========================================
# When building GAMMA based programs using the GAMMA script gamma, the
# programs will be linked to the GAMMA library. This all works quite
# simply if linking is done to the static library, libgamma.a, because
# the GAMMA script tell the compiler to look in the Linux binaries
# directory for the GAMMA library. Once the program has been compiled
# into an executable it contains the library functionality (at the
# expense of large executable file size) and no longer needs the static
# library. Hence the executable program can be moved anywhere on the
# system and run without troubles. They have the library embedded in
# them.	All that is required for all this to work is that the gamma
# script (or the compilation command) used to build GAMMA based
# programs know where the static library is located. This can be either
# a GAMMA directory (e.g /E/gamma/gamma-4.2.0/i686-pc-linux-gnu) or a
# Linux system directory (e.g. /usr/lib)
#
# The same is NOT true when GAMMA based programs link to the shared
# library, libgamma.so. THIS IS THE DEFAULT LINKAGE IN THIS GAMMA
# SCRIPT. In this instance the executable program will not contain
# the functions of the library (this keeps the executable file size
# small). Rather it just keeps a record that it needs a library with a
# particular name (the library soname). When the executable is run it
# will attempt to use the library having the name (soname) specified
# when the executable was built and linked. As a consequence, the
# system must know where to find the shared library when the executable
# is run. Since the executable may be moved around on the system, one
# must tell the system where the library exists.
#
# To accomplish this, libgamma.so needs to be put in a standard
# directory where shared libraries are kept. For Linux this is usually
# /usr/local/lib
# Furthermore, libgamma.so should simply exist as a symbolic link to
# the library that has the soname that the executable needs. Even more
# complicated, the true library naming should use GAMMA version
# numbering in it (libgamma.so.#.#.#, the last number optional). It is
# this latter file that is built when GAMMA itself is compiled. So, the
# actual physical library file will be /usr/lib/libgamma.so.#.#(.#),
# the file with the soname will be /usr/lib/libgamma.so.# and exist as
# a symbolic link to the physical library, and the basic library
# /usr/lib/libgamma.so will be a symbolic link to the library having
# the soname name, /usr/lib/libgamma.so.#. Whew. In any event, the
# GAMMA directory (e.g /E/gamma/gamma-4.2.0/i686-pc-linux-gnu) is NOT
# used in this case to keep the shared object libraries on the system,
# although they may be built in such a directory.
#
# 3.) The GAMMA headers (*.h)
# ==========================================
# When building GAMMA based programs using the GAMMA script gamma, the
# programs will include the GAMMA source headers. Whether linking to
# the static or shared library, the headers must exist in a known
# location. On some systems this can be the GAMMA source directory
# (e.g /E/gamma/gamma-4.2.0/src) or it can be a Linux system directory
# (e.g. /usr/include/gamma-4.2.0). It is the latter directory that is
# used in typical Linux system GAMMA binary installations.
#
# To summarize, this target (install) does the following:
#
# a.) Determines the GAMMA version directory (this is GAMMA_PATH)
# b.) Converts gamma.in into gamma where GAMMA_PATH is set properly.
# c.) Copies gamma to INSTALLDIR (e.g. /usr/local/bin).
# d.) Copies the GAMMA shared library to LIBINSTALLDIR (e.g. /usr/lib)
#     This is the full library name, libgamma.so.#.#.#
# e.) Runs LDCONFIG (/sbin/ldconfig) to set up a symbolic link between
#     the built library & a file with the soname. That is (in lib dir)
#                libgamma.so.# --> libgamma.so.#.#.#
# f.) Makes a final symbolic link between the common library name &
#     the soname library. Again this is done in LIBINSTALLDIR
#                    libgamma.so --> libgamma.so.#
#
#########################################################################

GSCRIPT = gamma

####################################################################
#
# 	                     Targets For The GAMMA Headers
#
#####################################################################

installhdrs: force
	@if (test -d ${HDRINSTALLDIR}) \
	then echo The GAMMA Headers Directory Exists. ${HDRINSTALLDIR}; \
	else echo Creating GAMMA Headers Directory. ${HDRINSTALLDIR}; \
	     mkdir ${HDRINSTALLDIR}; \
	     for i in ${MODULEDIRS}; do \
                echo  "* Making Subdirectory $$i"; \
                mkdir ${HDRINSTALLDIR}/$$i; \
	        done; \
	fi
	@${MAKE} ${HDRINSTALLDIR}/gamma.h
	@${MAKE} ${HDRINSTALLDIR}/GamGen.h
	@${MAKE} ${INS_HDRS}

##########################################################################
#
#                           Targets For The GAMMA Script
#
# gamma:		  Construct gamma script from gamma.in
# ${BININSTALLDIR}/gamma: Copy local gamma script to script in BININSTALLDIR
# script:		  Same as gamma
# installscript:	  Same as ${BININSTALLDIR}/gamma
#
##########################################################################

SEDSUB_VERSION = 's|VERSION_NUMBER_GOES_HERE|${GAMMAVERSION}|g'

${BININSTALLDIR}/gamma: gamma force
	@echo Copying The GAMMA Script To The Bin Directory ${BININSTALLDIR}
	@cp -f ./${GSCRIPT} ${BININSTALLDIR}
	@perl -i -p -e ${SEDSUB_VERSION} ${BININSTALLDIR}/${GSCRIPT}
	@chmod a+x ${BININSTALLDIR}/${GSCRIPT}

script: gamma
installscript: ${BININSTALLDIR}/gamma

#######################################################################
#
# 	                     Targets For The GAMMA Library
#
#######################################################################

${LIBINSTALLDIR}/${SO_NAMEXXX}: ${BINDIR}/${SO_NAMEXXX}
	@echo Copying The GAMMA Shared Object Library To The Lib Directory ${LIBINSTALLDIR}
	@cp -f ${BINDIR}/${SO_NAMEXXX} ${LIBINSTALLDIR}

installsolib: force
	@if (test -d ${LIBINSTALLDIRV}) \
        then echo The System GAMMA Library Directory Exists. ${LIBINSTALLDIRV}; \
        else echo Creating GAMMA Library Directory On System. ${LIBINSTALLDIRV}; \
             echo * mkdir ${LIBINSTALLDIRV}; \
             mkdir ${LIBINSTALLDIRV}; \
        fi
	@if (test -d ${LIBINSTALLDIR}) \
        then echo The System GAMMA Library Directory Exists. ${LIBINSTALLDIR}; \
        else echo Creating GAMMA Library Directory On System. ${LIBINSTALLDIR}; \
             echo * mkdir ${LIBINSTALLDIR}; \
             mkdir ${LIBINSTALLDIR}; \
        fi
	@echo " * Copying Shared Library To ${LIBINSTALLDIR}"
	@cp -f ${BINDIR}/${SO_NAMEXXX} ${LIBINSTALLDIR}
	@echo "Making Links For Shared Library"
	@rm -f ${LIBINSTALLDIR}/${SO_NAME}
	@rm -f ${LIBINSTALLDIR}/${SO_NAMEX}
	@rm -f ${LIBINSTALLDIR}/${SO_NAMEXX}
	@ln -s ${LIBINSTALLDIR}/${SO_NAMEXXX} ${LIBINSTALLDIR}/${SO_NAMEXX}
	@ln -s ${LIBINSTALLDIR}/${SO_NAMEXXX} ${LIBINSTALLDIR}/${SO_NAMEX}
	@ln -s ${LIBINSTALLDIR}/${SO_NAMEXXX} ${LIBINSTALLDIR}/${SO_NAME}

${LIBINSTALLDIR}/${LIB_NAME}: ${BINDIR}/${LIB_NAME}
	@echo Copying The GAMMA Static Object Library To The Lib Directory ${LIBINSTALLDIR}
	@cp -f ${BINDIR}/${LIB_NAME} ${LIBINSTALLDIR}

installlib: force
	@if (test -d ${LIBINSTALLDIRV}) \
	then echo The System GAMMA Library Directory Exists. ${LIBINSTALLDIRV}; \
	else echo Creating GAMMA Library Directory On System. ${LIBINSTALLDIRV}; \
             mkdir ${LIBINSTALLDIRV}; \
	fi
	@if (test -d ${LIBINSTALLDIR}) \
	then echo The System GAMMA Library Directory Exists. ${LIBINSTALLDIR}; \
        else echo Creating GAMMA Library Directory On System. ${LIBINSTALLDIR}; \
             mkdir ${LIBINSTALLDIR}; \
	fi
	@echo " * Copying Static Library To ${LIBINSTALLDIR}"
	@cp -f ${BINDIR}/${LIB_NAME} ${LIBINSTALLDIR}

######################################################################
#
# 	                    Targets For The GAMMA Installation
#
######################################################################

install: force
	@echo
	@echo "           Installing GAMMA"
	@echo
	@echo " * GAMMA Script Installation Directory Is ${BININSTALLDIR}"
	@echo " * GAMMA Binary Library Installation Directory Is ${LIBINSTALLDIR}"
	@echo " * GAMMA Headers Installation Directory Is ${HDRINSTALLDIR}"
	@${MAKE} installscript
	@${MAKE} installhdrs
	@${MAKE} installsolib
	@${MAKE} installlib
	@echo " * Finished. Compile Command gamma Should Now Work"
	@echo


# PYSGDIST is the directory where files built for distribution should reside.
# They'll be used by PyGamma's setup.py.
PYSGDIST := $(ROOTDIR)/pygamma/dist_staging

# OSINFO, PYTHONINFO, and UNDERLINE help to construct the build notes (which
# are just informative and only need to be human-readable, not
# programmatically readable).
OSINFO := ${shell uname -a}
PYTHONINFO := ${shell ${PYTHON} -c "import sys; print(sys.version)"}
UNDERLINE := "------------------------------------------"

${PYSGDIST}: force
	@if (test ! -d ./$(PYSGDIST)) then \
          mkdir $(PYSGDIST); \
          chmod ug+w $(PYSGDIST); \
	fi

pysgdist: ${PYSGDIST} ${BINDIR} ${BINDIR}/$(PYS_NAME)
	@# Copy files to the dist dir
	@cp -f  ${BINDIR}/pygamma.py   $(PYSGDIST)
	@cp -f  ${BINDIR}/$(PYS_NAME)  $(PYSGDIST)
	@# Write a tiny Python file that contains the PyGAMMA version
	@echo "__version__ = '$(GAMMAVERSION)'" > $(PYSGDIST)/pygamma_version.py
	@# Copy files to the distribution staging area and construct the build notes
	@# Note that the first line below creates build_notes.txt, the
	@# second and subsequent lines append to it.
	@echo "OS Info (uname -a)"        >  $(PYSGDIST)/build_notes.txt
	@echo $(UNDERLINE)                >> $(PYSGDIST)/build_notes.txt
	@echo "${OSINFO}"                 >> $(PYSGDIST)/build_notes.txt
	@echo "Compiler Info (${CXX} -v)" >> $(PYSGDIST)/build_notes.txt
	@echo $(UNDERLINE)                >> $(PYSGDIST)/build_notes.txt
	@echo "${shell ${CXX} -v 2>&1}"   >> $(PYSGDIST)/build_notes.txt
	@echo "Python Info (sys.version)" >> $(PYSGDIST)/build_notes.txt
	@echo $(UNDERLINE)                >> $(PYSGDIST)/build_notes.txt
	@echo "${PYTHONINFO}"             >> $(PYSGDIST)/build_notes.txt

########################################################################
#
#                                Clean
#
# This is the typical Gnuish clean target. It will remove any libraries
# (dynamic and static), all object files, and the test program(s)
#
########################################################################

clean:
	@rm -f ${LIB_OBJS}
	@rm -f ${SO_OBJS}
	@rm -f $(PYSG_OBJS) $(PYSG_WRAP_SRCS) $(PYSG_WRAP_OBJS) $(PYSGAMMA_CHK)
	@rm -f $(ALLTESTOBJFILES)
	@rm -f *.asc
	@rm -f *.gnu

strip:
	@strip ${BINDIR}/${LIB_NAME}
	@strip ${BINDIR}/${SO_NAME}
	@strip ${LIBINSTALLDIR}/${SO_NAMEXXX}

testclean:
	@rm -f $(ALLTESTOBJFILES)

distclean:
	@${MAKE} clean
	@rm -f ${BINDIR}/${LIB_NAME}
	@rm -f $(TEST_PGM)
	@rm -f ${BINDIR}/${SO_NAME}
	@rm -f ${BINDIR}/${SO_NAMEX}
	@rm -f ${BINDIR}/${SO_NAMEXX}
	@rm -f ${BINDIR}/${SO_NAMEXXX}
	@rm -f ${BINDIR}/${PYS_NAME}

force:
