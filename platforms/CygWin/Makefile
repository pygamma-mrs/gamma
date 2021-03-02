#########################################################################
#                                                                       #
#          GAMMA Primal Makefile: CygWin ==> (Py)GAMMA Library		#
#                                                                       #
# This is a primitive hand made Makefile to explicitly show how one     #
# builds and uses GAMMA as either a static library or DLL and how one	#
# builds and uses PyGAMMA. It assumes that one is using the GCC C++	#
# compiler provided with (& amidst) the Cygwin environment. At least,	#
# the binaries will be placed in a CygWin related directory. So, run	#
# this makefile from a CygWin bash shell so GCC, make, etc. are from	#
# Cygwin. As for PyGAMMA, one must have Python (Cygwin NOT Windows)	#
# as well as Boost.Python on the system.				#
#                                                                       #
# --------------------------------------------------------------------- #
#                                                                       #
#                             Quick Info                                #
#                                                                       #
#  1. Type "make" to build GAMMA DLL library. 		  	#
#  2. Type "make demo" to compile and run the GAMMA demo program(s).    #
#  3. Type "make test" to compile and run the GAMMA test program(s).    #
#  4. Type "make clean" to remove all object files.			#
#  5. Type "make install" make this build your default GAMMA.		#
#  6. The sources are assumed in SRCDIR == ../../src			#
#  7. The object are be produced in BINDIR ==  ../../i686-pc-cygwin	#
#  8. The library files are produced in BINDIR ==  ../../i686-pc-cygwin	#
#  9. Any executable files (test,demo) are produced in this directory.	#
# 10. The important static library file is libgamma.a 			#
# 11. The important DLL library files are gamma.dll and libgamma.dll.a	#
# 12. The important PyGMAMA DLL library file is pygamma.dll.		#
# 13. Produces the following: a.) libgamma.a - The GAMMA static library	#
#			      b.) gamma.dll  - The GAMMA DLL		#
#			      c.) libgamma.dll.a - GAMMA imports lib.   #
#                             d.) pygamma.dll - PyGAMMA DLL             #
#                             e.) DemoGAMMA.exe - The demo program	#
#                             f.) TestGAMMA.exe - The test program	#
#                                                                       #
# --------------------------------------------------------------------- #
#                                                                       #
#                                 Notes					#
#                                                                       #
# This is NOT a standalone makefile. It includes an additional file	#
# ../../make/MakeMods which is part of standard GAMMA. This latter file #
# defines which files are part of the GAMMA source codes. Changes in    #
# MakeMods will evoke changes in the build herein.                      #
#									#
# all:									#
# make all: These both run the default target all (first target). It    #
#           will compile all source codes in GAMMA as defined in file	#
#           MakeMods into object files. It will then take the object	#
#           files and convert them into a static library named		#
#           libgamma.a. These same object files will then be combined	#
#           into a dynamic library (DLL) with its assocatied imports    #
#           library, gamma.dll and libgamma.dll.a respectively. If we   #
#           are able to find Python and Boost.Python, then the PyGAMMA  #
#           DLL is also built, pygamma.dll. All libraries and object	#
#           files will be made in the directory specified by BINDIR.    #
#									#
# make clean: This will delete all object files from any earlier build. #
#             It is important to run make clean since the object files  #
#             can use lots of disk space as well as interfere with      #
#             other GAMMA builds (not part of this Makefile)		#
#									#
# make demo : This will build and run a set of programs that demo the   #
#             the modules provided with GAMMA.				#
#                                                                       #
# make test: This will build and run a set of programs that test the	#
#             the modules provided with GAMMA.				#
#                                                                       #
#########################################################################

#########################################################################
#                                                                       #
#                Semi-Standard Gnuish Makefile Variables                #
#									#
# CXX         = C++ compiler (used in default compiler rules)		#
# CPPFLAGS    = C preprocessor flags (used in default C/C++ rules)	#
# CXXFLAGS    = C++ compilation flags (used in default C++ builds)	#
#									#
# The default rule in make that uses the variables CXX, CPPFLAGS and	#
# CXXFLAGS to generate object files from .cc files is as follows:	#
#									#
# %.o : %.cc								#
#        $(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@			#
#									#
#########################################################################

DEBUG    = -g #-O2

CC       = gcc
CFLAGS   = $(DEBUG)

CXX       = c++
CPPFLAGS  = -I../../src
CXXFLAGS  = $(DEBUG)

LDFLAGS  = 
exeext   = .exe
defext   = .def

AR       = ar
ARFLAGS  = r
RANLIB   = ranlib

CPYGFLAGS   = -DUSE_DL_IMPORT -DBOOST_PYTHON_DYNAMIC_LIB -DPYGAMMA=1 -g -O0 -fno-inline \
              -I../../src -isystem ${BOOSTINCDIR} -isystem ${PYINCDIR}
#CPYGFLAGS   = -I${BOOSTINCDIR} -I${PYINCDIR} -DPYGAMMA=1
PYG_LDFLAGS = -L${BOOSTLIBDIR} -L${PYLIBDIR} -L${PYLIB2DIR} -lpython2.2.dll -lboost_python
PYGLDLIBS   = 

NPD      = --no-print-directory

#########################################################################
#									#
#	           Names Of Various Targets To Build			#
#									#
# These names are associated with the library we plan to build.		#
#									#
# Static Library Variables						#
# ========================						#
# LIB_NAME    = Name of library to be generated/targeted.  (libgamma.a)	#
#									#
# Dynamic Library (DLL) Variables					#
# ===============================					#
# DLL_NAME    = Name of DLL to be generated/targeted.      (gamma.dll)  #
# DLL_EXP_LIB = Name of exports library generated.    (libgamma.dll.a)  #
# DLL_EXP_DEF = Name of definiton file to use.             (gamma.def)  #
#									#
# PyGAMMA Dynamic Library (DLL) Variables				#
# =======================================				#
# PYG_NAME    = Name of PyGAMAM DLL to be generated       (pygamma.dll) #
# PYG_EXP_LIB = Name of exports library generated.   (libpygamma.dll.a) #
# PYINCDIR    = Header include directory for Python			#
# PYLIBDIR    = Library include directory for Python			#
# PYLIBDIR2   = Library include directory for Python			#
# BOOSTINCDIR = Header include directory for Boost.Python		#
# BOOSTLIBDIR = Library include directory for Boost.Python		#
#									#
#########################################################################

LIB_NAME    = libgamma.a

DLL_NAME    = gamma.dll
DLL_EXP_LIB = libgamma.dll.a
DLL_EXP_DEF = gamma.def

PYG_NAME    = pygamma.dll
PYG_EXP_LIB = libpygamma.dll.a
PYINCDIR    = C:/cygwin/usr/include/python2.2
PYLIBDIR    = C:/cygwin/lib/python2.2
PYDLLDIR    = C:/cygwin/lib/python2.2
PYLIB2DIR   = C:/cygwin/lib/python2.2/config
BOOSTINCDIR = "C:/Program Files/Boost/boost_1_29_0"
BOOSTLIBDIR = "C:/Program Files/Boost/boost_1_29_0/libs/python/build/bin/boost_python.dll/gcc/release/inlining-on/runtime-link-dynamic" 

#########################################################################
#                                                                       #
#                       Gamma Related Variables                         #
#                                                                       #
#                                                                       #
# These variables are used to specify the GAMMA sources. Exactly what   #
# sources are associated with any one GAMMA module is contained in the  #
# external makefile MakeMods. From these we can construct the names of  #
# module source files and object files. These are directory specific.   #
#                                                                       #
# MakeMods    = File listing all current GAMMA source files.            #
# SRCDIR      = Base directory of GAMMA source codes.                   #
# BINDIR      = Base directory of Cygwin binaries.                      #
# LIB_SRCS    = All code files to be put into the GAMMA static library.	#
# LIB_OBJS    = All object files to be put in the GAMMA static library.	#
# DLL_SRCS    = All code files to be put into the GAMMA DLL.		#
# DLL_OBJS    = All object files to be put in the GAMMA DLL.		#
#                                                                       #
# The variable ALLCCFILES should contain all GAMMA source .cc files.    #
# The variable ALLOFILES should contain all GAMMA object .o files.	#
# The source file names will have the source sub-directory as a prefix. #
#     but not the full path relative to this directory. 		# 
#                                                                       #
# The included makefile supplies NO targets or rules of any kind.	#
#                                                                       #
#########################################################################

include ../../make/MakeMods
SRCDIR   =  ../../src
BINDIR   = ../../i686-pc-cygwin

LIB_SRCS = $(ALLCCFILES)
LIB_OBJS = $(addsuffix .o, $(addprefix $(BINDIR)/, $(ALLFILES)))

DLL_SRCS = $(ALLCCFILES)
DLL_OBJS = $(addsuffix .o, $(addprefix $(BINDIR)/, $(ALLFILES)))

PYG_SRCS =  $(ALLCCFILES)
PYG_SRCS += $(addsuffix .cc, $(addprefix $(SRCDIR)/, PyGAMMAIF))
PYG_OBJS =  $(addsuffix .po, $(addprefix $(BINDIR)/, $(ALLFILES)))
PYG_OBJS += $(addsuffix .po, $(addprefix $(BINDIR)/, PyGAMMAIF))

#########################################################################
#									#
#                           GAMMA Module Directories			#
#									#
# These definitions are required for the pattern rules (later in this   #
# file) that dictate how source code is compiled. Since we like to keep #
# GAMMA code in several different subdirectories yet compile these into #
# object code in the same directory, these defintions are required to   #
# link specific source code directory names to the single directory we  #
# use for binary output (BINDIR).					#  
#									#
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

#########################################################################
#									#
#              Build GAMMA Static Library & DLL Pattern Rules		#
#									#
# These rules specify how we build object file (*.o) from C++ code	#
# files (*.cc). Gnu make has a built in implicit pattern rule to do	#
# this task:								#
#									#
# %.o : %.cc								#
#        $(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@			#
#									#
# But that does not suffice for good reason. The GAMMA source files	#
# reside in module directories off of the source directory SRCDIR.	#
# Second, we want all the object files to be built in the binaries 	#
# directory, BINDIR. Thus we must build ${BINDIR}/*.o from a code file	#
# ${SRCDIR}/${MODDIR}/*.cc. The main GAMMA  makefile handles this by	#
# recursing through each module directory and letting it decide how to	#
# build itself. Here we will define a rule for each module directory.   #
#									#
# For CygWin builds we need no special flags to distinguish what is	#
# needed when compiling files for a static library versus a DLL.	#
#									#
#########################################################################

${BINDIR}/%.o : ${SRCDIR}/${BASICS}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${MATRIX}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${HSPACE}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${LEVEL1}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${LSPACE}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${LEVEL2}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${BLOCH}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${GAMIO}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${BWR}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${PULSES}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${INTRANK2}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${FLOQUET}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${GRADIENTS}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${MULTISYS}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${ESRLIB}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${DEPREC}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

${BINDIR}/%.o : ${SRCDIR}/${TESTING}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

#########################################################################
#									#
#                     Build PyGAMMA DLL Pattern Rules			#
#									#
# These rules specify how we build object file (*.po) from C++ code	#
# files (*.cc). Gnu make has a built in implicit pattern rule to do	#
# this task:								#
#									#
# %.o : %.cc								#
#        $(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@			#
#									#
# But that does not suffice for good reasons. The GAMMA source files	#
# reside in module directories off of the source directory SRCDIR.	#
# Second, we want all the object files to be built in the binaries 	#
# directory, BINDIR. Thus we must build ${BINDIR}/*.po from a code file	#
# ${SRCDIR}/${MODDIR}/*.cc. The main GAMMA  makefile handles this by	#
# recursing through each module directory and letting it decide how to	#
# build itself. Here we will define a rule for each module directory.   #
# Also we are using the suffix .po rather than .o in order to clearly   #
# distinguish binary files associated with PyGAMMA from those built for #
# for GAMMA. Lastly, and this is the reason the binaries are different, #
# we need to in include the preprocessor defintion PYGAMMA=1 so that    #
# additonal sources for PyGAMMA are included when the source files are  #
# compiled. This is done with the flag CPYGFLAGS, and the result makes  #
# the *.po files dependent upon the Python and Boost.Python libraries.  #
#									#
#########################################################################

${BINDIR}/%.po : ${SRCDIR}/${BASICS}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) ${CPYGFLAGS} $< -o $@

${BINDIR}/%.po : ${SRCDIR}/${MATRIX}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) ${CPYGFLAGS} $< -o $@

${BINDIR}/%.po : ${SRCDIR}/${HSPACE}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) ${CPYGFLAGS} $< -o $@

${BINDIR}/%.po : ${SRCDIR}/${LEVEL1}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) ${CPYGFLAGS} $< -o $@

${BINDIR}/%.po : ${SRCDIR}/${LSPACE}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) ${CPYGFLAGS} $< -o $@

${BINDIR}/%.po : ${SRCDIR}/${LEVEL2}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) ${CPYGFLAGS} $< -o $@

${BINDIR}/%.po : ${SRCDIR}/${BLOCH}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) ${CPYGFLAGS} $< -o $@

${BINDIR}/%.po : ${SRCDIR}/${GAMIO}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) ${CPYGFLAGS} $< -o $@

${BINDIR}/%.po : ${SRCDIR}/${BWR}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) ${CPYGFLAGS} $< -o $@

${BINDIR}/%.po : ${SRCDIR}/${PULSES}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) ${CPYGFLAGS} $< -o $@

${BINDIR}/%.po : ${SRCDIR}/${INTRANK2}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) ${CPYGFLAGS} $< -o $@

${BINDIR}/%.po : ${SRCDIR}/${FLOQUET}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) ${CPYGFLAGS} $< -o $@

${BINDIR}/%.po : ${SRCDIR}/${GRADIENTS}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) ${CPYGFLAGS} $< -o $@

${BINDIR}/%.po : ${SRCDIR}/${MULTISYS}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) ${CPYGFLAGS} $< -o $@

${BINDIR}/%.po : ${SRCDIR}/${ESRLIB}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) ${CPYGFLAGS} $< -o $@

${BINDIR}/%.po : ${SRCDIR}/${DEPREC}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) ${CPYGFLAGS} $< -o $@

${BINDIR}/%.po : ${SRCDIR}/${TESTING}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) ${CPYGFLAGS} $< -o $@

${BINDIR}/%.po : ${SRCDIR}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) ${CPYGFLAGS} $< -o $@

#########################################################################
#									#
#                              Main Targets				#
#									#
# These targets produce the GAMMA and PyGAMMA libraries.		#
#									#
# Supplied Targets:							#
# =================							#
#									#
# all: - Builds the GAMMA library as a DLL, gamma.dll			#
# lib: - Builds only the GAMMA static library, libgamma.a		#
# dll: - Builds only the GAMMA dynamic library, gamma.dll		#
# pyg: - Builds only the PyGAMMA dynamic library, pygamma.dll           #
#									#
# Variables Needed:							#
# =================							#
# BINDIR   = CygWin binaries directory (../../i686-pc-cygwin)		#
# LIB_NAME = GAMMA static library name (libgamma.a)			#
# DLL_NAME = GAMMA dynamic library name (gamma.dll)			#
# PYG_NAME = PyGAMMA dynamic library name (pygamma.dll)			#
# DLL_NAME = GAMMA dynamic library name (gamma.dll)			#
#									#
#########################################################################

all: $(BINDIR)/$(DLL_NAME)

lib: $(BINDIR)/$(LIB_NAME)
dll: $(BINDIR)/$(DLL_NAME)
pyg: $(BINDIR)/$(PYG_NAME)


${LIB_NAME}: ${BINDIR}/${LIB_NAME}
${DLL_NAME}: ${BINDIR}/${DLL_NAME}
${PYG_NAME}: ${BINDIR}/${PYG_NAME}

#########################################################################
#									#
#                       Library Related Variables			#
#									#
# These variables are used exclusively when building the static lib. 	#
#									#
# LIB_SRCS    = All code files to be put into the library.		#
# LIB_OBJS    = All object files to be put into the library.		#
# LIB_NAME    = Name of library to be generated/targeted.  (libgamma.a)	#
# LIB_OBJS    = Name of object files used in library.			#
#									#
# These variables are used exclusively when building the DLL. 		#
#									#
# DLL_LDFLAGS = Flags used in DLL build.     (-Wl,--export-all-symbols) #
# DLL_LDLIBS  = Any libraries the DLL to be build depends upon.         #
# DLL_SRCS    = All code files to be put into the DLL.                  #
# DLL_OBJS    = All object files to be put into the DLL. Note that they #
#               are one and the same as LIB_OBJ since we use no special #
#               compiler flags to generate them over the static library #
#									#
# These variables are used exclusively when building the PyGAMMA DLL.	#
#									#
# PYG_SRCS    = All code files to be put into the PyGAMMA library.	#
# PYG_OBJS    = All object files to be put into the PYGAMMA library.	#
#               These will have special names because they are built    #
#               with special flags, in particular _PYGAMMA being set.   #
#									#
# These variables are used generally. 					#
#									#
# BINDIR      = Where we will ultimately place or binaries. This is	#
#               set to replace a CygWin build.				#
# NPD         = Turns of directory name printing when recursing through #
#               other directories.                 			#
#									#
#########################################################################

DLL_LDFLAGS = -Wl,--export-all-symbols -Wl,--enable-auto-import
DLL_LDLIBS  =

#########################################################################
#									#
#                      Build The GAMMA Static Library 			#
#									#
# LIB_SRCS    = All code files to be put into the library.		#
# LIB_OBJS    = All object files to be put into the library.		#
# LIB_NAME    = Name of library to be generated/targeted.  (libgamma.a)	#
# LIB_OBJS    = Name of object files used in library.			#
#									#
# The static library build depends on the object files in LIB_OBJS. If  #
# or when these are up to date then the target uses ar to form the 	#
# library from the object files and then ranlib to adjust the library.	#
#									#
# The LIB_OBJS variable is set from LIB_SRCS earlier in this file.	#
# These are apparently implicitly dependent upon the source codes in	#
# LIB_SRCS (i.e. the .cc file with the same base name as the .o file)   #
# & have an implicit pattern rule as to how the .cc files are compiled	#
# into .o files. This pattern rule, as specified in the documentation	#
# for GNU make, is							#
#									#
# %.o : %.cc								#
#        $(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@			#
#									#
# This implicit rule can be replaced by defining a new pattern rule	#
# with the same target and prerequisites, but different commands. 	#
# If there are no commands following the rule, the rule is cancelled.	#
# When the new rule is defined, the built-in one is replaced. The new	#
# rule's position in the sequence of implicit rules is determined by 	#
# where you write the new rule and affects how they are chained.	#
#									#
#########################################################################

${BINDIR}/$(LIB_NAME): ${BINDIR} ${LIB_OBJS}
	@echo
	@echo "Making Static Library $(LIB_NAME) From Object Files"
	@echo
	${AR} ${ARFLAGS} ${BINDIR}/$(LIB_NAME) $(LIB_OBJS)
	@echo
	@echo "Adjusting Static Library $(LIB_NAME) For Random Access"
	@echo
	${RANLIB} ${BINDIR}/$(LIB_NAME)

#########################################################################
#                                                                       #
#                Build The GAMMA DLL And The Exports Library            #
#                                                                       #
# DLL_NAME    = Name of DLL to be generated/targeted.      (gamma.dll)  #
# DLL_EXP_LIB = Name of exports library generated.    (libgamma.dll.a)  #
# DLL_OBJS    = Name of object files used in DLL.                       #
# DLL_LDFLAGS = Flags used in DLL build.     (-Wl,--export-all-symbols) #
# DLL_LDLIBS  = Libraries DLL depends upon.                             #
#                                                                       #
# The DLL build depends on the object files in DLL_OBJS. If or when	#
# these are up to date then the target uses g++ to form the DLL from 	#
# the object files using appropriate linking commands.			#
#									#
# The DLL_OBJS variable is set from DLL_SRCS earlier in this file.	#
# These are one and the same as LIB_OBJS for the build so the two are	#
# identical. As with LIB_OBJS, the objects files are implicitly 	#
# dependent upon the source codes in DLL_SRCS (i.e. the .cc file with	#
# the same base name as the .o file) & have an implicit pattern rule as	#
# to how the .cc files are compiled into .o files. This pattern rule, 	#
# as specified in the documentation for GNU make, is			#
#									#
# %.o : %.cc								#
#        $(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@			#
#									#
# Since we need nothing additonal, this is quite adequate.		#
#									#
#########################################################################

$(BINDIR)/$(DLL_NAME) $(BINDIR)/$(DLL_EXP_LIB): $(DLL_OBJS)
	@echo
	@echo "Making DLL $(DLL_NAME) & Imports Library $(DLL_EXP_LIB) From Object Files"
	@echo
	${CXX} -shared -Wl,--out-implib=${BINDIR}/$(DLL_EXP_LIB) \
	-o ${BINDIR}/$(DLL_NAME) $(DLL_OBJS) $(DLL_LDFLAGS) $(DLL_LDLIBS)

#########################################################################
#                                                                       #
#                       Build The PyGAMMA DLL           		#
#                                                                       #
# The PyGAMMA build depends on the object files in PYG_OBJS. If or when	#
# these are up to date then the target uses g++ to form the DLL from	#
# the object files using g++ with appropriate linking commands.		#
#									#
# The PYG_OBJS variable is set from PYG_SRCS earlier in this file.	#
# These are NOT one and the same as LIB_OBJS nor DLL_OBJS. This is 	#
# because the sources are compiled into object files with the variable	#
# PYGAMMA defined. This causes the build to include addtional sources	#
# and headers required for PyGAMMA. The actual source code files from	#
# which PYG_OBJS are made, PYG_SRCS, is one and the same as LIB_SRCS &	#
# DLL_SRCS, but the PYGAMMA variable forces the preprocessor to use	#
# code commented out in the latter two.					#
#									#
# Normally objects files are implicitly dependent upon the source codes #
# That is, the files in DLL_OBJS depend upon DLL_SRCS (i.e. the .cc 	#
# file with the same base name as the .o file) & have an implicit 	#
# pattern rule as to how the .cc files are compiled into .o files. 	#
# But since PyGAMMA demands are different, the files in PYG_OBJS are	#
# set to end with .po rather than .o so that the implicit pattern rule 	#
# is not invoked. Rather we have earlier defined our own pattern rule	#
# to handle the compilation of .cc to .po. The default is		#
#									#
# %.o : %.cc								#
#        $(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@			#
#                                                                       #
# and the one we will invoke in the following target is			#
#                                                                       #
# %.po : %.cc								#
#        $(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(CPYGFLAGS) $< -o $@	#
#									#
# Supplied Targets:							#
# =================							#
#									#
# PyGAMMAChk:            Checks that the appropriate Boost.Python and 	#
#                        Python directories (header and library) exist.	#
# $(BINDIR)/$(PYG_NAME): Builds the PyGAMMA DLL in the Cygwin binaries  #
#                        directory. This will update all PyGAMMA object #
#			 files by compiling source codes if needed.     #
#                                                                       #
# Variables Needed:							#
# =================							#
# PYG_NAME    = Name of DLL to be generated/targeted.    (pygamma.dll)  #
# PYINCDIR    = Python header includes directory.                       #
# PYLIBDIR    = Python library directory.				#
# BOOSTINCDIR = Boost.Python header includes directory.			#
# BOOSTLIBDIR = Boost.Python library directory.				#
# BINDIR      = Cygwin binaries directory.				#
# PYG_OBJS    = PyGAMMA object files.					#
# PYG_EXP_LIB = PyGAMMA exports library.				#
# PYG_LDFLAGS = PyGAMMA linker flags.					#
# PYG_LDLIBS  = PyGAMMA additonal libraries.				#
#                                                                       #
#########################################################################

PyGAMMAChk:
	@echo
	@echo "Making PyGAMMA DLL $(PYG_NAME) From Pythonized Object Files"
	@echo
	@echo
	@echo
	@echo "* Checking For The Existence Of CygWin Built Python Directory......."
	@if (test -d ./${PYINCDIR}) \
	then echo "  Good. The Cygwin Python Includes Directory Exists. ${PYINCDIR}"; \
	else echo "  Bad. Cannot Find The Cygwin Includes Python Directory! ${PYINCDIR}"; \
	     echo "  Sorry, Unable To Build PyGAMMA Without Cygwin Python Installed!"; \
	     echo; \
             false; \
	fi
	@echo
	@echo "* Checking For The Existence Of CygWin Python Library Directory......."
	@if (test -d ./${PYLIBDIR}) \
	then echo "  Good. The Cygwin Python Library Directory Exists. ${PYLIBDIR}"; \
	else echo "  Cannot Find The Cygwin Python Library Directory! ${PYLIBDIR}"; \
	     echo "  Sorry, Unable To Build PyGAMMA Without Cygwin Python Installed!"; \
	     echo; \
             false; \
	fi
	@echo
	@echo "* Checking For The Existence Of CygWin Boost.Python Includes Directory......."
	@if (test -d ./${BOOSTINCDIR}) \
	then echo "  Good. The Cygwin Boost Includes Directory Exists. ${BOOSTINCDIR}"; \
	else echo "  Cannot Find The Cygwin Boost Includes Directory! ${BOOSTINCDIR}"; \
	     echo "  Sorry, Unable To Build PyGAMMA Without Boost.Python!"; \
	     echo; \
             false; \
	fi
	@echo
	@echo "* Checking For The Existence Of CygWin Boost.Python Library Directory......."
	@if (test -d ./${BOOSTLIBDIR}) \
	then echo "  Good. The Cygwin Boost Library Directory Exists."; \
	else echo "  Cannot Find The Cygwin Built Boost Library Directory!"; \
	     echo "  Make Sure You First Build Boost.Python With Cygwin."; \
	     echo "  Sorry, Unable To Build PyGAMMA Without Boost.Python!"; \
	     echo; \
             false; \
	fi
	@echo "  ${BOOSTLIBDIR}"
	@echo
	@echo "* Now Attempting Compile GAMMA Source Into PyGAMMA Object Code"
	@echo
	@echo


$(BINDIR)/$(PYG_NAME): PyGAMMAChk $(PYG_OBJS)
	@echo
	@echo "Making PyGAMMA DLL $(PYG_NAME) From Pythonized Object Files"
	@echo
	${CXX} -Wl,--export-all-symbols -Wl,--exclude-symbols,_bss_end__:_bss_start__:_data_end__:_data_start__ \
        -Wl,--out-implib=${BINDIR}/$(PYG_EXP_LIB) -g -shared -o ${BINDIR}/$(PYG_NAME) \
        $(PYG_OBJS) $(PYG_LDFLAGS) $(PYG_LDLIBS)

# sosi

#########################################################################
#									#
#                   Insure CygWin Binaries Directory Exists		#
#									#
#########################################################################

${BINDIR}: force
	@if (test ! -d ./$(BINDIR)) \
        then \
          echo "First Constructing CygWin Binaries Directory For GAMMA"; \
          echo "Making Directory ./$(BINDIR)"; \
          mkdir $(BINDIR); \
          chmod ug+w $(BINDIR); \
	fi

force:

#########################################################################
#									#
#                        Build And Run The Demos			#
#									#
# The purpose of this target is to run the GAMMA demo program(s). The   #
# demo programs painstakingly walk the user through each GAMMA module   #
# so that they can get a feel of what is available.                     #
#									#
# The program will be named DemoGAMMA and should be found in the GAMMA  #
# demo directory in the modules sub-directory. The target should always	#
# run the demo program. It will insure that the GAMMA library (DLL) is  #
# up to date, that the GAMMA library exists in the PATH, and that the   #
# demo executable is up to date.					#
#									#
# $(exeext)   = Executable extension on system.                  (.exe) #
# CXX         = C++ compiler name                                 (c++)	#
# CXXFLAGS    = C++ compiler flags                             $(DEBUG)	#
# LDFLAGS     = Linker flags (automatically used)			#
#									#
# DEMO_DIR     = Directory where demo program sources are		#
# DEMO_PGM    = Name of produced demo executable			#
# DEMO_SRC    = Full name of demo program source file			#
# DEMO_OBJ    = Full name of demo program object file			#
#									#
# -o $@       : Sets output executable name to match the name of 	#
#               the target (TestGAMMA.exe)				#
# -L          : Look in local directory for libraries			#
# -lgamma    : Link to library lgamma					#
#									#
#########################################################################

DEMO_DIR      = ../../demo/modules
DEMO_NAME    = DemoGAMMA
DEMO_PGM     = ${DEMO_NAME}${exeext}
DEMO_SRC     = $(DEMO_DIR)/${DEMO_NAME}.cc
DEMO_PY      = $(DEMO_DIR)/${DEMO_NAME}.py
DEMO_OBJ     = ./${DEMO_NAME}.o
DEMO_FLAGS   = -D_GAMMADEMO=1

DEMO_OBJ_S   = ./${DEMO_NAME}.static.o
DEMO_PGM_S   = ${DEMO_NAME}.static${exeext}

PYTHON       = python

# ---------------- Demo Targets With Default (DLL) Library -------------------
#              (Use DLL If Present, If Not Use Static Library)

%.o : ${DEMO_DIR}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -DUSEDLL=1 $< -o $@

demo:	
	@${MAKE} $(NPD) ${DEMO_PGM}
	@strip ./${DEMO_PGM}
	${DEMO_PGM} $(DEMO_DIR)

${DEMO_PGM}: ${DEMO_OBJ}
	$(CXX) -o $@ ${DEMO_OBJ} ${DEMO_FLAGS} $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -DUSEDLL=1 -L${BINDIR}/ -lgamma

# --------------------- Demo Targets With Static Library ---------------------

%.static.o : ${DEMO_DIR}/%.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

demostatic:	
	@${MAKE} $(NPD) ${DEMO_PGM_S}
	@strip ./${DEMO_PGM}
	${DEMO_PGM_S} $(DEMO_DIR)

${DEMO_PGM_S}: ${DEMO_OBJ_S}
	$(CXX) -o $@ ${DEMO_OBJ_S} ${DEMO_FLAGS} $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -L${BINDIR}/ -static -lgamma

# -------------------- Demo Targets With PyGAMMA Library ---------------------

demopython:	
	cp -f ${BINDIR}/$(PYG_NAME) .
	${PYTHON} ${DEMO_PY}

#########################################################################
#									#
#                     Build And Run The Tests				#
#									#
# The purpose of this target is to run the GAMMA test program(s). The   #
# test programs painstakingly run through each GAMMA module and do some	#
# rather extensive checking for validity.				#
#									#
# The program will be named TestGAMMA and should be found in the GAMMA  #
# test directory. This target should always run the test program. It	#
# will insure that the GAMMA library (DLL) is up to date, that the 	#
# GAMMA library exists in the PATH, and that the test executable is up 	#
# to date.								#
#									#
# $(exeext)   = Executable extension on system.                  (.exe) #
# CXX         = C++ compiler name                                 (c++)	#
# CXXFLAGS    = C++ compiler flags                             $(DEBUG)	#
# LDFLAGS     = Linker flags						#
#									#
# -o $@       : Sets output executable name to match the name of 	#
#               the target (TestGAMMA.exe)				#
# -L          : Look in local directory for libraries			#
# -lgamma    : Link to library lgamma					#
#									#
#########################################################################

TESTDIR      = ../../test
TEST_PGM     = $(TESTDIR)/TestGAMMA${exeext}
TEST_SRC     = $(TESTDIR)/TestGAMMA.cc
TEST_OBJ     = $(TESTDIR)/TestGAMMA.o

test:	
	@${MAKE} $(NPD) ${TEST_PGM}
	${TEST_PGM}

${TEST_PGM}: ${TEST_OBJ}
	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) ${TEST_OBJ} -L${BINDIR}/ -lgamma 


#########################################################################
#									#
#                            TestGAMMA.o 				#
#									#
# This is the target for building an object file from the test file.	#
#									#
# CXX         = C++ compiler name                                 (c++)	#
# CPPFLAGS    = C++ compiler flags                        (-I../../src)	#
# CXXFLAGS    = C++ compiler flags                             $(DEBUG)	#
#									#
#########################################################################

TestGAMMA.o: %o: %cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

#########################################################################
#									#
#                           Installation 				#
#									#
# This is the typical Gnuish install target. It will copy the GAMMA DLL #
# (gamma.dll) and the GAMMA script (gamma) into the installation 	#
# directory. Of course`,' first the DLL must be built and the GAMMA     #
# script configured. Since we do not run ./configure for this primal    #
# Makefile build, the script must be configured directly herein. All    #
# that entails is setting the variable GAMMA_PATH with the script. This #
# is done using sed on the precursor file gamma.in. So, this target     #
# does the following:                                                   #
#									#
# 1.) Determines the GAMMA version directory (this is GAMMA_PATH)       #
# 2.) Converts gamma.in into gamma where GAMMA_PATH is set properly.    #
# 3.) Copies gamma.dll from the binaries directory to the install dir.  #
# 4.) Copies gamma from this directory into the install directory.      #
#									#
# The default installation directory is /usr/local/bin.			#
#									#
#########################################################################

INSTALLBINDIR = /usr/local/bin
INSTALLLIBDIR = /usr/lib
GSRCDIR       = /cygdrive/c/gamma/gamma-4.1.0/src
GSCRIPT    = gamma
GSCRIPTIN  = gamma.in
CURRENTDIR = ${shell pwd}
GAMMA_PATH = ${shell cd ../..; pwd; cd primal/CygWin}
SEDSUB     = 's|CONFIGUREDGAMMAPATH|${GAMMA_PATH}|g'

script: force
	@echo " * The GAMMA Path Is Set to ${GAMMA_PATH}"
	@echo " * Setting GAMMA Script To Know GAMMA Path"
	@sed ${SEDSUB} ${GSCRIPTIN}>${GSCRIPT}

install: force
	@echo
	@echo "           Installing CygWin GAMMA"
	@echo
	@echo " * Binaries Installation Directory Is ${INSTALLBINDIR}"
	@echo " * Copying GAMMA DLL ${DLL_NAME} From ${BINDIR}"
	@cp ${BINDIR}/${DLL_NAME} ${INSTALLBINDIR}
	@echo " * Libraries Installation Directory Is ${INSTALLLIBDIR}"
	@echo " * Copying GAMMA DLL Imports Library ${DLL_EXP_LIB} From ${BINDIR}"
	@cp ${BINDIR}/${DLL_EXP_LIB} ${INSTALLLIBDIR}
	@echo " * The GAMMA Path Is Set to ${GAMMA_PATH}"
	@echo " * Setting GAMMA Script To Know GAMMA Path"
	@sed ${SEDSUB} ${GSCRIPTIN}>${GSCRIPT}
	@echo " * Copying GAMMA script ${GSCRIPT}"
	@cp -f ./${GSCRIPT} ${INSTALLBINDIR}

pyinstall: force
	@echo " * Copying PyGAMMA DLL ${PyG_NAME} From ${BINDIR}"
	@echo "   Into CygWin Python DLL directory ${PYDLLDIR}"
	@cp ${BINDIR}/${PYG_NAME} ${PYDLLDIR}

#########################################################################
#									#
#                                Clean 					#
#									#
# This is the typical Gnuish clean target. It will remove any libraries #
# (both dll and .a), any definition files, all object files, and the    #
# build of the test program(s)                                          #
#									#
#########################################################################

clean:
	@rm -f $(LIB_OBJS)
	@rm -f $(PYG_OBJS)
	@rm -f $(TEST_OBJ)
	@rm -f *.mif
	@rm -f *.asc
	@rm -f *.gnu
	@rm -f *.exe
	@rm -f *.o
	@rm -f *.~*
	@rm -f *stack*
	@rm -f $(DEMO_OBJ)
	@rm -f $(DEMO_PGM)
	@rm -f gamma

pyclean:
	@rm -f $(PYG_OBJS)

testclean:
	@rm -f $(TEST_PGM)
	@rm -f $(TEST_OBJ)

democlean:
	@rm -f $(DEMO_OBJ)
	@rm -f $(DEMO_PGM)
	@rm -f *.mif
	@rm -f *.asc
	@rm -f *.gnu
	@rm -f *.exe
	@rm -f *.o
	@rm -f *.~*
	@rm -f *stack*

distclean:
	@rm -f $(LIB_OBJS)
	@rm -f ${LIB_NAME}
	@rm -f $(DEMO_OBJ)
	@rm -f $(DEMO_PGM)
	@rm -f $(TEST_OBJ)
	@rm -f $(TEST_PGM)

dllclean:
	@rm -f ./${BINDIR}/*.dll
	@rm -f ./${BINDIR}/*.dll.*

crap:
	@echo "Here Is LIB_SRCS " ${LIB_SRCS}
	@echo "Here Is LIB_OBJS " ${LIB_OBJS}
	@echo "Here Is PYG_SRCS " ${PYG_SRCS}
	@echo "Here Is PYG_OBJS " ${PYG_OBJS}

#################################################################################
#										#
#                              GAMMA Distributions				#
#										#
# What follows has nothing to do with any builds of GAMMA. This stuff hass been #
# included in this Makefile only so that GAMMA distributions will include the   #
# distributed files in this directory. Do not worry about this for any builds   #
# and or runs of GAMMA programs.						#
#										#
#################################################################################

THISDIR       = CygWin
PRIMALPLAT    = "using CygWin development environment"
PRIMALMAKEDIR = ../make
MAINMAKEDIR   = ../../make
SUBDIRS       =
HDRFILES      =
CCFILES       =
HTMLFILES     = index.html
IMGFILES      = 
INCFILES      = Makefile CygWinHeap.reg gamma.in cleanup.bat
DOSCLEANFILES = Makefile gamma.in
UPDATEFILES   = ${HTMLFILES}
DISTSUBDIRS   =
MODLABEL      = ${THISDIR}

include ${PRIMALMAKEDIR}/MakeDists
include ${PRIMALMAKEDIR}/CleanDir
include ${PRIMALMAKEDIR}/MakeHelp
include ${MAINMAKEDIR}/MakeGam
include ${MAINMAKEDIR}/Distribute
include ${MAINMAKEDIR}/Update
include ${MAINMAKEDIR}/DieDOSDie
include ${MAINMAKEDIR}/CopyDelDir
include ${MAINMAKEDIR}/TarZipDir

#################################################################################
#										#
#                      GAMMA CygWin Binary Distribution				#
#										#
# What follows has nothing to do with any builds of GAMMA either. This stuff is #
# included in this Makefile only so that GAMMA binary distributions for users	#
# that run GAMMA with Cygwin and GCC (MinGW option) can be readily made. Do not	#
# worry about this for any builds and or runs of GAMMA programs.		#
#										#
#################################################################################

GVER          = 4.1.0
BINDISTFILES  = ${INSTALLBINDIR}/${GSCRIPT} ${INSTALLBINDIR}/${DLL_NAME}
BINDISTFILES += ${INSTALLLIBDIR}/${DLL_EXP_LIB}
BINDISTFILES += $(addprefix $(GSRCDIR)/, $(ALLHFILES))
BINDISTFILES += $(addprefix $(GSRCDIR)/, gamma.h)
BINDISTFILES += $(addprefix $(GSRCDIR)/, GamGen.h)
BINDISTFILES += ./Readme.txt
BINTARFILE    = gamma-${GVER}.tar
BINZIPFILE    = ${BINTARFILE}.gz
BINFTPDIR     = ../../ftp/gamma-${GVER}/cygwin

./${BINTARFILE}: install
	@echo " * Generating Tar File Of Binaries ${BINTAR}"
	@rm -f ${BINTARFILE}
	@tar -cvf ${BINTARFILE} ${BINDISTFILES} 

./${BINZIPFILE}: ./${BINTARFILE}
	@echo " * Compressing Tar File Using GnuZip To Produce ${BINZIPFILE}"
	@rm -f ${BINZIPFILE}
	@gzip ${BINTARFILE}
	@echo " * Removing The Original Tar File"
	@rm -f ${BINTARFILE}

${BINFTPDIR}/${BINZIPFILE}: ./${BINZIPFILE}
	@echo " * Moving GAMMA Binaary Distrutionb To FTP Directory"
	@echo " * FTP Directory Is ${BINFTPDIR}"
	@mv -f ./${BINZIPFILE} ${BINFTPDIR}/${BINZIPFILE}

bindist: ${BINFTPDIR}/${BINZIPFILE}
	
#binfiles: force
#	tar -cvf gamma-${GVER}.tar ${BINDISTFILES} 
#	gzip gamma-${GVER}.tar
#	rm -f gamma-${GVER}.tar
#
#	@echo ${INSTALLBINDIR}/${GSCRIPT}
#	@echo ${INSTALLBINDIR}/${DLL_NAME}
#	@for i in ${ALLHFILES}; do \
#          echo ${GSRCDIR}/$$i; \
#  	done
#	install -d .${INSTALLBINDIR}
#	install -d .${GSRCDIR}
#	@for i in ${ALLMODDIRS}; do \
#          install -d .${GSRCDIR}/$$i; \
#  	done
#	install ${INSTALLBINDIR}/${GSCRIPT} .${INSTALLBINDIR}/${GSCRIPT}
#	install ${INSTALLBINDIR}/${DLL_NAME} .${INSTALLBINDIR}/${DLL_NAME}
#	install ${INSTALLLIBDIR}/${DLL_NAME} .${INSTALLLIBDIR}/${DLL_EXP_LIB}
#	@for i in ${ALLHFILES}; do \
#          install ${GSRCDIR}/$$i .${GSRCDIR}/$$i; \
#  	done
#	tar -cvf gamma-4.1.0.tar ${GSRCDIR} ${INSTALLBINDIR}
#	gzip  gamma-4.1.0.tar
#	rm -rf .${GSRCDIR}
#	rm -rf . ${INSTALLBINDIR}

force:
