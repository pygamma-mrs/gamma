/* GamGen.h *****************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      GAMMA Generics				Interface		**
**                                                                      **
**      Copyright (c) 2002						**
**      Scott A. Smith 							**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**      email: gamma@magnet.fsu.edu					**
**      www: http://gamma.magnet.fsu.edu				**
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
**  This file contains definitions that can be system specific. In the	**
**  future this header may actually be altered by the configuration.	**
**  For now it is static.						**
**                                                                      **
*************************************************************************/
 
#ifndef GAMGENERIC_H__			// Is this header known?
#define GAMGENERIC_H__			// If not then define & include

// ____________________________________________________________________________ 
// A                     Microsoft MSVC++ Specifics
// ____________________________________________________________________________

/* The macro _MSC_VER is set by Visual C++ during compilations. It is this
   macro that is mainly used to test for use of MSVC++ and to set up any code
   specifics needed for it. An additional macro was set up by myself, MSVCDLL,
   to handle function declspec decorations needed to build DLLs in windows for
   some systems. This is independent of _MSC_VER, the former is used in all
   MSVC++ compilations whereas the latter is used anytime declspec stuff is
   needed when building DLLs in a Windows environment. Since we minimize the
   amount of declspec usage in GAMMA code, it is ususally only when linking
   to any GAMMA DLL that this is used. DLLs built using Cygwin or MinGW ports
   of GCC to windows do not need the declspec decorations either. Don't know
   about Borland, Intel, and any others.
*/
 
// ----------------------------------------------------------------------------
//                The _declspec Stuff For MSVC++ DLL Builds
// ----------------------------------------------------------------------------

/* This variable, MSVCDLL, will be used in GAMMA to handle any __declspec
   declarations required by Microsoft when building DLLs using Visual C++.
   This Microsoftism is to qualify everything (functions, constants, etc) with
   the __declspec(dllexport) when building a DLL and with __declspec(dllimport)
   when using a DLL. Of course, this adds a lot of garbage to the code headers
   which is Microsoft specific. If you are not building a DLL on a Windows
   PC (mainly with MSVC++) then these declspec quailfiers are not required.
   
   GAMMA attempts to avoid trashing its headers with these goofy _declspec
   declarations everywhere. As a consequence, things can be quite painful when
   it comes to DLL building! Here is what is done.

   Microsoft sez that __declspec(dllexport) is not required if the entity one
   is exporting is listed in a .DEF file. The .DEF file is specified when the
   DLL is built & this provides information used in making the static library
   (foo.dll.lib) associated with the DLL (foo.dll).So, instead of making use
   of __declspec(dllexport) in GAMMA, there exists .DEF files that contain all
   of the functions and constants we want exported in GAMMA DLL(s). This is a
   big pain too, because they are partially made by hand. But so be it if we
   can avoid MS code specific additions.

   Use of _declspec(dllimport) is a different story. Although use of a .DEF
   file will completely eliminate the need for _declspec(dllexport), Microsoft
   demands that _declspec(dllimport) be used for static constants. In fact,
   MS recommends that it be used for functions as well, but that is not
   mandatory. To minimize MS specific code, we use _declspec(dllimport) as
   needed on static constants only. The definition of MSVCDLL will handle
   both cases, but is only used in the declaration of all exported contants,
   It is only needed when linking to an MSVC++ built DLL of GAMMA. Hence its
   left empty for the vast majority of GAMMA compilations, including DLLs made
   with Cygwin/GCC, MinGW, and even the PyGAMMA DLL under MSVC++. That is, do
   NOT define USEDLL under any builds except those that are linking to the
   GAMMA DLL under MSVC++. Note that the exporting of these constants for
   MSVC++ is done through use of a .DEF file.                               */

#ifndef MSVCDLL
# if defined(MAKEDLL)
#   define MSVCDLL __declspec(dllexport)
# elif defined(USEDLL)
#   define MSVCDLL __declspec(dllimport)
# else
#   define MSVCDLL
# endif
#endif

#ifndef MSVCDLC
# if defined(MAKEDLL)  && defined(_MSC_VER)
#   define MSVCDLC __declspec(dllexport)
# elif defined(USEDLL) && defined(_MSC_VER)
#   define MSVCDLC __declspec(dllimport)
# else
#   define MSVCDLC
# endif
#endif
 
// ----------------------------------------------------------------------------
//             General Declarations For Any MSVC++ Builds
// ----------------------------------------------------------------------------

//   Note that _MSC_VER is automatically set by the Visual C++ compiler!

// One of the big pains about MSVC++ (as of Version 6) is that long function
// names producing warning messages. I hear this is fixed in later versions,
// but with the use of templated classes one cannot help but get these damn
// warnings. Issuing the pragma below is supposed to take care of that....
// It gets a lot of them, but not all.

#if defined(_MSC_VER)				// If using MSVC++ then
#  pragma warning (disable : 4786)		// kill STL namelength warnings
#endif

// ____________________________________________________________________________ 
// B                         Gnu GCC Compiler
// ____________________________________________________________________________
 
// ----------------------------------------------------------------------------
//                       GCC Version Evaluation
// ----------------------------------------------------------------------------

/* There are several problems one might deem as GCC growing pains. These are
 * troubles that have arisen as GCC becomes more and more ANSI stardardized.
 * Consequently, it is nice to know the GCC version number during GAMMA builds
 * so that such version problems can be handled cleanly.  These next lines use
 * the * GCC defined macros __GNUC__, __GNUC_MINOR__, & __GNUC_PATCHLEVEL__.
 * These are self defining and would be 3, 2, 0 respectively for GCC 3.2.0   */

 
/* One problem is in the std::string.compare function. In GCC 2.9.x it seems
   that the argument order is switched over that in GCC 3.2.x. Other compilers
   (MSVC++, Sun Pro,...) have theirs defined the same as 3.2.x, so I assume
   that it is the standard. I am not sure about 3.0.x and 3.1.x (I hear they
   are bad anyway). The only time I have thus far dealt by 3.1.x is on OSX.2
   and that made for all kinds of troubles, mostly from the way Apple did
   their port.                                                               */
 
// ----------------------------------------------------------------------------
//                GCC Pragma Interface and Implementation Usage
// ----------------------------------------------------------------------------

/* Another problem is that there are problems in GAMMA code due to use of
   #pragma interface and #pragma implementation in GCC 3.x.x. Perhaps the GCC
   documentation is out of date on this subject, but if pragmas are used then
   there are multiple "undefined reference" errors when any programs are linked
   to the compiled GAMMA library (or libraries). So, I put an ifdef on all of
   the GNuish pragma (implementation, interface) statments so that they are
   used only if 1.) __GNUG__ is defined and 2.) GAMPRAGMA is defined. Here I
   will set GAMPRAGMA only if we are using GCC version 2.x.x.                */


#if __GNUC__ >= 4
#include <cstdlib>
#endif


// ____________________________________________________________________________ 
// C                   Sun SunPro C++ Compiler (SPARC)
// ____________________________________________________________________________

/*
#if defined(__SUNPRO_CC)			// Section for Sun C++ compiler
#endif						// End of Sun C++ block
*/

// ____________________________________________________________________________ 
// D                      Borland C++ Compiler
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                        Compiler Warnings Removal
// ----------------------------------------------------------------------------

/* The free command line tools (C++ compiler) provided by Borland produces way
   too many warnings for me, at least for current GAMMA building purposes. So
   these pragmas will remove many of the most annoying ones for now. I might
   remove some of these later and patch the GAMMA sources if there is ever
   time for it, but for now this is easier....                               */

#if defined(__BORLANDC__)			// If using Borland C++ then
#  pragma warn -8057				//   kill unused variable warnings
#  pragma warn -8004				//   kill unused var assignment warn
#  pragma warn -8066				//   kill unreachable code warnings
#  pragma warn -8027				// * kill inline expansion warnings
#  pragma warn -8022				// * hide virtual function warnings
#  pragma warn -8070				// * kill function return warnings
#endif						// End Borland C++ specifics

// ----------------------------------------------------------------------------
//                              Bessel Functions
// ----------------------------------------------------------------------------

/* Apparently the Borland compiler does not provide Bessel functions. GAMMA
   uses one of these, j0, someplace. This function should be in the math
   library but it isn't... forcing me to define one here.                    */

#if defined(__BORLANDC__)			// If using Borland C++ then
						// we must define j0
#include <math.h>
double j0(double x)
  { 
  double ax,z; 
  double xx,y,ans,ans1,ans2;  
  if ((ax=fabs(x)) < 8.0)
    { 
	y=x*x; 
	ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7 
	    +y*(-11214424.18+y*(77392.33017+y*(-184.9052456))))); 
	ans2=57568490411.0+y*(1029532985.0+y*(9494680.718 
            +y*(59272.64853+y*(267.8532712+y*1.0)))); 
	ans=ans1/ans2; 
    }
  else
    { 
	z=8.0/ax; 
	y=z*z; 
	xx=ax-0.785398164; 
	ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4 
	    +y*(-0.2073370639e-5+y*0.2093887211e-6))); 
	ans2 = -0.1562499995e-1+y*(0.1430488765e-3 
	    +y*(-0.6911147651e-5+y*(0.7621095161e-6 
	    -y*0.934935152e-7))); 
	ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2); 
    } 
  return ans; 
  } 
#endif						// End Borland C++ specifics

// ____________________________________________________________________________ 
// Z                          All Platforms
// ____________________________________________________________________________

/*       defines things that overcome problems with specific platforms       */

// ----------------------------------------------------------------------------
//                   Generic Max and Min Functions
// ----------------------------------------------------------------------------

/* Normally there is a generic max and min function specified in the ANSI C++
   header <cmath>. This turns out not to be true in MSVC++ (at least as of
   version 6.0 which defines the functions __max and __min instead in <stdlib.h>
   in order to bypass naming problems with other similar functions it provides.
   Also, it is claimed (but I have not verified this, richtig Matthias?) that
   the Sun Pro C++ compiler on SPARC suffers the same problem. So I just declare
   my own but call them gmax and gmin to avoid all of the potential cross
   platform problems with max and min... and these are so simple...           */

#ifndef __GMAXMIN				//   If we have no max/min
#  define __GMAXMIN 1				//      flag that we do now
#  define gmin(a,b) (((a)<(b)) ? (a):(b))	//      define max and min
#  define gmax(a,b) (((a)>(b)) ? (a):(b))
#endif						//   End of max/min block

#ifndef HUGE					// Sometimes HUGE isn't known
  #include <limits>
  #define HUGE std::numeric_limits<double>::infinity();

//#  define HUGE HUGE_VAL				// but HUGE_VAL is.  We'll make
#endif

#endif								// GamGen.h

