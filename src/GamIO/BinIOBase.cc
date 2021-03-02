/* BinIOBase.cc *************************************************-*-c++-*-
**									**
**	                       G A M M A			 	**
**								 	**
**	Binary Input/Output Base                   Implementation	**
**								 	**
**	Copyright (c) 1999					 	**
**	Scott Smith 	          				 	**
**      Dr. Scott A. Smith                                              **
**      1800 E. Paul Dirac Drive                                        **
**      National High Magnetic Field Laboratory                         **
**      Tallahassee FL 32306 USA                                        **
**						 			**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**							 		**
**  Description							 	**
**								 	**
** These are simple functions which let users perform byte swapping	**
** of data when read from/written to machines that store their values	**
** in a different byte order.						**
**							 		**
*************************************************************************/

#ifndef   GBinIOB_cc_				// Is file already included?
#  define GBinIOB_cc_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#  endif

#include <GamIO/BinIOBase.h>			// Inlcude the header


#if defined(_MSC_VER) || defined(__MINGW32__)	// If using MSVC or MinGW                              // and if using MSVC++ then I guess
 #include <direct.h>				// Include mkdir function
#else
  #include <sys/stat.h>			// Include POSIX mkdir function
#endif

/*
union longchars
  {
  long longval;
  char chars[4];
  };

union shortchars
  {
  short shortval;
  char chars[2];
  };

union intchars
  {
  int intval;
  char chars[2];
  };
*/


// This patches the difference between GCC & MSVC++ mkdir functions. It
// appears as if maybe MSVC is more ANSI than GCC here? Oddly, MinGW is
// the same as MSVC but it is based on GCC 2.95.x... perhaps a later
// version that that I am using with Cygwin & Linux.

#if defined(__MINGW32__)
  int MakeADir(const std::string& dname, int no)
    { return mkdir(dname.c_str()); }
#elif defined(_MSC_VER)
  int MakeADir(const std::string& dname, int no)
    {
		return _mkdir(dname.c_str());
	  }
#else
  int MakeADir(const std::string& dname, int no)
    { return mkdir(dname.c_str(),  no); }
#endif


void Swap(int32_t& i)
  {
  longchars tmpj, tmpi;
  tmpi.longval=i;
  tmpj.chars[0] = tmpi.chars[3];
  tmpj.chars[1] = tmpi.chars[2];
  tmpj.chars[2] = tmpi.chars[1];
  tmpj.chars[3] = tmpi.chars[0];
  i = tmpj.longval;
  }

void Swap(short& i)
  {
  shortchars tmpj, tmpi;
  tmpi.shortval=i;
  tmpj.chars[0] = tmpi.chars[1];
  tmpj.chars[1] = tmpi.chars[0];
  i = tmpj.shortval;
  }

void Swap(double& d)
  {
  doublechars tmpj, tmpi;
  tmpi.doubleval=d;
  tmpj.chars[0] = tmpi.chars[7];
  tmpj.chars[1] = tmpi.chars[6];
  tmpj.chars[2] = tmpi.chars[5];
  tmpj.chars[3] = tmpi.chars[4];
  tmpj.chars[4] = tmpi.chars[3];
  tmpj.chars[5] = tmpi.chars[2];
  tmpj.chars[6] = tmpi.chars[1];
  tmpj.chars[7] = tmpi.chars[0];
  d = tmpj.doubleval;
  }

bool WeRBigEnd()
  {
  int x = 1;
  if(*(char *)&x == 1) return false;
  else                 return true;
  }

/* This function exists for changes in GCC from version 2.9x to 3.0.x.
   Earlier versions took an integer for the mode (the open mode was an
   enumeration in class ios) but newer version do not. Here I convert
   from integer to ios_base.  The conversion is

        in  = 1 (0x01)  out       = 2  (0x02)   ate      = 4   (0x03)
        app = 8 (0x08)  trunc     = 16 (0x10)   nocreate = 32  (0x20)
                        noreplace = 64 (0x40)   binary   = 128 (0x80)

   Other current C++ compilers do not seem to have a problems with this.
   In particular MSVC++ Versions 5 & 6. So, what I do here is to check
   if we are dealing with a 3.x.x version of GCC. If we are I will use
   the function below to convert to _Ios_Openmode. If not, the function
   will do nothing but return the same integer!                     */

#if (__GNUG__ > 2 )
  const std::ios_base::openmode Int2Mode(int mode)
    {
    return(std::ios_base::openmode(mode));
    }
#else
  int Int2Mode(int mode) { return mode; }
#endif

#endif							 // BinIoBase.cc
