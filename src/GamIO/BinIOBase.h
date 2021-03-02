/* BinIOBase.h **************************************************-*-c++-*-
**									**
**	                       G A M M A			 	**
**								 	**
**	Binary Input/Output Base                   Interface		**
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

#ifndef   GBinIOB_h_			// Is file already included?
#  define GBinIOB_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif


#ifdef _MSC_VER
#  if (_MSC_VER == 1500)
#    include <ms_stdint.h>      // #include "ms_stdint.h"
#  elif (_MSC_VER > 1910)		// first VS 2017 v15.x was 1910
#    include <stdint.h>
#  endif
#else
#  include <stdint.h>
#endif


#include <string>				// Include libstdc++ strings
#include <iostream>				// Know ios stuff

union longchars
  {
  int32_t longval;
  char chars[4];
  };

union doublechars
  {
  double doubleval;
  char chars[8];
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

void Swap(int32_t& i);
void Swap(double& d);
void Swap(short& i);
bool WeRBigEnd();


// This patches the difference between GCC & MSVC++ mkdir functions

int MakeADir(const std::string& dname, int no);

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
   const std::ios_base::openmode Int2Mode(int mode);
#else
   int Int2Mode(int mode);
#endif

#endif							 // BinIoBase.h
