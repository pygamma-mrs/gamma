/* StringCut.h **************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      String Manipulation Supplement 		      Interface		**
**      Copyright (c) 1991, 1992, 1993					**
**      Tilo Levante, Scott Smith					**
**      Eidgenoessische Technische Hochschule				**
**      Labor fur physikalische Chemie					**
**      8092 Zurich / Switzerland					**
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description								**
**                                                                      **
**  This module defines several string manipulation functions commonly	**
**  used in GAMMA to parse input parameters.  Mostly, these allow for	**
**  the clipping of strings in ways typically needed by GAMMA.  Since	**
**  the ANSI C++ library no longer contains Regex, this code doesn't	**
**  use it either (too bad).  Furthermore, it is based on the standard  **
**  "string" class rather than the older GNU libg++ "String" class.	**
**                                                                      **
*************************************************************************/

#ifndef   StringCut_h_			// Is this file already included?
#  define StringCut_h_ 1		// If no, then remember it 
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// If so, this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <string>			// Include libstdc++ strings
#include <fstream>			// Include file streams (ifstream)

// ____________________________________________________________________________
// A                        STRING PARSING FUNCTIONS
// ____________________________________________________________________________

/* It would be nice to get regex working clean in C++ but with the new ANSI
   compilers I don't know how it works anymore.  So, I had to bite the bullet
   and define some primitive string parsing routines for GAMMA. One day I'll
   revisit regex and replace all of this stuff with something more general. 
	
	// Input  	Sinp      : An input string
	// Output	cutSinp   : Input string modified by being parsed
	//			    often just having beginning "blanks"
	//               xwhite   : Flag to signal white space removal
	// Note			  : String Sinp is altered here too
	//	                    so that Sinp --> cutSinp + Sinp

  Function                         Constructed Array
 ----------    -----------------------------------------------------------------
  cutWhite	Removes white space from input string start { RXwhite }
  cutString	Remove first isolated from string start input string
 cutParBlks	Removes "( blanks" from string start { Regex("([ \t]*") }
cutBlksXBlks	Removes "blanksXblanks" at string start {Regex("[ \t]*X[ \t]*")}
  cutDouble     Removes a double from string start { Regex RXdouble }
   cutInt       Removes an integer from string start { Regex RXint }          */

MSVCDLL std::string cutWhite(std::string&     Sinp);
MSVCDLL std::string cutString(std::string&    Sinp,                       bool xwhite=true);
MSVCDLL std::string cutParBlks(std::string&   Sinp); 
MSVCDLL std::string cutBlksXBlks(std::string& Sinp, const std::string& X, bool xwhite=true);
MSVCDLL std::string cutDouble(std::string&    Sinp,                       bool xwhite=true);
MSVCDLL std::string cutInt(std::string&       Sinp,                       bool xwhite=true);

// ____________________________________________________________________________
// B                INTEGER TO STRING CONVERSION FUNCTIONS
// ____________________________________________________________________________

/* Some C++ compilers used to contain the dec function (which is defined in C)
   that will convert from integer to string.  However, such is not the case for
   Microsoft Visual C++. So I have removed all uses of dec in the GAMMA code
   and replaced them with the function Gdec. These are defined below.        */

MSVCDLL std::string Gitoa(int i);
MSVCDLL std::string Gdec(int i);
MSVCDLL std::string Gdec2(long li);
MSVCDLL std::string Gdec(const std::string& fmt, int i);
MSVCDLL std::string Gdec(int i, int digs);
 
        // Input                i       : An integer
        // Output               Si      : String representation of i
        // Note                         : Older GNU libg++ had dec function
        //                                so this replaces it
 
// ____________________________________________________________________________
// C                NUMERIC TO STRING CONVERSION FUNCTIONS
// ____________________________________________________________________________

/* The standard C library had "form" functions to do the task of converting
   various numbers into strings. This was included in GNU's outdated libg++
   but is not defined (AFAIK) in libstdc++. It still seems to be viable in all
   GCC flavors but it is missing from MSVC++. So I have removed their use from
   the GAMMA code and replaced them with Gform functions defined below.      */

MSVCDLL std::string Gform(const std::string& fmt, double d);
MSVCDLL std::string Gform(const std::string& fmt, int i);

// ____________________________________________________________________________
// D                        LINE READING FUNCTIONS
// ____________________________________________________________________________

MSVCDLL bool Greadline(std::ifstream& s, std::string& Sout, char terminator='\n');

        // Input		s	: An input stream
	//			Sout	: An output string
        // Output		TF	: True if line read
        // Note				: Older GNU libg++ had readline
	//				  function so this replaces it

// ____________________________________________________________________________
// E                         STRING ALIGNMENT FUNCTIONS
// ____________________________________________________________________________

MSVCDLL std::string   CenterString(const std::string& str, int width=80);
MSVCDLL std::ostream& CenterString(std::ostream& ostr, const std::string& str, int width=80);

#endif						                      // StringCut.h
