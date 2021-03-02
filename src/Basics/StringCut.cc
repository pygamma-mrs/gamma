/* StringCut.cc *************************************************-*-c++-*-
**              	                                                **
**                             G A M M A                                **
**              	                                                **
**      String Manipulation Supplement 		      Implementation	**
**      Copyright (c) 1991, 1992, 1993					**
**      Tilo Levante, Scott Smith	                                **
**      Eidgenoessische Technische Hochschule                           **
**      Labor fur physikalische Chemie                                  **
**      8092 Zurich / Switzerland                                       **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**									**
**  Description								**
**									**
**  This module defines several string manipulation functions commoly	**
**  used in GAMMA to parse input parameters.  Mostly, these allow for	**
**  the clipping of strings in ways typically needed by GAMMA.  Since	**
**  the ANSI C++ library no longer contains Regex, this code doesn't	**
**  use it either (too bad).  Furthermore, it is based on the standard	**
**  "string" class rather than the older GNU libg++ "String" class.	**
**									**
*************************************************************************/

#ifndef   StringCut_cc_			// Is this file already included?
#  define StringCut_cc_ 1		// If no, then remember it 
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// If so, this is the implementation
#  endif

#include <GamGen.h>			// Include OS specific stuff
#include <Basics/StringCut.h>		// Include the header file
#include <string>			// Include C++ standard strings
#include <cstdio>			// Include function sprintf

#ifdef __BORLANDC__
  using std::sprintf;
#endif

using std::string;			// Using libstdc++ strings

        // Input                Sinp    : An input string
        // Output               cutSinp : String containing any of the
        //                                "blanks" removed from Sinp beginning
        // Note                         : String Sinp is altered here too
        //                                so that Sinp --> cutSinp + Sinp
        // Note                         : Mimicks using  RXwhite
 
string cutWhite(string& Sinp)
  {
  if(!Sinp.length()) return string("");	// Bail of no input string
  string RXwhite(" \t");			// Our mock Regex
  string cutSinp;				// String for white clipped out
  int i = Sinp.find_first_not_of(RXwhite);      // Index of 1st non-blank
  if(i<0)					// If the string is all blanks
    { 						// then set it to an empty
    cutSinp = Sinp; 				// string & return the original 
    Sinp = string("");
    }
  else
    {
    cutSinp = Sinp.substr(0,i);	 		// Copy initial blanks
    Sinp = Sinp.substr(i); 			// Clip initial blanks
    }
  return cutSinp;                               // Return blanks string
  }

        // Input                Sinp    : An input string
	//			xwhite  : Remove white space
        // Output               cutSinp : String containing any of the
        //                                "nonblanks" removed from Sinp start
        // Note                         : String Sinp is altered here too
        //                                so that Sinp --> (B)+cutSinp(B)+Sinp
 
string cutString(string& Sinp, bool xwhite)
  {
  if(!Sinp.length()) return string("");	// Bail of no input string
  if(xwhite) cutWhite(Sinp);			// Remove any white space
  string RXwhite(" \t");			// Our mock space Regex
  int i = Sinp.find_first_of(RXwhite);		// Index of 1st blank
  if(i < 0)					// If no white space
    {						// all of input is a single
    string cutSinp(Sinp);			// string
    Sinp = string("");
    return cutSinp;
    }
  string cutSinp = Sinp.substr(0,i); 	// Copy input to 1st blank
  Sinp = Sinp.substr(i);                        // Clip initial string
  if(xwhite) cutWhite(Sinp);			// Remove any white space
  return cutSinp;				// Return string
  }

        // Input                Sinp    : An input string
        // Output               cutSinp : String containing any of the
        //                                "( blanks" removed from Sinp beginning
        // Note                         : String Sinp is altered here too
        //                                so that Sinp --> cutSinp + Sinp
        // Note                         : Mimicks Regex RXParBlks("([ \t]*")
 
string cutParBlks(string& Sinp)
  {
  if(!Sinp.length()) return string("");	// Bail of no input string
  string cutSinp;				// This will be the return
  cutWhite(Sinp);                               // Remove initial blanks
  if(Sinp[0] == '(')                            // See if "(" is in Sinp
    {                                           // If so, put it into cutSinp
    cutSinp = string("(");			// and then clip it and any
    Sinp = Sinp.substr(1);                      // additional blanks
    cutSinp += cutWhite(Sinp);
    }
  return cutSinp;
  }

        // Input                Sinp    : An input string
	//			X       : String to parse out
	//			xwhite  : Remove white space
        // Output               cutSinp : String containing any of the
        //                                "blanksXblanks" from Sinp start
        // Note                         : String Sinp is altered here too
        //                                so that Sinp --> cutSinp + Sinp
	// Note				: Like Regex("[ \t]*X[ \t]*")
	// Notes			: Seems the compare function will
	//				  return FALSE when strings match!
 
string cutBlksXBlks(string& Sinp, const string& X, bool xwhite)
  {
  if(!Sinp.length()) return string("");	// Bail of no input string
  string cutSinp;				// This will be the return
  if(xwhite) cutWhite(Sinp);			// Remove initial blanks
#if defined(_MSC_VER) || defined(__SUNPRO_CC) || defined(__GNUC__) || defined(__BORLANDC__)
  if(!Sinp.compare(0,X.length(),X))		// See if X at start in Sinp
#else
 if(!Sinp.compare(X,0,X.length()))		// See if X at start in Sinp
#endif
    {                                           // If so, put it into cutSinp
    cutSinp = X;				// and then clip it and any
    Sinp = Sinp.substr(X.length());             // additional blanks
    if(xwhite) cutSinp += cutWhite(Sinp);	// Clip white at string end
    }
  return cutSinp;
  }

string cutDouble(string& Sinp, bool xwhite)
 
        // Input                Sinp    : An input string
	//			xwhite  : Remove white space
        // Output               cutSinp : String containing a double as
        //                                "blanksDoubleblanks" from Sinp start
        // Note                         : String Sinp is altered here too
        //                                  Sinp --> (B) + cutSinp + (B) + Sinp
	// Note				: Like RXdouble == 
	//				  "-?\\(\\([0-9]+\\.[0-9]*\\)\\|
        //                                 \\([0-9]+\\)\\|\\(\\.[0-9]+\\)\\)
        //                                 \\([eE][---+]?[0-9]+\\)?"

/* Here are some doubles that this should handle:

   1234    -1234     0.01234     -0.01234    .01234   -.01234 
   123e4   -123e4       123.e4      -123.e4    123e+4    123e-4
   123E4   -123E4	123.E4      -123.E4    123E+4    123E-4
   12.3e4  -12.3e4      12.3e-4     -12.3e-4
   12.3E4  -12.3E4      12.3E-4     -12.3E-4
   1234.   -1234.                                                             */

  {
  if(!Sinp.length()) return string("");	// Bail of no input string
  if(xwhite) cutWhite(Sinp);                    // Remove initial blanks
  string Sin(Sinp);			// Working copy of Sinp
  string cutSinp;				// Workout output strings

// All Doubles Will Start With Either An Integer Or A "-.", The Latter
// Only Occurs When The Value Is Set As -.#.  We'll First Clip Off The
// Integer (or Sign) From The Input String.

  string Sout = cutInt(Sin,0);		// First look for an integer
  string RXnum("0123456789");		// Set of integer constructs
  if(!Sout.length())				// If no integer was there
    {						// the next SHOULD be -.# 
    if(!Sin.find("-.") &&			// Check if we start with -.#
                 Sin.find_first_of(RXnum) == 2)
      {
      Sout = "-";				//   Set output for "-"
      Sin = Sin.substr(1);			//   Clip - from input
      }
    else if(Sin[0] == '.' &&			// It is possible just .#
              Sin.find_first_of(RXnum) == 1)	// used to set the double
      ;
    else					// No int nor -.# start so
      return cutSinp;				// we can't parse any double
    }

// Now That We Have Removed The Initial Integer (Or Sign), The Double
// Can Have a Decimal Or An Exponent.  We First Check To Make Sure The
// Next Character Is In {., e, E}, If Not The Integer is The Double

  int peE = Sin.find_first_of(".eE");		// Get index of .,e,orE
  if(peE != 0)					// If none exists we are done
    {						// as the double set by int
    Sinp = Sin;					//	Copy remaining string
    if(xwhite) cutWhite(Sinp);			//	Remove any white space
    return Sout;				//	Return the integer
    }

// Next We'll Check To See If There's a Decimal Value In The Double
// At This Point We Know Input Starts With Either ".", "e", or "E"

  if(Sin[0] == '.')				// Any decimal will begin as .
    {
    Sout += ".";				//   Add "." to output
    Sin = Sin.substr(1);			//   Clip . from input
    if(Sin[0] == '-')				//   Cannot have "-" next!
      {						//   if so, we are done
      Sinp = Sin;				//   Copy remaining string
      if(xwhite) cutWhite(Sinp);		//   Remove any white space
      return Sout;				//   Return the integer
      }
    Sout += cutInt(Sin,0);			// Cut any integer out
    }

// Next We'll Check To See If There's Any Exponential Value In The Double
// If There Is We Know Input Starts With Either "e" or "E"

  if(Sin.find_first_of("eE"))			// No exponent if not e or E
    {
    Sinp = Sin;					//   Copy remaining string
    if(xwhite) cutWhite(Sinp);			//   Remove any white space
    return Sout;				//   Return the integer
    }

  string etype = Sin.substr(0,1);		// Store "e" or "E"
  Sin = Sin.substr(1);				// Clip "e" or "E" from input
  if(!Sin.find_first_of("+-"))			// Exponent can have a sign
    {						// If a sign exists then
    etype += Sin.substr(0,1);			// add it to etype and remove
    Sin = Sin.substr(1);			// it form our input string
    }
  string sint = cutInt(Sin,0);			// Try for exponent value
  if(!sint.length())				// If no exponent, we were
    {						// done
    Sinp = etype + Sin;				//   Copy remaining string
    return Sout;				//   Return the integer
    }
  else						// We found the exponent
    {
    Sout += etype + sint; 
    Sinp = Sin;					//   Copy remaining string
    if(xwhite) cutWhite(Sinp);			//   Remove any white space
    cutSinp = Sout;				//   Set to return the integer
    }
  return cutSinp;                               // Return empty string
  }

        // Input                Sinp    : An input string
	//			xwhite  : Remove white space
        // Output               cutSinp : String containing an integer as
        //                                "blanksIntblanks" from Sinp start
        // Note                         : String Sinp altered here too so that
        //                                  Sinp --> (B) + cutSinp + (B) + Sinp
	// Note				: Like RXint == Regex("-?[0-9]+")
	// Note				: Integer -0 remains as -0

string cutInt(string& Sinp, bool xwhite)
  {
  if(!Sinp.length()) return string("");	// Bail of no input string
  if(xwhite) cutWhite(Sinp);			// Remove initial blanks
  string cutSinp("");			// Output string 
  string RXint("-0123456789");		// Our mock Regex
  string RXnum("0123456789");		// Just the integers
  if(!Sinp.find_first_of(RXint))		// Possible 1st char is int?
    {
    int ifi = Sinp.find_first_of(RXnum);	// See where first digit is
    if(ifi > 1) return cutSinp;			// if past col 2, can't be int 
    ifi = Sinp.find_first_not_of(RXnum, 1);	// Index of 1st non-int char
    if(ifi<0)					// Here if all part of int
      {						//      Everything is #'s!
      cutSinp = Sinp; 				//	Full string is int
      Sinp = string("");			//	Rest is empty
      }
    else if(ifi == 0 && Sinp[0]=='-')		// Possible we only have -
      return cutSinp;				// 	Leave all alone
    else					// Here if only part is int
      {
      cutSinp = Sinp.substr(0,ifi);		// Clip integer to cutSinp
      Sinp = Sinp.substr(ifi);			// Clip integer from Sinp
      if(xwhite) cutWhite(Sinp);		// Remove initial blanks
      }
    }
  return cutSinp;                               // Return empty string
  }

// ____________________________________________________________________________
// B                INTEGER TO STRING CONVERSION FUNCTIONS
// ____________________________________________________________________________

/* Some C++ compilers used to contain the dec function (which is defined in C)
   that will convert from integer to string.  However, such is not the case for
   Microsoft Visual C++. So I have removed all uses of dec in the GAMMA code
   and replaced them with the function Gdec. These are defined below.        */

string Gitoa(int i) { return Gdec(i); }
string Gdec(int i)
  { 
	char buffer[81];  

#ifdef _MSC_VER
	sprintf_s(buffer, 80, "%d", i); 
#else
	sprintf(buffer, "%d", i); 
#endif

	return string(buffer); 
	}
string Gdec2(long li)
  { 
	char buffer[161]; 

#ifdef _MSC_VER
	sprintf_s(buffer, 160, "%ld", li); 
#else
	sprintf(buffer, "%ld", li); 
#endif

	return string(buffer); 
	}
string Gdec(const string& fmt, int i)
  { 
	char buffer[81]; 

#ifdef _MSC_VER
	sprintf_s(buffer, 80, fmt.c_str(), i);
#else
	sprintf(buffer, fmt.c_str(), i); 
#endif

	return string(buffer); 
	}
string Gdec(int i, int digs)
  { string fmt = string("%") + Gdec(digs) + string("d"); return Gdec(fmt, i); }

// ____________________________________________________________________________
// C                NUMERIC TO STRING CONVERSION FUNCTIONS
// ____________________________________________________________________________

/* The standard C library had "form" functions to do the task of converting
   various numbers into strings. This was included in GNU's outdated libg++
   but is not defined (AFAIK) in libstdc++. It still seems to be viable in all
   GCC flavors but it is missing from MSVC++. So I have removed their use from
   the GAMMA code and replaced them with Gform functions defined below.      */

string Gform(const string& fmt, double d)
  { 
	char buffer[81]; 

#ifdef _MSC_VER
	sprintf_s(buffer, 80, fmt.c_str(), d); 
#else
	sprintf(buffer, fmt.c_str(), d); 
#endif

	return string(buffer); 
	}
string Gform(const string& fmt, int i)
  { 
	char buffer[81]; 

#ifdef _MSC_VER
	sprintf_s(buffer, 80, fmt.c_str(), i); 
#else
	sprintf(buffer, fmt.c_str(), i); 
#endif

	return string(buffer); 
	}

// ____________________________________________________________________________
// D                        LINE READING FUNCTIONS
// ____________________________________________________________________________

        // Input                s	: An input stream
        //			Sout	: An output string
        // Output               TF      : True if line read
        // Note                         : Older GNU libg++ had readline
        //                                function so this replaces it

bool Greadline(std::ifstream& s, string& Sout, char terminator)
  {
  char buf[81];				// Need a character buffer
  s.getline(buf, 81, terminator);   	// Read next (2nd) line
  if(s.eof()) return false;		// Exit if no data
  Sout = string(buf);		// Set buffer to string
  return true;				// Flag we were successful
  }

/*
int Greadline(std::istream& s, string& Sout, char terminator)

        // Input                s	: An input stream
	//			Sout	: An output string
        // Output               TF      : True if line read
        // Note                         : Older GNU libg++ had readline
	//				  function so this replaces it

  {
  char buf[200];				// Need a character buffer
  s.getline(buf, 200, terminator);		// Read next (2nd) line
  if(s.eof()) return 0;				//   Exit if no data
  Sout = string(buf);			// Set buffer to string
  return 1;
  }
*/

/*
int readline(istream& s, string& Sout, char terminator, int discard)
  {
  if(!s.ipfx(0)) return 0;
  int ch;
  int i = 0;
  Sout.rep = Sresize(Sout.rep, 80);
  register streambuf *sb = s.rdbuf();
  while((ch = sb->sbumpc()) != EOF)
    { 
    if (ch != terminator || !discard)
    {
      if (i >= Sout.rep->sz - 1)
        Sout.rep = Sresize(Sout.rep, i+1);
      Sout.rep->s[i++] = ch;
    }
    if (ch == terminator)
      break;
  }
  Sout.rep->s[i] = 0;
  Sout.rep->len = i;
  if (ch == EOF) s.clear(ios::eofbit|s.rdstate());
  return i;
  }   

        // Input                Sinp    : An input string
        // Output               i 	: Integer value of string
        // Note                         : I'm sick of C++ library changes
	//				  where atoi, dec, etc seem to disappear
	//				  so here is one possible conversion
int string2int(const string& Sinp)
  {
  int i;
  istringstream(Sinp) >> i;
  return i;
  }
*/

// ____________________________________________________________________________
// E                         STRING ALIGNMENT FUNCTIONS
// ____________________________________________________________________________

string CenterString(const string& str, int width)
  {
  int sl = str.length();
  int len = (width-sl)/2;
  if(len>0) return string(len, ' ') 
                 + str 
                 + string(width-len-sl, ' ');
  else      return str;
  }

std::ostream& CenterString(std::ostream& ostr, const string& str, int width)
  {
  int sl = str.length();
  int len = (width-sl)/2;
  if(len>0) 
    ostr << string(len, ' ')
         << str 
         << string(width-len-sl, ' ');
  else
    ostr << str
         << string(width-sl, ' ');
  return ostr;
  }


#endif						           // StringCut.cc
