/* ML5Hdr.cc ****************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A                           **
**                                                                      **
**	MATLAB Header (MAT Version 5)              Implementation	**
**                                                                      **
**      Scott Smith                                                     **
**      Dr. Scott A. Smith                                              **
**      1800 E. Paul Dirac Drive                                        **
**      National High Magnetic Field Laboratory                         **
**      Tallahassee FL 32306 USA                                        **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
** This MATLAB class provides functions for the input and output of     **
** headers from/to MAT version 5 files.  MATLAB is a commercial		**
** product and can be purchaced from The MathWorks, Inc.  This  module  **
** deals with binary MAT files version 5 circa 1999.                    **
**                                                                      **
*************************************************************************/

#ifndef _ML5HDR_CC_				// Is file already included?
#  define _ML5HDR_CC_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#  endif

#include <GamIO/ML5Hdr.h>			// Include the header file
#include <GamIO/BinIOBase.h>			// Include byte ordering stuff
#include <Basics/Gutils.h>                      // Include GAMMA errors
#include <string>				// Include libstdc++ strings
#include <fstream>				// Include libstdc++ file streams
#include <iostream>				// Include input output streams
#include <cstdio>				// Inlcude standard io stuff
#include <ctime>				// Include time/date functions

#if !defined(_MSC_VER) && !defined(__MINGW32__) // MSVC++ & MinGW do not use this
# include <sys/utsname.h>                       // Include platform information
#endif

using std::string;				// Using libstdc++ strings
using std::fstream;				// Using libstdc++ file streams
using std::ios;					// Using libstdc++ file type settings
using std::ostream;				// Using libstdc++ output streams

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                      MATLAB MAT 5 Header Error Handling
// ____________________________________________________________________________
 
/*              Input           ML5H    : MAT version 5 header
                                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
		Output		none    : Error message output
                                          Execution stopped (if fatal)       */  
 
void MatLab5Hdr::MLH5error(int eidx, int noret) const
  {                                                     
  string hdr("MATLAB MAT V5 Header");
  switch(eidx)
    {
    case 0: GAMMAerror(hdr, 0, noret); break;   // Program Aborting        (0)
    case 7: GAMMAerror(hdr, "End of File Reached Before Data Found", noret);
      break;                                     //                        (7)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }

 
void MatLab5Hdr::MLH5error(int eidx, const string& pname, int noret) const
  {
  string hdr("MATLAB MAT V5 Header");
  string msg;
  switch(eidx)
    {
    case 1:  GAMMAerror(hdr,   1, pname, noret); break;	// File Problems   (1)
    case 2:  GAMMAerror(hdr,   2, pname, noret); break;	// !Read SingPar   (2)
    default: GAMMAerror(hdr, -11, pname, noret); break;	// Unknown Error   (-1)
    }
  }
     
 
volatile void MatLab5Hdr::MLH5fatality(int eidx) const
  {
  MLH5error(eidx, 1);
  if(eidx) MLH5error(0);
  GAMMAfatal();					// Clean exit from program
  }
 
// ____________________________________________________________________________
// ii                   MATLAB MAT 5 Header Byte Ordering
// ____________________________________________________________________________
 
/* These functions deal with the header data storage endian flags.  Mostly they
   rely on the function WeRBigEnd found in the BinIOBase module which tests the
   current computer system for its "endianness".  In contrast, the header value
   BigEndIn is set when the header is read in from an external file so that we
   can know (in conjunction with WeRBigEnd) whether we need to do byte swapping
   on the values read.
 
	   Input		ML5H    : MAT version 5 header (this)
           Output               int     : True if computer big-endian
					  False if computer little-endian
				void    : If computer is big endian, MI is
					  set to {M,I} whereas if it is small
					  endian then MI is set to {I,M}.  It
					  is the order of the two bytes of MI
					  which indicate the "endianness" of
					  stored data in MAT Ver. 5 files.   */

void MatLab5Hdr::SetEndian()    { (WeRBigEnd())?SetBigEnd():SetLittleEnd(); }
void MatLab5Hdr::SetLittleEnd() { M='I'; I='M'; }
void MatLab5Hdr::SetBigEnd()    { M='M'; I='I'; }
int  MatLab5Hdr::BigEndian()    { if(M=='M') return 1; return 0; }


// ---------------------------------------------------------------------------- 
// ---------------------------- PUBLIC FUNCTIONS ------------------------------ 
// ---------------------------------------------------------------------------- 

// ____________________________________________________________________________
// A                      MATLAB MAT 5 Header Constructors
// ____________________________________________________________________________

/* These functions deal with constructing a MATLAB Version 5 MAT header in
   memory.  The default construdtor will create an appropriate header based
   on the current computer platform.                                         */

MatLab5Hdr::MatLab5Hdr()
  {
  string s("MATLAB 5.0 MAT-file");		// Begin string for textfiled

  s += string(", Platform: ");			// Add in platform info
#if defined(_MSC_VER) || defined(__SUNPRO_CC) || defined(__MINGW32__)
  s += string("Windows");
#else
  struct utsname name;			// Structure for info
  uname(&name);							// Fill up structure
  s += name.sysname;				// Add system name
#endif



  s += string(", Created on: ");		// Add in the time
  time_t now;					// Here is a time value
  time(&now);					// Set time value

#if defined(_MSC_VER)
	char timeStr[31];
  ctime_s(timeStr, 30, &now);				// Add as ctime format
	s += timeStr;
#else
  s += ctime(&now);		
#endif
  
	if(s[s.length()-1] == '\n')							// Remove return if ctime
  s = s.substr(0, s.length()-1);					// puts one on the end
  s += string(", Generated by GAMMA");		// Made with GAMMA

#ifdef _MSC_VER
	strcpy_s(TextField,124,s.c_str());
#else
	s.copy(TextField,124);  // Copy string s to TextField
#endif

  char empty(' ');
  for(int i=s.length(); i<124; i++)		// Zero rest of textfield
    TextField[i] = empty;
  Ver = 0x0100;					// Set the version
//  Ver = 1;					// Set the version
  SetEndian();					// Determine byte order
  size = 128;					// Total number of bytes
  BigEndIn = -1;				// No input endian type
  }

MatLab5Hdr::MatLab5Hdr(const MatLab5Hdr& MLH1)
  {
  for(int i=0; i<124; i++)			// Copy the entire textfield
    TextField[i] = MLH1.TextField[i];
  Ver = MLH1.Ver;				// Copy the MAT version
  M = MLH1.M;					// Copy the endian flags
  I = MLH1.I;
  BigEndIn = MLH1.BigEndIn;			// Copy input endian type
  }

// ____________________________________________________________________________
// B                    MATLAB MAT 5 Header Access Functions
// ____________________________________________________________________________

void MatLab5Hdr::skip(fstream& fp) { fp.seekg(size,ios::cur); }

// ____________________________________________________________________________
// C                MATLAB MAT 5 Header Binary Output Functions
// ____________________________________________________________________________
 
/* These are the functions which will output a header in MATLAB ".mat" binary
   format, version 5.  Note that in a valid MATLAB MAT file, this header will
   appear only once and reside at the beginning of the file.  Also note that
   this function should really output the default constructed header ONLY.
   But since this function should be private, used only by the friend class
   MatLabFile, we need only be carefull in that class.  If not there would be
   problems writing a header after reading it from a computer with a different
   byte storage order.

		Input           ML5H    : MAT version 5 header (this)
           			fp	: Output file stream
                                warn     : Warning output level
                                                0 = no warnings
                                                1 = warnings
                                               >1 = fatal warnings
		Output		void    : Output file stream is modified
					  by having the ML5H written to it
		Note			: No care is taken to insure the
					  header is written tothe top of the
					  file (as it should be in MATLAB)   */
 
int MatLab5Hdr::write(fstream& fp, int warn) const
  {
  for(int i=0; i<124 && fp.is_open(); i++)		// First output 124 char bytes
    fp.write(&(TextField[i]), sizeof(char));	// as the header text field
  if(fp.is_open())					// Continue writing as long
    {						// as output stream is OK
//    short V = Ver;
//    if(WeRBigEnd() == BigEndIn) Swap(V);	// Byte swap Ver if needed
    fp.write((char*)&Ver, sizeof(short));		// Write the two byte version
    fp.write((char*)&M,   sizeof(char));		// Write out endian indicator
    fp.write((char*)&I,   sizeof(char));		// (two bytes MI) 
    }
  if(!fp)					// Check that output stream
    {						// is still OK.  If not,
    if(warn)					// we'll output a warning
      {						// if needed and even die
      MLH5error(1, 1);				// if the error is fatal
      if(warn>1) MLH5fatality(11);
      else       MLH5error(11);
      }
    return 0;
    }
  return 1;
  } 

 
// ____________________________________________________________________________
// D                MATLAB MAT 5 Header Binary Input Functions
// ____________________________________________________________________________

/* These are the functions which will input a header in MATLAB ".mat" binary
   format, version 5.  Note that in a valid MATLAB V5 MAT file this header will
   appear only once and reside at the beginning of the file.  These functions
   can not alter the MAT file, they can only alter where the file pointer is.
   Upon a read, any data will be byte swapped as required by computer byte
   ordering differences.  Upon header output, any header text input will be
   altered accordingly.                                                      
 
                Input           ML5H    : MAT version 5 header (this)
                                fp      : Input file stream
                                warn     : Warning output level
                                                0 = no warnings
                                                1 = warnings
                                               >1 = fatal warnings
                Output          void    : ML5H is set from values found
                                          in the input stream
                Note                    : No care is taken to insure the
                                          header is read from the top of the
                                          file (as it should be in MATLAB)   */

int MatLab5Hdr::read(fstream& fp, int warn)
  {
  for(int i=0; i<124; i++)			// First input 124 char bytes
    { 						// as the header text field
    if(!fp || fp.eof())				// Insure filestream is open
      {						// Here if bad filestream
      if(warn)					// We'll output a warning
        {					// if needed and even die
        MLH5error(1, 1);			// if the error is fatal
        if(warn>1) MLH5fatality(11);
        else       MLH5error(11);
        }
      return 0;
      }
    fp.read(&(TextField[i]), sizeof(char));
    }
  fp.read((char*)&Ver, sizeof(short));			// Read the two byte version
  fp.read((char*)&M, sizeof(char));			// Read out endian indicator
  fp.read((char*)&I, sizeof(char));
  (M=='M')?BigEndIn=1:BigEndIn=0;		// Set input endian flag 
  if(int(WeRBigEnd()) == BigEndIn) Swap(Ver);	// Byte swap Ver if needed
  return 1;
  } 

// ____________________________________________________________________________
// E                  MATLAB MAT 5 Header ASCII Output Functions
// ____________________________________________________________________________
 

void MatLab5Hdr::print(ostream& ostr, int hpf) const

	// Input 		ML5H    : MAT version 5 header (this)
	// 			ostr	: An output stream
	//			hpf	: Header print flag
	// Output		ostr	: Output stream with MATLAB MAT file
	//				  header version 5 info placed into it
	// Note				: This information exists only at the
	//				  beginning of the MAT file

  {
  if(hpf) ostr << "\n\t\tHeader";	// Header Title
  ostr << "\n\t\t  Text:         ";	// 124 byte text field
  for(int i=0, nz=0, no=0; i<123; i++)	// Output it, leaving last
    {					// one off (a zero)
    if(TextField[i] == ',')		// A comma triggers new line
      {
      ostr << "\n\t\t               ";
      no = 0;
      }
    else if(TextField[i] == ' ')	// Skip zeros if needed
      nz++; 
    else
      {
      ostr << string(nz, ' ');
      ostr << TextField[i];
      no += nz + 1;			// Total output so far
      if(no > 55)
        {
        ostr << "\n\t\t               ";
        no = 0;
        }
      nz = 0;
      }
    }
short V = 0x0100;
if(Ver == V)
  ostr << "\n\t\t  Version:      " << 1;	// Version
else
  ostr << "\n\t\t  Version:      " << Ver;	// Version
  ostr << "\n\t\t  Endian:       " << M		// Endian Indicator
       << " " << I;
  if(M == 'I' && I == 'M') 
       ostr << " (stored little endian)";
  else ostr << " (stored big endian)";
  }


ostream& operator<< (ostream& ostr, const MatLab5Hdr& MLH5)

        // Input                ostr	: An output stream
	// 			ML5H    : MAT version 5 header
        // Output               ostr	: The output stream modified by
        //                                the header parameters

  { MLH5.print(ostr); return ostr; }


#endif							// ML5Hdr.cc
