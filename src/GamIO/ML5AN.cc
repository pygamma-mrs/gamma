/* ML5AN.cc *****************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A                           **
**                                                                      **
**     MATLAB Array Name (MAT Version 5)        Implementation		**
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
** This MATLAB class provides functions for the input and output of a   **
** MATLAB array name sub-element in MAT version 5 files.  MATLAB 	**
** is a commercial product and can be purchaced from The MathWorks, Inc.**
** This module deals with binary MAT files version 5 circa 1999.      	**
**                                                                      **
*************************************************************************/

#ifndef _ML5AN_CC_				// Is file already included?
#define _ML5AN_CC_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)					// Using the GNU compiler?
#    pragma implementation				// This is the implementation
#endif

#include <GamIO/ML5AN.h>			// Include the header file
#include <GamIO/ML5SubE.h>			// Include base class header
#include <Basics/Gutils.h>                      // Include GAMMA errors
#include <string>				// Include libstdc++ strings
#include <stdio.h>


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i          MATLAB MAT 5 Array Name Sub-Element Error Handling
// ____________________________________________________________________________

/*              Input           MLAN5   : MAT version 5 sub-element
                                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */  

void MatLab5AN::MLAN5error(int eidx, int noret)
  {                                                     
  std::string hdr("MATLAB MAT V5 Array Name Sub-Element");
  switch(eidx)
    {
    case 0: GAMMAerror(hdr, 0, noret); break;   // Program Aborting        (0)
    case 1: GAMMAerror(hdr, "End of File Reached Before Data Found", noret);
      break;                                     //                        (1)
    case 10: GAMMAerror(hdr, "Sorry, Cannot Read Myself", noret); 	// (10)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }

 
void MatLab5AN::MLAN5error(int eidx, const std::string& pname, int noret)
  {
  std::string hdr("MATLAB MAT V5 Array Name Sub-Element");
  std::string msg;
  switch(eidx)
    {
    case 1:  GAMMAerror(hdr,   1, pname, noret); break;	// File Problems   (1)
    case 2:  GAMMAerror(hdr,   2, pname, noret); break;	// !Read SingPar   (2)
    default: GAMMAerror(hdr, -11, pname, noret); break;	// Unknown Error   (-1)
    }
  }
     
 
volatile void MatLab5AN::MLAN5fatality(int eidx)
  {
  MLAN5error(eidx, 1);
  if(eidx) MLAN5error(0);
  GAMMAfatal();					// Clean exit from program
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------


// ____________________________________________________________________________
// A                  MATLAB MAT 5 Array Name Sub-Element Constructors
// ____________________________________________________________________________


MatLab5AN::MatLab5AN() : MatLab5SE() { NC = 0; MxName = ""; }
MatLab5AN::MatLab5AN(const std::string& name) : MatLab5SE()
  {
  NC = name.length();			// Store the number of chars
  MxName = name;			// Store the array name
  MLTag = MatLab5Tag(1,NC,(NC<5)?1:0);
  }

MatLab5AN::~MatLab5AN() { }
 

// ____________________________________________________________________________
// B            MATLAB MAT 5 Array Name Sub-Element Access Functions
// ____________________________________________________________________________
   
/* These functions return the array name of the sub-element.  Without any
   arguments this will return the currently listed name.  If given a file
   stream, it is assumed that the file points to the ed to the */

std::string MatLab5AN::Name() { return MxName; }

// ____________________________________________________________________________
// C            MATLAB MAT 5 Array Name Sub-Element Binary Output Functions
// ____________________________________________________________________________
   
/* These are the functions which will output this sub-element in MATLAB ".mat"
   binary format, version 5.  MAT sub-elements are output as a pair {tag,data}.
   In this case, the output can be compressed if the array name is four chars
   or less. This we will do if it is possible.

   Remember, the MATLAB MAT version 5 tag for this sub-element contains
   8-bytes with two important settings: 1.) data type & 2.) # of bytes.
   The data type is placed in the first 4 bytes and the number of bytes
   placed into the last 4 bytes for non-compressed data elements.  For
   compressed data elements the 1st 2 bytes contains the # of bytes,
   the 2nd 2 bytes contain the data type, and the last 4 bytes are the
   characters for the array name.  In this case the data type is written
   as miINT8 (i.e. value 1)
 
                Input           ML5AN	: MAT version 5 array name SE (this)
                                fp      : Output file stream
				aname   : Array name
                Output          void    : Output file stream is modified
                                          by having the ML5AN written to it
                Note                    : No care is taken to insure the
                                          tag is written.                    */
 
                                 
int MatLab5AN::write(std::fstream& fp) const
  {
  MLTag.write(fp);				// First write tag to file
  char c;					// Char for name character
  int i;
  for(i=0; i<NC; i++)				// Output the characters of
    {						// the array name
    c = MxName[i];
    fp.write(&c, sizeof(char));
    }
  int left = NC%8;				// Amount to fill 8 byte 
  if(NC <= 4) left = NC%4;			// boundary of sub-element
  if(left)					// This varies whether the
    {						// tag is compressed or not
    left = 8-left;				// If not compressed span 8, if
    if(NC <= 4) left = left-4;			// compressed span 4 (4 in tag)
    c = ' ';					// We'll use empty to pad till
    for(i=0; i<left; i++)			// we reach the boundary
      fp.write(&c, sizeof(char));
    }
  return 1;
  }


int MatLab5AN::write(std::fstream& fp, const std::string& aname) const
  {
  int len = aname.length();			// Length of string
  if(len < 1) return 0;				// Nothing if no name
  if(len < 5)					// If short name write
    { 						// compressed sub-element tag
    MatLab5Tag T(1, len, 1);			// Construct appropriate tag
    T.write(fp);				// Write tag to file
    }
  else						// If name longer than 4 chars
    {						// write normal sub-element
    int32_t TB = 1;				// MATLAB data type miDOUBLE
    fp.write((char*)&TB, sizeof(int32_t));		// Output the data type  (TAG)
    TB = len;					// Set number of bytes
    fp.write((char*)&TB, sizeof(int32_t));		// Output # of characters(TAG)
    }						// (assumes 1 char = 1 byte)
  char c;
  int i;
  for(i=0; i<len; i++)				// Output the characters of
    {						// the array name
    c = aname[i];
    fp.write(&c, sizeof(char));
    }
  int left = len%8;				// Amount to fill 8 byte 
  if(len < 4) left = len%4;			// boundary of sub-element
  if(left)
    {
    left = 8-left;
    if(len < 4) left = left-4;
    c = ' ';					// (we'll use empty) till we
    for(i=0; i<left; i++)			// reach the boundary
      fp.write(&c, sizeof(char));
    }
  return 1;
  }

/* These are the functions which will return the number of bytes that are  
   written upon output this sub-element in MATLAB ".mat" binary format, V5   */ 
  
int MatLab5AN::Size(const std::string& name) const

  { 
  int NB = MLTag.Size();			// Set for tag length (8 bytes)
  int nl = name.length();			// Number of characters
  if(nl>4)					// If >4 chars, NOT compressed
    {						// The characters are in the
    NB += nl;					// data portion which must
    int left = nl%8;				// align on an 8-byte boundary 
    if(left) NB += 8-left;
    }	
  return NB;
  }



// ____________________________________________________________________________
// D          MATLAB MAT 5 Array Name Sub-Element Binary Input Functions
// ____________________________________________________________________________

/* These are the functions which will input an Array Name sub-element in MATLAB
   ".mat" binary format, version 5.  Note that in a valid MATLAB MAT file, a 
   "Array Name Sub-Element" is a combination of Tag + Data
 
                Input           ML5AN	: MAT version 5 array name (this)
                                fp      : Input file stream
				bigend	: Data storage endian type
                                warn	: Warning output level
                                            0 = no warnings
                                            1 = warnings
                                            >1 = fatal warnings
                Output          void    : ML5AN is set from values found
                                          in the input stream
                Note                    : No care is taken to insure the
                                          array name is read from the top of a
                                          data element (as it should be in
					  MATLAB)                            */

int MatLab5AN::read(std::fstream& fp, int bigend, int warn)
  {
  if(!MatLab5SE::read(fp, bigend, (warn)?0:1))	// Read sub-element
    {                                           // If this has problems, then
    if(warn)                                    // we must deal with them
      {                                         // appropriately now
      if(warn==1) MLAN5error(10, 1);
      else        MLAN5fatality(10);
      }  
    return 0;
    }
  NC = MLTag.Bytes();				// Set number of characters
  MxName = "";					// Begin with empty name
  if(MLTag.Compressed)
    for(int i=0; i<NC; i++)			// Loop over name characters
      MxName += MLTag.EightBytes[i+4];		// and set the string for 'em
  else
    for(int i=0; i<NC; i++)			// Loop over name characters
      MxName += MLData[i];			// and set the string for 'em
  return 1;
  }

 
// ____________________________________________________________________________
// E                MATLAB MAT 5 Array Name Sub-Element ASCII Output Functions
// ____________________________________________________________________________

                                                                                
void MatLab5AN::print(std::ostream& ostr) const

        // Input                ML5AN	: MAT version 5 sub-element (this)
        //                      ostr    : An output stream
        // Output               ostr    : Output stream with MATLAB MAT file
        //                                sub-element version 5 info placed
	//				  into it
        // Note                         : This information exists at the
        //                                beginning of each sub-element

  {
  ostr << "\n\t\tSE Array Name   ";		// Print the header info
  MLTag.print(ostr, 0);				// Output 8-byte tag
  ostr << "\n\t\t  Name Chars:   " << NC; 	// # of characters in name
  ostr << "\n\t\t  Array Name:   " << MxName; 	// Matrix name
  }
 
std::ostream& operator<< (std::ostream& ostr, const MatLab5AN& ML5AN)
 
        // Input                ostr    : An output stream
        //                      ML5AN	: MAT version 5 sub-element
        // Output               ostr    : The output stream modified by
        //                                the sub-element parameters

  { ML5AN.print(ostr); return ostr; }


#endif							// ML5AN.cc
