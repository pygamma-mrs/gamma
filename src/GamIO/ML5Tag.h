/* ML5Tag.h *****************************************************-*-c++-*-
**									**
**	                            G A M M A			 	**
**								 	**
**	MATLAB Tag (MAT Version 5)                     Interface	**
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
** This MATLAB class provides functions for the input and output of	**
** tags associated with a data element in a MATLAB MAT version 5 file.  **
** MATLAB is a commercial product and can be purchaced from The 	**
** MathWorks, Inc. This module deals with binary MAT files version 5 	**
** circa 1999. 								**
**							 		**
*************************************************************************/

#ifndef _GML5TAG_h_				// Is file already included?
#define _GML5TAG_h_				// If no, then remember it

#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// This is the interface 
#  endif

#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <string>				// Know libstdc++ strings
#include <fstream>				// Know libstdc++ filestreams
#include <GamIO/BinIOBase.h>			// Know binary byte stuff
 
class MatLab5Tag
  {
  char EightBytes[8];				// 8 byte field
  int  DType;					// Data type
  int  NBytes;					// Number of bytes
  int  Compressed;				// Flag if compressed input

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                      MATLAB MAT 5 Tag Error Handling
// ____________________________________________________________________________

/* 		Input 		ML5T    : MAT version 5 tag	
				eidx    : Error index
               			pname   : string in message 
				noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */
void MLT5error(int eidx, int noret=0);
void MLT5error(int eidx, const std::string& pname, int noret=0); 
volatile void MLT5fatality(int eidx);

int IsBigEndian() const
  {
  int x = 1;
  if(*(char *)&x == 1) return 0;
  else                 return 1;
  }

public:

friend class MatLab5SE;				// Data Sub-Element Full Access
friend class MatLab5AF;				// Array Flags Full Access
friend class MatLab5DA;				// Dimensions Array Full Access
friend class MatLab5AN;				// Array Name Full Access
friend class MatLab5Re;				// Real Array Full Access
friend class MatLab5Im;				// Imag Array Full Access
friend class MatLab5DE;				// Data Element Full Access
 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                    MATLAB MAT 5 Tag Constructors
// ____________________________________________________________________________

 
MatLab5Tag();
MatLab5Tag(int TY, int NB, int CMP=0);
MatLab5Tag(const MatLab5Tag& T);
void operator=(const MatLab5Tag& T);
 
// ____________________________________________________________________________
// B                    MATLAB MAT 5 Tag Access Functions
// ____________________________________________________________________________


int Type(int LE=0) const;

        // Input                ML5T    : MAT version 5 tag (this)
        //                      LE      : Little Endian flag (default = big)
        // Output               type    : An integer for the data type
 
 
std::string DataType() const;
 
        // Input                ML5T    : MAT version 5 tag (this)
        // Output               string  : String labeling the type


std::string DataSymbol() const;

        // Input                ML5T    : MAT version 5 tag (this)
        // Output               string  : String labeling the type
        //                                according to MATLAB symbol 
 
int Bytes() const;

        // Input                ML5T    : MAT version 5 tag (this)
        // Output               NBytes  : Number of bytes


void SetCompressed(int bigend);

        // Input                ML5T    : MAT version 5 tag (this)
	//			bigend  : Data storage endian type
        // Output               T/F     : True if tag is compressed
        // Note                         : Returned value is # of bytes 
        //                                if true 


// ____________________________________________________________________________
// C                MATLAB MAT 5 Tag Binary Output Functions
// ____________________________________________________________________________

   
/* These are the functions which will output a tag in MATLAB ".mat" binary
   format, version 5.  Note that in a valid MATLAB MAT file, this tag will
   be associated with a "Data Element".
 
                Input           ML5H    : MAT version 5 tag (this)
                                fp      : Output file stream
                Output          void    : Output file stream is modified
                                          by having the ML5H written to it
                Note                    : No care is taken to insure the
                                          tag is written tothe top of the
                                          file (as it should be in MATLAB)   */
 
int write(std::fstream& fp) const;
 
/* These are the functions which will return the number of bytes that are
   written upon output a Tag in MATLAB ".mat" binary format, Ver. 5. Note
   that all Tags are 8 bytes long (even those of compressed data elements)   */

int Size() const;
       

// ____________________________________________________________________________
// D                MATLAB MAT 5 Tag Binary Input Functions
// ____________________________________________________________________________

/* These are the functions which will input a tag in MATLAB ".mat" binary       
   format, version 5.  Note that in a valid MATLAB MAT file, this tag will
   associated with a "Data Element", a combination of Tag + Data
 
                Input           ML5T    : MAT version 5 tag (this)
                                fp      : Input file stream
				bigend  : Data storage endian type
                                warn	: Warning output level
                                                0 = no warnings
                                                1 = warnings
                                               >1 = fatal warnings
                Output          void    : ML5T is set from values found
                                          in the input stream
                Note                    : No care is taken to insure the
                                          tag is read from the top of a data
                                          element (as it should be in MATLAB)*/

int read(std::fstream& fp, int bigend, int warn=1);

 
// ____________________________________________________________________________
// E                  MATLAB MAT 5 Tag ASCII Output Functions
// ____________________________________________________________________________

                                                                                
void print(std::ostream& ostr, int hf=1, int bg=0) const;
 
        // Input                ML5T    : MAT version 5 tag (this)
        //                      ostr    : An output stream
	//			hf      : Flag to print header
	//			bf      : Flag to print tag chars
        // Output               ostr    : Output stream with MATLAB MAT file
        //                                tag version 5 info placed into it
        // Note                         : This information exists at the
        //                                beginning of each data element
 
friend std::ostream& operator<< (std::ostream& ostr, const MatLab5Tag& MLT5);
 
        // Input                ostr    : An output stream
        //                      ML5T    : MAT version 5 tag
        // Output               ostr    : The output stream modified by
        //                                the tag parameters
  };


#endif							 // _GML5TAG_h_	
