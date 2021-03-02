/* ML4Tag.h *****************************************************-*-c++-*-
**									**
**	                            G A M M A			 	**
**								 	**
**	MATLAB Tag (MAT Version 4)                     Interface	**
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
** tags associated with a data element in a MATLAB MAT version 4 file.  **
** MATLAB is a commercial product and can be purchaced from The 	**
** MathWorks, Inc. This module deals with binary MAT files version 4 	**
** circa 1999. 								**
**								 	**
** A data element is a { tag, data } pairing, and each version 4 MAT    **
** file is a series of data elements.  Herein we treat the tag. Unlike  **
** MAT version 5 files, MAT version 4 files contain no main (initial)   **
** header. Thus, in older MATLAB nomeclature, the MAT version 4 tag is  **
** equivalent to a header - where each data element has its own header. **
**								 	**
** MAT version 4 tags/headers each contain 20 bytes, 5 int32_t (4 byte)	**
** integers. The first integer, type, is the only non-obvious value.	**
** When taken as an integer of the form MOPT, where M is the value in 	**
** the thousands place, O the value in the hundreds, etc., the table	**
** below applies (from MATLAB External Interface Library documentation) **
**								 	**
**    M:binary format    O:empty         P:data format        T:mx type **
** --------------------  -------  --------------------------  --------- **
** 0:IEEE little endian  0:fixed  0:double precision(64-bit)  0:full 	**
** 1:IEEE big endian              1:single precision(32-bit)  1:test    **
** 2:VAX D-float                  2:32-bit signed integer     2:sparse  **
** 3:VAX G-float                  3:16-bit signed integer		**
** 4:Cray                         4:16-bit unsigned integer		**
**							 		**
*************************************************************************/

#ifndef   GML4TAG_h_				// Is file already included?
#  define GML4TAG_h_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// This is the interface
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


#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <string>				// Know libstdc++ strings
#include <fstream>				// Know libstdc++ filestreams
#include <Matrix/matrix.h>			// Know GAMMA matrices

class MatLab4Tag
  {
  int32_t type;                                    // Computer type
  int32_t mrows;                                   // Number of data rows
  int32_t ncols;                                   // Number of data columns
  int32_t cmplxf;                                  // Flag complex (0=real, 1=cmp)
  int32_t nlen;                                    // Length of name for data
  int MLM; 					// Endian flag
  int MLO;					// Extra flag
  int MLP;					// Data precision flag
  int MLT;					// Matrx type flag



// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                      MATLAB MAT 4 Tag Error Handling
// ____________________________________________________________________________

/* 		Input 		ML4T    : MAT version 4 tag
				eidx    : Error index
               			pname   : string in message
				noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */

void MLT4error(int eidx, int noret=0);
void MLT4error(int eidx, const std::string& pname, int noret=0);
volatile void MLT4fatality(int eidx);


// ____________________________________________________________________________
// ii              MATLAB MAT 4 Tag Binary/Byte Functions
// ____________________________________________________________________________

void MOPT2Type();
void Type2MOPT();

public:

  friend class MatLab4DE;			// Allow data elements access
  friend class MatLabFile;			// Allow MATLAB files access

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                    MATLAB MAT 4 Tag Constructors
// ____________________________________________________________________________


MatLab4Tag();
MatLab4Tag(const std::string& N, const matrix& mx, int cf=1);
MatLab4Tag(const MatLab4Tag& T);
void operator=(const MatLab4Tag& T);

// ____________________________________________________________________________
// B                    MATLAB MAT 4 Tag Access Functions
// ____________________________________________________________________________


std::string MType() const;

        // Input                ML4T    : MAT version 4 tag (this)
        // Output               string  : String labeling the byte
        //                                ordering of stored data

std::string PType() const;

        // Input                ML4T    : MAT version 4 tag (this)
        // Output               string  : String labeling the type
        //                                of data that is stored


std::string TType() const;

        // Input                ML4T    : MAT version 4 tag (this)
        // Output               string  : String labeling the type
        //                                of matrix that is stored


int ElemBytes() const;

        // Input                ML4T    : MAT version 4 tag (this)
        // Output               int     : Number of bytes in an array element


int Bytes() const;

        // Input                ML4T    : MAT version 4 tag (this)
        // Output               NBytes  : Number of bytes


// ____________________________________________________________________________
// C                MATLAB MAT 4 Tag Binary Output Functions
// ____________________________________________________________________________


/* These are the functions which will output a tag in MATLAB ".mat" binary
   format, version 4.  Note that in a valid MATLAB MAT file, this tag will
   be associated with a "Data Element".

                Input           ML4H    : MAT version 4 tag (this)
                                fp      : Output file stream
                Output          void    : Output file stream is modified
                                          by having the ML4H written to it
                Note                    : No care is taken to insure the
                                          tag is written before data in the
                                          file (as it should be in MATLAB)   */

int write(std::fstream& fp) const;

/* These are the functions which will return the number of bytes that are
   written upon output a Tag in MATLAB ".mat" binary format, Ver. 5. Note
   that all Tags are 8 bytes long (even those of compressed data elements)   */

int Size() const;


// ____________________________________________________________________________
// D                MATLAB MAT 4 Tag Binary Input Functions
// ____________________________________________________________________________

/* These are the functions which will input a tag in MATLAB ".mat" binary
   format, version 4.  Note that in a valid MATLAB MAT file, this tag will
   associated with a "Data Element", a combination of Tag + Data

                Input           ML4T    : MAT version 4 tag (this)
                                fp      : Input file stream
                                warn	: Warning output level
                                                0 = no warnings
                                                1 = warnings
                                               >1 = fatal warnings
                Output          void    : ML4T is set from values found
                                          in the input stream
                Note                    : No care is taken to insure the
                                          tag is read from the top of a data
                                          element (as it should be in MATLAB)*/

int read(std::fstream& fp, int warn=1);


// ____________________________________________________________________________
// E                  MATLAB MAT 4 Tag ASCII Output Functions
// ____________________________________________________________________________


void print(std::ostream& ostr, int hf=1) const;

        // Input                ML4T    : MAT version 4 tag (this)
        //                      ostr    : An output stream
	//			hf      : Flag to print header
        // Output               ostr    : Output stream with MATLAB MAT file
        //                                tag version 4 info placed into it
        // Note                         : This information exists at the
        //                                beginning of each data element


friend std::ostream& operator<< (std::ostream& ostr, const MatLab4Tag& MLT4);

        // Input                ostr    : An output stream
        //                      ML4T    : MAT version 4 tag
        // Output               ostr    : The output stream modified by
        //                                the tag parameters
  };


#endif							 // ML4Tag.h
