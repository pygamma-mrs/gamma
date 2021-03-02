/* ML5DA.h ******************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A                           **
**                                                                      **
**      MATLAB Dimensions Array (MAT Version 5)        Interface	**
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
** MATLAB array flags sub-element in MAT version 5 files.  MATLAB is a  **
** commercial product and can be purchaced from The MathWorks, Inc.     **
** This  module  deals with binary MAT files version 5 circa 1999.      **
**                                                                      **
*************************************************************************/

#ifndef   _GML5DA_h_				// Is file already included?
#define   _GML5DA_h_  1				// If no, then remember it

#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// This is the interface 
#  endif

#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <GamIO/ML5SubE.h>			// Include base class
#include <Matrix/matrix.h>                      // Include GAMMA matrices
#include <Matrix/row_vector.h>                  // Include GAMMA row vectors
#include <Matrix/col_vector.h>                  // Include GAMMA column vectors
#include <fstream>				// Include libstdc++ filestreams
#include <string>				// Include libstdc++ strings

class MatLab5DA : MatLab5SE
  {
  int    ND;					// Number of dimensions
  int*   Sizes;					// Array of dimensions
 
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// i             MATLAB MAT 5 Dimensions Array Sub-Element Error Handling
// ____________________________________________________________________________
 
/*              Input           MLDA5   : MAT version 5 sub-element
                                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */  
 
void MLDA5error(int eidx, int noret=1);
void MLDA5error(int eidx, const std::string& pname, int noret=1);
volatile void MLDA5fatality(int eidx);


public:

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

 
// ____________________________________________________________________________
// A           MATLAB MAT 5 Dimensions Array Sub-Element Constructors
// ____________________________________________________________________________
 

MatLab5DA();
MatLab5DA(const matrix& mx);

        // Input                ML5DA   : MAT version 5 dimensions array (this)
        //                      mx      : A GAMMA matrix
        // Output               this    : Constructed appropriately
        //                                for the dimensions array sub-element
        //                                of a GAMMA array
        // Note                         : The tag is set to be of type
        //                                5 (i.e. MiINT32) and the size
        //                                is set to eight bytes long
        // Note                           The data is as long row, long cols
        //                                & output size is eight bytes as is
        //                                indicated in the preceding tag

virtual ~MatLab5DA();

// ____________________________________________________________________________
// B         MATLAB MAT 5 Dimensions Array Sub-Element Access Functions
// ____________________________________________________________________________

int Dims()     const;
int Dim(int i) const;
 
// ____________________________________________________________________________
// C     MATLAB MAT 5 Dimensions Array Sub-Element Binary Output Functions
// ____________________________________________________________________________

/* These are the functions which will output this sub-element in MATLAB ".mat"
   binary format, version 5.  This is done with a pair {tag,data}.  Both the
   tag & data will be written directly from these functions, bypassing the
   associated class structures which have no knowledge of GAMMA array types.

   Remember, the MATLAB MAT version 5 tag for this sub-element contains
   8-bytes with two important settings: 1.) data type & 2.) # of bytes.
   The data type is placed in the first 4 bytes and the number of bytes
   placed into the last 4 bytes (for non-compressed data elements which we
   will always output...).  In this case, we output the elements as 32-bit
   integers, so the data type should be MATLAB type miINT32, i.e. value 5.

                Input           ML5DA   : MAT version 5 dimensions array (this)
                                fp      : Output file stream
                                mx      : GAMMA matrix
                                rv      : GAMMA row vector
                                cv      : GAMMA column vector
                Output          void    : Output file stream is modified
                                          by having the ML5DA written to it
                Note                    : No care taken herein to insure the
                                          sub-element is properly written.   */
  
int write(std::fstream& fp)                       const;
int write(std::fstream& fp, const matrix& mx)     const;
int write(std::fstream& fp, const row_vector& rv) const;
int write(std::fstream& fp, const col_vector& cv) const;

  
/* These are the functions which will return the number of bytes that are  
   written upon output this sub-element in MATLAB ".mat" binary format, V5   */ 
  
int Size(const matrix& mx) const;


// ____________________________________________________________________________
// D          MATLAB MAT 5 Dimensions Array Sub-Element Binary Input Functions
// ____________________________________________________________________________

/* These are the functions which will input an array flags sub-element in 
   MATLAB ".mat" binary format, version 5.  Note that in a valid MATLAB MAT
   file, an "Dimensions Array Sub-Element" is a combination of Tag + Data that has
   an 8 byte tag and no data.
 
                Input           ML5DA 	: MAT version 5 tag (this)
                                fp      : Input file stream
				bigend	: Data storage endian type
                                warn	: Warning output level
                                                0 = no warnings
                                                1 = warnings
                                               >1 = fatal warnings
                Output          void    : ML5DA is set from values found
                                          in the input stream
                Note                    : No care is taken to insure the
                                          tag is read from the top of a data
                                          element (as it should be in MATLAB)*/

int read(std::fstream& fp, int bigend, int warn=1);

 
// ____________________________________________________________________________
// E      MATLAB MAT 5 Dimensions Array Sub-Element ASCII Output Functions
// ____________________________________________________________________________

                                                                                
void print(std::ostream& ostr) const;
 
        // Input                ML5DA	: MAT version 5 array flags (this)
        //                      ostr    : An output stream
        // Output               ostr    : Output stream with MATLAB MAT file
        //                                array flags version 5 info placed
	//				  into it
        // Note                         : This information exists at the
        //                                beginning of each array flags
 
friend std::ostream& operator<< (std::ostream& ostr, const MatLab5DA& ML5DA);
 
        // Input                ostr    : An output stream
        //                      ML5DA	: MAT version 5 array flags
        // Output               ostr    : The output stream modified by
        //                                the array flags parameters
  };


#endif							 // _GML5DA_h_
