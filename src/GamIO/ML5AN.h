/* ML5AN.h ******************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A                           **
**                                                                      **
**      MATLAB Array Name (MAT Version 5)  	      Interface		**
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
** MATLAB array name sub-element in MAT version 5 files.  MATLAB is a  	**
** commercial product and can be purchaced from The MathWorks, Inc.     **
** This  module  deals with binary MAT files version 5 circa 1999.      **
**                                                                      **
*************************************************************************/

#ifndef GML5AN_H_				// Is file already included?
#  define GML5AN_H_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// This is the interface 
#  endif

#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <GamIO/ML5SubE.h>			// Include base class
#include <GamIO/ML5Tag.h>			// Include MAT 5 tags
#include <string>

class MatLab5AN : MatLab5SE
  {
  int    NC;					// Number of characters
  std::string MxName;				// Matrix name
 
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// i             MATLAB MAT 5 Array Name Sub-Element Error Handling
// ____________________________________________________________________________
 
/*              Input           MLAN5   : MAT version 5 sub-element
                                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */  
 
void MLAN5error(int eidx, int noret=1);
void MLAN5error(int eidx, const std::string& pname, int noret=1);
volatile void MLAN5fatality(int eidx);


public:

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

 
// ____________________________________________________________________________
// A           MATLAB MAT 5 Array Name Sub-Element Constructors
// ____________________________________________________________________________
 

MatLab5AN();
MatLab5AN(const std::string& name);
virtual ~MatLab5AN();

// ____________________________________________________________________________
// B            MATLAB MAT 5 Array Name Sub-Element Access Functions
// ____________________________________________________________________________

std::string Name();
 
// ____________________________________________________________________________
// C        MATLAB MAT 5 Array Name Sub-Element Binary Output Functions
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
   characters for the array name.

                Input           ML5AN   : MAT version 5 array name SE (this)
                                fp      : Output file stream
                                aname   : Array name
                Output          void    : Output file stream is modified
                                          by having the ML5AN written to it
                Note                    : No care is taken to insure the
                                          tag is written.                    */  
 
   
int write(std::fstream& fp) const;
int write(std::fstream& fp, const std::string& aname) const;

/* These are the functions which will return the number of bytes that are
   written upon output this sub-element in MATLAB ".mat" binary format, V5   */

int Size(const std::string& name) const;



// ____________________________________________________________________________
// D          MATLAB MAT 5 Array Name Sub-Element Binary Input Functions
// ____________________________________________________________________________

/* These are the functions which will input an array name sub-element in 
   MATLAB ".mat" binary format, version 5.  Note that in a valid MATLAB MAT
   file, an "Array Name Sub-Element" is a combination of Tag + Data that has
   an 8 byte tag and no data.
 
                Input           ML5AN	: MAT version 5 tag (this)
                                fp      : Input file stream
				bigedn	: Data storage endian type
                                warn	: Warning output level
                                                0 = no warnings
                                                1 = warnings
                                               >1 = fatal warnings
                Output          void    : ML5AN is set from values found
                                          in the input stream
                Note                    : No care is taken to insure the
                                          tag is read from the top of a data
                                          element (as it should be in MATLAB)*/

int read(std::fstream& fp, int bigend, int warn=1);

 
// ____________________________________________________________________________
// E      MATLAB MAT 5 Array Name Sub-Element ASCII Output Functions
// ____________________________________________________________________________

                                                                                
void print(std::ostream& ostr) const;
 
        // Input                ML5AN	: MAT version 5 array name (this)
        //                      ostr    : An output stream
        // Output               ostr    : Output stream with MATLAB MAT file
        //                                array name version 5 info placed
	//				  into it
        // Note                         : This information exists at the
        //                                beginning of each array name
 
friend std::ostream& operator<< (std::ostream& ostr, const MatLab5AN& ML5AN);
 
        // Input                ostr    : An output stream
        //                      ML5AN	: MAT version 5 array name
        // Output               ostr    : The output stream modified by
        //                                the array name parameters
  };


#endif							 // ML5AN.h
