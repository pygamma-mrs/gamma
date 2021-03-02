/* ML5AF.h ******************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A                           **
**                                                                      **
**      MATLAB Array Flags (MAT Version 5)        Interface		**
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

#ifndef GML5AF_h_				// Is file already included?
#  define GML5AF_h_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// This is the interface 
#  endif

#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <GamIO/ML5SubE.h>			// Include base class
#include <Matrix/matrix.h>			// Include GAMMA matrices
//#include <bitset>				// Include libstdc++ bitsets

class MatLab5AF : MatLab5SE
  {
  int  Fs;					// 3rd byte as integer
  int  Class;					// 4th byte as integer

 
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// i             MATLAB MAT 5 Array Flags Sub-Element Error Handling
// ____________________________________________________________________________
 
/*              Input           MLAF5   : MAT version 5 sub-element
                                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */  
 
void MLAF5error(int eidx, int noret=1);
void MLAF5error(int eidx, const std::string& pname, int noret=1);
volatile void MLAF5fatality(int eidx);


public:

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

 
// ____________________________________________________________________________
// A                  MATLAB MAT 5 Array Flags Sub-Element Constructors
// ____________________________________________________________________________
 

MatLab5AF();
MatLab5AF(int cmplx, int C, int global=0,int logical=0);
MatLab5AF(const matrix& mx, int cmplx);
 
        // Input                ML5T    : MAT version 5 tag (this)
        //                      mx      : A GAMMA matrix
        //                      cmplx   : Flag for whether real vs complex
        // Output               this    : Constructed appropriately
        //                                for the array flags sub-element
        //                                of a GAMMA array
        // Note                         : The tag is set to be of type
        //                                6 (i.e. MiUINT32) and the size 
        //                                is set to eight bytes long
        // Note                           The data is set to be of type
        //                                6 (i.e. mxDOUBLE_CLASS) and the
        //                                output size is eight bytes as is
        //                                indicated in the preceding tag 

// ____________________________________________________________________________
// B                MATLAB MAT 5 Array Flags Sub-Element Access Functions
// ____________________________________________________________________________

 
std::string SClass() const;
 
        // Input                ML5T    : MAT version 5 tag (this) 
        // Output               string  : String labeling the type 

std::string Symbol() const;

        // Input                ML5T    : MAT version 5 tag (this)
        // Output               string  : String labeling the type


int IsComplex() const;
 
        // Input                ML5AF   : MAT version 5 array flags SE (this)
        // Output               TF      : True if complex flag set
 
/***********************************************************************
  sosi - will want to replace this when libstdc++ bitset class works!
         or perhaps just get rid of it entirely!
***********************************************************************/
 
int IsFlag(int flg) const;

   
// ____________________________________________________________________________
// C        MATLAB MAT 5 Array Flags Sub-Element Binary Output Functions
// ____________________________________________________________________________

   
/* These are the functions which will output an array flags sub-element in
   MATLAB ".mat" binary format, version 5. 
 
                Input           ML5AF	: MAT version 5 array flags (this)
                                fp      : Output file stream
                Output          void    : Output file stream is modified
                                          by having the ML5AF written to it
                Note                    : No care is taken to insure where
					  the tag is written in the file.    */
 
int write(std::fstream& fp) const;

 
/* These are the functions which will return the number of bytes that are 
   written upon output this sub-element in MATLAB ".mat" binary format, V5   */ 
 
int Size() const;

// ____________________________________________________________________________
// D          MATLAB MAT 5 Array Flags Sub-Element Binary Input Functions
// ____________________________________________________________________________

/* These are the functions which will input an array flags sub-element in 
   MATLAB ".mat" binary format, version 5.  Note that in a valid MATLAB MAT
   file, an "Array Flags Sub-Element" is a combination of Tag + Data that has
   an 8 byte tag and no data.
 
                Input           ML5AF	: MAT version 5 tag (this)
                                fp      : Input file stream
				bigend	: Data storage endian type
                                warn    : Warning output level
                                                0 = no warnings
                                                1 = warnings
                                               >1 = fatal warnings
                Output          void    : ML5AF is set from values found
                                          in the input stream
                Note                    : No care is taken to insure the
                                          tag is read from the top of a data
                                          element (as it should be in MATLAB)*/

virtual int read(std::fstream& fp, int bigend, int warn=1);

 
// ____________________________________________________________________________
// E         MATLAB MAT 5 Array Flags Sub-Element ASCII Output Functions
// ____________________________________________________________________________

                                                                                
void print(std::ostream& ostr) const;
 
        // Input                ML5AF	: MAT version 5 array flags (this)
        //                      ostr    : An output stream
        // Output               ostr    : Output stream with MATLAB MAT file
        //                                array flags version 5 info placed
	//				  into it
        // Note                         : This information exists at the
        //                                beginning of each array flags
 
friend std::ostream& operator<< (std::ostream& ostr, const MatLab5AF& ML5AF);
 
        // Input                ostr    : An output stream
        //                      ML5AF	: MAT version 5 array flags
        // Output               ostr    : The output stream modified by
        //                                the array flags parameters
  };


#endif							 // ML5AF.h
