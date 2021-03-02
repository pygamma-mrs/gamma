/* ML5DElem.h ***************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A                           **
**                                                                      **
**      MATLAB Data Element (MAT Version 5)             Interface	**
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
** MATLAB data element in MAT version 5 files.  MATLAB is a commerical  **
** product and can be purchaced from The MathWorks, Inc.     		**
** This  module  deals with binary MAT files version 5 circa 1999.      **
**                                                                      **
*************************************************************************/

#ifndef   _GML5DELEM_h_	  			// Is file already included?
#define   _GML5DELEM_h_   			// If no, then remember it

#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// This is the interface 
#  endif

#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <GamIO/ML5Tag.h>			// Include main MAT tags
#include <GamIO/ML5SubE.h>			// Include MAT sub-elements
#include <GamIO/ML5AF.h>			// Include MAT array flags SE
#include <GamIO/ML5DA.h>			// Include MAT dim. array SE
#include <GamIO/ML5AN.h>			// Include MAT array name SE
#include <GamIO/ML5Reals.h>			// Include MAT reals SE
#include <GamIO/ML5Imags.h>			// Include MAT imags SE

class MatLab5DE
  {
  MatLab5Tag Tag;				// The main data element tag
  MatLab5AF  AF;				// The array flags sub-element
  MatLab5DA  DA;				// The dim. array sub-element
  MatLab5AN  AN;				// The array name sub-element
  MatLab5Re  RData;				// The reals sub-element
  MatLab5Im  IData;				// The imags sub-element
  int        BigEndIn;				// Flag for big endian input
  matrix     MX;				// Array of data for output
 
  friend class MatLabFile;			// Allow MatLabFile full access

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// i                MATLAB MAT 5 Data Element Error Handling
// ____________________________________________________________________________
 
/*              Input           MLDE5   : MAT version 5 data element
                                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */  
 
void MLDE5error(int eidx, int noret=1);
void MLDE5error(int eidx, const std::string& pname, int noret=1);
volatile void MLDE5fatality(int eidx);


//public:	/* Keeping Everything Private, Allow MatLabFile Use!!! */

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

 
// ____________________________________________________________________________
// A           MATLAB MAT 5 Data Element Constructors
// ____________________________________________________________________________
 

MatLab5DE();
MatLab5DE(const matrix& mx, const std::string& Name, int cmplx=1);
virtual ~MatLab5DE();

// ____________________________________________________________________________
// B            MATLAB MAT 5 Data Element Access Functions
// ____________________________________________________________________________
 
void whos(std::ostream& ostr, std::fstream& fp, int bigend);
 
        // Input                ML5DE   : MAT version 5 data element (this)
        //                      ostr    : An output stream
        //                      fp      : A file stream
	//			bigend	: Data storage endian type
        // Output               ostr    : Output stream with MATLAB MAT file
        //                                data element version 5 info placed
        //                                into it in MATLAB "whos" fashion
        // Note                         : It is assumed that fp points to the
        //                                beginning of a data element and that
        //                                it is open for reading.
        // Note                         : The file stream position will be
        //                                advanced to the end of the element
 
int    Size    (std::fstream& fp, int bigend);
std::string Name    (std::fstream& fp, int bigend);
void   skip    (std::fstream& fp, int bigend);


// ____________________________________________________________________________
// C            MATLAB MAT 5 Data Element Binary Output Functions
// ____________________________________________________________________________

/* These are the functions which will output a data element in MATLAB ".mat"
   binary format, version 5. This is done with a pair {tag,data}.  However,
   the data itself can be broken up into sub-elements, each of which will have
   the {tag,data} structure as well.  For matrix output the data is actually
   5 sub-elements of the following:

          Sub-Element     GAMMA Class    Data Type        Number of Bytes
     -------------------  -----------  ------------   ------------------------
     1. Array Flags        MatLab5AF   miUINT32 (6)   2*miUINT32  (8 bytes)
     2. Dimensions Array   MatLab5DA   miINT32  (5)   ND*miINT32  (ND*4 bytes)
     3. Array Name         MatLab5AN   miINT8   (1)   NC*miINT8   (NC bytes)
     4. Real Part          MatLab5Re   miDOUBLE (9)   NE*miDOUBLE (NE*8 bytes)
     5. Imaginary Part     MatLab5Im   miDOUBLE (9)   NE*miDOUBLE (NE*8 bytes)

   Thus our array output will consist of a Tag followed by the 5 Sub-Elements.
 
                Input           ML5DE   : MAT version 5 data element (this)
                                fp      : Output file stream
                Output          void    : Output file stream is modified
                                          by having the ML5DE written to it
                Note                    : No care is taken to insure the
                                          data element is written &
                                          written at data element start.      */ 
 
int write(std::fstream& fp) const;
int write(std::fstream& fp,const matrix& mx,    const std::string& AN) const;
int write(std::fstream& fp,const row_vector& rv,const std::string& AN) const;
int write(std::fstream& fp,const col_vector& cv,const std::string& AN) const;
 
 
/* These are the functions which will return the number of bytes that are
   written upon output a data element in MATLAB ".mat" binary format, Ver. 5 */
 
int DataSize(const matrix& mx,const std::string& name,int cmplx=0) const;
int Size(const matrix& mx, const std::string& name, int cmplx=0) const;

// ____________________________________________________________________________
// D          MATLAB MAT 5 Data Element Binary Input Functions
// ____________________________________________________________________________

/* These are the functions which will input an data element data element in 
   MATLAB ".mat" binary format, version 5.  Note that in a valid MATLAB MAT
   file, an "Data Element" is a combination of Tag + Data that has
   an 8 byte tag and no data.
 
                Input           ML5DE 	: MAT version 5 tag (this)
                                fp      : Input file stream
				bigend	: Data storage endian type
                                warn 	: Warning output level
                                                0 = no warnings
                                                1 = warnings
                                               >1 = fatal warnings
                Output          void    : ML5DE is set from values found
                                          in the input stream
                Note                    : No care is taken to insure the
                                          tag is read from the top of a data
                                          element (as it should be in MATLAB)*/

int read(std::fstream& fp, int bigend, int warn=1);
matrix GetMatrix(std::fstream& fp, int bigend, int warn=1);

 
// ____________________________________________________________________________
// E            MATLAB MAT 5 Data Element ASCII Output Functions
// ____________________________________________________________________________

                                                                                
void print(std::ostream& ostr) const;
 
        // Input                ML5DE	: MAT version 5 data element (this)
        //                      ostr    : An output stream
        // Output               ostr    : Output stream with MATLAB MAT file
        //                                data element version 5 info placed
	//				  into it
        // Note                         : This information exists at the
        //                                beginning of each data element
 
friend std::ostream& operator<< (std::ostream& ostr, const MatLab5DE& ML5DE);
 
        // Input                ostr    : An output stream
        //                      ML5DE	: MAT version 5 data element
        // Output               ostr    : The output stream modified by
        //                                the data element parameters

// ____________________________________________________________________________
// F       MATLAB MAT 5 Data Element ASCII Output Auxiliary Functions
// ____________________________________________________________________________

/* These functions are just used for testing purposes.  The mirror the binary
   write functions, but instead of writing in binary to a filestream the write
   in ASCII to an output stream.  Thus it is easy to see most of what is
   being written, at least in principle.....                                 */
                                            
void print(std::ostream& ostr,
                           const matrix& mx, const std::string& N, int cmplx) const;
};


#endif							 // _GML5DELEM_h_
