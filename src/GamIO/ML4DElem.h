/* ML4DElem.h ***************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A                           **
**                                                                      **
**      MATLAB Data Element (MAT Version 4)             Interface	**
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
** MATLAB data element in MAT version 4 files.  MATLAB is a commerical  **
** product and can be purchaced from The MathWorks, Inc.     		**
** This  module  deals with binary MAT files version 4 circa 1999.      **
**                                                                      **
*************************************************************************/

#ifndef   _GML4DELEM_h_				// Is file already included?
#define   _GML4DELEM_h_ 				// If no, then remember it

#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// This is the interface 
#  endif

#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <GamIO/ML4Tag.h>			// Include MAT V4 tags/headers

class MatLab4DE
  {
  MatLab4Tag Tag;				// The main data element tag
  std::string     Name;				// Element name
  matrix     Data;				// Element data
 
  friend class MatLabFile;			// Allow MatLabFile full access

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// i                MATLAB MAT 4 Data Element Error Handling
// ____________________________________________________________________________
 
/*              Input           MLDE4   : MAT version 4 data element
                                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */  
 
void MLDE4error(int eidx, int noret=1);
void MLDE4error(int eidx, const std::string& pname, int noret=1);
volatile void MLDE4fatality(int eidx);


// ____________________________________________________________________________
// ii       MATLAB MAT 4 Data Element Specialized Read/Write Functions
// ____________________________________________________________________________
 
/* These functons perform binary I/O of the "data" components of an MATLAB
   MAT version 4 data element. There are two components in each element, a
   name (series of characters ending in \0) followed by data.  The functions
   below will read/write these directly from a file.  They do NOT have any
   regard over file positioning and utilize the Tag/Header information.  Thus
   the Tag values should be proper and the name/data will be related to the
   string Name and the matrix Data.  Upon read, the data elements will
   perform the required conversion(s) as needed.

           Input                fp       : A file stream
           Output               void/int : The data array and array name
                                           are written/read to/from the file
	  				   stream in MATLAB MAT V4 format    */

int WriteMx(std::fstream& fp) const;
int ReadMx(std::fstream& fp, int warn=1);
int WriteName(std::fstream& fp) const;
int ReadName(std::fstream& fp, int warn=1);


//public:	/* Keeping Everything Private, Allow MatLabFile Use!!! */

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

 
// ____________________________________________________________________________
// A           MATLAB MAT 4 Data Element Constructors
// ____________________________________________________________________________
 

MatLab4DE();
MatLab4DE(const std::string& Name, const matrix& mx, int cmplx=1);
virtual ~MatLab4DE();

// ____________________________________________________________________________
// B            MATLAB MAT 4 Data Element Access Functions
// ____________________________________________________________________________
 
void whos(std::ostream& ostr, std::fstream& fp);
 
        // Input                ML4DE   : MAT version 4 data element (this)
        //                      ostr    : An output stream
        //                      fp      : A file stream
        // Output               ostr    : Output stream with MATLAB MAT file
        //                                data element version 4 info placed
        //                                into it in MATLAB "whos" fashion
        // Note                         : It is assumed that fp points to the
        //                                beginning of a data element and that
        //                                it is open for reading.
        // Note                         : The file stream position will be
        //                                advanced to the end of the element
 
int  Size(std::fstream& fp);
void skip(std::fstream& fp);


// ____________________________________________________________________________
// C            MATLAB MAT 4 Data Element Binary Output Functions
// ____________________________________________________________________________

/* These are the functions which will output a data element in MATLAB ".mat"
   binary format, version 4. This is done with a pair {tag,data}.  However,
   the data itself can be broken up into sub-elements, each of which will have
   the {tag,data} structure as well.  For matrix output the data is actually
   4 sub-elements of the following:

          Sub-Element     GAMMA Class    Data Type        Number of Bytes
     -------------------  -----------  ------------   ------------------------
     1. Array Flags        MatLab5AF   miUINT32 (6)   2*miUINT32  (8 bytes)
     2. Dimensions Array   MatLab5DA   miINT32  (5)   ND*miINT32  (ND*4 bytes)
     3. Array Name         MatLab5AN   miINT8   (1)   NC*miINT8   (NC bytes)
     4. Real Part          MatLab5Re   miDOUBLE (9)   NE*miDOUBLE (NE*8 bytes)
     5. Imaginary Part     MatLab5Im   miDOUBLE (9)   NE*miDOUBLE (NE*8 bytes)

   Thus our array output will consist of a Tag followed by the 5 Sub-Elements.
 
                Input           ML4DE   : MAT version 4 data element (this)
                                fp      : Output file stream
                Output          void    : Output file stream is modified
                                          by having the ML4DE written to it
                Note                    : No care is taken to insure the
                                          data element is written &
                                          written at data element start.      */ 
 
int write(std::fstream& fp) const;
 
 
/* These are the functions which will return the number of bytes that are
   written upon output a data element in MATLAB ".mat" binary format, Ver. 4 */
 
int Size(const matrix& mx, const std::string& name, int cmplx=0) const;

// ____________________________________________________________________________
// D          MATLAB MAT 4 Data Element Binary Input Functions
// ____________________________________________________________________________

/* These are the functions which will input an data element data element in 
   MATLAB ".mat" binary format, version 5.  Note that in a valid MATLAB MAT
   file, an "Data Element" is a combination of Tag + Data that has
   an 8 byte tag and no data.
 
        	Input		ML4DE	: MAT version 4 data element (this)
                                fp      : Input file stream
                                warn	: Warning output level
                                                0 = no warnings
                                                1 = warnings
                                               >1 = fatal warnings
                Output          void    : ML4DE is set from values found
                                          in the input stream
                Note                    : No care is taken to insure the
                                          tag is read from the top of a data
                                          element (as it should be in MATLAB)*/

int read(std::fstream& fp, int warn=1);

 
// ____________________________________________________________________________
// E            MATLAB MAT 4 Data Element ASCII Output Functions
// ____________________________________________________________________________


void print(std::ostream& ostr) const;
 
        // Input                ML4DE	: MAT version 4 data element (this)
        //                      ostr    : An output stream
        // Output               ostr    : Output stream with MATLAB MAT file
        //                                data element version 4 info placed
	//				  into it
        // Note                         : This information exists at the
        //                                beginning of each data element
 
friend std::ostream& operator<< (std::ostream& ostr, const MatLab4DE& ML4DE);
 
        // Input                ostr    : An output stream
        //                      ML4DE	: MAT version 4 data element
        // Output               ostr    : The output stream modified by
        //                                the data element parameters
};


#endif							 // _GML4DELEM_h_
