/* ML4DElem.cc **************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**     MATLAB Data Element (MAT Version 4)        Implementation	**
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
** MATLAB data element in MAT version 4 files.  MATLAB is a commercial  **
** product and can be purchaced from The MathWorks, Inc.		**
** This module deals with binary MAT files version 4 circa 1999.      	**
**                                                                      **
*************************************************************************/

#ifndef   GML4DE_CC_				// Is file already included?
#  define GML4DE_CC_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#  endif

#include <GamIO/ML4DElem.h>			// Include the header file
#include <GamIO/ML4Tag.h>			// Include MAT V4 tags/headers
#include <Basics/Gutils.h>                      // Include GAMMA errors
#include <Matrix/matrix.h>			// Include GAMMA matrices
#include <string>				// Include libstdc++ strings
#include <fstream>				// Include file streams
#include <Basics/StringCut.h>		  // Include Gdec function
//#include <stream.h>				// Include dec function


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i              MATLAB MAT 4 Data Element Error Handling
// ____________________________________________________________________________

/*              Input           MLDE4   : MAT version 4 data element
                                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */  

void MatLab4DE::MLDE4error(int eidx, int noret)
  {                                                     
  std::string hdr("MATLAB MAT V4 Data Element");
  switch(eidx)
    {
    case 0: GAMMAerror(hdr, 0, noret); break;   // Program Aborting        (0)
    case 1: GAMMAerror(hdr, "End of File Reached Before Data Found", noret);
      break;                                     //                        (1)
    case 10: GAMMAerror(hdr, "Cannot Read Data Element Header", noret);	// (10)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }

 
void MatLab4DE::MLDE4error(int eidx, const std::string& pname, int noret)
  {
  std::string hdr("MATLAB MAT V4 Data Element");
  std::string msg;
  switch(eidx)
    {
    case 1:  GAMMAerror(hdr,   1, pname, noret); break;	// File Problems   (1)
    case 2:  GAMMAerror(hdr,   2, pname, noret); break;	// !Read SingPar   (2)
    default: GAMMAerror(hdr, -11, pname, noret); break;	// Unknown Error   (-1)
    }
  }
     
 
volatile void MatLab4DE::MLDE4fatality(int eidx)
  {
  MLDE4error(eidx, 1);
  if(eidx) MLDE4error(0);
  GAMMAfatal();					// Clean exit from program
  }
 

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

int MatLab4DE::WriteMx(std::fstream& fp) const
  {
  double d;
  int i,j;
  for(j=0; j<Tag.ncols; j++)		// First output the real part
    for(i=0; i<Tag.mrows; i++) 		// We output the matrix reals
      { 				// column by column
      d = Data.getRe(i,j);
      fp.write((char*)&d, sizeof(double));
      } 
  if(Tag.cmplxf)			// Next output the imaginary part
    {					// if they are desired
    for(j=0; j<Tag.ncols; j++) 		// We output the matrix imaginaries
      for(i=0; i<Tag.mrows; i++)	// column by column
        {
        d = Data.getIm(i,j);
        fp.write((char*)&d, sizeof(double));
        } 
    }
  return 1;
  }

 
int MatLab4DE::ReadMx(std::fstream& fp, int warn)
  {
  Data = matrix(Tag.mrows, Tag.ncols);	// Set data as an empty array
  double d;				// Temp double
  int i,j;				// Temp indices
  for(j=0; j<Tag.ncols; j++)		// MATLAB array is stored column
    { 					// by column, reals first so this
    for(i=0; i<Tag.mrows; i++)		// is the order in which it is read
      { 
      fp.read((char*)&d, sizeof(double));
// sosi - may need byte swapping & conversion here
      Data.put(d,i,j);
      } 
    } 
  if(Tag.cmplxf)			// If array is complex, next comes
    {					// imaginary elements, column
    for(j=0; j<Tag.ncols; j++)		// by column
      { 
      for(i=0; i<Tag.mrows; i++)
        { 
        fp.read((char*)&d, sizeof(double));
// sosi - may need byte swapping & conversion here
        Data.put(Data.get(i,j)+complex(0,d),i,j);
        } 
      }
    } 
// sosi warning still off
  if(warn) return 1;
  return 1;
  }


int MatLab4DE::WriteName(std::fstream& fp) const
  { fp.write(Name.c_str(), Tag.nlen*sizeof(char));  return 1; }

int MatLab4DE::ReadName(std::fstream& fp, int warn)
  {
  char nc[80];				// An array of chars for name
  fp.read((char*)&nc, Tag.nlen*sizeof(char));	// Read in the characters
  Name = nc;				// Store as a string
// sosi warning still off
  if(warn) return 1;
  return 1;
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                  MATLAB MAT 4 Data Element Constructors
// ____________________________________________________________________________


MatLab4DE::MatLab4DE() { }
MatLab4DE::MatLab4DE(const std::string& N, const matrix& mx, int cmplx)
          :Tag(N,mx,cmplx), Name(N), Data(mx) { }
MatLab4DE::~MatLab4DE() { }
 

// ____________________________________________________________________________
// B                MATLAB MAT 4 Data Element Access Functions
// ____________________________________________________________________________


void MatLab4DE::whos(std::ostream& ostr, std::fstream& fp)

        // Input                ML4DE	: MAT version 4 data element (this)
        //                      ostr    : An output stream
	//			fp	: A file stream
        // Output               ostr    : Output stream with MATLAB MAT file
        //                                data element version 4 info placed
	//				  into it in MATLAB "whos" fashion
        // Note                         : It is assumed that fp points to the
        //                                beginning of a data element and that
	//				  it is open for reading.
	// Note				: The file stream position will be
	//				  advanced to the end of the element

  {
  if(!fp) return;			// No output if no filestream
  int pos = fp.tellp();			// Get current file position
  int TF = Tag.read(fp, 1);		// Read data element tag
  TF *= ReadName(fp, 1);		// Read array name sub-element
  ostr << "\n" << std::string(12, ' ');	// Start new line & output margin
  int len = Name.length();		// Get array name length
  if(len > 10)				// Output the array name
    ostr << Name.substr(0,7) << "...";	// within the 9 columns in whos
  else					// output allocated for it
    ostr << Name << std::string(10-len, ' ');
  ostr << "  ";				// Output a spacer
  std::string dimstr = Gdec(int(Tag.mrows))	// String for dimension output
                + std::string(" x ")
                + Gdec(int(Tag.ncols));
  len = dimstr.length();		// Get string length
  ostr << dimstr << std::string(10-len, ' ');// Output dimensions
  int NB = 11111;			// Get the bytes
// sosi still need byte count set right
  if(Tag.cmplxf) NB +=2;		// Twice a many if complex
  std::string byt = Gdec(NB);			// Set bytes in string form
  len = byt.length();			// Length of bytes string
  ostr << "  ";				// Output a spacer
  ostr << byt << std::string(10-len, ' ');	// Output the bytes
  ostr << "  ";				// Output a spacer
//  ostr << AF.SClass();			// Ouput the MATLAB class type
// sosi still need class set right
  fp.seekp(pos);			// Set original current file position
  }


int MatLab4DE::Size(std::fstream& fp)
  {
  if(!fp) return 0;			// No output if no filestream
  int pos = fp.tellp();			// Get current file position
  Tag.read(fp, 1);			// Read data element tag
  fp.seekp(pos);			// Set original current file position
  return Tag.Bytes();			// Return specified size in tag
  }


void MatLab4DE::skip(std::fstream& fp)
  {
  if(!fp) return;			// No output if no filestream
  Tag.read(fp, 1);			// Read data element tag
// sosi still need to skip more
  }					// jump the file pointer after it



// ____________________________________________________________________________
// C             MATLAB MAT 4 Data Element Binary Output Functions
// ____________________________________________________________________________

/* These are the functions which will output a data element in MATLAB ".mat"
   binary format, version 4. This is done with a pair {tag,data}.  
 
                Input           ML4DE	: MAT version 4 data element (this)
                                fp      : Output file stream
                Output          void    : Output file stream is modified
                                          by having the ML4DE written to it
                Note                    : No care is taken to insure the
                                          data element is written &
					  written at data element start.      */

int MatLab4DE::write(std::fstream& fp) const
  {
  int TF = Tag.write(fp);		// Write main data element tag
  TF *= WriteName(fp);			// Write the data element name
  TF *= WriteMx(fp);			// Write the data
  return TF;
  }

/* These are the functions which will return the number of bytes that are
   written upon output a data element in MATLAB ".mat" binary format, Ver. 4 */

int MatLab4DE::Size(const matrix& mx, const std::string& name, int cmplx) const
  {
  int NB  = Tag.Size();			// Bytes from data element tag
  NB += Tag.nlen*sizeof(char);		// Bytes from data element name
// sosi need to add data length
  return NB;
  }

// ____________________________________________________________________________
// D             MATLAB MAT 4 Data Element Binary Input Functions
// ____________________________________________________________________________

/* These are the functions which will input a data element in MATLAB ".mat"
   binary format, version 4.  Note that in a valid MATLAB MAT file a 
   "Data Element" is a combination of Tag + Data.
 
                Input           ML4DE   : MAT version 4 dim. array (this)
                                fp      : Input file stream
                                warn	: Warning output level
						0 = no warnings
						1 = warnings
                                               >1 = fatal warnings
                Output          void    : ML4DE is set from values found
                                          in the input stream
                Note                    : No care is taken to insure the
                                          dim. array is read from the top of
 					  a data element (as it should be 
					  in MATLAB)                         */

int MatLab4DE::read(std::fstream& fp, int warn)
  {
  int TF = Tag.read(fp, (warn)?1:0);	// Read data element tag
  TF *= ReadName(fp, (warn)?1:0);	// Read data element name
  TF *= ReadMx(fp, (warn)?1:0);		// Read data element data
  if(TF)
    {					// If this has problems, then
    if(warn)				// we must deal with them
      {					// appropriately now
      if(warn==1) MLDE4error(10, 1);
      else        MLDE4fatality(10);
      }
    return 0;
    }
  return 1;
  }

 
// ____________________________________________________________________________
// E            MATLAB MAT 4 Data Element ASCII Output Functions
// ____________________________________________________________________________

                                                                                
void MatLab4DE::print(std::ostream& ostr) const

        // Input                ML4DE	: MAT version 4 data element (this)
        //                      ostr    : An output stream
        // Output               ostr    : Output stream with MATLAB MAT file
        //                                data element version 4 info placed
	//				  into it

  {
  ostr << "\n\t\tData Element";			// Output a header
  Tag.print(ostr, 0);				// Print data element tag
  ostr << "\n\t\t  Name:          " << Name;	// Number of data columns
  ostr.flush();
  }
 
std::ostream& operator<< (std::ostream& ostr, const MatLab4DE& ML4DE)
 
        // Input                ostr    : An output stream
        //                      ML4DE	: MAT version 4 data element
        // Output               ostr    : The output stream modified by
        //                                the data element parameters

  { ML4DE.print(ostr); return ostr; }

#endif							// ML4DElem.cc
