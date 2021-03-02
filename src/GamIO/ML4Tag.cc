/* ML4Tag.cc ****************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A                           **
**                                                                      **
**	MATLAB Tag (MAT Version 4)                 Implementation	**
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
** tags associated with a data element in a MATLAB MAT version 4 file.  **
** MATLAB is a commercial product and can be purchaced from The         **
** MathWorks, Inc. This module deals with binary MAT files version 4    **
** circa 1999.                                                          **
**                                                                      **
** A data element is a { tag, data } pairing, and each version 4 MAT    **
** file is a series of data elements.  Herein we treat the tag. Unlike  **
** MAT version 4 files, MAT version 4 files contain no main (initial)   **
** header. Thus, in older MATLAB nomeclature, the MAT version 4 tag is  **
** equivalent to a header - where each data element has its own header. **
**                                                                      **
** MAT version 4 tags/headers each contain 20 bytes, 5 int32_t (4 byte)    **
** integers. The first integer, type, is the only non-obvious value.    **
** When taken as an integer of the form MOPT, where M is the value in   **
** the thousands place, O the value in the hundreds, etc., the table    **
** below applies (from MATLAB External Interface Library documentation) **
**                                                                      **
**    M:binary format    O:empty         P:data format        T:mx type **
** --------------------  -------  --------------------------  --------- **
** 0:IEEE little endian  0:fixed  0:double precision(64-bit)  0:full    **
** 1:IEEE big endian              1:single precision(32-bit)  1:test    **
** 2:VAX D-float                  2:32-bit signed integer     2:sparse  **
** 3:VAX G-float                  3:16-bit signed integer               **
** 4:Cray                         4:16-bit unsigned integer             **
**                                                                      **
*************************************************************************/

#ifndef   GML4TAG_CC_				// Is file already included?
#  define GML4TAG_CC_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#  endif

#include <GamIO/ML4Tag.h>			// Include the header file
#include <GamIO/BinIOBase.h>			// Include byte swapping
#include <Basics/Gutils.h>                      // Include GAMMA errors
#include <string>				// Include libstdc++ strings
#include <stdio.h>
#include <time.h>				// Include time/date functions

#if !defined(_MSC_VER) && !defined(__MINGW32__) // MSVC++ & MinGW do not use this
# include <sys/utsname.h>			// Include platform information
#endif

#include <Matrix/matrix.h>			// Include GAMMA matrices


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                      MATLAB MAT 4 Tag Error Handling
// ____________________________________________________________________________
 
/*              Input           ML4T    : MAT version 4 header
                                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
		Output		none    : Error message output
                                          Execution stopped (if fatal)       */  
 
void MatLab4Tag::MLT4error(int eidx, int noret)
  {                                                     
  std::string hdr("MATLAB MAT V4 Tag");
  switch(eidx)
    {
    case 0: GAMMAerror(hdr, eidx, noret); break;   // Program Aborting        (0)
    case 1: GAMMAerror(hdr, eidx, noret); break;   // Input FileStream Bad    (1)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }

 
void MatLab4Tag::MLT4error(int eidx, const std::string& pname, int noret)
  {
  std::string hdr("MATLAB MAT V4 Tag");
  std::string msg;
  switch(eidx)
    {
    case 1:  GAMMAerror(hdr,   1, pname, noret); break;	// File Problems   (1)
    case 2:  GAMMAerror(hdr,   2, pname, noret); break;	// !Read SingPar   (2)
    default: GAMMAerror(hdr, -11, pname, noret); break;	// Unknown Error   (-1)
    }
  }
     
 
volatile void MatLab4Tag::MLT4fatality(int eidx)
  {
  MLT4error(eidx, 1);
  if(eidx) MLT4error(0);
  GAMMAfatal();					// Clean exit from program
  }
 
// ____________________________________________________________________________
// ii  		   MATLAB MAT 4 Tag Binary/Byte Functions
// ____________________________________________________________________________
 
void MatLab4Tag::MOPT2Type()
  { type = int32_t(1000*MLM + 100*MLO + 10*MLP + MLT); }

void MatLab4Tag::Type2MOPT()
  { 
  MLM = (type - type%1000)/1000;
  int32_t tmpx = type - MLM;
  MLO = (tmpx - tmpx%100)/100;
  tmpx -= MLO; 
  MLP = (tmpx - tmpx%10)/10;
  MLT = tmpx - MLP;
  }
 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                      MATLAB MAT 4 Tag Constructors
// ____________________________________________________________________________


MatLab4Tag::MatLab4Tag()
  {
  int x = 1;					// First set endian flag
  if(*(char *)&x == 1) MLM = 0;			// If here, little endian
  else                 MLM = 1;			// else,    big endian
  MLO    = 0;					// Extra flag (always 0)
  MLP    = 0; 					// Default precision (0=double)
  MLT    = 0;					// Default array type (matrix)
  MOPT2Type();					// Set computer type ala MOPT
  mrows  = 0;					// Number of array rows
  ncols  = 0;					// Number of array columns
  cmplxf = 0;					// Flag complex (0=real, 1=cmp)
  nlen   = 0;                                   // Length of name for data
  }

MatLab4Tag::MatLab4Tag(const std::string& N, const matrix& mx, int cf)
  {
  int x = 1;					// First set endian flag
  if(*(char *)&x == 1) MLM = 0;			// If here, little endian
  else                 MLM = 1;			// else,    big endian
  MLO    = 0;					// Extra flag (always 0)
  MLP    = 0; 					// Default precision (0=double)
  MLT    = 0;					// Default array type (matrix)
  MOPT2Type();					// Set computer type ala MOPT
  mrows  = mx.rows();				// Number of array rows
  ncols  = mx.cols();				// Number of array columns
  cmplxf = cf;					// Flag complex (0=real, 1=cmp)
  nlen   = N.length()+1; 			// Length of name for data
  }

MatLab4Tag::MatLab4Tag(const MatLab4Tag& T)
  {
  type   = T.type;				// Copy Computer type (none)
  mrows  = T.mrows;				// Copy Number of array rows
  ncols  = T.ncols;				// Copy number of array columns
  cmplxf = T.cmplxf;				// Copy flag complex
  nlen   = T.nlen;				// Copy length of name
  MLM    = T.MLM;				// Copy binary format
  MLO    = T.MLO;				// Copy extra flag
  MLP    = T.MLP;				// Copy element precision
  MLT    = T.MLT;				// Copy matrix type
  }

void MatLab4Tag::operator=(const MatLab4Tag& T)
  {
  type   = T.type;				// Copy Computer type (none)
  mrows  = T.mrows;				// Copy Number of array rows
  ncols  = T.ncols;				// Copy number of array columns
  cmplxf = T.cmplxf;				// Copy flag complex
  nlen   = T.nlen;				// Copy length of name
  MLM    = T.MLM;				// Copy binary format
  MLO    = T.MLO;				// Copy extra flag
  MLP    = T.MLP;				// Copy element precision
  MLT    = T.MLT;				// Copy matrix type
  }

// ____________________________________________________________________________
// B                    MATLAB MAT 4 Tag Access Functions
// ____________________________________________________________________________


std::string MatLab4Tag::MType() const

        // Input 		ML4T    : MAT version 4 tag (this)
        // Output               string  : String labeling the byte
        //                                ordering of stored data
  {
  switch(MLM)
    {
    case -1: return std::string("none"); break;
    case  0: return std::string("IEEE Little Endian"); break;
    case  1: return std::string("IEEE Big Endian");    break;
    case  2: return std::string("VAX D-float");        break;
    case  3: return std::string("VAX G-float");        break;
    case  4: return std::string("Cray");               break;
    }
  return std::string("Unknown");
  }


std::string MatLab4Tag::PType() const

        // Input 		ML4T    : MAT version 4 tag (this)
        // Output               string  : String labeling the type
	//				  of data that is stored
  {
  switch(MLP)
    {
    case -1: return std::string("none"); break;
    case  0: return std::string("64 bit, double precision"); break;
    case  1: return std::string("32 bit, single precision"); break;
    case  2: return std::string("32 bit, signed integer");   break;
    case  3: return std::string("16 bit, signed integer");   break;
    case  4: return std::string("16 bit, unsigned integer"); break;
    case  5: return std::string("8 bit,  unsigned integer"); break;
    }
  return std::string("Unknown");
  }


std::string MatLab4Tag::TType() const

        // Input 		ML4T    : MAT version 4 tag (this)
        // Output               string  : String labeling the type
	//				  of matrix that is stored
  {
  switch(MLT)
    {
    case -1: return std::string("none"); break;
    case  0: return std::string("Full Matrix");   break;
    case  1: return std::string("Text Matrix");   break;
    case  2: return std::string("Sparse Matrix"); break;
    }
  return std::string("Unknown");
  }


int MatLab4Tag::ElemBytes() const

        // Input 		ML4T    : MAT version 4 tag (this)
        // Output               int     : Number of bytes in an array element
  {
  switch(MLP)
    {
    case  0: return 8; break;
    case  1: return 4; break;
    case  2: return 4; break;
    case  3: return 2; break;
    case  4: return 2; break;
    case  5: return 1; break;
    default: return 0; break;
    }
  return 0;
  }


int MatLab4Tag::Bytes() const { return mrows*ncols*ElemBytes(); }

        // Input 		ML4T    : MAT version 4 tag (this)
        // Output               NBytes	: Number of bytes (in data)
	//				: This will be equivalent to 
	//				  (mrows)x(cols)x(sizeof data)


// ____________________________________________________________________________
// C                  MATLAB MAT 4 Tag Binary Output Functions
// ____________________________________________________________________________

/* These are the functions which will output a tag in MATLAB ".mat" binary
   format, version 4.  Note that in a valid MATLAB MAT file, this tag will
   be associated with a "Data Element" or a "Sub-Element".
 
                Input           ML4H    : MAT version 4 tag (this)
                                fp      : Output file stream
                Output          void    : Output file stream is modified
                                          by having the ML4H written to it
                Note                    : No care is taken to insure the
                                          tag is written tothe top of the
                                          file (as it should be in MATLAB)   */

int MatLab4Tag::write(std::fstream& fp) const                                       
  {
  fp.write((char*)&type,   sizeof(int32_t));	// Write in computer type
  fp.write((char*)&mrows,  sizeof(int32_t));	// Write in array rows
  fp.write((char*)&ncols,  sizeof(int32_t));	// Write in array columns
  fp.write((char*)&cmplxf, sizeof(int32_t));	// Write in real/complex flag
  fp.write((char*)&nlen,   sizeof(int32_t));	// Write in name length
  return 1;
  }


int MatLab4Tag::Size() const { return 20; }	// Bytes used in a file


// ____________________________________________________________________________
// D                MATLAB MAT 4 Tag Binary Input Functions
// ____________________________________________________________________________
 
/* These are the functions which will input a tag in MATLAB ".mat" binary
   format, version 4.  Note that in a valid MATLAB MAT file, this tag will
   associated with a "Data Element", a combination of Tag + Data 

		Input           ML4T    : MAT version 4 tag (this)
           			fp	: Input file stream
                                warn	: Warning output level
                                                0 = no warnings
                                                1 = warnings
                                               >1 = fatal warnings
		Output		void    : ML4T is set from values found
					  in the input stream
		Note			: No care is taken to insure the
					  tag is read from the top of a data
					  element (as it should be in MATLAB)*/
 
int MatLab4Tag::read(std::fstream& fp, int warn)
  {
  fp.read((char*)&type,   sizeof(int32_t));	// Read in computer type
  fp.read((char*)&mrows,  sizeof(int32_t));	// Read in array rows
  fp.read((char*)&ncols,  sizeof(int32_t));	// Read in array columns
  fp.read((char*)&cmplxf, sizeof(int32_t));	// Read in real/complex flag
  fp.read((char*)&nlen,   sizeof(int32_t));	// Read in name length
  if(!fp)
    {
    MLT4error(1, 1);				// Problems with filestream
    if(warn > 1) MLT4fatality(3);		// Can't read header
    else         MLT4error(3,1);
    return 0;
    }
  Type2MOPT();					// Set MOPT from type
  if((!MLM && WeRBigEnd()) ||( MLM && !WeRBigEnd()))// If read data is little
    {						// endian yet we are big
    Swap(mrows);				// endian (or vice-versa)
    Swap(ncols);				// we must do byte swapping
    Swap(cmplxf);
    Swap(nlen);
    }
  return 1;
  } 


// ____________________________________________________________________________
// E                  MATLAB MAT 4 Tag ASCII Output Functions
// ____________________________________________________________________________
 

void MatLab4Tag::print(std::ostream& ostr, int hf) const

	// Input 		ML4T    : MAT version 4 tag (this)
	// 			ostr	: An output stream
	//			hf 	: Flag to print header   
	//                      bf      : Flag to print tag chars
	// Output		ostr	: Output stream with MATLAB MAT file
	//				  tag version 4 info placed into it
	// Note				: This information exists at the
	//				  beginning of each data element

  {
  if(hf) ostr << "\n\t\tHeader";
  ostr << "\n\t\t  Rows:          " << mrows;	// Number of data rows
  ostr << "\n\t\t  Columns:       " << ncols;	// Number of data columns
  ostr << "\n\t\t  R/C Flag:      "; 		// Real/Complex flag
  if(cmplxf) ostr << "Complex";
  else       ostr << "Real";
  ostr << "\n\t\t  Binary Format: " << PType();	// The binary format type
  ostr << "\n\t\t  Byte Ordering: " << MType();	// The number byte order
  ostr << "\n\t\t  Data Type:     " << TType();	// The data array type
  ostr << "\n\t\t  Name Len:      " << nlen;	// Variable name length
  ostr.flush();
  }


std::ostream& operator<< (std::ostream& ostr, const MatLab4Tag& MLT4)

        // Input                ostr	: An output stream
	// 			ML4T    : MAT version 4 tag
        // Output               ostr	: The output stream modified by
        //                                the tag parameters

  { MLT4.print(ostr); return ostr; }


#endif							// ML4Tag.cc
