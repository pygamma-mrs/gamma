/* ML5Imags.cc **************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A                           **
**                                                                      **
**     MATLAB Imaginaries Array (MAT Version 5)     Implementation	**
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
** MATLAB imags array sub-element in MAT version 5 files.  MATLAB 	**
** is a commercial product and can be purchaced from The MathWorks, Inc.**
** This module deals with binary MAT files version 5 circa 1999.      	**
**                                                                      **
*************************************************************************/

#ifndef _ML5IM_CC_				// Is file already included?
#define _ML5IM_CC_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)					// Using the GNU compiler?
#    pragma implementation				// This is the implementation
#endif

#include <GamIO/ML5Imags.h>			// Include the header file
#include <GamIO/ML5SubE.h>			// Include base class header
#include <GamIO/ML5DA.h>			// Include dimensions array SE
#include <Basics/Gutils.h>                      // Include GAMMA errors
#include <string>				// Include libstdc++ strings
#include <stdio.h>


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i          MATLAB MAT 5 Imags Array Sub-Element Error Handling
// ____________________________________________________________________________

/*              Input           MLIM5   : MAT version 5 sub-element
                                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */  

void MatLab5Im::MLIM5error(int eidx, int noret)
  {                                                     
  std::string hdr("MATLAB MAT V5 Imags Array Sub-Element");
  switch(eidx)
    {
    case 0: GAMMAerror(hdr, 0, noret); break;   // Program Aborting        (0)
    case 1: GAMMAerror(hdr, "End of File Imached Before Data Found", noret);
      break;                                     //                        (1)
    case 10: GAMMAerror(hdr, "Cannot Imad Data Element Header", noret); // (10)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }

 
void MatLab5Im::MLIM5error(int eidx, const std::string& pname, int noret)
  {
  std::string hdr("MATLAB MAT V5 Imags Array Sub-Element");
  std::string msg;
  switch(eidx)
    {
    case 1:  GAMMAerror(hdr,   1, pname, noret); break;	// File Problems   (1)
    case 2:  GAMMAerror(hdr,   2, pname, noret); break;	// !Imad SingPar   (2)
    default: GAMMAerror(hdr, -11, pname, noret); break;	// Unknown Error   (-1)
    }
  }
     
 
volatile void MatLab5Im::MLIM5fatality(int eidx)
  {
  MLIM5error(eidx, 1);
  if(eidx) MLIM5error(0);
  GAMMAfatal();					// Clean exit from program
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------


// ____________________________________________________________________________
// A                  MATLAB MAT 5 Imags Array Sub-Element Constructors
// ____________________________________________________________________________


MatLab5Im::MatLab5Im() : MatLab5SE() { }
 
MatLab5Im::MatLab5Im(const matrix& mx) : MatLab5SE()
  { MLTag = MatLab5Tag(9, mx.rows()*mx.cols()); }

MatLab5Im::~MatLab5Im() { }
 

// ____________________________________________________________________________
// B          MATLAB MAT 5 Imags Array Sub-Element Access Functions
// ____________________________________________________________________________
   

// ____________________________________________________________________________
// C       MATLAB MAT 5 Imags Array Sub-Element Binary Output Functions
// ____________________________________________________________________________

 
/* These are the functions which will output the sub-element in MATLAB ".mat" 
   binary format, version 5.  This is done with a pair {tag,data}.  Both the
   tag & data will be written directly from these functions, bypassing the 
   associated class structures which have no knowledge of GAMMA array types. 
  
   Remember, the MATLAB MAT version 5 tag for this sub-element contains 
   8-bytes with two important settings: 1.) data type & 2.) # of bytes. 
   The data type is placed in the first 4 bytes and the number of bytes 
   placed into the last 4 bytes (for non-compressed data elements which we 
   will always output...).  In this case, we output the elements as doubles
   so that the data type should be MATLAB type miDOUBLE, i.e. value 9. 
  
                Input           ML5Im   : MAT version 5 imags SE (this)
                                fp      : Output file stream 
                                mx      : GAMMA matrix 
                                rv      : GAMMA row vector 
                                cv      : GAMMA column vector 
                Output          void    : Output file stream is modified 
                                          by having the ML5Im written to it 
                Note                    : No care taken herein to insure the 
                                          sub-element is properly written.   */ 
  
int MatLab5Im::write(std::fstream& fp, const matrix& mx) const
  {
  int32_t TB = 9;                                  // MATLAB data type miDOUBLE
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the data type  (TAG)
  TB = mx.rows()*mx.cols()*sizeof(double);      // Number of bytes to output
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the # of bytes (TAG)
  int i,j;
  double d;
  for(j=0; j<mx.cols(); j++)                    // Output the imag elements of
    for(i=0; i<mx.rows(); i++)                  // the matrix column by column
      {
      d = mx.getIm(i,j);
      fp.write((char*)&d, sizeof(double));
      }
  return 1;
  }

  
int MatLab5Im::write(std::fstream& fp, const row_vector& rv) const                      
  {
  int32_t TB = 9;                                  // MATLAB data type miDOUBLE
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the data type  (TAG)
  TB = rv.size()*sizeof(double);                // Number of bytes to output
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the # of bytes (TAG)
  double d;
  for(int i=0; i<rv.size(); i++)		// Output the vector elements
    {
    d = rv.getIm(i);
    fp.write((char*)&d, sizeof(double));
    } 
  return 1;
  } 
 
  
int MatLab5Im::write(std::fstream& fp, const col_vector& cv) const           
  {
  int32_t TB = 9;                                  // MATLAB data type miDOUBLE
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the data type  (TAG)
  TB = cv.size()*sizeof(double);                // Number of bytes to output
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the # of bytes (TAG)
  double d;
  for(int i=0; i<cv.size(); i++)		// Output the vector elements
    {
    d = cv.getIm(i);
    fp.write((char*)&d, sizeof(double));
    }
  return 1;
  }
 

/* These are the functions which will return the number of bytes that are
   written upon output the sub-element in MATLAB ".mat" binary format, V.5
   For the Imaginary data, we will the the Tag (8 bytes), followed by the
   array/vector elements (NE*double), followed by any additional padding to
   fill the sub-element up to an eight byte boundary.                        */

int MatLab5Im::Size(const matrix& mx) const
  {
  int NB =  MLTag.Size();
  NB += mx.rows()*mx.cols()*sizeof(double);
  return NB;
  }


   
// ____________________________________________________________________________
// D         MATLAB MAT 5 Imags Array Sub-Element Binary Input Functions
// ____________________________________________________________________________

/* These are the functions which will input an Imags Array sub-element in MATLAB
   ".mat" binary format, version 5.  Note that in a valid MATLAB MAT file, a 
   "Imags Array Sub-Element" is a combination of Tag + Data
 
                Input           ML5Im	: MAT version 5 imags array (this)
                                fp      : Input file stream
                                warn    : Warning output level
                                                0 = no warnings
                                                1 = warnings
                                               >1 = fatal warnings
                Output          void    : ML5Im is set from values found
                                          in the input stream
                Note                    : No care is taken to insure the
                                          imags array is read from the top of a
                                          data element (as it should be in
					  MATLAB)                            */


int MatLab5Im::fread(std::fstream& fp, int bigend, int warn)
  {
  if(!MLTag.read(fp, bigend, (warn)?1:0)) 	// Read sub-element tag
    {                                           // If this has problems, then
    if(warn)                                    // we must deal with them
      {                                         // appropriately now
      if(warn==1) MLIM5error(10, 1);
      else        MLIM5fatality(10);
      }  
    return 0;
    }
  if(!MLTag.Compressed)				// If the data element isn't
    fp.seekp(MLTag.Bytes(), std::ios::cur);		// compressed, skip to data end
  return 1;					// (if compressed, at end now)
  }

int MatLab5Im::read(std::fstream& fp, int bigend, matrix& mx, int warn)
  {
  if(!MLTag.read(fp, bigend, (warn)?1:0)) 	// Read sub-element tag
    {                                           // If this has problems, then
    if(warn)                                    // we must deal with them
      {                                         // appropriately now
      if(warn==1) MLIM5error(10, 1);
      else        MLIM5fatality(10);
      }  
    return 0;
    }
  int nr = mx.rows();				// Get number of rows
  int nc = mx.cols();;				// Get number of columns
  int swap = 0;
  if(bigend != int(WeRBigEnd())) swap++;
  double d;
  char c;
  int i, j, nb;
  for(j=0, nb=0; j<nc; j++)
    {
    for(i=0; i<nr; i++, nb++)
      {
      switch(MLTag.DType)
        {
        case 2:  fp.read(&c,sizeof(char)); d=double(c); break;	// miUINT8
	case 9:							// miDOUBLE
        default: fp.read((char*)&d,sizeof(double)); 
                 if(swap) Swap(d);			break;	// GAMMA 
        }
      mx.put(mx.get(i,j)+complex(0,d),i,j);
      }
    }
  return 1;
  }

 
 
// ____________________________________________________________________________
// E                MATLAB MAT 5 Imags Array Sub-Element ASCII Output Functions
// ____________________________________________________________________________

                                                                                
void MatLab5Im::print(std::ostream& ostr) const

        // Input                ML5Im	: MAT version 5 sub-element (this)
        //                      ostr    : An output stream
        // Output               ostr    : Output stream with MATLAB MAT file
        //                                sub-element version 5 info placed
	//				  into it
        // Note                         : This information exists at the
        //                                beginning of each sub-element

  {
  ostr << "\n\t\tSE Imags Array  ";		// Print the header info
  MLTag.print(ostr, 0, 0);			// Output 8-byte tag
  }
 
std::ostream& operator<< (std::ostream& ostr, const MatLab5Im& ML5Im)
 
        // Input                ostr    : An output stream
        //                      ML5Im	: MAT version 5 sub-element
        // Output               ostr    : The output stream modified by
        //                                the sub-element parameters

  { ML5Im.print(ostr); return ostr; }


#endif							// ML5Imags.cc
