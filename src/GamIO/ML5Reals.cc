/* ML5Reals.cc **************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A                           **
**                                                                      **
**     MATLAB Reals Array (MAT Version 5)        Implementation		**
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
** MATLAB reals array sub-element in MAT version 5 files.  MATLAB 	**
** is a commercial product and can be purchaced from The MathWorks, Inc.**
** This module deals with binary MAT files version 5 circa 1999.      	**
**                                                                      **
*************************************************************************/

#ifndef _ML5REALS_CC_				// Is file already included?
#define _ML5REALS_CC_ 1				// If no, then remember it

#  if defined(GAMPRAGMA)					// Using the GNU compiler?
#    pragma implementation				// This is the implementation
#endif

#include <GamIO/ML5Reals.h>			// Include the header file
#include <GamIO/ML5SubE.h>			// Include base class header
#include <GamIO/ML5DA.h>			// Include dimensions array SE
#include <Basics/Gutils.h>                      // Include GAMMA errors
#include <string>				// Include libstdc++ strings
#include <Matrix/matrix.h>			// Include GAMMA matrices

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i          MATLAB MAT 5 Reals Array Sub-Element Error Handling
// ____________________________________________________________________________

/*              Input           MLRE5   : MAT version 5 sub-element
                                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */  

void MatLab5Re::MLRE5error(int eidx, int noret)
  {                                                     
  std::string hdr("MATLAB MAT V5 Reals Array Sub-Element");
  switch(eidx)
    {
    case 0: GAMMAerror(hdr, 0, noret); break;   // Program Aborting        (0)
    case 1: GAMMAerror(hdr, "End of File Reached Before Data Found", noret);
      break;                                     //                        (1)
    case 10: GAMMAerror(hdr, "Cannot Read Data Element Header", noret); // (10)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }

 
void MatLab5Re::MLRE5error(int eidx, const std::string& pname, int noret)
  {
  std::string hdr("MATLAB MAT V5 Reals Array Sub-Element");
  std::string msg;
  switch(eidx)
    {
    case 1:  GAMMAerror(hdr,   1, pname, noret); break;	// File Problems   (1)
    case 2:  GAMMAerror(hdr,   2, pname, noret); break;	// !Read SingPar   (2)
    default: GAMMAerror(hdr, -11, pname, noret); break;	// Unknown Error   (-1)
    }
  }
     
 
volatile void MatLab5Re::MLRE5fatality(int eidx)
  {
  MLRE5error(eidx, 1);
  if(eidx) MLRE5error(0);
  GAMMAfatal();					// Clean exit from program
  }
 

// ____________________________________________________________________________
// A                  MATLAB MAT 5 Reals Array Sub-Element Constructors
// ____________________________________________________________________________


MatLab5Re::MatLab5Re() : MatLab5SE() { }

MatLab5Re::MatLab5Re(const matrix& mx) : MatLab5SE()
  { MLTag = MatLab5Tag(9, mx.rows()*mx.cols()); }

MatLab5Re::~MatLab5Re() { }
 

// ____________________________________________________________________________
// B                MATLAB MAT 5 Reals Array Sub-Element Access Functions
// ____________________________________________________________________________
   

// ____________________________________________________________________________
// C     MATLAB MAT 5 Reals Array Sub-Element Binary Output Functions
// ____________________________________________________________________________

/* These are the functions which will output this sub-element in MATLAB ".mat"
   binary format, version 5.  This is done with a pair {tag,data}.  Both the
   tag & data will be written directly from these functions, bypassing the
   associated class structures which have no knowledge of GAMMA array types.
   
   Remember, the MATLAB MAT version 5 tag for this sub-element contains
   8-bytes with two important settings: 1.) data type & 2.) # of bytes.
   The data type is placed in the first 4 bytes and the number of bytes
   placed into the last 4 bytes (for non-compressed data elements which we
   will always output...).  In this case, we output the elements as doubles
   so that the data type should be MATLAB type miDOUBLE, i.e. value 9.
 
                Input           ML5Re	: MAT version 5 sub-element (this)
                                fp      : Output file stream
				mx      : GAMMA matrix
				rv	: GAMMA row vector
				cv	: GAMMA column vector
                Output          void    : Output file stream is modified
                                          by having the ML5Re written to it
		Note			: Since we write doubles (8 bytes)
					  per element we always end on an 8 
					  byte boundary.
                Note                    : No care taken herein to insure the
                                          sub-element is properly written.   */
 
int MatLab5Re::write(std::fstream& fp, const matrix& mx) const	
  {
  int32_t TB = 9;					// MATLAB data type miDOUBLE
  fp.write((char*)&TB, sizeof(int32_t));			// Output the data type  (TAG)
  TB = mx.rows()*mx.cols()*sizeof(double);	// Number of mx bytes to output
  fp.write((char*)&TB, sizeof(int32_t));			// Output the # of bytes (TAG)
  int i,j;
  double d;
  for(j=0; j<mx.cols(); j++)			// Output the real elements of
    for(i=0; i<mx.rows(); i++) 			// the matrix column by column
      {
      d = mx.getRe(i,j);
      fp.write((char*)&d, sizeof(double));
      }
  return 1;
  }
  
int MatLab5Re::write(std::fstream& fp, const row_vector& rv) const                      
  {
  int32_t TB = 9;                                  // MATLAB data type miDOUBLE
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the data type  (TAG)
  TB = rv.size()*sizeof(double);                // Number of rv bytes to output
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the # of bytes (TAG)
  double d;
  for(int i=0; i<rv.size(); i++)		// Output the vector elements
    {
    d = rv.getRe(i);
    fp.write((char*)&d, sizeof(double));
    }
  return 1;
  } 
 
  
int MatLab5Re::write(std::fstream& fp, const col_vector& cv) const           
  {
  int32_t TB = 9;                                  // MATLAB data type miDOUBLE
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the data type  (TAG)
  TB = cv.size()*sizeof(double);                // Number of cv bytes to output
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the # of bytes (TAG)
  double d;
  for(int i=0; i<cv.size(); i++) 		// Output the vector elements
    {
    d = cv.getRe(i);
    fp.write((char*)&d, sizeof(double));
    }
  return 1;
  }


/* These are the functions which will return the number of bytes that are
   written upon output the sub-element in MATLAB ".mat" binary format, V.5   */

int MatLab5Re::Size(const matrix& mx) const
  {
  int NB =  MLTag.Size();
  NB += mx.rows()*mx.cols()*sizeof(double);
  return NB;
  }

// ____________________________________________________________________________
// D       MATLAB MAT 5 Reals Array Sub-Element Binary Input Functions
// ____________________________________________________________________________

/* These are the functions which will input a reals array sub-element in MATLAB
   ".mat" binary format, version 5.  In a valid MATLAB MAT file, the reals 
   array sub-element is a combination of {tag,data}, and both will be read in
   these funcitons.  
   Remember, the MATLAB MAT version 5 tag for this sub-element contains
   8-bytes with two important settings: 1.) data type & 2.) # of bytes.
   The data type resides in the first 4 bytes and the number of bytes
   resides in the last 4 bytes (for non-compressed data elements, which should
   be the case for arrays). 
   It would be splendid if the array elements themselves were stored as doubles
   (MATLAB type miDOUBLE, i.e. value 9), but such will only be certain if the
   MAT file was generated from GAMMA. Instead, we may have to type cast the
   values and byte swap as needed.
   Finally, class MatLab5Tag does handle reading of the tag but class 
   MatLab5SE does NOT handle reading of the data.  That is because here we want
   the data placed into a GAMMA array, and that same array we will wish to 
   share with the imaginaries array sub-element if necessary.  So we handle the
   data input directly.
 
                Input           ML5Re	: MAT version 5 sub-element (this)
                                fp      : Output file stream
				mx      : GAMMA matrix
				rv	: GAMMA row vector
				cv	: GAMMA column vector
                Output          void    : Output file stream is modified
                                          by having the ML5Re written to it
                Note                    : No care taken herein to insure the
 
                Input           ML5Re	: MAT version 5 reals array (this)
                                fp      : Input file stream
                                warn    : Warning output level
                                                0 = no warnings
                                                1 = warnings
                                               >1 = fatal warnings
                Output          void    : ML5Re is set from values found
                                          in the input stream
                Note                    : No care is taken to insure the
                                          reals array is read from the top of a
                                          data element (as it should be in
					  MATLAB)                            */

int MatLab5Re::fread(std::fstream& fp, int bigend, int warn)
  {
  if(!MLTag.read(fp, bigend, (warn)?1:0)) 	// Read sub-element tag
    {                                           // If this has problems, then
    if(warn)                                    // we must deal with them
      {                                         // appropriately now
      if(warn==1) MLRE5error(10, 1);
      else        MLRE5fatality(10);
      }  
    return 0;
    }
  if(!MLTag.Compressed)				// If the data element isn't
    fp.seekp(MLTag.Bytes(), std::ios::cur);		// compressed, skip to data end
  return 1;					// (if compressed, at end now)
  }


matrix MatLab5Re::read(std::fstream& fp, int bigend, const MatLab5DA& DA, int warn)
  {
  if(!MLTag.read(fp, bigend, (warn)?1:0)) 	// Read sub-element tag
    {                                           // If this has problems, then
    if(warn)                                    // we must deal with them
      {                                         // appropriately now
      if(warn==1) MLRE5error(10, 1);
      else        MLRE5fatality(10);
      }  
    return 0;
    }
  int nr = DA.Dim(0);				// Get number of rows
  int nc = DA.Dim(1);				// Get number of columns
  matrix mx(nr, nc);				// Construct starting array
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
      mx.put(d,i,j);
      }
    }
  return mx;
  }

 
// ____________________________________________________________________________
// E        MATLAB MAT 5 Reals Array Sub-Element ASCII Output Functions
// ____________________________________________________________________________

                                                                                
void MatLab5Re::print(std::ostream& ostr) const

        // Input                ML5Re	: MAT version 5 sub-element (this)
        //                      ostr    : An output stream
        // Output               ostr    : Output stream with MATLAB MAT file
        //                                sub-element version 5 info placed
	//				  into it
        // Note                         : This information exists at the
        //                                beginning of each sub-element

 {
 ostr << "\n\t\tSE Reals Array  ";		// Print the header info
 MLTag.print(ostr, 0, 0);			// Output 8-byte tag
// ostr << "\n\t\t  Name Chars:   " << NC; 	// # of characters in name
// ostr << "\n\t\t  Reals Array:   " << MxName; 	// Matrix name
// if(NC && data)
//   {
//   ostr << "\n\t\t  Reals Array:   " ; 		// Output the reals array
//   for(int i=0; i<NC; i++) ostr << data[i];
//   }
 }
 
std::ostream& operator<< (std::ostream& ostr, const MatLab5Re& ML5Re)
 
        // Input                ostr    : An output stream
        //                      ML5Re	: MAT version 5 sub-element
        // Output               ostr    : The output stream modified by
        //                                the sub-element parameters

  { ML5Re.print(ostr); return ostr; }

#endif							// _ML5REALS_CC_
