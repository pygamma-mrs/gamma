/* ML5DA.cc *****************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A                           **
**                                                                      **
**     MATLAB Dimensions Array (MAT Version 5)        Implementation	**
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
** MATLAB dimensions array sub-element in MAT version 5 files.  MATLAB 	**
** is a commercial product and can be purchaced from The MathWorks, Inc.**
** This module deals with binary MAT files version 5 circa 1999.      	**
**                                                                      **
*************************************************************************/

#ifndef _ML5DA_CC_				// Is file already included?
#define _ML5DA_CC_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)					// Using the GNU compiler?
#    pragma implementation				// This is the implementation
#endif

#include <GamIO/ML5DA.h>			// Include the header file
#include <GamIO/ML5SubE.h>			// Include base class header
#include <Basics/Gutils.h>                      // Include GAMMA errors
#include <string>				// Include libstdc++ strings
#include <stdio.h>


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i          MATLAB MAT 5 Dimensions Array Sub-Element Error Handling
// ____________________________________________________________________________

/*              Input           MLDA5   : MAT version 5 sub-element
                                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */  

void MatLab5DA::MLDA5error(int eidx, int noret)
  {                                                     
  std::string hdr("MATLAB MAT V5 Dimensions Array Sub-Element");
  switch(eidx)
    {
    case 0: GAMMAerror(hdr, 0, noret); break;   // Program Aborting        (0)
    case 1: GAMMAerror(hdr, "End of File Reached Before Data Found", noret);
      break;                                     //                        (1)
    case 10: GAMMAerror(hdr, "Cannot Read Data Element Header", noret);	 //(10)
    case 11: GAMMAerror(hdr, "Unreasonable Number of Dimensions",noret); //(11)
    case 12: GAMMAerror(hdr, "Unreasonable Array Dimension", noret);     //(12)
    case 20: GAMMAerror(hdr, "Error During Input From File", noret);     //(20)

    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }

 
void MatLab5DA::MLDA5error(int eidx, const std::string& pname, int noret)
  {
  std::string hdr("MATLAB MAT V5 Dimensions Array Sub-Element");
  std::string msg;
  switch(eidx)
    {
    case 1:  GAMMAerror(hdr,   1, pname, noret); break;	// File Problems   (1)
    case 2:  GAMMAerror(hdr,   2, pname, noret); break;	// !Read SingPar   (2)
    default: GAMMAerror(hdr, -11, pname, noret); break;	// Unknown Error   (-1)
    }
  }
     
 
volatile void MatLab5DA::MLDA5fatality(int eidx)
  {
  MLDA5error(eidx, 1);
  if(eidx) MLDA5error(0);
  GAMMAfatal();					// Clean exit from program
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------


// ____________________________________________________________________________
// A                  MATLAB MAT 5 Dimensions Array Sub-Element Constructors
// ____________________________________________________________________________


MatLab5DA::MatLab5DA() : MatLab5SE() { ND = 0; Sizes = NULL; }
MatLab5DA::MatLab5DA(const matrix& mx) : MatLab5SE()

	// Input		ML5DA   : MAT version 5 dimensions array (this)
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

  {
  MLTag  = MatLab5Tag(5,8);		// Set up the tag
  ND = 2;
  Sizes = new int[ND];
  Sizes[0] = mx.rows();
  Sizes[1] = mx.cols();
  }

MatLab5DA::~MatLab5DA() { if(Sizes) delete [] Sizes; }
 

// ____________________________________________________________________________
// B        MATLAB MAT 5 Dimensions Array Sub-Element Access Functions
// ____________________________________________________________________________
   
int MatLab5DA::Dims()     const { return ND; }
int MatLab5DA::Dim(int i) const
  { 
  if(i<0 || i>ND-1) return 0;
  return Sizes[i];
  }

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
                                 

int MatLab5DA::write(std::fstream& fp) const
  {
  int TF = MLTag.write(fp);			// Write the Tag (uncompressed)
  int32_t LI;					// Use long integer for output
  for(int i=0; i<ND; i++)			// Now write out all the array
    {						// dimensions as long integers
    LI = Sizes[i];
    fp.write((char*)&LI, sizeof(int32_t));
    }
  return TF;
  }


int MatLab5DA::write(std::fstream& fp, const matrix& mx) const
  {
  int32_t TB = 5;                                  // MATLAB data type miDOUBLE
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the data type  (TAG)
  TB = 2*sizeof(int32_t);				// Number of bytes to output
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the # of bytes (TAG)
  TB = mx.rows();				// The 1st mx dimension
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the 1st mx dimension
  TB = mx.cols();				// The 2nd mx dimension
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the 2nd mx dimension
  return 1;
  }


int MatLab5DA::write(std::fstream& fp, const row_vector& rv) const
  {
  MLTag.write(fp);                              // First write the tag
  int32_t TB = 5;                                  // MATLAB data type miDOUBLE
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the data type  (TAG)
  TB = 2*sizeof(int32_t);				// Number of bytes to output
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the # of bytes (TAG)
  TB = 1;					// The vector "row" dimension
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the row dimension
  TB = rv.size();				// The vector "column" dimen.
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the vector size
  return 1;
  }


int MatLab5DA::write(std::fstream& fp, const col_vector& cv) const
  {
  MLTag.write(fp);                              // First write the tag
  int32_t TB = 5;                                  // MATLAB data type miDOUBLE
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the data type  (TAG)
  TB = 2*sizeof(int32_t);				// Number of bytes to output
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the # of bytes (TAG)
  TB = cv.size();				// The vector "row" dimension
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the vector size
  TB = 1;					// The vector "column" dimen.
  fp.write((char*)&TB, sizeof(int32_t));                  // Output the col dimension
  return 1;
  }

 
/* These are the functions which will return the number of bytes that are 
   written upon output this sub-element in MATLAB ".mat" binary format, V5   */ 
 
int MatLab5DA::Size(const matrix& mx) const { return 4*sizeof(int32_t); }


// ____________________________________________________________________________
// D     MATLAB MAT 5 Dimensions Array Sub-Element Binary Input Functions
// ____________________________________________________________________________

/* These are the functions which will input a sub-element in MATLAB ".mat"
   binary format, version 5.  Note that in a valid MATLAB MAT file, a 
   "Dimensions Array Sub-Element" is a combination of Tag + Data
 
                Input           ML5DA   : MAT version 5 dim. array (this)
                                fp      : Input file stream
                                warn     : Warning output level
                                                0 = no warnings
                                                1 = warnings
                                               >1 = fatal warnings
                Output          void    : ML5DA is set from values found
                                          in the input stream
                Note                    : No care is taken to insure the
                                          dim. array is read from the top of
 					  a data element (as it should be 
					  in MATLAB)                         */

int MatLab5DA::read(std::fstream& fp, int bigend, int warn)
  {
int MaxDims=5;
int MaxDim=524288;
  if(!MatLab5SE::read(fp, bigend, (warn)?1:0))	// Read sub-element
    {						// If this has problems, then
    if(warn)					// we must deal with them
      {						// appropriately now
      if(warn==1) MLDA5error(10, 1);
      else        MLDA5fatality(10);
      }
    return 0;
    }
  ND = MLTag.Bytes()/4;				// Set number of dimensions
  if(ND < 2 || ND>MaxDims)
    {
    if(warn)
      {
      MLDA5error(11, 1);
      if(warn > 1) MLDA5fatality(20);
      else         MLDA5error(20);
      }
    }
  if(Sizes) delete [] Sizes;			// Delete dimensions array
  Sizes = new int[ND];				// Set new dimensions array
  longchars LC;					// Read long integers
  int swap = 0;
  if(bigend != int(WeRBigEnd())) swap++;
  for(int i=0; i<ND; i++)
    { 
    LC.chars[0] = MLData[0+4*i];
    LC.chars[1] = MLData[1+4*i];
    LC.chars[2] = MLData[2+4*i];	
    LC.chars[3] = MLData[3+4*i];
    if(swap) Swap(LC.longval);
    Sizes[i] = int(LC.longval);
    if(Sizes[i]<1 || Sizes[i]>MaxDim)
      {
      if(warn)
        {
        MLDA5error(12, 1);
        if(warn > 1) MLDA5fatality(20);
        else         MLDA5error(20);
        }
      }
    } 
  return 1;
  }

 
// ____________________________________________________________________________
// E      MATLAB MAT 5 Dimensions Array Sub-Element ASCII Output Functions
// ____________________________________________________________________________

                                                                                
void MatLab5DA::print(std::ostream& ostr) const

        // Input                ML5DA	: MAT version 5 sub-element (this)
        //                      ostr    : An output stream
        // Output               ostr    : Output stream with MATLAB MAT file
        //                                sub-element version 5 info placed
	//				  into it
        // Note                         : This information exists at the
        //                                beginning of each sub-element

 {
 ostr << "\n\t\tSE Dimen. Array ";		// Output a header
 MLTag.print(ostr, 0);				// Print sub-element tag
 ostr << "\n\t\t  Dimensions:   " << ND; 	// # of Dimensions
 if(ND)
   {
//   longchars LC;
   ostr << "   ("; 
   for(int i=0; i<ND; i++)
     { 						// We clip out the dimensions
/*
     LC.chars[0] = MLData[0+4*i];		// from the data.
     LC.chars[1] = MLData[1+4*i];
     LC.chars[2] = MLData[2+4*i];	
     LC.chars[3] = MLData[3+4*i];
     ostr << LC.longval;
*/
     ostr << Sizes[i];
     if(i<ND-1) ostr << " x ";
     }
   ostr << ")"; 
   }
 ostr.flush();
 }
 
std::ostream& operator<< (std::ostream& ostr, const MatLab5DA& ML5DA)
 
        // Input                ostr    : An output stream
        //                      ML5DA	: MAT version 5 sub-element
        // Output               ostr    : The output stream modified by
        //                                the sub-element parameters

  { ML5DA.print(ostr); return ostr; }


#endif							// ML5DA.cc
