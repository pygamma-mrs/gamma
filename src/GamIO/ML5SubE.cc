/* ML5SubE.cc ***************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A                           **
**                                                                      **
**	MATLAB Sub-Element (MAT Version 5)        Implementation	**
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
** MATLAB sub-elements in MAT version 5 files.  MATLAB is a commercial	**
** product and can be purchaced from The MathWorks, Inc.  This  module  **
** deals with binary MAT files version 5 circa 1999.                    **
**                                                                      **
*************************************************************************/

#ifndef _ML5SUBE_CC_				// Is file already included?
#define _ML5SUBE_CC_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)					// Using the GNU compiler?
#    pragma implementation				// This is the implementation
#endif

#include <GamIO/ML5SubE.h>			// Include the header file
#include <Basics/Gutils.h>                      // Include GAMMA errors
#include <string>				// Include libstdc++ strings
#include <stdio.h>


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                 MATLAB MAT 5 Sub-Element Error Handling
// ____________________________________________________________________________

/*              Input           MLSE5   : MAT version 5 sub-element
                                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */  

void MatLab5SE::MLSE5error(int eidx, int noret)
  {                                                     
  std::string hdr("MATLAB MAT V5 Sub-Element");
  switch(eidx)
    {
    case 0: GAMMAerror(hdr, 0, noret); break;   // Program Aborting        (0)
    case 1: GAMMAerror(hdr, "End of File Reached Before Data Found", noret);
      break;                                     //                        (1)
    case 10: GAMMAerror(hdr, "Cannot Read Data Element Header", noret); // (10)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }

 
void MatLab5SE::MLSE5error(int eidx, const std::string& pname, int noret)
  {
  std::string hdr("MATLAB MAT V5 Sub-Element");
  std::string msg;
  switch(eidx)
    {
    case 1:  GAMMAerror(hdr,   1, pname, noret); break;	// File Problems   (1)
    case 2:  GAMMAerror(hdr,   2, pname, noret); break;	// !Read SingPar   (2)
    default: GAMMAerror(hdr, -11, pname, noret); break;	// Unknown Error   (-1)
    }
  }
     
 
volatile void MatLab5SE::MLSE5fatality(int eidx)
  {
  MLSE5error(eidx, 1);
  if(eidx) MLSE5error(0);
  GAMMAfatal();					// Clean exit from program
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------


// ____________________________________________________________________________
// A                  MATLAB MAT 5 Sub-Element Constructors
// ____________________________________________________________________________


MatLab5SE::MatLab5SE() { MLData = NULL; }

MatLab5SE::~MatLab5SE() { if(MLData) delete [] MLData; }


 
// ____________________________________________________________________________
// B                MATLAB MAT 5 Sub-Element Access Functions
// ____________________________________________________________________________
   

// ____________________________________________________________________________
// C            MATLAB MAT 5 Sub-Element Binary Output Functions
// ____________________________________________________________________________

   
/* These are the functions which will output a sub-element in MATLAB ".mat"
   binary format, version 5.  This is done with a pair {tag,data}.  The data
   MUST be output in a format which matches that specified in the tag. 
 
                Input           ML5SE	: MAT version 5 sub-element (this)
                                fp      : Output file stream
                Output          void    : Output file stream is modified
                                          by having the ML5SE written to it
                Note                    : No care is taken to insure the
                                          sub-element is written and written
					  at the beginning of an element.    */
 
int MatLab5SE::write(std::fstream& fp) const
  {
  MLTag.write(fp);				// First write the tag
//  for(int i=0; i<124; i++)			// Now we output the data
//    fp.write(&(TextField[i]), sizeof(char));	// as the header text field
  return 1;
  } 

// ____________________________________________________________________________
// D               MATLAB MAT 5 Sub-Element Binary Input Functions
// ____________________________________________________________________________

/* These are the functions which will input a sub-element in MATLAB ".mat"
   binary format, version 5.  Note that in a valid MATLAB MAT file, a 
   "Sub-Element" is a combination of Tag + Data
 
                Input           ML5SE	: MAT version 5 sub-element (this)
                                fp      : Input file stream
				bigend  : Stored data endian type
				warn    : Warning output level
                                                0 = no warnings
                                                1 = warnings
                                               >1 = fatal warnings
                Output          void    : ML5SE is set from values found
                                          in the input stream
                Note                    : No care is taken to insure the
                                          sub-element is read at the start
					  of the element (as it should be
 					  in MATLAB)                         */

int MatLab5SE::read(std::fstream& fp, int bigend, int warn)
  {
  if(!MLTag.read(fp, bigend, (warn)?1:0))	// Read sub-element tag (only)
    {                                           // If this has problems, then
    if(warn)                                    // we must deal with them
      {                                         // appropriately now
      if(warn==1) MLSE5error(10, 1);
      else        MLSE5fatality(10);
      }  
    return 0;
    }
  if(MLData) delete [] MLData;			// Delete any current data
  MLData = NULL;				// Insure it is empty
  int nb = MLTag.NBytes;			// Get total number of bytes
  if(!MLTag.Compressed && nb)			// Read data only if tag is
    { 						// not compressed & has bytes
    MLData = new char[nb];			//   Allocate new data array
    for(int i=0; i<nb; i++)			//   Read in the data array
      fp.read(&(MLData[i]), sizeof(char));
    int left = nb%8;				//   End on 8 byte boundary
    if(left)
      {
      left = 8 - left;	
      char x;
      for(int k=0; k<left; k++)
        fp.read(&x, sizeof(char));
      }
//    if(left) fp.seekp(left, ios::cur);	//   so skip to this point
    }
  else 						// If it is compressed, only
    {						// the last 4 tag bytes are
    MLData = new char[4];			// actually data. Just set
    for(int i=0; i<4; i++)			// those herein
      MLData[i] = MLTag.EightBytes[4+i];	// (& hope MLBytes <= 4 )
    }
  return 1;
  }

 
// ____________________________________________________________________________
// E                MATLAB MAT 5 Sub-Element ASCII Output Functions
// ____________________________________________________________________________

                                                                                
//void MatLab5SE::print(ostream& ostr) const

        // Input                ML5SE	: MAT version 5 sub-element (this)
        //                      ostr    : An output stream
        // Output               ostr    : Output stream with MATLAB MAT file
        //                                sub-element version 5 info placed
	//				  into it
        // Note                         : This information exists at the
        //                                beginning of each sub-element

// {
// ostr << "\n\t\tSub-Element:  ";				// 8 byte field
//   for(int i=0; i<8; i++) ostr << MLTag.EightBytes[i];	// Output all bytes
// }
 
//ostream& operator<< (ostream& ostr, const MatLab5SE& ML5SE)
 
        // Input                ostr    : An output stream
        //                      ML5SE	: MAT version 5 sub-element
        // Output               ostr    : The output stream modified by
        //                                the sub-element parameters

//  { ML5SE.print(ostr); return ostr; }


#endif							// ML5SubE.cc
