/* ML5SubE.h ****************************************************-*-c++-*-
**									**
**	                            G A M M A			 	**
**								 	**
**	MATLAB Sub-Element (MAT Version 5)             Interface	**
**								 	**
**	Copyright (c) 1999					 	**
**	Scott Smith 	          				 	**
**      Dr. Scott A. Smith                                              **
**      1800 E. Paul Dirac Drive                                        **
**      National High Magnetic Field Laboratory                         **
**      Tallahassee FL 32306 USA                                        **
**						 			**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**							 		**
**  Description							 	**
**								 	**
** This MATLAB class provides functions for the input and output of	**
** Sub-Elements of a MATLAB MAT version 5 file.  MATLAB is a commercial	**
** product and can be purchaced from The MathWorks, Inc.  This  module	**
** deals with binary MAT files version 5 circa 1999. 			**
**							 		**
*************************************************************************/

#ifndef   _GML5SUBE_H_				// Is file already included?
#define   _GML5SUBE_H_ 				// If no, then remember it

#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// This is the interface 
#  endif

#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <GamIO/ML5Tag.h>			// Include ML5 Tags

class MatLab5SE
  {
  MatLab5Tag MLTag;				// Sub-element Tag
  char*      MLData;				// Sub-element Data

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                 MATLAB MAT 5 Data Sub-Element Error Handling
// ____________________________________________________________________________

/*              Input           MLSE5   : MAT version 5 data sub-element
                                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */  
void MLSE5error(int eidx, int noret=0);
void MLSE5error(int eidx, const std::string& pname, int noret=0);
volatile void MLSE5fatality(int eidx);

public:

friend class MatLab5AF;                         // Array Flags Full Access
friend class MatLab5DA;                         // Dimensions Array Full Access
friend class MatLab5AN;                         // Array Name Full Access
friend class MatLab5Re;				// Reals Array Full Access
friend class MatLab5Im;				// Imags Array Full Access


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------


// ____________________________________________________________________________
// A                  MATLAB MAT 5 Data Sub-Element Constructors
// ____________________________________________________________________________


MatLab5SE();
virtual ~MatLab5SE();

 
// ____________________________________________________________________________
// B                MATLAB MAT 5 Data Sub-Element Access Functions
// ____________________________________________________________________________

 
//string Class() const;
 
        // Input                ML5T    : MAT version 5 tag (this) 
        // Output               string  : String labeling the type 
 

   
// ____________________________________________________________________________
// C            MATLAB MAT 5 Data Sub-Element Binary Output Functions
// ____________________________________________________________________________

   
/* These are the functions which will output a data sub-element in MATLAB ".mat"
   binary format, version 5. 
 
                Input           ML5SE	: MAT version 5 data sub-element (this)
                                fp      : Output file stream
                Output          void    : Output file stream is modified
                                          by having the ML5SE written to it
                Note                    : No care is taken to insure the
                                          tag is written.                    */
 
virtual int write(std::fstream& fp) const;


// ____________________________________________________________________________
// D               MATLAB MAT 5 Data Sub-Element Binary Input Functions
// ____________________________________________________________________________

/* These are the functions which will input a data sub-element in MATLAB ".mat"
   binary format, version 5.  Note that in a valid MATLAB MAT file, a 
   "Data Sub-Element" is a combination of Tag + Data
 
                Input           ML5SE	: MAT version 5 tag (this)
                                fp      : Input file stream
				bigend  : Stored data endian type
                                warn    : Warning output level
                                                0 = no warnings
                                                1 = warnings
                                               >1 = fatal warnings
                Output          void    : ML5SE is set from values found
                                          in the input stream
                Note                    : No care is taken to insure the
                                          tag is read from the top of a data
                                          element (as it should be in MATLAB)*/

virtual int read(std::fstream& fp, int bigend, int warn=1);

 
// ____________________________________________________________________________
// E                MATLAB MAT 5 Data Sub-Element ASCII Output Functions
// ____________________________________________________________________________

                                                                                
//virtual void print(ostream& ostr) const;
 
        // Input                ML5SE	: MAT version 5 data sub-element (this)
        //                      ostr    : An output stream
        // Output               ostr    : Output stream with MATLAB MAT file
        //                                data sub-element version 5 info placed
	//				  into it
        // Note                         : This information exists at the
        //                                beginning of each data sub-element
 
//friend ostream& operator<< (ostream& ostr, const MatLab5SE& ML5SE);
 
        // Input                ostr    : An output stream
        //                      ML5SE	: MAT version 5 data sub-element
        // Output               ostr    : The output stream modified by
        //                                the data sub-element parameters
  };


#endif							 // _GML5SUBE_H_
