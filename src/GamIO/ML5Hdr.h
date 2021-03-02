/* ML5Hdr.h *****************************************************-*-c++-*-
**									**
**	                            G A M M A			 	**
**								 	**
**	MATLAB Header (MAT Version 5)                  Interface	**
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
** This MATLAB module provides functions for the input and output of	**
** MATLAB headers for MAT version 5 files..  MATLAB is a commercial	**
** product and can be purchaced from The MathWorks, Inc.  This  module	**
** deals with binary MAT files version 5 circa 1999. 			**
**							 		**
*************************************************************************/

#ifndef   GML5HDR_h_				// Is file already included?
#  define GML5HDR_h_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// This is the interface 
#  endif

#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <string>				// Know about libstdc++ strings
#include <fstream>				// Know about filestreams

class MatLab5Hdr

// ----------------------------------------------------------------------------
// ------------------------------- STRUCTURE ----------------------------------
// ----------------------------------------------------------------------------

/* MATLAB version 5 MAT files begin with a 128 byte header.  The first
   124 bytes are text whereas the last 4 bytes are split into two parts.
   Bytes 125-126 contain the MATLAB version and bytes 127-128 contain
   an endian indicator(MI) for byte swapping.                          */
 
  {
  char TextField[124];                          // 124 byte text field
  short Ver; 	                                // 2 short integers (4 bytes)
  char M, I;					// Endian indicator
  int size;					// Total bytes in header
  int BigEndIn;					// Big endian type read
 
  friend class MatLabFile;			// Allow MatLabFile full access
 
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                      MATLAB MAT 5 Header Error Handling
// ____________________________________________________________________________

/* 		Input 		ML5H    : MAT version 5 header	
				eidx    : Error index
               			pname   : string in message 
				noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */

void MLH5error(int eidx, int noret=0) const;
void MLH5error(int eidx, const std::string& pname, int noret=0) const; 
volatile void MLH5fatality(int eidx) const;


// ____________________________________________________________________________
// ii                   MATLAB MAT 5 Header Byte Ordering
// ____________________________________________________________________________

/* These functions deal with the header data storage endian flags.  Mostly they
   rely on the function WeRBigEnd found in the BinIOBase module which tests the
   current computer system for its "endianness".  In contrast, the header value
   BigEndIn is set when the header is read in from an external file so that we
   can know (in conjunction with WeRBigEnd) whether we need to do byte swapping
   on the values read.

           Input                ML5H    : MAT version 5 header (this)
           Output               void	: Sets the ML5H flags that indicate how
					  this computer maintains its byte 
					  order of dava values.
                                void    : If computer is big endian, MI is
                                          set to {M,I} whereas if it is small
                                          endian then MI is set to {I,M}.  It
                                          is the order of the two bytes of MI
                                          which indicate the "endianness" of
                                          stored data in MAT Ver. 5 files.   */

void SetEndian();			// Set endian flags in ML5H
void SetLittleEnd();		// Set little endian in ML5H
void SetBigEnd();			// Set big endian in ML5H
int  BigEndian();			// See if big endian input

// public:	/* Keep Everything Private Except for Friend Classes Above!! */
 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                    MATLAB MAT 5 Header Constructors
// ____________________________________________________________________________

 
MatLab5Hdr();
MatLab5Hdr(const MatLab5Hdr& MLH1);

// ____________________________________________________________________________
// B                    MATLAB MAT 5 Header Access Functions
// ____________________________________________________________________________
 
void skip(std::fstream& fp);
   
// ____________________________________________________________________________
// C                MATLAB MAT 5 Header Binary Output Functions
// ____________________________________________________________________________

/* These are the functions which will output a header in MATLAB ".mat" binary   
   format, version 5.  Note that in a valid MATLAB MAT file, this header will
   appear only once and reside at the beginning of the file.
 
                Input           ML5H    : MAT version 5 header (this)
                                fp      : Output file stream
                                warn     : Warning output level
                                                0 = no warnings
                                                1 = warnings
                                               >1 = fatal warnings
                Output          void    : Output file stream is modified
                                          by having the ML5H written to it
                Note                    : No care is taken to insure the
                                          header is written tothe top of the
                                          file (as it should be in MATLAB)   */

int write(std::fstream& fp, int warn=1) const;

// ____________________________________________________________________________
// D                MATLAB MAT 5 Header Binary Input Functions
// ____________________________________________________________________________

/* These are the functions which will input a header in MATLAB ".mat" binary    
   format, version 5.  Note that in a valid MATLAB V5 MAT file this header will
   appear only once and reside at the beginning of the file.  These functions
   can not alter the MAT file, they can only alter where the file pointer is.
   Upon a read, any data will be byte swapped as required by computer byte
   ordering differences.  Upon header output, any header text input will be
   altered accordingly.                                                      
 
                Input           ML5H    : MAT version 5 header (this)
                                fp      : Input file stream
                                warn     : Warning output level
                                                0 = no warnings
                                                1 = warnings
                                               >1 = fatal warnings
                Output          void    : ML5H is set from values found
                                          in the input stream
                Note                    : No care is taken to insure the
                                          header is read from the top of the
                                          file (as it should be in MATLAB)   */

int read(std::fstream& fp, int warn=1);

                                                                               
// ____________________________________________________________________________
// E                  MATLAB MAT 5 Header ASCII Output Functions
// ____________________________________________________________________________

/* These functions
*/
                                                                                
void print(std::ostream& ostr, int hpf=1) const;
 
        // Input                ML5H    : MAT version 5 header (this)
        //                      ostr    : An output stream
	//			hpf	: Header print flag
        // Output               ostr    : Output stream with MATLAB MAT file
        //                                header version 5 info placed into it
        // Note                         : This information exists only at the
        //                                beginning of the MAT file
 
friend std::ostream& operator<< (std::ostream& ostr, const MatLab5Hdr& MLH5);
 
        // Input                ostr    : An output stream
        //                      ML5H    : MAT version 5 header
        // Output               ostr    : The output stream modified by
        //                                the header parameters
};


#endif							 // ML5Hdr.h
