/* Gascii.h *****************************************************-*-c++-*-*
**									**
**	                          G A M M A 				**
**									**
**	GAMMA ASCII File Functions			Interface	**
**									**
**	Copyright (c) 1997					 	**
**	Dr. Scott A. Smith					 	**
**	1800 E. Paul Dirac Drive					**
**	National High Magnetic Field Laboratory				**
**	Tallahassee FL 32306 USA					**
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
**  Description							 	**
**									**
**  This module deals with I/O of ASCII files in GAMMA.  Although GAMMA	**
**  data types often come with the ability to read themselves from set	**
**  parameter files (ASCII files in GAMMA parameter format), often one	**
**  needs to input a set of values (int, double, string, ....) that	**
**  control program function but have nothing to do with GAMMA objects.	**
**  These functions simplify that process by allowing flexible parsing	**
**  of ASCII files.							**
**								 	**
*************************************************************************/

#ifndef   Gamascii_h_			// Is file already included?
#  define Gamascii_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface		    	// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <string>			// Know about libstdc++ strings
#include <Matrix/row_vector.h>		// Know about row_vectors

// __________________________________________________________________________________
// i                        ASCII File I/O Error Functions
// __________________________________________________________________________________
 

void ASCIIerr(int eidx, int noret=0);

        // Input                eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message

 
void ASCIIdie(int eidx);
  
        // Input                eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message


void ASCIIerr(int eidx, const std::string& pname, int noret=0);
 
        // Input                eidx    : Flag for error type
        //                      pname   : string included in message
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message
 

// __________________________________________________________________________________
// A                        ASCII File Reading & Parsing
// __________________________________________________________________________________
 
// ----------------------------------------------------------------------------------
//                     I/O To Obtain strings From An ASCII File
// ----------------------------------------------------------------------------------


MSVCDLL std::string* Getstrings(const std::string& filename, int& N);

        // Input        filename: ASCII input file name
        //		        N	    : Number of strings
        // Return       Svect   : A string vector is filled with names
        //                        found in the file "filename"
        //			              Value of N is set to # strings found
        // Note                 : Assumes file contains 1 name per
        //                        line

 
// ----------------------------------------------------------------------------------
//                    I/O To Obtain Doubles From An ASCII File
// ----------------------------------------------------------------------------------
  

MSVCDLL double* GetDoubles(const std::string& filename, int& N);
 
        // Input        filename: ASCII input file name
        //		        N      	: Number of doubles 
        // Return       Dvect   : A double vector is filled with values
        //                        found in the file "filename"
        //			              Value of N is set to # doubles found
        // Note                 : Assumes file contains 1 value per line

 
MSVCDLL row_vector GetDoubles(const std::string& filename);
 
        // Input        filename: ASCII input file name
        // Return       data    : A vector is filled with values
        //                        found in the file "filename"
 


 
MSVCDLL void vxread(const std::string& filename, row_vector& vx);
 
	// Input           filename : ASCII input file name 
	//                 vx       : A row vector 
	// Return          none     : vx is filled with data in 
	//                            file "filename" 

// __________________________________________________________________________________
// G                       ASCII File Auxiliary Functions
// __________________________________________________________________________________


MSVCDLL void ASCIItell(const std::string& filename, int verbose=0);

	// Input           filename : ASCII input file name
	//		   vebose   : Flag how much info to output
	// Return          none     : The ASCII file is probed to determine
	//                            what type of information it contains

MSVCDLL int LineCount(const std::string& filename);
 
        // Input           filename	: ASCII input file name   
        // Return          int		: The nubmer of lines in the ASCII file
		//                            named filename

#endif                                                             // Gascii.h

