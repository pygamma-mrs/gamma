/*************************************************************************
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**	Felix		           		      Interface		**
**                                                                      **
**      Copyright (c) 1991, 1992                                        **
**      Scott Smith							**
**      Eidgenoessische Technische Hochschule                           **
**      Labor fuer physikalische Chemie                                 **
**      8092 Zurich / Switzerland                                       **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** The Felix module provides functions for input and output of GAMMA	**
** data to/from the program Felix from Molecular Simulations (MSI), 	**
** nee Biosym, nee Hare Research, Inc.					**
**                                                                      **
*************************************************************************/

#ifndef _Felix_h_			// Is file already included?
#define _Felix_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma interface			// This is the interface
#endif

#include <string>			// Include libstdc++ strings
#include <Deprecated/block_1D.h>	// Include 1D data blocks
#include <Deprecated/block_2D.h>	// Include 2D data blocks
#include <Matrix/matrix.h>		// Include matrices
//#include <iomanip.h>			// Functions setw & setprecision


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
//     (Not Really Private, But We Don't Advertise That These Exist....)
 
// ____________________________________________________________________________
// i                         Felix ERROR HANDLING
// ____________________________________________________________________________
 
 
void FelixErr(int eidx, int noret=0);

        // Input 		eidx    : Error index
        //                      noret   : Flag for return (0=linefeed)
        // Output               none    : Error message
 

volatile void FelixFatal(int eidx);
 
        // Input 		eidx    : Error index
        // Output               none	: Stop execution, output error


void FelixErr(int eidx, const std::string& pname, int noret=0);
 
        // Input 		eidx    : Error index
        //                      pname   : Additional string for error
        //                      noret   : Flag for return (0=linefeed)
        // Output               none    : Error message

 
// ____________________________________________________________________________
// ii                   Felix Serial I/O Base Functions
// ____________________________________________________________________________


void WriteFelixBlk(std::ofstream& fp,const row_vector& data,int nc,int nr,int rc);
 
        // Input                fp      : Output file pointer
        //                      data    : 1D data to be output
        //                      nc      : Number of complex values
        //                      nr      : Number of real values
        //                      rc      : flag for real versus complex
        //                                0:real, <0:imag, >0:complex

void WriteFelixBlk(std::ofstream& fp, const matrix& data, int rc=1);

        // Input                fp      : Output file pointer
        //                      data    : 2D data to be output
        //                      rc      : flag for real versus complex
        //                                0:real, <0:imag, >0:complex
        // Note                         : No error checking herein
 
void WriteFelixPSet(std::ofstream& fp, const ParameterSet &pset, int cpx, int rc=1);
   
        // Input                fp      : A file handle
        //                      pset    : Parameter set
        //                      cpx     : Number of COMPLEX points in each
        //                                line of Felix .dat file
        //                      rc      : Real (0) versus Complex(1) flag
        // Output               none    : The input parameter set is scanned
        //                                for Felix parameters.  Those found
        //                                are written out in Felix .dat header
        //                                format to the file fp.
        // Note                         : Currently only 32 Felix parameters
        //                                are recognized (from a parameter set)
        //                                as seen in the function Fp_anac.
     

// ____________________________________________________________________________
// iii                 Felix Serial I/O Base Input Functions
// ____________________________________________________________________________
  
row_vector ReadFelixBlk(std::fstream& fp, int rc, int FOR=1);

        // Input                fp      : Pointer to an open file
        //                      rc      : Flag for real vs. complex data
        //                      FOR     : Flag for reading 2 FORTRAN integers
        // Output               BLK     : 1D-data block is returned filled
        //                                from data in the file fp
        //                                presumed in Felix ".dat" format
        // Note                         : Assumes there is no header present
        // Note                         : Reads the FORTRAN integers at the
        //                                start and end of the block


// ____________________________________________________________________________
// A                       Felix .dat Output Functions
// ____________________________________________________________________________
 
/* These are the functions which will output a "serial" file in Felix fomat.
   A serial file is essentialy written as block after block after block....
   without any of the Felix matrix structure.  These are typically unprocessed
   data, although they don't have to be.  In a 2D case, the date is normally
   run through a Felix macro which reads the lines successively and processes
   the result into a spectrum for display.

           Input                filename : Output file name
                                fp       : Output file stream
                                BLK      : Output data block
                                rc       : Flag for real versus complex output
                                          0:real, <0:imag, >0:complex (default)
                                           default of 1 == complex
	  			reset	 : Flag to reset spectrum count
           Output               none     : Disk file in Felix serial format
                                           is either produced or added to
           Note                          : ".dat" format is equivalent for
                                           both Felix and its predecessor
                                           FTNMR.  Hence the file FTNMR_dat
                                           is equivalent to this.            */  

void Felix(const std::string& filename, const row_vector& data, int rc=1);
void Felix(const std::string& filename, const block_1D&   BLK,  int rc=1);
void Felix(std::ofstream &fp,           const row_vector& data, int rc=1, int rs=0);
void Felix(std::ofstream &fp,           const block_1D&   BLK,  int rc=1, int rs=0);
void Felix(const std::string& filename, const matrix&     data, int rc=1);
void Felix(const std::string& filename, const block_2D&   BLK,  int rc=1);
 
// ____________________________________________________________________________
// B                       Felix .dat Input Functions
// ____________________________________________________________________________
 
/* These are the functions which will read a "serial" file in Felix fomat.
   A serial file is essentialy stored as block after block after block....
   without any of the Felix matrix structure.  These are typically unprocessed
   data, although they don't have to be.  In a 2D case, the date is normally
   run through a Felix macro which reads the lines successively and processes
   the result into a spectrum for display.
 
           Input                filename : Input file name
                                fp       : Input file stream
				block    : Block index desired
				rowi     : Initial block index desired
				rows     : Number of blocks desired
	  Output		BLK      : 1D or 2D data block containing
					   block or blocks read in from
					   specified file or filestream
          Note				 : File assumed Felix serial!        */

block_1D Felix_1D(const std::string& filename, int block=0);
block_2D Felix_2D(const std::string& filename, int rowi=0, int rows=0);

void Felix_header(const std::string& filename);

	// Input		filename : Input file name
	// Output		none     : Read header in Felix dat file.


block_2D Felix_d_cat(const std::string& filename, int IOout = 0);

	// Input		filename : Felix concatonated dat file
	//			IOout	 : Flag for Felix reading to std. output
	// Output		BLK	 : 2D-data block is returned filled
	//				   from data in the file "filename"
	//				   presumed in Felix concatonated
	//				   ".dat" format
	// Note				 : Currently limited to 4K headers+blocks


// ____________________________________________________________________________
//                       Felix .mat FILE WRITING
// ____________________________________________________________________________


//void Felix_mat(const std::string& filename, block_2D &BLK, int rc=1);

	// Input		filename : output file name
	// 			BLK	 : 2D-data block
	//			rc       : flag for real versus complex
	//				   default of 1 == complex
	// Output		none     : disk file filename is produced
	//				 : in Felix ".mat" format


//void Felix_mat(std::ofstream &fp, block_1D &BLK, int rc=1, int reset=0);

	// Input		fp       : file handle
	//			BLK	 : 1D data block
	//			rc       : flag for real versus complex
	//			reset	 : Flag to reset spectrum count
	//				   default of 1 == complex
	// Output		none     : data in BLK is added to file
	//				   fp in Felix .mat readable format


// ____________________________________________________________________________
//                           Felix .mat FILE READING
// ____________________________________________________________________________


//block_2D Felix_mat(const std::string& filename);

	// Input		filename : input file name
	// Output		BLK	 : 2D-data block read from 
	// 				   disk file which exists
	//				   in Felix ".mat" format

//block_1D Felix_mat_1D(const std::string& filename, int block=0);

	// Input		filename : Input file name
	// 			block    : Data block desired
	// Output		BLK	 : 1D-data block is returned filled
	//				   from data in the file "filename"
	//				   presumed in Felix ".mat" format


//block_2D Felix_mat_2D(const std::string& filename, int rowi=0, int rows=0);

	// Input		filename : output file name
	// Output		BLK	 : 2D-data block read from 
	// 				   disk file filename which
	//				   exists in Felix ".mat" format


// ____________________________________________________________________________
//                          Felix AUXILIARY FUNCTIONS
// ____________________________________________________________________________


//void Felix_mat_header(const std::string& filename, int verbose=0);

	// Input		filename : Felix file name
	//			vebose   : Flag for how detailed
	// Output		cout     : output stream containing
	//				   the header information
	//				   read from the Felix file.


// ____________________________________________________________________________
//                    Felix HEADER AUXILIARY FUNCTIONS
// ____________________________________________________________________________


std::string Fp_String(int par);
  
	// Input	       par    :	Parameter number
	// Output	       string :	A string which defines the
	//				Felix .dat file parameter


std::string Fp_anac(int par);
  
	// Input	       par    :	Parameter number
	// Output	       string :	A string which is an anacronym
	//				for the Felix .dat file parameter


std::string Fp_String_qual(int pos, float parval);
  
	// Input	       par    :	Parameter number
	//		       parval : Parameter value
	// Output	       string :	A string which units and/or
	//				specifics of Felix .dat file parameter

#endif								// Felix.h
