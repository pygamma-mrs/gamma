/* Felix.cc *****************************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      Felix                                         Implementation    **
**                                                                      **
**      Copyright (c) 1991, 1992                                        **
**      Scott Smith                                                     **
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
** The Felix module provides functions for input and output of GAMMA    **
** data to/from the program Felix from Molecular Simulations (MSI),     **
** nee Biosym, nee Hare Research, Inc.                                  **
**                                                                      **
*************************************************************************/

#ifndef _Felix_cc_				// Is this file included?
#define _Felix_cc_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)					// Using the GNU compiler?
#    pragma implementation				// this is the implementation
#endif

#include <GamIO/Felix.h>			// Inlcude the header file
#include <Basics/Gutils.h>			// Include GAMMA error handling
#include <Deprecated/block_1D.h>		// Include 1D data blocks
#include <Deprecated/block_2D.h>		// Include 2D data blocks
#include <Basics/StringCut.h>			// Include Gdec,Gform (dec,form)
   
union fparam					// Used for Felix parameter I/O
  {
  char fc;					// A single character
  int fi;					// A single integer
  float ff;					// A single float
  };


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

//     (Not Really Private, But We Don't Advertise That These Exist....)

// ____________________________________________________________________________
// i                         Felix ERROR HANDLING
// ____________________________________________________________________________


void FelixErr(int eidx, int noret)
     
        // Input                eidx    : Error index
        //                      noret   : Flag for return (0=linefeed)
        // Output               none    : Error message


/* The following error messages use the defaults set in the Gutils package

                Case                          Error Message                
 
                (0)                     Program Aborting.....
                (1)                     Problems With Input File Stream
                (2)                     Problems With Output File Stream
                default                 Unknown Error                        */

  {
  string hdr("Felix");
  switch (eidx)
    {
    case 3:							// (3)
      cout << "Spectrum Number Beyond Amount Present in File";
      break;
    case 6:							// (6)
      cout << "Cannot Open Felix .dat File for Writing";
      break;
    case 7:							// (7)
      cout << "Problems Reading Data From Felix .dat File....";
      break;
    case 8:							// (8)
      cout << "Trouble Reading Felix .dat File, Size = 0";
      break;
    case 11:							// (11)
      cout << "Felix Block Size Read as Zero";
      break;
    case 12:							// (12)
      cout << "End of File Reached Before Data Found";
      break;
    case 17:							// (17)
      cout << "Cannot Open Felix .mat File for Reading";
      break;
    case 18:							// (18)
      cout << "Error Reading Felix .mat File";
      break;
    case 9:							// (9)
      cout << "Cannot Handle Number of Dimensions in Felix .mat File";
      break;
    case 19:							// (19)
      cout << "Accessing Row Outside Felix Matrix Dimension";
      break;
    case 20:							// (20)
      cout << "Cannot Figure Data Point Size for Header";
      break;
    default:
      cout << "Unkown Error, number " << eidx;
    }
  cout << "\n";
  }


volatile void FelixFatal(int eidx)

	// Input		eidx  : error flag
	// Output		none  : Stops Execution & Error Message Output

  {
  FelixErr(eidx, 1);				// First output error message
  if(eidx) FelixErr(0);				// Output its fatal
  GAMMAfatal();					// Clean exit from program
  }


void FelixErr(int eidx, const string& pname, int noret)
 
        // Input                eidx    : Error index
        //                      pname   : Additional string for error 
        //                      noret   : Flag for return (0=linefeed)
        // Output               none    : Error message

/* The following error messages use the defaults set in the Gutils package

                Case                          Error Message        
 
                (1)                     Problems With File 
                default                 Unknown Error                        */

  {                                                     
  string hdr("Felix");
  switch(eidx)
    {
    case 1:  GAMMAerror(hdr, 1, pname, noret);  break; // File Problems  (1)
    default: GAMMAerror(hdr, -1, pname, noret); break; // Unknown Error  (-1)
    }
  }
     
// ____________________________________________________________________________
// ii                  Felix Serial I/O Base Output Functions
// ____________________________________________________________________________

/* These are core routines for writing to a Felix serial file.  The do NO
   error checking or anything sophisticated.  
*/

void WriteFelixBlk(ofstream& fp, const row_vector& data, int nc, int nr, int rc)
 
        // Input                fp      : Output file pointer
        //                      data    : 1D data to be output
	//			nc      : Number of complex values
	//			nr      : Number of real values
        //                      rc      : flag for real versus complex
        //                                0:real, <0:imag, >0:complex
        // Note                         : No error checking herein 
	// Note				: Row vector function write is used.
	//				  It writes out floats: Re, Im.
	// Note				: The size is written preceeding and
	//				   following the data in keeping with
	//				   the FORTRAN (used by Felix) I/O.

  {
  int m = nr*sizeof(float)			// Integer for FORTRAN blk size
                          + sizeof(int);	// 2*reals floating + 1 int
  fp.write((char*)&m, sizeof(int));			// Total storage (for FORTRAN)
  fp.write((char*)&nc, sizeof(int));			// Output number complex pts.
  if(rc)					// Output entire complex vector
    write(fp, data);				// (vector typed): r,i,r,i....
  else
    {						// Output only reals, done
    float r;					// 1 at a time: r,r,r,r,r....
    for(int i=0; i<nc; i++)
      { 
      r = float(data.getRe(i));
      fp.write((char*)&r, sizeof(float));
      } 
    } 
  fp.write((char*)&m, sizeof(int));			// Again total needed (FORTRAN)
  }


void WriteFelixBlk(ofstream& fp, const matrix& data, int rc)
 
        // Input                fp      : Output file pointer
        //                      data    : 2D data to be output
        //                      rc      : flag for real versus complex
        //                                0:real, <0:imag, >0:complex
        // Note                         : No error checking herein 

  {
  int SPECNUM = data.rows();			// No. rows in data (blocks)
  int SIZE    = data.cols();			// No. cols in data (blocksize)
  float r;
  int cmplx = SIZE;				// Set cmplx,real values count,
  if(rc==0) cmplx /= 2; 			// stored static in this module
  int reals = 2*cmplx;
  int i,j,m   = reals*sizeof(float)		// Integer for FORTRAN blk size
                          + sizeof(int);	// 2*reals floating + 1 int
  for(i=0; i<SPECNUM; i++)
    {
    fp.write((char*)&m, sizeof(int));			// Total storage (for FORTRAN)
    fp.write((char*)&cmplx, sizeof(int));		// Output number complex pts.
    if(rc)
      { 					// Output entire complex vector
      for(j=0; j<cmplx; j++) 			// pairwise: r,i,r,i,r....
        { 
        r = float(data.getRe(i,j));
        fp.write((char*)&r, sizeof(float));
        r = float(data.getIm(i,j));
        fp.write((char*)&r, sizeof(float));
        } 
      }
    else
      {						// Output only reals, done
      for(j=0; j<cmplx; j++)	 		// 1 at a time: r,r,r,r,r....
        { 
        r = float(data.getRe(i,j));
        fp.write((char*)&r, sizeof(float));
        } 
      } 
    fp.write((char*)&m, sizeof(int));			// Again total needed (FORTRAN)
    }
  }


void WriteFelixPSet(ofstream& fp, const ParameterSet& pset, int cpx, int rc)

	// Input		fp      : A file handle
	//			pset    : Parameter set
	//			cpx     : Number of COMPLEX points in each
	//				  line of Felix .dat file
	//			rc	: Real (0) versus Complex(1) flag
	// Output		none    : The input parameter set is scanned
	//				  for Felix parameters.  Those found
	//				  are written out in Felix .dat header
	//				  format to the file fp.
	// Note				: Currently only 32 Felix parameters
	//				  are recognized (from a parameter set)
	//				  as seen in the function Fp_anac.

  {
  int csize = 16;			// Currently 16 complex parameters written
  int rsize = 2*csize;			// Equivalently, 32 real parameters written
  float *header;
  header = new float[rsize];			// Allocate an array of floats for a header
  for(int i=0; i<rsize; i++)		// Initially zero all header parameters
    header[i] = 0.0;
  int m = sizeof(float)*rsize		// Get total header memory needed (for FORTRAN)
		       + sizeof(int);	// size floating + 1 int (for FORTRAN)
  fp.write((char*)&m, sizeof(int));		// Output total storage needed (for FORTRAN)
  csize *= -1;				// Output -1*csize (neg. # of complex parameters)
  fp.write((char*)&csize, sizeof(int));	// as this flags that parameters exist 

  rc = 1;				// We'll write things as complex
  cpx = abs(cpx);			// Reset Number of COMPLEX points positive

  union fparam fpar;			// This is a general Felix parameter
  double pdatad = 0;			// This will be used to read double parameters
  string param;				// Use this for Felix parameter names

//        Begin Filling in Header with p-set Felix Parameters

  for(int k=0; k<rsize; k++)		// Loop over all Felix serial parametrs
    {
    param = Fp_anac(k);			// Get parameter k's anacronym
    switch(k)				// Then write it accordingly
      {
      case 0:				// Set size of data in complex points
        if(cpx > 0) fpar.fi = cpx; 	// 1. cpx exists, set it as point size
// sosi
//        else 				// 2. try to read Felix parameter
//          pset.get(param, fpar.fi, 0);
        if(cpx <=0) FelixErr(20); 	// 3. fatal error if can't get size
        header[k] = fpar.ff;
        break;
      case 1:				// Set data type
        if(rc==0 || rc==1) fpar.fi=rc; 	// 1. rc exists, set it data type flag
// sosi
//        else pset.get(param,fpar.fi,0); // 2. try to read Felix parameter
        if(rc<0 || rc>1) fpar.fi = 1; 	// 3. assume complex points
        header[k] = fpar.ff;
        break;
      case 2: case 3: case 4:		// Read the parameters from
// sosi
//      case 5:				// the p_set as integers
//        if(pset.get(param , fpar.fi))
//          header[k] = fpar.ff;
//        break;
      case 6: case 7: case 8:		// Read these Felix parameters
      case 9: case 10: case 11:		// from the parameter set as
      case 12: case 13: case 14:	// double precision numbers
      case 15: case 16: case 17:
      case 18: case 19: case 20:
      case 21: case 22: case 23:
      case 24: case 25: case 26:
      case 27: case 28: case 29:
      case 30: case 31:
      default:
// sosi
//        if(pset.get(param , pdatad))
          header[k] = float(pdatad);
        break;
      }
    }
  fp.write((char*)header,rsize*sizeof(float));	// Output header to file
  fp.write((char*)&m, sizeof(int));
  delete [] header;		// Output total storage again (for FORTRAN)
  return;
  }

     
// ____________________________________________________________________________
// iii                 Felix Serial I/O Base Input Functions
// ____________________________________________________________________________
  
row_vector ReadFelixBlk(fstream& fp, int rc, int FOR)

	// Input		fp	: Pointer to an open file
	// 			rc      : Flag for real vs. complex data
	//			FOR	: Flag for reading 2 FORTRAN integers
	// Output		BLK	: 1D-data block is returned filled
	//				  from data in the file fp
	//				  presumed in Felix ".dat" format
	// Note				: Assumes there is no header present
	// Note				: Reads the FORTRAN integers at the
	//				  start and end of the block

  {
  float fr, fi;
  int cmplx = 0;
  int reals = 0;
  int i = 0;
  if(FOR) fp.read((char*)&i, sizeof(int));		// Storage taken (for FORTRAN)
  fp.read ((char*)&cmplx, sizeof(int));		// Read # of complex points
  if(!rc) reals = 2*cmplx; 			// Set the number of reals
  else reals = cmplx; 				// Double if reals, not complex
  row_vector data(reals);			// Declare data for output
  if(!rc)					// Do this if reading reals
    for(int j=0; j<reals; j++)			// Reals in 1 at a time
      { 
      fp.read((char*)&fr, sizeof(float));
      data.put(fr,j);
      } 
  else						// Do this if reading complex
    for(int j=0; j<cmplx; j++)			// Points read in Re, Im
      { 
      fp.read((char*)&fr, sizeof(float));
      fp.read((char*)&fi, sizeof(float));
      data.put(complex(fr, fi), j);
      } 
  if(FOR) fp.read((char*)&i, sizeof(int));		// Storage taken (for FORTRAN)
  return data;
  }


// ____________________________________________________________________________
// A                       Felix .dat Output Functions
// ____________________________________________________________________________

/* These are the functions which will output a "serial" file in Felix fomat.
   A serial file is essentialy written as block after block after block....
   without any of the Felix matrix structure.  These are typically unprocessed
   data, although they don't have to be.  In a 2D case, the date is normally
   run through a Felix macro which reads the lines successively and processes
   the result into a spectrum for display. 

           Input		filename : Output file name
           			fp       : Output file stream
           			BLK      : Output data 1D/2D block
           			data     : Output data vector/matrix
	  			rc       : Flag for real versus complex output
                                          0:real, <0:imag, >0:complex (default)
	  				   default of 1 == complex
	   Output		none     : Disk file in Felix serial format
                                           is either produced (filename input)
					   or added (file stream input)
	   Note				 : ".dat" format is equivalent for
	  				   both Felix and its predecessor
	  				   FTNMR.  Hence the file FTNMR_dat
	  				   is equivalent to this.            */

void Felix(const string& filename, const row_vector& data, int rc)
  {
  ofstream fp;					// File handle
  fp.open(filename.c_str(), ios::out);		// Open file filename
  if(!fp.is_open()) 				// Insure file is OK
    {
    FelixErr(1, filename, 1);			//   Problems with file
    FelixFatal(6);				//   Can't write to file
    }
  Felix(fp, data, rc, 1);			// Just use overload function
  fp.close();					// Close file filename
  }

void Felix(const string& filename, const block_1D& BLK, int rc)
  {
  ofstream fp;					// File handle
  fp.open(filename.c_str(), ios::out);		// Open file filename
  if(!fp.is_open()) 				// Insure file is OK
    {
    FelixErr(1, filename, 1);			//   Problems with file
    FelixFatal(6);				//   Can't write to file
    }
  Felix(fp, BLK, rc, 1);			// Just use overload function
  fp.close();					// Close file filename
  }

void Felix(ofstream &fp, const row_vector& data, int rc, int reset)
  {
  static int SPECNUM = 0;		// Track number of 1D spectra
  static int cmplx   = 0;		// Track number of cmplex values
  static int reals   = 0;		// Track number of real values
  if(reset)  SPECNUM = 0; 		// Reset spectrum counter
  if(SPECNUM == 0)			// If first data out to file then
    {					// we set static parameters for
    cmplx = data.elements();		// cmplx & real values count, then
    if(rc==0) cmplx /= 2; 		// output the header info before
    reals = 2*cmplx;			// dumping any data into file
    }
  WriteFelixBlk(fp,data,rc,cmplx,reals);// Output vector of points
  SPECNUM++;				// Increment spectrum count
  }

void Felix(ofstream &fp, const block_1D& BLK, int rc, int reset)
  {
  static int SPECNUM = 0;		// Track number of 1D spectra
  static int cmplx   = 0;		// Track number of cmplex values
  static int reals   = 0;		// Track number of real values
  if(reset)  SPECNUM = 0; 		// Reset spectrum counter
  if(SPECNUM == 0)			// If first data out to file then
    {					// we set static parameters for
    cmplx = BLK.elements();		// cmplx & real values count, then
    if(rc==0) cmplx /= 2; 		// output the header info before
    reals = 2*cmplx;			// dumping any data into file
    WriteFelixPSet(fp, BLK, cmplx, rc);
    }
  WriteFelixBlk(fp,BLK,rc,cmplx,reals);	// Output vector of points
  SPECNUM++;				// Increment spectrum count
  }


void Felix(const string& filename, const matrix& data, int rc)
  {
  ofstream fp;					// File handle
  fp.open(filename.c_str(), ios::out);		// Open file filename
  if(!fp.is_open()) 				// Insure file is OK
    {
    FelixErr(1, filename, 1);			//   Problems with file
    FelixFatal(6);				//   Can't write to file
    }
  WriteFelixBlk(fp, data, rc);			// Output data Felix .dat
  fp.close();					// Close the file
  } 


void Felix(const string& filename, const block_2D& BLK, int rc)
  {
  ofstream fp;					// File handle
  fp.open(filename.c_str(), ios::out);		// Open file filename
  if(!fp.is_open()) 				// Insure file is OK
    {
    FelixErr(1, filename, 1);			//   Problems with file
    FelixFatal(6);				//   Can't write to file
    }
  WriteFelixBlk(fp, BLK, rc);			// Output data Felix .dat
  fp.close();					// Close the file
  } 

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
          Output                BLK      : 1D or 2D data block containing
                                           block or blocks read in from
                                           specified file or filestream
          Note                           : File assumed Felix serial!        */  


block_1D Felix_1D(const string& filename, int block)

	// Input		filename : Output file name
	// 			block    : Data block desired
	// Output		BLK	 : 1D-data block is returned filled
	//				   from data in the file "filename"
	//				   presumed in Felix ".dat" format

  {
  if(block < 0) block = 0; 			// Insure a positive block
  fstream fp;					// File handle for input
  fp.open(filename.c_str(), ios::in);		// Open file filename
  if(!fp.is_open()) 				// Insure input file is OK
    {
    FelixErr(1, filename, 1);			//   Problems with file
    FelixFatal(7);				//   Problems reading .dat file
    }
  int i, ncmplx;
  fp.read((char*)&i, sizeof(int));			// Storage taken (for FORTRAN)
  fp.read((char*)&ncmplx, sizeof(int));		// Read # of complex points
  if(!ncmplx)					// Insure this is OK so far
    {
    FelixErr(11,1);				// Problems reading points
    FelixFatal(7);		 		// Problems reading .dat file
    }
  int is_header = 0;				// Flag for header existence
  int cmplx = 1;				// Flag for real vs. complex
  if(ncmplx < 0)				// If ncmplx negative, header
    {						// is present ....
    ncmplx = abs(ncmplx);			// This number of header values
    int *header;
    header = new int[2*ncmplx];			// Allocate space for header
    fp.read((char*)header, 2*ncmplx*sizeof(int));	// Read in the header values
    fp.read((char*)&i, sizeof(int));			// Storage taken (for FORTRAN)
    cmplx = header[1];				// This is # of complex points
    is_header = 1;				// Flag header 
    delete [] header;
    }
  else						// If no header, then set
    fp.seekp(0, ios::beg);			// the file back to beginning

  int j = 0;					// This is the block index
  int fpspot = fp.tellp();			// Get file location
  int blksize = 3*sizeof(int)			// Compute entire block size, 
		+ 2*ncmplx*sizeof(float);	// 3 integers & 2*size floats
  while(j < block)				// Loop block(s) we don't want
    {						// to read (by moving fp)
    fpspot += blksize;				// Where next block begins
    fp.seekp(fpspot, ios::beg);			// Skip pointer to next block
    j++;					// This is block we are at now
    }
  block_1D BLK;					// Declare BLK for output
// sosi - out until next time through
//    BLK = Felix_parameters(header);		// Set block parameters
  BLK = ReadFelixBlk(fp, cmplx);		// Read in data for block
  return BLK;
  }

  
block_2D Felix_2D(const string& filename, int rowi, int rows)

	// Input		filename : Output file name
	// 			rowi     : Index of first data block desired
	// 			rows     : Number of data blocks desired
	// Output		BLK	 : 2D-data block is returned filled
	//				   from data in the file "filename"
	//				   presumed in Felix ".dat" format

  {
  if(rowi < 0)					// Insure a positive block
    rowi = 0;
  fstream fp;					// File handle

  fp.open(filename.c_str(), ios::in);		// Open file filename	
  if(!fp.is_open())
    FelixFatal(17);

  fp.seekp(0,ios::end); 		        // Go to the end of the file
  int fsize = fp.tellp();			// Get total bytes in the file
  fp.seekp(0,ios::beg);				// Back to beginning of file

  int i, ncmplx;
  fp.read((char*)&i, sizeof(int));			// Storage taken (for FORTRAN)
  fp.read((char*)&ncmplx, sizeof(int));		// Read # of complex points
  if(ncmplx == 0)				// Insure points present
    FelixFatal(11);
  int is_header = 0;				// Flag for header existence
  int cmplx = 1;				// Flag for real vs. complex
  if(ncmplx < 0)				// Assumed complex if no
    {						// parameters in Felix file
    ncmplx = abs(ncmplx);
    int *header;
    header = new int[2*ncmplx];
    fp.read((char*)header, 2*ncmplx*sizeof(int));
    fp.read((char*)&i, sizeof(int));			// Storage taken (for FORTRAN)
    cmplx = header[1];
    is_header = 1;
    delete [] header;
    }
  else						// If no header, then set
    fp.seekp(0, ios::beg);			// the file back to beginning

// sosi - modified by tilo here in alpha version

  int fpspot = fp.tellp();			// Get file location
  int blksize = 3*sizeof(int)			// Compute entire block size, 
		+ 2*ncmplx*sizeof(float);	// 3 integers & 2*size floats
  int trows = (fsize-fpspot)/blksize;		// Number of rows in the file
  if(rowi > trows) FelixFatal(12);
  if(rowi+rows > trows) FelixFatal(12);
  if(rows==0) rows=trows-rowi;

  fpspot += rowi*blksize;
  fp.seekp(fpspot, ios::beg);			// Skip a block

  block_2D BLK(rows, ncmplx);			// Declare 2D-BLK for output
// Need a check for reaching the end of file
//    BLK = Felix_parameters(header);		// Attach any parameters
  int j = 0;
  while(j < rows)
    {
    BLK.put_block(j,0,ReadFelixBlk(fp, cmplx));
    j++;
    }
  fp.close();					// Close file filename
  return BLK;
  }


void Felix_header(const string& fname)

	// Input		fname    : Input file name
	// Output		none     : Read header in Felix dat file.

  {
  union fparam fpar;			// A general Felix parameter
  int size = 0;				// Integer for size
  fstream fp;				// Declare a pointer to a file
  fp.open(fname.c_str(),ios::in);	// Open the file for reading
  fp.read((char*)&size, sizeof(int));		// Read intial integer (for FORTRAN)
//  cout << "\n\t     " << itoa(size,10,8)
  cout << "\n\t     " << Gdec(size,8)
       << "    \t" << Fp_String(-2);
  fp.read((char*)&size, sizeof(int));		// Read size of block
  if(size==0)				// Error if size is zero
    FelixErr(8);
  else if(size >0)
    cout << "\nThere Were No Parameters Found in this Felix .dat File."
	 << "\nThe Number of Complex Points = " << size << "\n";
  else
    {
//    cout << "\n\t     " << itoa(size,10,8)
    cout << "\n\t     " << Gdec(size,8)
         << "    \t" << Fp_String(-1);
    size = 2*abs(size);			// Set number of reals
    float *header;
    header = new float[size];			// Array for header
    fp.read((char*)header, size*sizeof(int));	// Read Felix header
    for(int k=0; k<size; k++)
      {
      fpar.ff = header[k];
//      cout << "\n\t" << itoa(k+1,10,3) << ". ";
      cout << "\n\t" << Gdec(k+1,3) << ". ";
      switch(k)
        {
        case 0: case 1: case 2:		// Write as an integer
        case 3: case 4: case 5:
//          cout << itoa(fpar.fi,10,8)
          cout << Gdec(fpar.fi,8)
               << "    ";
          break;
        case 6: case 7: case 8:		// Write as a float
        case 9: case 10: case 11:
        case 12: case 13: case 14:
        case 15: case 16: case 17:
        case 18: case 19: case 20:
        case 21: case 22: case 23:
        case 24: case 25: case 26:
        case 27: case 28: case 29:
        case 30: case 31:
        default:
          cout << Gform("%10.2f", fpar.ff);
          break;
        }
      cout << "\t" << Fp_String(k) << Fp_String_qual(k, header[k]);
      delete [] header;
      }
    fp.read((char*)&size, sizeof(int));	// Read last integer (for FORTRAN)
//    cout << "\n\t     " << itoa(size,10,8)
    cout << "\n\t     " << Gdec(size,8)
         << "    \t" << Fp_String(-2);
    cout << "\n";
    }
  return;
  }


 block_2D Felix_d_cat(const string& filename, int IOout)

	// Input		filename : Felix concatonated .dat file
	//			IOout	 : Flag for Felix reading to std. output
	// Output		BLK	 : 2D-data block is returned filled
	//				   from data in the file "filename"
	//				   presumed in Felix concatonated
	//				   ".dat" format
	// Note				 : Currently limited to 4K headers+blocks

  {
  const int limit=4096;				// Change this if > 4K blocks
  fstream fp;					// File handle
  fp.open(filename.c_str(), ios::in);		// Open Felix file filename
  if(!fp.is_open()) FelixFatal(7); 		// Check that file is O.K.

		  // First Scan File for Information

  fp.clear();
  int IObad = 0;				// Flag for good I/O
  int endoffile = 0;				// Flag for reaching end of Felix file
  int blocksize = 0;
  int i=0, ncmplx=0, icnt=0, fpspot=0;
  int storage[limit];				// Declare array for block sizes
  int sizes[limit];				// Declare array for data sizes
  int cmplx = 1;				// Flag for real vs. complex data
  while(!IObad && !endoffile)
    {
    fp.read((char*)&i, sizeof(int));			// Storage taken (for FORTRAN)
    IObad = fp.fail();				// Test file status
    endoffile = fp.eof();			// Test end of file status
    if(!IObad && !endoffile)
      {
      fpspot = fp.tellp();			// Store beginning of block point
      fp.read((char*)&ncmplx, sizeof(int));		// Read # of complex points stored
      storage[icnt] = ncmplx;			// Storage allocated to data this block
      if(ncmplx < 0)				// If a header, get size & type of points
        {					// actually stored in following blocks
        fp.read((char*)&i, sizeof(int));		// 1st Header Parameter, # points
        fp.read((char*)&cmplx, sizeof(int));		// 2st Header Parameter, real vs. complex
        if(!cmplx)
          sizes[icnt] = -i;			// Store real data sizes negative 
        else
          sizes[icnt] = i;			// Store complex data sizes positive
        }
      else if(icnt)				// If data, use size reported last header
        sizes[icnt] = sizes[icnt-1];
      else					// If data at start, take storage size
        sizes[icnt] = ncmplx;
      blocksize = 2*abs(ncmplx)*sizeof(float)	// Block is complex points and 2 integers
		               + 2*sizeof(int);	// at start and end put in by FORTRAN I/O
      fpspot += blocksize;			// Compute end of block file point
      fp.seekp(fpspot, ios::beg);		// Skip file to next Felix block
      icnt++;
      IObad = fp.fail();			// Test file status
      endoffile = fp.eof();			// Test end of file status
      }
    else if(endoffile)
      if(IOout)
        cout << "\n\tEnd of Felix File Reached, " << icnt << " Blocks Read (1)";
    else
      if(IOout)
        cout << "\n\tProblems Reading Concatonated Felix File(1)?";
    }

		  // Now Figure Out What Is In The Felix File

  int size = -32768;				// Sets max block to 32K
  int bsize = 0;
  int hsize = 0;
  int nheaders = 0;
  int ndatablks = 0;
  int k;
  for(k=0; k<icnt; k++)			// Search for the biggest stored block
    {						// and the largest data set
    i = storage[k];
    if(i < 0)
      {
      nheaders++;
      if(hsize > i) hsize = i;
      }
    if(i > 0)
      {
      ndatablks++;
      if(bsize < i) bsize = i;
      i = sizes[k];
      if(size < i) size = i;
      }
    }
  if(IOout)
    {
    cout << "\n\t" << nheaders << " Felix Headers Read, Largest Uses " << abs(hsize) << " Complex Points";
    cout << "\n\t" << ndatablks  << " Felix Data Blocks Read, Largest Uses " << bsize << " Complex Points";
    cout << "\n\tLargest Felix Data Block is " << abs(size);
    if(size > 0) cout << " Complex";
    else if(size < 0) cout << " Real";
    cout << " Points";
    cout << "\n";
    }
		  // Now Re-Read Felix File, Store Data into Array

  fp.close();					// Close the Felix file
  fp.open(filename.c_str(), ios::in);		// Re-open Felix file filename
  if(!fp.is_open()) FelixFatal(7); 		// Check that file is O.K.
  size = abs(size);
  block_2D BLK(ndatablks, size);		// Construct matrix for data blocks
  float fr, fi;					// For individual data points
  IObad = 0;					// Reset flag for good I/O
  endoffile = 0;				// Reset flag for reaching end of Felix file
  int nreals = 0;				// Number of real points
  int bcnt = 0;
  if(IOout)
    cout << "\n\tReading Felix Concatonated File Into " << ndatablks << " by " << size << " Array";
  for(k=0; k<icnt && !IObad && !endoffile; k++)
    {
    fp.read((char*)&i, sizeof(int));			// Storage taken (for FORTRAN)
    fp.read((char*)&ncmplx, sizeof(int));		// Read # of complex points stored
    if(ncmplx != storage[k])			// Check versus previous read
      cout << "\nTrouble with block storage size";
    IObad = fp.fail();				// Test file status
    endoffile = fp.eof();			// Test end of file status
    if(ncmplx<0 && !IObad && !endoffile)	// This is a header, get real vs. complex
      {
      fp.read((char*)&i, sizeof(int));			// 1st Header Parameter, # complex points
      fp.read((char*)&cmplx, sizeof(int));		// 2st Header Parameter, real vs. complex
      if(!cmplx)				// Real data, set # of real points
        nreals = i;
      fpspot = fp.tellp();			// Get file location
      blocksize = 2*abs(ncmplx)*sizeof(float)	// Block is all data + last integer that
				 + sizeof(int); // was put in by FORTRAN I/O
      blocksize -= 2*sizeof(float);		// Except the two integers just read
      fpspot += blocksize;			// Figure where the next block starts then
      fp.seekp(fpspot, ios::beg);		// Skip rest of Felix header
      }
    else if(ncmplx>0 && !IObad && !endoffile)	// Data, use real vs. complex last header
      {
      fpspot = fp.tellp();				// Get initial location
      blocksize = 2*abs(ncmplx)*sizeof(float);		// Actual storage block uses
      if(IOout)
        cout << "\n\tData block " << bcnt+1;
      int j;
      if(!cmplx)					// Real data
        {
        for(j=0; j<nreals && !fp.fail(); j++)	// Points in 1 at a time
          { 
          fp.read((char*)&fr, sizeof(float));
          BLK.put(fr, bcnt, j);
          IObad = fp.fail();
          } 
        if(IOout)
          cout << ": " << j << " Real Points";
        }
      else						// Complex data
        {
        for(j=0; j<ncmplx && !fp.fail(); j++)	// Points read in Re, Im
          { 
          fp.read((char*)&fr, sizeof(float));
          fp.read((char*)&fi, sizeof(float));
          BLK.put(complex(fr, fi), bcnt, j);
          IObad = fp.fail();
          } 
        if(IOout)
          cout << ": " << j << " Complex Points";
        }
      fp.seekp(fpspot+blocksize, ios::beg);		// Skip to end of block storage
      fp.read((char*)&i, sizeof(int));				// Storage taken (for FORTRAN)
      bcnt++;
      }
    else if(endoffile)
      if(IOout)
        cout << "\n\tEnd of Felix File Reached, " << icnt << " Blocks Read(2)";
    else
      if(IOout)
        cout << "\n\tProblems Reading Concatonated Felix File(2)?";
    }
  fp.close();						// Close the Felix file
  cout << "\n";
  return BLK;
  }


// ______________________________________________________________________
//                        Felix .mat FILE WRITING
// ______________________________________________________________________

/*

void Felix_mat_pset(fstream& fp, const ParameterSet& pset,
                                            int ndim, int *sizes, int rc)

	// Input		fp      : A file handle
	// 			pset    : Parameter set
	//			sizes	: Sizes in points of each dimension
	//			rc	: Real (0) versus Complex(1) flag
	// Output		none    : pset is scanned for Felix parameters
	//				  those found are put into header
	// Note				: size

  {
  int blks[ndim];			// Array of blocks sizes
  int i;
  for(i=0; i<ndim; i++)			// This is just 64 point chunks
    blks[i] = sizes[i]/64;		// with twice a many in dim2 if
  if((ndim > 1) && (rc)) blks[1] *= 2; 	// data is complex

  int dbrick=1;				// Determine the number of
  for(i=0; i<ndim; i++)			// bricks in the data
    dbrick *= blks[i];

  int header[4096];			// Set up Felix .mat header
  for(i=0; i<4096; i++)			// Zero all header parameters
    header[i]=0;
					// Begin filling in parameters
  header[0] = ndim;			// Number of dimensions
  if(rc)
    header[1] = 2;			// Data Type, Felix has 2 = complex
  else					//                      1 = real
    header[1] = 1;
  header[3] = dbrick+1;			// Total bricks present
//
//			Writing Vector Values
//	    Begin at 20, ndim Written for Each Parameter
//
  int indx = 20;
  for(i=0; i<ndim; i++)
    {					// Unknown Parameter
    header[indx] = 1;			// Assumed 1 for each dimension
    indx++;				// 2D indices: 20, 21
    }
  for(i=0; i<ndim; i++)
    {					// Points in each dimension
    header[indx] = sizes[i];		// Irreguardless of real vs cmplx,
    indx++;				// these are the matrix dimensions
    }					// 2D indices: 22,23
  for(i=0; i<ndim; i++)
    {					// Blocks in each dimension
    header[indx] = blks[i];		// Just points/64. 2nd dim.
    indx++;				// is doubled if complex
    }					// 2D indices: 24,25
  for(i=0; i<ndim; i++)
    {					// Unknow Parameter
    header[indx] = blks[0];		// Seems 1st is 1 and rest are
    if(i==0)				// blocks in dim1
      header[indx] = 1;			// 2D indices: 26,27
    indx++;
    }
  for(i=0; i<ndim; i++)
    {					// Unknown Parameter
    header[indx] = 64;			// Seems to be points/block
    if(i==1 && rc)			// with 1/2 2nd dim. if complex
      header[indx] = 32;		// 2D indices; 28,29
    indx++;
    }
  for(i=0; i<ndim; i++)
    {					// Unknow Parameter
    header[indx] = 64;			// Seems 1st is 1 and rest are
    if(i==0)				// 64
      header[indx] = 1;			// 2D indices: 30,31
    indx++;
    }
  for(i=0; i<ndim; i++)
    {					// Unknown Parameter
    header[indx] = 1140457472;		// Seems to be this number always
    indx++;				// 2D indices: 32,33
    }
  for(i=0; i<ndim; i++)
    {					// Referencing (Unknown)
    header[indx] = 1;			// Seems to be always 1
    indx++;				// 2D indices: 40,41
    }
  header[100] = 1;

  header[101] = dbrick+1;		// Segments of bricks

  header[111] = 5;			// Length of full filename 
  header[121] = int('G');		// Actual filename
  header[122] = int('A');		// GAMMA put in here
  header[123] = int('M');		// Felix, when constructing
  header[124] = int('M');		// a matrix puts the full
  header[125] = int('A');		// filename with path here

  indx = 220;
  char a, b;
  a = 'D';
  for(i=0; i<ndim; i++)			// Fill in ASCII numbers for
    {					// D1, D2, ...., D8 referencing
    b = char(i);
    header[indx] = int(a);
    indx++;
    header[indx] = int(b);
    indx++;
    }

  header[4096] = 1;			// Tell Felix its a nice mat file

					// Output the header
  fp.write(header, 4096*sizeof(int));	// Write the header
  return;
// sosi
//  if(pset.status() == 0) indx = 0;	// Compiler likes this used
  }


// Tilo - Again doesn't work because get_block does not return block_1D

void Felix_mat(const string& filename, block_2D& BLK, int rc)

	// Input		filename : output file name
	//			BLK	 : 2D data block
	//			rc       : flag for real versus complex
	//				   1:complex, 0:real
	// Output		none     : disk file filename is produced
	//				 : in Felix .dat readable format

  {
  ofstream fp;				// Declare a pointer to a file
  fp.open(filename.c_str(),ios::out);	// Open the file for writing
  if(!fp.is_open()) FelixFatal(17);
  int SPECNUM = BLK.rows();		// Number of rows in data block
  int SIZE = BLK.cols();		// Number of cols in data block
  block_1D BLK1D(SIZE);
  BLK1D = BLK.get_block(0,0,1,SIZE);
  Felix_mat(fp, BLK1D, rc, 1);
  for(int i=1; i<SPECNUM; i++)		// Write successive data blocks
    {					// in Felix .dat format row by row
    BLK1D = BLK.get_block(i,0,1,SIZE);
    Felix_mat(fp, BLK1D, rc);
    }
  fp.close();				// Close the file
  } 
*/


/*
void Felix_mat(ofstream& fp, block_1D& BLK, int rc, int reset)

	// Input		fp       : file handle
	//			BLK	 : 1D data block
	//			rc       : flag for real versus complex
	//				   1:complex, 0:real
	//			reset	 : Flag for resetting spectrum count
	// Output		none     : data in BLK is added to file
	//				   fp in Felix .dat readable format
	// Note				 : The function write(file,BLK) is
	//				   supplied in class matrix but used
	//				   with derived class block_1D.  It
	//				   writes out floats: Re, Im.
	// Note				 : The size is written preceeding and
	//				   following the data in keeping with
	//				   the FORTRAN (used by Felix) I/O.

  {
  static int SPECNUM = 0;		// Track number of 1D spectra
  static int cmplx = 0;			// Track number of cmplx
  static int reals = 0;			// Track number of reals
  if(reset) SPECNUM = 0; 		// Reset spectrum counter
  if(SPECNUM == 0)			// First data out to file !
    {					// Set static parameters
    cmplx = BLK.elements();		// cmplx and reals then
    if(rc==0)				// output the header info
      cmplx /= 2;
    reals = 2*cmplx;
//    int ndim = 2;
    int sizes[2];
    int nrows = cmplx;
//    int ncols = 128;
    sizes[0] = nrows;
    sizes[1] = BLK.elements();
// sosi
//    Felix_mat_pset(fp, BLK, ndim, sizes, rc);
    }
//else					// Multiple data out to file !
//    {					// Insure the sizes are
//					// compatible
//    }

  int m = reals*sizeof(float)		// Integer for FORTRAN blk size
                          + sizeof(int);// 2*reals floating + 1 int
  fp.write(&m, sizeof(int));		// Total storage (for FORTRAN)
  fp.write(&cmplx, sizeof(int));	// Output the number complex pts.
  if(rc)				// Output entire complex BLK
    write(fp, (row_vector)BLK);		// (vector typed): r,i,r,i....
  else
    {					// Output only reals, done
    float r;				// 1 at a time: r,r,r,r,r....
    for(int i=0; i<cmplx; i++)
      { 
      r = float(Re(BLK(i)));
      fp.write(&r, sizeof(float));
      } 
    } 
  fp.write(&m, sizeof(int));		// Again total needed (FORTRAN)
  SPECNUM++;				// Increment spectrum count
  }

*/

// ______________________________________________________________________
//                        Felix .mat FILE READING
// ______________________________________________________________________


/*
block_1D Felix_mat_BLK(fstream &fp, int pts, int cmplx, int pset)

	// Input		fp	 : Pointer to an open Felix .mat file
	// 			pts      : Number of points to return 
	// 			cmplx    : Flag for real vs. complex data
	//			pset	 : Flag if parameters need including
	// Output		BLK	 : 1D-data block is returned filled
	//				   from data in the file fp
	//				   presumed in Felix ".mat" format
	// Note				 : The value of pts must be a
	//				   positive integer multiple of 64
	// Note				 : Currently only reads blocks
	//				   parallel to the Felix dim1

  {
  block_1D BLK(pts);			// Declare BLK for output
  int j = 0;
  int k = 0;
  float r = 0;
  int blks = pts/64;			// Number of Felix 64 pt blocks 
int fpspot = fp.tellp();		// Get file location
cout << "\n\tData Requested Takes " << blks << " Felix 64 point blocks";
cout << "\n\tTry a multiple of 64 Points for skipping: ";
int imult64;
cin >> imult64;
cout << "\n";
  while(blks)
    {
    for(j=0; j<64; j++)			// Reals in 64 points at a time
      { 				// This is a block inside of Felix
      fp.read(&r, sizeof(float));
      BLK(k+j) = r;
      } 
cout << "\n\t64 Points read into real values";
fp.seekp(fpspot+(imult64*64*sizeof(float)));			// Get file location
fpspot = fp.tellp();
cout << "\n\tSkipping " << imult64*64 << " points";
// NEED TO JUMP THE FILE POINTER HERE YOU IDIOT!
    if(cmplx)
      {  
      for(j=0; j<64; j++)		// Imaginaries in 64 points a time
        { 
        fp.read(&r, sizeof(float));
        BLK(k+j) += complex(0, r);
        } 
cout << "\n\t64 Points read into imaginary values";
fp.seekp(fpspot+(imult64*64*sizeof(float)));			// Get file location
fpspot = fp.tellp();
cout << "\n\tSkipping " << imult64*64 << " points";
      } 
    k += 64;				// Increment point index
    blks--;				// Decrement block index
    } 
  return BLK;
  k = pset;				// Compiler likes this used
  }


block_1D Felix_mat_1D(const string& filename, int row)

	// Input		filename : Output file name
	// 			row      : Felix matrix row desired [0, n)
	// Output		BLK	 : 1D-data block is returned filled
	//				   from data in the file "filename"
	//				   presumed in Felix ".mat" format

  {
  const int Fmatblock = 64*sizeof(float);	// Size of a Felix .mat block (64 pts)
  const int Fmatbrick = 4096*sizeof(float);	// Size of a Felix .mat brick (4K pts)
  fstream fp;					// File handle
  fp.open(filename.c_str(), ios::in);		// Open Felix .mat file filename
  if(!fp.is_open()) FelixFatal(17); 		// Insure the matrix file is open
  int header[4096];				// Allocate Felix .mat header array
  fp.read(header, 4096*sizeof(int));		// Read it in from the .mat file
  int nd = header[0];				// Number of matrix dimensions
  if(nd > 2)
    {
    cout << "\n\tThe routine Felix_mat_1D is valid only for 2D matrices!";
    cout << "\n\tSorry, can't handle a " << nd << "D matrix";
    cout << "\n\tAborting.......";
    FelixFatal(0);
    }
  int nrows = header[20+nd+1];			// Number of matrix rows
  int ncols = header[20+nd];			// Number of matrix columns
  if((row < 0) || (row >= nrows))		// Insure row actually exists
    FelixFatal(19);
  int cmplx = (header[1]) - 1;			// Set real vs. complex matrix flag
  int rowsize = ncols*sizeof(float);		// Compute entire row size, 
  if(cmplx)					// Twice a large if complex data
    rowsize *= 2;
  int rpbr = header[20+4*nd+1];		 	// Rows stored in each Felix brick
  int brkinc = ncols/64;			// Bricks to complete rpbr rows
  int ibrks = brkinc*int(row/rpbr);		// Bricks preceding this rows data
  int iblks = row%rpbr;				// Blocks preceding this rows data
  int fpspot = fp.tellp();			// Get file location (after header)
  fpspot += iblks*Fmatblock + ibrks*Fmatbrick;	// Move to data start for this row
  fp.seekp(fpspot, ios::beg);			// Set file to point at row start					
  block_1D BLK(ncols);				// Declare BLK for output
  int j = 0;
  int k = 0, kr = 0, ki = 0;
  float r = 0;
  while(k < ncols)				// Begin reading in the row
    {
    fpspot = fp.tellp();			// Get this blocks file location
    for(j=0; j<64 && kr<ncols; j++)		// Read a Felix block of reals
      {
      fp.read(&r, sizeof(float));
      BLK(k+j) = r;
      kr++;
      } 
    fpspot += rpbr*Fmatblock;
    fp.seekp(fpspot, ios::beg);			// Skip to next Felix block this row
    if(cmplx)
      {  
      fpspot = fp.tellp();			// Get this blocks file location
      for(j=0; j<64 && ki<ncols; j++)		// Read a Felix block of imaginaries
        { 
        fp.read(&r, sizeof(float));
        BLK(k+j) += complex(0, r);
        ki++;
        } 
      fpspot += rpbr*Fmatblock;
      fp.seekp(fpspot, ios::beg);		// Skip to next Felix block this row
      } 
    k += 64;					// Increment row point index
    } 
  return BLK;
  }


block_1D Felix_mat_1Dx(const string& filename, int row)

	// Input		filename : Output file name
	// 			row      : Felix matrix row desired [0, n)
	// Output		BLK	 : 1D-data block is returned filled
	//				   from data in the file "filename"
	//				   presumed in Felix ".mat" format

  {
  const int Fmatblock = 64*sizeof(float);	// Size of a Felix .mat block (64 pts)
  const int Fmatbrick = 4096*sizeof(float);	// Size of a Felix .mat brick (4K pts)
  fstream fp;					// File handle
  fp.open(filename.c_str(), ios::in);		// Open Felix .mat file filename
  if(!fp.is_open())				// Insure the matrix file is open
    FelixFatal(17);
  int header[4096];				// Allocate Felix .mat header array
  fp.read(header, 4096*sizeof(int));		// Read it in from the .mat file
  int nd = header[0];				// Number of matrix dimensions
  if(nd > 2)
    {
    cout << "\n\tThe routine Felix_mat_1D is valid only for 2D matrices!";
    cout << "\n\tSorry, can't handle a " << nd << "D matrix";
    cout << "\n\tAborting.......";
    FelixFatal(0);
    }
//  int nrows = header[20+nd];			// Number of matrix rows
//  int ncols = header[20+nd+1];			// Number of matrix columns
  int nrows = header[20+nd+1];			// Number of matrix rows
  int ncols = header[20+nd];			// Number of matrix columns
  if((row < 0) || (row >= nrows))		// Insure row actually exists
    FelixFatal(19);
  int cmplx = (header[1]) - 1;			// Set real vs. complex matrix flag
  int rowsize = ncols*sizeof(float);		// Compute entire row size, 
  if(cmplx)					// Twice a large if complex data
    rowsize *= 2;
  int rpbr = header[20+4*nd+1];		 	// Rows stored in each Felix brick
  int brkinc = ncols/64;			// Bricks to complete rpbr rows
  int ibrks = brkinc*int(row/rpbr);		// Bricks preceding this rows data
  int iblks = row%rpbr;				// Blocks preceding this rows data
cout << "\n\tNumber of Matrix Rows = " << nrows;
cout << "\n\tNumber of Matrix Columns = " << ncols;
cout << "\n\tRows per brick = " << rpbr;
cout << "\n\tBricks to Complete" << rpbr << " Rows = " << brkinc;
cout << "\n\tBricks preceding this row = " << ibrks;
cout << "\n\tBlocks preceding this row = " << iblks;
cout << "\n";
  int fpspot = fp.tellp();			// Get file location (after header)
cout << "\n\tFile location after header " << fpspot;
  fpspot += iblks*Fmatblock + ibrks*Fmatbrick;	// Move to data start for this row
cout << "\n\tSkipping " << ibrks << " Bricks, Hopefully " << rpbr*ibrks
     << " Complete Rows (+" << ibrks*Fmatbrick << ")";
cout << "\n\tSkipping " << iblks << " Blocks, Hopefully "
     << "to the Start of Felix Row " << row+1 << " (+" << iblks*Fmatblock << ")";
cout << "\n\tTotal Blocks Skipped = " << 64*ibrks+iblks;
cout << "\n";
cout << "\n\tFile location after skip " << fpspot << "\n";
//cout << "\n\n\tRows per brick = " << rpbr << ", Input a new value: ";
//cin >> rpbr;
  fp.seekp(fpspot, ios::beg);			// Set file to point at row start					
  block_1D BLK(ncols);				// Declare BLK for output
  int j = 0;
  int k = 0, kr = 0, ki = 0;
  float r = 0;
  while(k < ncols)				// Begin reading in the row
    {
    fpspot = fp.tellp();			// Get this blocks file location
    for(j=0; j<64 && kr<ncols; j++)		// Read a Felix block of reals
      {
      fp.read(&r, sizeof(float));
      BLK(k+j) = r;
      kr++;
      } 
cout << "\n\tSkipping " << rpbr << " Felix 64 pt blocks to next read (+" << rpbr*Fmatblock << ")"; 
    fpspot += rpbr*Fmatblock;
    fp.seekp(fpspot, ios::beg);			// Skip to next Felix block this row
cout << "\n\tFile location after skip: " << fpspot << "\n";
    if(cmplx)
      {  
      fpspot = fp.tellp();			// Get this blocks file location
      for(j=0; j<64 && ki<ncols; j++)		// Read a Felix block of imaginaries
        { 
        fp.read(&r, sizeof(float));
        BLK(k+j) += complex(0, r);
        ki++;
        } 
cout << "\n\tSkipping " << rpbr << " Felix 64 pt blocks to next read"; 
      fpspot += rpbr*Fmatblock;
      fp.seekp(fpspot, ios::beg);		// Skip to next Felix block this row
cout << "\n\tFile location after skip: " << fpspot << "\n";
      } 
    k += 64;					// Increment row point index
    } 
  return BLK;
  }


block_2D Felix_mat_2D(const string& filename, int rowi, int rows)

	// Input		filename : output file name
	// Output		BLK	 : 2D-data block read from 
	// 				   disk file filename which
	//				   exists in Felix ".mat" format

{
fstream fp;				// Declare a pointer to a file
fp.open(filename.c_str(), ios::in);	// Open the file for reading
if(!fp.is_open()) FelixFatal(17);
int header[4096];			// Read Felix .mat header
fp.read(header,4096*sizeof(int));
if(header[0] >2 || header[0]<1)		// Check number of dimensions
  FelixFatal(9);
int ncols = header[22];			// Retrieve number of points
int nrows = 1;
if(header[0] == 2)
  nrows = header[23];
block_2D BLK(nrows, ncols);		// Set up 2D block for output
return BLK;
rowi = rows;				// Compiler likes this used
}

*/

// ______________________________________________________________________
//                   Felix FILE AUXILIARY ROUTINES
// ______________________________________________________________________

/*

void Felix_mat_header(const string& filename, int verbose)

	// Input		filename : Felix file name
	//			vebose   : Flag for how detailed
	// Output		cout     : output stream containing
	//				   the header information
	//				   read from the Felix file.

{
  int header[4096];				// Array of 4096 integers
  fstream fp;					// Declare a pointer to a file
  fp.open(filename.c_str(), ios::in);		// Open the file for reading
  if(!fp.is_open()) FelixFatal(17); 		// Insure file opened O.K.
  fp.read(header, 4096*sizeof(int));		// Read Felix .mat header
  int ndim = header[0];				// Store # of matrix dimensions
  if(ndim <0 || ndim>8)				// Check for reasonable value
    {
    FelixErr(17);
    cout << "\nThe Number of Dimensions Read is " << ndim;
    cout << "\nResetting to 2, Hopefully for Some Logical Output";
    }
  int prnt = 0;					// Used for a print flag

// 0 ************** Begin Outputting Scalar Parameters *************** 19 //

    cout << "\n\t\tFelix Header Scalar Parameters\n";
    cout << "\n     1. " << header[0] << "\tNumber of Matrix Dimensions";
    cout << "\n     2. " << header[1] << "\tMatrix Data Type Flag";
    cout << "\n     3. " << header[2] << "\tTotal Number of Matrix Bricks (4K Point Chunks)";
    cout << "\n     4. " << header[3] << "\tNumber of Blocks per Brick";
    int kk;
    for(kk=4; kk<9 && verbose; kk++)
      cout << "\n     " << kk+1 << ". " << header[kk] << "\tUnknown";
    for(kk=9; kk<20 && verbose>1; kk++)
      cout << "\n    " << kk+1 << ". " << header[kk] << "\tUnknown";

// 20 ************* Begin Outputting Vector Parameters ************** 100 //

    cout << "\n\n\t\tFelix Matrix Header Vector Parameters\n";
    int vindex=20;
    int i;
    for(i=0; i<ndim; i++)
      {
      cout << "\n    " << vindex+1 << ". " << header[vindex];
      cout << "\tMatrix Data Type Dimension " << i+1;
      vindex++;
      }
    cout << "\n";
    for(i=0; i<ndim; i++)
      {
      cout << "\n    " << vindex+1 << ". " << header[vindex];
      cout << "\tPoint Size Dimension " << i+1;
      vindex++;
      }
    cout << "\n";
    for(i=0; i<ndim; i++)
      {
      cout << "\n    " << vindex+1 << ". " << header[vindex];
      cout << "\tNumber of Blocks Dimension " << i+1;
      vindex++;
      }
    cout << "\n";
    for(i=0; i<ndim; i++)
      {
      cout << "\n    " << vindex+1 << ". " << header[vindex];
      cout << "\tUnknown Dimension " << i+1;
      vindex++;
      }
    cout << "\n";
    for(i=0; i<ndim; i++)
      {
      cout << "\n    " << vindex+1 << ". " << header[vindex];
      cout << "\tBlock Size Dimension " << i+1;
      vindex++;
      }
    cout << "\n";
    for(i=0; i<ndim; i++)
      {
      cout << "\n    " << vindex+1 << ". " << header[vindex];
      cout << "\tUnknown Dimension " << i+1;
      vindex++;
      }
    cout << "\n";
    for(i=0; i<ndim; i++)
      {
      cout << "\n    " << vindex+1 << ". " << header[vindex];
      cout << "\tUnknown Dimension " << i+1;
      vindex++;
      }
    cout << "\n";

// 101 ************* Begin Outputting Other Parameters ************* 4094 //

    if(verbose)
      {
      cout << "\n\n\t\tFelix Matrix Header Other Parameters\n";
      cout << "\n   102. " << header[101] << "\tSegments of Bricks";
  
      cout << "\n   112. " << header[111] << "\tCharacters in Filename";
      vindex = 121;
      for(i=0; i<header[111]; i++)
        {
        cout << "\n   " << vindex+1 << ". " << header[vindex]
             << "\tCharacter in Filename: " << char(header[vindex]);
        vindex++;
        }

      for(kk=vindex; kk<999; kk++)
        {
        prnt = verbose - 9;
        if(header[kk]) prnt++;
        if(prnt>0)
        cout << "\n   " << kk+1 << ". " << header[kk] << "\tUnknown";
        }
      for(kk=999; kk<4094; kk++)
        {
        prnt = verbose - 9;
        if(header[kk]) prnt++;
        if(prnt>0)
        cout << "\n  " << kk+1 << ". " << header[kk] << "\tUnknown";
        }
      }

// 4095 ***************** Output Status Flag *********************** 4095 //

  cout << "\n\n\t\tFelix Matrix Header Status Flag\n";
  cout << "\n  4096. " << header[4095] << "\tHeader Status Flag";
  cout << "\n";
  return;
  }
*/


/*
ParameterSet Felix_p_set(int *header)

	// Input		header  : Array of Felix .dat parameters
	// Output		pset    : Parameter set with Felix parameters

  {
  ParameterSet pset;
  string pname;
  //string pdatac;
  int pdatai=0;
  double pdatad=0;
  string pstate;
  int i;
  for(i=0; i<32; i++)
    {
    pname = Fp_anac(i);			// Get the Felix .dat parameter anacronym
    pstate = Fp_String(i); 		// Get the Felix parameter comment
    switch(i)
      {
      case 0: case 1: case 2:
      case 3: case 4: case 5:
        pdatai = header[i];		// Add the parameter as an integer 
        pset.add(pname, pdatai, pstate); 
        break;
      case 6: case 7: case 8:
      case 9: case 10: case 11:
      case 12: case 13: case 14:
      case 15: case 16: case 17:
        pdatad = header[i];		// Add the parameter as a double
        pset.add(pname, pdatad, pstate); 
        break;
      case 18:
        pdatai = header[i];		// Add the parameter as an integer 
        pset.add(pname, pdatai, pstate); 
        break;
      case 19:
        pdatad = header[i];		// Add the parameter as a double
        pset.add(pname, pdatad, pstate); 
        break;
      case 20:
        pdatai = header[i];		// Add the parameter as an integer 
        pset.add(pname, pdatai, pstate); 
        break;
      case 21: case 22: case 23:
      case 24:
        pdatad = header[i];		// Add the parameter as a double
        pset.add(pname, pdatad, pstate); 
        break;
      case 25:
        pdatai = header[i];		// Add the parameter as an integer 
        pset.add(pname, pdatai, pstate); 
        break;
      case 26: case 27: case 28:
      case 29: case 30: case 31:
      default:
        pdatad = header[i];		// Add the parameter as a double
        pset.add(pname, pdatad, pstate); 
        break;
      }
    }
  return pset;
  }
*/


// ____________________________________________________________________________
//                    Felix HEADER AUXILIARY FUNCTIONS
// ____________________________________________________________________________


string Fp_String(int par)
  
	// Input	       par    :	Parameter number
	// Output	       string :	A string which defines the
	//				Felix .dat file parameter

  {
  switch(par)
    {
    case -2: return string("FORTRAN Integer"); 				break;
    case -1: return string("Storage Size in Complex Points");		break;
    case  0: return string("Number of Complex Points"); 		break;
    case  1: return string("Data Type"); 				break;
    case  2: return string("Transform State"); 				break;
    case  3: return string("W2 Axis Type"); 				break;
    case  4: return string("W2-W1 Axis Equivalence"); 			break;
    case  5: return string("W1 Axis Type"); 				break;
    case  6: return string("Empty"); 					break;
    case  7: return string("Empty"); 					break;
    case  8: return string("Empty"); 					break;
    case  9: return string("Empty"); 					break;
    case 10: return string("Pointer to Comments"); 			break;
    case 11: return string("Length of Comments"); 			break;
    case 12: return string("Pointer to Raw Header");	 		break;
    case 13: return string("Length of Raw Header"); 			break;
    case 14: return string("Spectrometer Type"); 			break;
    case 15: return string("Unknown");		 			break;
    case 16: return string("W2 Spectral Width (Hz.)"); 			break;
    case 17: return string("W2 Spectrometer Frequency (MHz.)");		break;
    case 18: return string("W2 Reference Point"); 			break;
    case 19: return string("W2 Reference Frequency (Hz.)");	 	break;
    case 20: return string("W2 Phase Correction Pivot Point");		break;
    case 21: return string("W2 Zero Order Phase Correction (Deg.)");	break;
    case 22: return string("W2 First Order Phase Correction"); 		break;
    case 23: return string("W1 Spectral Width (Hz.)"); 			break;
    case 24: return string("W1 Spectrometer Frequency (MHz.)");		break;
    case 25: return string("W1 Reference Point"); 			break;
    case 26: return string("W1 Reference Frequency (Hz.)"); 		break;
    case 27: return string("W1 Phase Correction Pivot Point");		break;
    case 28: return string("W1 Zero Order Phase Correction (Deg.)"); 	break;
    case 29: return string("W1 First Order Phase Correction (Deg.)");	break;
    case 30: return string("Reserved"); 				break;
    case 31: return string("Reserved"); 				break;
    default: return string("Unknown"); 					break;
    }
  return string("Unknown");
  }


string Fp_anac(int par)
  
	// Input	       par    :	Parameter number
	// Output	       string :	A string which is an anacronym
	//				for the Felix .dat file parameter

  {
  switch(par)
    {
    case  0: return string("DataPts");	break;	// Number of Complex Data Points
    case  1: return string("DataTyp");	break;	// Data Type
    case  2: return string("Trans1");	break;	// Transform Status Flag
    case  3: return string("AxTyp1");	break;	// W2 Axis Type          (1st Dim)
    case  4: return string("Axneq");	break;  // W2 - W1 Equivalence
    case  5: return string("AxTyp2");	break;	// W1 Axis Type          (2nd Dim)
    case  6: return string("Empty1"); 	break;
    case  7: return string("Empty2"); 	break;
    case  8: return string("Empty3"); 	break;
    case  9: return string("Empty4"); 	break;
    case 10: return string("P2Cmnt");	break;	// Pointer to Comments
    case 11: return string("LCmnts"); 	break;	// Length of Comments
    case 12: return string("P2Hder");	break;	// Pointer to Header
    case 13: return string("LHeadr"); 	break;	// Length of Header
    case 14: return string("Spectrm"); 	break;	// Spectrometer Type
    case 15: return string("Unknown"); 	break;
    case 16: return string("SW1"); 	break;	// W2 Spectral Width      (1st Dim)
    case 17: return string("SF1");	break;	// W2 Spectrometer Freq   (1st Dim)
    case 18: return string("RefPt1");	break;	// W2 Reference Point     (1st Dim)
    case 19: return string("RefFr1");	break;	// W2 Reference Frequency (1st Dim)
    case 20: return string("PivPt1");	break;	// W2 Phasing Pivot Point (1st Dim)	
    case 21: return string("ZerPhs1");	break;	// W2 Zero Order Phase    (1st Dim)	
    case 22: return string("FstPhs1");	break;	// W2 First Order Phase   (1st Dim)	
    case 23: return string("SW2");	break;	// W1 Spectral Width	  (2nd Dim)
    case 24: return string("SF2");	break;  // W1 Spectrometer Freq   (2nd Dim)
    case 25: return string("PivPt2");	break;	// W1 Reference Point	  (2nd Dim)
    case 26: return string("RefFr2");	break;	// W1 Reference Frequency (2nd Dim)
    case 27: return string("PivPt2");	break;	// W1 Phasing Pivot Point (2nd Dim)
    case 28: return string("ZerPhs2");	break;	// W1 Zero Order Phasing  (2nd Dim)
    case 29: return string("FstPhs2");	break;	// W1 First Order Phasing (2nd Dim)
    case 30: return string("Reserv1");	break;	// Felix Reserved
    case 31: return string("Reserv2");	break;	// Felix Reserved
    default: return string("Unknown");	break;
    }
  return string("Unknown");
  }


string Fp_String_qual(int pos, float parval)
  
	// Input	       par    :	Parameter number
	//		       parval : Parameter value
	// Output	       string :	A string which units and/or
	//				specifics of Felix .dat file parameter

  {
  string pstring;
  union fparam fpar;			// A general Felix parameter
  fpar.ff = parval;
  switch(pos)
    {
    case 1: 				// Type of data
      if(fpar.fi)
        pstring = ": Complex";
      else
        pstring = ": Real"; 
      break;
    case 2: 				// Transform State
      if(fpar.fi)
        pstring = ": Spectrum";
      else
        pstring = ": FID"; 
      break;
    case 3: 				// W2 Axis Type 
    case 5: 				// W1 Axis Type 
      if(fpar.fi == 0)
        pstring = ": None";
      else if(fpar.fi == 1)
        pstring = ": Points";
      else if(fpar.fi == 2)
        pstring = ": Hertz";
      else if(fpar.fi == 3)
        pstring = ": PPM";
      else if(fpar.fi == 4)
        pstring = ": Sec";
      else if(fpar.fi == 5)
        pstring = ": Cm-1";
      else
        pstring = ": Unknown";
        break;
    default:
      pstring = "";
      break;
    }
  return pstring;
  }

#endif							// Felix.cc
