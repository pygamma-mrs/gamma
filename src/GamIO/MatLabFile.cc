/* MatLabFile.cc ************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A                           **
**                                                                      **
**      MATLAB File                                Implementation	**
**                                                                      **
**      Copyright (c) 1990                                              **
**      Scott Smith                                                     **
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
** The MatLabFile class facilitates reading and writing of MATLAB MAT   **
** binary files. MATLAB is a commercial product and can be purchased    **
** from The MathWorks, Inc.                                             **
**                                                                      **
*************************************************************************/

#ifndef _MLFILE_cc_			// Is this file already included?
#define _MLFILE_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler
#    pragma implementation
#endif

#if defined(_MSC_VER)			// If using MSVC++ then we
 #pragma warning (disable : 4786)       // Kill STL namelength warnings
#endif 

#include <GamIO/MatLabFile.h>			// Include interface
#include <Basics/Gutils.h>			// Include GAMMA errors
#include <Basics/StringCut.h>			// Include GAMMA Gdec function
#include <GamIO/ML5Hdr.h>			// Include MAT V5 headers
#include <GamIO/ML5DElem.h>			// Include MAT V5 data elements
#include <GamIO/ML4DElem.h>			// Include MAT V4 data elements
#include <GamIO/BinIOBase.h>			// Include Int2Mode function
#include <string>				// Include libstdc++ strings

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ----------------------------------------------------------------------------
// i                         MatLab File Error Handling
// ____________________________________________________________________________

/*              Input           MLF	: MatLab MAT file
                                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */  

void MatLabFile::MLFerror(int eidx, int noret) const
  {
  std::string hdr("MATLAB MAT File");
  switch(eidx)
    {
    case 0: GAMMAerror(hdr, 0, noret); break;   // Program Aborting        (0)
    case 10: GAMMAerror(hdr,"Cannot Read Data Element Header", noret);   //(10)
      break;
    case 20: GAMMAerror(hdr,"Cannot Open File In ASCII Mode", noret);	 //(20)
      break;
    case 21: GAMMAerror(hdr,"These Files Must Be Binary", noret);	 //(21)
      break;
    case 22: GAMMAerror(hdr,"Open For Input & Output Disallowed", noret);//(22)
      break;
    case 23: GAMMAerror(hdr,"Open As Either Input Or Output", noret);	 //(23)
      break;
    case 24: GAMMAerror(hdr,"Open At End Of File Is Disallowed", noret); //(24)
      break;
    case 25: GAMMAerror(hdr,"Open Without Using ios::ate", noret);	 //(25)
      break;
    case 26: GAMMAerror(hdr,"Open File For Append Is Disallowed", noret);//(26)
      break;
    case 27: GAMMAerror(hdr,"Open Without Using ios::app", noret);	 //(27)
      break;
    default: GAMMAerror(hdr,eidx, noret); break;// Usually Unknown Error  (-1)
      break;
    }
  }
     

void MatLabFile::MLFerror(int eidx, const std::string& pname, int noret) const
  {
  std::string hdr("MATLAB MAT File");
  switch(eidx)
    {
    case 1:  GAMMAerror(hdr,   1, pname, noret); break; // File Problems   (1)
    case 2:  GAMMAerror(hdr,   2, pname, noret); break; // !Read SingPar   (2)
    case 10: GAMMAerror(hdr,"Can't Read Array "+pname, noret); break;    //(10)
    default: GAMMAerror(hdr,  -1, pname, noret); break; // Unknown Error   (-1)
    }
  }
     
 
volatile void MatLabFile::MLFfatality(int eidx) const
  {
  MLFerror(eidx, 1);
  if(eidx) MLFerror(0);
  GAMMAfatal();					// Clean exit from program
  }
     
 
volatile void MatLabFile::MLFfatality(int eidx, const std::string& pname) const
  {
  MLFerror(eidx, pname, 1);
  if(eidx) MLFerror(0);
  GAMMAfatal();					// Clean exit from program
  }



// ----------------------------------------------------------------------------
// ii                       Mode Checking Functions
// ----------------------------------------------------------------------------

/* These functions check the mode specified when opening a MATLAB MAT file.
   We do not allow all of the possible modes availble to C++ filestreams.  In
   particular, we do not allow users to append to an existing file (since this
   would be bad if the existing data was written with a different byte order).
   Additonally, there seems to be some discrepancy between compilers as to
   what modes are available to fstream, so we just stick to the basics.

	in  = 1 (0x01)	out       = 2  (0x02)	ate      = 4   (0x03)
	app = 8	(0x08)	trunc     = 16 (0x10)	nocreate = 32  (0x20)
			noreplace = 64 (0x40) 	binary   = 128 (0x80)
   
   Use of noreplace & nocreate maybe didn't make it to ANSI standard...      */

// sosi - what about trunc?

void MatLabFile::CheckMode(int mode)
  {
  if(!(mode & std::ios::binary))
    {
    MLFerror(20, 1);			// File can't be ASCII
    MLFfatality(21);			// File must be binary
    }
  if((mode & std::ios::out) && (mode & std::ios::in))
    {
    MLFerror(22, 1);			// File can't be input & output
    MLFfatality(23);			// Open file as either input or output
    }
  if(mode & std::ios::ate)
    {
    MLFerror(24, 1);			// File can't open at end
    MLFfatality(25);			// Must be opened without ios::ate
    }
  if(mode & std::ios::app)
    {
    MLFerror(26, 1);			// File can't open in append mode
    MLFfatality(27);			// Must be opened with ios:app
    }
  if(mode & std::ios::out)     wxr = -1;	// This flag says we are outputting
  else if(mode & std::ios::in) wxr = 1;	// This flag says we are inputting
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                    MatLabFile CONSTRUCTORS, DESTRUCTOR
// ____________________________________________________________________________


MatLabFile::MatLabFile() { MATVerOut=5; fsize=0; wxr=0; fname=""; }

MatLabFile::MatLabFile(const std::string& name, int mode)
           :fp(name.c_str(), Int2Mode(mode))

        // Input              	MLF	: MATLAB MAT file
	//                      name	: External filename
	//			mode	: I/O mode
	//				   binary = I/O in binary
	//				   in     = open for reading
	//				   out    = open for writing
        // Output                MLF	: MATLAB MAT file created/opened/etc.
	// Note				: There is NO data written/read in
	//				  from this routine!  Truncation,
	//				  etc. may delete the file though

  {
  if(!fp.good())			// Insure MATLAB File is O.K.
    {
    if(wxr>1)				// Errors if open for reading
      {
      MLFerror(1,1);			//   Problems with input filestream
      MLFerror(1,name,1);		//   Problems with file
      MLFfatality(14);			//   Can't read a damn thing
      }
    else
      {
      MLFerror(2,1);			//   Problems with output filestream
      MLFerror(1,name,1);		//   Problems with file
      MLFfatality(15);			//   Can't write a damn thing
      }
    }
  CheckMode(mode);			// Check if mode is acceptable
  MATVerOut = 5;			// Set default MAT output version
  fname = name;				// Set the MAT filename
  fp.seekp(0,std::ios::end);			// Go to the end of the file
  fsize = fp.tellp();			// See how many bytes are in the file
  fp.seekp(0);				// Set things back to the start
  }

MatLabFile::MatLabFile(const MatLabFile& MLF)
  {
//sosi
//  fp        = MLF.fp;			// Input/Output file stream
  fname     = MLF.fname;		// Filename associated with MatLabFile
  fsize     = MLF.fsize;		// Size of the file (in bytes)
  ML5Hdr    = MLF.ML5Hdr;		// MatLab MAT Version 5 main header
  MATVerOut = MLF.MATVerOut;		// MatLab MAT Version (for output)
  wxr       = MLF.wxr;			// Flag for write/nothing/read (-1,0,1);
  }

MatLabFile::~MatLabFile() {}				


MatLabFile& MatLabFile::operator=(const MatLabFile& MLF)
  {
  if(this == &MLF) return *this;	// Avoid self copying
//sosi
//  fp        = MLF.fp;			// Input/Output file stream
  fname     = MLF.fname;		// Filename associated with MatLabFile
  fsize     = MLF.fsize;		// Size of the file (in bytes)
  ML5Hdr    = MLF.ML5Hdr;		// MatLab MAT Version 5 main header
  MATVerOut = MLF.MATVerOut;		// MatLab MAT Version (for output)
  wxr       = MLF.wxr;			// Flag for write/nothing/read (-1,0,1);
  return *this;
  }


void MatLabFile::close()
  {
  fp.close();				// Close the file stream
  fsize = 0;				// The file size is zero noe
  wxr = 0;				// We are not read nor write
  }
 

void MatLabFile::open(int mode)
  {
  fp.close();				// Close the file stream
  if(!fname.length())			// We cannot open the file without
    {					// an external file name 
    }
  fp.open(fname.c_str(),Int2Mode(mode));// Open the file stream
  if(!fp.good())			// Insure MATLAB File is O.K.
    {
    if(wxr>1)				// Errors if open for reading
      {
      MLFerror(1,1);			//   Problems with input filestream
      MLFerror(1,fname,1);		//   Problems with file
      MLFfatality(14);			//   Can't read a damn thing
      }
    else
      {
      MLFerror(2,1);			//   Problems with output filestream
      MLFerror(1,fname,1);		//   Problems with file
      MLFfatality(15);			//   Can't write a damn thing
      }
    }
  CheckMode(mode);			// Check if mode is acceptable
  fp.seekp(0,std::ios::end);			// Go to the end of the file
  fsize = fp.tellp();			// See how many bytes are in the file
  fp.seekp(0);				// Set things back to the start
  }

 
// ____________________________________________________________________________
// B                        MatLab File Access Functions
// ____________________________________________________________________________

/* These functions allow access to the MATLAB MAT version.  If the file is 
   open for reading the version returned reflects the file on disk.  If the
   file is either unopened or opened for writing the version refects that
   which is written. Setting the version affect only the output format.      */ 

int  MatLabFile::Version() const      { return MATVerOut; }
void MatLabFile::Version(int V)       { MATVerOut = V; }
int  MatLabFile::Version(const std::string& filename)
  {
  MatLabFile MLF(filename);		// Open the MatLab file
  MLF.fp.seekp(0);			// Move to file start
  char X;				// Done by charactor input
  int ver = 5;				// Assume version 5 format
  for(int i=0; i<4 && ver==5; i++)	// Look in 1st 4 bytes of file
    {					// and if any of them are zero then
    MLF.fp.read(&X, sizeof(char));	// we are dealing with a version 4
    if(!int(X)) ver = 4;		// MAT file (we only know V4 & V5)
    }
  return ver;
  }


int  MatLabFile::Version(std::fstream& fp) 
  {
  if(!fp) return -1;			// Return crazy if no file open
  if(!fsize) return -1;			// Return crazy if file size is zero
  int pos =  fp.tellp();		// Get current position
  fp.seekp(0);				// Move to file start
  char X;				// Done by charactor input
  int ver = 5;				// Assume version 5 format
  for(int i=0; i<4 && ver==5; i++)	// Look in 1st 4 bytes of file
    {					// and if any of them are zero then
    fp.read(&X, sizeof(char));		// we are dealing with a version 4
    if(!int(X)) ver = 4;		// MAT file (we only know V4 & V5)
    }
  fp.seekp(pos);			// Reset original file position
  return ver;
  }
 
// ____________________________________________________________________________
// C                   MatLab File Auxiliary Functions
// ____________________________________________________________________________


void MatLabFile::whos(const std::string& filename, std::ostream& ostr)
  {
  MatLabFile MLF(filename);		// Open the MatLab file
  MLF.whos(ostr);			// Use function overload
  }

void MatLabFile::whos(std::ostream& ostr)

        // Input              	MLF	: MATLAB MAT file
	//                      ostr	: An output stream
	// Output		void	: This adds information about the
	//				  variables in the MAT file into
	//				  the output stream.
	// Note				: If the file is empty then no
	//				  output is produced
	// Note				: Current file position unaltered

  {
  if(!fp) return;			// Do nothing if no file open
  if(!fsize) return;			// No output if file size is zero
  std::string margin = std::string(12, ' ');
  std::string spcr = std::string(2, ' ');
  std::string ul10 = std::string(10, '-'); 
  ostr << "\n" << margin
       << "   Name   " << spcr
       << "   Size   " << spcr
       << "   Bytes  " << spcr
       << "          Class";
    ostr << "\n" << margin
         << ul10 << spcr
         << ul10 << spcr
         << ul10 << spcr
         << ul10 << ul10 << "----";
  int pos =  fp.tellp();		// Get current position
  fp.seekp(0);				// Move to file start
  if(Version(fp) == 4) 			// No main header in version 4 files
    {					// so we begin by outputting the
    MatLab4DE ML4DE;			// A version 4 MAT data element
    int bytesleft = fsize;
    while(fp.is_open() && bytesleft)
      {
      ML4DE.whos(ostr, fp);		// Output line of "whos" info
      ML4DE.skip(fp);			// Skip to next data element
      bytesleft = fsize - fp.tellp();
      }
    }
  else					// For version 5 MAT files we must
    {					// skip the header and then start
    int bytesleft = fsize;
    MatLab5Hdr MH;			// Set up temporary header
    MH.read(fp);			// Read the header from the file start
    MatLab5DE ML5DE;			// A version 5 MAT data element
    bytesleft = fsize - fp.tellp();
    bool bigend = (MH.BigEndIn)?true:false;
    while(fp.is_open() && bytesleft>0)
      {
      ML5DE.whos(ostr,fp,bigend);	// Output line of "whos" info
      ML5DE.skip(fp, bigend);	// Skip to next data element
      bytesleft = fsize - fp.tellp();
      }
    }
  fp.seekp(pos);			// Reset original file position
  }


void MatLabFile::details(const std::string& filename, std::ostream& ostr)
  {
  MatLabFile MLF(filename);		// Open the MatLab file
  MLF.details(ostr);			// Use overload function
  }


void MatLabFile::details(std::ostream& ostr)
  {
  if(!fp) return;			// Do nothing if no file open
  if(!fsize) return;			// No output if file size is zero
  int pos =  fp.tellp();		// Get current position
  if(Version(fp) == 4) 			// No main header in version 4 files
    {					// so we begin by outputting the
    fp.seekp(0);			// Move to file start
    MatLab4DE ML4DE;			// A version 4 MAT data element
    int bytesleft = fsize;
    while(fp.is_open() && bytesleft)
      {
      ML4DE.whos(ostr, fp);		// Output line of "whos" info
      ML4DE.skip(fp);			// Skip to next data element
      bytesleft = fsize - fp.tellp();
      }
    }
  else					// For version 5 MAT files we must
    {					// skip the header and then start
    fp.seekp(0);			// Move to file start
    int bytesleft = fsize;
    MatLab5Hdr MH;			// Set up temporary header
    MH.read(fp);			// Read the header from the file start
    bytesleft = fsize - fp.tellp();
    MatLab5DE ML5DE;			// A version 5 MAT data element
    bool bigend = (MH.BigEndIn)?true:false;
    while(fp.is_open() && bytesleft>0)
      {
      ML5DE.read(fp,bigend);	// Read in the data element
      ML5DE.print(ostr);		// Output its gory details
      bytesleft = fsize - fp.tellp();
      }
    }
  fp.seekp(pos);			// Reset original file position
  }


void MatLabFile::Header(const std::string& filename, std::ostream& ostr)
  {
  MatLabFile MLF(filename);		// Open the MatLab file
  if(MLF.Version(MLF.fp) == 4) return;	// No main header in version 4 files
  MatLab5Hdr MH;			// Set up temporary header
  MH.read(MLF.fp);			// Read the header
  MH.print(ostr);			// Output main header
  }


// ____________________________________________________________________________
// D                  MatLab File (Binary) Read Functions
// ____________________________________________________________________________
         
/* These functions will access a MATLAB file and read one of its array data
   elements into a GAMMA array. Since MATLAB arrays are not distinguished
   from vectors, the read always produces matrices as well.  Users must do
   a type cast to produce the data in a row or column vector.                */

matrix MATLAB(const std::string& filename, const std::string& name, int warn) 
  {
  MatLabFile MLF(filename, std::ios::binary|std::ios::in);// Open a MATLAB File
  matrix mx = MLF.GetMatrix(name, warn);	// Return matrix name
  if(!mx.rows() && warn)
    {
    MLF.MLFerror(1, filename, 1);		// Problems with filename
    if(warn > 1) MLF.MLFfatality(10,name);	// Can't read requested array
    }
  return mx;
  }

matrix MatLabFile::GetMatrix(const std::string& name, int warn)
  {
  matrix mx;				// For our return matrix
  if(!fp || !fsize)			// Insure the file is OK
    {
    if(warn)
      {
      MLFerror(1, 1);			// Problems with input filestream
      if(warn > 1) MLFfatality(10,name);// Can't read requested array
      else         MLFerror(10, 10); 
      }
    return mx;
    }
  int pos =  fp.tellp();		// Get current position
  fp.seekp(0);				// Move to file start
  if(Version(fp) == 4) 			// No main header in version 4 files
    {					// so we begin by outputting the
// sosi this all needs work
    MatLab4DE ML4DE;			// A version 4 MAT data element
    int bytesleft = fp.tellp();
    while(fp.is_open() && bytesleft)
      {
//      ML4DE.whos(ostr, fp);		// Output line of "whos" info
      ML4DE.skip(fp);			// Skip to next data element
      bytesleft = fsize - fp.tellp();
      }
    }
  else
    {
    MatLab5Hdr MH;			// Set up temporary header
    MH.read(fp);			// Read the header from the file start
    MatLab5DE ML5DE;			// A version 5 MAT data element
    std::string DEname;			// Data element name
    int bytesleft = fp.tellp();		// This many bytes left in file
    bool bigend = (MH.BigEndIn)?true:false;		// Flag for input endian type
    while(fp.is_open() && bytesleft)
      {
      DEname = ML5DE.Name(fp,bigend);	// Get name of current data eleemnt
      if(name != DEname)		// If not the element being sought
        ML5DE.skip(fp, bigend);		// then skip to next data element
      else
        {				// This is the element being sought
        mx = ML5DE.GetMatrix(fp,bigend);// Read in the matrix
        fp.seekp(pos);			// Reset original file position
        return mx;
        }
      bytesleft = fsize - fp.tellp();	// Determine bytes left in the file
      }
    if(warn)
      {
      if(warn > 1) MLFfatality(10,name);// Can't read requested array
      else         MLFerror(10, 10); 
      }
    }
  fp.seekp(pos);			// Reset original file position
  return mx;				// This will be a dummy matrix
  }

// ____________________________________________________________________________
// E                   MatLab File (Binary) Write Functions
// ____________________________________________________________________________
 
/* These functions will access/create a MATLAB file and write a data element
   from a GAMMA array.                                                       */

int MATLAB(const std::string& fileout, const std::string& name,
                                            const matrix& mx, int rc, int warn) 
  {
  MatLabFile MLF(fileout,std::ios::binary|std::ios::out);	// Open a MATLAB File
  int TF = MLF.write(name, mx, rc, warn?1:0);	// Write matrix with set name
  (MLF.fp).close();				// Close the file
  return TF;
  }

int MatLabFile::write(const std::string& N, const matrix& mx, int rc, int warn)
  {
  int TF = 1;
  if(MATVerOut == 5)
    {
    if(!fsize) TF *= ML5Hdr.write(fp, warn?1:0);// Write main header if start
    MatLab5DE DE(mx,N,rc);			// Construct a data element
    DE.write(fp);				// Write the data element
    fsize = fp.tellp();				// Reset the file size
    }
  else
    {
    //    MatLab5DE DE(mx,N,rc);			// Construct a data element
    //    DE.write(fp, mx, N);			// Write the matrix
    // sosi
    return 0;
    }
  return 1;
  }


int MATLAB(const std::string& fileout, const std::string& name,
                                               const coord_vec& cvec, int warn) 
  {
  MatLabFile MLF(fileout,std::ios::binary|std::ios::out);	// Open a MATLAB File
  int TF = MLF.write(name, cvec, warn?1:0);	// Write cvec with set name(s)
  (MLF.fp).close();				// Close the file
  return TF;
  }

int MatLabFile::write(const std::string& N, const coord_vec& cvec, int warn)
  {
  int TF = 1;
  if(MATVerOut == 5)
    {
    if(!fsize) TF *= ML5Hdr.write(fp, warn?1:0);// Write main header if start
    int npts = cvec.size();			// Get vector size
    row_vector X(npts),Y(npts),Z(npts);
    for(int i=0; i<npts; i++)
      {
      X.put(cvec.x(i), i);
      Y.put(cvec.y(i), i);
      Z.put(cvec.z(i), i);
      }
    MatLab5DE DEX(X,N+Gdec(0),0);		// Construct a data element
    DEX.write(fp);				// Write the data element
    MatLab5DE DEY(Y,N+Gdec(1),0);		// Construct a data element
    DEY.write(fp);				// Write the data element
    MatLab5DE DEZ(Z,N+Gdec(2),0);		// Construct a data element
    DEZ.write(fp);				// Write the data element
    fsize = fp.tellp();				// Reset the file size
    }
  else
    {
//    MatLab5DE DE(mx,N,rc);			// Construct a data element
//    DE.write(fp, mx, N);			// Write the matrix
// sosi
return 0;
    }
  return 1;
  }

// ____________________________________________________________________________
// F                     MatLab File (ASCII) Print Functions
// ____________________________________________________________________________
 
/* These functions are used to print out the contents of MATLAB MAT files.
   They will NOT alter the contents of such files and will close the file if
   it is open for writing.

	Input			MLF     : MAT version 5 file (this)
 				ostr    : An output stream
				full    : Flag for full output
	Output			ostr    : Output stream with MATLAB MAT file
                                          base info placed into it           */

void MatLabFile::print(const std::string& filename, std::ostream& ostr, int full)
  {
  MatLabFile MLF(filename);		// Open the MatLab file
  MLF.print(ostr, full);			// Use function overload
  }


void MatLabFile::print(std::ostream& ostr, int full)
  {
  ostr << std::string(34,' ') 			// Print a title
       << "MATLAB Binary MAT File\n";
  ostr << "\n\t\tMAT File Name   ";		// Print the file name
  if(!fname.length()) ostr << "None";
  else                ostr << fname;
  ostr << "\n\t\tMAT File Size   " 		// Print the file size
       << fsize << " Bytes";
  if(!fp)					// Print the default file
    {						// version if no file exists
    ostr << "\n\t\tMAT Status      Not Open";
    return;
    }
  int pos =  fp.tellp();			// Get current position
  int vers = Version(fp);			// Get MAT version 
  ostr << "\n\t\tMAT Version     " << vers;	// Output the MAT version
  if(vers == 5)					// MAT Version 5 Output
    {
    if(fsize)					// If there is data in the
      {						// file we always use that
      fp.seekp(0);				// Start at file beginning
      ML5Hdr.read(fp,1);			// Read version 5 header 
      }
    std::cout << "\n\t\tMAT Header";
    ML5Hdr.print(ostr, 0);			// Output main header
    std::cout << "\n";				// Add a blank line
    if(full) details(ostr);
    else     whos(ostr);			// Do MATLAB "whos" output
    }
  else if(vers == 4)				// MAT Version 4 Output
    {
    std::cout << "\n";				// Add a blank line
    whos(ostr);					// Do MATLAB "whos" output
    }
  fp.seekp(pos);				// Back to original position
  }


std::ostream& operator<< (std::ostream& ostr, MatLabFile& MLF)
  { MLF.print(ostr); return ostr; }


// ____________________________________________________________________________
// G              MatLab File (ASCII) Auxiliary Print Functions
// ____________________________________________________________________________
 
/* These functions are just used for testing purposes.  The mirror the binary
   write functions, but instead of writing in binary to a filestream the write
   in ASCII to an output stream.  Thus it is easy to see most of what is
   being written, at least in principle.....  The function overload using
   class matrix should also suffice for GAMMA row_vector and col_vector.     */

void MatLabFile::print(std::ostream& ostr,
                           const matrix& mx, const std::string& N, int cmplx) const 
  {
  if(MATVerOut == 5)
    {
    ML5Hdr.print(ostr);				// Output main header
    MatLab5DE DE(mx,N,cmplx);			// Construct mx data element
    DE.print(ostr);				// Output the data element 
    }
  else
    {
    MatLab4DE DE(N,mx,cmplx);			// Output main header
    DE.print(ostr);				// Output the data element 
    }
  }

#endif							// MatLabFile.cc
