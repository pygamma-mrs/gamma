/* XWinSer.h ****************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      XWinSer                                    Interface		**
**                                                                      **
**      Copyright (c) 1999                                              **
**      Scott Smith                                                     **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**                                                                      **
**      $Header:  $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**      Description                                                     **
**                                                                      **
** The XWin* files provide an interface to Bruker XWinNMR (uxnmr) data  **
** sets. This class embodies a Bruker binary time domain file, usually  **
** named ser, associated with a multi-dimension NMR experiment. As is	**
** implied, such serial files simply contain stacked FIDs. However, the	**
** data points are stored as 32-bit integers ordered re,im,re,im,......	**
** Important aspects are 1.) How many points are in each block, 	**
** 2.) What is the byte order of stored data is,  and 3.) How many	**
** blocks are stored in the file. Such information is normally found	**
** in associated parameter files, acqus and acqu2s. Below is a typical	**
** Bruker directory structure associated with a 2D acquisiton.		**
**									**
**                          __ acqu, acqu2 (changable parameter files)	**
**			   / 						**
**                        /___ acqus, acqu2s (static parameter files)	**
**			 /						**
**  expname -- expnum --< ---- ser (binary data)			**
**			 \						**
**			  \___ pdata -- 1 -- proc, proc2, procs, ....   **
**									**
** This class would then handle the file ser.  It will NOT deal with	**
** the directory hierarchy shown above nor with any associated ASCII	**
** parameter files.  That is left up to higher level classes supporting **
** XWinNMR. 								**
**                                                                      **
** Note that this is a DISK entity, not a row_vector. That is, when	**
** this is instructed to write a ser file it must either be given a	**
** complete data matrix (which probably uses too much memory) or a	**
** series of row vectors which are written to disk as provided. The 	**
** class will also provide some basic info as to how the files is,	**
** or should be, written. Similarly, it can read any block of the 	**
** serial file (or the entire array.)					**
**                                                                      **
*************************************************************************/

#ifndef   XWinSer_H_ 				// Is file already included?
#  define XWinSer_H_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// This is the implementation
#  endif

#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <fstream>                              // File streams from libstdc++
#include <string>				// Include libstdc++ strings
#include <Matrix/row_vector.h>			// Know about row vectors
#include <Matrix/matrix.h>			// Know about matrices

class XWinSer

// ----------------------------------------------------------------------------
// ------------------------------- STRUCTURE ----------------------------------
// ----------------------------------------------------------------------------
 
  {
  std::fstream ffp;		// Input/Output binary file stream
  std::string  ffname;		// Bruker Ser/Serial file name
  bool	       fbigend;		// Flag for big- vs. little-endian (this arch)
  bool         fbyteordin;	// Flag for big- vs. little-endian input
  int          ftotpts;		// Total points per Ser (real + imag)    
  int          fpadding;	// Bytes needed for padding to 256 pt boundary
  int          ffsize;		// Our total file size (bytes)
 
// ----------------------------------------------------------------------------
// ----------------------------- PRIVATE FUNCTIONS ----------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                      XWinNMR Ser File Error Handling
// ____________________________________________________________________________

/* These functions take care of any errors encountered when reading, writing,   
   and setting parameters in Bruker acquisition parameter files.
 
	Input		SerPar  : XWinNMR acqusition parameters (this)
			eidx    : Error index
                        pname   : Additional error message
                        noret   : Flag for linefeed (0=linefeed)
        Output          void    : An error message is output                 */

         void XWinSererror(int    eidx, int noret=0) const;
         void XWinSererror(int    eidx, const std::string& pname, int noret=0) const;
volatile void XWinSerfatality(int eidx) const;
volatile void XWinSerfatality(int eidx, const std::string& pname) const;


// ----------------------------------------------------------------------------
// ii                      Mode Checking Functions
// ----------------------------------------------------------------------------

/* These functions check the mode specified when opening a Bruker binary file.
   We do not allow all of the possible modes availble to C++ filestreams.  In
   particular, we do not allow users to append to an existing file (since this
   would be bad if the existing data was written with a different byte order).
   Additonally, there seems to be some discrepancy between compilers as to
   what modes are available to fstream, so we just stick to the basics.

        in  = 1 (0x01)  out       = 2  (0x02)   ate      = 4   (0x03)
        app = 8 (0x08)  trunc     = 16 (0x10)   nocreate = 32  (0x20)
                        noreplace = 64 (0x40)   binary   = 128 (0x80)

   Use of noreplace & nocreate maybe didn't make it to ANSI standard...      */  

void CheckMode(int mode);
int  CheckSize(int warn=2);

 
// ----------------------------------------------------------------------------
// iii                 File Padding & Boundary Functions
// ----------------------------------------------------------------------------
 
/* Regardless of the block size, a Bruker serial file (ser) enforces the rule
   that each block must begin on a 1024 byte block boundary (256 points).  So,
   when writing serial files we will pad the output file with zeros until the
   boundary is reached. Similarly, we'll check that all blocks are all
   written into the file beginning at such a boundary.                       */  
 
void SetPadding();
bool CheckBoundary();
void SkipPadding();
void AddPadding();

// ----------------------------------------------------------------------------
// iv                    File Default Settings
// ----------------------------------------------------------------------------

void SetDefaults();


public:
// ____________________________________________________________________________ 
// A                  XWinSer File Constructors, Destructor
// ____________________________________________________________________________

/* These are the constructors of the class handling Bruker XWinNMR serial
   binary files.  The default name of such files is ser.  The output
   byte order is automatically set depending upon the computer architecture.

        //                      name    : External filename 
        //                      mode    : I/O mode 
        //                                 app    = append 
        //                                 ate    = open & seek end (at end) 
        //                                 binary = I/O in binary
        //                                 in     = open for reading 
        //                                 out    = open for writing 
        //                                 trunc  = truncate file at 0 length 
*/

XWinSer();
XWinSer(const std::string& name, int TD, bool byteord,
                                      int mode = std::ios::binary|std::ios::in);
XWinSer(const XWinSer& XWF);
virtual ~XWinSer();                                                             
virtual XWinSer& operator=(const XWinSer& XWF);


// ____________________________________________________________________________
// B                   XWinNMR Serial File Access Functions
// ____________________________________________________________________________
 
/* These functions allow users to get some simple intformation regarding the
   contents of the Bruker data serial file.                                  */

	std::string sername() const;		// File name
virtual int         size()    const;		// No. complex points
	int         TDS()     const;		// No. total points
virtual bool        order()   const;		// Data byte order
virtual int         bytes()   const;		// File size in bytes
virtual int         blocks()  const;		// No. Sers
virtual int         pad()     const;		// No. padding bytes

// ____________________________________________________________________________
// C                       XWinSer Input Functions
// ____________________________________________________________________________
 
/* These functions will read in a single block or a series of blocks from the
   Bruker serial file. The block points are read in as re,im,re,im,...., each
   point stored as a 32-bit integer.  It is very important that we know the
   number of points in each block and the byte order used in storing the data.
   This information must be supplied, either when calling the function(s) or
   when constructing objects of this class type.  Such information is usually
   read in from corresponding ASCII parameter files, typically named acqus 
   and acqu2s.

     Function   Arguments              Action                       Return
     --------  ------------  -----------------------------------  ----------
     readFID   fn,TD,BO,I    Read ith block in serial/FID file    row_vector
   + readFID   idx           Read block of index idx (-1=current) row_vector
     readFIDs  fn,TD,BO,I,J  Read J FIDs starting from Ith one      matrix
   + readFIDs  I,J           Read J FIDs starting from Ith one      matrix
     readSer   fn,TD,bytord  Read entire serial file                matrix
   + readSer                 Read entire serial file                matrix

 Above X=Currently Inactive, +=File Must Be Properly Opened For Reading      */

virtual row_vector readFID(const std::string& F,int TD,int BO=-1,int I=-1);
row_vector readFID(int indx=-1);
matrix     readFIDs(const std::string& name, int TD, int BO, int I, int J);
matrix     readFIDs(int I, int J);
matrix     readSer(const std::string& name, int TD, int byteord=-1);
matrix     readSer();
row_vector readSlice(const std::string& fin, int TD, int BO=-1, int I=-1);
row_vector readSlice(int idx);
row_vector readSlices(const std::string& F, int TD, int BO, int I, int J);
row_vector readSlices(int idx, int NB);


// ____________________________________________________________________________
// D                       XWinSer Output Functions
// ____________________________________________________________________________

/* These functions will write a serial file in Bruker XWinNMR format. This can
   be done either by directly writing a matrix into such a file, or by writing
   a series of row vectors. Note that usually there are parameter file to go
   along with such data files, normally called acqus and acqu2s.             */

virtual int write(const std::string& F, const row_vector& data, int warn=2);
virtual int write(const row_vector& data, int warn=2);
virtual int write(const std::string& F, const matrix& data, int warn=2);
virtual int write(const matrix& data, int warn=2);

// ____________________________________________________________________________
// E                    XWinSer Standard Output Functions
// ____________________________________________________________________________
 
/* These functions allow users to have a quick peek at the "fid" contents.
   They do not do any manipulations to the file whatsoever, they will only
   report on what in known about the data values and file sizes.             */ 

std::ostream& print(std::ostream& ostr, int full=0, int hdr=1);
friend  std::ostream& operator<< (std::ostream& O, XWinSer& F);

// ____________________________________________________________________________
// F                  XWinSer File Handling Functions
// ____________________________________________________________________________

bool open(const std::string& name, int TD, bool byteord,
                         int mode = std::ios::binary|std::ios::in, int warn=2);


};

#endif 								// XWinSer.h
