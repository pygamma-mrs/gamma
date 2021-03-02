/* XWinFid.h ****************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      XWinFid                                    Interface		**
**                                                                      **
**      Copyright (c) 1999                                              **
**      Scott Smith                                                     **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**      Description                                                     **
**                                                                      **
** The XWin* files provide an interface to Bruker XWinNMR (uxnmr) data  **
** sets. This class embodies a Bruker binary data file associated, fid,	**
** associated with a 1D acquisition. There is not anything fancy about	**
** such files.  They contains 32-bit integers ordered re,im,re,im,....	**
** Two important aspects are 1.) How many points are in the FID, and	**
** 2.) What is the byte order of the stored data. That information is	**
** normally contained in an associated parameter file, acqus. Below	**
** is a typical Bruker directory structure associated with a 1D		**
** acquisiton.								**
**                          __ acqu  (changable parameter file)		**
**			   / 						**
**                        /___ acqus (static parameter file)		**
**			 /						**
**  expname -- expnum --< ---- fid (binary data)			**
**			 \						**
**			  \___ pdata -- 1 -- proc, procs, meta, outd	**
**									**
** This class would then handle the file "fid".  It will NOT deal with	**
** the directory hierarchy shown above nor with any associated ASCII	**
** parameter files.  That is left up to higher level classes supporting **
** XWinNMR. 								**
**                                                                      **
** Note that this is a DISK entity, not a row_vector. That is, when	**
** this is instructed to write an FID it must be given a row_vector &	**
** some basic info as to how it should be written.  Similarly, it can	**
** read the file fid into a row_vector if needed.			**
**                                                                      **
*************************************************************************/

#ifndef   XWinFid_H_ 				// Is file already included?
#  define XWinFid_H_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// This is the implementation
#  endif

#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <fstream>                              // File streams from libstdc++
#include <string>				// Include libstdc++ strings
#include <Matrix/row_vector.h>			// Know about row vectors

class XWinFid

// ----------------------------------------------------------------------------
// ------------------------------- STRUCTURE ----------------------------------
// ----------------------------------------------------------------------------
 
  {
  std::fstream    ffp;		// Input/Output binary file stream
  std::string     ffname;		// Bruker FID/Serial file name
  bool	     fbigend;		// Flag for big- vs. little-endian (this arch)
  bool       fbyteordin;	// Flag for big- vs. little-endian input
  int        ftotpts;		// Total points per FID (real + imag)    
  int        fpadding;		// Bytes needed for padding to 256 pt boundary
  int        ffsize;		// Our total file size (bytes)
  row_vector fdata;		// The fid data points
 
// ----------------------------------------------------------------------------
// ----------------------------- PRIVATE FUNCTIONS ----------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                      XWinNMR Fid File Error Handling
// ____________________________________________________________________________

/* These functions take care of any errors encountered when reading, writing,   
   and setting parameters in Bruker acquisition parameter files.
 
	Input		FidPar  : XWinNMR acqusition parameters (this)
			eidx    : Error index
                        pname   : Additional error message
                        noret   : Flag for linefeed (0=linefeed)
        Output          void    : An error message is output                 */

         void XWinFiderror(int    eidx, int noret=0) const;
         void XWinFiderror(int    eidx, const std::string& pname, int noret=0) const;
volatile void XWinFidfatality(int eidx) const;
volatile void XWinFidfatality(int eidx, const std::string& pname) const;


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
 
/* Regardless of the FID size, a Bruker serial file (ser) enforces the rule
   that each FID must begin on a 1024 byte block boundary (256 points).  So,
   when writing serial files we will pad the output file with zeros until the
   boundary is reached. Similarly, we will check that all FIDs are all written
   written into the file beginning at such a boundary.                       */  
 
void SetPadding();
bool CheckBoundary();
void SkipPadding();
void AddPadding();

public:
// ____________________________________________________________________________ 
// A                  XWinFid File Constructors, Destructor
// ____________________________________________________________________________

/* These are the constructors of the class handling Bruker XWinNMR 1D
   acquisition binary files.  The default name such files is "fid". The
   output byte order is automatically set for the computer architecture.

        //                      name    : External filename 
	//			vx	: Data vector
	//			TD	: Total points (real + imag)
	//			byteord : Input byte order                   */
 
        XWinFid();
        XWinFid(const std::string& name, const row_vector& vx);
        XWinFid(const std::string& name, int TD, bool byteord);
        XWinFid(const XWinFid& XWF);
virtual ~XWinFid();
void    operator= (const XWinFid& XWF);

// ____________________________________________________________________________
// B                   XWinNMR Fid File Access Functions
// ____________________________________________________________________________
 
/* These functions allow users to get some simple intformation regarding the
   contents of the Bruker data acquisition file.                             */
 
        std::string     fidname()   const;		// File name
virtual int        size()      const;		// No. complex points
        int        TDF()       const;		// No. total points
virtual bool       order()     const;		// Data byte order
virtual int        bytes()     const;		// File size in bytes
virtual int        blocks()    const;		// No. FIDs
virtual int        pad()       const;		// No. padding bytes
virtual row_vector data()      const;		// Data points

// ____________________________________________________________________________
// C                       XWinFid Input Functions
// ____________________________________________________________________________
 
/* This function will read in a binary fid file from a Bruker 1D acqusition
   in XWinNMR.  Typically this is just done during construction of the class.
   However, it is convenient to have an empty class constructor which can be
   used in conjunction with read/write to perform the same feat.  The Bruker
   file "fid" contains values that are read in as re,im,re,im,...., each point
   stored as a 32-bit integer.  It is very important that we know the number 
   of points in an FID and the byte order used in storing the data.  This
   information must be supplied, either when calling the read function or when
   constructing objects of this class type.  Such information is usually stored
   in corresponding ASCII parameter file, typically named acqu and acqus. Such
   details should be attended to in higher classes utilizing these functions.*/ 
 
virtual bool read(const std::string& name, bool bord, 
                                                         int TD=0, int warn=2);

// ____________________________________________________________________________
// D                       XWinFid Output Functions
// ____________________________________________________________________________

/* These functions will write an FID or series of FIDs to a file in Bruker 
   XWinNMR format. This can be done either by directly writing a row vector 
   or matrix into such a file, or by writing a series of row vectors. Note that
   ususally there is a parameter file associated with such data files, called
   acqus in the case of a single FID file or called acqu2s for a serial file  
   file of multiple FIDs.the                                                 */

virtual int write(const std::string& F, const row_vector& data, int wrn=2);

// ____________________________________________________________________________
// E                    XWinFid Standard Output Functions
// ____________________________________________________________________________
 
 /* These functions allow users to have a quick peek at the "fid" contents.
    They do not do any manipulations to the file whatsoever, they will only
    report on what in known about the data values and file sizes.             */ 

virtual std::ostream& print(std::ostream& ostr, int full=0, int hdr=1) const;
friend  std::ostream& operator<< (std::ostream& O, const XWinFid& F);
};

#endif 								// XWinFid.h
