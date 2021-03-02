/* XWinSpec.h ****************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      XWinSpec                                    Interface		**
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
** sets. This class embodies two Bruker binary data files, 1r & 1i,	**
** associated with a spectrum generated during processing of a 1D	**
** acquisition. There is not anything fancy about such files, they	**
** just contain a series of 32-bit integers. Two important aspects are	**
** 1.) How many points they contain , and 2.) The data byte order. That	**
** information is normally contained in an associated parameter file,	**
** procs. Below	is a typical Bruker directory structure associated with	**
** a 1D	acquisiton that has been processed.				**
**                                                                      **
**                          __ acqu  (changable parameter file)		**
**			   / 						**
**                        /___ acqus (static parameter file)		**
**			 /						**
**  expname -- expnum --< ---- fid (binary data)			**
**			 \						**
**			  \___ pdata -- 1 -- proc, procs, meta, 	**
**                                           1r,   1i,    outd		**
**									**
** This class would then handle the files "1r" & "1i".  It will NOT	**
** deal with the directory hierarchy shown above nor with any of the 	**
** associated ASCII parameter files.  That is left up to higher level 	**
** classes supporting XWinNMR. 						**
**                                                                      **
** Note that this is a DISK entity, not a row_vector. That is, when	**
** this is instructed to write a spectrum it must be given a row_vector	**
** & some basic info as to how it should be written.  Similarly, it can	**
** read the two files into a row_vector if needed.			**
**                                                                      **
*************************************************************************/

#ifndef   XWinSpec_H_ 				// Is file already included?
#  define XWinSpec_H_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// This is the implementation
#  endif
 
#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <fstream>                              // File streams from libstdc++
#include <string>				// Include libstdc++ strings
#include <Matrix/row_vector.h>			// Know about row vectors
#include <Matrix/matrix.h>			// Know about matrices

class XWinSpec

// ----------------------------------------------------------------------------
// ------------------------------- STRUCTURE ----------------------------------
// ----------------------------------------------------------------------------
 
  {
  std::fstream    sfpre;	// Input/Output binary file stream (reals)
  std::fstream    sfpim;	// Input/Output binary file stream (imags)
  std::string     sfname;	// Bruker Spectrum base file name
  bool	     sbigend;		// Flag for big- vs. little-endian (this arch)
  bool       sbyteordin;	// Flag for big- vs. little-endian input
  int        stotpts;		// Total points per spectrum (real + imag)    
  int        spadding;		// Bytes needed for padding to 256 pt boundary
  int        sfsize;		// Our total file size (bytes)
  row_vector sdata;		// The spectrum data points
 
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

         void XWSerror(int eidx, int noret=0) const;
         void XWSerror(int eidx, const std::string& pname, int noret=0) const;
volatile void XWSfatal(int eidx) const;
volatile void XWSfatal(int eidx, const std::string& pname) const;


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
// A                  XWinSpec File Constructors, Destructor
// ____________________________________________________________________________

/* These are the constructors of the class handling Bruker XWinNMR 1D
   acquisition binary files.  The default name such files is "fid". The
   output byte order is automatically set for the computer architecture.

        //                      name    : External filename 
	//			vx	: Data vector
	//			SI	: Total points (real + imag)
	//			byteord : Input byte order                   */
 
        XWinSpec();
        XWinSpec(const std::string& name, const row_vector& vx);
        XWinSpec(const std::string& name, int SI, bool byteord);
        XWinSpec(const XWinSpec& XWF);
virtual ~XWinSpec();
void    operator= (const XWinSpec& XWF);

// ____________________________________________________________________________
// B                   XWinNMR Fid File Access Functions
// ____________________________________________________________________________
 
/* These functions allow users to get some simple intformation regarding the
   contents of the Bruker data acquisition file.                             */
 
        std::string     specname()  const;	// Base file name
        std::string     specrname() const;	// Reals file name
        std::string     speciname() const;	// Imagins. file name
virtual int             size()      const;	// No. complex points
        int             SSI()       const;	// No. total points
virtual bool            order()     const;	// Data byte order
virtual int             bytes()     const;	// File size in bytes
virtual int             blocks()    const;	// No. FIDs
virtual int             pad()       const;	// No. padding bytes
virtual row_vector      data()      const;	// Data points

// ____________________________________________________________________________
// C                       XWinSpec Input Functions
// ____________________________________________________________________________
 
/* This function will read in a binary spectrum from a Bruker 1D data set
   in XWinNMR.  Typically this is just done during construction of the class.
   However, it is convenient to have an empty class constructor which can be
   used in conjunction with read/write to perform the same feat.  The Bruker
   files "1r" and "1i" contain real and imaginary values respectively. Each point
   stored as a 32-bit integer.  It is very important that we know the number 
   of points in an FID and the byte order used in storing the data.  This
   information must be supplied, either when calling the read function or when
   constructing objects of this class type.  Such information is usually stored
   in corresponding ASCII parameter file, typically named acqu and acqus. Such
   details should be attended to in higher classes utilizing these functions.*/ 
 
virtual bool read(const std::string& name, bool bord, 
                                                          int SI=0, int wrn=2);

// ____________________________________________________________________________
// D                       XWinSpec Output Functions
// ____________________________________________________________________________

/* These functions will write a binary spectrum to (two) files in Bruker 
   XWinNMR format. This is done by directly writing a row vector to the files.
   Note that ususally there is a parameter file associated with such files,
   called procss in the case of a single spectrum.                           */

virtual int write(const std::string& F, const row_vector& dat, int wrn=2);

// ____________________________________________________________________________
// E                    XWinSpec Standard Output Functions
// ____________________________________________________________________________
 
 /* These functions allow users to have a quick peek at the "fid" contents.
    They do not do any manipulations to the file whatsoever, they will only
    report on what in known about the data values and file sizes.             */ 

virtual std::ostream& print(std::ostream& ostr, int full=0, int hdr=1) const;
friend  std::ostream& operator<< (std::ostream& O, const XWinSpec& F);
};

#endif 								// XWinSpec.h
