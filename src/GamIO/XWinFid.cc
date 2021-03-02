/* XWinFid.cc ***************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A                           **
**                                                                      **
**      XWinFid						Implementation	**
**                                                                      **
**      Copyright (c) 1999                                              **
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
**			  \___ pdata -- 1 -- proc, procs, meta		**
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

#ifndef _XWinFid_cc_			// Is this file already included?
#  define _XWinFid_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler
#    pragma implementation
#  endif

#include <GamIO/XWinFid.h>		// Include interface
#include <GamIO/BinIOBase.h>		// Include binary I/O functions
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <Matrix/row_vector.h>		// Include GAMMA row vectors
#include <stdlib.h>

using std::string;			// Using libstdc++ strings
using std::ostream;			// Using libstdc++ output streams
using std::ios;				// Using libstdc++ file type settings

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                      XWinNMR Fid File Error Handling
// ____________________________________________________________________________

/* These functions take care of any errors encountered when reading, writing,
   and setting parameters in Bruker acquisition parameter files.

        Input           FidPar  : XWinNMR acqusition parameters (this)
                        eidx    : Error index
                        pname   : Additional error message
                        noret   : Flag for linefeed (0=linefeed)
        Output          void    : An error message is output                 */  


void XWinFid::XWinFiderror(int eidx, int noret) const
  {
  string hdr("Bruker Fid/Serial File");
  switch(eidx)
    {
    case 20: GAMMAerror(hdr,"ASCII Mode Disallowed",     noret);break;	//(20)
    case 21: GAMMAerror(hdr,"The Type MUST Be Binary",   noret);break;	//(21)
    case 22: GAMMAerror(hdr,"Open For Input & Output!",  noret);break;	//(22)
    case 23: GAMMAerror(hdr,"Do ios::in Or ios::iout",   noret);break;	//(23)
    case 24: GAMMAerror(hdr,"Open At End Disallowed",    noret);break;	//(24)
    case 25: GAMMAerror(hdr,"Do Not Open With ios::ate", noret);break;	//(25)
    case 26: GAMMAerror(hdr,"Open For Append Disallowed",noret);break;	//(26)
    case 27: GAMMAerror(hdr,"Do Not Open With ios::app", noret);break;	//(27)
    case 28: GAMMAerror(hdr,"Attached File Stream Bad",  noret);break;	//(28)
    case 29: GAMMAerror(hdr,"Cannot Perform Any I/O",    noret);break;	//(29)
    case 30: GAMMAerror(hdr,"Can't Set Block Boundary?", noret);break;	//(30)
    case 31: GAMMAerror(hdr,"Can't Check Block Boundary",noret);break;	//(31)
    case 40: GAMMAerror(hdr,"Read Start Off Boundary",   noret);break;	//(40)
    case 41: GAMMAerror(hdr,"Problems Reading FID",      noret);break;	//(41)
    case 42: GAMMAerror(hdr,"Write Start Off Boundary",  noret);break;	//(42)
    case 43: GAMMAerror(hdr,"Problems Writing FID",      noret);break;	//(43)
    case 44: GAMMAerror(hdr,"Problems With File Size",   noret);break;	//(44)
    case 45: GAMMAerror(hdr,"1 FID Bigger Than File!",   noret);break;	//(45)
    case 46: GAMMAerror(hdr,"Block Size (TD) Too Big?",  noret);break;	//(46)
    default: GAMMAerror(hdr,eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }
     

void XWinFid::XWinFiderror(int eidx, const string& pname, int noret) const
  {
  string hdr("Bruker Fid/Serial File");
  switch(eidx)
    {
    case 1:  GAMMAerror(hdr,   1, pname, noret); break; // File Problems   (1)
    case 10: GAMMAerror(hdr,"Can't Read Array "+pname, noret); break;    //(10)
    default: GAMMAerror(hdr,  -1, pname, noret); break; // Unknown Error   (-1)
    }
  }
     
volatile void XWinFid::XWinFidfatality(int eidx) const
  {
  XWinFiderror(eidx, 1);
  if(eidx) XWinFiderror(0);
  GAMMAfatal();					// Clean exit from program
  }
     
volatile void XWinFid::XWinFidfatality(int eidx, const string& pname) const
  {
  XWinFiderror(eidx, pname, 1);
  if(eidx) XWinFiderror(0);
  GAMMAfatal();					// Clean exit from program
  }

// ----------------------------------------------------------------------------
// ii                      Mode & Size Checking Functions
// ----------------------------------------------------------------------------

/* These functions check the mode specified when opening a XWinFid MAT file.
   We do not allow all of the possible modes availble to C++ filestreams.  In
   particular, we do not allow users to append to an existing file (since this
   would be bad if the existing data was written with a different byte order).
   Additonally, there seems to be some discrepancy between compilers as to
   what modes are available to fstream, so we just stick to the basics.

	in  = 1 (0x01)	out       = 2  (0x02)	ate      = 4   (0x03)
	app = 8	(0x08)	trunc     = 16 (0x10)	nocreate = 32  (0x20)
			noreplace = 64 (0x40) 	binary   = 128 (0x80)
   
   Use of noreplace & nocreate maybe didn't make it to ANSI standard...      */

void XWinFid::CheckMode(int mode)
  {
  if(!(mode & ios::binary))
    {
    XWinFiderror(20, 1);			// File can't be ASCII
    XWinFidfatality(21);			// File must be binary
    }
  if((mode & ios::out) && (mode & ios::in))
    {
    XWinFiderror(22, 1);			// File can't be input & output
    XWinFidfatality(23);			// Open file as either input or output
    }
  if(mode & ios::ate)
    {
    XWinFiderror(24, 1);			// File can't open at end
    XWinFidfatality(25);			// Must be opened without ios::ate
    }
  if(mode & ios::app)
    {
    XWinFiderror(26, 1);			// File can't open in append mode
    XWinFidfatality(27);			// Must be opened with ios:app
    }
  }

int XWinFid::CheckSize(int warn)
  {
  int fidbytes = 4*ftotpts;		// Number of bytes per FID
  if(fidbytes > ffsize)			// Insure 1 FID is NOT bigger than 
    {					// the entire file size
    if(warn)
      {
      XWinFiderror(44, 1);		// File size problems
      if(warn > 1) XWinFidfatality(45); 	// Block size too big
      else         XWinFiderror(45, 1);
      }
    return 0;
    }
  if(ffsize > fidbytes)			// For a serial file, the number of
    {					// bytes in the file must exactly
    fidbytes +=  fpadding; 		// match a multiple of padded FIDs
    int delb = ffsize;			// This is bytes in file
    while(delb > 0) delb -= fidbytes;	// Subtract off buffered Fids
    if(delb < 0)			// Should now be exactly 0!
      {
      if(warn)
        {
        XWinFiderror(44, 1);		// File size problems
        if(warn > 1) XWinFidfatality(45);// Block size too big
        else         XWinFiderror(45, 1);
        }
      return 0;
      }
    }
  return 1;
  }
 
// ----------------------------------------------------------------------------
// iii                 File Padding & Boundary Functions
// ----------------------------------------------------------------------------
 
/* Regardless of the FID size, a Bruker serial file (ser) enforces the rule     
   that each FID must begin on a 1024 byte block boundary (256 points).  So,
   when writing serial files we will pad the output file with zeros until the
   boundary is reached. Similarly, we will check that all FIDs are all written
   written into the file beginning at such a boundary.                       */

void XWinFid::SetPadding()
  {
  fpadding = 0;				// Assume there is no padding
  if(ftotpts == 0) return;		// No padding if no points
  fpadding = ftotpts; 			// Begin will points per FID 
  while(fpadding >= 256) fpadding-=256;	// Clip chunks of 256 points
  if(fpadding < 0) XWinFidfatality(30);	// Quit if can't find boundary
  if(fpadding > 0)			// If padding needed, set it to number
    fpadding = 4*(256-fpadding);	// of bytes to reach 256 point boundary
  }

bool XWinFid::CheckBoundary()
  {
  int pos = ffp.tellp();		// Get current file position
  if(!pos) return 1;			// Start of file is boundary
  if(!ftotpts) XWinFidfatality(31);	// Quit if can't check boundary
  int blkbytes = 4*ftotpts+fpadding;	// Total bytes per "padded" FID
  while(pos >= blkbytes) pos-=blkbytes;	// Clip off blocks
  if(!pos) return 1;			// This is zero if a boundary
  if(pos < 0) XWinFidfatality(31);	// Problems if negative
  return 0; 				// Here if NOT a boundary
  }

void XWinFid::SkipPadding() { ffp.seekp(fpadding, ios::cur); }
  
void XWinFid::AddPadding()
  {
  char pc = ' ';			// Padding charactor
  for(int i=0; i<fpadding; i++) 	// Write blanks, assuming there
    ffp.write(&pc, sizeof(char));	// is only 1 byte per character
  }
 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                    XWinFid CONSTRUCTORS, DESTRUCTOR
// ____________________________________________________________________________

/* Note that size this class should always match what is on the disk, there
   are no access functions for directly setting class elements.  Thus, the
   convenient empty constructor below will either need an assigment to fill
   it or use of a read/write function.                                       */

XWinFid::XWinFid()
  { 
  ffname     = string("fid");		// Set a default file name
  fbigend    = WeRBigEnd();		// Set output byte order (our arch)
  fbyteordin = 0;			// Set default input byte order
  ftotpts    = 0;			// Set no FID point size
  fpadding   = 0;			// Set no block padding
  ffsize     = 0;			// Set file size to 0
  }

XWinFid::XWinFid(const string& name, const row_vector& data)
        :ffp(name.c_str(), ios::binary|ios::out)
  {
  if(!ffp.good())			// Insure Bruker File is O.K.
    {
    XWinFiderror(28,1);			//   Problems with filestream
    XWinFiderror(1,name,1);		//   Problems with file
    XWinFidfatality(29);		//   Can't write a damn thing
    }
  ffname     = name;			// Set the Bruker filename
  fbigend    = WeRBigEnd();		// Set output byte order (our arch)
  fbyteordin = fbigend;			// Set input byte order (our arch)
  fdata      = data;			// Set our data points
  ftotpts    = 2*data.size();		// Set points per fid (re + im)
  SetPadding();				// Set padding for data boundary
  ffp.seekp(0);				// Set things at file start
  int npts   = fdata.size();		// Number of points (2*TD)
  int32_t rval, ival;			// For point reals,imags
  for(int i=0; i<npts; i++)		// Loop over all points
    {
    rval = int32_t(data.getRe(i));		//   Real point as int
    ival = int32_t(data.getIm(i));		//   Imag point as int
    ffp.write((char*)&rval, sizeof(int32_t));	//   Write real point
    ffp.write((char*)&ival, sizeof(int32_t));	//   Write imag point
    }
  AddPadding();				// Pad to end of block
  ffp.seekp(0,ios::end);		// Go to the end of the file
  ffsize = ffp.tellp();			// Get total bytes in the file
  ffp.close();				// Close the file
  }

XWinFid::XWinFid(const string& name, int TD, bool byteord)
        :ffp(name.c_str(), ios::binary|ios::in)
  {
  if(!ffp.good())			// Insure Bruker File is O.K.
    {
    XWinFiderror(28,1);			//   Problems with filestream
    XWinFiderror(1,name,1);		//   Problems with file
    XWinFidfatality(29);		//   Can't read a damn thing
    }
  ffname     = name;			// Set the Bruker filename
  fbigend    = WeRBigEnd();		// Set output byte order (our arch)
  fbyteordin = byteord;			// Set input byte order (Bruker)
  ftotpts    = TD;			// Set points per fid (re + im)
  SetPadding();				// Set any padding for data boundary
  ffp.seekp(0,ios::end);		// Go to the end of the file
  ffsize = ffp.tellp();			// See how many bytes are in the file
  ffp.seekp(0);				// Set things back to the start
  if(!CheckSize(1))			// Insure file size O.K.
    XWinFidfatality(46);		// for writing, TD too large else
  int swapon = 0;			// Assume no byte swapping
  if(fbigend != fbyteordin) swapon = 1;	// Need to swap if mismatch
  fdata = row_vector(ftotpts/2);	// Array for fid data
  int32_t ptre, ptim;			// These will be input values
  for(int i=0; i<ftotpts/2; i++)	// Loop over fid points
    {  
    ffp.read((char*)&ptre,sizeof(int32_t));//   Read a (real) point
    ffp.read((char*)&ptim,sizeof(int32_t));//   Read a (imag) point
    if(swapon) {Swap(ptre); Swap(ptim);}//   Byte swap if needed
    fdata.put(complex(ptre,ptim), i);	//   Store point
    }  
  ffp.close();				// Close the file
  }

XWinFid::XWinFid(const XWinFid& XWF)
  
        // Input              	this 	: XWinNMR Binary FID file
        // 			XWF	: XWinNMR Binary FID file
	// Output		none	: this set identical to XWF

  {
//ffp        = XWF.ffp;           	// Always opened/closed anew
  ffname     = XWF.ffname;		// Copy the file name
  fbigend    = XWF.fbigend;		// Copy our endian type
  fbyteordin = XWF.fbyteordin;		// Copy file endian type
  ftotpts    = XWF.ftotpts;		// Copy total bytes (re+im)
  fpadding   = XWF.fpadding;		// Copy needed padding
  ffsize     = XWF.ffsize;		// Copy total file size
  fdata      = XWF.fdata;		// Copy data vector
  }

XWinFid::~XWinFid() {}
  
        // Input              	XWF	: XWinNMR Binary FID file
	// Output		none	: XWF is destructed
	// Note			  	: Of course, this doesn't touch
	//				  the external file at all

void XWinFid::operator= (const XWinFid& XWF)
  
        // Input              	this 	: XWinNMR Binary FID file
        // 			XWF	: XWinNMR Binary FID file
	// Output		none	: XWF is copied to this

  {
//ffp        = XWF.ffp;           	// Always opened/closed anew
  ffname     = XWF.ffname;		// Copy the file name
  fbigend    = XWF.fbigend;		// Copy our endian type
  fbyteordin = XWF.fbyteordin;		// Copy file endian type
  ftotpts    = XWF.ftotpts;		// Copy total bytes (re+im)
  fpadding   = XWF.fpadding;		// Copy needed padding
  ffsize     = XWF.ffsize;		// Copy total file size
  fdata      = XWF.fdata;		// Copy data vector
  }

// ____________________________________________________________________________
// B                   XWinNMR Fid File Access Functions
// ____________________________________________________________________________

/* These functions allow users to get some simple information regarding the
   contents of the Bruker binary data acqusition file, fid.                  */

string     XWinFid::fidname()  const { return ffname; }
int        XWinFid::size()     const { return ftotpts/2; }
int        XWinFid::TDF()      const { return ftotpts; }
int        XWinFid::bytes()    const { return ffsize; }
int        XWinFid::pad()      const { return fpadding; }
bool       XWinFid::order()    const { return fbyteordin; }
row_vector XWinFid::data()     const { return fdata; }

int XWinFid::blocks() const
  {
  int npts = ffsize/4;				// Number 32-bit ints in file
  int ptsperfid = ftotpts + fpadding/4;		// Points used by each FID
  if(!ptsperfid) return 0;			// Insure no size overflow
  return npts/ptsperfid;			// Number of FIDs in the file
  }

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
   constructing objects of this class type.  Such information is usually stroed
   in corresponding ASCII parameter file, typically named acqu and acqus. Such
   details should be attended to in higher classes utilizing these functions.*/ 
 
bool XWinFid::read(const string& name, bool byteord, int TD, int warn)
  {
  ffp.open(name.c_str(), ios::binary|ios::in);
  if(!ffp.good())			// Insure Bruker File is O.K.
    {
    if(warn)
      {
      XWinFiderror(28,1);		//   Problems with filestream
      XWinFiderror(1,name,1);		//   Problems with file
      if(warn>1) XWinFidfatality(29);	//   Can't read a damn thing
      else       XWinFiderror(1,name);
      }
    return false;
    }
  ffname     = name;			// Set the Bruker filename
  fbigend    = WeRBigEnd();		// Set output byte order (our arch)
  fbyteordin = byteord;			// Set input byte order (Bruker)
  ffp.seekp(0,ios::end);		// Go to the end of the file
  ffsize = ffp.tellp();			// See how many bytes are in the file
  ffp.seekp(0);				// Set things back to the start
  if(TD>0) ftotpts= TD;			// Set points per fid (re + im)
  else     ftotpts=ffsize/sizeof(int32_t);	// Set points per fid (re + im)
  SetPadding();				// Set any padding for data boundary
  if(!CheckSize(1))			// Insure file size O.K. if open
    XWinFidfatality(46);		// for writing, TD too large else
  int swapon = 0;			// Assume no byte swapping
  if(fbigend != fbyteordin) swapon = 1;	// Need to swap if mismatch
  fdata = row_vector(ftotpts/2);	// Array for fid data
  int32_t ptre, ptim;			// These will be input values
  for(int i=0; i<ftotpts/2; i++)	// Loop over fid points
    {  
    ffp.read((char*)&ptre,sizeof(int32_t));//   Read a (real) point
    ffp.read((char*)&ptim,sizeof(int32_t));//   Read a (imag) point
    if(swapon) {Swap(ptre); Swap(ptim);}//   Byte swap if needed
    fdata.put(complex(ptre,ptim),i);	//   Store point
    }  
  ffp.close();				// Close the file
  return true;				// Return successfull
  }
 
 
// ____________________________________________________________________________
// D                       XWinFid Output Functions
// ____________________________________________________________________________
 
/* This function will write a binary fid file for a Bruker 1D acqusition
   in XWinNMR.  Typically this is just done during construction of the class.
   However, it is convenient to have an empty class constructor which can be
   used in conjunction with read/write to perform the same feat.  The Bruker
   file "fid" contains values that are read in as re,im,re,im,...., each point
   stored as a 32-bit integer. The fid file should be consistent with the
   corresponding ASCII parameter files, acqu and acqus. However, such details
   should be attended to in higher classes utilizing these functions.        */ 
 
int XWinFid::write(const string& name, const row_vector& vx, int warn)
  {
  ffp.open(name.c_str(), ios::binary|ios::out);
  if(!ffp.good())			// Insure Bruker File is O.K.
    {
    if(warn)
      {
      XWinFiderror(28,1);		//   Problems with filestream
      XWinFiderror(1,name,1);		//   Problems with file
      if(warn>2) XWinFidfatality(29);	//   Can't write a damn thing
      else       XWinFiderror(29);
      }
    return 0;
    }
  ffname     = name;			// Set the Bruker filename
  fbigend    = WeRBigEnd();		// Set output byte order (our arch)
  fbyteordin = fbigend;			// Set input byte order (our arch)
  fdata      = vx;			// Set our data points
  ftotpts    = 2*vx.size();		// Set points per fid (re + im)
  SetPadding();				// Set padding for data boundary
  ffp.seekp(0);				// Set things at file start
  int npts   = fdata.size();		// Number of points (2*TD)
  int32_t rval, ival;			// For point reals,imags
  for(int i=0; i<npts; i++)		// Loop over all points
    {
    rval = int32_t(vx.getRe(i));		//   Real point as int
    ival = int32_t(vx.getIm(i));		//   Imag point as int
    ffp.write((char*)&rval, sizeof(int32_t));	//   Write real point
    ffp.write((char*)&ival, sizeof(int32_t));	//   Write imag point
    }
  AddPadding();				// Pad to end of block
  ffp.seekp(0,ios::end);		// Go to the end of the file
  ffsize = ffp.tellp();			// Get total bytes in the file
  ffp.close();				// Close the file
  return 1;
  }
 
// ____________________________________________________________________________
// E                    XWinFid Standard Output Functions
// ____________________________________________________________________________
 
/* These functions allow users to have a quick peek at the "fid" contents.
   They do not do any manipulations to the file whatsoever, they will only
   report on what in known about the data values and file sizes.             */ 

ostream& XWinFid::print(ostream& ostr, int full, int hdr) const
  {
  if(hdr)
    ostr << "\n\n\t\tXWinNMR Binary FID File " << ffname << "\n";
  ostr << "\n\t\tFile Name:              " << ffname;
  ostr << "\n\t\tTotal Points (Re+Im):   " << ftotpts;
  ostr << "\n\t\tTotal Bytes:            " << ffsize;
  if(fpadding)
    ostr << "\n\t\tAdded Padding:          " << fpadding;
  ostr << "\n\t\tFile Byte Ordering:     ";
  if(!fbyteordin) ostr << "Little Endian";
  else            ostr << "Big Endian";
  double imax=-1.e-15, imin=1.e15;
  double rmax=-1.e-15, rmin=1.e15;
  double nmax=0,       nmin=1.e15;
  double n, vi, vr;
  for(int i=0; i<fdata.size(); i++)
    {
    n  = norm(fdata.get(i));
    vi = fdata.getIm(i);
    vr = fdata.getRe(i);
    if(n  > nmax) nmax = n;
    if(n  < nmin) nmin = n;
    if(vi > imax) imax = vi;
    if(vi < imin) imin = vi;
    if(vr < rmin) rmin = vr;
    if(vr > rmax) rmax = vr;
    }
  ostr << "\n\t\tData Norm Max/Min:      "
       << nmax << "/" << nmin;
  ostr << "\n\t\tData Real Max/Min:      "
       << rmax << "/" << rmin;
  ostr << "\n\t\tData Imaginary Max/Min: "
       << imax << "/" << imin;
  if(full)
    {
    ostr << "\n\t\tData Points:            ";
    for(int k=0; k<abs(full); k++)
      ostr << "\n\t\t\t" << k << ". " << fdata.get(k);
    }
  return ostr;
  }

ostream& operator<< (ostream& O, const XWinFid& F) { F.print(O); return O; }

#endif						// XWinFid.cc
