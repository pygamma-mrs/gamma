/* XWinSpec.cc ***************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A                           **
**                                                                      **
**      XWinSpec					Implementation	**
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
** sets. This class embodies two Bruker binary data files, 1r & 1i,     **
** associated with a spectrum generated during processing of a 1D       **
** acquisition. There is not anything fancy about such files, they      **
** just contain a series of 32-bit integers. Two important aspects are  **
** 1.) How many points they contain , and 2.) The data byte order. That **
** information is normally contained in an associated parameter file,   **
** procs. Below is a typical Bruker directory structure associated with **
** a 1D acquisiton that has been processed.                             **
**                                                                      **
**                          __ acqu  (changable parameter file)         **
**                         /                                            **
**                        /___ acqus (static parameter file)            **
**                       /                                              **
**  expname -- expnum --< ---- fid (binary data)                        **
**                       \                                              **
**                        \___ pdata -- 1 -- proc, procs, meta,         **
**                                           1r,   1i,    outd          **
**                                                                      **
** This class would then handle the files "1r" & "1i".  It will NOT     **
** deal with the directory hierarchy shown above nor with any of the    **
** associated ASCII parameter files.  That is left up to higher level   **
** classes supporting XWinNMR.                                          **
**                                                                      **
** Note that this is a DISK entity, not a row_vector. That is, when     **
** this is instructed to write a spectrum it must be given a row_vector **
** & some basic info as to how it should be written.  Similarly, it can **
** read the two files into a row_vector if needed.                      **
**                                                                      **
*************************************************************************/

#ifndef   XWinSpec_cc_			// Is this file already included?
#  define XWinSpec_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler
#    pragma implementation
#  endif

#include <GamIO/XWinSpec.h>		// Include interface
#include <GamIO/BinIOBase.h>		// Include binary I/O functions
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <Matrix/row_vector.h>		// Include GAMMA row vectors
#include <stdlib.h>

using std::string;			// Using libstdc++ strings
using std::ios;
using std::ostream;			// Using libstdc++ output streams

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                      XWinNMR Spectrum File Error Handling
// ____________________________________________________________________________

/* These functions take care of any errors encountered when reading, writing,
   and setting parameters in Bruker acquisition parameter files.

        Input           FidPar  : XWinNMR acqusition parameters (this)
                        eidx    : Error index
                        pname   : Additional error message
                        noret   : Flag for linefeed (0=linefeed)
        Output          void    : An error message is output                 */  


void XWinSpec::XWSerror(int eidx, int noret) const
  {
  string hdr("Bruker 1D Spectrum Files");
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
    case 32: GAMMAerror(hdr,"Real.-Imag. Sizes Differ?", noret);break;	//(32)
    case 40: GAMMAerror(hdr,"Read Start Off Boundary",   noret);break;	//(40)
    case 41: GAMMAerror(hdr,"Problems Reading Spectrum",      noret);break;	//(41)
    case 42: GAMMAerror(hdr,"Write Start Off Boundary",  noret);break;	//(42)
    case 43: GAMMAerror(hdr,"Problems Writing Spectrum",      noret);break;	//(43)
    case 44: GAMMAerror(hdr,"Problems With File Size",   noret);break;	//(44)
    case 45: GAMMAerror(hdr,"1 Spectrum Bigger Than File!",   noret);break;	//(45)
    case 46: GAMMAerror(hdr,"Block Size (SI) Too Big?",  noret);break;	//(46)
    default: GAMMAerror(hdr,eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }
     

void XWinSpec::XWSerror(int eidx, const string& pname, int noret) const
  {
  string hdr("Bruker 1D Spectrum Files");
  switch(eidx)
    {
    case 1:  GAMMAerror(hdr,   1, pname, noret); break; // File Problems   (1)
    case 10: GAMMAerror(hdr,"Can't Read Array "+pname, noret); break;    //(10)
    default: GAMMAerror(hdr,  -1, pname, noret); break; // Unknown Error   (-1)
    }
  }
     
volatile void XWinSpec::XWSfatal(int eidx) const
  {
  XWSerror(eidx, 1);
  if(eidx) XWSerror(0);
  GAMMAfatal();					// Clean exit from program
  }
     
volatile void XWinSpec::XWSfatal(int eidx, const string& pname) const
  {
  XWSerror(eidx, pname, 1);
  if(eidx) XWSerror(0);
  GAMMAfatal();					// Clean exit from program
  }

// ----------------------------------------------------------------------------
// ii                      Mode & Size Checking Functions
// ----------------------------------------------------------------------------

/* These functions check mode specified when opening a XWinSpec binary file.
   We do not allow all of the possible modes availble to C++ filestreams.  In
   particular, we do not allow users to append to an existing file (since this
   would be bad if the existing data was written with a different byte order).
   Additonally, there seems to be some discrepancy between compilers as to
   what modes are available to fstream, so we just stick to the basics.

	in  = 1 (0x01)	out       = 2  (0x02)	ate      = 4   (0x03)
	app = 8	(0x08)	trunc     = 16 (0x10)	nocreate = 32  (0x20)
			noreplace = 64 (0x40) 	binary   = 128 (0x80)
   
   Use of noreplace & nocreate maybe didn't make it to ANSI standard...      */

void XWinSpec::CheckMode(int mode)
  {
  if(!(mode & ios::binary))
    {
    XWSerror(20, 1);			// File can't be ASCII
    XWSfatal(21);			// File must be binary
    }
  if((mode & ios::out) && (mode & ios::in))
    {
    XWSerror(22, 1);			// File can't be input & output
    XWSfatal(23);			// Open file as either input or output
    }
  if(mode & ios::ate)
    {
    XWSerror(24, 1);			// File can't open at end
    XWSfatal(25);			// Must be opened without ios::ate
    }
  if(mode & ios::app)
    {
    XWSerror(26, 1);			// File can't open in append mode
    XWSfatal(27);			// Must be opened with ios:app
    }
  }

int XWinSpec::CheckSize(int warn)
  {
  int spec = 4*stotpts;			// Bytes per Spectrum (Re or Im)
  if(spec > sfsize)			// Insure 1 Spectrum is NOT bigger than 
    {					// the entire file size
    if(warn)
      {
      XWSerror(44, 1);			// File size problems
      if(warn > 1) XWSfatal(45); 	// Block size too big
      else         XWSerror(45, 1);
      }
    return 0;
    }
  if(sfsize > spec)			// For a serial file, the number of
    {					// bytes in the file must exactly
    spec +=  spadding; 			// match a multiple of padded Spectrums
    int delb = sfsize;			// This is bytes in file
    while(delb > 0) delb -= spec;	// Subtract off buffered Fids
    if(delb < 0)			// Should now be exactly 0!
      {
      if(warn)
        {
        XWSerror(44, 1);		// File size problems
        if(warn > 1) XWSfatal(45);// Block size too big
        else         XWSerror(45, 1);
        }
      return 0;
      }
    }
  return 1;
  }
 
// ----------------------------------------------------------------------------
// iii                 File Padding & Boundary Functions
// ----------------------------------------------------------------------------
 
/* Regardless of the spectrum size, a Bruker binary file (1r & 1i) enforces the 
   rule that each file must begin on a 1024 byte block boundary (256 points).  
   So, when writing these data files we will pad the output file with zeros
   until the boundary is reached. Similarly, we will check that any spectra 
   written into the file beginning at such a boundary.                       */

void XWinSpec::SetPadding()
  {
  spadding = 0;				// Assume there is no padding
  if(stotpts == 0) return;		// No padding if no points
  spadding = stotpts; 			// Begin will points per Spectrum 
  while(spadding >= 256) spadding-=256;	// Clip chunks of 256 points
  if(spadding < 0) XWSfatal(30);	// Quit if can't find boundary
  if(spadding > 0)			// If padding needed, set it to number
    spadding = 4*(256-spadding);	// of bytes to reach 256 point boundary
  }

bool XWinSpec::CheckBoundary()
  {
  int pos = sfpre.tellp();		// Get current file position
  if(!pos) return 1;			// Start of file is boundary
  if(!stotpts) XWSfatal(31);		// Quit if can't check boundary
  int blkbytes = 4*stotpts+spadding;	// Total bytes per "padded" Spectrum
  while(pos >= blkbytes) pos-=blkbytes;	// Clip off blocks
  if(!pos) return 1;			// This is zero if a boundary
  if(pos < 0) XWSfatal(31);		// Problems if negative
  return 0; 				// Here if NOT a boundary
  }

void XWinSpec::SkipPadding()
  {
  sfpre.seekp(spadding, ios::cur);
  sfpim.seekp(spadding, ios::cur);
  }
  
void XWinSpec::AddPadding()
  {
  char pc = ' ';			// Padding charactor
  for(int i=0; i<spadding; i++) 	// Write blanks, assuming there
    { 					// is only 1 byte per character
    sfpre.write(&pc, sizeof(char));
    sfpim.write(&pc, sizeof(char));
    }
  }
 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                    XWinSpec CONSTRUCTORS, DESTRUCTOR
// ____________________________________________________________________________

/* Note that size this class should always match what is on the disk, there
   are no access functions for directly setting class elements.  Thus, the
   convenient empty constructor below will either need an assigment to fill
   it or use of a read/write function.                                       */

XWinSpec::XWinSpec()
  { 
  sfname     = string("1");		// Set a default base file name
  sbigend    = WeRBigEnd();		// Set output byte order (our arch)
  sbyteordin = 0;			// Set default input byte order
  stotpts    = 0;			// Set no spectrum point size
  spadding   = 0;			// Set no block padding
  sfsize     = 0;			// Set file size to 0
  }

XWinSpec::XWinSpec(const string& basename, const row_vector& data)
        :sfpre((basename+"r").c_str(), ios::binary|ios::out),
         sfpim((basename+"i").c_str(), ios::binary|ios::out)
  {
  if(!sfpre.good())			// Insure Bruker File is O.K.
    {
    XWSerror(28,1);			//   Problems with filestream
    XWSerror(1,basename+"r",1);		//   Problems with file
    XWSfatal(29);			//   Can't write a damn thing
    }
  if(!sfpim.good())			// Insure Bruker File is O.K.
    {
    XWSerror(28,1);			//   Problems with filestream
    XWSerror(1,basename+"i",1);		//   Problems with file
    XWSfatal(29);			//   Can't write a damn thing
    }
  sfname     = basename;		// Set the Bruker base filename
  sbigend    = WeRBigEnd();		// Set output byte order (our arch)
  sbyteordin = sbigend;			// Set input byte order (our arch)
  sdata      = data;			// Set our data points
  stotpts    = data.size();		// Set points per spectrum
  SetPadding();				// Set padding for data boundary
  sfpre.seekp(0);			// Set things at real file start
  sfpim.seekp(0);			// Set things at imag file start
  int npts   = sdata.size();		// Number of points (SI)
  int32_t rval, ival;			// For point reals,imags
  for(int i=0; i<npts; i++)		// Loop over all points
    {
    rval = int32_t(data.getRe(i));		//   Real point as int
    ival = int32_t(data.getIm(i));		//   Imag point as int
    sfpre.write((char*)&rval, sizeof(int32_t));	//   Write real point
    sfpim.write((char*)&ival, sizeof(int32_t));	//   Write imag point
    }
  AddPadding();				// Pad to end of block
  sfpre.seekp(0,ios::end);		// Go to the end of the file
  sfsize = sfpre.tellp();		// Get total bytes in the file
  sfpre.close();			// Close the reals file
  sfpim.close();			// Close the imaginaries file
  }

XWinSpec::XWinSpec(const string& name, int SI, bool byteord)
        :sfpre((name+"r").c_str(), ios::binary|ios::in),
         sfpim((name+"i").c_str(), ios::binary|ios::in)
  {
  if(!sfpre.good())			// Insure Bruker File is O.K.
    {
    XWSerror(28,1);			//   Problems with filestream
    XWSerror(1,name+"r",1);		//   Problems with file
    XWSfatal(29);			//   Can't read a damn thing
    }
  if(!sfpim.good())			// Insure Bruker File is O.K.
    {
    XWSerror(28,1);			//   Problems with filestream
    XWSerror(1,name+"i",1);		//   Problems with file
    XWSfatal(29);			//   Can't read a damn thing
    }
  sfname     = name;			// Set the Bruker base filename
  sbigend    = WeRBigEnd();		// Set output byte order (our arch)
  sbyteordin = byteord;			// Set input byte order (Bruker)
  stotpts    = SI;			// Set points per spectrum (re or im)
  SetPadding();				// Set any padding for data boundary
  sfpre.seekp(0,ios::end);		// Go to the end of the file
  sfsize = sfpre.tellp();		// See how many bytes are in the file
  sfpre.seekp(0);			// Set things back to the start
  if(!CheckSize(1))			// Insure file size O.K.
    XWSfatal(46);			// for writing, SI too large else
  int swapon = 0;			// Assume no byte swapping
  if(sbigend != sbyteordin) swapon = 1;	// Need to swap if mismatch
  sdata = row_vector(stotpts/2);	// Array for spectrum data
  int32_t ptre, ptim;			// These will be input values
  for(int i=0; i<stotpts/2; i++)	// Loop over spectrum points
    {  
    sfpre.read((char*)&ptre,sizeof(int32_t)); 	//   Read a (real) point
    sfpim.read((char*)&ptim,sizeof(int32_t)); 	//   Read a (imag) point
    if(swapon) {Swap(ptre); Swap(ptim);}//   Byte swap if needed
    sdata.put(complex(ptre,ptim), i);	//   Store point
    }  
  sfpre.close();			// Close the reals file
  sfpim.close();			// Close the imags file
  }

XWinSpec::XWinSpec(const XWinSpec& XWF)
  
        // Input              	this 	: XWinNMR Binary Spectrum file
        // 			XWF	: XWinNMR Binary Spectrum file
	// Output		none	: this set identical to XWF

  {
//sfpre      = XWF.sfpre;           	// Always opened/closed anew
//sfpim      = XWF.sfpim;           	// Always opened/closed anew
  sfname     = XWF.sfname;		// Copy the file name
  sbigend    = XWF.sbigend;		// Copy our endian type
  sbyteordin = XWF.sbyteordin;		// Copy file endian type
  stotpts    = XWF.stotpts;		// Copy total bytes (re+im)
  spadding   = XWF.spadding;		// Copy needed padding
  sfsize     = XWF.sfsize;		// Copy total file size
  sdata      = XWF.sdata;		// Copy data vector
  }

XWinSpec::~XWinSpec() {}
  
        // Input              	XWF	: XWinNMR Binary Spectrum file
	// Output		none	: XWF is destructed
	// Note			  	: Of course, this doesn't touch
	//				  the external file at all

void XWinSpec::operator= (const XWinSpec& XWF)
  
        // Input              	this 	: XWinNMR Binary Spectrum file
        // 			XWF	: XWinNMR Binary Spectrum file
	// Output		none	: XWF is copied to this

  {
//sfpre      = XWF.sfpre;           	// Always opened/closed anew
//sfpim      = XWF.sfpim;           	// Always opened/closed anew
  sfname     = XWF.sfname;		// Copy the file name
  sbigend    = XWF.sbigend;		// Copy our endian type
  sbyteordin = XWF.sbyteordin;		// Copy file endian type
  stotpts    = XWF.stotpts;		// Copy total bytes (re+im)
  spadding   = XWF.spadding;		// Copy needed padding
  sfsize     = XWF.sfsize;		// Copy total file size
  sdata      = XWF.sdata;		// Copy data vector
  }

// ____________________________________________________________________________
// B                   XWinNMR Fid File Access Functions
// ____________________________________________________________________________

/* These functions allow users to get some simple information regarding the
   contents of the Bruker binary data acqusition files, 1r & 1i.             */

string     XWinSpec::specname()  const { return sfname; }
string     XWinSpec::specrname() const { return sfname + "r"; }
string     XWinSpec::speciname() const { return sfname + "i"; }
int        XWinSpec::size()      const { return stotpts; }
int        XWinSpec::SSI()       const { return stotpts; }
int        XWinSpec::bytes()     const { return sfsize; }
int        XWinSpec::pad()       const { return spadding; }
bool       XWinSpec::order()     const { return sbyteordin; }
row_vector XWinSpec::data()      const { return sdata; }

int XWinSpec::blocks() const
  {
  int npts = sfsize/4;				// Number 32-bit ints in file
  int ptsperspec = stotpts + spadding/4;	// Points used by each spectrum
  if(!ptsperspec) return 0;			// Insure no size overflow
  return npts/ptsperspec;			// Number of Spectra in the file
  }

// ____________________________________________________________________________
// C                       XWinSpec Input Functions
// ____________________________________________________________________________

/* This function will read in a binary spectrum from a Bruker 1D data set
   in XWinNMR.  Typically this is just done during construction of the class.
   However, it is convenient to have an empty class constructor which can be
   used in conjunction with read/write to perform the same feat.  The Bruker
   files "1r" and "1s" contain values that are read in sequentially, each point
   stored as a 32-bit integer.  It is very important that we know the number 
   of points in a spectrum and the byte order used in storing the data.  This
   information must be supplied, either when calling the read function or when
   constructing objects of this class type.  Such information is usually stored
   in corresponding ASCII parameter file, typically named proc and procs. Such
   details should be attended to in higher classes utilizing these functions.*/ 
 
bool XWinSpec::read(const string& name, bool byteord, int SI, int warn)
  {
  string rname = name + "r";
  string iname = name + "i";
  sfpre.open(rname.c_str(), ios::binary|ios::in);
  sfpim.open(iname.c_str(), ios::binary|ios::in);
  if(!sfpre.good())			// Insure Bruker File is O.K.
    {
    if(warn)
      {
      XWSerror(28,1);			//   Problems with filestream
      XWSerror(1,rname,1);		//   Problems with reals file
      if(warn>1) XWSfatal(29);		//   Can't read a damn thing
      else       XWSerror(1,rname);
      }
    return false;
    }
  if(!sfpim.good())			// Insure Bruker File is O.K.
    {
    if(warn)
      {
      XWSerror(28,1);			//   Problems with filestream
      XWSerror(1,iname,1);		//   Problems with imaginaries file
      if(warn>1) XWSfatal(29);		//   Can't read a damn thing
      else       XWSerror(1,iname);
      }
    return false;
    }
  sfname     = name;			// Set the Bruker filename
  sbigend    = WeRBigEnd();		// Set output byte order (our arch)
  sbyteordin = byteord;			// Set input byte order (Bruker)
  sfpre.seekp(0,ios::end);		// Go to the end of the file
  sfsize = sfpre.tellp();		// See how many bytes are in the file
  sfpre.seekp(0);			// Set things back to the start
  sfpim.seekp(0,ios::end);		// Go to the end of the file
  if(sfsize != sfpim.tellp())		// Insure file sizes match up
    {
    if(warn)
      {
      XWSerror(32,1);			//   Problems with file sizes
      if(warn>1) XWSfatal(29);		//   Can't read a damn thing
      else       XWSerror(1,name);
      }
    return false;
    }
  sfpim.seekp(0);			// Set things back to the start
  if(SI>0) stotpts= SI;			// Set points per spectrum (re or im)
  else     stotpts=sfsize/4;		// Set points per spectrum (re or im)
  SetPadding();				// Set any padding for data boundary
  if(!CheckSize(1))			// Insure file size O.K. if open
    XWSfatal(46);			// for writing, SI too large else
  int swapon = 0;			// Assume no byte swapping
  if(sbigend != sbyteordin) swapon = 1;	// Need to swap if mismatch
  sdata = row_vector(stotpts);		// Array for spectrum data
  int32_t ptre, ptim;			// These will be input values
  for(int i=0; i<stotpts; i++)		// Loop over spectrum points
    {  
    sfpre.read((char*)&ptre,sizeof(int32_t)); 	//   Read a (real) point
    sfpim.read((char*)&ptim,sizeof(int32_t)); 	//   Read a (imag) point
    if(swapon) {Swap(ptre); Swap(ptim);}//   Byte swap if needed
    sdata.put(complex(ptre,ptim),i);	//   Store point
    }  
  sfpre.close();			// Close the reals file
  sfpim.close();			// Close the imaginaries file
  return true;				// Return successfull
  }
 
 
// ____________________________________________________________________________
// D                       XWinSpec Output Functions
// ____________________________________________________________________________
 
/* This function will write a binary spectrum for a Bruker 1D acqusition
   in XWinNMR.  Typically this is just done during construction of the class.
   However, it is convenient to have an empty class constructor which can be
   used in conjunction with read/write to perform the same feat.  The Bruker
   files "1r" & "1i" contain real and imaginary points reapectively, each point
   stored as a 32-bit integer. The spectrum files should be consistent with the
   corresponding ASCII parameter files, proc and procs. However, such details
   should be attended to in higher classes utilizing these functions.        */ 
 
int XWinSpec::write(const string& name, const row_vector& vx, int warn)
  {
  sfpre.open((name+"r").c_str(), ios::binary|ios::out);
  sfpim.open((name+"i").c_str(), ios::binary|ios::out);
  if(!sfpre.good())			// Insure Bruker File is O.K.
    {
    if(warn)
      {
      XWSerror(28,1);			//   Problems with filestream
      XWSerror(1,name+"r",1);		//   Problems with file
      if(warn>2) XWSfatal(29);		//   Can't write a damn thing
      else       XWSerror(29);
      }
    return 0;
    }
  if(!sfpim.good())			// Insure Bruker File is O.K.
    {
    if(warn)
      {
      XWSerror(28,1);			//   Problems with filestream
      XWSerror(1,name+"i",1);		//   Problems with file
      if(warn>2) XWSfatal(29);		//   Can't write a damn thing
      else       XWSerror(29);
      }
    return 0;
    }
  sfname     = name;			// Set the Bruker filename
  sbigend    = WeRBigEnd();		// Set output byte order (our arch)
  sbyteordin = sbigend;			// Set input byte order (our arch)
  sdata      = vx;			// Set our data points
  stotpts    = vx.size();		// Set points per spectrum (re or im)
  SetPadding();				// Set padding for data boundary
  sfpre.seekp(0);			// Set things at file start
  int npts   = sdata.size();		// Number of points (SI)
  int32_t rval, ival;			// For point reals,imags
  for(int i=0; i<npts; i++)		// Loop over all points
    {
    rval = int32_t(vx.getRe(i));		//   Real point as int
    ival = int32_t(vx.getIm(i));		//   Imag point as int
    sfpre.write((char*)&rval, sizeof(int32_t));	//   Write real point
    sfpim.write((char*)&ival, sizeof(int32_t));	//   Write imag point
    }
  AddPadding();				// Pad to end of block
  sfpre.seekp(0,ios::end);		// Go to the end of the file
  sfsize = sfpre.tellp();		// Get total bytes in the file
  sfpre.close();			// Close the reals file
  sfpim.close();			// Close the imaginaries file
  return 1;
  }
 
// ____________________________________________________________________________
// E                    XWinSpec Standard Output Functions
// ____________________________________________________________________________
 
/* These functions allow users to have a quick peek at the "spectrum" contents.
   They do not do any manipulations to the files whatsoever, they will only
   report on what in known about the data values and file sizes.             */ 

ostream& XWinSpec::print(ostream& ostr, int full, int hdr) const
  {
  if(hdr)
    ostr << "\n\n\t\tXWinNMR Binary Spectrum File " << sfname << "\n";
  ostr << "\n\t\tFile Name:              " << sfname;
  ostr << "\n\t\tReals File:             " << sfname + "r";
  ostr << "\n\t\tImaginaries File:       " << sfname + "r";
  ostr << "\n\t\tTotal Points (Re/Im):   " << stotpts;
  ostr << "\n\t\tTotal Bytes:            " << sfsize;
  if(spadding)
    ostr << "\n\t\tAdded Padding:          " << spadding;
  ostr << "\n\t\tFile Byte Ordering:     ";
  if(!sbyteordin) ostr << "Little Endian";
  else            ostr << "Big Endian";
  double imax=-1.e-15, imin=1.e15;
  double rmax=-1.e-15, rmin=1.e15;
  double nmax=0,       nmin=1.e15;
  double n, vi, vr;
  for(int i=0; i<sdata.size(); i++)
    {
    n  = norm(sdata.get(i));
    vi = sdata.getIm(i);
    vr = sdata.getRe(i);
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
      ostr << "\n\t\t\t" << k << ". " << sdata.get(k);
    }
  return ostr;
  }

ostream& operator<< (ostream& O, const XWinSpec& F) { F.print(O); return O; }

#endif						// XWinSpec.cc
