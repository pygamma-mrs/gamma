/* XWinSer.cc ***************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A                           **
**                                                                      **
**      XWinSer						Implementation	**
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

#ifndef _XWinSer_cc_			// Is this file already included?
#define _XWinSer_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler
#    pragma implementation
#endif

#include <GamIO/XWinSer.h>			// Include interface
#include <GamIO/BinIOBase.h>			// Include binary I/O functions
#include <Basics/Gutils.h>                      // Include GAMMA errors

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                      XWinNMR Serial File Error Handling
// ____________________________________________________________________________

/* These functions take care of any errors encountered when reading, writing,
   and setting parameters in Bruker acquisition parameter files.

        Input           FidPar  : XWinNMR acqusition parameters (this)
                        eidx    : Error index
                        pname   : Additional error message
                        noret   : Flag for linefeed (0=linefeed)
        Output          void    : An error message is output                 */  


void XWinSer::XWinSererror(int eidx, int noret) const
  {
  std::string hdr("Bruker Serial File");
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
     

void XWinSer::XWinSererror(int eidx, const std::string& pname, int noret) const
  {
  std::string hdr("Bruker Serial File");
  switch(eidx)
    {
    case 1:  GAMMAerror(hdr,   1, pname, noret); break; // File Problems   (1)
    case 10: GAMMAerror(hdr,"Can't Read Array "+pname, noret); break;    //(10)
    default: GAMMAerror(hdr,  -1, pname, noret); break; // Unknown Error   (-1)
    }
  }
     
volatile void XWinSer::XWinSerfatality(int eidx) const
  {
  XWinSererror(eidx, 1);
  if(eidx) XWinSererror(0);
  GAMMAfatal();					// Clean exit from program
  }
     
volatile void XWinSer::XWinSerfatality(int eidx, const std::string& pname) const
  {
  XWinSererror(eidx, pname, 1);
  if(eidx) XWinSererror(0);
   GAMMAfatal();					// Clean exit from program
  }

// ----------------------------------------------------------------------------
// ii                      Mode & Size Checking Functions
// ----------------------------------------------------------------------------

/* These functions check the mode specified when opening a XWinSer MAT file.
   We do not allow all of the possible modes availble to C++ filestreams.  In
   particular, we do not allow users to append to an existing file (since this
   would be bad if the existing data was written with a different byte order).
   Additonally, there seems to be some discrepancy between compilers as to
   what modes are available to fstream, so we just stick to the basics.

	in  = 1 (0x01)	out       = 2  (0x02)	ate      = 4   (0x03)
	app = 8	(0x08)	trunc     = 16 (0x10)	nocreate = 32  (0x20)
			noreplace = 64 (0x40) 	binary   = 128 (0x80)
   
   Use of noreplace & nocreate maybe didn't make it to ANSI standard...      */

void XWinSer::CheckMode(int mode)
  {
  if(!(mode & std::ios::binary))
    {
    XWinSererror(20, 1);			// File can't be ASCII
    XWinSerfatality(21);			// File must be binary
    }
  if((mode & std::ios::out) && (mode & std::ios::in))
    {
    XWinSererror(22, 1);			// File can't be input & output
    XWinSerfatality(23);			// Open file as either input or output
    }
  if(mode & std::ios::ate)
    {
    XWinSererror(24, 1);			// File can't open at end
    XWinSerfatality(25);			// Must be opened without std::ios::ate
    }
  if(mode & std::ios::app)
    {
    XWinSererror(26, 1);			// File can't open in append mode
    XWinSerfatality(27);			// Must be opened with ios:app
    }
  }

int XWinSer::CheckSize(int warn)
  {
  int fidbytes = 4*ftotpts;		// Number of bytes per FID
  if(fidbytes > ffsize)			// Insure 1 FID is NOT bigger than 
    {					// the entire file size
    if(warn)
      {
      XWinSererror(44, 1);		// File size problems
      if(warn > 1) XWinSerfatality(45); 	// Block size too big
      else         XWinSererror(45, 1);
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
        XWinSererror(44, 1);		// File size problems
        if(warn > 1) XWinSerfatality(45);// Block size too big
        else         XWinSererror(45, 1);
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

void XWinSer::SetPadding()
  {
  fpadding = 0;				// Assume there is no padding
  if(ftotpts == 0) return;		// No padding if no points
  fpadding = ftotpts; 			// Begin will points per FID 
  while(fpadding >= 256) fpadding-=256;	// Clip chunks of 256 points
  if(fpadding < 0) XWinSerfatality(30);	// Quit if can't find boundary
  if(fpadding > 0)			// If padding needed, set it to number
    fpadding = 4*(256-fpadding);	// of bytes to reach 256 point boundary
  }

bool XWinSer::CheckBoundary()
  {
  int pos = ffp.tellp();		// Get current file position
  if(!pos) return 1;			// Start of file is boundary
  if(!ftotpts) XWinSerfatality(31);	// Quit if can't check boundary
  int blkbytes = 4*ftotpts+fpadding;	// Total bytes per "padded" FID
  while(pos >= blkbytes) pos-=blkbytes;	// Clip off blocks
  if(!pos) return 1;			// This is zero if a boundary
  if(pos < 0) XWinSerfatality(31);	// Problems if negative
  return 0; 				// Here if NOT a boundary
  }

void XWinSer::SkipPadding() { ffp.seekp(fpadding, std::ios::cur); }
  
void XWinSer::AddPadding()
  {
  char pc = ' ';			// Padding charactor
  for(int i=0; i<fpadding; i++) 	// Write blanks, assuming there
    ffp.write(&pc, sizeof(char));	// is only 1 byte per character
  }
 
// ----------------------------------------------------------------------------
// iv                    File Default Settings
// ----------------------------------------------------------------------------

void XWinSer::SetDefaults()
  {
  ffname     = std::string("");	// Bruker Ser/Serial file name empty
  fbigend    = WeRBigEnd();	// Set output byte order (our arch)
  fbyteordin = 1; 		// Set input byte order big (as from SGI)
  ftotpts    = 0;		// Total points per Ser (real + imag)    
  fpadding   = 0;		// Bytes needed for padding to 256 pt boundary
  ffsize     = 0;		// Our total file size (bytes)
  if(ffp.is_open())		// Close binary file stream if open
    ffp.close();     
  }


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                    XWinSer CONSTRUCTORS, DESTRUCTOR
// ____________________________________________________________________________


XWinSer::XWinSer() { SetDefaults(); }

XWinSer::XWinSer(const std::string& name, int TD, bool byteord, int mode)
       :ffp(name.c_str(), Int2Mode(mode))
  {
  if(!ffp.good())			// Insure Bruker File is O.K.
    {
    XWinSererror(28,1);			//   Problems with filestream
    XWinSererror(1,name,1);		//   Problems with file
    XWinSerfatality(29);		//   Can't read a damn thing
    }
  CheckMode(mode);			// Check if mode is acceptable
  ffname     = name;			// Set the Bruker filename
  fbigend    = WeRBigEnd();		// Set output byte order (our arch)
  fbyteordin = byteord;			// Set input byte order (Bruker)
  ftotpts    = TD;			// Set points per fid (re + im)
  SetPadding();				// Set any padding for data boundary
  ffp.seekp(0,std::ios::end);		// Go to the end of the file
  ffsize = ffp.tellp();			// See how many bytes are in the file
  ffp.seekp(0);				// Set things back to the start
  if((mode & std::ios::in) && !CheckSize(1))	// Insure file size O.K. if open
    XWinSerfatality(46);			// for writing, TD too large else
  }

XWinSer::XWinSer(const XWinSer& XWF) 
  {
// sosi - this is a problem for constant XWF
//  ffp        = XWF.ffp;			// Input/Output binary file stream
  ffname     = XWF.ffname;		// Bruker Ser/Serial file name
  fbigend    = XWF.fbigend;		// Flag for big- vs. little-endian (this arch)
  fbyteordin = XWF.fbyteordin;		// Flag for big- vs. little-endian input
  ftotpts    = XWF.ftotpts;		// Total points per Ser (real + imag)
  fpadding   = XWF.fpadding;		// Bytes needed for padding to 256 pt boundary
  ffsize     = XWF.ffsize;		// Our total file size (bytes)
  }

XWinSer::~XWinSer() {}
  
        // Input              	XWF	: XWinNMR Binary FID file
	// Output		none	: XWF is destructed
	// Note			  	: Of course, this doesn't touch
	//				  the external file at all

XWinSer& XWinSer::operator=(const XWinSer& XWF)
  {
  if(this == &XWF) return *this;	// Avoid self copying
// sosi - this is a problem for constant XWF
//  ffp        = XWF.ffp;			// Input/Output binary file stream
  ffname     = XWF.ffname;		// Bruker Ser/Serial file name
  fbigend    = XWF.fbigend;		// Flag for big- vs. little-endian (this arch)
  fbyteordin = XWF.fbyteordin;		// Flag for big- vs. little-endian input
  ftotpts    = XWF.ftotpts;		// Total points per Ser (real + imag)
  fpadding   = XWF.fpadding;		// Bytes needed for padding to 256 pt boundary
  ffsize     = XWF.ffsize;		// Our total file size (bytes)
  return *this;
  }

// ____________________________________________________________________________
// B                 XWinNMR Fid File Auxiliary Functions
// ____________________________________________________________________________

/* These functions allow users to get some simple intformation regarding the
   contents of the Bruker serial acquisition file.                           */

std::string XWinSer::sername()  const { return ffname; }
int    XWinSer::size()     const { return ftotpts/2; }
int    XWinSer::TDS()      const { return ftotpts; }
int    XWinSer::bytes()    const { return ffsize; }
int    XWinSer::pad()      const { return fpadding; }
bool   XWinSer::order()    const { return fbyteordin; }

int XWinSer::blocks() const
  {
  int npts = ffsize/4;				// Number 32-bit ints in file
  int ptsperfid = ftotpts + fpadding/4;		// Points used by each FID
  return npts/ptsperfid;			// Number of FIDs in the file
  }

// ____________________________________________________________________________
// C                         XWinSer Input Functions
// ____________________________________________________________________________

/* These functions will read in a single FID from the Bruker acqusition file
   and or a series of FID's from a Bruker serial file. The FID values are read
   in as re,im,re,im,...., each point stored as a 32-bit integer.  It is very
   important that we know the number of points in an FID and the byte order
   used in storing the data.  This information must be supplied, either when
   calling the function(s) or when constructing objects of this class type.
   Such information is usually read in from a corresponding ASCII parameter
   file, typically named acqus in the case of a sinlge FID file (fid) or 
   called acqu2s in case of a stacked FID serial file (ser). Such details
   should be attended to in higher classes utilizing these functions.

   Function   Arguments              Action                       Return
   --------  ------------  -----------------------------------  ----------
   readFID   fn,TD,BO,I    Read Ith block in serial/FID file    row_vector
 + readFID   idx           Read block of index idx(-1=current)  row_vector
   readFIDs  fn,TD,BO,I,J  Read a set of J FID's from Ith one     matrix
 + readFIDs  I,NB          Read a set of NB FID's from Ith one    matrix
   readSer   fn,TD,bytord  Read entire serial file                matrix
 + readSer                 Read entire serial file                matrix
 
 Above X=Currently Inactive, +=File Must Be Properly Opened For Reading      */
 

row_vector XWinSer::readFID(const std::string& fin, int TD, int bord, int idx)
  {
  bool byteord = (bord)?true:false;
  if(bord == -1) bord = WeRBigEnd();	// Our byte order if not set
  XWinSer XWA(fin,TD,byteord,std::ios::binary|std::ios::in);	// Open Bruker serial file
  return XWA.readFID(idx);			// Use function overload
  }
 
row_vector XWinSer::readFID(int idx)
  {
  if(idx >= 0)					// If we want a specific
    {						// block, then adjust pointer
    int blocksize = 4*ftotpts+fpadding;		// This is block size in bytes
    ffp.seekp(idx*blocksize, std::ios::beg);		// Skip idx blocks in file
    }
  if(!CheckBoundary())				// If we are not on a boundary 
    {						// then there is trouble
    XWinSererror(40, 1);			//   Read off block boundary
    XWinSerfatality(41);			//   Bad read of FID
    }
  int swappts = 0;				// Assume byte no swapping
  if(fbigend != fbyteordin) swappts = 1;	// Need to swap if mismatch
  row_vector data(ftotpts/2);			// Array for fid data
  int32_t ptre, ptim;				// These will be input values
  for(int i=0; i<ftotpts/2; i++)		// Loop over fid points
    {  
    ffp.read((char*)&ptre,sizeof(int32_t)); 	// Read a (real) point
    ffp.read((char*)&ptim,sizeof(int32_t)); 	// Read a (imag) point
    if(swappts) { Swap(ptre); Swap(ptim); }	// Byte swap if needed
    data.put(complex(ptre,ptim), i);		// Store point
// sosi - check for end of file here & issue warning
    }  
  SkipPadding();				// Set file to next boundary
  return data;					// Return the FID as row_vector
  }

matrix XWinSer::readFIDs(const std::string& fin, int TD, int BO, int I, int J)
  {
   bool byteord = (BO)?true:false;
  if(BO == -1) byteord = WeRBigEnd();		// Our byte order if not set
  XWinSer XA(fin,TD,byteord,std::ios::binary|std::ios::in);	// Open Bruker serial file
  return XA.readFIDs(I, J);			// Use function overload
  }

matrix XWinSer::readFIDs(int idx, int NB)
  {
  int swappts = 0;				// Assume byte no swapping
  if(fbigend != fbyteordin) swappts = 1;	// Need to swap if mismatch
  int npts = ftotpts/2;				// Get number of complex pts
  matrix data(NB, npts);			// Array for fid data
  int blocksize = 4*ftotpts+fpadding;		// This is block size in bytes
  ffp.seekp(idx*blocksize, std::ios::beg);		// Skip idx blocks in file
  int32_t ptre, ptim;				// These will be input values
  int i,j;					// Element indices
  for(i=0; i<NB; i++)				// Loop over fid blocks
    {
    for(j=0; j<npts; j++)			// Loop over fid points
      {  
      ffp.read((char*)&ptre,sizeof(int32_t)); 	// Read a (real) point
      ffp.read((char*)&ptim,sizeof(int32_t)); 	// Read a (imag) point
      if(swappts) { Swap(ptre); Swap(ptim); }	// Byte swap if needed
      data.put(complex(ptre,ptim), i,j);	// Store point
      }  
    SkipPadding();				// Set file to next boundary
    }
  return data;
  }
 
matrix XWinSer::readSer(const std::string& fin, int TD, int bytord)
  {
  bool byteord = (bytord)?true:false;
  if(bytord == -1) byteord = WeRBigEnd();	// Our byte order if not set
  XWinSer XA(fin,TD,byteord,std::ios::binary|std::ios::in);	// Open Bruker serial file
  return XA.readSer();				// Use function overload
  }

matrix XWinSer::readSer()
  {
  int swappts = 0;				// Assume byte no swapping
  if(fbigend != fbyteordin) swappts = 1;	// Need to swap if mismatch
  int nblks = blocks();				// Get number of FIDs
  int npts = ftotpts/2;				// Get number of complex pts
  ffp.seekp(std::ios::beg);				// Go to file beginning
  int32_t ptre, ptim;				// These will be input values
  int i,j;					// Looping variables
  matrix data(nblks, npts);			// Array for fid data
  for(i=0; i<nblks; i++)			// Loop over fid blocks
    {
    for(j=0; j<npts; j++)			// Loop over fid points
      {  
      ffp.read((char*)&ptre,sizeof(int32_t)); 	// Read a (real) point
      ffp.read((char*)&ptim,sizeof(int32_t)); 	// Read a (imag) point
      if(swappts) { Swap(ptre); Swap(ptim); }	// Byte swap if needed
      data.put(complex(ptre,ptim),i,j);		// Store point
      }  
    SkipPadding();				// Set file to next boundary
    }
  return data;					// Return the FID as matrix
  }

row_vector XWinSer::readSlice(const std::string& fin, int TD, int bord, int idx)
  {
  bool byteord = (bord)?true:false;
  if(bord == -1) byteord = WeRBigEnd();		// Our byte order if not set
  XWinSer XWA(fin,TD,byteord,std::ios::binary|std::ios::in);	// Open Bruker serial file
  return XWA.readSlice(idx);			// Use function overload
  }
 
row_vector XWinSer::readSlice(int idx)
  {
  int swappts = 0;				// Assume byte no swapping
  if(fbigend != fbyteordin) swappts = 1;	// Need to swap if mismatch
  int nblks = blocks();				// Get number of FIDs
  int blocksize = 4*ftotpts+fpadding;		// This is block size in bytes
  ffp.seekp(8*idx, std::ios::beg);			// Skip idx (complex) points
  row_vector data(nblks);			// Array for fid data
  int32_t ptre, ptim;				// These will be input values
  for(int i=0; i<nblks; i++)			// Loop over fid points
    {  
    ffp.read((char*)&ptre,sizeof(int32_t)); 	// Read a (real) point
    ffp.read((char*)&ptim,sizeof(int32_t)); 	// Read a (imag) point
    if(swappts) { Swap(ptre); Swap(ptim); }	// Byte swap if needed
    data.put(complex(ptre,ptim),i);		// Store point
    ffp.seekp(blocksize-8, std::ios::cur);		// Skip to next block, next pt
    }  
  return data;					// Return the FID as row_vector
  }

row_vector XWinSer::readSlices(const std::string& F, int TD, int BO, int I, int J)
  {
  bool byteord = (BO)?true:false;
  if(BO == -1) byteord = WeRBigEnd();		// Our byte order if not set
  XWinSer XWA(F,TD,byteord,std::ios::binary|std::ios::in);	// Open Bruker serial file
  return XWA.readSlices(I,J);			// Use function overload
  }
 
row_vector XWinSer::readSlices(int idx, int NB)
  {
  int swappts = 0;				// Assume byte no swapping
  if(fbigend != fbyteordin) swappts = 1;	// Need to swap if mismatch
  int nblks = blocks();				// Get number of FIDs
  int blocksize = 4*ftotpts+fpadding;		// This is block size in bytes
  blocksize -= 4*2*NB;				// Subtract complex pts/block
  ffp.seekp(8*idx, std::ios::beg);			// Skip idx (complex) points
  int i, j;					// Row, column point indices
  int32_t ptre, ptim;				// These will be input values
  matrix data(nblks, NB);			// Array for slice data
  for(j=0; j<nblks; j++)			// Loop over fid blocks
    for(i=0; i<NB; i++)				// Loop over (complex) points
      {  
      ffp.read((char*)&ptre,sizeof(int32_t)); 	// Read a (real) point
      ffp.read((char*)&ptim,sizeof(int32_t));	// Read a (imag) point
      if(swappts) { Swap(ptre); Swap(ptim); }	// Byte swap if needed
      data.put(complex(ptre,ptim),i,j);		// Store point
      ffp.seekp(blocksize, std::ios::cur);		// Skip to next block, next pt
      }  
  return data;					// Return slices as matrix
  }

// ____________________________________________________________________________
// D                       XWinSer Output Functions
// ____________________________________________________________________________
 
/* These functions will write an FID or series of FIDs to a file in Bruker
   XWinNMR format. This can be done either by directly writing a row vector
   or matrix into such a file, or by writing a series of row vectors. Note that
   ususally there is a parameter file associated with such data files, called
   acqus in the case of a single FID file or called acqu2s for a serial file
   file of multiple FIDs. NOTE: The input data from GAMMA (row_vector or
   matrix) contains double precision numbers which are converted to long
   integers. Any numbers of magnidude less than 1 will become zero.  Thus
   it is very wise to SCALE THE INPUT DATA (~10^6 - ~10^8) so that one will
   not loose any informaton prior to output in Bruker format               .*/

int XWinSer::write(const std::string& fout, const row_vector& data, int warn)
  {
  bool bord = WeRBigEnd();			// Our byte order 
  int TD = 2*data.size();			// FID point size
  XWinSer XA(fout,TD,bord,std::ios::binary|std::ios::out);// Open Bruker serial file
  int TF = XA.write(data);			// Write vector to file
  ffp.close();					// Close the file
  return TF;
  }

int XWinSer::write(const row_vector& vx, int warn)
  {
  if(!CheckBoundary())				// If we are not on a boundary 
    {						// then there is trouble
    if(warn)
      {
      XWinSererror(42, 1);			//   Write off block boundary
      if(warn>1) XWinSerfatality(43);		//   Bad write of FID
      else       XWinSererror(43,1);
      }
    return 0;
    }
  int npts = vx.size();				// Number of points (2*TD)
  int32_t rval, ival;				// For point reals,imags
  for(int i=0; i<npts; i++)			// Loop over all points
    {
    rval = int32_t(vx.getRe(i));			//   Real point as int
    ival = int32_t(vx.getIm(i));			//   Imag point as int
    ffp.write((char*)&rval, sizeof(int32_t));		//   Write real point
    ffp.write((char*)&ival, sizeof(int32_t));		//   Write imag point
    }
  AddPadding();					//   Pad to end of block
  return 1;
  }

int XWinSer::write(const std::string& fout, const matrix& data, int warn)
  {
  bool bord = WeRBigEnd();			// Our byte order 
  int TD = 2*data.cols();			// FID point size
  XWinSer XA(fout,TD,bord,std::ios::binary|std::ios::out);	// Open Bruker serial file
  int TF = XA.write(data);			// Write matrix to file
  ffp.close();					// Close the file
  return TF;
  }

int XWinSer::write(const matrix& data, int warn)
  {
  if(!CheckBoundary())				// If we are not on a boundary 
    {						// then there is trouble
    if(warn)
      {
      XWinSererror(42, 1);			//   Write off block boundary
      if(warn>1) XWinSerfatality(43);		//   Bad write of FID
      else       XWinSererror(43,1);
      }
    return 0;
    }
  int npts  = data.cols();			// Number of points (2*TD)
  int nfids = data.rows();			// Number of FIDs
  int i,j;					// Loop indexing variables
  int32_t rval, ival;				// For point reals,imags
  for(i=0; i<nfids; i++)			// Loop over FIDs
    {
    for(j=0; j<npts; j++)			// Loop over points
      {
      rval = int32_t(data.getRe(i,j));		// Real point as int
      ival = int32_t(data.getIm(i,j));		// Imag point as int
      ffp.write((char*)&rval, sizeof(int32_t));		// Write real point
      ffp.write((char*)&ival, sizeof(int32_t));		// Write imag point
      }
    AddPadding();				//   Pad to end of block
    }
  return 1;
  }
 
// ____________________________________________________________________________
// E                    XWinSer Standard Output Functions
// ____________________________________________________________________________
 
/* These functions allow users to have a quick peek at the "ser" contents.
   They do not do any manipulations to the file whatsoever, they will only
   report on what in known about the data values and file sizes.             */ 

std::ostream& XWinSer::print(std::ostream& ostr, int full, int hdr)
  {
  if(hdr)
    ostr << "\n\n\t\tXWinNMR Binary Serial File " << ffname << "\n";
  ostr << "\n\t\tFile Name:              " << ffname;
  ostr << "\n\t\tBlock Points (Re+Im):   " << ftotpts;
  ostr << "\n\t\tNumber of Blocks:       " << blocks();
  ostr << "\n\t\tTotal Bytes:            " << ffsize;
  if(fpadding)
    ostr << "\n\t\tAdded Padding:          " << fpadding;
  else
    ostr << "\n\t\tAdded Padding:          none";
  ostr << "\n\t\tFile Byte Ordering:     ";
  if(!fbyteordin) ostr << "Little Endian";
  else            ostr << "Big Endian";
  if(full)
    {
    double imax=-1.e-15, imin=1.e15;
    double rmax=-1.e-15, rmin=1.e15;
    double nmax=0,       nmin=1.e15;
    double n, vi, vr;
    complex z;
    int swappts = 0;                              // Assume byte no swapping
    if(fbigend != fbyteordin) swappts = 1;        // Need to swap if mismatch
    int nblks = blocks();                         // Get number of FIDs
    int npts = ftotpts/2;                         // Get number of complex pts
    ffp.seekp(std::ios::beg);                          // Go to file beginning
    int32_t ptre, ptim;                              // These will be input values
    int i,j;                                      // Looping variables
    for(i=0; i<nblks; i++)                        // Loop over fid blocks
      {
      for(j=0; j<npts; j++)                       // Loop over fid points
        {
	ffp.read((char*)&ptre,sizeof(int32_t));             // Read a (real) point
	ffp.read((char*)&ptim,sizeof(int32_t));             // Read a (imag) point
	if(swappts) { Swap(ptre); Swap(ptim); }   // Byte swap if needed
	z = complex(ptre,ptim); 	          // Store point
        n  = norm(z);
        vi = Im(z);
        vr = Re(z);
        if(n  > nmax) nmax = n;
        if(n  < nmin) nmin = n;
        if(vi > imax) imax = vi;
        if(vi < imin) imin = vi;
        if(vr < rmin) rmin = vr;
        if(vr > rmax) rmax = vr;
	} 
      SkipPadding();                              // Set file to next boundary
      }
    ostr << "\n\t\tData Norm Max/Min:      "
         << nmax << "/" << nmin;
    ostr << "\n\t\tData Real Max/Min:      "
         << rmax << "/" << rmin;
    ostr << "\n\t\tData Imaginary Max/Min: "
         << imax << "/" << imin;
    }
  return ostr;
  }

std::ostream& operator<< (std::ostream& O, XWinSer& F) { F.print(O); return O; }

// ____________________________________________________________________________
// F                  XWinSer File Handling Functions
// ____________________________________________________________________________

bool XWinSer::open(const std::string& name, int TD, bool byteord, int mode, int warn)
  {
  SetDefaults();			// Insure all is empty
  ffp.open(name.c_str(),Int2Mode(mode));// Try and open file
  if(!ffp.good())			// Insure Bruker File is O.K.
    {
    if(warn)
      {
      XWinSererror(28,1);		//   Problems with filestream
      XWinSererror(1,name,1);		//   Problems with file
      if(warn>1) XWinSerfatality(29);	//   Can't read a damn thing
      else       XWinSererror(29);
      }
    return false;
    }
  CheckMode(mode);			// Check if mode is acceptable
  ffname     = name;			// Set the Bruker filename
  fbigend    = WeRBigEnd();		// Set output byte order (our arch)
  fbyteordin = byteord;			// Set input byte order (Bruker)
  ftotpts    = TD;			// Set points per fid (re + im)
  SetPadding();				// Set any padding for data boundary
  ffp.seekp(0,std::ios::end);		// Go to the end of the file
  ffsize = ffp.tellp();			// See how many bytes are in the file
  ffp.seekp(0);				// Set things back to the start
  if((mode & std::ios::in) && !CheckSize(1))	// Insure file size O.K. if open
    XWinSerfatality(46);		// for writing, TD too large else
  return true;
  }
 

#endif						// XWinSer.cc
