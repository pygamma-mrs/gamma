/* XWin2D.cc ****************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A                           **
**                                                                      **
**      XWin2D				     Implementation		**
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
** The XW* files provide an interface to Bruker XWinNMR (uxnmr) data    **
** sets. This class embodies a Bruker data set for a 2D acquisition.    **
** Each 2D acquisition generates several files (binary and ASCII) 	**
** spanning several directories.  This class is intended to handle 	**
** this structure.  Below is a typical Bruker directory structure 	**
** associated with a 2D acquisiton.					**
**									**
**                          __ acqu, acqu2 (changable parameter files)	**
**			   / 						**
**                        /___ acqus, acqu2s (static parameter files)	**
**			 /						**
**  expname -- expnum --< ---- ser (binary data)			**
**			 \						**
**			  \___ pdata -- 1 -- proc, proc2, procs, proc2s **
**									**
** This class will handle the directory hierarchy shown above as	**
** well as most of the files therein. Thus, given a base directory 	**
** for an experiment/simulation output (expname) & an experiment	**
** number (expnum) this class will generate the subdirectories shown,	**
** the acquisition parameter (acqu*) and binary (ser) files as well as  **
** the processing files (proc*, meta).  It is NOT intended to generated **
** processed 2D/ND data files - that should be done using the XWinNMR   **
** software.  The class will also allow for import of 2D data sets.     **
** The serial file (ser), & to a more limited extent processed spectra, **
** can be read and placed into GAMMA row_vectors and matrices. The      **
** associated parameters are also provided for.                         **
**									**
*************************************************************************/

#ifndef   XWin2D_cc_			// This file already included?
#  define XWin2D_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler
#    pragma implementation
#  endif

#include <string>			// Include libstdc++ strings
#include <stdlib.h>
#include <GamIO/XWin2D.h>		// Include interface
#include <GamIO/XWinAcqus.h>		// Include acqus file interface
#include <GamIO/XWinAcqu2s.h>		// Include acqu2s file interface
#include <GamIO/XWinSer.h>		// Include serial file interface
#include <GamIO/XWinOutd.h>		// Include outd file interface
#include <GamIO/XWinProcs.h>		// Include procs file interface
#include <GamIO/BinIOBase.h>		// Include binary I/O functions
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <Basics/StringCut.h>		// Include GAMMA Gdec function
#include <Basics/Gconstants.h>          // Include HZ2RaD & GAMMA1H
#include <Basics/Isotope.h>		// Include gamma (gyrmag.) function
#include <sys/stat.h>			// Include POSIX mkdir function
#if defined(_MSC_VER)                   // If not using MSVC++ then we                              // and if using MSVC++ then I guess
 #include <direct.h>			//   Include mkdir function
#endif

using std::string;			// Using libstdc++ strings
using std::ostream;			// Using libstdc++ output streams
using std::ios;				// Using libstdc++ file type settings

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// This patches the difference between GCC & MSVC++ mkdir functions
/*
#if defined(_MSC_VER) || defined(__SUNPRO_CC)  
  int MakeADir(const string& dname, int no)
    { return mkdir(dname.c_str()); }
#else
  int MakeADir(const string& dname, int no)
    { return mkdir(dname.c_str(),  no); }
#endif
*/

// ____________________________________________________________________________ 
// i                   XWinNMR 2D Data Set Error Handling
// ____________________________________________________________________________
 
/* These functions take care of any errors encountered when reading, writing,
   and setting parameters in Bruker acquisition data sets.
 
        Input           eidx    : Error index
                        pname   : Additional error message
                        noret   : Flag for linefeed (0=linefeed)
        Output          void    : An error message is output                 */  

void XWin2D::XWin2Derror(int eidx, int noret) const                       
  {
  string hdr("XWinNMR 2D Data Set");
  string msg;
  switch (eidx)
    {
    case 20: GAMMAerror(hdr,"Cannot Read/Write Data Set",    noret); break; // (20)
    case 21: GAMMAerror(hdr,"Cannot Read Full Data Set",     noret); break; // (21)
    case 22: GAMMAerror(hdr,"Cannot Generate Full Data Set", noret); break; // (22)
    case 23: GAMMAerror(hdr,"Cannot Write Block To Disk",    noret); break; // (23)
    case 24: GAMMAerror(hdr,"Block Size Inconsistent",       noret); break; // (24)
    case 25: GAMMAerror(hdr,"Accessed FID Out of Range",     noret); break; // (25)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }  


void XWin2D::XWin2Derror(int eidx, const string& pname, int noret) const
  {
  string hdr("XWinNMR 2D Data Set");
  string msg;
  switch(eidx)
    {
    case 5: msg = string("Bad Use Of ") + pname + string(" Function ");
             GAMMAerror(hdr,msg+pname,noret);  break;                  // (5)
    case 21:msg = string("Problems With Files in Directory ");
             GAMMAerror(hdr,msg+pname,noret);  break;                  // (21)
    case 22:msg = string("Problems With Use of Directory ");
             GAMMAerror(hdr,msg+pname,noret);  break;                  // (22)
    case 23:msg = string("Problems Reading Parameter File ");
             GAMMAerror(hdr,msg+pname,noret);  break;                  // (23)
    case 25:msg = string("FID Range Is [0,") + pname + string("]");
             GAMMAerror(hdr,msg,noret);  break;                        // (25)
    default: GAMMAerror(hdr, eidx, pname, noret); break;// Unknown Error  (-1)
    }
  }  
 
volatile void XWin2D::XWin2Dfatality(int eidx) const
  {
  XWin2Derror(eidx, 1);			// Normal non-fatal error
  if(eidx) XWin2Derror(0);		// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

volatile void XWin2D::XWin2Dfatality(int eidx, const string& pname) const
  {
  XWin2Derror(eidx, pname, 1);		// Normal non-fatal error
  if(eidx) XWin2Derror(0);		// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________ 
// ii            XWinNMR 2D Data Set File & Directory I/O Handling
// ____________________________________________________________________________
 
int XWin2D::CheckDir(int TF, int warn, const string& dout) const
  {
  if(!TF) return 1;				// Good mkdir returns 0
  string cmd = "cd " + dout;			// Perhaps file is from
  if(!system(cmd.c_str())) return 1;		// pre-exisiting dir.
  if(warn)
    {
    XWin2Derror(21, dout, 1);			// Problems with directory
    if(warn > 1) XWin2Dfatality(22);		// Can't generate data set
    else         XWin2Derror(22,1);
    }
  return TF;
  }


int XWin2D::CheckWrite(int TF, int warn, const string& dout) const
  {
  if(!TF)
    {
    if(warn)
      {
      XWin2Derror(21, dout, 1);			// Problems with files
      if(warn > 1) XWin2Dfatality(22);		// Can't generate data set
      else         XWin2Derror(22,1);
      }
    }
  return TF;
  }

// ____________________________________________________________________________ 
// iii                   XWinNMR 2D Data Set Setup Functions
// ____________________________________________________________________________

/* These functions are used to quickly set up values associated with the 2D
   data set.  In particular they handle setting up a standard Bruker directory
   structure and the file names in the data set.  Another function is used to
   insure the that the parameter sets are self-consistent. This should be 
   called before the data set (or at least the parameter files) are output.  */


void XWin2D::SetNames()
  {
  NAIdir  = dname + "/" + Gdec(Aidx);		// Where to find acq. data
  Nacqu   = NAIdir + "/acqu";			// Parameter file (acqu)
  Nacqus  = NAIdir + "/acqus";			// Parameter file (acqus)
  Nacqu2  = NAIdir + "/acqu2";			// Parameter file (acqu2)
  Nacqu2s = NAIdir + "/acqu2s";			// Parameter file (acqu2)
  Nser    = NAIdir + "/ser";			// Data      file (ser)
  NPdir   = NAIdir + "/pdata";			// Base processing dir
  NPIdir  = NPdir  + "/" + Gdec(Pidx);		// Where to find proc. data
  Nproc   = NPIdir + "/proc";			// Parameter file (proc)
  Nprocs  = NPIdir + "/procs";			// Parameter file (procs)
  Nproc2  = NPIdir + "/proc2";			// Parameter file (proc2)
  Nproc2s = NPIdir + "/proc2s";			// Parameter file (proc2s)
  Nmeta   = NPIdir + "/meta";			// Parameter file (meta)
  Noutd   = NPIdir + "/outd";			// Parameter file (outd)
  }

// sosi - To placate MSVC++'s different definition of mkdir, I have made the
//        function MakeADir (in BinIOBase.cc) that will call either Microsofts
//        or GCC's version.

int XWin2D::MakeDirs(int warn)
  {
  int TF = 1;
//  TF *= mkdir(dname.c_str(),  493);		// Try & make base directory
  TF *= MakeADir(dname,  493);		// Try & make base directory  CheckDir(TF, warn, dname);			// Insure directory OK
//  TF *= mkdir(NAIdir.c_str(), 493);		// Try & make acq directory
  TF *= MakeADir(NAIdir,  493);
  CheckDir(TF, warn, NAIdir);			// Insure directory OK
//  TF *= mkdir(NPdir.c_str(),  493);		// Try & make process dir.
  TF *= MakeADir(NPdir,  493);
  CheckDir(TF, warn, NPdir);			// Insure directory OK
//  TF *= mkdir(NPIdir.c_str(), 493);		// Try & make process dir.
  TF *= MakeADir(NPIdir,  493);
  CheckDir(TF, warn, NPIdir);			// Insure directory OK
  return TF;
  }


int XWin2D::ReadPars(int warn)
  {
  int TF  = Acqs.readAPar(Nacqus, warn?1:0);	// Read in acqus file
  if(!TF)
  { if(warn>1) XWin2Dfatality(23, Nacqus);	// Troubles with NAcqus
    else if(warn) XWin2Derror(23,Nacqus,1);

      TF *= Acq2s.readAPar(Nacqu2s, warn?1:0);	// Read in acqu2s file
  }
  if(!TF)
  { if(warn>1) XWin2Dfatality(23, Nacqu2s);	// Troubles with NAcqu2s
    else if(warn) XWin2Derror(23,Nacqu2s,1);

      TF *= Procs.readPPar(Nprocs, warn?1:0);	// Read in procs file
  }
  if(!TF)
  { if(warn>1) XWin2Dfatality(23, Nprocs);	// Troubles with procs 
    else if(warn) XWin2Derror(23,Nprocs,1);

      TF *= Proc2s.readPPar(Nproc2s,warn?1:0);	// Read in the proc2s file 
  }
  if(!TF)
  { if(warn>1) XWin2Dfatality(23, Nproc2s);	// Troubles with proc2s 
    else if(warn) XWin2Derror(23,Nproc2s,1);

  }
  if(TF) SetConsistent();
  return TF;
  }

/* This function checks the consistency between acqus, acqu2s, procs & proc2s.
   When these files are read (at least earlier XWinNMR) acqu2s may be missing
   some values that are assumed taken from acqus. Similarly when writing 
   such files GAMMA must insure this is the case.                            */

void XWin2D::SetConsistent()
  {
//                First Set t2-t1 Acqusition Parameters Consistent

  Acq2s.BF1(Acqs.BF2());			// Set F1 Om1 = F2 Om2  
  Acq2s.BF2(Acqs.BF1());			// Set F1 Om2 = F2 Om1  
  Acq2s.SFO1(Acqs.SFO2());			// Set F1 SFO1 = F2 SFO2  
  Acq2s.SFO2(Acqs.SFO1());			// Set F1 SFO2 = F2 SFO1  
  Acq2s.XW_IN(1, 1.0/Acq2s.SW_h());		// Insure IN0 = 1/SW

// sosi
//  int i;
//  for(i=3; i<=8; i++)				// Set NUCLEI to empty
//    NUC(i, string(""));				// on channels 3-8 acqu(s)
//  for(i=1; i<=8; i++)				// Set all NUCLEI to empty
//    Acq2s.NUC(i, string(""));			// on channels 1-8 acqu2s
//  Acq2s.NUC(1, NUC(2));				// Make sure 
//                Next Set t2-f2 Processing Parameters Consistent


  Procs.SF(Acqs.SFO1());			// Insure W2 SF matches
  Procs.SW_p(Acqs.SW_h());			// Insure W2 SW matches
  Procs.OFFSET(Acqs.O1()+Acqs.SW()/2);		// Insure W2 offset matches

//                Next Set t1-f1 Processing Parameters Consistent

  Proc2s.SF(Acq2s.SFO1());			// Insure W1 SF matches
  Proc2s.SW_p(Acq2s.SW_h());			// Insure W1 SW matches
  Proc2s.OFFSET(Acq2s.O1()+Acq2s.SW()/2);	// Insure W1 offset matches
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A         XWinNMR 2D Data Set Constructors, Destructor, Assignment
// ____________________________________________________________________________

/* These are the constructors of the class handling Bruker XWinNMR 2D data
   sets. There are several files associated with such a data set. Those
   associated with the acquisition are the four ASCII parameter files
  (acqu, acqus, acqu2, and acqu2s) and the binary serial data file (ser).
   In the processed data sub-directory, pdata will reside more files - whether 
   or NOT the data has actually been processed.  These will be the six
   parameter files (proc, procs, proc2, proc2s, meta, and outd).  If the data
   has been processed there may also exist more binary files (e.g. 2rr) if the
   data has been processed.  GAMMA has individual XWin* classes for handling
   these individual files, so this class has the job of keeping track of
   directories and making sure all these are appropriate files are dealt with.

   Note that these constructors do NO writing! Furthermore, they will only read
   the ASCII parameter files and NOT the larger binary files.  Writing and/or
   binary reading is left to other functions.                                */


XWin2D::XWin2D()
  {
  dname = string("");				// No directory name
  oldMeta = 0;					// Assume new meta format
  Aidx = 1;					// Default acq. directory
  Pidx = 1;					// Default proc. directory
  Acqs   = XWinAcqus();				// Default acqus file
  Procs  = XWinProcs();				// Default procs file
  Acq2s  = XWinAcqu2s();			// Default acqu2s file
  Proc2s = XWinProc2s();			// Default proc2s file
  }

XWin2D::XWin2D(const string& name, int mode, int expno, int procno)
  {
  dname  = name;				// Set Bruker directory name
  oldMeta = 0;					// Assume new meta format
  Aidx   = expno;				// Set the experiment number
  Pidx   = procno;				// Set the processing number 
  SetNames();					// Set up all file names
  bool TF = true;
  if(mode & ios::in)				// If open for reading then
    {
    TF = Acqs.readAPar(Nacqus, 1);		// Read in acqus file
    if(!TF) XWin2Dfatality(20);			// Bail if unsuccessful
    TF = Acq2s.readAPar(Nacqu2s, 1);		// Read in acqu2s file
    if(!TF) XWin2Dfatality(20);			// Bail if unsuccessful
    int N    = Acqs.TD();			// This is the block size
    bool BO = Acqs.BYTORDA()?true:false;		// This is the byte order
    int mode = ios::in|ios::binary; 		// This is the input mode
    int warn = 1;				// This is the warning level
    TF = Ser.open(Nser,N,BO,mode,warn);	// Open the serial file
    if(!TF) XWin2Dfatality(20);			// Bail if unsuccessful
    TF = Procs.readPPar(Nprocs,1);		// Read in the procs file 
    if(!TF) XWin2Dfatality(20);			// Bail if unsuccessful
    TF = Proc2s.readPPar(Nproc2s,1);		// Read in the proc2s file 
    if(!TF) XWin2Dfatality(20);			// Bail if unsuccessful
    }
  else if(mode & ios::out)
    {
    MakeDirs(); 				// If writing, create dirs
    Acqs   = XWinAcqus();			// Set up the acqu(s) pars.
    Procs  = XWinProcs();			// Set up the proc(s) pars.
    Acq2s  = XWinAcqu2s();			// Set up the acqu2(s) pars.
    Proc2s = XWinProc2s();			// Set up the proc2(s) pars.
    }
  }

XWin2D::XWin2D(const XWin2D& XW2D)
  {
  oldMeta = XW2D.oldMeta;
  dname   = XW2D.dname;
  Aidx    = XW2D.Aidx;
  Pidx    = XW2D.Pidx;
  NAIdir  = XW2D.NAIdir;
  NPdir   = XW2D.NPdir;
  NPIdir  = XW2D.NPIdir;
  Nacqu   = XW2D.Nacqu;
  Nacqus  = XW2D.Nacqus;
  Nacqu2  = XW2D.Nacqu2;
  Nacqu2s = XW2D.Nacqu2s;
  Nser    = XW2D.Nser;
  Nproc   = XW2D.Nproc;
  Nprocs  = XW2D.Nprocs;
  Nproc2  = XW2D.Nproc2;
  Nproc2s = XW2D.Nproc2s;
  Nmeta   = XW2D.Nmeta;
  Noutd   = XW2D.Noutd;
  Acqs    = XW2D.Acqs;
  Procs   = XW2D.Procs;
  Ser     = XW2D.Ser;
  Acq2s   = XW2D.Acq2s;
  Proc2s  = XW2D.Proc2s;
  }

XWin2D::~XWin2D() {}

XWin2D& XWin2D::operator=(const XWin2D& XW2D)
  {
  if(this == &XW2D) return *this;
  oldMeta = XW2D.oldMeta;
  dname   = XW2D.dname;
  Aidx    = XW2D.Aidx;
  Pidx    = XW2D.Pidx;
  NAIdir  = XW2D.NAIdir;
  NPdir   = XW2D.NPdir;
  NPIdir  = XW2D.NPIdir;
  Nacqu   = XW2D.Nacqu;
  Nacqus  = XW2D.Nacqus;
  Nacqu2  = XW2D.Nacqu2;
  Nacqu2s = XW2D.Nacqu2s;
  Nser    = XW2D.Nser;
  Nproc   = XW2D.Nproc;
  Nprocs  = XW2D.Nprocs;
  Nproc2  = XW2D.Nproc2;
  Nproc2s = XW2D.Nproc2s;
  Nmeta   = XW2D.Nmeta;
  Noutd   = XW2D.Noutd;
  Acqs    = XW2D.Acqs;
  Procs   = XW2D.Procs;
  Ser     = XW2D.Ser;
  Acq2s   = XW2D.Acq2s;
  Proc2s  = XW2D.Proc2s;
  return *this;
  }


// ____________________________________________________________________________
// B                    XWin2D Parameter Access Functions
// ____________________________________________________________________________
 
/* These functions allow direct access to some of the more important parameters
   that each 2D data set should know.  This class contains objects of many
   types associated with the Bruker 2D data set files. Since they are quite
   similar we do not derive XWin2D from these, since this might cause conflicts
   (that may or may not be sorted out by virtual tables). Rather we rewrite the
   functions herein to match the lower classes and route to the proper
   dimension when necessary.....                                             */
 
/* --------------------------------------------------------------------------  
                          Functions To Retrieve Parameters
   --------------------------------------------------------------------------*/

/*                   NOT REALLY INHERITED FROM CLASS XWinAcqPar
                        (Applies to Files acqu(s) & acqu2(s)                 */

string XWin2D::acqname(int d)      const 	// Acq. Par. File name (short)
   { if(d==1) return Acq2s.acqname(); return Acqs.acqname(); }
double XWin2D::AQ(int d)         const        // Acquisition length (sec)
   { if(d==1) return Acq2s.AQ();      return Acqs.AQ(); }
int    XWin2D::AQ_mod(int d)     const        // Acquisition mode
   { if(d==1) return Acq2s.AQ_mod();  return Acqs.AQ_mod(); }
double XWin2D::BF1(int d)        const        // Base Spectrometer freq.
   { if(d==1) return Acq2s.BF1();     return Acqs.BF1(); }
double XWin2D::BF2(int d)        const        // Base Spectrometer freq.
   { if(d==1) return Acq2s.BF2();     return Acqs.BF2(); }
int    XWin2D::BYTORDA(int d)    const        // Binary byte order
   { if(d==1) return Acq2s.BYTORDA(); return Acqs.BYTORDA(); }
int    XWin2D::DSc(int d)         const        // Number of dummy scans
   { if(d==1) return Acq2s.DSc();      return Acqs.DSc(); }
string XWin2D::EXP(int d)        const        // Experiment name
   { if(d==1) return Acq2s.EXP();     return Acqs.EXP(); }
double XWin2D::XW_IN(int i,int d)   const        // Dwell time
   { if(d==1) return Acq2s.XW_IN(i);     return Acqs.XW_IN(i); }
string XWin2D::NAME(int d)       const        // Full File Name
   { if(d==1) return Acq2s.NAME();    return Acqs.NAME(); }
int    XWin2D::NS(int d)         const        // Number of scans
   { if(d==1) return Acq2s.NS();      return Acqs.NS(); }
string XWin2D::NUC(int i,int d)  const	// Nucleus for a channel
   { if(d==1) return Acq2s.NUC(i);    return Acqs.NUC(i); }
string XWin2D::NUCLEUS(int d)    const        // Base nucleus
   { if(d==1) return Acq2s.NUCLEUS(); return Acqs.NUCLEUS(); }
double XWin2D::O1(int d)         const        // Offset freq.
   { if(d==1) return Acq2s.O1();      return Acqs.O1(); }
double XWin2D::O2(int d)         const        // Offset freq.
   { if(d==1) return Acq2s.O2();      return Acqs.O2(); }
//int    XWin2D::PARMODE(int d=0)    const        // Acquisiiton dimension
//   { if(d==1) return Acq2s.PARMODE(); return Acqs.PARMODE(); }
string XWin2D::PULPROG(int d)    const        // Pulse program
   { if(d==1) return Acq2s.PULPROG(); return Acqs.PULPROG(); }
double XWin2D::SFO1(int d)       const        // Spectrometer freq.
   { if(d==1) return Acq2s.SFO1();    return Acqs.SFO1(); }
double XWin2D::SFO2(int d)       const        // Spectrometer freq.
   { if(d==1) return Acq2s.SFO2();    return Acqs.SFO2(); }
double XWin2D::SFO3(int d)       const        // Spectrometer freq.
   { if(d==1) return Acq2s.SFO3();    return Acqs.SFO3(); }
string XWin2D::SOLVENT(int d)    const        // Solvent
   { if(d==1) return Acq2s.SOLVENT(); return Acqs.SOLVENT(); }
double XWin2D::SW(int d)         const        // Spectral width (PPM)
   { if(d==1) return Acq2s.SW();      return Acqs.SW(); }
double XWin2D::SW_h(int d)       const        // Spectral width (Hz)
   { if(d==1) return Acq2s.SW_h();    return Acqs.SW_h(); }
//int    XWin2D::TD(int d)         const        // Total points
//   { if(d==1) return Acq2s.TD();      return Acqs.TD(); }
double XWin2D::TE(int d)         const        // Sample temperature
   { if(d==1) return Acq2s.TE();      return Acqs.TE(); }

/*                   NOT REALLY INHERITED FROM CLASS XWinSer
                             (Applies to File ser)                           */

string XWin2D::sername() const { return Ser.sername(); }// File name
int    XWin2D::TDS()     const { return Ser.TDS(); }	// No. total points
 
/*                  NOT REALLY INHERITED FROM CLASS XWinProcPar
                           (Applies To Procs and Proc2s)                    */

string XWin2D::parname(int d)   const		// ASCII File name
   { if(d==1) return Proc2s.parname();  return Procs.parname(); }
int    XWin2D::BYTORDP(int d)   const		// Binary byte order
   { if(d==1) return Proc2s.BYTORDP();  return Procs.BYTORDP(); }
int    XWin2D::FT_mod(int d)    const		// How FFT is performed
   { if(d==1) return Proc2s.FT_mod();   return Procs.FT_mod(); }
double XWin2D::LB(int d)        const		// Line Broadening
   { if(d==1) return Proc2s.LB();       return Procs.LB(); }
int    XWin2D::MC2(int d)       const		// FT type on t1
   { if(d==1) return Proc2s.MC2();      return Procs.MC2(); }
double XWin2D::OFFSET(int d)    const		// Spectrum offset
   { if(d==1) return Proc2s.OFFSET();   return Procs.OFFSET(); }
double XWin2D::PHC0(int d)      const		// Zero order phase
   { if(d==1) return Proc2s.PHC0();     return Procs.PHC0(); }
double XWin2D::PHC1(int d)      const		// 1st order phase
   { if(d==1) return Proc2s.PHC1();     return Procs.PHC1(); }
int XWin2D::PH_mod(int d)      const		// Zero order phase
   { if(d==1) return Proc2s.PH_mod();   return Procs.PH_mod(); }
int XWin2D::REVERSE(int d)   const		// Plot spectrum reverse
   { if(d==1) return Proc2s.REVERSE();  return Procs.REVERSE(); }
double XWin2D::SF(int d)        const		// Spectrometer frequency
   { if(d==1) return Proc2s.SF();       return Procs.SF(); }
int    XWin2D::SI(int d)        const		// Data size (re+im)
   { if(d==1) return Proc2s.SI();       return Procs.SI(); }
int    XWin2D::SSB(int d)       const		// Sine bell
   { if(d==1) return Proc2s.SSB();      return Procs.SSB(); }
int    XWin2D::STSI(int d)      const		// Strip size
   { if(d==1) return Proc2s.STSI();     return Procs.STSI(); }
int    XWin2D::STSR(int d)      const		// Strip start
   { if(d==1) return Proc2s.STSR();     return Procs.STSR(); }
double XWin2D::SW_p(int d)      const		// Spectral width (PPM)
   { if(d==1) return Proc2s.SW_p();     return Procs.SW_p(); }
double XWin2D::TDeff(int d)     const		// Effective FFT size
   { if(d==1) return Proc2s.TDeff();    return Procs.TDeff(); }
int    XWin2D::WDW(int d)       const		// Window function
   { if(d==1) return Proc2s.WDW();      return Procs.WDW(); }

/*                   ACCESS FUNCTIONS ON THIS LEVEL ONLY                     */

string XWin2D::name()  const { return dname;  }		// Directory base name 
double XWin2D::Field() const { return Acqs.field(); }	// Bo Field in Teslas
 
/* --------------------------------------------------------------------------  
                          Functions To Specify Parameters
   --------------------------------------------------------------------------*/

/*                NOT REALLY INHERITED FROM CLASS XWinAcqPar
                         (Applies To Acqs and Acq2s)                         */

void XWin2D::AQ_mod(int aqmo, int d)
  { Acq2s.AQ_mod(aqmo); Acqs.AQ_mod(aqmo); }
/*
void XWin2D::BF1(double bf, int d) 			// Not allowed
  {
  if(d==1) { Acq2s.BF1(bf);  Acqs.BF2(bf); }
  else     { Acqs.BF1(bf);   Acq2s.BF2(bf); }
  }
void XWin2D::BF2(double bf, int d)			// Not allowed
  {
  if(d==1) { Acq2s.BF2(bf);  Acqs.BF1(bf); }
  else     { Acqs.BF2(bf);   Acq2s.BF1(bf); }
  }
void XWin2D::BYTORDA(int bo, int d)			// Not allowed
  { Acq2s.BYTORDA(bo); Acqs.BYTORDA(bo); }
*/
int  XWin2D::D(int idx, double tsec, int d, int warn)
  { return Acq2s.D(idx,tsec,warn); return Acqs.D(idx,tsec,warn); }
void XWin2D::DSc(int ds, int d)    
   { if(d==1) Acq2s.DSc(ds); else Acqs.DSc(ds); }
void XWin2D::EXP(const string& exp, int d)
  { Acq2s.EXP(exp); Acqs.EXP(exp); }
void XWin2D::XW_IN(int channel, double in, int d) 
  { Acq2s.XW_IN(channel,in); Acqs.XW_IN(channel,in); }
void XWin2D::O1(double of, int d)   
  {
  if(d==1) { Acq2s.O1(d);  Acqs.O2(d); }
  else     { Acqs.O1(d);   Acq2s.O2(d); }
  }
void XWin2D::O2(double of, int d)
  {
  if(d==1) { Acq2s.O2(d);  Acqs.O1(d); }
  else     { Acqs.O2(d);   Acq2s.O1(d); }
  }
//void XWin2D::NAME(const string& name, int d)   
void XWin2D::NS(int ns, int d)
   { if(d==1) Acq2s.NS(ns); else Acqs.NS(ns); }
//void XWin2D::NUC(int i, const string& N, int d)
void XWin2D::NUCLEI(int channel, const string& I, double O, int warn)
  {
  Acqs.NUCLEI(channel, I, O, warn?1:0); 
  if(channel == 1)
    Acq2s.NUCLEI(channel+1, I, O, warn?1:0); 
  if(channel == 2)
    Acq2s.NUCLEI(channel-1, I, O, warn?1:0); 
  }
//void XWin2D::NUCLEUS(const string& I, int d)
//   { if(d==1) Acq2s.NUCLEUS(I); else Acqs.NUCLEUS(I); }
int  XWin2D::P(int idx, double tp, int d, int warn)
  { return Acq2s.P(idx,tp); return Acqs.P(idx,tp); }
//void XWin2D::PARMODE(int pm, int d)
//  { Acq2s.PARMODE(pm); Acqs.PARMODE(pm); }
void XWin2D::PULPROG(const string& P, int d)
  { if(d==1) Acq2s.PULPROG(P); else Acqs.PULPROG(P); }
void XWin2D::SFO1(double sf, int d)
  {
  if(d==1) { Acq2s.SFO1(sf);  Acqs.SFO2(sf); }
  else     { Acqs.SFO1(sf);   Acq2s.SFO2(sf); Procs.SF(sf); }
  }
void XWin2D::SFO2(double sf, int d)
  {
  if(d==1) { Acq2s.SFO2(sf);  Acqs.SFO1(sf); Proc2s.SF(sf); }
  else     { Acqs.SFO2(sf);   Acq2s.SFO1(sf); }
  }
//void XWin2D::SFO3(double sf, int d)
//void XWin2D::SFO(double sf, int i, int d)
void XWin2D::SOLVENT(const string& S, int d)
  { Acq2s.SOLVENT(S); Acqs.SOLVENT(S); }
void XWin2D::SW(double sw, int d)
  {
  if(d)
    {
    Acq2s.SW(sw);			// Set Sweep width (ppm)
    Acq2s.XW_IN(1, 1.0/Acq2s.SW_h());	// Adjust IN1 time (sec)
    Proc2s.SW_p(Acq2s.SW_h());		// Set processing SW (Hz)
    }
  if(d<1)
    {
    Acqs.SW(sw);			// Set sweep width (ppm)
    Acqs.XW_IN(1, 1.0/Acqs.SW_h());	// Adjust IN1 time (sec)
    Procs.SW_p(Acqs.SW_h()); 		// Set processing SW (Hz)
    }
  }
void XWin2D::SW_h(double sw, int d)
  {
  if(d)
    { 
    Acq2s.SW_h(sw);			// Set Sweep width (Hz)
    Acq2s.XW_IN(1, 1.0/sw);		// Adjust IN1 time (sec)
    Proc2s.SW_p(sw);			// Set processing SW (Hz)
    }
  if(d<1)
    { 
    Acqs.SW_h(sw);			// Set sweep width (Hz)
    Acqs.XW_IN(1, 1.0/sw);		// Adjust IN1 time (sec)
    Procs.SW_p(sw); }			// Set processing SW (Hz)
    }

void XWin2D::TE(double te, int d)
  { Acq2s.TE(te); Acqs.TE(te); }
//void XWin2D::TD(int npts, int d)
//   { if(d==1) Acq2s.TD(npts); else Acqs.TD(npts); }

/*                  NOT REALLY INHERITED FROM CLASS XWinProcPar
                           (Applies To Procs and Proc2s)                    */
 
//void XWin2D::BYTORDP(int bo, int d)        // Set binary byte order
void XWin2D::GB(int gb, int d)
  { if(d) { Proc2s.GB(gb); } if(d<1) { Procs.GB(gb); } }
void XWin2D::FT_mod(int ft, int d)
  { if(d) Proc2s.FT_mod(ft); if(d<1) Procs.FT_mod(ft); }
void XWin2D::FT_mod(const string& ft, int d)
  { if(d==1) { Proc2s.FT_mod(ft); } else { Procs.FT_mod(ft); } }
void XWin2D::LB(int lb, int d)
  { if(d) { Proc2s.LB(lb); } if(d<1){ Procs.LB(lb); } }
void XWin2D::PHC0(double ph0, int d)
  { if(d) { Proc2s.PHC0(ph0); } if(d<1){ Procs.PHC0(ph0); } }
void XWin2D::PHC1(double ph1, int d)
  { if(d) { Proc2s.PHC1(ph1); } if(d<1){ Procs.PHC1(ph1); } }
void XWin2D::PH_mod(const string& p, int d)
  { if(d) { Proc2s.PH_mod(p); } if(d<1){ Procs.PH_mod(p); } }
void XWin2D::MC2(int mc, int d)           { Procs.MC2(mc); Proc2s.MC2(mc); }
void XWin2D::MC2(const string& mc, int d) { Procs.MC2(mc); Proc2s.MC2(mc); }
void XWin2D::REVERSE(int yn, int d)
  { if(d) { Proc2s.REVERSE(yn); } if(d<1) { Procs.REVERSE(yn); } }
//void XWin2D::PPARMOD(int pm, int d)        // Set data dimension
//void XWin2D::SI(int si, int d)             // Set data size (re+im)
//void XWin2D::SF(double SF, int d)            // Set spectrometer freq.
void XWin2D::SSB(int sb, int d)
  { if(d) { Proc2s.SSB(sb); } if(d<1) { Procs.SSB(sb); } }
//void XWin2D::STSI(int sb, int d)          { Procs.x Proc2s.STSI(sb);  }
//void XWin2D::STSR(int sr, int d)          { Procs.x Proc2s.STSR(sr);  }
//void XWin2D::SW_p(double swp, int d)      { Procs.x Proc2s.SW_p(swp); }
void XWin2D::WDW(int wd, int d)
  { if(d) { Proc2s.WDW(wd); } if(d<1) { Procs.WDW(wd); } }
void XWin2D::WDW(const string& wd, int d)
  { if(d) { Proc2s.WDW(wd); } if(d<1) { Procs.WDW(wd); } }

/*                   ACCESS FUNCTIONS ON THIS LEVEL ONLY                     */

void XWin2D::Field(double Bo)
  {
  Acqs.field(fabs(Bo));
  Acq2s.field(fabs(Bo)); 
  }

void XWin2D::Field(double Om, const string& I)
  {
  Isotope X(I);
  double Bo = Om*1.e10*HZ2RAD/X.gamma();
  Acqs.field(fabs(Bo));
  Acq2s.field(fabs(Bo));
  }

void XWin2D::OldMeta(int om) { om?oldMeta=1:oldMeta=0; }

/*
void XWin2D::NUC(int i, int j, string S)
  {
  if(i==0)				// Spectrometer frequency
    {					// in the t2/w2 dimension
    XWinAcqus::NUC(j, S);		//   Set in acqu/acqus
    }
  else
    {					// in the t1/w1 dimension
    Acq2s.NUC(j, S);			//   Set in acqu2/acqu2s
    }
  }

void XWin2D::SFO1(int i, double sf)
  {
  if(i==0)				// Spectrometer frequency
    {					// in the t2/w2 dimension
    XWinAcqus::SFO1(sf);		//   Set in acqu/acqus
    XWinProcs::SF(sf);			//   Set in proc/procs
    }
  else					// Spectrometer frequency
    {					// in the t1/w1 dimension
    Acq2s.SFO1(sf);			//   Set in acqu2/acqu2s
    Proc2s.SF(sf);			//   Set in proc2/proc2s
    }
  }

*/

// ____________________________________________________________________________
// C                 XWinNMR 2D Data Set Input Functions
// ____________________________________________________________________________

/* These functions will read part, or the entire, Bruker XWinNMR 2D Data Set. 

                           __ acqu, acqu2 (changable parameter files)
			  / 
                         /___ acqus, acqu2s (static parameter files)
			/
   expname -- expnum --< ---- ser (binary data)
			\
			 \--- pdata -- 1 -- proc, proc2, procs, proc2s, meta

   Here expname is the base directory containing the data set. The number
   expnum is the experiment number an will default to 1. The acquisition
   parameters are taken from the file acqus and the acquisition binary
   points from ser.  The processed data, if accessed, will come from the
   parameter files procs, meta, etc. and binary data files 1r and 1i.

                       INHERITED FROM CLASS XWinSer

   Function   Arguments              Action                       Return
   --------  ------------  -----------------------------------  ----------
   readFID   fn,TD,bytord  Read 1st block in serial/FID file    row_vector
   readFID   fn,TD,bord,i  Read ith block in serial/FID file    row_vector
 + readFID   idx           Read block of index idx(-1=current)  row_vector
   readFIDs  fn,TD,bytord  Read entire serial file                matrix
 + readFIDs                Read entire serial file                matrix
		    
	     Above +=File Must Be Properly Opened For Reading

 row_vector XWinSer::readFID(const string& fin, int TD, int BO=-1, int idx=-1)
 row_vector XWinSer::readFID(int idx=-1)
 matrix     XWinSer::readFIDs(const string& fin, int TD, int bytord)
 matrix     XWinSer::readFIDs()                                              */

row_vector XWin2D::readFID(const string& dirin, int idx, int ai, int pi, int warn)
  {
  dname = dirin; 
  Aidx = ai;
  Pidx = pi;
  SetNames();						// Set file & dirs
  int    TF  = Acqs.readAPar(Nacqus, warn?1:0);		// Read in acqus file
  if(TF) TF *= Acq2s.readAPar(Nacqu2s, warn?1:0);	// Read in acqu2s file
  if(!TF)
    {							// then output a 
    if(warn)						// warning and bail
      {
      XWin2Derror(21, 1);				// Can't read data set
      if(warn > 1) XWin2Dfatality(21, dname);		// Troubles with dname
      else         XWin2Derror(21, dname, 1);
      }
    }
  return readFID(idx, warn); 
  }

row_vector XWin2D::readFID(int idx, int warn)
  {
  int SBO = Acqs.BYTORDA();				// Get binary byte order
  int STD = Acqs.TD();					// Get block size (re+im)
  int SNB = Acq2s.TD();					// Get # of blocks (re+im)
  if(idx >= SNB)					// Insure range is proper
    {
    XWin2Derror(25, 1);					//   Accessed FID bad
    if(warn > 1) XWin2Dfatality(25, Gdec(SNB-1));	//   FID range
    else         XWin2Derror(25, Gdec(SNB-1));
    }
  row_vector data;					// Array for data
  data = Ser.readFID(Nser,STD,SBO,idx);			// Read FID from ser file
  return row_vector();
  return data;
  }
 
// ____________________________________________________________________________
// D                 XWinNMR 2D Data Set Output Functions
// ____________________________________________________________________________

/* These functions will write an entire Bruker XWinNMR 2D Data Set. 

                           __ acqu, acqu2  (changable parameter files)
			  / 
                         /___ acqus, acqu2s (static parameter file)
			/
     expname -- Aidx --< ---- ser (binary data)
			\
			 \--- pdata -- Pidx -- proc, procs, proc2, proc2s

   Here expname is the base directory specified for output. The numbers Aidx
   and Pidx default to 1.  The produced files are all ASCII except ser. The
   files acqu and acqus will be equivalent as will acqu2-acqu2s, proc-procs,
   and proc2-proc2s.  All files in pdata, since this class will NOT do any
   OUTPUT of frequency domain 2D data, will be defaults based on values in
   acqus and acqu2s. 

   Function  Arguments                         Result
   ________  _________  _______________________________________________________
    write    dir,mx,i     XWinNMR 2D data set produced in dir with data mx 

   Note that GAMMA handles only some of the ASCII parameters one might want
   set in these data sets (# dimensions, # blocks, block size,...). Thus users
   should set additional parameters prior using these functions if values other
   than defaults are desired.  Suggested parameters one might set are:

   Parmeter File  Parameter                    Settings
   _____________  _________  __________________________________________________
    acqu,acqus     AQ_mod    0=qf=one channel, 1=qsim=quadrature simultaneous
			     2=qseq=quadrature sequential, 3=digital quadrature
                   BF1,BF2   Base frequencies in MHz
                   O1,O2     Offset frequencies in Hz
                  SFO{1,2,3} Spectrometer frequencies in MHz
		  SW, SW_h   Spectral Width (PPM, Hz)
									     */

int XWin2D::write(const string& Bdir, const matrix& data, int warn)

  { 
  dname = Bdir;						// Set base directory
  SetNames();						// Set file & dirs
  MakeDirs(); 						// Create directories
  Acqs.PARMODE(1);					// Flag data is 2D
  Acqs.NAME(Nacqus);					// Set full file name
  Acqs.TD(data.cols()*2);				// Set block size re+im
  Procs.PPARMOD(1);					// Flag data is 2D
  Procs.SI(data.cols()*2);				// Set block size re+im
  Acq2s.TD(data.rows());				// Set # of blocks
  Proc2s.SI(data.rows());				// Set # of blocks
  SetConsistent();					// Insure parfiles match
  int TF = Acqs.writeAPar(Nacqus, warn?1:0);		// Write the parameters
  CheckWrite(TF, warn, Nacqus);				// Insure output OK
  TF = Acqs.writeAPar(Nacqu, warn?1:0);			// Write the parameters
  CheckWrite(TF, warn, Nacqu);				// Insure output OK
  TF = Acq2s.writeAPar(Nacqu2s, warn?1:0);		// Write the parameters
  CheckWrite(TF, warn, Nacqu2s);			// Insure output OK
  TF = Acq2s.writeAPar(Nacqu2, warn?1:0);		// Write the parameters
  CheckWrite(TF, warn, Nacqu2);				// Insure output OK
  TF = Ser.write(Nser, data, warn?1:0);			// Write binary data
  CheckWrite(TF, warn, Nser);				// Insure output OK
  XWinMeta XWM;						// XWinNMR meta file
  XWM.OldFlag(oldMeta);					// Set meta format flag
  TF = XWM.write(Nmeta, warn?1:0);			// Output meta file
  CheckWrite(TF, warn, Nmeta);				// Insure output OK
  XWinOutd XWO;						// XWinNMR outd file
  TF = XWO.write(Noutd, warn?1:0);			// Output outd file
  CheckWrite(TF, warn, Noutd);				// Insure output OK
  TF = Procs.writePPar(Nprocs, warn?1:0);		// Output procs file
  CheckWrite(TF, warn, Nprocs);				// Insure output OK
  TF = Procs.writePPar(Nproc, warn?1:0);		// Output proc file
  CheckWrite(TF, warn, Nproc);				// Insure output OK
  TF = Proc2s.writePPar(Nproc2s, warn?1:0);		// Output proc2s file
  CheckWrite(TF, warn, Nproc2s);			// Insure output OK
  TF = Proc2s.writePPar(Nproc2, warn?1:0);		// Output proc2 file
  CheckWrite(TF, warn, Nproc2);				// Insure output OK
  return TF;
  } 


/*
int XWin2D::write(const row_vector& data, int warn)
  { 
// 1. Check That File is Open
// 2. Check That Parameter Files Written
  if(data.size() != size())
    {
    if(warn)
      {
      XWin2Derror(23, 1);				// Cannot write block
      if(warn>1) XWin2Dfatality(24);			// Bad block size
      else       XWin2Derror(24);
      }
    return 0;
    }
  int TF = Ser.write(data, warn?1:0);			// Write the block
  int blks = Acq2s.TD();				// Current blocks
  Acq2s.TD(blks+1);					// Updata block count
  TF = Acq2s.write(Nacqu2s, warn?1:0);			// Update acqu2s file
  CheckWrite(TF, warn, Nacqu2s);			// Insure output OK
  TF = Acq2s.write(Nacqu2, warn?1:0);			// Update acqu2 file
  CheckWrite(TF, warn, Nacqu2);				// Insure output OK
  return TF;
  } 
*/
 
// ____________________________________________________________________________
// E                 XWinNMR 2D Data Set ASCII Output Support
// ____________________________________________________________________________

/* These functions allow for a quick output of the data set contents. They
   don't have anything to do with output while running XWinNMR, rather users
   may just glance at acquisition parameters & values or store then in a 
   small file if needed.  The two function "dpa" and "dpp" mimick the output
   that would be seen in XWinNMR.                                            */
 
ostream& XWin2D::dpa(ostream& ostr, const string& dirin)
  {
  dname = dirin; 					// Copy base directory
  SetNames();						// Set file & dirs
  ReadPars();						// Read parameter files
  return dpa(ostr);					// Display acq pars
  }

ostream& XWin2D::dpa(ostream& ostr) const
  {
  Acqs.dpa(ostr);
  Acq2s.dpa(ostr);
  return ostr;
  }

ostream& XWin2D::dpp(ostream& ostr, const string& dirin)
  {
  dname = dirin; 					// Copy base directory
  SetNames();						// Set file & dirs
  ReadPars();						// Read parameter files
  return dpp(ostr);					// Display proc pars
  }

ostream& XWin2D::dpp(ostream& ostr) const
  {
  Procs.dpp(ostr,BF1(0));				// Display t2/w2 ppars
  Proc2s.dpp(ostr,BF1(1));				// Display t1/w1 ppars
  return ostr;
  }

// ostream& XWin2D:;printPset(ostr) const;		// XWinAcqus Inherited

ostream& XWin2D::print(ostream& ostr, int full) const
  {
  ostr << "\n" << string(21, ' ')
       << "Bruker XWinNMR 2D Data Set\n";
  ostr << "\n\t\tData Set Directory:         ";
  if(!dname.length()) ostr << "Unspecified";
  else                ostr << dname;
  ostr << "\n\t\tASCII Parameter Files:      ";
  if(!dname.length()) ostr << "acqus";
  else                ostr << dname << "/" << Aidx << "/acqus";
  ostr << "\n\t\t                            ";
  if(!dname.length()) ostr << "acqu2s";
  else                ostr << dname << "/" << Aidx << "/acq2us";
  ostr << "\n\t\tBinary Data File:           ";
  if(!dname.length()) ostr << "ser";
  else                ostr << dname + "/" << Aidx << "/ser";
  if(full)
    {
    Acqs.print(ostr, full-1, 0);
    ostr << "\n\t\tNumber of Blocks:           "
	 << Acq2s.TD();
    }
  ostr << "\n";
  return ostr;
  }

ostream& operator<<(ostream& ostr, const XWin2D& XW2D)
  { XW2D.print(ostr); return ostr; }

// ____________________________________________________________________________
// F              XWinNMR 2D Data Set Interactive Functions
// ____________________________________________________________________________

/*

string XWin2D::ask_read(int argc, char* argv[], int& argn, int idx)

	// Input                XW2D : An XWinNMR 2D data set
	//                      argc    : Number of arguments
	//                      argv    : Vector of argc arguments
	//                      argn    : Argument index
	//			idx	: Flag to request expt #
	//			           idx < 1: Ask for it
	//			           idx >=1: Use idx itself def (1)
	// Output               dirname : The parameter argn of array argc
	//                                is used to supply a directory name
	//                                from which the 2D data set resides
	//                                If the argument argn is not in argv,
	//                                the user is asked to supply a name
	//				  and experiment number
	//                                The directory is returned
	// Note                         : The directory should contain an
	//				  XWinNMR 2D data set.

  {
  string dirname;				// Name of base directory
  query_parameter(argc, argv, argn,             // Get directory from command
   "\n\tXWinNMR 2D Directory? ", dirname); 	// Or ask for it
  if(idx < 1)
    {
    argn++;
    query_parameter(argc, argv, argn,		// Get directory index
   "\n\tXWinNMR 2D Sub-Directory Index? ", idx);// Or ask for it
    }
  readFID(dirname, idx);			// Read data set from filename
  return dirname;				// Give back the directory
  }
*/



#endif							// XWin2D.cc
