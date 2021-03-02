/* XWin1D.cc ****************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A                           **
**                                                                      **
**      XWin1D				     Implementation		**
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
** sets. This class embodies a Bruker data set for a 1D acquisition.    **
** Each 1D acquisition generates several files (binary and ASCII) 	**
** spanning several directories.  This class is intended to handle 	**
** this structure.  Below is a typical Bruker directory structure 	**
** associated with a 1D acquisiton.					**
**									**
**                          __ acqu  (changable parameter file)		**
**			   / 						**
**                        /___ acqus (static parameter file)		**
**			 /						**
**  expname -- expnum --< ---- fid (binary data)			**
**			 \						**
**			  \___ pdata -- 1 -- proc, procs, meta		**
**									**
** This class will handle the directory hierarchy shown above as	**
** well as most of the files therein. Thus, given a base directory 	**
** for an experiment/simulation output (expname) & an experiment	**
** number (expnum) this class will generate the subdirectories shown,	**
** the acquisition parameter (acqu, acqus) and binary (fid) files	**
** as well as the processing files (proc, procs,meta).  It is NOT 	**
** intended to generated processed 1D data files (1i, 1r) - that 	**
** should be done using the XWinNMR software.  The class will also	**
** allow for import of 1D data sets.  The fid file and a processed	**
** spectrum can be read and placed into a GAMMA row_vector. Associated	**
** parameters are also provided for.					**
**									**
*************************************************************************/

#ifndef   XWin1D_cc_			// This file already included?
#  define XWin1D_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler
#    pragma implementation		// then this is the implementation
#  endif

#include <string>			// Include libstdc++ strings
#include <iostream>			// Include input output streams
#include <stdlib.h>
#include <GamIO/XWin1D.h>		// Include interface
#include <GamIO/XWinAcqus.h>		// Include acqus file interface
#include <GamIO/XWinFid.h>		// Include fid file interface
#include <GamIO/XWinMeta.h>		// Include meta file interface
#include <GamIO/XWinOutd.h>		// Include outd file interface
#include <GamIO/XWinProcs.h>		// Include procs file interface
#include <GamIO/XWinSpec.h>		// Include spec files interface
#include <GamIO/BinIOBase.h>		// Include binary I/O functions
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <Basics/StringCut.h>		// Include GAMMA Gdec function
#include <sys/stat.h>			// Include POSIX mkdir function
#if defined(_MSC_VER)                   // If not using MSVC++ then we                              // and if using MSVC++ then I guess
 #include <direct.h>			// Include mkdir function
#endif

using std::string;			// Using libstdc++ strings
using std::ostream;			// Using libstdc++ output streams


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
// i                   XWinNMR 1D Data Set Error Handling
// ____________________________________________________________________________
 
/* These functions take care of any errors encountered when reading, writing,
   and setting parameters in Bruker acquisition data sets.
 
        Input           eidx    : Error index
                        pname   : Additional error message
                        noret   : Flag for linefeed (0=linefeed)
        Output          void    : An error message is output                 */  

void XWin1D::XWin1Derror(int eidx, int noret) const                       
  {
  string hdr("XWinNMR 1D Data Set");
  string msg;
  switch (eidx)
    {
    case 21: GAMMAerror(hdr,"Cannot Read Full Data Set",     noret); break; // (21)
    case 22: GAMMAerror(hdr,"Cannot Generate Full Data Set", noret); break; // (22)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }  


void XWin1D::XWin1Derror(int eidx, const string& pname, int noret) const
  {
  string hdr("XWinNMR 1D Data Set");
  string msg;
  switch(eidx)
    {
    case 5: msg = string("Bad Use Of ") + pname + string(" Function ");
             GAMMAerror(hdr,msg+pname,noret);  break;                  // (5)
    case 21:msg = string("Problems With Files in Directory ");
             GAMMAerror(hdr,msg+pname,noret);  break;                  // (21)
    case 22:msg = string("Problems With Use of Directory ");
             GAMMAerror(hdr,msg+pname,noret);  break;                  // (22)
    default: GAMMAerror(hdr, eidx, pname, noret); break;// Unknown Error  (-1)
    }
  }  
 

volatile void XWin1D::XWin1Dfatality(int eidx) const
  {
  XWin1Derror(eidx, 1);			// Normal non-fatal error
  if(eidx) XWin1Derror(0);		// Program aborting error
  GAMMAfatal();					// Clean exit from program			// Quit program
  }
 

volatile void XWin1D::XWin1Dfatality(int eidx, const string& pname) const
  {
  XWin1Derror(eidx, pname, 1);		// Normal non-fatal error
  if(eidx) XWin1Derror(0);		// Program aborting error
  GAMMAfatal();					// Clean exit from program			// Quit program
  }

// ____________________________________________________________________________ 
// ii                   XWinNMR 1D Data Set Error Handling
// ____________________________________________________________________________
 
int XWin1D::CheckDir(int TF,   int warn, const string& dout) const
  {
  if(!TF) return 1;				// Good mkdir returns 0
  string cmd = "cd " + dout;			// Perhaps file is from
  if(!system(cmd.c_str())) return 1;		// pre-exisiting dir.
  if(warn)
    {
    XWin1Derror(22, dout, 1);			// Problems with directory
    if(warn > 1) XWin1Dfatality(22);		// Can't generate data set
    else         XWin1Derror(22,1);
    }
  return TF;
  }


int XWin1D::CheckWrite(int TF, int warn, const string& dout) const
  {
  if(!TF)
    {
    if(warn)
      {
      XWin1Derror(21, dout, 1);		// Problems with files
      if(warn > 1) XWin1Dfatality(22);	// Can't generate data set
      else         XWin1Derror(22,1);
      }
    }
  return TF;
  }

// ____________________________________________________________________________
// iii                   XWinNMR 1D Data Set Setup Functions
// ____________________________________________________________________________
 
/* These functions are used to quickly set up values associated with the 1D
   data set.  In particular they handle setting up a standard Bruker directory
   structure and the file names in the data set.  Another function is used to
   insure the that the parameter sets are self-consistent. This should be
   called before the data set (or at least the parameter files) are output.  */

void XWin1D::SetNames()
  {
  NAIdir  = dname + "/" + Gdec(Aidx);           // Where to find acq. data
  Nacqu   = NAIdir + "/acqu";                   // Parameter file (acqu)
  Nacqus  = NAIdir + "/acqus";                  // Parameter file (acqus)
  Nfid    = NAIdir + "/fid";                    // Data      file (fid)
  NPdir   = NAIdir + "/pdata";                  // Base processing dir
  NPIdir  = NPdir  + "/" + Gdec(Pidx);          // Where to find proc. data
  Nproc   = NPIdir + "/proc";                   // Parameter file (proc)
  Nprocs  = NPIdir + "/procs";                  // Parameter file (procs)
  Nmeta   = NPIdir + "/meta";                   // Parameter file (meta)
  Noutd   = NPIdir + "/outd";                   // Parameter file (outd)
  Nspec   = NPIdir + "/1"; 			// Base processed data name
  }
 
// sosi - To placate MSVC++'s different definition of mkdir, I have made the
//        function MakeADir (in BinIOBase.cc) that will call either Microsofts
//        or GCC's version.

int XWin1D::MakeDirs(int warn)
  {
  int TF = 1;
//  TF *= mkdir(dname.c_str(),  493);             // Try & make base directory
  TF *= MakeADir(dname,493);
  CheckDir(TF, warn, dname);                    // Insure directory OK
//  TF *= mkdir(NAIdir.c_str(), 493);             // Try & make acq directory
  TF *= MakeADir(NAIdir, 493);
  CheckDir(TF, warn, NAIdir);                   // Insure directory OK
//  TF *= mkdir(NPdir.c_str(),  493);             // Try & make process dir.
  TF *= MakeADir(NPdir, 493);
  CheckDir(TF, warn, NPdir);                    // Insure directory OK
//  TF *= mkdir(NPIdir.c_str(), 493);             // Try & make process dir.
  TF *= MakeADir(NPIdir,493);
  CheckDir(TF, warn, NPIdir);                   // Insure directory OK
  return TF;
  }

/*
int XWin1D::ReadPars(int warn)
  {
  int TF = 1;
// sosi
  return TF;
  }
*/

void XWin1D::SetConsistent()
  {
  BYTORDP(BYTORDA());			// Set matching byte order
  SF(SFO1());				// Set matching frequecies
  SW_p(SW());				// Set matching spectral widths
  }


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A               XWinNMR 1D Data Set Constructors, Destructor
// ____________________________________________________________________________


XWin1D::XWin1D() : XWinAcqus(), XWinFid(), XWinProcs(), XWinSpec()
  { dname = string(""); Aidx = 1; Pidx = 1;}

XWin1D::~XWin1D() {}

// ____________________________________________________________________________
// B                    XWinNMR Fid File Access Functions
// ____________________________________________________________________________

/*                      Functions To Obtain Parameters

   These functions allow users to get some simple information regarding the
   contents of the Bruker data acqusition file.

                        INHERITED FROM CLASS XWinAcqus
                        ==============================
 
int    XWin1D::acqname()   const { return _parfile;          }
double XWin1D::AQ()        const { return _TD/(2*_SW*_SFO1); }
int    XWin1D::AQ_mod()    const { return _AQ_mod;           }
double XWin1D::BF1()       const { return _BF[0];            }
double XWin1D::BF2()       const { return _BF[1];            }
int    XWin1D::BYTORDA()   const { return _BYTORDA           }
int    XWin1D::DS()        const { return _DS;               }
string XWin1D::EXP()       const { return _EXP;              }
double XWin1D::IN(int i)   const { return _IN[i];            }
string XWin1D::NAME()      const { return _NAME;             }
int    XWin1D::NS()        const { return _NS;               }
string XWin1D::NUC(int i)  const { return _NUC[i];           }
string XWin1D::NUCLEUS()   const { return _NUCLEUS;          }
double XWin1D::O1()        const { return _O[0];             }
double XWin1D::O2()        const { return _O[1];             }
int    XWin1D::PARMODE()   const { return _PARMODE;          }
string XWin1D::PULPROG()   const { return _PULPROG;          }
double XWin1D::SFO1()      const { return _SFO[0];           }
double XWin1D::SFO2()      const { return _SFO[1];           }
double XWin1D::SFO3()      const { return _SFO[2];           }
string XWin1D::SOLVENT()   const { return _SOLVENT;          }
double XWin1D::SW()        const { return _SW;               }
double XWin1D::SW_h()      const { return _SW/_SFO[0];       }
double XWin1D::TE()        const { return _TE;               }
int    XWin1D::TD()        const { return _TD;               }

                        INHERITED FROM CLASS XWinFid
                        ============================

string     XWinFid::fidname()  const { return ffname; }
int        XWinFid::size()     const { return ftotpts/2; }
int        XWinFid::TDF()      const { return ftotpts; }
int        XWinFid::bytes()    const { return ffsize; }
int        XWinFid::pad()      const { return fpadding; }
bool       XWinFid::order()    const { return fbyteordin; }
row_vector XWinFid::data()     const { return fdata; }

                        INHERITED FROM CLASS XWinSpec
                        =============================
 
string     XWinSpec::specname()  const { return sfname; }
string     XWinSpec::specrname() const { return sfname + "r"; }
string     XWinSpec::speciname() const { return sfname + "i"; }
int        XWinSpec::size()      const { return stotpts; }
int        XWinSpec::SSI()       const { return stotpts; }
int        XWinSpec::bytes()     const { return sfsize; }
int        XWinSpec::pad()       const { return spadding; }
bool       XWinSpec::order()     const { return sbyteordin; }
row_vector XWinSpec::data()      const { return sdata; }


                       INHERITED FROM CLASS XWinProcs
                       ==============================

string XWinProcs::parname()   const { return parfile;  } // File name
int    XWinProcs::BYTORDP()   const { return _BYTORDP; } // Byte order
int    XWinProcs::FT_mod()    const { return _FT_mod;  } // FFT mode
double XWinProcs::LB()        const { return _LB;      } // Line Broadening
int    XWinProcs::MC2()       const { return _MC2;     } //
double XWinProcs::OFFSET()    const { return _OFFSET;  } // Spectrum offset
double XWinProcs::PHC0()      const { return _PHC0;    } // Zero order phase
double XWinProcs::PHC1()      const { return _PHC1;    } // 1st order phase
string XWinProcs::REVERSE()   const { return _REVERSE; } // Spectrum reverse
double XWinProcs::SF()        const { return _SF;      } // Spectrom. freq.
int    XWinProcs::SI()        const { return _SI;      } // Size (re+im)
int    XWinProcs::SSB()       const { return _SSB;     } // Sine bell
int    XWinProcs::STSI()      const { return _STSI;    } // Sine bell
int    XWinProcs::STSR()      const { return _STSR;    } // Sine bell
double XWinProcs::SW_p()      const { return _SW_p;    } // Spec. width(ppm)
double XWinProcs::TDeff()     const { return _TDeff;   } // Effective FID size
int    XWinProcs::WDW()       const { return _WDW;     } // Window function


                        INHERITED FROM CLASS XWinAcqus
                        ==============================
 
void XWin1D::BYTORDP(int bo            { bo?_BYTORDP=1:_BYTORDP=0;}
void XWin1D::AQ_mod(int aqmo)          { _AQ_mod  = aqmo;        }
int  XWin1D::D(int i, double t, int w) { return SetDelay(i,t,w); }
void XWin1D::DS(int ds)                { _DS      = ds;          }
void XWin1D::EXP(const string& exp)    { _EXP     = exp;         }
void XWin1D::NAME(const string& nm)    { _NAME    = fulldir+nm;  }
void XWin1D::NS(int ns)                { _NS      = ns;          }
void XWin1D::NUCLEUS(const string& I)  { _NUCLEUS = I;           }
void XWin1D::NUCLEI(int CH, const string& I, double OF, int warn)
int  XWin1D::P(int i, double t, int w) { return SetPulse(i,t,w); }
void XWin1D::PULPROG(const string& P)  { _PULPROG = P;           }
void XWin1D::SFO1(double sf)           { _SFO[0]  = sf;          }
void XWin1D::SFO2(double sf)           { _SFO[1]  = sf;          }
void XWin1D::SFO3(double sf)           { _SFO[2]  = sf;          }
void XWin1D::SFO(double sf, int i)     { _SFO[i]  = sf;          }
void XWin1D::SOLVENT(const string& S)  { _SOLVENT = S;           }
void XWin1D::SW(double sw)             { _SW      = sw;          }
void XWin1D::SW_h(double sw)           { _SW      = sw/_SFO[0];  }
void XWin1D::TD(int td)                { _TD      = td;          }
void XWin1D::TE(double te)             { _TE      = te;          }

                       INHERITED FROM CLASS XWinProcs
                       ==============================
 
void XWinProcs::MC2(int mc)             { _MC2     = mc; }
void XWinProcs::REVERSE(int yn)         { yn?_REVERSE="yes":_REVERSE="no"; }
void XWinProcs::PPARMOD(int pm)         { _PPARMOD = pm; }
void XWinProcs::SI(int si)              { _SI      = si; }
void XWinProcs::SF(double SF)           { _SF      = SF; SetOffset(); }
void XWinProcs::SSB(int sb)             { _SSB     = sb; }
void XWinProcs::STSI(int st)            { _STSI    = st; }
void XWinProcs::STSR(int sr)            { _STSR    = sr; }
void XWinProcs::SW_p(double swp)        { _SW_p    = swp; SetOffset(); }
void XWinProcs::WDW(int wd)             { _WDW     = wd; }                   */  

string     XWin1D::name()     const { return dname; } 
row_vector XWin1D::FID()      const { return XWinFid::data(); }
row_vector XWin1D::Spectrum() const { return XWinSpec::data(); }

// ____________________________________________________________________________
// C                 XWinNMR 1D Data Set Input Functions
// ____________________________________________________________________________

/* These functions will read an entire Bruker XWinNMR 1D Data Set. 

                           __ acqu  (changable parameter file)
			  / 
                         /___ acqus (static parameter file)
			/
   expname -- expnum --< ---- fid (binary data)
			\
			 \--- pdata -- 1 -- proc, procs, 1r, 1i, meta, outd

   Here expname is the base directory containing the data set. The number
   expnum is the experiment number an will defaul to 1. The acquisition
   parameters are taken from the file acqus and the acquisition binary
   points from fid.  The processed data, if accessed, will come from the
   parameter files procs, meta, etc. and binary data files 1r and 1i.

   Function                               Purpose
  ----------        --------------------------------------------------------

    readFID         Read in data set FID for class object.  This is done
		    special for both parameters in ASCII files and for
		    binary data. Since Bruker parameter format differs
		    from GAMMA parameter format all their parameter files
		    must be parsed appropriately. Similarly, their binary
		    data format is specific to XWinNMR and must be read in
		    accordingly. Parameters that relate to processed 1D data
		    are ignored (i.e. we only deal with acqus & fid).
  parsePSet         Converts parameters in internal pset to specific values
   getPar           Returns parameter found in the parameter set herein 
                    This function is inherited from base class XWinPSet.     */

int XWin1D::read(const string& dirin, int aidx, int pidx, int warn)
  {
  dname = dirin;
  Aidx = aidx;
  Pidx = pidx;
  return read(warn);
  }

int XWin1D::read(int wn)
  {
  SetNames();						// Set all filenames
  bool   TF = XWinAcqus::readAPar(Nacqus,       wn?1:0);// Read in acqus file
  bool bytorda = (BYTORDA())?true:false;  		// FID data byte order
  if(TF) TF = XWinFid::read(Nfid,  bytorda,TD(),wn?1:0);// Read in fid if OK
  if(TF) TF = XWinProcs::readPPar(Nprocs,       wn?1:0);// Read in procs file
  bool bytordp = (BYTORDP())?true:false;		// Spec data byte order
  if(TF) TF = XWinSpec::read(Nspec,bytordp,SI(),wn?1:0);// Read in spec if OK
  if(!TF)						// If we can't read the
    {							// acqus,procs, or spec 
    if(wn)						// files, output warning
      { 						// and bail if desired
      XWin1Derror(21, 1);				// Can't read data set
      if(wn > 1) XWin1Dfatality(21, dname);		// Troubles with dname
      else       XWin1Derror(21, dname, 1);
      return -1;
      }
    return -1;
    }
  return (TF)?1:0;
  }

int XWin1D::readFID(const string& dirin, int idx, int warn)
  { dname = dirin; Aidx = idx; Pidx = 1; return readFID(warn); }

int XWin1D::readFID(int warn)
  {
  SetNames();			// Set all filenames
  bool byteord = (BYTORDA())?true:false;	
  bool   TF=XWinAcqus::readAPar(Nacqus, warn?1:0);	// Read in acqus file
  if(TF) TF=XWinFid::read(Nfid,byteord,TD(),warn?1:0);// Read in fid if OK
  if(!TF)						// If we can't read the
    {							// acqus or fid file
    if(warn)						// output warning & bail
      {
      XWin1Derror(21, 1);				// Can't read data set
      if(warn > 1) XWin1Dfatality(21, dname);		// Troubles with dname
      else         XWin1Derror(21, dname, 1);
      return -1;
      }
    return -1;
    }
  return (TF)?1:0;
  }

int XWin1D::readSpectrum(const string& dirin, int aidx, int pidx, int warn)
  { dname = dirin; Aidx = aidx; Pidx = pidx; return readSpectrum(warn); }

int XWin1D::readSpectrum(int wn)
  {
  SetNames();						// Set all filenames
  bool bytord = (BYTORDP())?true:false;	
  bool   TF=XWinAcqus::readAPar(Nacqus, wn?1:0);	// Read in acqus file
  if(TF) TF=XWinProcs::readPPar(Nprocs, wn?1:0);	// Read in procs file
  if(TF) TF=XWinSpec::read(Nspec,bytord,SI(),wn?1:0);// Read in spec if OK
  if(!TF)						// If we can't read the
    {							// acqus,procs, or spec 
    if(wn)						// files, output warning
      { 						// and bail if desired
      XWin1Derror(21, 1);				// Can't read data set
      if(wn > 1) XWin1Dfatality(21, dname);		// Troubles with dname
      else       XWin1Derror(21, dname, 1);
      return -1;
      }
    return -1;
    }
  return (TF)?1:0;
  }

// ____________________________________________________________________________
// D                 XWinNMR 1D Data Set Output Functions
// ____________________________________________________________________________

/* These functions will write an entire Bruker XWinNMR 1D Data Set. 

                           __ acqu  (changable parameter file)
			  / 
                         /___ acqus (static parameter file)
			/
   expname -- expnum --< ---- fid (binary data)
			\
			 \___ pdata -- 1 -- proc, procs, meta

   Here expname is that specified output.  The number expnum
   will default to 1.  The produced files are all ASCII except fid. The
   files acqu and acqus will be equivalent as will proc and procs.  All
   files in pdata, since this class will NOT do any OUTPUT of frequency
   domain 1D data, will be defaults based on values in acqus. 

   Function                               Purpose
 ____________________________________________________________________________

    read            Read in parameter set for class object.  This is done */


//       Input           XW1D	: XWinNMR 1D data set
//                       Bdir	: Base directory of 1D data set
//                       data	: 1D FID data vector
//			 warn	: Flag for warning levels
//       Output          XW1D	: An 1D data set is written to Bdir
//			          This includes Bdir/idx/{acqu, acqus, fid}
//				  & Bdir/idx/pdata/ix2/{proc,procs,meta,outd}
//	 Note			: The user MUST specify several variables
//				  that are essential to the data set prior

int XWin1D::write(const string& Bdir, const row_vector& data, int warn)
  { 
  dname = Bdir;						// Set base directory
  return write(data, warn);				// Use overload
  }

int XWin1D::write(const row_vector& data, int warn)
  {
  SetNames();						// Set file & dirs
  MakeDirs();						// Create directories
  PARMODE(0);						// Flag data is 1D
  TD(data.cols()*2);					// Set block size
  int TF = XWinAcqus::writeAPar(Nacqus, warn?1:0);	// Write params  (acqus)
  CheckWrite(TF, warn, Nacqus);                         // Insure output OK
  TF = XWinAcqus::writeAPar(Nacqu, warn?1:0);  		// Write params  (acqu)
  CheckWrite(TF, warn, Nacqu);                          // Insure output OK
  TF = XWinFid::write(Nfid, data, warn?1:0);		// Write bin. data (fid)
  CheckWrite(TF, warn, Nfid);                           // Insure output OK
  XWinOutd XWO;                                         // XWinNMR outd file
  TF = XWO.write(Noutd, warn?1:0);                      // Output outd file  
  CheckWrite(TF, warn, Noutd);                          // Insure output OK
  XWinMeta XWM;						// XWinNMR meta file
  TF = XWM.write(Nmeta, warn?1:0);			// Output meta file
  CheckWrite(TF, warn, Nmeta);				// Insure output OK
  XWinProcs XWPS;					// XWinNMR procs file
  TF = XWPS.writePPar(Nprocs, warn?1:0);		// Output procs file
  CheckWrite(TF, warn, Nprocs);				// Insure output OK
  TF = XWPS.writePPar(Nproc, warn?1:0);			// Output proc file
  CheckWrite(TF, warn, Nproc);				// Insure output OK
  return TF;
  } 

/*
int XWin1D::write(const row_vector& data, int warn=2)
  { 
  string fout = dname + string("/acqus");		// ASCII parameter file
  int TF = XWinAcqus::write(fout, data, warn?1:0);	// Write the parameters
  fout = dname + string("/fid");			// Binary data file
  TF = XWinFid::write(fout, data, warn?1:0);		// Write binary data
  } 
*/
 
// ____________________________________________________________________________
// E                 XWinNMR 1D Data Set ASCII Output Support
// ____________________________________________________________________________

/* These functions allow for a quick output of the data set contents. They
   don't have anything to do with output while running XWinNMR, rather users
   may just glance at fid parameters & values or store then in a small file. */
 
// ostream& XWin1D:;printPset(ostr) const;		// XWinAcqus Inherited

ostream& XWin1D::print(ostream& ostr, int full) const
  {
  ostr << "\n" << string(21, ' ')
       << "Bruker XWinNMR Single 1D Data Set\n";
  ostr << "\n\t\tData Set Directory:          ";
  if(!dname.length()) ostr << "Unspecified";
  else                ostr << dname;

  ostr << "\n\t\tAcquisition Parameter File:  ";
  if(!dname.length()) ostr << "acqus";
  else                ostr << Nacqus;
  ostr << "\n\t\tBinary FID Data File:        ";
  if(!dname.length()) ostr << "fid";
  else                ostr << Nfid;
  ostr << "\n\t\tProcessing Parameter File:   ";
  if(!dname.length()) ostr << "procs";
  else                ostr << Nprocs;
  ostr << "\n\t\tBinary Processed Data Files: ";
  if(!dname.length()) ostr << "1r, 1i";
  else                ostr << Nspec << "{r,i}";
  if(full)
    {
    ostr << "\n\t  Acquisition Info";
    XWinAcqus::print(ostr, full-1, 0);
    ostr << "\n\t  Processing Info";
    XWinProcs::print(ostr, full-1, 0);
    }
  ostr << "\n";
  return ostr;
  }

ostream& operator<<(ostream& ostr, const XWin1D& XW1D)
  { XW1D.print(ostr); return ostr; }

// ____________________________________________________________________________
// F              XWinNMR 1D Data Set Interactive Functions
// ____________________________________________________________________________


string XWin1D::ask_read(int argc, char* argv[], int& argn, int aidx, int pidx)

	// Input                XW1D : An XWinNMR 1D data set
	//                      argc    : Number of arguments
	//                      argv    : Vector of argc arguments
	//                      argn    : Argument index
	//			aidx	: Flag to request experiment #
	//			           aidx < 1: Ask for it
	//			           aidx >=1: Use aidx itself def (1)
	//			pidx	: Flag to request processing #
	//			           pidx < 1: Ask for it
	//			           pidx >=1: Use pidx itself def (1)
	// Output               dirname : The parameter argn of array argc
	//                                is used to supply a directory name
	//                                from which the 1d data set resides
	//                                If the argument argn is not in argv,
	//                                the user is asked to supply a name
	//				  and experiment number
	//                                The directory is returned
	// Note                         : The directory should contain an
	//				  XWinNMR 1D data set.

  {
  string dirname;				// Name of base directory
  query_parameter(argc, argv, argn,             // Get directory from command
       "\n\tXWinNMR 1D Directory? ", dirname); 	// Or ask for it
  if(aidx<1 || aidx>5000)
    {
    argn++;
    query_parameter(argc, argv, argn,
       "\n\tXWinNMR 1D Experiment SubDirectory Index? ", aidx);
    }
  if(pidx<1 || pidx>5000)
    {
    argn++;
    query_parameter(argc, argv, argn,
       "\n\tXWinNMR 1D Processing SubDirectory Index? ", pidx);
    }
  read(dirname, aidx, pidx);			// Read data set from filename
  return dirname;				// Give back the directory
  }



#endif							// XWin1D.cc
