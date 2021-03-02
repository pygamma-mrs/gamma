/* XWinProcPar.cc ***********************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      XWinProcPar                                  Implementation	**
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
** sets. This class embodies a Bruker parameter file, procs, that which **
** to control access to NMR processing parameters within XWinNMR. This	**
** is another ASCII file that has nothing to do with GAMMA and will	**
** only be used	for output. The file is usually dir/1/pdata/1/procs	**
** where dir is the directory base that contains XWinNMR data.  Note	**
** that without	the existence of a proc file XWinNMR will complain when	**
** looking at the data, even if it is just the associated fid.  One may	**
** still be able to plot the fid, but then one cannot display status 	**
** parameters under the "Output" menu in XWinNMR.  Below is a typical	**
** Bruker directory structure in a 1D acquisiton case.			**
**									**
**                          __ Acqu  (changable parameter file)		**
**			   / 						**
**                        /___ acqus (static parameter file)		**
**			 /						**
**  expname -- expnum --< ---- fid (binary data)			**
**			 \						**
**			  \___ pdata -- 1 -- proc, procs, Procs		**
**									**
*************************************************************************/

#ifndef   XWinProcPar_CC_		// Is file already included?
#  define XWinProcPar_CC_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <GamIO/XWinProcPar.h>		// Include the interface
#include <GamIO/BinIOBase.h>            // Include binary IO functions
#include <GamIO/XWinPSet.h>		// Include Bruker parameter sets
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <Basics/StringCut.h>		// Include GAMMA string cutting
#include <string>			// Include libstdc++ strings
#ifndef _MSC_VER                        // If not using MSVC++ then we
 #include <sys/time.h>			// Include time and date access
 #include <unistd.h>			// Include POSIX getcwd function
#else                                   // and if using MSVC++ then I guess
 #include <time.h>			// Include time and date access
 #include <direct.h>			// Include getcwd function
#endif

using std::string;			// Using libstdc++ strings
using std::ofstream;			// Using libstdc++ output file streams
using std::ostream;			// Using libstdc++ output streams

// ----------------------------------------------------------------------------
// ----------------------------- PRIVATE FUNCTIONS ----------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                      XWinNMR Procs File Error Handling
// ____________________________________________________________________________
 
/* These functions take care of any errors encountered when reading, writing,
   and setting parameters in Bruker processing parameter files.

        Input		eidx    : Error index
        		pname   : Additional error message
        		noret	: Flag for linefeed (0=linefeed)
        Output		void    : An error message is output                 */
 
void XWinProcPar::XWPPerror(int eidx, int noret) const
  {
  string hdr("XWinNMR Process. Param. File");
  string msg;
  switch (eidx)
    {
    case 25: GAMMAerror(hdr,"Cannot Write Parameters To File",noret);break; // (25)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }  
    
void XWinProcPar::XWPPerror(int eidx, const string& pname, int noret) const
  {                                                                             
  string hdr("XWinNMR Process. Param. File");
  string msg;
  switch(eidx)
    {
    case 21: msg = string("Cannot Write To ") + pname;
             GAMMAerror(hdr,msg+pname,noret);  break;			// (21)
    case 22: msg = string("Cannot Set ") + pname + string(" Value");
             GAMMAerror(hdr,msg,noret);  break;				// (22)
    case 23: msg = string("Valid Settings Are ") + pname;
             GAMMAerror(hdr,msg,noret);  break;				// (23)
    case 24: msg = string("Setting Parameter Default To ") + pname;
             GAMMAerror(hdr,msg,noret);  break;  			// (24)
    default: GAMMAerror(hdr, eidx, pname, noret); break;// Unknown Error  (-1)
    }
  }

volatile void XWinProcPar::XWPPfatal(int eidx) const
  {                                                                 
  XWPPerror(eidx, 1);			// Normal non-fatal error
  if(eidx) XWPPerror(0);		// Program aborting error
  GAMMAfatal();				// Clean exit from program
  }

volatile void XWinProcPar::XWPPfatal(int eidx, const string& pname) const
  {                                                                 
  XWPPerror(eidx, pname, 1);		// Normal non-fatal error
  if(eidx) XWPPerror(0);		// Program aborting error
  GAMMAfatal();				// Clean exit from program
  }
 
// ____________________________________________________________________________
// ii                 XWinNMR Procs File Default Parmameters
// ____________________________________________________________________________

void XWinProcPar::SetDefaults(const string& fname)
  {
  parfile   = fname;
  _TITLE    = "Parameter File,  Version xwinnmr1.1  GAMMA";
  _JCAMPDX  = 5.0;
  _DATATYPE = "Parameter Values";
  _ORIGIN   = "UXNMR, Bruker Analytische Messtechnik GmbH";
  _OWNER    = "GAMMA";
  time_t longtime;              // Need a time structure
  longtime = NULL;	//

#ifdef _MSC_VER
  struct tm newtime;		// For setting current date
  localtime_s(&newtime, &longtime);
	char date[31];
	asctime_s(date,30,&newtime);
  _DATE = date;
#else
  struct tm *ptr;		// For setting current date
  ptr = localtime(&longtime);
  _DATE = asctime(ptr);	
#endif

  _ABSF1   = 0;
  _ABSF2   = 0;
  _ABSG    = 0;
  _ABSL    = 0;
  _XALPHA   = 0;
  _AQORDER = 0;
  _ASSFAC  = 0;
  _ASSFACI = 0;
  _ASSFACX = 0;
  _ASSWID  = 0;
  _AUNMP   = "<>";		// For AU
  _AZFE    = 0;			// Intetral ends
  _AZFW    = 0;			// Multiplet integral
  _BCFW    = 0;			// Unknown
  _BC_mod  = 0;			// BC modify
  _BYTORDP = WeRBigEnd();	// Byte ordering
  _COROFFS = 0;			// Unknown
  _DATMOD  = 1;			// Add spectra
  _DC      = 2;			// For spectra add
  _DFILT   = "<>";		// Filter file
  _DTYPP   = 0;			// Unknown
  _FCOR    = 0.0;		// 1st pt multiply
  _FTSIZE  = 0;			// Unknown
  _FT_mod  = 0;			// FT aquire mod
  _GAMMA   = 1;			// Quad correct. 
  _GB      = 0;			// Gauss. multiply
  _INTBC   = 1;			// Integration param.
  _INTSCL  = 1;			// Integral compare
  _ISEN    = 128;		// Small integral
  _LB      = 0.3;		// Line broadening
  _LEV0    = 0;			// Contouring param
  _LPBIN   = 0;			// Linear prediction
  _MAXI    = 1000;		// Peak picking
  _MC2     = 0;			// 2D tranform type
  _MEAN    = 0;			// Unknown
  _ME_mod  = 0;			// Linear prediciton
  _MI      = 0;			// Peak picking
  _NCOEF   = 0;			// LP coefficients
  _NC_proc = 0;			// 2D processing
  _NLEV    = 6;			// Contour levels
  _NOISF1  = 0;			// SN ratio
  _NOISF2  = 0;			// SN ratio
  _NSP     = 0;			// Left shift
  _NTH_PI  = 0;			// Unknown
  _NZP     = 0;			// For zp
  _OFFSET  = 2000.0;		// 1st pt shift (Hz)
  __PC      = 1;		// Peak search
  _PHC0    = 0.0;		// Zero order phase
  _PHC1    = 0.0;		// 1st order phase
  _PH_mod  = 0;			// How to do phasing
  _PKNL    = 0;			// 5th order phase
  _PPARMOD = 0;			// Parameter display
  _PSCAL   = 4;			// Peak picking
  _PSIGN   = 0;			// Unknown
  _REVERSE = 0;			// Spectrum Reverse
  _SF      = 400.0;		// Spect. Freq.
  _SI      = 32768;		// Size
  _SIGF1   = 0;			// Signal/Noise
  _SIGF2   = 0;			// Signal/Noise
  _SINO    = 0;			// Signal/Noise
  _SIOLD   = 0;			// Unknown
  _SREGLST = "<>";		// Unknown
  _SSB     = 2;			// Sin Window
  _STSI    = 0;			// 2D storage
  _STSR    = 0;			// 2D storage
  _SW_p    = 10.0*_SF;		// Sweep Width (Hz)
  _SYMM    = 0;			// Symmetrization
  _S_DEV   = 0;			// Std. Dev.
  _TDeff   = 0;			// Effective FFT size
  _TDoff   = 0;			// FFT offset
  _TI      = "<>";		// Unknown
  _TILT    = "no";		// Tilt flag (2D)
  _TM1     = 0;			// Unknown
  _TM2     = 0;			// Unknown
  _TOPLEV  = 0;			// Contour top level
  _USERP1  = "<user>";		// User Permission?
  _USERP2  = "<user>";		// Unknown
  _USERP3  = "<user>";		// Unknown
  _USERP4  = "<user>";		// Unknown
  _USERP5  = "<user>";		// Unknown
  _WDW     = 1;			// Window function
  _XDIM    = 64;		// Submatrix size
  _YMAX_p  = 0;			// 2D Max
  _YMIN_p  = 0;			// 2D Min
//sosi  soffset  = 0;			// Spectrum Offset
  }

void XWinProcPar::SetDefaults1(const string& fname)
  { SetDefaults(fname); }	// Use generic defaults

void XWinProcPar::SetDefaults2(const string& fname)
  {
  SetDefaults(fname);		// Use generic defaults
  _MC2  = 5;			// Transform type (echo-antiecho)
  _SSB  = 2;			// Sine window
  _STSI = 0;			// 2D storage
  _STSR = 0;			// 2D storage
  _WDW  = 2;			// Window function
  }

void XWinProcPar::Copy(const XWinProcPar& XWPP)
  {
  parfile   = XWPP.parfile;
  _TITLE    = XWPP._TITLE;
  _JCAMPDX  = XWPP._JCAMPDX;
  _DATATYPE = XWPP._DATATYPE;
  _ORIGIN   = XWPP._ORIGIN;
  _OWNER    = XWPP._OWNER;
  _DATE     = XWPP._DATE;
  _ABSF1    = XWPP._ABSF1;
  _ABSF2    = XWPP._ABSF2;
  _ABSG     = XWPP._ABSG;
  _ABSL     = XWPP._ABSL;
  _XALPHA    = XWPP._XALPHA;
  _AQORDER  = XWPP._AQORDER;
  _ASSFAC   = XWPP._ASSFAC;
  _ASSFACI  = XWPP._ASSFACI;
  _ASSFACX  = XWPP._ASSFACX;
  _ASSWID   = XWPP._ASSWID;
  _AUNMP    = XWPP._AUNMP;	// For AU
  _AZFE     = XWPP._AZFE;	// Intetral ends
  _AZFW     = XWPP._AZFW;	// Multiplet integral
  _BCFW     = XWPP._BCFW;	// Unknown
  _BC_mod   = XWPP._BC_mod;	// BC modify
  _BYTORDP  = XWPP._BYTORDP;	// Byte ordering
  _COROFFS  = XWPP._COROFFS;	// Unknown
  _DATMOD   = XWPP._DATMOD;	// Add spectra
  _DC       = XWPP._DC;		// For spectra add
  _DFILT    = XWPP._DFILT;	// Filter file
  _DTYPP    = XWPP._DTYPP;	// Unknown
  _FCOR     = XWPP._FCOR;	// 1st pt multiply
  _FTSIZE   = XWPP._FTSIZE;	// Unknown
  _FT_mod   = XWPP._FT_mod;	// FT aquire mod
  _GAMMA    = XWPP._GAMMA;	// Quad correct. 
  _GB       = XWPP._GB;		// Gauss. multiply
  _INTBC    = XWPP._INTBC;	// Integration param.
  _INTSCL   = XWPP._INTSCL;	// Integral compare
  _ISEN     = XWPP._ISEN;	// Small integral
  _LB       = XWPP._LB;		// Line broadening
  _LEV0     = XWPP._LEV0;	// Contouring param
  _LPBIN    = XWPP._LPBIN;	// Linear prediction
  _MAXI     = XWPP._MAXI;	// Peak picking
  _MC2      = XWPP._MC2;	// 2D tranform type
  _MEAN     = XWPP._MEAN;	// Unknown
  _ME_mod   = XWPP._ME_mod;	// Linear prediciton
  _MI       = XWPP._MI;		// Peak picking
  _NCOEF    = XWPP._NCOEF;	// LP coefficients
  _NC_proc  = XWPP._NC_proc;	// 2D processing
  _NLEV     = XWPP._NLEV;	// Contour levels
  _NOISF1   = XWPP._NOISF1;	// SN ratio
  _NOISF2   = XWPP._NOISF2;	// SN ratio
  _NSP      = XWPP._NSP;	// Left shift
  _NTH_PI   = XWPP._NTH_PI;	// Unknown
  _NZP      = XWPP._NZP;	// For zp
  _OFFSET   = XWPP._OFFSET;	// 1st pt shift (Hz)
  __PC      = XWPP.__PC;	// Peak search
  _PHC0     = XWPP._PHC0;	// Zero order phase
  _PHC1     = XWPP._PHC1;	// 1st order phase
  _PH_mod   = XWPP._PH_mod;	// How to do phasing
  _PKNL     = XWPP._PKNL;	// 5th order phase
  _PPARMOD  = XWPP._PPARMOD;	// Parameter display
  _PSCAL    = XWPP._PSCAL;	// Peak picking
  _PSIGN    = XWPP._PSIGN;	// Unknown
  _REVERSE  = XWPP._REVERSE;	// Spectrum Reverse
  _SF       = XWPP._SF;		// Spect. Freq.
  _SI       = XWPP._SI;		// Size
  _SIGF1    = XWPP._SIGF1;	// Signal/Noise
  _SIGF2    = XWPP._SIGF2;	// Signal/Noise
  _SINO     = XWPP._SINO;	// Signal/Noise
  _SIOLD    = XWPP._SIOLD;	// Unknown
  _SREGLST  = XWPP._SREGLST;	// Unknown
  _SSB      = XWPP._SSB;	// Sin Window
  _STSI     = XWPP._STSI;	// 2D storage
  _STSR     = XWPP._STSR;	// 2D storage
  _SW_p     = XWPP._SW_p;	// Sweep Width (ppm)
  _SYMM     = XWPP._SYMM;	// Symmetrization
  _S_DEV    = XWPP._S_DEV;	// Std. Dev.
  _TDeff    = XWPP._TDeff;	// Effective FFT size
  _TDoff    = XWPP._TDoff;	// FFT offset
  _TI       = XWPP._TI;		// Unknown
  _TILT     = XWPP._TILT;	// Tilt flag (2D)
  _TM1      = XWPP._TM1;	// Unknown
  _TM2      = XWPP._TM2;	// Unknown
  _TOPLEV   = XWPP._TOPLEV;	// Contour top level
  _USERP1   = XWPP._USERP1;	// User Permission?
  _USERP2   = XWPP._USERP2;	// Unknown
  _USERP3   = XWPP._USERP3;	// Unknown
  _USERP4   = XWPP._USERP4;	// Unknown
  _USERP5   = XWPP._USERP5;	// Unknown
  _WDW      = XWPP._WDW;	// Window function
  _XDIM     = XWPP._XDIM;	// Submatrix size
  _YMAX_p   = XWPP._YMAX_p;
  _YMIN_p   = XWPP._YMIN_p;
  }


// ____________________________________________________________________________
// iii                  XWinNMR ProcPar Inter-Related Parameters
// ____________________________________________________________________________
 
/* Some of the parameters in an XWinNMR processing parameter files dependent on
   on another.  Thus when one parameter is specified or altered, on or more
   other parameters must be adjusted too.  These functions are meant to take
   care of such details.                                                     */
 
/* This sets the value of _OFFSET in the parameter set.  In this case _OFFSET
   is maintained in PPM and is the value of the 1st point in the spectrum. 
   If the center of the spectrum is NOT offset, this this value needs to be
   simply _SW_p*_SF/2.0 since (conveniently) _SW_p is stored in Hz...        */

void XWinProcPar::SetOffset()		// Sets OFFSET, 1st point shift (Hz)
{ _OFFSET = _SW_p*_SF/2.0; }	// OFFSET = spectrum offset + SW/2
//  { _OFFSET = soffset+_SW_p*_SF/2.0; }	// OFFSET = spectrum offset + SW/2
 

 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A          XWinProcPar Parameter File Constructors, Destructor
// ____________________________________________________________________________
 
/* These are constructors of the class handling Bruker XWinNMR 1D processing
   parameter files.  This doesn't do anything in particular, it is the write
   functions that perform the work.                                          */

XWinProcPar::XWinProcPar() : XWinPSet() { SetDefaults1("procs");}

XWinProcPar::XWinProcPar(const string& fname, int type) : XWinPSet(fname)
  {
  switch(type)
    {
    case 1:
    default: SetDefaults1(fname); break;        // Set defaults for procs
    case 2:  SetDefaults2(fname); break;        // Set defaults for proc2s
    }  
  }   


XWinProcPar::XWinProcPar(const XWinProcPar& XWP) : XWinPSet(XWP) { Copy(XWP); }

XWinProcPar::~XWinProcPar()  { }

void XWinProcPar::operator= (const XWinProcPar& XWP)
  { XWinPSet::operator= (XWP); Copy(XWP); }

// ____________________________________________________________________________
// B                  XWinProcPar Parameter Access Functions
// ____________________________________________________________________________

/* These functions allow users to set any more important parameters that
   these files contain.  The two primary values are the spectrum point size
   (SI) and the data byte order (BYTORDP).                                   */

string XWinProcPar::parname()   const { return parfile;  } // File name
int    XWinProcPar::BYTORDP()   const { return _BYTORDP; } // Byte order
int    XWinProcPar::FT_mod()    const { return _FT_mod;  } // FFT mode
double XWinProcPar::LB()        const { return _LB;      } // Line Broadening
int    XWinProcPar::MC2()       const { return _MC2;     } //
double XWinProcPar::OFFSET()    const { return _OFFSET;  } // Spectrum offset
double XWinProcPar::PHC0()      const { return _PHC0;    } // Zero order phase
double XWinProcPar::PHC1()      const { return _PHC1;    } // 1st order phase
int    XWinProcPar::PH_mod()    const { return _PH_mod;  } // Type of phasing
int    XWinProcPar::REVERSE()   const { return _REVERSE; } // Spectrum reverse
double XWinProcPar::SF()        const { return _SF;      } // Spectrom. freq.
int    XWinProcPar::SI()        const { return _SI;      } // Size (re+im)
int    XWinProcPar::SSB()       const { return _SSB;     } // Sine bell
int    XWinProcPar::STSI()      const { return _STSI;    } // Sine bell
int    XWinProcPar::STSR()      const { return _STSR;    } // Sine bell
double XWinProcPar::SW_p()      const { return _SW_p;    } // Spec. width(ppm)
double XWinProcPar::TDeff()     const { return _TDeff;   } // Effective FID size
int    XWinProcPar::WDW()       const { return _WDW;     } // Window function

void XWinProcPar::BYTORDP(int bo)         { bo?_BYTORDP=1:_BYTORDP=0; }
void XWinProcPar::GB(int gb)              { _GB      = gb; }
void XWinProcPar::FT_mod(int ft)          { _FT_mod  = ft; }
void XWinProcPar::FT_mod(const string& ft)
  { 
       if(ft == string("NO")  || ft == string("no"))  _FT_mod = 0; 
  else if(ft == string("ISR") || ft == string("isr")) _FT_mod = 1; 
  else if(ft == string("IQC") || ft == string("iqc")) _FT_mod = 2; 
  else if(ft == string("IQR") || ft == string("iqr")) _FT_mod = 3; 
//  else if(ft == string("FSC") || ft == string("fsc")) _FT_mod = 5; 
//  else if(ft == string("FSR") || ft == string("fsr")) _FT_mod = 4; 
  else if(ft == string("FSC") || ft == string("fsc")) _FT_mod = 4; 
  else if(ft == string("FSR") || ft == string("fsr")) _FT_mod = 5; 
  else if(ft == string("FQC") || ft == string("fqc")) _FT_mod = 6; 
  else if(ft == string("FQR") || ft == string("fqr")) _FT_mod = 7; 
  else if(ft == string("ISC") || ft == string("isc")) _FT_mod = 8; 
  else 
    {
    XWPPerror(22, "FT_mod", 1);
    XWPPerror(23, "no,isr,iqc,iqr,fsc,fsr,fqc,fqr,isc", 1); 
    XWPPerror(24, "fqr");  _FT_mod = 7;
    }
  }
void XWinProcPar::LB(double lb)           { _LB      = lb; }
void XWinProcPar::MC2(int mc)             { _MC2     = mc; }
void XWinProcPar::MC2(const string& mc)
  { 
       if(mc == string("QF")            || mc == string("qf"))            _MC2 = 0; 
  else if(mc == string("QSEQ")          || mc == string("qseq"))          _MC2 = 1; 
  else if(mc == string("TPPI")          || mc == string("tppi"))          _MC2 = 2; 
  else if(mc == string("States")        || mc == string("states"))        _MC2 = 3; 
  else if(mc == string("States-TPPI")   || mc == string("states-tppi"))   _MC2 = 4; 
  else if(mc == string("Echo-Antiecho") || mc == string("echo-antiecho")) _MC2 = 5; 
  else 
    {
    XWPPerror(22, "MC2", 1);
    XWPPerror(23, "QF, QSEQ, TPPI, States, States-TPPI, Echo-Antiecho", 1); 
    XWPPerror(24, "TPPI");  _MC2 = 2;
    }
  }
void XWinProcPar::OFFSET(double off)      { _OFFSET  = off;}
void XWinProcPar::PHC0(double ph0)        { _PHC0    = ph0;}
void XWinProcPar::PHC1(double ph1)        { _PHC1    = ph1;}
void XWinProcPar::PH_mod(int  phm)        { _PH_mod  = phm;}
void XWinProcPar::PH_mod(const string& p)
  { 
       if(p == string("no")  || p == string("NO")) _PH_mod = 0; 
  else if(p == string("pk")  || p == string("PK")) _PH_mod = 1; 
  else if(p == string("mc")  || p == string("MC")) _PH_mod = 2; 
  else if(p == string("ps")  || p == string("PS")) _PH_mod = 3; 
  else 
    {
    XWPPerror(22, "PH_mod", 1);
    XWPPerror(23, "no, pk, mc, ps", 1); 
    XWPPerror(24, "no");  _PH_mod = 0;
    }
  }

void XWinProcPar::REVERSE(int yn)         { yn?_REVERSE=1:_REVERSE=0;    }
void XWinProcPar::PPARMOD(int pm)         { _PPARMOD = pm; }
void XWinProcPar::SI(int si)              { _SI      = si; }
void XWinProcPar::SF(double SF)           { _SF      = SF; SetOffset();  }
void XWinProcPar::SSB(int sb)             { _SSB     = sb; }
void XWinProcPar::STSI(int st)            { _STSI    = st; }
void XWinProcPar::STSR(int sr)            { _STSR    = sr; }
void XWinProcPar::SW_p(double swp)        { _SW_p    = swp; SetOffset(); }
void XWinProcPar::WDW(int wd)             { _WDW     = wd; }
void XWinProcPar::WDW(const string& wd)
  { 
       if(wd == string("NO")    || wd == string("no"))    _WDW = 0; 
  else if(wd == string("EM")    || wd == string("em"))    _WDW = 1; 
  else if(wd == string("GM")    || wd == string("gm"))    _WDW = 2; 
  else if(wd == string("SINE")  || wd == string("sine"))  _WDW = 3; 
  else if(wd == string("QSINE") || wd == string("qsine")) _WDW = 4; 
  else if(wd == string("TRAP")  || wd == string("trap"))  _WDW = 5; 
  else if(wd == string("USER")  || wd == string("user"))  _WDW = 6; 
  else if(wd == string("SINC")  || wd == string("sinc"))  _WDW = 7; 
  else if(wd == string("QSINC") || wd == string("qsinc")) _WDW = 8; 
  else if(wd == string("TRAF")  || wd == string("traf"))  _WDW = 9; 
  else if(wd == string("TRAF5") || wd == string("traf5")) _WDW = 10; 
  else 
    {
    XWPPerror(22, "WDW", 1);
    string tmp("NO,EM,GM,SINE,QSINE");
    tmp += string("\n") + string(30, ' ')
        + string("TRAP,USER,SINC,QSINC,TRAF,TRAF5");
    XWPPerror(23, tmp, 1);
    XWPPerror(24, "EM");  _WDW = 1;
    }
  }

// ____________________________________________________________________________
// C                       XWinProcPar Input Functions
// ____________________________________________________________________________

/* These functions will read in the processing parameters from an XWinNMR
   parameter file, typically named proc, procs, proc2, proc2s,.....  By design,
   the Bruker parameter file is initially read into a GAMMA parameter set so
   that ALL parameters in the file are stored (class XWinPSet).  Subsequently,
   the parameters in the parameter set are parsed (herein) to obtain values of
   consequence to GAMMA and these are explicitly maintained variables in this
   class.

      Function                               Purpose
    ____________________________________________________________________________

       read              Read in parameter set for class object.  This is done
		         special since Bruker format is NOT in GAMMA parameter
		         format so it must be parsed appropriately.
      parsePSet          Converts parameters in internal pset to specific values
       getPar            Returns parameter found in the parameter set herein
			 This function is inherited from base class XWinPSet. */

bool XWinProcPar::readPPar(const string& filein, int warn)
  {
  parfile = filein;				// Set internal filename
  bool TF = readPPar(warn);			// Read in all parameters
  if(!TF && warn)                               // Output errors if trouble
    {
    XWPPerror(1, filein, 1);		// Problems with filein
    if(warn > 1) XWPPerror(3);		// Can't construct from pset
    else         XWPPerror(3,1);
    }
  return TF;
  }

bool XWinProcPar::readPPar(int warn)
  {
  bool TF = XWinPSet::readPSet(parfile, warn);	// Use base class to read
  if(!TF) return TF;                            // Fail if we cannot read it
  return parsePSet(warn);			// Parse our essential params
  }

bool XWinProcPar::parsePSet(int warn)
  {
  string pvalue;				// String value for # points
  getPar("TITLE",     _TITLE); 			// Try & get the title
  getPar("JCAMPDX",   _JCAMPDX); 		// Try & get JCAMP version
  getPar("DATATYPE",  _DATATYPE); 		// Try & get DATATYPE 
  getPar("ORIGIN",    _ORIGIN);			// Try & get the origin
  getPar("OWNER",     _OWNER);			// Try & get the owner
  getPar("DATE",      _DATE); 			// Try & get the date
  getPar("ABSF1",     _ABSF1);
  getPar("ABSF2",     _ABSF2);
  getPar("ABSG",      _ABSG);
  getPar("ABSL",      _ABSL);
  getPar("ALPHA",     _XALPHA);
  getPar("AQORDER",   _AQORDER);
  getPar("ASSFAC",    _ASSFAC);
  getPar("ASSFACI",   _ASSFACI);
  getPar("ASSFACX",   _ASSFACX);
  getPar("ASSWID",    _ASSWID);
  getPar("AUNMP",     _AUNMP);
  getPar("AZFE",      _AZFE);
  getPar("AZFW",      _AZFW);
  getPar("BCFW",      _BCFW);
  getPar("BC_mod",    _BC_mod);
  if(!getPar("BYTORDP",_BYTORDP)) _BYTORDP = 0;
  getPar("COROFFS",   _COROFFS);
  getPar("DATMOD",   _DATMOD);
  getPar("DC",       _DC);
  getPar("DFILT",   _DFILT);
  getPar("DTYPP",   _DTYPP);
  getPar("FCOR",    _FCOR);
  getPar("FTSIZE",  _FTSIZE);
  getPar("FT_mod",  _FT_mod);
  getPar("GAMMA",   _GAMMA);
  getPar("GB",      _GB);
  getPar("INTBC",   _INTBC);
  getPar("INTSCL",  _INTSCL);
  getPar("ISEN",    _ISEN);
  getPar("LB",      _LB);
  getPar("LEV0",    _LEV0);
  getPar("LPBIN",   _LPBIN);
  getPar("MAXI",    _MAXI);
  getPar("MC2",     _MC2);
  getPar("MEAN",    _MEAN);
  getPar("ME_mod",  _ME_mod);
  getPar("MI",      _MI);
  getPar("NCOEF",   _NCOEF);
  getPar("NC_proc", _NC_proc);
  getPar("NLEV",    _NLEV);
  getPar("NOISF1",  _NOISF1);
  getPar("NOISF2",  _NOISF2);
  getPar("NSP",     _NSP);
  getPar("NTH_PI",  _NTH_PI);
  getPar("NZP",     _NZP);
  getPar("OFFSET",  _OFFSET);
  getPar("PC",      __PC);
  getPar("PHC0",    _PHC0);
  getPar("PHC1",    _PHC1);
  getPar("PH_mod",  _PH_mod);
  getPar("PKNL",    _PKNL);
  getPar("PPARMOD", _PPARMOD);
  getPar("PSCAL",   _PSCAL);
  getPar("PSIGN",   _PSIGN);
  getPar("REVERSE", _REVERSE);
  if(!getPar("SF",  _SF,23,1)) XWPPfatal(2, "SF");	// Try to get spect. freq.
  if(!getPar("SI",  _SI,23,1)) XWPPfatal(2, "SI");	// Try to get spect. freq.
  getPar("SIGF1",   _SIGF1);
  getPar("SIGF2",   _SIGF2);
  getPar("SINO",    _SINO);
  getPar("SIOLD",   _SIOLD);
  getPar("SREGLST", _SREGLST);
  getPar("SSB",     _SSB);
  getPar("STSI",    _STSI);
  getPar("STSR",    _STSR);
  if(!getPar("SW_p",_SW_p,24,1)) XWPPfatal(2, "SW");	// Try to get spect. width
  getPar("SYMM",     _SYMM);
  getPar("S_DEV",    _S_DEV);
  getPar("TDeff",    _TDeff);
  getPar("TDoff",    _TDoff);
  getPar("TI",       _TI);
  getPar("TILT",     _TILT);
  getPar("TM1",      _TM1);
  getPar("TM2",      _TM2);
  getPar("TOPLEV",   _TOPLEV);
  getPar("USERP2",   _USERP2);
  getPar("USERP3",   _USERP3);
  getPar("USERP4",   _USERP4);
  getPar("USERP5",   _USERP5);
  getPar("WDW",      _WDW);
  getPar("XDIM",     _XDIM);
  getPar("YMAX_p",   _YMAX_p);
  getPar("YMIN_p",   _YMIN_p);
  return true;
  }

//int XWinProcPar::getPar(const string& n,string& v,int i, int w) const;

// ____________________________________________________________________________
// D                         XWinProcPar Output Functions
// ____________________________________________________________________________
 
/* These functions allow for output of NMR parameters directly into a Bruker
   XWinNMR ASCII parameter file (procs).                                     */

int XWinProcPar::writePPar(const string& dname, int warn)
  { parfile = dname; return writePPar(warn); }

int XWinProcPar::writePPar(int warn) const
  {
  ofstream ofstr(parfile.c_str());		// Open ASCII file for output
  if(!ofstr.good())				// If file bad then exit
     {
     if(warn)
       {
       XWPPerror(2, 1);				// Output filestream problems
       XWPPerror(25,1);				// Can't write parameters out
       if(warn==1) XWPPerror(21,parfile,1);	// Can't write parameters out
       else        XWPPfatal(21,parfile);	// with out without die
       } 
     return 0;
     }  
  string nn("##");			// Bruker ## parameter line start
  string nns("##$");			// Bruker ##$ parameter line start
  string ss("$$");			// Bruker $$ parameter line start
  ofstr << nn << "TITLE= "    << _TITLE << "\n";
  ofstr << nn << "JCAMPDX= "  << _JCAMPDX << "\n";
  ofstr << nn << "DATATYPE= " << _DATATYPE << "\n";
  if(_ORIGIN != "")
    ofstr << nn << "ORIGIN= "   << _ORIGIN << "\n";
  if(_OWNER != "")
    ofstr << nn << "OWNER= "    << _OWNER << "\n";
  time_t longtime;                              // Need a time structure
  longtime = NULL;

#ifdef _MSC_VER
  struct tm newtime;		// For setting current date
  localtime_s(&newtime, &longtime);
	char date[31];
	asctime_s(date,30,&newtime);
	ofstr << ss << " " << date;
#else
  struct tm *ptr;
  ptr = localtime(&longtime);
  ofstr << ss << " " << asctime(ptr);		// Output current date & time
#endif

  string cwd;

#ifdef _MSC_VER
	char cwdtmp[130];
  _getcwd(cwdtmp, 128);		// Current working directory
	cwd = string(cwdtmp);
#else
  cwd = string(getcwd(NULL, 128));		// Current working directory
#endif

  ofstr << ss << " " << cwd			// Output full file name
        << "/" << parfile << "\n";
  ofstr << nns << "ABSF1= " << _ABSF1   << "\n";// Left limit for bc (ppm)
  ofstr << nns << "ABSF2= " << _ABSF2   << "\n";// Right limit for bc (ppm)
  ofstr << nns << "ABSG= "  << _ABSG    << "\n";// Polynomical degree for bc
  ofstr << nns << "ABSL= "  << _ABSL    << "\n";// Integral detection level
  ofstr << nns << "ALPHA= " << _XALPHA  << "\n";// Quadrature correction const.
  ofstr << nns << "AQORDER= "<<_AQORDER << "\n";// 3D processing order
  ofstr << nns << "ASSFAC= "<< _ASSFAC  << "\n";// Vertical scaling flag
  ofstr << nns << "ASSFACI= "<<_ASSFACI << "\n";// Undocumented
  ofstr << nns << "ASSFACX= "<<_ASSFACX << "\n";// Undocumented
  ofstr << nns << "ASSWID= "<< _ASSWID  << "\n";// Search interval for ASSFAC
  ofstr << nns << "AUNMP= " << _AUNMP   << "\n";// For AU
  ofstr << nns << "AZFE= "  << _AZFE    << "\n";// Integration ends (ppm)
  ofstr << nns << "AZFW= "  << _AZFW    << "\n";// For multiplet integration
  ofstr << nns << "BCFW= "  << _BCFW    << "\n";// Undocumented
  ofstr << nns << "BC_mod= "<< _BC_mod  << "\n";// Baseline correct modifiy
						// no          - do nothing
						// single/quad - subtract constant
						// spol/qpol   - polynomial subtr.
						// sfil/qfil   - Bax/Marion filter
  ofstr << nns << "BYTORDP= "<<_BYTORDP << "\n";// Byte ord (1=big,0=little endian)
  ofstr << nns << "COROFFS= "<< _COROFFS<< "\n";// Undocumented
  ofstr << nns << "DATMOD= " << _DATMOD << "\n";// Add spectra (0=raw, 1=processed)
  ofstr << nns << "DC= "     << _DC     << "\n";// Used when adding spectra
  ofstr << nns << "DFILT= "  << _DFILT  << "\n";// Filter coefficients file
  ofstr << nns << "DTYPP= "  << _DTYPP  << "\n";// Undocumented
  ofstr << nns << "FCOR= "   << _FCOR   << "\n";// 1st point multiply factor
  ofstr << nns << "FTSIZE= " << _FTSIZE << "\n";// Undocumented
  ofstr << nns << "FT_mod= " << _FT_mod << "\n";// How FT's are performed
  ofstr << nns << "GAMMA= "  << _GAMMA  << "\n";// Quadrature correction const.
  ofstr << nns << "GB= "     << _GB     << "\n";// Used in Gaussian multiply
  ofstr << nns << "INTBC= "  << _INTBC  << "\n";// Integration parameter
  ofstr << nns << "INTSCL= " << _INTSCL << "\n";// Integral comparison parameter
  ofstr << nns << "ISEN= "   << _ISEN   << "\n";// Small integral discard 
  ofstr << nns << "LB= "     << _LB     << "\n";// Line broadening
  ofstr << nns << "LEV0= "   << _LEV0   << "\n";// Contouring parameter
  ofstr << nns << "LPBIN= "  << _LPBIN  << "\n";// Linear prediciton value
  ofstr << nns << "MAXI= "   << _MAXI   << "\n";// Peak picking max intensity
  ofstr << nns << "MC2= "    << _MC2    << "\n";// 2D transform type
  ofstr << nns << "MEAN= "   << _MEAN   << "\n";//
  ofstr << nns << "ME_mod= " << _ME_mod << "\n";// Linear prediciton flag
  ofstr << nns << "MI= "     << _MI     << "\n";// Peak picking min intensity
  ofstr << nns << "NCOEF= "  << _NCOEF  << "\n";// Coefficients linear prediction
  ofstr << nns << "NC_proc= "<< _NC_proc<< "\n";// For 2D processing check
  ofstr << nns << "NLEV= "   << _NLEV   << "\n";// Contour levels
  ofstr << nns << "NOISF1= " << _NOISF1 << "\n";// For calculating S/N ratios
  ofstr << nns << "NOISF2= " << _NOISF2 << "\n";// For calculating S/N ratios
  ofstr << nns << "NSP= "    << _NSP    << "\n";// Left shift (points)
  ofstr << nns << "NTH_PI= " << _NTH_PI << "\n";// Undocumented
  ofstr << nns << "NZP= "    << _NZP    << "\n";// For zp command
  ofstr << nns << "OFFSET= " << _OFFSET << "\n";// Shift of 1st point (ppm)
  ofstr << nns << "PC= "     << __PC    << "\n";// Peak search sensitivity
  ofstr << nns << "PHC0= "   << _PHC0   << "\n";// Zero order phase
  ofstr << nns << "PHC1= "   << _PHC1   << "\n";// First order phase
  ofstr << nns << "PH_mod= " << _PH_mod << "\n";// How to do phasing
  ofstr << nns << "PKNL= "   << _PKNL   << "\n";// 5th order phase correction
  ofstr << nns << "PPARMOD= "<< _PPARMOD<< "\n";// For parameter display
  ofstr << nns << "PSCAL= "  << _PSCAL  << "\n";// For peak picking
  ofstr << nns << "PSIGN= "  << _PSIGN  << "\n";// Undocumented
  ofstr << nns << "REVERSE= "<< _REVERSE<< "\n";// Reverse spectrum
  ofstr << nns << "SF= "     << _SF     << "\n";// Spectrometer freq.
  ofstr << nns << "SI= "     << _SI     << "\n";// Size (TD -> power of 2)
  ofstr << nns << "SIGF1= "  << _SIGF1  << "\n";// For S/N ratio
  ofstr << nns << "SIGF2= "  << _SIGF2  << "\n";// For S/N ratio
  ofstr << nns << "SINO= "   << _SINO   << "\n";// S/N ratio
  ofstr << nns << "SIOLD= "  << _SIOLD  << "\n";// Undocumented
  ofstr << nns << "SREGLST= "<< _SREGLST<< "\n";// Undocumented
  ofstr << nns << "SSB= "    << _SSB    << "\n";// For Sin window functions
  ofstr << nns << "STSI= "   << _STSI   << "\n";// 2D matrix storage 
  ofstr << nns << "STSR= "   << _STSR   << "\n";// 2D matrix storage
  ofstr << nns << "SW_p= "   << _SW_p   << "\n";// Sweep Width (ppm)
  ofstr << nns << "SYMM= "   << _SYMM   << "\n";// Symmetrization
  ofstr << nns << "S_DEV= "  << _S_DEV  << "\n";// 2D standard dev.
  ofstr << nns << "TDeff= "  << _TDeff  << "\n";// Effective size for FFT
  ofstr << nns << "TDoff= "  << _TDoff  << "\n";// FFT offset
  ofstr << nns << "TI= "     << _TI     << "\n";// Undocumented
  ofstr << nns << "TILT= "   << _TILT   << "\n";// Tilt flag for 2D spectra
  ofstr << nns << "TM1= "    << _TM1    << "\n";// Undocumented
  ofstr << nns << "TM2= "    << _TM2    << "\n";// Undocumented
  ofstr << nns << "TOPLEV= " << _TOPLEV << "\n";// Contouring top level
  ofstr << nns << "USERP1= " << _USERP1 << "\n";// User permission files?
  ofstr << nns << "USERP2= " << _USERP2 << "\n";
  ofstr << nns << "USERP3= " << _USERP3 << "\n";
  ofstr << nns << "USERP4= " << _USERP4 << "\n";
  ofstr << nns << "USERP5= " << _USERP5 << "\n";
  ofstr << nns << "WDW= "    << _WDW    << "\n";// Window function
  ofstr << nns << "XDIM= "   << _XDIM   << "\n";// Sub-matrix size
  ofstr << nns << "YMAX_p= " << _YMAX_p << "\n";// Max of 2D spectrum
  ofstr << nns << "YMIN_p= " << _YMIN_p << "\n";// Min of 2D spectrum
  ofstr << nn  << "END= "    <<            "\n";
  return 1;
  }
 
// ____________________________________________________________________________
// E                 XWinProcPar Parameter Access Functions
// ____________________________________________________________________________
 
/* These functions allow access to the Bruker XWinNMR parameter set, not just
  those parameters which are directly handled within this class. 
 
                       INHERITED FROM BASE CLASS XWinPSet 
 
int XWinProcPar::getPar(const string& pn,int& val,   int id,int warn) const
int XWinProcPar::getPar(const string& pn,double& val,int id,int warn) const
int XWinProcPar::getPar(const string& pn,string& val,int id,int warn) const
ParameterSet XWinProcPar::getPSet() const                                    */

 
// ____________________________________________________________________________
// F                 XWinProcPar Standard Output Functions
// ____________________________________________________________________________
 
// ostream& XWinProcPar::printPset(ostream& O) const		INHERITED

ostream& XWinProcPar::print(ostream& ostr, int full, int hdr) const
  {
  string marg(16, ' ');
  string ls = string("\n") + marg;
  if(hdr)
    ostr << ls << "Procs Parameter File:       " << parfile;
  XWinPSet::print(ostr, 0, 0);
  if(!full)
    {
    ostr << ls << "Owner:                      " << _OWNER;
    ostr << ls << "FT Mode:                    " << _FT_mod;
    ostr << ls << "Line Broadening:            " << _LB;
    ostr << ls << "Zeroth Order Phase:         " << _PHC0;
    ostr << ls << "First Order Phase:          " << _PHC1;
    ostr << ls << "Spectrum Reverse:           " << _REVERSE;
    ostr << ls << "Spectrometer Frequency:     " << _SF                << " MHz";
    ostr << ls << "Real Data Point Size:       " << _SI                << " Pts";
    ostr << ls << "Spectral Width:             " << _SW_p              << " PPM";
    ostr << ls << "Effective FFT Size:         " << _TDeff             << " Pts";
    }
  else
    {
    ostr << ls << "Title:                      " << _TITLE;
    ostr << ls << "JCAMP Version:              " << _JCAMPDX;
    ostr << ls << "Data Type:                  " << _DATATYPE;
    if(_ORIGIN != "")
      ostr << ls << "File Origin:              " << _ORIGIN;
    if(_OWNER != "")
      ostr << ls << "File Owner:                 " << _OWNER;
    ostr << ls << "File Date:                  " << _DATE;
    ostr << ls << "FT Mode:                    " << _FT_mod;
    ostr << ls << "Line Broadening:            " << _LB;
    ostr << ls << "Zeroth Order Phase:         " << _PHC0;
    ostr << ls << "First Order Phase:          " << _PHC1;
    ostr << ls << "Spectrum Reverse:           " << _REVERSE;
    ostr << ls << "Spectrometer Frequency:     " << _SF                << " MHz";
    ostr << ls << "Real Data Point Size:       " << _SI                << " Pts";
    ostr << ls << "Spectral Width:             " << _SW_p              << " PPM";
    ostr << ls << "Effective FFT Size:         " << _TDeff             << " Pts";
    }
  return ostr;
  }  

ostream& operator<< (ostream& O, const XWinProcPar& P) { P.print(O); return O; }

// ____________________________________________________________________________
// G                      XWinProcPar Auxiliary Functions
// ____________________________________________________________________________

 /* These are just a group of functions that return strings indicating what
    the Bruker parameter values mean.  I just add things as I learn them here
    so that GAMMA output can remind me what all of these parameters are...   */
 
string XWinProcPar::AQORDERS() const
  {
  string ret;
  switch(_AQORDER)
    {
    case 0:  ret = "3-2-1"; break;
    case 1:  ret = "3-1-2"; break;
    default: ret = "Unknown";
    }
  return ret;
  }

string XWinProcPar::BC_modS()   const
  {
  string ret;
  switch(_BC_mod)
    {
    case 0:  ret = "no"; break;
    case 1:  ret = "single"; break;
    case 2:  ret = "quad"; break;
    case 3:  ret = "spol"; break;
    case 4:  ret = "qpol"; break;
    default: ret = "Unknown";
    }
  return ret;
  }

string XWinProcPar::BYTORDPS() const
  {
  string ret;
  switch(_BYTORDP)
    {
    case 0:  ret = "Little Endian"; break;
    case 1:  ret = "Big Endian";    break;
    default: ret = "Default";
    }
  return ret;
  }

string XWinProcPar::DATMODS()   const
  {
  string ret;
  switch(_DATMOD)
    {
    case 0:  ret = "raw"; break;
    case 1:  ret = "proc"; break;
    default: ret = "Unknown";
    }
  return ret;
  }

string XWinProcPar::FT_modS()   const
  {
  string ret;
  switch(_FT_mod)
    {
    case 0:  ret = "no"; break;
    case 1:  ret = "isr"; break;
    case 2:  ret = "iqc"; break;
    case 3:  ret = "iqr"; break;
    case 4:  ret = "fsc"; break;
    case 5:  ret = "fsr"; break;
    case 6:  ret = "fqc"; break;
    case 7:  ret = "fqr"; break;
    case 8:  ret = "isc"; break;
    default: ret = "Unknown";
    }
  return ret;
  }

string XWinProcPar::INTBCS() const
  {
  string ret;
  switch(_INTBC)
    {
    case 0:  ret = "no"; break;
    case 1:  ret = "yes"; break;
    default: ret = "Unknown";
    }
  return ret;
  }

string XWinProcPar::MC2S() const
  {
  string ret;
  switch(_MC2)
    {
    case 0:  ret = "QF";            break;
    case 1:  ret = "QSEQ";          break;
    case 2:  ret = "TPPI";          break;
    case 3:  ret = "States";        break;
    case 4:  ret = "States-TPPI";   break;
    case 5:  ret = "Echo-Antiecho"; break;
    default: ret = "Unknown";
    }
  return ret;
  }

string XWinProcPar::ME_modS()   const
  {
  string ret;
  switch(_ME_mod)
    {
    case 0:  ret = "no"; break;
    case 1:  ret = "LPbr,LPbc"; break;
    default: ret = "Unknown";
    }
  return ret;
  }

string XWinProcPar::PH_modS()   const
  {
  string ret;
  switch(_PH_mod)
    {
    case 0:  ret = "no"; break;
    case 1:  ret = "pk"; break;
    case 2:  ret = "mc"; break;
    case 3:  ret = "ps"; break;
    default: ret = "Unknown";
    }
  return ret;
  }

string XWinProcPar::PKNLS() const
  {
  string ret;
  switch(_PKNL)
    {
    case 0:  ret = "FALSE";          break;
    case 1:  ret = "TRUE";          break;
    default: ret = "Unknown";
    }
  return ret;
  }

string XWinProcPar::PPARMODS() const
  {
  string ret;
  switch(_PPARMOD)
    {
    case 0:  ret = "1D";          break;
    case 1:  ret = "2D";          break;
    case 2:  ret = "3D";          break;
    default: ret = "Unknown";
    }
  return ret;
  }

string XWinProcPar::PSCALS() const
  {
  string ret;
  switch(_PSCAL)
    {
    case 0:  ret = "global";   break;
    case 1:  ret = "preg";     break;
    case 2:  ret = "ireg";     break;
    case 3:  ret = "pireg";    break;
    case 4:  ret = "sreg";     break;
    case 5:  ret = "psreg";    break;
    case 6:  ret = "noise";    break;
    default: ret = "Unknown";
    }
  return ret;
  }

string XWinProcPar::PSIGNS() const
  {
  string ret;
  switch(_PSIGN)
    {
    case 0:  ret = "pos.";         break;
    case 1:  ret = "neg.";         break;
    default: ret = "Unknown";
    }
  return ret;
  }

string XWinProcPar::REVERSES() const
  {
  string ret;
  switch(_REVERSE)
    {
    case 0:  ret = "FALSE";          break;
    case 1:  ret = "TRUE";          break;
    default: ret = "Unknown";
    }
  return ret;
  }

string XWinProcPar::SYMMS() const
  {
  string ret;
  switch(_SYMM)
    {
    case 0:  ret = "no";          break;
    case 1:  ret = "yes";         break;
    default: ret = "Unknown";
    }
  return ret;
  }

string XWinProcPar::WDWS()   const
  {
  string ret;
  switch(_WDW)
    {
    case 0:  ret = "no"; break;
    case 1:  ret = "EM"; break;
    case 2:  ret = "GM"; break;
    case 3:  ret = "SINE"; break;
    case 4:  ret = "QSINE"; break;
    case 5:  ret = "TRAP"; break;
    default: ret = "Unknown";
    }
  return ret;
  }


#endif							// XWinProcPar.cc
