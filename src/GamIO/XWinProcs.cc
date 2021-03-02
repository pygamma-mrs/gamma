/* XWinProcs.cc **************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      XWinProcs                                  Implementation	**
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

#ifndef   XWinProcs_CC_			// Is file already included?
#  define XWinProcs_CC_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <GamIO/XWinProcs.h>		// Include the interface
#include <GamIO/XWinProcPar.h>		// Include the base class interface
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
#endif
#include <cmath>			// Include fabs function

using std::string;			// Using libstdc++ strings
using std::ostream;			// Using libstdc++ output streams


// ----------------------------------------------------------------------------
// ----------------------------- PRIVATE FUNCTIONS ----------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                      XWinNMR Procs File Error Handling
// ____________________________________________________________________________
 
/* These functions take care of any errors encountered when reading, writing,
   and setting parameters in Bruker acquisition parameter files.

        Input		eidx    : Error index
        		pname   : Additional error message
        		noret	: Flag for linefeed (0=linefeed)
        Output		void    : An error message is output                 */
 
void XWinProcs::XWinProcserror(int eidx, int noret) const
  {
  string hdr("XWinNMR Procs Parameter File");
  string msg;
  switch (eidx)
    {
    case 25: GAMMAerror(hdr,"Cannot Write Parameters To File",noret);break; // (25)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }  

    
void XWinProcs::XWinProcserror(int eidx, const string& pname, int noret) const
  {                                                                             
  string hdr("XWinNMR Procs Parameter File");
  string msg;
  switch(eidx)
    {
    case 21:msg = string("Cannot Write To ") + pname;
             GAMMAerror(hdr,msg+pname,noret);  break;                  // (21)
    default: GAMMAerror(hdr, eidx, pname, noret); break;// Unknown Error  (-1)
    }
  }

volatile void XWinProcs::XWinProcsfatality(int eidx) const
  {                                                                 
  XWinProcserror(eidx, 1);			// Normal non-fatal error
  if(eidx) XWinProcserror(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

volatile void XWinProcs::XWinProcsfatality(int eidx, const string& pname) const
  {                                                                 
  XWinProcserror(eidx, pname, 1);		// Normal non-fatal error
  if(eidx) XWinProcserror(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }
 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A          XWinProcs Parameter File Constructors, Destructor
// ____________________________________________________________________________
 
/* These are constructors of the class handling Bruker XWinNMR 1D processing
   parameter files.  This doesn't do anything in particular, it is the write
   functions that perform the work.                                          */

     XWinProcs::XWinProcs()                     : XWinProcPar("procs", 1) { }
     XWinProcs::XWinProcs(const string& name)   : XWinProcPar(name, 1)    { }
     XWinProcs::XWinProcs(const XWinProcs& XWP) : XWinProcPar(XWP)        { }
     XWinProcs::~XWinProcs()                                              { }
void XWinProcs::operator= (const XWinProcs& X) { XWinProcPar::operator= (X); }

// ____________________________________________________________________________
// B                  XWinProcs Parameter Access Functions
// ____________________________________________________________________________

/* These functions allow users to set some of the basic parameters that
   these files contain.  Remember, GAMMA doesn't intend to output processed
   data sets to XWinNMR. It only generates such files because the XWinNMR
   software demands it whenever a 1D acquisition is performed - or whenever
   a simulated FID is output in their format.  Input is quite different, as
   GAMMA does support the ability to read in a processed 1D spectrum in
   XWinNMR format as well as access all associated parameters.

                             INHERITED FROM XWinProcPar
                             ==========================
	  
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

// ____________________________________________________________________________
// C                       XWinProcs Input Functions
// ____________________________________________________________________________

/* These functions will read in the 1D processing parameters from an XWinNMR
   parameter file, typically named procs.  By design, the Bruker parameter
   file is initially read into a GAMMA parameter set so that ALL parameters in
   the file are stored (class XWinPSet).  Subsequently, the parameters in the
   parameter set are parsed (herein) to obtain values of consequence to GAMMA
   and these are explicitly maintained variables in this class.

      Function                               Purpose
    ____________________________________________________________________________

       read              Read in parameter set for class object.  This is done
		         special since Bruker format is NOT in GAMMA parameter
		         format so it must be parsed appropriately.
      parsePSet          Converts parameters in internal pset to specific values

                         INHERITED FROM CLASS XWinProcPar
                         ================================

bool XWinProcs::read(const string& filein, int warn)
bool XWinProcs::read(int warn)
bool XWinProcs::parsePSet(int warn)                                          */

// ____________________________________________________________________________
// D                         XWinProcs Output Functions
// ____________________________________________________________________________
 
/* These functions allow for output of NMR parameters directly into a Bruker
   XWinNMR ASCII parameter file (procs).

                         INHERITED FROM CLASS XWinProcPar
                         ================================

int XWinProcs::write(const string& dname, int warn)
int XWinProcs::write(int warn) const                                         */
 
// ____________________________________________________________________________
// E                    XWinProds Standard Output Functions
// ____________________________________________________________________________

/* These functions allow for a quick output of the parameter file contents.
   They don't have anything to do with output while running XWinNMR, rather
   users can just glance at procs parameters or store them in a small file.  */

// ostream& XWinProcs::printPset(ostream& O) const		INHERITED

/*
ostream& XWinProcs::print(ostream& ostr, int full, int hdr) const
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

ostream& operator<< (ostream& O, const XWinProcs& P) { P.print(O); return O; }

*/

// ____________________________________________________________________________
// F                    XProcPar BrukerLike Output Functions
// ____________________________________________________________________________

/* These functions perform ASCII output of processing parameters.  The format
   here is similar to issuing a "dpp" command in XWinNMR or using the menu
   choice "Output/Display status pars./Processing only".  It isn't as complete
   as the Bruker output because I'm bored and its often a mystery how they
   get some of their values.....  Most of the values output are taken from the
   parameter set which holds all the values which were read in from the
   Bruker parameter file.  It makes no sense to use these functions if on
   hasn't first read such a file since all the parameters will be empty.  A
   couple of values are used directly from the class when parameters need to
   be calculated.						             */

ostream& XWinProcs::dpp(ostream& ostr, double BF0) const
  {
  if(!BF0) BF0 = _SF;				// procs has NO BF0, yet it
  double SRval = (_SF-BF0)*1.e6;		// figures out SR like this!!
  double HZPPval = _SW_p/_SI;
  ostr << "\n" << string(29, ' ') << "F2 - Processing Parameters";
  ostr << "\n" << string(79, '=') << "\n\n";
  bru(ostr, "SI",      SI(),       "",                            0);
  bru(ostr, "PPARMOD", PPARMODS(), "",                            1);
  bru(ostr, "SF",      _SF,        "MHz",     bruform(_SF,8),     0);
  bru(ostr, "OFFSET",  _OFFSET,    "ppm",     bruform(_OFFSET,3), 1);
  bru(ostr, "SR",      SRval,      "Hz",      bruform(SRval, 2),  0);
  bru(ostr, "HzpPT",   HZPPval,    "Hz",      bruform(HZPPval,6), 1);
  bru(ostr, "SW_p",    _SW_p,      "Hz",      bruform(_SW_p,2),   0);
  bru(ostr, "XDIM",    _XDIM,      "",                            1);
  bru(ostr, "WDW",     WDWS(),     "",                            0);
  bru(ostr, "SSB",     _SSB,       "",                            1);
  bru(ostr, "LB",      _LB ,       "Hz",      bruform(_LB,2),     0); 
  bru(ostr, "GB",      _GB,        "",                            1);
  bru(ostr, "PH_mod",  PH_modS(),  "",                            0);
  bru(ostr, "PKNL",    PKNLS(),    "",                            1);
  bru(ostr, "PHC0",    _PHC0,      "degrees", bruform(_PHC0,3),   0);
  bru(ostr, "PHC1",    _PHC1,      "degrees", bruform(_PHC1,3),   1);
  bru(ostr, "BC_mod",  BC_modS(),  "",                            0);
  bru(ostr, "BCFW",    _BCFW,      "ppm",     bruform(_BCFW,3),   1);
  bru(ostr, "FT_mod",  FT_modS(),  "",                            0);
  bru(ostr, "FCOR",    _FCOR,      "",        bruform(_FCOR,1),   1);
  bru(ostr, "ME_mod",  ME_modS(),  "",                            0);
  bru(ostr, "COROFFS", _COROFFS,   "Hz",      bruform(_COROFFS,2),1);
  bru(ostr, "NCOEF",   _NCOEF,     "",                            0);
  bru(ostr, "LPBIN",   _LPBIN,     "",                            1);
  bru(ostr, "ABSF1",   _ABSF1,     "ppm",     bruform(_ABSF1,3),  0);
  bru(ostr, "ABSF2",   _ABSF2,     "ppm",     bruform(_ABSF2,3),  1);
  bru(ostr, "ABSG",    _ABSG,      "",                            0);
  bru(ostr, "ABSL",    _ABSL,      "",                            1);
  bru(ostr, "AZFE",    _AZFE,      "ppm",     bruform(_AZFE,3),   0);
  bru(ostr, "AZFW",    _AZFW,      "ppm",     bruform(_AZFW,3),   1);
  bru(ostr, "TDeff",   _TDeff,     "",                            0);
  bru(ostr, "TDoff",   _TDoff,     "",                            1);
  bru(ostr, "STSR",    _STSR,      "",                            0);
  bru(ostr, "STSI",    _STSI,      "",                            1);
  bru(ostr, "PSCAL",   PSCALS(),   "",                            0);
  bru(ostr, "SREGLST", _SREGLST,   "",                            1);
  bru(ostr, "PC",      __PC,        "",       bruform(__PC,2),    0);
  bru(ostr, "PSIGN",   PSIGNS(),   "",                            1);
  bru(ostr, "MI",      _MI,        "cm",     bruform(_MI,2),      0);
  bru(ostr, "MAXI",    _MAXI,      "cm",     bruform(_MAXI,2),    1);
  bru(ostr, "INTBC",   INTBCS(),   "",                            0);
//  bru(ostr, "INTSCL",  _INTSCL,    "",       "%12.5e",            1);
  bru(ostr, "INTSCL",  _INTSCL,    "",                            1);
  bru(ostr, "ISEN",    _ISEN,      "",                            0);
  bru(ostr, "REVERSE", REVERSES(), "",                            1);
  bru(ostr, "AUNMP",   _AUNMP,     "",                            0);
  bru(ostr, "SINO",    _SINO,      "",       bruform(_SINO,2),    1);
  bru(ostr, "NOISF1",  _NOISF1,    "ppm",    bruform(_NOISF1,3),  0);
  bru(ostr, "NOISF2",  _NOISF2,    "ppm",    bruform(_NOISF2,3),  1);
  bru(ostr, "SIGF1",   _SIGF1,     "ppm",    bruform(_SIGF1,3),   0);
  bru(ostr, "SIGF2",   _SIGF2,     "ppm",    bruform(_SIGF2,3),   1);
  bru(ostr, "ASSFAC",  _ASSFAC,    "",                            0);
  bru(ostr, "ASSFACI", _ASSFACI,   "",                            1);
  bru(ostr, "ASSFACX", _ASSFACX,   "",                            0);
  bru(ostr, "ASSWID",  _ASSWID,    "",                            1);
  bru(ostr, "DATMOD",  DATMODS(),  "",                            0);
  bru(ostr, "DC",      _DC,        "",                            1);
  bru(ostr, "NSP",     _NSP,       "",                            0);
  bru(ostr, "NZP",     _NZP,       "",                            1);
  bru(ostr, "DFILT",   _DFILT,     "",                            0);
  bru(ostr, "TI",      _TI,        "",                            1);
  bru(ostr, "TM1",     _TM1,       "",                            0);
  bru(ostr, "TM2",     _TM2,       "",                            1);
  bru(ostr, "ALPHA",   _XALPHA,    "",                            0);
  bru(ostr, "GAMMA",   _GAMMA,     "",                            1);
  bru(ostr, "YMAX_p",  _YMAX_p,    "",                            0);
  bru(ostr, "YMIN_p",  _YMIN_p,    "",                            1);
  bru(ostr, "NC_proc", _NC_proc,   "",                            0);
  bru(ostr, "MEAN",    _MEAN,      "",        bruform(_MEAN,3),   1);
  bru(ostr, "S_DEV",   _S_DEV,     "",        bruform(_S_DEV,3),  0);
  bru(ostr, "AQORDER", AQORDERS(), "",                            1);
  bru(ostr, "NLEV",    _NLEV,      "",                            0);
  bru(ostr, "LEV0",    _LEV0,      "%",       bruform(_LEV0,2),   1);
  bru(ostr, "TOPLEV",  _TOPLEV,    "%",       bruform(_TOPLEV,2), 0);
  ostr << "\n";
  return ostr;
  }

#endif							// XWinProcs.cc
