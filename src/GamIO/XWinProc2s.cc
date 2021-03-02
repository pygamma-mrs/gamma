/* XWinProc2s.cc **************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      XWinProc2s                                  Implementation	**
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
** sets. This class embodies a Bruker parameter file proc2s, that which **
** to control access to NMR processing parameters within XWinNMR. This	**
** is another ASCII file that has nothing to do with GAMMA and will	**
** only be used	for output. The file is usually dir/1/pdata/1/proc2s	**
** where dir is the directory base that contains XWinNMR data.  Note	**
** that without	the existence of a proc2 file XWinNMR complains when	**
** looking at the data, even if it is just the associated serial file. 	**
** One may still be able process the data (albeit unlikely), but then   **
** one cannot display status parameters under the "Output" menu in 	**
** XWinNMR.  Below is a typical	Bruker directory structure in a 2D	**
** (serial) acquisition case.						**
**									**
**                          __ Acqu  (changable parameter file)		**
**			   / 						**
**                        /___ acqus (static parameter file)		**
**			 /						**
**  expname -- expnum --< ---- fid (binary data)			**
**			 \						**
**			  \___ pdata -- 1 -- proc(s),proc2(s),meta,...  **
**									**
*************************************************************************/

#ifndef   XWinProc2s_CC_		// Is file already included?
#  define XWinProc2s_CC_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <GamIO/XWinProc2s.h>		// Include the interface
#include <GamIO/XWinProcPar.h>		// Include processing parameter sets
#include <GamIO/BinIOBase.h>            // Include binary IO functions
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <string>			// Include libstdc++ strings
#ifndef _MSC_VER                        // If not using MSVC++ then we
 #include <sys/time.h>			// Include time and date access
 #include <unistd.h>			// Include POSIX getcwd function
#else                                   // and if using MSVC++ then I guess
 #include <time.h>			// Include time and date access
#endif

using std::string;			// Using libstdc++ strings
using std::ostream;			// Using libstdc++ output streams

// ----------------------------------------------------------------------------
// ----------------------------- PRIVATE FUNCTIONS ----------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                      XWinNMR Proc2s File Error Handling
// ____________________________________________________________________________
 
/* These functions take care of any errors encountered when reading, writing,
   and setting parameters in Bruker acquisition parameter files.

        Input		eidx    : Error index
        		pname   : Additional error message
        		noret	: Flag for linefeed (0=linefeed)
        Output		void    : An error message is output                 */
 
void XWinProc2s::XWinPP2error(int eidx, int noret) const
  {
  string hdr("XWinNMR Proc2s Parameter File");
  string msg;
  switch (eidx)
    {
    case 25: GAMMAerror(hdr,"Cannot Write Parameters To File",noret);break; // (25)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }  

    
void XWinProc2s::XWinPP2error(int eidx, const string& pname, int noret) const
  {                                                                             
  string hdr("XWinNMR Proc2s Parameter File");
  string msg;
  switch(eidx)
    {
    case 21:msg = string("Cannot Write To ") + pname;
             GAMMAerror(hdr,msg+pname,noret);  break;                  // (21)
    default: GAMMAerror(hdr, eidx, pname, noret); break;// Unknown Error  (-1)
    }
  }

volatile void XWinProc2s::XWinPP2fatal(int eidx) const
  {                                                                 
  XWinPP2error(eidx, 1);			// Normal non-fatal error
  if(eidx) XWinPP2error(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

volatile void XWinProc2s::XWinPP2fatal(int eidx, const string& pname) const
  {                                                                 
  XWinPP2error(eidx, pname, 1);		// Normal non-fatal error
  if(eidx) XWinPP2error(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }
 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A          XWinProc2s Parameter File Constructors, Destructor
// ____________________________________________________________________________
 
/* These are constructors of the class handling Bruker XWinNMR 2D processing
   parameter files.  This doesn't do anything in particular, it is the write
   functions that perform the work.                                          */

     XWinProc2s::XWinProc2s()                   : XWinProcPar("proc2s", 2) { }
     XWinProc2s::XWinProc2s(const string& name) : XWinProcPar(name, 2)     { }
     XWinProc2s::XWinProc2s(const XWinProc2s& X): XWinProcPar(X)           { }
     XWinProc2s::~XWinProc2s()                                             { }
void XWinProc2s::operator= (const XWinProc2s& X) {XWinProcPar::operator= (X);}

// ____________________________________________________________________________
// B                  XWinProc2s Parameter Access Functions
// ____________________________________________________________________________

/* These functions allow users to set some of the basic parameters that
   these files contain.  Remember, GAMMA doesn't intend to output processed
   data sets to XWinNMR. It only generates such files because the XWinNMR
   software demands it whenever a 2D acquisition is performed - or whenever
   a simulated 2D data set is output in their format. Input is quite different, 
   as GAMMA does support the ability to read in a processed 2D spectra in
   XWinNMR format as well as access all associated parameters.
 
                             INHERITED FROM XWinProcPar
                             ==========================

string XWinProc2s::parname()   const { return parfile;  } // File name
int    XWinProc2s::BYTORDP()   const { return _BYTORDP; } // Byte order
int    XWinProc2s::FT_mod()    const { return _FT_mod;  } // FFT mode
double XWinProc2s::LB()        const { return _LB;      } // Line Broadening
int    XWinProc2s::MC2()       const { return _MC2;     } //
double XWinProc2s::OFFSET()    const { return _OFFSET;  } // Spectrum offset
double XWinProc2s::PHC0()      const { return _PHC0;    } // Zero order phase
double XWinProc2s::PHC1()      const { return _PHC1;    } // 1st order phase
string XWinProc2s::REVERSE()   const { return _REVERSE; } // Spectrum reverse
double XWinProc2s::SF()        const { return _SF;      } // Spectrom. freq.
int    XWinProc2s::SI()        const { return _SI;      } // Size (re+im)
int    XWinProc2s::SSB()       const { return _SSB;     } // Sine bell
int    XWinProc2s::STSI()      const { return _STSI;    } // Sine bell
int    XWinProc2s::STSR()      const { return _STSR;    } // Sine bell
double XWinProc2s::SW_p()      const { return _SW_p;    } // Spec. width(ppm)
double XWinProc2s::TDeff()     const { return _TDeff;   } // Effective FID size
int    XWinProc2s::WDW()       const { return _WDW;     } // Window function
 
void XWinProc2s::MC2(int mc)             { _MC2     = mc; }
void XWinProc2s::REVERSE(int yn)         { yn?_REVERSE="yes":_REVERSE="no"; }
void XWinProc2s::PPARMOD(int pm)         { _PPARMOD = pm; }
void XWinProc2s::SI(int si)              { _SI      = si; }
void XWinProc2s::SF(double SF)           { _SF      = SF; SetOffset(); }
void XWinProc2s::SSB(int sb)             { _SSB     = sb; }
void XWinProc2s::STSI(int st)            { _STSI    = st; }
void XWinProc2s::STSR(int sr)            { _STSR    = sr; }
void XWinProc2s::SW_p(double swp)        { _SW_p    = swp; SetOffset(); }
void XWinProc2s::WDW(int wd)             { _WDW     = wd; }                   */  

// ____________________________________________________________________________
// C                       XWinProc2s Input Functions
// ____________________________________________________________________________

/* These functions will read in the 2D processing parameters from an XWinNMR
   parameter file, typically named proc2s.  By design, the Bruker parameter
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
 
bool XWinProc2s::read(const string& filein, int warn)
bool XWinProc2s::read(int warn)
bool XWinProc2s::parsePSet(int warn)                                         */

// ____________________________________________________________________________
// D                         XWinProc2s Output Functions
// ____________________________________________________________________________
 
/* These functions allow for output of NMR parameters directly into a Bruker
   XWinNMR ASCII parameter file (proc2s).

                         INHERITED FROM CLASS XWinProcPar
                         ================================
 
int XWinProc2s::write(const string& dname, int warn)
int XWinProc2s::write(int warn) const                                         */


// ____________________________________________________________________________
// E                    XWinProc2s Standard Output Functions
// ____________________________________________________________________________

/* These functions allow for a quick output of the parameter file contents.
   They don't have anything to do with output while running XWinNMR, rather
   users can just glance at proc2s parameters or store then in a small file. */

// ostream& XWinProc2s::printPset(ostream& O) const		INHERITED

/*
ostream& XWinProc2s::print(ostream& ostr, int full, int hdr) const
  {
  string marg(16, ' ');
  string ls = string("\n") + marg;
  if(hdr)
    ostr << ls << "Proc2s Parameter File:       " << prname;
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

ostream& operator<< (ostream& O, const XWinProc2s& P) { P.print(O); return O; }
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
   be calculated.                                                            */

ostream& XWinProc2s::dpp(ostream& ostr, double BF0) const
  {
  if(!BF0) BF0 = _SF;
  double SRval = (_SF-BF0)*1.e6;
  double HZPPval = _SW_p/_SI;
  ostr << "\n" << string(29, ' ') << "F1 - Processing Parameters";
  ostr << "\n" << string(79, '=') << "\n\n";
  bru(ostr, "SI",      SI(),       "",                            0);
  bru(ostr, "MC2",     MC2S(),     "",                            1);
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
  bru(ostr, "TDeff",   _TDeff,     "",                            0);
  bru(ostr, "TDoff",   _TDoff,     "",                            1);
  bru(ostr, "STSR",    _STSR,      "",                            0);
  bru(ostr, "STSI",    _STSI,      "",                            1);
  bru(ostr, "BC_mod",  BC_modS(),  "",                            0);
  bru(ostr, "PH_mod",  PH_modS(),  "",                            1);
  bru(ostr, "PHC0",    _PHC0,      "degrees", bruform(_PHC0,3),   0);
  bru(ostr, "PHC1",    _PHC1,      "degrees", bruform(_PHC1,3),   1);
  bru(ostr, "ABSF1",   _ABSF1,     "ppm",     bruform(_ABSF1,3),  0);
  bru(ostr, "ABSF2",   _ABSF2,     "ppm",     bruform(_ABSF2,3),  1);
  bru(ostr, "ABSG",    _ABSG,      "",                            0);
  bru(ostr, "ABSL",    _ABSL,      "",                            1);
  bru(ostr, "FT_mod",  FT_modS(),  "",                            0);
  bru(ostr, "REVERSE", REVERSES(), "",                            1);
  bru(ostr, "SYMM",    SYMMS(),    "",                            0);
  bru(ostr, "TILT",    _TILT,      "",                            1);
  bru(ostr, "TM1",     _TM1,       "",                            0);
  bru(ostr, "TM2",     _TM2,       "",                            1);
  bru(ostr, "ME_mod",  ME_modS(),  "",                            0);
  bru(ostr, "NCOEF",   _NCOEF,     "",                            1);
  bru(ostr, "LPBIN",   _LPBIN,     "",                            0);
  ostr << "\n";
  return ostr;
  }




#endif							// XWinProc2s.cc
