/* XWinAcqus.cc *************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      XWinAcqus                                  Implementation	**
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
** sets. This class embodies a Bruker parameter file, acqus, which will **
** contain parameters detailing an acquisition.  When such a file is    **
** read in the class will contain all important acquisition parameters, **
** such as TD (# points) and SW (spectral width PPM), and make them     **
** accessible to outside functions.  When such a file is output this    **
** class will insure the file acqus is in Bruker XWinNMR format (ASCII) **
** When used in conjunction with other XWinNMR related GAMMA classes,   **
** reading and writing of XWinNMR FID's and serial files should be      **
** trivial.								**
**                                                                      **
** A typical Bruker directory structure in a 1D acquisiton case:	**
**									**
**                          __ acqu  (changable parameter file)		**
**			   / 						**
**                        /___ acqus (static parameter file)		**
**			 /						**
**  expname -- expnum --< ---- fid (binary data)			**
**			 \						**
**			  \___ pdata -- 1 -- proc, procs, meta		**
**									**
** Note: Of critical importance to these data sets is the byte order	**
**       in which the binary data is stored.  The order is specified	**
**       in the ASCII parameter file acqus as BYTORDA. If its value 	**
**       is zero then the data is little endian (e.g. written from an	**
**       Intel based machine).  If its value is 1 then the data is	**
**       big endian (e.g. written on an SGI or SPARC). Here is the fun	**
**       part: Bruker leaves this parameter out of the acqus file	**
**       that part of their "exam1d" data set.  One might assume then	**
**       that the data is big endian, say written from their beloved	**
**       SGI based consoles.  But NO, it is little endian.  Thus I will	**
**       herein assume that if there is no BYTORDA specified in acqus	**
**       then the data is little endian.				**
**									**
*************************************************************************/

#ifndef _XWinAcqus_CC_			// Is file already included?
#  define _XWinAcqus_CC_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <GamIO/XWinAcqus.h>		// Include the interface
#include <GamIO/XWinPSet.h>		// Include Bruker parameter parse
#include <Basics/ParamSet.h>		// Include parameters sets
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <Basics/Gconstants.h>		// Include GAMMA constants (HZ2RAD)
#include <Basics/StringCut.h>		// Include GAMMA string parsing
#include <Basics/Isotope.h>		// Include GAMMA spin isotopes
#include <GamIO/BinIOBase.h>		// Include binary IO functions
#include <string>			// Include libstdc++ strings
#ifndef _MSC_VER                        // If not using MSVC++ then we
 #include <sys/time.h>			// Include time and date access
 #include <unistd.h>			// Include POSIX getcwd function
#else                                   // and if using MSVC++ then I guess
 #include <ctime>			// Include time and date access
#endif

using std::string;			// Using libstdc++ strings
using std::ostream;			// Using libstdc++ output streams

// ----------------------------------------------------------------------------
// ----------------------------- PRIVATE FUNCTIONS ----------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                      XWinNMR Acqus File Error Handling
// ____________________________________________________________________________
 
/* These functions take care of any errors encountered when reading, writing,
   and setting parameters in Bruker acquisition parameter files.

        Input		eidx    : Error index
        		pname   : Additional error message
        		noret	: Flag for linefeed (0=linefeed)
        Output		void    : An error message is output                 */
 
void XWinAcqus::XWinAcquserror(int eidx, int noret) const
  {
  string hdr("XWinNMR 1D Acqus Parameter File");
  switch (eidx)
    {
    case 19: GAMMAerror(hdr,"Cannot Get Acquisition Params", noret); break; // (19)
    case 20: GAMMAerror(hdr,"Cannot Parse Bruker Parameters",noret); break; // (20)
    case 21: GAMMAerror(hdr,"Cannot Determine 1D Data Size",noret);  break; // (21)
    case 23: GAMMAerror(hdr,"Cannot Find Spectrometer Freq.",noret); break; // (23)
    case 24: GAMMAerror(hdr,"Cannot Determine Spectral Width",noret);break; // (24)
    case 25: GAMMAerror(hdr,"Cannot Write Parameters To File",noret);break; // (25)
    case 26: GAMMAerror(hdr,"Delay Time Index Must Be [0,31]",noret);break; // (26)
    case 27: GAMMAerror(hdr,"Cannot Have Negative Delay Time",noret);break; // (27)
    case 28: GAMMAerror(hdr,"Pulse Length Index Outside [0,31]",noret);break;//(28)
    case 29: GAMMAerror(hdr,"Cannot Have Negative Pulse Length",noret);break;//(29)
    case 30: GAMMAerror(hdr,"Channels Must Be Labeled [1,8]", noret);break; // (30)
    case 31: GAMMAerror(hdr,"No Spectrometer Field Specified",noret);break; // (31)
    case 32: GAMMAerror(hdr,"Set Spectrometer Field Using B0",noret);break; // (32)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }  

void XWinAcqus::XWinAcquserror(int eidx, const string& pname, int noret) const
  {                                                                             
  string hdr("XWinNMR 1D Acquisition Parameter File");
  string msg;
  switch(eidx)
    {
    case 20:msg = string("Cannot Find ");
             GAMMAerror(hdr,msg+pname,noret);  break;                  // (20)
    case 21:msg = string("Cannot Write To ");
             GAMMAerror(hdr,msg+pname,noret);  break;                  // (21)
    case 22:msg = string("Cannot Set Parameter ");
             GAMMAerror(hdr,msg+pname,noret);  break;                  // (22)
    case 23:msg = string("Cannot Set Channel For Nucleus ");
             GAMMAerror(hdr,msg+pname,noret);  break;                  // (23)
    case 24:msg = string("Problems With Unknown Nucleus ");
             GAMMAerror(hdr,msg+pname,noret);  break;                  // (24)
    default: GAMMAerror(hdr, eidx, pname, noret); break;// Unknown Error  (-1)
    }
  }

volatile void XWinAcqus::XWinAcqusfatality(int eidx) const
  {                                                                 
  XWinAcquserror(eidx, 1);			// Normal non-fatal error
  if(eidx) XWinAcquserror(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

volatile void XWinAcqus::XWinAcqusfatality(int eidx, const string& pname) const
  {                                                                 
  XWinAcquserror(eidx, pname, 1);		// Normal non-fatal error
  if(eidx) XWinAcquserror(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }
 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A             XWinAcqus Parameter File Constructors, Destructor
// ____________________________________________________________________________

/* These are the constructors of the class handling Bruker XWinNMR acquisition
   parameter files.  This doesn't do anything in particular, it is the read
   and write functions that perform the work.  Thus, we only have a default
   constructor specified.  The reading and writing of the associated ASCII
   parameter file is done in one step, so we don't need anything complex.    */

     XWinAcqus::XWinAcqus()                     : XWinAcqPar("acqus", 1) { }
     XWinAcqus::XWinAcqus(const string& name)   : XWinAcqPar(name, 1)    { }
     XWinAcqus::XWinAcqus(const XWinAcqus& XWA) : XWinAcqPar(XWA)        { }
     XWinAcqus::~XWinAcqus()                                             { }
void XWinAcqus::operator= (const XWinAcqus& XWA) {XWinAcqPar::operator=(XWA);}

// ____________________________________________________________________________
// B                  XWinAcqus Parameter Access Functions
// ____________________________________________________________________________

/* These functions allow direct access to some of the more important parameters
   that each acqusition file should know.  The two primary values are the fid
   point size (TD) and the data byte order (BYTORDA).
 
                        INHERITED FROM CLASS XWinAcqPar
                        ===============================
 
double XWinAcqus::field()     const { return _Bo;               }
int    XWinAcqus::acqname()   const { return parfile;           }
double XWinAcqus::AQ()        const { return _TD/(2*_SW*_SFO1); }
int    XWinAcqus::AQ_mod()    const { return _AQ_mod;           }
double XWinAcqus::BF1()       const { return _BF[0];            }
double XWinAcqus::BF2()       const { return _BF[1];            }
int    XWinAcqus::BYTORDA()   const { return _BYTORDA           }
int    XWinAcqus::DSc()        const { return _DSc;               } 
string XWinAcqus::EXP()       const { return _EXP;              } 
double XWinAcqus::IN(int i)   const { return _IN[i];            }
string XWinAcqus::NAME()      const { return _NAME;             }
int    XWinAcqus::NS()        const { return _NS;               } 
string XWinAcqus::NUC(int i)  const { return _NUC[i];           }
string XWinAcqus::NUCLEUS()   const { return _NUCLEUS;          }
double XWinAcqus::O1()        const { return _O[0];             }
double XWinAcqus::O2()        const { return _O[1];             }
int    XWinAcqus::PARMODE()   const { return _PARMODE;          }
string XWinAcqus::PULPROG()   const { return _PULPROG;          }
double XWinAcqus::SFO1()      const { return _SFO[0];           }
double XWinAcqus::SFO2()      const { return _SFO[1];           }
double XWinAcqus::SFO3()      const { return _SFO[2];           }
string XWinAcqus::SOLVENT()   const { return _SOLVENT;          }
double XWinAcqus::SW()        const { return _SW;               }
double XWinAcqus::SW_h()      const { return _SW/_SFO[0];       }
double XWinAcqus::TE()        const { return _TE;               }
int    XWinAcqus::TD()        const { return _TD;               }

void XWinAcqus::field(double Bo)          { _Bo      = Bo;          }
void XWinAcqus::AQ_mod(int aqmo)          { _AQ_mod  = aqmo;        }   
int  XWinAcqus::D(int i, double t, int w) { return SetDelay(i,t,w); }
void XWinAcqus::DSc(int ds)                { _DSc      = ds;          }   
void XWinAcqus::EXP(const string& exp)    { _EXP     = exp;         }   
void XWinAcqus::NAME(const string& nm)    { _NAME    = fulldir+nm;  }
void XWinAcqus::NS(int ns)                { _NS      = ns;          }   
void XWinAcqus::NUCLEUS(const string& I)  { _NUCLEUS = I;           }   
void XWinAcqus::NUCLEI(int CH, const string& I, double OF, int warn)
int  XWinAcqus::P(int i, double t, int w) { return SetPulse(i,t,w); }
void XWinAcqus::PULPROG(const string& P)  { _PULPROG = P;           }   
void XWinAcqus::SFO1(double sf)           { _SFO[0]  = sf;          }   
void XWinAcqus::SFO2(double sf)           { _SFO[1]  = sf;          }   
void XWinAcqus::SFO3(double sf)           { _SFO[2]  = sf;          }   
void XWinAcqus::SFO(double sf, int i)     { _SFO[i]  = sf;          }   
void XWinAcqus::SOLVENT(const string& S)  { _SOLVENT = S;           }   
void XWinAcqus::SW(double sw)             { _SW      = sw;          }   
void XWinAcqus::SW_h(double sw)           { _SW      = sw/_SFO[0];  }
void XWinAcqus::TD(int td)                { _TD      = td;          }   
void XWinAcqus::TE(double te)             { _TE      = te;          }        */

// ____________________________________________________________________________
// C                       XWinAcqus Input Functions
// ____________________________________________________________________________

/* These functions will read in the acquisition parameters from an XWinNMR
   parameter file, typically named acqus.  By design, the Bruker parameter
   file is initially read into a GAMMA parameter set so that ALL parameters in
   the file are stored.  Subsequently, the parameters in the parameter set are
   parsed to obtain values of consequence to GAMMA and these are explicitly
   maintained variables in the class.

      Function                               Purpose
   ____________________________________________________________________________

     parsePSet          Parses a Bruker parameter set (XWinPSet) for parameters
                        of consequence to an acquisition.
       read             Read in parameter set (class XWinPset) and parse out
                        parameters of consequence to an acquisition.

bool XWinAcqus::readAPar(const string& filein, int warn)
bool XWinAcqus::readAPar(int warn)
bool XWinAcqus::parsePSet(int warn) 					     */
 
// ____________________________________________________________________________
// D                         XWinAcqus Output Functions
// ____________________________________________________________________________
 
/* These function allow for output of NMR parameters directly into a Bruker
   XWinNMR ASCII parameter file (acqus).  We don't write everything that
   Bruker does since users would have to set all of them if so.

                     INHERITED FROM CLASS XWinAcqPar

int XWinAcqus::writeAPar(const string& name, int warn)
int XWinAcqus::writeAPar(int warn=2) const                                   */

// ____________________________________________________________________________
// E                    XWinAcqus Standard Output Functions
// ____________________________________________________________________________

/* These function allow for a quick output of the parameter file contents.
   They don't have anything to do with output while running XWinNMR, rather
   users can just glance at acqus parameters or store then in a small file.  */

// ostream& XWinAcqus::printPset(ostream& O) const                INHERITED

ostream& XWinAcqus::print(ostream& ostr, int full, int hdr) const
  {
  string marg(16, ' ');
  string ls = string("\n") + marg;
  if(hdr)
    ostr << ls << "Acquisition Parameter File: " << parfile;
  string byord = "little";			// Byte order label
  if(_BYTORDA) byord = string("big");
  ostr << ls << "Data Byte Ordering:         ";
  ostr << byord             << " endian";
  ostr << ls << "Number of Points (Re+Im):   " << _TD;
  ostr << ls << "Spectrometer Frequency:     " << _SFO[0]           << " MHz";
  ostr << ls << "Spectral Width:             " << _SW               << " PPM";
  ostr << ls << "                            " << _SW*_SFO[0]       << " Hz";
  ostr << ls << "Dwell Time:                 " << 1./(_SW*_SFO[0])  << " sec";
  ostr << ls << "FID Length:                 " << _TD/(_SW*_SFO[0]) << " sec";
  return ostr;
  }  

ostream& operator<< (ostream& O, const XWinAcqus& A) { A.print(O); return O; }

// ____________________________________________________________________________
// F                    XWinAcqus BrukerLike Output Functions
// ____________________________________________________________________________
 
/* These functions perform ASCII output of acquisiton parameters.  The format
   here is similar to issuing a "dpa" command in XWinNMR or using the menu
   choice "Output/Display status pars./Acquisition only".  It isn't as complete
   as the Bruker output because I'm bored and its often a mystery how they
   get some of their values.....  Most of the values output are taken from the
   parameter set which holds all the values which were read in from the
   Bruker parameter file.  It makes no sense to use these functions if on
   hasn't first read such a file since all the parameters will be empty.  A
   couple of values are used directly from the class when parameters need to
   be calculated.						             */


ostream& XWinAcqus::dpa(ostream& ostr) const
  {
  ostr << "\n" << string(29, ' ') << "Acquisition Parameters";
  ostr << "\n" << string(79, '=') << "\n\n";
	// *** _DATE has been redefined as a time_t.  
	// Cast as long on output to make sure no loss of data.
  bru(ostr, "Date_",    static_cast<long>(_DATE),         "",        0);
  bru(ostr, "Time",     "NOCH NICHT",         "",                    1);
// int tf = pset.getString("Time",sval2);
//  bru(ostr, "Time",     _TIME,         "",                           1);
  bru(ostr, "PULPROG",  _PULPROG,      "",                           0);
  bru(ostr, "AQ_mod",   AQ_modS(),     "",                           1);
  bru(ostr, "TD",       TD(),          "",                           0);
  bru(ostr, "PARMODE",  PARMODES(),    "",                           1);
  bru(ostr, "NS",       _NS,           "",                           0);
  bru(ostr, "DS",       _DSc,           "",                           1);
  bru(ostr, "D",        "** Array **", "sec",                        0);
  bru(ostr, "P",        "** Array **", "usec",                       1);
  bru(ostr, "SW",       _SW,           "ppm",    bruform(_SW,4),     0);
  bru(ostr, "SWH",      SW_h(),        "Hz",     bruform(SW_h(),3),  1);
  bru(ostr, "FIDRES",   FIDRES(),      "Hz",     bruform(FIDRES(),6),0);
  bru(ostr, "FW",       _FW,           "Hz",     bruform(_SW,2),     1);
  bru(ostr, "AQ",       AQ(),          "sec",    bruform(_SW,7),     0);
  bru(ostr, "RG",       _RG,           "",                           1);
  bru(ostr, "DW",       DW(),          "usec",   bruform(DW(),3),    0);
  bru(ostr, "DIGTYP",   DIGTYPS(),     "",                           1);
  bru(ostr, "DE",       _DE,           "usec",   bruform(_DE,2),     0);
  bru(ostr, "DR",       _DR,           "",                           1);
  bru(ostr, "PHP",      _PHP,          "",                           0);
  bru(ostr, "QNP",      _QNP,          "",                           1);
  bru(ostr, "HL1",      _HL1,          "db",                         0);
  bru(ostr, "HL2",      _HL2,          "db",                         1);
  bru(ostr, "HL3",      _HL3,          "db",                         0);
  bru(ostr, "HL4",      _HL4,          "db",                         1);
  bru(ostr, "",         "",            "",                           0);
  bru(ostr, "NUCLEUS",  _NUCLEUS,      "",                           1);
  bru(ostr, "O1",       _O[0],         "Hz",     bruform(_O[0],2),   0);
  bru(ostr, "SFO1",     _SFO[0],       "MHz",    bruform(_SFO[0],7), 1);
  bru(ostr, "",         "",            "",                           0);
  bru(ostr, "BF1",      _BF[0],        "MHz",    bruform(_BF[0],7),  1);
  bru(ostr, "",         "",            "",                           0);
  bru(ostr, "DECNUC",   _DECNUC,       "",                           1);
  bru(ostr, "O2",       _O[1],         "Hz",     bruform(_O[1],2),   0);
  bru(ostr, "SFO2",     _SFO[1],       "MHz",    bruform(_SFO[1],7), 1);
  bru(ostr, "",         "",            "",                           0);
  bru(ostr, "BF2",      _BF[1],        "MHz",    bruform(_BF[1],7),  1);
  bru(ostr, "",         "",            "",                           0);
  bru(ostr, "DECBNUC",  _DECBNUC,      "",                           1);
  bru(ostr, "O3",       _O[2],         "Hz",     bruform(_O[2],2),   0);
  bru(ostr, "SFO3",     _SFO[2],       "MHz",    bruform(_SFO[2],7), 1);
  bru(ostr, "",         "",            "",                           0);
  bru(ostr, "BF3",      _BF[2],        "MHz",    bruform(_BF[2],7),  1);
  bru(ostr, "O4",       _O[3],         "",       bruform(_O[3],2),   0);
  bru(ostr, "SFO4",     _SFO[3],       "Hz",     bruform(_SFO[3],2), 1);
  bru(ostr, "",         "",            "",                           0);
  bru(ostr, "BF4",      _BF[3],        "MHz",    bruform(_BF[3],7),  1);
  bru(ostr, "TL",       "** Array **", "dB",                         0);
  bru(ostr, "DL",       "** Array **", "dB",                         1);
  bru(ostr, "DBL",      "** Array **", "dB",                         0);
  bru(ostr, "CNST",     "** Array **", "",                           1);
  bru(ostr, "CPDPRG",   _CPDPRG,       "",                           0);
  bru(ostr, "CPDPRGT",  _CPDPRGT,      "",                           1);
  bru(ostr, "CPDPRGB",  _CPDPRGB,      "",                           0);
  bru(ostr, "CPDPRG4",  _CPDPRGN[3],   "",                           1);
  bru(ostr, "GRDPROG",  _GRDPROG,      "",                           0);
  bru(ostr, "LOCNUC",   _LOCNUC,       "",                           1);
  bru(ostr, "SOLVENT",  _SOLVENT,      "",                           0);
  bru(ostr, "PROBHD",   _PROBHD,       "",                           1);
  bru(ostr, "",         "",            "",                           0);
  bru(ostr, "INSTRUM",  _INSTRUM,      "",                           1);
  bru(ostr, "RO",       _RO,           "Hz",                         0);
  bru(ostr, "TE",       _TE,           "K",      bruform(_TE,1),     1);
  bru(ostr, "NBL",      _NBL,          "",                           0);
  bru(ostr, "V9",       _V9,           "%",      bruform(_V9,2),     1);
  bru(ostr, "WBST",     _WBST,         "",                           0);
  bru(ostr, "WBSW",     _WBSW,         "MHZ",    bruform(_WBSW,7),   1);
  bru(ostr, "AUNM",     _AUNM,         "",                           0);
  bru(ostr, "POWMOD",   POWMODS(),     "",                           1);
  bru(ostr, "HPPRGN",   HPPRGNS(),     "",                           0);
  bru(ostr, "PRGAIN",   PRGAINS(),     "",                           1);
  bru(ostr, "IN",       "** Array **", "sec",                        0);
  bru(ostr, "INP",      "** Array **", "used",                       1);
  bru(ostr, "L",        "** Array **", "",                           0);
  bru(ostr, "SEOUT",    SEOUTS(),      "",                           1);
  bru(ostr, "S",        "** Array **", "dB",                         0);
  bru(ostr, "FS",       "** Array **", "dB",                         1);
  bru(ostr, "PH_ref",   _PH_ref,       "degree", bruform(_PH_ref,3), 0);
  bru(ostr, "PHCOR",    "** Array **", "degree",                     1);
  bru(ostr, "ROUTWD1",  "** Array **", "",                           0);
  bru(ostr, "ROUTWD2",  "** Array **", "",                           1);
  bru(ostr, "TP",       "** Array **", "dB",                         0);
  bru(ostr, "TPOFFS",   "** Array **", "Hz",                         1);
  bru(ostr, "TPNAME0",  _TPNAME[0],    "",                           0);
  bru(ostr, "TPNAME1",  _TPNAME[1],    "",                           1);
  bru(ostr, "TPNAME2",  _TPNAME[2],    "",                           0);
  bru(ostr, "TPNAME3",  _TPNAME[3],    "",                           1);
  bru(ostr, "TPNAME4",  _TPNAME[4],    "",                           0);
  bru(ostr, "TPNAME5",  _TPNAME[5],    "",                           1);
  bru(ostr, "TPNAME6",  _TPNAME[6],    "",                           0);
  bru(ostr, "TPNAME7",  _TPNAME[7],    "",                           1);
  bru(ostr, "DP",       "** Array **", "dB",                         0);
  bru(ostr, "DPOFFS",   "** Array **", "Hz",                         1);
  bru(ostr, "DPNAME0",  _DPNAME[0],    "",                           0);
  bru(ostr, "DPNAME1",  _DPNAME[1],    "",                           1);
  bru(ostr, "DPNAME2",  _DPNAME[2],    "",                           0);
  bru(ostr, "DPNAME3",  _DPNAME[3],    "",                           1);
  bru(ostr, "DPNAME4",  _DPNAME[4],    "",                           0);
  bru(ostr, "DPNAME5",  _DPNAME[5],    "",                           1);
  bru(ostr, "DPNAME6",  _DPNAME[6],    "",                           0);
  bru(ostr, "DPNAME7",  _DPNAME[7],    "",                           1);
  bru(ostr, "DBP",      "** Array **", "dB",                         0);
  bru(ostr, "DBPOFFS",  "** Array **", "Hz",                         1);
  bru(ostr, "DBPNAM0",  _DBPNAM[0],    "",                           0);
  bru(ostr, "DBPNAM1",  _DBPNAM[1],    "",                           1);
  bru(ostr, "DBPNAM2",  _DBPNAM[2],    "",                           0);
  bru(ostr, "DBPNAM3",  _DBPNAM[3],    "",                           1);
  bru(ostr, "DBPNAM4",  _DBPNAM[4],    "",                           0);
  bru(ostr, "DBPNAM5",  _DBPNAM[5],    "",                           1);
  bru(ostr, "DBPNAM6",  _DBPNAM[5],    "",                           0);
  bru(ostr, "DBPNAM7",  _DBPNAM[7],    "",                           1);
  bru(ostr, "FL1",      _FL1,          "dB",                         0);
  bru(ostr, "FL2",      _FL2,          "dB",                         1);
  bru(ostr, "FL3",      _FL3,          "dB",                         0);
  bru(ostr, "FL4",      _FL4,          "dB",                         1);
  bru(ostr, "XL",       _XL,           "",                           0);
  bru(ostr, "YL",       _YL,           "",                           1);
  bru(ostr, "ZL1",      _ZL1,          "dB",     bruform(_ZL1,2),    0);
  bru(ostr, "ZL2",      _ZL2,          "dB",     bruform(_ZL2,2),    1);
  bru(ostr, "ZL3",      _ZL3,          "dB",     bruform(_ZL3,2),    0);
  bru(ostr, "ZL4",      _ZL4,          "dB",     bruform(_ZL4,2),    1);
  bru(ostr, "DSLIST",   _DSLIST,       "",                           0);
  bru(ostr, "F1LIST",   _F1LIST,       "",                           1);
  bru(ostr, "F2LIST",   _F2LIST,       "",                           0);
  bru(ostr, "F3LIST",   _F3LIST,       "",                           1);
  bru(ostr, "VCLIST",   _VCLIST,       "",                           0);
  bru(ostr, "VDLIST",   _VDLIST,       "",                           1);
  bru(ostr, "VPLIST",   _VPLIST,       "",                           0);
  bru(ostr, "VTLIST",   _VTLIST,       "",                           1);
  bru(ostr, "FOV",      _FOV,          "cm",                         0);
  bru(ostr, "DECSTAT",  DECSTATS(),    "",                           1);
  bru(ostr, "YMAX_a",   _YMAX_a,       "",                           0);
  bru(ostr, "YMIN_a",   _YMIN_a,       "",                           1);
  bru(ostr, "NC",       _NC,           "",                           0);
  bru(ostr, "RD",       _RD,           "sec",    bruform(_RD,7),     1);
  bru(ostr, "VD",       _VD,           "sec",    bruform(_VD,8),     0);
  bru(ostr, "PW",       _PW,           "usec",   bruform(_PW,1),     1);
  bru(ostr, "PAPS",     PAPSS(),        "",                          0);
  ostr << "\n";
  return ostr;
  }


// ____________________________________________________________________________
// G                    XWinAcqus GAMMA Output Functions
// ____________________________________________________________________________
 
/* This function is used to quickly tell GAMMA the important values in acqus */

ostream& XWinAcqus::printGAMMA(ostream& ostr) const

  {
  ostr << "\n\n\t\tXWinNMR 1D acqus Parameters For GAMMA\n";
  ostr << "\n\t\tFile Name:                " << parfile;
  ostr << "\n\t\tFile Byte Ordering:       ";
  if(!_BYTORDA) ostr << "Little Endian";
  else          ostr << "Big Endian";
  ostr << "\n\t\tBlock Size (TD/2):        " << _TD/2       << " Pts";
  ostr << "\n\t\tOmega (SFO1):             " << _SFO[0]     << " MHz";
  ostr << "\n\t\tSweep Width (SW):         " << _SW         << " PPM";
  ostr << "\n\t\t                          " << _SW*_SFO[0] << " Hz";
  ostr << "\n\t\tDwell Time (DW)           " << 1.e6/(_SW*_SFO[0])  << " usec";
  ostr << "\n\t\tAcquisition Length (AQ)   " << AQ()        << " sec";
  ostr << "\n\t\tAcquisition Mode(AQ_mod): " << _AQ_mod;
  ostr << "\n\t\tBase Nucleus (_NUCLEUS):  " << _NUCLEUS;
  return ostr;
  }

#endif							// XWinAcqus.cc
