/* XWinAcqu2s.cc *************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      XWinAcqu2s                                  Implementation	**
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
** sets. This class embodies a Bruker parameter file, Acqu2s, which 	**
** contains parameters detailing an acquisition.  When such a file is   **
** read in the class will contain all important acquisition parameters, **
** such as TD (# points) and SW (spectral width PPM), and make them     **
** accessible to outside functions.  When such a file is output this    **
** class will insure the file acqu2s is in Bruker XWinNMR format (ASCII)**
** When used in conjunction with other XWinNMR related GAMMA classes,   **
** reading and writing of XWinNMR 2D (and ND) serial files should be    **
** trivial.								**
**                                                                      **
** A typical Bruker directory structure in a 2D acquisiton case:	**
**									**
**                          __ acqu, acqu2 (changable parameter file)	**
**			   / 						**
**                        /___ acqus acqu2s (static parameter file)	**
**			 /						**
**  expname -- expnum --< ---- ser (binary data)			**
**			 \						**
**			  \___ pdata -- 1 -- proc, proc2, procs, proc2s **
**					     meta, outd,  2rr, 2ii      **
**									**
** Note: Of critical importance to these data sets is the byte order	**
**       in which the binary data is stored.  The order is specified	**
**       in the ASCII parameter file acqu2s as BYTORDA. If its value 	**
**       is zero then the data is little endian (e.g. written from an	**
**       Intel based machine).  If its value is 1 then the data is	**
**       big endian (e.g. written on an SGI or SPARC). 			**
**									**
*************************************************************************/

#ifndef   XWinAcqu2s_CC_		// Is file already included?
#  define XWinAcqu2s_CC_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <GamIO/XWinAcqu2s.h>		// Include the interface
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
// i                      XWinNMR Acqu2s File Error Handling
// ____________________________________________________________________________
 
/* These functions take care of any errors encountered when reading, writing,
   and setting parameters in Bruker acquisition parameter files.

        Input		eidx    : Error index
        		pname   : Additional error message
        		noret	: Flag for linefeed (0=linefeed)
        Output		void    : An error message is output                 */
 
void XWinAcqu2s::XWinAcqu2serror(int eidx, int noret) const
  {
  string hdr("XWinNMR Acqu2s Parameter File");
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

void XWinAcqu2s::XWinAcqu2serror(int eidx, const string& pname, int noret) const
  {                                                                             
  string hdr("XWinNMR 2D Acquisition Parameter File");
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

volatile void XWinAcqu2s::XWinAcqu2sfatality(int eidx) const
  {                                                                 
  XWinAcqu2serror(eidx, 1);			// Normal non-fatal error
  if(eidx) XWinAcqu2serror(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

volatile void XWinAcqu2s::XWinAcqu2sfatality(int eidx, const string& pname) const
  {                                                                 
  XWinAcqu2serror(eidx, pname, 1);		// Normal non-fatal error
  if(eidx) XWinAcqu2serror(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }
 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A             XWinAcqu2s Parameter File Constructors, Destructor
// ____________________________________________________________________________

/* These are the constructors of the class handling Bruker XWinNMR acquisition
   parameter files.  This doesn't do anything in particular, it is the read
   and write functions that perform the work.  Thus, we only have a default
   constructor specified.  The reading and writing of the associated ASCII
   parameter file is done in one step, so we don't need anything complex.    */

     XWinAcqu2s::XWinAcqu2s()                      :XWinAcqPar("acqu2s", 2) { }
     XWinAcqu2s::XWinAcqu2s(const string& name)    :XWinAcqPar(name, 2)     { }
     XWinAcqu2s::XWinAcqu2s(const XWinAcqu2s& XWA) :XWinAcqPar(XWA)         { }
     XWinAcqu2s::~XWinAcqu2s()                                              { }
void XWinAcqu2s::operator= (const XWinAcqu2s& XWA){XWinAcqPar::operator=(XWA);}

// ____________________________________________________________________________
// B                  XWinAcqu2s Parameter Access Functions
// ____________________________________________________________________________

/* These functions allow direct access to some of the more important parameters
   that each Acqu2sition file should know.  The two primary values are the fid
   point size (TD) and the data byte order (BYTORDA).

                        INHERITED FROM CLASS XWinAcqPar
                        ===============================

int    XWinAcqu2s::filename() const { return _parfile;           }
double XWinAcqu2s::field()    const { return _Bo;                }
double XWinAcqu2s::AQ()       const { return _TD/(2*_SW*_SFO1);  }
int    XWinAcqu2s::AQ_mod()   const { return _AQ_mod;            }
double XWinAcqu2s::BF1()      const { return _BF[0];             }
double XWinAcqu2s::BF2()      const { return _BF[1];             }
int    XWinAcqu2s::BYTORDA()  const { return _BYTORDA            }
int    XWinAcqu2s::DS()       const { return _DS;                }
string XWinAcqu2s::EXP()      const { return _EXP;               }
double XWinAcqu2s::IN(int i)  const { return _IN[i];             }
string XWinAcqu2s::NAME()     const { return _NAME;              }
int    XWinAcqu2s::NS()       const { return _NS;                }
string XWinAcqu2s::NUC(int i) const { return _NUC[i];            }
string XWinAcqu2s::NUCLEUS()  const { return _NUCLEUS;           }
double XWinAcqu2s::O1()       const { return _O[0];              }
double XWinAcqu2s::O2()       const { return _O[1];              }
int    XWinAcqu2s::PARMODE()  const { return _PARMODE;           }
string XWinAcqu2s::PULPROG()  const { return _PULPROG;           }
double XWinAcqu2s::SFO1()     const { return _SFO[0];            }
double XWinAcqu2s::SFO2()     const { return _SFO[1];            }
double XWinAcqu2s::SFO3()     const { return _SFO[2];            }
string XWinAcqu2s::SOLVENT()  const { return _SOLVENT;           }
double XWinAcqu2s::SW()       const { return _SW;                }
double XWinAcqu2s::SW_h()     const { return _SW/_SFO[0];        }
double XWinAcqu2s::TE()       const { return _TE;                }
int    XWinAcqu2s::TD()       const { return _TD;                }

void XWinAcqu2s::field(double BoT)         { Bo = fabs(BoT);         }
void XWinAcqu2s::AQ_mod(int aqmo)          { _AQ_mod = aqmo;         }
int  XWinAcqu2s::D(int i, double t, int w) { return SetDelay(i,t,w); }
void XWinAcqu2s::DS(int ds)                { _DS     = ds;           }
void XWinAcqu2s::EXP(const string& exp)    { _EXP    = exp;          }
void XWinAcqu2s::NS(int ns)                { _NS     = ns;           }
void XWinAcqu2s::NUCLEUS(const string& I)  { _NUCLEUS = I;           }
void XWinAcqu2s::NUCLEI(int CH, const string& I, double OF, int warn)
int  XWinAcqu2s::P(int i, double t, int w) { return SetPulse(i,t,w); }
void XWinAcqu2s::PULPROG(const string& P)  { _PULPROG = P;           }
void XWinAcqu2s::SFO1(double sf)           { _SFO[0]  = sf;          }
void XWinAcqu2s::SFO2(double sf)           { _SFO[1]  = sf;          }
void XWinAcqu2s::SFO3(double sf)           { _SFO[2]  = sf;          }
void XWinAcqu2s::SFO(double sf, int i)     { _SFO[i]  = sf;          }
void XWinAcqu2s::SOLVENT(const string& S)  { _SOLVENT = S;           }
void XWinAcqu2s::SW_h(double sw)           { _SW      = sw/_SFO[0];  }
void XWinAcqu2s::SW(double sw)             { _SW     = sw;           }
void XWinAcqu2s::SWHZ(double sw)           { _SW     = sw/_SFO[0];   }
void XWinAcqu2s::TD(int td)                { _TD     = td;           }
void XWinAcqu2s::TE(double te)             { _TE     = te;           }       */

// ____________________________________________________________________________
// C                       XWinAcqu2s Input Functions
// ____________________________________________________________________________

/* These functions will read in the acquisition parameters from an XWinNMR
   parameter file, typically named acqu2s.  By design, the Bruker parameter
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
 
bool XWinAcqu2s::readAPar(const string& filein, int warn)
bool XWinAcqu2s::readAPar(int warn)
bool XWinAcqu2s::parsePSet(int warn)                                         */

// ____________________________________________________________________________
// D                         XWinAcqu2s Output Functions
// ____________________________________________________________________________
 
/* These function allow for output of NMR parameters directly into a Bruker
   XWinNMR ASCII parameter file (Acqu2s).  We don't write everything that
   Bruker does since users would have to set all of them if so.

                     INHERITED FROM CLASS XWinAcqPar

int XWinAcqu2s::writeAPar(const string& name, int warn)
int XWinAcqu2s::writeAPar(int warn=2) const                                  */

// ____________________________________________________________________________
// E                    XWinAcqu2s Standard Output Functions
// ____________________________________________________________________________

/* These function allow for a quick output of the parameter file contents.
   They don't have anything to do with output while running XWinNMR, rather
   users can just glance at Acqu2s parameters or store then in a small file.  */

// ostream& XWinAcqu2s::printPset(ostream& O) const                INHERITED

ostream& XWinAcqu2s::print(ostream& ostr, int full, int hdr) const
  {
  string marg(16, ' ');
  string ls = string("\n") + marg;
  if(hdr)
    ostr << ls << "Acquisition Parameter File: " << parfile;
  string byord = "little";			// Byte order label
  if(_BYTORDA) byord = string("big");
  ostr << ls << "Data Byte Ordering";
  ostr << ":         ";
  ostr << byord             << " endian";
  ostr << ls << "Number of Blocks:           " << _TD;
  ostr << ls << "Spectrometer Frequency:     " << _SFO[0]           << " MHz";
  ostr << ls << "Spectral Width:             " << _SW               << " PPM";
  ostr << ls << "                            " << _SW*_SFO[0]       << " Hz";
  ostr << ls << "Dwell Time:                 " << 1./(_SW*_SFO[0])  << " sec";
  ostr << ls << "FID Length:                 " << _TD/(_SW*_SFO[0]) << " sec";
  return ostr;
  }  

ostream& operator<< (ostream& O, const XWinAcqu2s& A) { A.print(O); return O; }

// ____________________________________________________________________________
// F                    XWinAcqu2s BrukerLike Output Functions
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


ostream& XWinAcqu2s::dpa(ostream& ostr) const
  {
  ostr << "\n" << string(27, ' ') << "F1 - Acquisition Parameters";
  ostr << "\n" << string(79, '=') << "\n\n";
  bru(ostr, "NDO",    _DR,           "",                           0);
  bru(ostr, "TD",     TD(),          "",                           1);
  bru(ostr, "SFO1",   _SFO[0],       "MHz",    bruform(_SFO[0],7), 0);
  bru(ostr, "FIDRES", FIDRES(),      "Hz",     bruform(FIDRES(),6),1);
  bru(ostr, "SW",     _SW,           "ppm",    bruform(_SW,4),     0); 
  bru(ostr, "SWH",    SW_h(),        "Hz",     bruform(SW_h(),3),  1);
  bru(ostr, "FOV",    _FOV,          "cm",     bruform(_FOV,2),    0);
  ostr << "\n";
  return ostr;
  }

// ____________________________________________________________________________
// G                    XWinAcqu2s GAMMA Output Functions
// ____________________________________________________________________________
 
/* This function is used to quickly tell GAMMA the important values in Acqu2s */

ostream& XWinAcqu2s::printGAMMA(ostream& ostr) const

  {
  ostr << "\n\n\t\tXWinNMR 1D Acqu2s Parameters For GAMMA\n";
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

#endif							// XWinAcqu2s.cc
