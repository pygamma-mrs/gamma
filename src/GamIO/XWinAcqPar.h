/* XWinAcqPar *************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      XWinAcqPar                                   Interface		**
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
** sets. This class embodies a "set" of Bruker acquisition parameters,	**
** parameters that might be found in their ASCII parameter files for	**
** an acquisiton: acqu, acqus, acqu2, acqu2s, ....  It is really a	**
** glorified structure rather than a class since everything is public.	**
** We allow higher XWinNMR classes that deal with acquisition parameter	**
** files to be our friend so they have direct access to all parameters.	**
**									**
*************************************************************************/

#ifndef XWinAcqPar_ 			// Is file already included?
#  define XWinAcqPar_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the implementation
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <string>			// Include libstdc++ std::strings
#include <Basics/ParamSet.h>		// Include GAMMA parameter sets
#include <GamIO/XWinPSet.h>             // Include Bruker parameter parsing

class XWinAcqPar: public XWinPSet

// ----------------------------------------------------------------------------
// ------------------------------- STRUCTURE ----------------------------------
// ----------------------------------------------------------------------------
 
  {
  friend class XWinAcqus;		// Allow class XWinAcqus full access
  friend class XWinAcqu2s;		// Allow class XWinAcqu2s full access
  double       _Bo;                     // Field strength in Teslas
  std::string  parfile;			// Parameter file name (base)
  std::string  _TITLE;                  // Parameter file title
  std::string  _JCAMPDX;                // JCAMP version (5.0)
  std::string  _DATATYPE;               // Type of file (Parameter Values)
  std::string  _ORIGIN;                 // Where the file came from
  std::string  _OWNER;                  // Who made the file
  std::string  _DAY;			// Long format day and date
  std::string  _NAME;                   // The full file name
  int	        _AQSEQ;			// Acquisition sequence
  int	        _AQ_mod;			// Acquisition mode
  std::string  _AUNM;                   // To start AU higher level
  double       _BF[8];			// Base frequencies
  int          _BYTORDA;                // Byte ordering
  int          _CFDGTYP;                // 
  int          _CFRGTYP;                // 
  double       _CNST[32];		// Constants
  std::string  _CPDPRG;                 // 
  std::string  _CPDPRGN[8];		//
  std::string  _CPDPRGB;                // 
  std::string  _CPDPRGT;                // 
  double       _D_[32];			// Delay times (sec)
//  double       _D[32];			// Delay times (sec)
  //long         _DATE;                   // When the file was made
	time_t	     _DATE;
  double       _DBL[8];			// 
  double       _DBP[8];			// 
  double       _DBP07;
  std::string  _DBPNAM[8];
  double       _DBPOAL[8];
  double       _DBPOFFS[8];
  double       _DE;
  std::string  _DECBNUC;
  double       _DECIM;
  std::string  _DECNUC;
  int          _DECSTAT;
  int          _DIGMOD;
  int          _DIGTYP;
  int          _DL[8];			// 
  int          _DP[8];			// 
  int          _DP07;			// 
  std::string  _DPNAME[8];
  double       _DPOAL[8];
  double       _DPOFFS[8];
  int          _DQDMODE;
  int          _DR;	
  int          _DSc;		// Number of dummy scans
  std::string       _DSLIST;	
  int          _DSPFIRM;
  int          _DSPFVS;
  int          _DTYPA;
  std::string       _EXP;		// Experiment performed
  std::string       _F1LIST;
  std::string       _F2LIST;
  std::string       _F3LIST;
  int          _FCUCHAN[10];
  int          _FL1;
  int          _FL2;
  int          _FL3;
  int          _FL4;
  int          _FOV;
  std::string  _FQNLIST[8];
  int          _FS[8];
  int          _FTLPGN;
  int          _FW;
  int          _GP031;
  std::string  _GPNAME[32];	// Gradient parameters
  int          _GPX[32];
  int          _GPY[32];
  int          _GPZ[32];
  std::string  _GRDPROG;
  int          _HGAIN[4];
  int          _HL1;
  int          _HL2;
  int          _HL3;
  int          _HL4;
  int          _HOLDER;
  int          _HPMOD[8];
  int          _HPPRGN;
  double       _IN[32];		// Delay increments
  double       _INP[32];	// Pulse increments
  std::string  _INSTRUM;
  int          __L[32];		// 
  int          _LFILTER;
  int          _LGAIN;
  int          _LOCKPOW;
  std::string  _LOCNUC;
  int          _LOCPHAS;
  std::string  _LOCSHFT;
  double       _LTIME;
  int          _MASR;
  std::string       _MASRLST;
  int          _NBL;
  int          _NC;
  int          _NS;		// Number of scans
  std::string  _NUC[8];		// Channel nuclei
  int          _NUCLEI;
  std::string  _NUCLEUS;	// Base nucleus
  double       _O[8];		// Channel frequency offsets
  double       _OBSCHAN[10];	// 
  int          ___OVERFLW;
  double       ___P[32];
  int	        _PAPS;		//
  int          _PARMODE;	//
  double       _PCPD[10];
  double       _PHCOR[32];
  int          _PHP;
  int          _PH_ref;
  double       _PL[32];
  int          _POWMOD;		//
  int          __PR;		//
  double       _PRECHAN[16];
  int          _PRGAIN;
  std::string  _PROBHD;
  std::string  _PULPROG;	// Pulse program
  int          _PW;		//
  int          _QNP;		//
  int          _RD;		//
  int          _RECPH;		//
  int          _RG;		//
  int          _RO;		//
  int          _ROUTWD1[24];	//
  int          _ROUTWD2[24];	//
  int          _RSEL[10];	//
  int          __S[8];		//
  int          _SEOUT;		//
  double       _SFO[8];		// Channel irradiation frequencies
  std::string  _SOLVENT;	// Solvent
// sosi - MSVC++ compiler didn't line _SP so switched to _XSP!!!
  double       _XSP[16];		//
  int          _SP07;		//
  std::string  _SPNAM[16];	//
  double       _SPOAL[16];	//
  double       _SPOFFS[16];	//
  double       _SW;		// Spectral width (PPM)
  double       _SWIBOX[16];	//
  double       _SW_h;		// Spectral width (Hz)
  int          _TD;		// Total number of points (real + imag)
  double       _TE;		// Temperature (K)
  int          _TL[8];
  int          _TP[8];
  int          _TP07;
  std::string  _TPNAME[8];
  double       _TPOAL[8];
  double       _TPOFFS[8];
  int          _TUNHIN;
  int          _TUNHOUT;
  int          _TUNXOUT;
  std::string  _USERA[5];
  int          _V9;
  std::string  _VALIST;
  std::string  _VCLIST;
  int          _VD;
  std::string  _VDLIST;
  std::string  _VPLIST;
  std::string  _VTLIST;
  int          _WBST;
  int          _WBSW;
  int          _XGAIN[4];
  int          _XL;
  int          _YL;
  int          _YMAX_a;
  int          _YMIN_a;
  int          _ZL1;
  int          _ZL2;
  int          _ZL3;
  int          _ZL4;

// ----------------------------------------------------------------------------
// ----------------------------- PRIVATE FUNCTIONS ----------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                      XWinNMR AcqPar File Error Handling
// ____________________________________________________________________________

/* These functions take care of any errors encountered when reading, writing,   
   and setting parameters in Bruker acquisition parameter files.
 
	Input		AcqPar  : XWinNMR Acqusition parameters (this)
			eidx    : Error index
                        pname   : Additional error message
                        noret   : Flag for linefeed (0=linefeed)
        Output          void    : An error message is output                 */

         void XWinAcqParerror(int    eidx, int noret=0) const;
         void XWinAcqParerror(int    eidx, const std::string& pname, int noret=0) const;
volatile void XWinAcqParfatality(int eidx) const;
volatile void XWinAcqParfatality(int eidx, const std::string& pname) const;

// ____________________________________________________________________________
// ii                 XWinNMR AcqPar File Default Parmameters
// ____________________________________________________________________________

void SetDefaults(const std::string& fname);
void SetDefaults1(const std::string& fname);
void SetDefaults2(const std::string& fname);
void Copy(const XWinAcqPar& XWP);
 
// ____________________________________________________________________________
// iii          XWinNMR AcqPar Parmameter Setting Functions
// ____________________________________________________________________________
/* These are functions that can mess up things so they are private.
 
  Function                                Explanation
  ---------- ------------------------------------------------------------------
  SetField   Sets the stored field strength Bo from the isotope _NUCLEUS. This
             does NOT influence any other acqusition parameters.  In fact, Bo
             is itself not an acqusition parameter directly and needed only for
             resetting channel nuclei { NUC, BF, SF, O }
  FieldReset This allows for resetting of the spectrometer field strength Bo.
             It will recurse through all channels and reset their base
             frequencies as well as their spectrometer frequencies, keeping
             the offset values the same.
*/ 
 
void CheckNuclei();
void SetField();
void FieldReset(double BoT);
void SetIN(int i, double inval);
bool SetDelay(int idx, double tsec, int warn=2);
bool SetNucleus(int channel, const std::string& I, double off, int w=2);
void SetO(int i, double offset);
bool SetPulse(int idx, double tusec, int warn=2);
void SetSW(double sw, int inHz=0);

public:
// ____________________________________________________________________________ 
// A             XWinAcqPar Parameter File Constructors, Destructor
// ____________________________________________________________________________
 
/* These are the constructors of the class handling Bruker XWinNMR acquisition
   parameter files.  This doesn't do anything in particular, it is the read
   and write functions that perform the work.  Thus, we only have a default
   constructor specified.  The reading and writing of the associated ASCII
   parameter file is done in one step, so we don't need anything complex.    */

        XWinAcqPar();
        XWinAcqPar(const std::string& fname, int type=1);
        XWinAcqPar(const XWinAcqPar& XWA2);
virtual ~XWinAcqPar();
void    operator= (const XWinAcqPar& XWA2);


// ____________________________________________________________________________
// B                 XWinAcqPar Parameter Access Functions
// ____________________________________________________________________________
 
 /* These functions allow direct access to important acquisition parameters.
    The two primary values are the point size (TD) and the data byte order
    (BYTORDA).                                                               */

std::string acqname()  const;	// File name (short)
double field()    const;	// Field strength (T)
double AQ()       const;	// Acquisition length
int    AQ_mod()   const;	// Acquisition mode
double BF1()      const;	// Base Spectrometer freq.
double BF2()      const;	// Base Spectrometer freq.
int    BYTORDA()  const;	// Binary byte order
int    DSc()       const;	// Number of dummy scans
double DW()       const;	// The dwell time (sec)
std::string EXP()      const;	// Experiment name
double FIDRES()   const;	// The FID resolution (Hz)
// sosi - IN must be some defined system value as this is trouble
//        in GCC 3.2. Changed IN to XW_IN
double XW_IN(int i)  const;	// Dwell time
std::string NAME()     const;	// Full file name
int    NS()       const;	// Number of scans
std::string NUC(int i) const;	// Nucleus for a channel
std::string NUCLEUS()  const;	// Base nucleus
double O1()       const;	// Offset freq.
double O2()       const;	// Offset freq.
int    PARMODE()  const;	// Acquisition dimension
std::string PULPROG()  const;	// Pulse program
double SFO1()     const;	// Spectrometer freq.
double SFO2()     const;	// Spectrometer freq.
double SFO3()     const;	// Spectrometer freq.
std::string SOLVENT()  const;	// Solvent
double SW()       const;	// Spectral width (PPM)
double SW_h()     const;	// Spectral width (Hz)
double TE()       const;	// Sample temperature
int    TD()       const;	// Total points

void field(double bo);
void AQ_mod(int aqmo);
void BF1(double bf);
void BF2(double bf);
void BYTORDA(int bo);
int  D(int idx, double tsec, int warn=2);
void DECBNUC(const std::string& I);
void DECNUC(const std::string& I);
void DSc(int ds);
void EXP(const std::string& exp);
// sosi GCC 3.2 has trouble with IN, must be system defined.
// switched name to XW_IN
void XW_IN(int i, double in);
void O1(double of);
void O2(double of);
void O(int i, double of);
void NAME(const std::string& nam);
void NS(int ns);
void NUC(int i, const std::string& N);
void NUCLEI(int channel, const std::string& I, double O, int warn=2);
void NUCLEUS(const std::string& I);
int  P(int idx, double tp, int warn=2);
void PARMODE(int pm);
void PULPROG(const std::string& P);
void SFO1(double sf);
void SFO2(double sf);
void SFO3(double sf);
void SFO(double sf, int i);
void SOLVENT(const std::string& S);
void SW(double sw);
void SW_h(double sw);
void TE(double te);
void TD(int npts);

// ____________________________________________________________________________
// C                       XWinAcqPar Input Functions
// ____________________________________________________________________________

/* These functions will read all parameters from a Bruker XWinNMR acquisition
   parameter file. The entire Bruker parameter file is simply read into a GAMMA
   parameter set (class XWinPSet) so that ALL parameters in the file are stored
   in this class. Were the Bruker files in GAMMA format then this would be done
   exclusively with class ParamSet, but since they are not we must handle the
   file parsing explicitly. This is relegated to the base class XWinPSet. In
   turn, specific parameters of consequnce to acquisitons (e.g. acquisition
   parameter files) are parsed within this class.

      Function                               Purpose
   ____________________________________________________________________________

      readPSet          Read in parameter set (class XWinPSet).  This is done
                        special since Bruker format is NOT in GAMMA parameter
                        format so it must be parsed appropriately.
     parsePSet		Parses a Bruker parameter set (XWinPSet) for parameters 
			of consequence to an acquisition.
       read             Read in parameter set (class XWinPset) and parse out
			parameters of consequence to an acquisition.

                       INHERITED FROM BASE CLASS XWinPSet

        bool readPSet(const std::string& filein, int warn=1);
        bool readPSet(int warn=1);                               */
        bool readAPar(const std::string& filein, int warn=2);
	bool readAPar(int warn=2);
virtual bool parsePSet(int defs=1, int warn=2);

// ____________________________________________________________________________
// D                       XWinAcqPar Output Functions
// ____________________________________________________________________________

/* These function allow for output of NMR parameters directly into a Bruker
   XWinNMR ASCII parameter file {acqu(s) & acqu2(s)}.  Note that the value
   of PARMODE will affect the output somewhat.                               */

int writeAPar(const std::string& name, int warn=2);
int writeAPar(int warn=2) const;

// ____________________________________________________________________________
// E                 XWinAcqPar Parameter Access Functions
// ____________________________________________________________________________

/* These functions allow access to the Bruker XWinNMR parameter set, not just
  those parameters which are directly handled within this class.

                       INHERITED FROM BASE CLASS XWinPSet

int getPar(const std::string& pn,int& val,   int id=0,int wrn=0) const;
int getPar(const std::string& pn,double& val,int id=0,int wrn=0) const;
int getPar(const std::string& pn,std::string& val,int id=0,int wrn=0) const;
ParameterSet getPSet() const;                                    */
 

// ____________________________________________________________________________
// F                      XWinAcqPar Auxiliary Functions
// ____________________________________________________________________________
 
 /* These are just a group of functions that return std::strings indicating what
    the Bruker parameter values mean.  I just add things as I learn them here
    so that GAMMA output can remind me what all of these parameters are...   */
 
std::string AQ_modS()   const;
std::string BYTORDAS()  const;
std::string DECSTATS()  const;
std::string DIGTYPS()   const;
std::string DSS()       const;
std::string HPPRGNS()   const;
std::string PAPSS()     const;
std::string PARMODES()  const;
std::string POWMODS()   const;
std::string PRGAINS()   const;
std::string SEOUTS()    const;
std::string SWS()       const;
};

#endif 								// XWinAcqPar
