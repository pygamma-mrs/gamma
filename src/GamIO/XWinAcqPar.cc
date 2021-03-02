/* XWinAcqPar.cc *************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      XWinAcqPar                                  Implementation	**
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

#ifndef   XWinAcqPar_CC_		// Is file already included?
#  define XWinAcqPar_CC_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <GamIO/XWinAcqPar.h>		// Include the interface
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <Basics/Gconstants.h>		// Include GAMMA constants
#include <Basics/StringCut.h>		// Include GAMMA string parsing
#include <Basics/Isotope.h>		// Include GAMMA spin isotopes
#include <GamIO/BinIOBase.h>		// Include binary IO functions
#include <string>			// Include libstdc++ strings
#include <ctime>			// Include time and date access
#ifndef _MSC_VER                        // If not using MSVC++ then we
 #include <sys/time.h>			// Include time and date access
 #include <unistd.h>			// Include POSIX getcwd function
#else                                   // and if using MSVC++ then I guess
 #include <direct.h>			// Include time and date access
#endif

using std::string;			// Using libstdc++ strings
using std::ofstream;			// Using libstdc++ output file streams

// ----------------------------------------------------------------------------
// ----------------------------- PRIVATE FUNCTIONS ----------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                      XWinNMR AcqPar File Error Handling
// ____________________________________________________________________________
 
/* These functions take care of any errors encountered when reading, writing,
   and setting parameters in Bruker acquisition parameter files.

        Input		eidx    : Error index
        		pname   : Additional error message
        		noret	: Flag for linefeed (0=linefeed)
        Output		void    : An error message is output                 */
 
void XWinAcqPar::XWinAcqParerror(int eidx, int noret) const
  {
  string hdr("XWinNMR AcqPar Parameter File");
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

void XWinAcqPar::XWinAcqParerror(int eidx, const string& pname, int noret) const
  {                                                                             
  string hdr("XWinNMR Acquisition Parameter File");
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

volatile void XWinAcqPar::XWinAcqParfatality(int eidx) const
  {                                                                 
  XWinAcqParerror(eidx, 1);			// Normal non-fatal error
  if(eidx) XWinAcqParerror(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

volatile void XWinAcqPar::XWinAcqParfatality(int eidx, const string& pname) const
  {                                                                 
  XWinAcqParerror(eidx, pname, 1);		// Normal non-fatal error
  if(eidx) XWinAcqParerror(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }
 
// ____________________________________________________________________________
// ii                   XWinNMR AcqPar Default Parmameters
// ____________________________________________________________________________
 

void XWinAcqPar::SetDefaults(const string& fname)
  {
  parfile   = fname;
  _TITLE    = "Parameter File,  Version xwinnmr1.1  GAMMA Generated";
  _JCAMPDX  = "5.0";
  _DATATYPE = "Parameter Values";
  _ORIGIN   = "UXNMR, Bruker Analytische Messtechnik GmbH";
  _OWNER    = "GAMMA";
  //struct tm *ptr;               	// For setting current date
  time_t longtime;              	// Need a time structure
  longtime  = NULL; 

#ifdef _MSC_VER
  struct tm newtime;
	errno_t err = localtime_s(&newtime, &longtime);
	char day[31];
	err = asctime_s(day, 30, &newtime);
	_DAY = day;
  string cwd;
  cwd = string(_getcwd(NULL, 128));
#else
  const struct tm * ptr = localtime(&longtime);
  _DAY      = asctime(ptr); 
  string cwd;
  cwd = string(getcwd(NULL, 128));
#endif
  
  _NAME = cwd + "/" + fname; 		// Full file name
  _AQSEQ    = 0;			// Acquire sequence in 2D/3D
  _AQ_mod   = 1;			// Acquisiton mode = qsim
  _AUNM     = "";			// Acquisiton program
//_BF 					// Values are in loop below
  _BYTORDA = WeRBigEnd();       	// Byte ordering
  _CFDGTYP = 0;				// Unknown
  _CFRGTYP = 0;				// Unknown
//_CNST 				// Values are in loop below
  _CPDPRG  = "";			// Unknown
//_CPDPRG#				// Values are in loop below
  _CPDPRGB = "";			// Unknown
  _CPDPRGT = "";			// Unknown
//_D					// Values are in loop below
  time(&_DATE); 			// Compact date value
//_DBL					// Values are in loop below
//_DBP					// Values are in loop below
  _DBP07   = 0;
//_DBPNAME				// Values are in loop below
//_DBPOAL				// Values are in loop below
//_DBPOFFS				// Values are in loop below
  _DE      = 0.0;			// 
  _DECBNUC = "off";			// Same as NUC3 in AMX/ARX
  _DECIM   = 1;				// 
  _DECNUC  = "1H";			// Same as NUC2 in AMX/ARX
  _DECSTAT = 4;				// 
  _DIGMOD  = 1;				// 
  _DIGTYP  = 1;				// 
//_DL					// Values are in loop below
//_DP					// Values are in loop below
  _DP07    = 0;				// 
//_DPNAME				// Values are in loop below
//_DPOAL				// Values are in loop below
//_DPOFFS				// Values are in loop below
  _DQDMODE = 0;				// 
  _DR      = 16;			// 
  _DSc      = 0;				// Number of dummy scans
  _DSLIST  = "SSSSSSSSSSSSSSS";		// 
  _DSPFIRM = 0;				// 
  _DSPFVS  = 0;				// 
  _DTYPA   = 0;				// 
  _EXP     = "gammasim";		// Experiment performed	
  _F1LIST  = "111111111111111";		// 
  _F2LIST  = "222222222222222";		// 
  _F3LIST  = "333333333333333";		// 
//_FCUCHAN				// Values are in loop below
  _FL1     = 1;				// 
  _FL2     = 83;			// 
  _FL3     = 83;			// 
  _FL4     = 83;			// 
  _FOV     = 20;			// 
//_FQNLISTN				// Values are in loop below
//_FS					// Values are in loop below
  _FTLPGN  = 0;				// 
  _FW      = 8000;			// 
  _GP031   = 0;				// 
//_GPNAM				// Values are in loop below
//_GPX					// Values are in loop below
//_GPY					// Values are in loop below
//_GPZ					// Values are in loop below
  _GRDPROG = "";			// 
//_HGAIN				// Values are in loop below
  _HL1     = 90;				// 
  _HL2     = 35;			// 
  _HL3     = 8;				// 
  _HL4     = 26;			// 
  _HOLDER  = 0;				// 
//_HPMOD				// Values are in loop below
  _HPPRGN  = 0;				// 
//_IN					// Values are in loop below
//_INP					// Values are in loop below
  _INSTRUM= "GAMMA";			// 
//_L					// Values are in loop below
  _LFILTER = 10;			// 
  _LGAIN   = -10;			// 
  _LOCKPOW = -20;			// 
  _LOCNUC  = "2H";			// 
  _LOCPHAS = 0;				// 
  _LOCSHFT = "no";			// 
  _LTIME   = 0.1;			// 
  _MASR    = 0;				// 
  _MASRLST = "masrlst";			// 
  _NBL     = 1;				// 
  _NC      = 1;				// 
  _NS      = 1;				// Number of scans
//_NUC					// Values are in loop below
  _NUCLEI  = 0;				// 
  _NUCLEUS = "1H";			// Same as NUC1 in AMX/ARX
//_O					// Values are in loop below
//_OBSCHAN				// Values are in loop below
//_OVERFLW				// Values are in loop below
//_P					// Values are in loop below
  _PAPS    = 2;
  _PARMODE = 0;
//_PCPD					// Values are in loop below
//_PHCOR				// Values are in loop below
  _PHP     = 1;
  _PH_ref  = 0;
//_PL					// Values are in loop below
  _POWMOD  = 0;
  __PR     = 1;
//_PRECHAN				// Values are in loop below
  _PROBHD  = "";			// 
  _PULPROG = "gamma";			// 
  _PW      = 0;
  _QNP     = 1;
  _RD      = 0;
  _RECPH   = 0;
  _RG      = 0;
  _RO      = 4;
//_ROUTWD1				// Values are in loop below
//_ROUTWD2				// Values are in loop below
//_RSEL					// Values are in loop below
//_S					// Values are in loop below
  _SEOUT   = 0;				// 
//_SF					// Values are in loop below
  _SOLVENT = "cdcl3";			// Solvent
//_SP					// Values are in loop below
  _SP07    = 0;				// 
//_SPNAME				// Values are in loop below
//_SPOAL				// Values are in loop below
//_SPOFFS				// Values are in loop below
  _SW      = 10;			// No spectral width (PPM)
//_SWIBOX				// Values are in loop below
  _SW_h    = 5000;			// No spectral width (Hz)
  _TD      = 1024;			// Set for 1K acqustion points
  _TE      = 293;			// Temperature in K
//_TL					// Values are in loop below
//_TP					// Values are in loop below
  _TP07    = 0;				// 
//_TPNAME				// Values are in loop below
//_TPOAL				// Values are in loop below
//_TPOFFS				// Values are in loop below
  _TUNHIN  = 0;				// 
  _TUNHOUT = 0;				// 
  _TUNXOUT = 0;				// 
//_USERA				// Values are in loop below
  _V9	   = 5;				// 
  _VALIST  = "valist";			// 
  _VCLIST  = "CCCCCCCCCCCCCCC";		// 
  _VD	   = 0;				// 
  _VDLIST  = "DDDDDDDDDDDDDDD";		// 
  _VPLIST  = "PPPPPPPPPPPPPPP";		// 
  _VTLIST  = "TTTTTTTTTTTTTTT";		// 
  _WBST    = 1024;			// 
  _WBSW    = 4;				// 
//_XGAIN				// Values are in loop below
  _XL      = 0;				// 
  _YL      = 0;				// 
  _YMAX_a  = 169433;			// 
  _YMIN_a  = -125068;			// 
  _ZL1     = 120;			// 
  _ZL2     = 120;			// 
  _ZL3     = 120;			// 
  _ZL4     = 120;			// 
  int i=0;
  for(i=0; i<4; i++)			// Parameters over 4 values
    {
    _HGAIN[i]   = 0;
    _XGAIN[i]   = 0;
    }
  for(i=0; i<5; i++)			// Parameters over 5 values
    {
    _USERA[i]   = "user";
    }
  for(i=0; i<8; i++)			// Parameters over 8 channels
    {
    _BF[i]      = 500.00;		// Set base frequencies to 500
    _CPDPRGN[i] = "mlev";
    _DBL[i]     = 120;			// Set decoupler stuff
    _DBP[i]     = 150;			// 
    _DBPNAM[i]  = "";			// 
    _DBPOAL[i]  = 0.5;			// 
    _DBPOFFS[i] = 0;			// 
    _DL[i]      = 120;			// 
    _DP[i]      = 150;			// 
    _DPNAME[i]  = "";			// 
    _DPOAL[i]   = 0.5;			// 
    _DPOFFS[i]  = 0;			// 
    _FQNLIST[i] = "freqlist";		// 
    _FS[i]      = 83;			// 
    _HPMOD[i]   = 0;			// 
    _NUC[i]     = "off";		// Set isotope types to off
    _O[i]       = 0;			// Set offset frequencies to 0
    __S[i]      = 83;			// 
    _SFO[i]     = _BF[i] - _O[i];	// Set irradiation frequencies
    _TL[i]      = 120;			// 
    _TP[i]      = 150;			// 
    _TPNAME[i]  = "";
    _TPOAL[i]   = 0.5;			// 
    _TPOFFS[i]  = 0;			// 
    }
  for(i=0; i<10; i++)			// Parameters over 10 values
    {
    _FCUCHAN[i] = 0;
    _OBSCHAN[i] = 0;
    _PCPD[i]    = 100;
    _RSEL[i]    = 0;
    }
  for(i=0; i<16; i++)			// Parameters over 16 values
    {
    _PRECHAN[i] = 0;
    _SPNAM[i]   = "";
    }
  for(i=0; i<24; i++)			// Parameters over 24 values
    {
    _ROUTWD1[i] = 0;
    _ROUTWD2[i] = 0;
    }
  for(i=0; i<32; i++)			// Parameters in 32 sized arrays
    {
    _CNST[i]   = 1;			// Set all constants to 1
    _D_[i]      = 0;			// Set all delay times to 0
    _GPNAME[i] = "";			// Set all gradient parameters to null
    _GPX[i]    = 0;			// 
    _GPY[i]    = 0;			// 
    _GPZ[i]    = 0;			// 
    _IN[i]     = 0;			// Set all # delay increments to 0
    _INP[i]    = 0;			// Set all pulse increment times to 0
    __L[i]     = 1;
    ___P[i]    = 1;
    _PHCOR[i]  = 0;
    _PL[i]     = 120;
    }
  _NUC[0] = _NUCLEUS;			// NUCLEUS is NUC1 in AMX/ARX
  _NUC[1] = _DECNUC;			// DECNUC  is NUC2 in AMX/ARX
  _NUC[2] = _DECBNUC;			// DECBNUC is NUC3 in AMX/ARX
  _IN[0] = 1/_SW_h;			// Set D0 increment
  SetField();				// Set the Bo field strength
  }
 

void XWinAcqPar::SetDefaults1(const string& fname)
  { SetDefaults(fname); }		// Use generic defaults

void XWinAcqPar::SetDefaults2(const string& fname)
  {
  SetDefaults(fname);			// Set generic defaults
  _DIGTYP  = 2;				// 
  _DP07    = 0;				// 
  _DR      = 2;				// May actually be ND?
  _FW      = 118;			// Filter width
  _HL2     = 90;			// 
  _HL3     = 90;			// 
  _HL4     = 90;			// 
  _INSTRUM= "";				// No instrument specified
  _LOCNUC  = "";			// No lock nucleus
  _PARMODE = 1;				// Aquisition dimension (1=2D)
  _TD      = 256;			// No t1 incrementation points
  }

void XWinAcqPar::Copy(const XWinAcqPar& XWAP)
  {
  parfile   = XWAP.parfile;		// Copy file name
  _Bo       = XWAP._Bo;			// Copy field srength
  _TITLE    = XWAP._TITLE;		// Copy file title
  _JCAMPDX  = XWAP._JCAMPDX;		// Copy JCAMP version
  _DATATYPE = XWAP._DATATYPE;		// Copy the data type
  _ORIGIN   = XWAP._ORIGIN;
  _OWNER    = XWAP._OWNER;
  _DAY      = XWAP._DAY;
  _NAME     = XWAP._NAME;
  _AQSEQ    = XWAP._AQSEQ;		// Copy
  _AQ_mod   = XWAP._AQ_mod;		// Copy acquisition mode
  _AUNM     = XWAP._AUNM;		// 
//_BF 					// Values are in loop below
  _BYTORDA  = XWAP._BYTORDA;		// Copy binary byte order
  _CFDGTYP  = XWAP._CFDGTYP;
  _CFRGTYP  = XWAP._CFRGTYP;
//_CNST 				// Values are in loop below
  _CPDPRG   = XWAP._CPDPRG;
//_CPDPRG#				// Values are in loop below
  _CPDPRGB  = XWAP._CPDPRGB;
  _CPDPRGT  = XWAP._CPDPRGT;
//_D					// Values are in loop below
//_DBL					// Values are in loop below
//_DBP					// Values are in loop below
  _DBP07    = XWAP._DBP07;
//_DBPNAME				// Values are in loop below
//_DBPOAL				// Values are in loop below
//_DBPOFFS				// Values are in loop below
  _DE       = XWAP._DE;
  _DECBNUC  = XWAP._DECBNUC;		// Copy 3rd nuc. AMX/ARX (NUC3)
  _DECIM    = XWAP._DECIM;	
  _DECNUC   = XWAP._DECNUC;		// Copy 2nd nuc. AMX/ARX (NUC2)
  _DECSTAT  = XWAP._DECSTAT;
  _DIGMOD   = XWAP._DIGMOD;
  _DIGTYP   = XWAP._DIGTYP;
//_DL					// Values are in loop below
//_DP					// Values are in loop below
  _DP07     = XWAP._DP07;
//_DPNAME				// Values are in loop below
//_DPOAL				// Values are in loop below
//_DPOFFS				// Values are in loop below
  _DQDMODE  = XWAP._DQDMODE;
  _DR       = XWAP._DR;
  _DSc       = XWAP._DSc;			// Copy dummy scans
  _DSLIST   = XWAP._DSLIST;
  _DSPFIRM  = XWAP._DSPFIRM;
  _DSPFVS   = XWAP._DSPFVS;
  _DTYPA    = XWAP._DTYPA;
  _EXP      = XWAP._EXP;
  _F1LIST   = XWAP._F1LIST;
  _F2LIST   = XWAP._F2LIST;
  _F3LIST   = XWAP._F3LIST;
//_FCUCHAN				// Values are in loop below
  _FL1      = XWAP._FL1;
  _FL2      = XWAP._FL2;
  _FL3      = XWAP._FL3;
  _FL4      = XWAP._FL4;
  _FOV      = XWAP._FOV;
//_FQNLISTN				// Values are in loop below
//_FS					// Values are in loop below
  _FTLPGN   = XWAP._FTLPGN;
  _FW       = XWAP._FW;
  _GP031    = XWAP._GP031;
//_GPNAM				// Values are in loop below
//_GPX					// Values are in loop below
//_GPY					// Values are in loop below
//_GPZ					// Values are in loop below
  _GRDPROG  = XWAP._GRDPROG;
//_HGAIN				// Values are in loop below
  _HL1      = XWAP._HL1;
  _HL2      = XWAP._HL2;
  _HL3      = XWAP._HL3;
  _HL4      = XWAP._HL4;
  _HOLDER   = XWAP._HOLDER;
//_HPMOD				// Values are in loop below
  _HPPRGN   = XWAP._HPPRGN;
//_IN					// Values are in loop below
//_INP					// Values are in loop below
  _INSTRUM  = XWAP._INSTRUM;
//_L					// Values are in loop below
  _LFILTER  = XWAP._LFILTER;
  _LGAIN    = XWAP._LGAIN;
  _LOCKPOW  = XWAP._LOCKPOW;
  _LOCNUC   = XWAP._LOCNUC;
  _LOCPHAS  = XWAP._LOCPHAS;
  _LOCSHFT  = XWAP._LOCSHFT;
  _LTIME    = XWAP._LTIME;
  _MASR     = XWAP._MASR;
  _MASRLST  = XWAP._MASRLST;
  _NBL      = XWAP._NBL;
  _NC       = XWAP._NC;
  _NS       = XWAP._NS;
//_NUC					// Values are in loop below
  _NUCLEI   = XWAP._NUCLEI;
  _NUCLEUS  = XWAP._NUCLEUS;
//_O					// Values are in loop below
//_OBSCHAN				// Values are in loop below
//_OVERFLW				// Values are in loop below
//_P					// Values are in loop below
  _PAPS     = XWAP._PAPS;
  _PARMODE  = XWAP._PARMODE;
//_PCPD					// Values are in loop below
//_PHCOR				// Values are in loop below
  _PHP      = XWAP._PHP;
  _PH_ref   = XWAP._PH_ref;
//_PL					// Values are in loop below
  _POWMOD   = XWAP._POWMOD;
  __PR      = XWAP.__PR;
//_PRECHAN				// Values are in loop below
  _PROBHD   = XWAP._PROBHD;
  _PULPROG  = XWAP._PULPROG;
  _PW       = XWAP._PW;
  _QNP      = XWAP._QNP;
  _RD       = XWAP._RD;
  _RECPH    = XWAP._RECPH;
  _RG       = XWAP._RG;
  _RO       = XWAP._RO;
//_ROUTWD1				// Values are in loop below
//_ROUTWD2				// Values are in loop below
//_RSEL					// Values are in loop below
//_S					// Values are in loop below
  _SEOUT    = XWAP._SEOUT;
//_SF					// Values are in loop below
  _SOLVENT  = XWAP._SOLVENT;
//_SP					// Values are in loop below
  _SP07     = XWAP._SP07;
//_SPNAM				// Values are in loop below
//_SPOAL				// Values are in loop below
//_SPOFFS				// Values are in loop below
  _SW       = XWAP._SW;			// Copy spectral width (ppm)
//_SWIBOX				// Values are in loop below
  _SW_h     = XWAP._SW_h;		// Copy spectral width (Hz)
  _TD       = XWAP._TD;
  _TE       = XWAP._TE;
//_TL					// Values are in loop below
//_TP					// Values are in loop below
  _TP07     = XWAP._TP07;
//_TPNAME				// Values are in loop below
//_TPOAL				// Values are in loop below
//_TPOFFS				// Values are in loop below
  _TUNHIN   = XWAP._TUNHIN;
  _TUNHOUT  = XWAP._TUNHOUT;
  _TUNXOUT  = XWAP._TUNXOUT;
//_USERA				// Values are in loop below
  _V9	    = XWAP._V9;
  _VALIST   = XWAP._VALIST;
  _VCLIST   = XWAP._VCLIST;
  _VD	    = XWAP._VD;
  _VDLIST   = XWAP._VDLIST;
  _VPLIST   = XWAP._VPLIST;
  _VTLIST   = XWAP._VTLIST;
  _WBST     = XWAP._WBST;
  _WBSW     = XWAP._WBSW;
//_XGAIN				// Values are in loop below
  _XL       = XWAP._XL;
  _YL       = XWAP._YL;
  _YMAX_a   = XWAP._YMAX_a;
  _YMIN_a   = XWAP._YMIN_a;
  _ZL1      = XWAP._ZL1;
  _ZL2      = XWAP._ZL2;
  _ZL3      = XWAP._ZL3;
  _ZL4      = XWAP._ZL4;
  int i=0;
  for(i=0; i<4; i++)			// Parameters over 4 values
    {
    _HGAIN[i]   = XWAP._HGAIN[i];
    _XGAIN[i]   = XWAP._XGAIN[i];
    }
  for(i=0; i<5; i++)			// Parameters over 5 values
    {
    _USERA[i]   = XWAP._USERA[i];
    }
  for(i=0; i<8; i++)			// Parameters over 8 channels
    {
    _BF[i]      = XWAP._BF[i];		// Set base frequencies to 500
    _CPDPRGN[i] = XWAP._CPDPRGN[i];
    _DBL[i]     = XWAP._DBL[i];		// Set decoupler stuff
    _DBP[i]     = XWAP._DBP[i];		// 
    _DBPNAM[i]  = XWAP._DBPNAM[i];	// 
    _DBPOAL[i]  = XWAP._DBPOAL[i];	// 
    _DBPOFFS[i] = XWAP._DBPOFFS[i];	// 
    _DL[i]      = XWAP._DL[i];		// 
    _DP[i]      = XWAP._DP[i];		// 
    _DPNAME[i]  = XWAP._DPNAME[i];	// 
    _DPOAL[i]   = XWAP._DPOAL[i];	// 
    _DPOFFS[i]  = XWAP._DPOFFS[i];	// 
    _FQNLIST[i] = XWAP._FQNLIST[i];	// 
    _FS[i]      = XWAP._FS[i];		// 
    _HPMOD[i]   = XWAP._HPMOD[i];	// 
    _NUC[i]     = XWAP._NUC[i];		// Set isotope types to off
    _O[i]       = XWAP._O[i];		// Set offset frequencies to 0
    __S[i]      = XWAP.__S[i];		// 
    _SFO[i]     = XWAP._SFO[i];		// Set irradiation frequencies to 0
    _TL[i]      = XWAP._TL[i];		// 
    _TP[i]      = XWAP._TP[i];		// 
    _TPNAME[i]  = XWAP._TPNAME[i];
    _TPOAL[i]   = XWAP._TPOAL[i];	// 
    _TPOFFS[i]  = XWAP._TPOFFS[i];	// 
    }
  for(i=0; i<10; i++)			// Parameters over 10 values
    {
    _FCUCHAN[i] = XWAP._FCUCHAN[i];
    _OBSCHAN[i] = XWAP._OBSCHAN[i];
    _PCPD[i]    = XWAP._PCPD[i];
    _RSEL[i]    = XWAP._RSEL[i];
    }
  for(i=0; i<16; i++)			// Parameters over 16 values
    {
    _PRECHAN[i] = XWAP._PRECHAN[i];
    _SPNAM[i]   = XWAP._SPNAM[i];
    }
  for(i=0; i<24; i++)			// Parameters over 24 values
    {
    _ROUTWD1[i] = XWAP._ROUTWD1[i];
    _ROUTWD2[i] = XWAP._ROUTWD2[i];
    }
  for(i=0; i<32; i++)			// Parameters in 32 sized arrays
    {
    _CNST[i]   = XWAP._CNST[i];		// Set all constants to 1
    _D_[i]     = XWAP._D_[i];		// Set all delay times to 0
    _GPNAME[i] = XWAP._GPNAME[i];	// Set all gradient parameters to null
    _GPX[i]    = XWAP._GPX[i];		// 
    _GPY[i]    = XWAP._GPY[i];		// 
    _GPZ[i]    = XWAP._GPZ[i];		// 
    _IN[i]     = XWAP._IN[i];		// Set all delay increment times to 0
    _INP[i]    = XWAP._INP[i];		// Set all pulse increment times to 0
    __L[i]     = XWAP.__L[i];
    ___P[i]    = XWAP.___P[i];
    _PHCOR[i]  = XWAP._PHCOR[i];
    _PL[i]     = XWAP._PL[i];
    }
  }

// ____________________________________________________________________________
// iii          XWinNMR AcqPar File Parmameter Setting Functions
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

void XWinAcqPar::CheckNuclei()
  {
  if(_NUCLEUS != _NUC[0])			// Insure NUCLEUS = NUC1 
    {
    if(Isotope::known(_NUCLEUS))	
      _NUC[0] = _NUCLEUS;
    else _NUCLEUS = _NUC[0];
    }
  if(_DECNUC != _NUC[1])			// Insure DECNUC = NUC2 
    {
    if(Isotope::known(_DECNUC))	
      _NUC[1] = _DECNUC;
    else _DECNUC = _NUC[1];
    }
  if(_DECBNUC != _NUC[2])			// Insure DECBNUC = NUC3 
    {
    if(Isotope::known(_DECBNUC))	
      _NUC[2] = _DECBNUC;
    else _DECBNUC = _NUC[2];
    }
  }

void XWinAcqPar::SetField()
  {
  Isotope X(_NUCLEUS);
  _Bo = _BF[0]*1.e6*HZ2RAD/X.gamma(); 
  }

void XWinAcqPar::FieldReset(double BoT)		// Reset Spectrometer Field (T)
  {  
  _Bo = BoT;					// Set the field strength
  string I;					// Possible nucleus symbol
  CheckNuclei();				// Insure nuclei consistent
  for(int i=0; i<8; i++)                        // Loop possible channels
    {
    I = _NUC[i];				//  Get isotope type
    if(Isotope::known(I))			//  If this is a valid type
      {
      Isotope X(I);
      _BF[i] = X.gamma()*_Bo/(1.e6*HZ2RAD);
      _SFO[i] = _BF[i] - _O[i];
      }
    }
  }

void XWinAcqPar::SetSW(double sw, int inHz)
  {
  if(inHz) { _SW_h = sw; _SW   = sw/_SFO[0]; }
  else     { _SW   = sw; _SW_h = sw*_SFO[0]; }
  if(_PARMODE == 1)			// For a 2D experiment
    { _IN[0] = 1.0/_SW_h*_DR; } 	// we must set IN0
//  if(_PARMODE == 2)			// For a 3D experiment
//    {
  }

void XWinAcqPar::SetIN(int i, double inval)
  {
  _IN[i] = inval;
/*
  if(_PARMODE == 1 && i==0)		// For a 2D experiment we
    {					// recalculate SW for the
    SW = 1.0/(_SFO[i]*_DR*inval;	// F1 dimension
    }
  if(_PARMODE == 2 && i==0)		// For a 3D experiment we
    {					// recalculate SW for the
    SW = 1.0/(_SFO[i]*_DR*inval;	// F1 dimension
*/
  }

void XWinAcqPar::SetO(int i, double offset)
  {
  _O[i] = offset;
  _SFO[i] = _BF[i] - _O[i];
  }


bool XWinAcqPar::SetNucleus(int channel, const string& I, double offset, int warn)

	// Input	channel	: Spectrometer channel [1, 8]
	//		I       : Isotope label
	// Output	void	: The XWin32 NUCLEI setup is performed on
	//			  the channel indicated.  The parameters
	//			        { NUC#, BF#, SFO#, O# } 
	//			  are set. This involves not only setting a
	//			  label but setting a basic frequency BF#
	//			  where # is the channel
	// Note			: The following is always maintained
	//				   SR = SF - BF1
	// Note			: Some parameters are identical in this
	//			  treatment as they are for AMX/ARX systems
	//			      NUCLEUS == NUC1 == _NUC[0]
	//			      DECNUC  == NUC2 == _NUC[1]
	//			      DECBNUC == NUC2 == _NUC[2]
  {
  if(!_Bo)
    {
    if(warn)
      {
      XWinAcqParerror(31, 1);				// No Bo Set
      XWinAcqParerror(32, 1);				// Set Bo with B0
      if(warn > 1) XWinAcqParfatality(23, I);		// Can't set I channel
      else         XWinAcqParerror(23, I, 1);
      }
    return false;
    }
  Isotope X;
  if(!X.exists(I))
    {
    if(warn)
      {
      XWinAcqParerror(22, I, 1);			// Can't set I channel
      if(warn > 1) XWinAcqParfatality(24, I);		// Unknown nucleus
      else         XWinAcqParerror(24, I, 1);
      }
    return false;
    }
  if(channel<1 || channel>8)
    {
    if(warn)
      {
      XWinAcqParerror(30, 1);				// Channesl [1,8] only
      string pname = string("BF") + Gdec(channel);
      if(warn > 1) XWinAcqParfatality(22, pname);	// Problems setting BF
      else         XWinAcqParerror(22, pname, 1);
      }
    return false;
    }
  Isotope II(I);					// For an isotope
  _NUC[channel-1] = I;					// Set nucleus
  _BF[channel-1]  = fabs(II.gamma())*_Bo*1.e-6/HZ2RAD;	// Set base frequency
  _O[channel-1]   = offset;				// Set offset frequency
  _SFO[channel-1] = _BF[channel-1] - offset;		// Set irrad. frequency
  if(channel == 1)      _NUCLEUS = I; 
  else if(channel == 2) _DECNUC = I;
  else if(channel == 3) _DECBNUC = I;
  return true;
  }

bool XWinAcqPar::SetDelay(int idx, double tsec, int warn)
  {
  if(idx<0 || idx>31)
    {
    if(warn)
      {
      XWinAcqParerror(26, 1);				// Delay [0,31] only
      string pname = string("D") + Gdec(idx);
      if(warn > 1) XWinAcqParfatality(22, pname);	// Problems setting D
      else         XWinAcqParerror(22, pname, 1);
      }
    return false;
    }
  if(tsec<0)
    {
    if(warn)
      {
      XWinAcqParerror(27, 1);				// Delay with - time
      string pname = string("D") + Gdec(idx);
      if(warn > 1) XWinAcqParfatality(22, pname);	// Problems setting D
      else         XWinAcqParerror(22, pname, 1);
      }
    return 0;
    }
  _D_[idx] = tsec;
  return true;
  }

bool XWinAcqPar::SetPulse(int idx, double tp, int warn)
  {
  if(idx<0 || idx>31)
    {
    if(warn)
      {
      XWinAcqParerror(28, 1);				// Pulse [0,31] only
      string pname = string("P") + Gdec(idx);
      if(warn > 1) XWinAcqParfatality(22, pname);	// Problems setting P
      else         XWinAcqParerror(22, pname, 1);
      }
    return false;
    }
  if(tp<0)
    {
    if(warn)
      {
      XWinAcqParerror(29, 1);				// Delay with - time
      string pname = string("P") + Gdec(idx);
      if(warn > 1) XWinAcqParfatality(22, pname);	// Problems setting P
      else         XWinAcqParerror(22, pname, 1);
      }
    return false;
    }
  ___P[idx] = tp;					// Set pulse length
  return true;
  }
 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A            XWinAcqPar Parameter File Constructors, Destructor
// ____________________________________________________________________________

/* These are the constructors of the class handling Bruker XWinNMR acquisition
   parameter files.  This doesn't do anything in particular, it is the read
   and write functions that perform the work.  Thus, we only have a default
   constructor specified.  The reading and writing of the associated ASCII
   parameter file is done in one step, so we don't need anything complex.    */

XWinAcqPar::XWinAcqPar() : XWinPSet() { SetDefaults1("acqus"); }

XWinAcqPar::XWinAcqPar(const string& fname, int type) : XWinPSet(fname)
  {
  switch(type)
    {
    case 1:
    default: SetDefaults1(fname); break;	// Set defaults for acqus
    case 2:  SetDefaults2(fname); break;	// Set defaults for acqu2s
    }
  }

XWinAcqPar::XWinAcqPar(const XWinAcqPar& XWAP) : XWinPSet(XWAP) { Copy(XWAP); }

XWinAcqPar::~XWinAcqPar() { }

void XWinAcqPar::operator= (const XWinAcqPar& XWAP)
  {
  XWinPSet::operator=(XWAP);			// Copy the parameter set
  Copy(XWAP);					// Copy the parameters
  }


// ____________________________________________________________________________
// B                 XWinAcqPar Parameter Access Functions
// ____________________________________________________________________________
 
/* These functions allow direct access to important acquisition parameters.
   The two primary values are the point size (TD) and the data byte order
   (BYTORDA).

  Function                                Explanation
  --------   ------------------------------------------------------------------
     AQ      For qf   = single channel, the time between pts is the dwell time.
             For qseq = quadra. sequential, time between pts is the dwell time.  
             For qsim = quadra. simultaneous, point pairs occur at the same
                        time but separated by 2*dwell time parameter DW
             For qdig = digital quadrature, uses DIGMOD.
     ND(i)   Number of delays which are increment in the dimension i.  This is
             then used to calculate the sweep widths in the F1 & F2 dimensions
             as
                    SF(F1) = 1/(SF01*ND0*IN0) SF(F2)=1/(SF01*ND10*IN10)
     IN(i)   Number of increments to take during a variable delay in dimension
             i. This is used to calculate the sweep width according the the 
             above equation.
     DW      Dwell time in microseconds.  Changing this parameter should alter
             SW as
                               SW = 10^6/[2.0*(0.05+DW)*SFO1]
             so that always
                                   DW = 10^6/[2.0*SW*SFO1]                    */


string XWinAcqPar::acqname()  const { return parfile;                     }
double XWinAcqPar::field()    const { return _Bo;                         }
double XWinAcqPar::AQ()       const
  { 
  double dws = 1.0/(2.0*_SW*_SFO[0]);			// dwell time (sec)
  double aqt, pts=double(_TD);
  switch(_AQ_mod)
    {
    default:
    case 0:					// qf=   single channel  
    case 1: aqt = pts*dws;       break;		// qseq= quadrature sequential (2)
    case 2: aqt = pts*dws; 	break;		// qsim= quad. simultaneous (1)**
    case 3: aqt = pts*dws;       break;		// qdig= digital quadrature
    }
  return aqt;
  }
int    XWinAcqPar::AQ_mod()   const { return _AQ_mod;                     }
double XWinAcqPar::BF1()      const { return _BF[0];                      }
double XWinAcqPar::BF2()      const { return _BF[1];                      }
int    XWinAcqPar::BYTORDA()  const { return _BYTORDA;                    }
int    XWinAcqPar::DSc()       const { return _DSc;                         }
double XWinAcqPar::DW()       const { return 1.e6/(2*_SW*_SFO[0]);        }
string XWinAcqPar::EXP()      const { return _EXP;                        }
double XWinAcqPar::FIDRES()   const { return _SW*_SFO[0]/_TD;             }
double XWinAcqPar::XW_IN(int i)  const { return _IN[i];                      }
string XWinAcqPar::NAME()     const { return _NAME;                       }
int    XWinAcqPar::NS()       const { return _NS;                         }
string XWinAcqPar::NUC(int i) const { return _NUC[i];                     }
string XWinAcqPar::NUCLEUS()  const { return _NUCLEUS;                    }
double XWinAcqPar::O1()       const { return _O[0];                       }
double XWinAcqPar::O2()       const { return _O[1];                       }
int    XWinAcqPar::PARMODE()  const { return _PARMODE;                    }
string XWinAcqPar::PULPROG()  const { return _PULPROG;                    }
double XWinAcqPar::SFO1()     const { return _SFO[0];                     }
double XWinAcqPar::SFO2()     const { return _SFO[1];                     }
double XWinAcqPar::SFO3()     const { return _SFO[2];                     }
string XWinAcqPar::SOLVENT()  const { return _SOLVENT;                    }
double XWinAcqPar::SW()       const { return _SW;                         }
double XWinAcqPar::SW_h()     const { return _SW_h;			  }
double XWinAcqPar::TE()       const { return _TE;                         }
int    XWinAcqPar::TD()       const { return _TD;                         }

void XWinAcqPar::field(double bo)          { FieldReset(bo);              }
void XWinAcqPar::AQ_mod(int aqmo)          { _AQ_mod  = aqmo;             }
void XWinAcqPar::BF1(double bf)            { _BF[0]   = bf;               }
void XWinAcqPar::BF2(double bf)            { _BF[1]   = bf;               }
void XWinAcqPar::BYTORDA(int bo)           { bo?_BYTORDA=1:_BYTORDA=0;    }
int  XWinAcqPar::D(int i, double t, int w) { return SetDelay(i,t,w);      }
void XWinAcqPar::DECBNUC(const string& I)  { NUC(2, I);                   }
void XWinAcqPar::DECNUC(const string& I)   { NUC(1, I);                   }
void XWinAcqPar::DSc(int ds)                { _DSc      = ds;               }
void XWinAcqPar::EXP(const string& exp)    { _EXP     = exp;              }
void XWinAcqPar::XW_IN(int i, double in)      { _IN[i]   = in;               }
void XWinAcqPar::O1(double of)             { SetO(0, of);                 }
void XWinAcqPar::O2(double of)             { SetO(1, of);                 }
void XWinAcqPar::O(int i, double of)       { SetO(i, of);                 }

void XWinAcqPar::NAME(const string& nm)    
  { 
#ifdef _MSC_VER
	_NAME   = string(_getcwd(NULL,128));
#else
	_NAME   = string(getcwd(NULL,128));
#endif
	_NAME += "/" + nm;          
  }
void XWinAcqPar::NS(int ns)                { _NS      = ns;               }
void XWinAcqPar::NUC(int i,const string& N){ _NUC[i]  = N;                }
void XWinAcqPar::NUCLEI(int CH, const string& I, double OF, int warn)
                                           { SetNucleus(CH, I, OF, warn); }
void XWinAcqPar::NUCLEUS(const string& I)  { NUC(0, I);                   }
int  XWinAcqPar::P(int i, double t, int w) { return SetPulse(i,t,w);      }
void XWinAcqPar::PARMODE(int pm)           { _PARMODE = pm;               }
void XWinAcqPar::PULPROG(const string& P)  { _PULPROG = P;                }
void XWinAcqPar::SFO1(double sf)           { SFO(sf,0);                   }
void XWinAcqPar::SFO2(double sf)           { SFO(sf,1);                   }
void XWinAcqPar::SFO3(double sf)           { SFO(sf,2);                   }
void XWinAcqPar::SFO(double sf, int i)     { _SFO[i]  = sf;               }
void XWinAcqPar::SOLVENT(const string& S)  { _SOLVENT = S;                }
void XWinAcqPar::SW(double sw)             { SetSW(sw);                   }
void XWinAcqPar::SW_h(double sw)           { SetSW(sw, 1);                }
void XWinAcqPar::TD(int td)                { _TD      = td;               }
void XWinAcqPar::TE(double te)             { _TE      = te;               }

// ____________________________________________________________________________
// C                       XWinAcqPar Input Functions
// ____________________________________________________________________________

/* These functions read in the acquisition parameters from an XWinNMR parameter
   file, typically named acqus or acqu2s.  By design, the Bruker parameter
   file is initially read into a GAMMA parameter set so that ALL parameters in
   the file are stored.  Subsequently, the parameters in the parameter set are
   parsed to obtain values of consequence to GAMMA and these are explicitly
   maintained variables in the class.                                        */

bool XWinAcqPar::readAPar(const string& filein, int warn)
  {
  parfile = filein;				// Set internal filename
  bool TF = (readAPar(warn))?true:false;           // Read in all parameters
  if(!TF && warn)                               // Output errors if trouble
    {
    XWinAcqParerror(1, filein, 1);              // Problems with filein
    if(warn > 1) XWinAcqParerror(3);            // Can't construct from pset
    else         XWinAcqParerror(3,1);
    }
  return TF;
  }

bool XWinAcqPar::readAPar(int warn)
  {
  bool TF = XWinPSet::readPSet(parfile, warn);	// Use base class to read
  if(!TF) return TF;                            // Fail if we can read it
  TF = parsePSet(warn);				// Parse our essential params
  if(TF) SetField();				// Set the Bo value
  return TF;
  }    


bool XWinAcqPar::parsePSet(int defs, int warn)
  {
  if(defs==2) 
		SetDefaults2(parfile);		// First set the default
  else        
		SetDefaults1(parfile);		// parameters (same name)

//                    First Try And Access All Single Parameters

  int i=0;
  getPar("TITLE",    _TITLE);			// Try & get file title
  getPar("JCAMPDX",  _JCAMPDX);			// Try for JCAMP version
  getPar("DATATYPE", _DATATYPE);		// Try & get file type
  getPar("ORIGIN",   _ORIGIN); 			// Try & get file ORIGIN
  getPar("OWNER",    _OWNER);			// Try & get file OWNER
  getPar("DAY",      _DAY);			// Try for date/time
  getPar("NAME",     _NAME);			// Try & get file name
  getPar("_AQSEQ",   _AQSEQ);   		// Acquisition sequence
  getPar("_AQ_mod",  _AQ_mod);   		// Try for acquis. mode
  getPar("AUNM",     _AUNM); 			// Try and get AUNM
						// BF Vals Parsed in Loop
  getPar("BYTORDA",  _BYTORDA);			// Try and get byte order *
  getPar("CFDGTYP",  _CFDGTYP);			// Try and get this
  getPar("CFRGTYP",  _CFRGTYP);			// Try and get this
						// CNST Values Parsed in Loop
  getPar("CPDPRG",   _CPDPRG);			// Try and get this
						// CPDPRG# Values Parsed in Loop
  getPar("CPDPRGB",  _CPDPRGB);			// Try and get this
  getPar("CPDPRGT",  _CPDPRGT);			// Try and get this
						// D# Delays Parsed in Loop
  getPar("DATE",     i); 
	//_DATE = long(i);	// Try and get file date
	_DATE = time_t(i);
						// DBL# values Parsed in Loop
						// DBP# values Parsed in Loop
  getPar("DBP07",    _DBP07);			// Try and get this
						// DBPNAM# values Parsed in Loop
						// DBPOAL# values Parsed in Loop
						// DBPOFFS# vals Parsed in Loop
  getPar("DE",       _DE);			// Try and get this
  getPar("DECBNUC",  _DECBNUC);			// Try and get this
  getPar("DECIM",    _DECIM);			// Try and get this
  getPar("DECNUC",   _DECNUC);			// Try and get this
  getPar("DECSTAT",  _DECSTAT);			// Try and get this
  getPar("DIGMOD",   _DIGMOD);			// Try and get this
  getPar("DIGTYP",   _DIGTYP);			// Try and get this
						// DL# values Parsed in Loop
						// DP# values Parsed in Loop
  getPar("DP07",     _DP07);			// Try and get this
						// DPNAME# values Parsed in Loop
						// DPOAL# values Parsed in Loop
						// DPOFFS# values Parsed in Loop
  getPar("DQDMODE",  _DQDMODE);			// Try and get this
  getPar("DR",       _DR);			// Try and get this
  getPar("DS",       _DSc);			// Try and get # dummy scans
  getPar("DSLIST",   _DSLIST);			// Try and get this
  getPar("DSPFIRM",  _DSPFIRM);			// Try and get this
  getPar("DSPFVS",   _DSPFVS);			// Try and get this
  getPar("DTYPA",    _DTYPA);			// Try and get this
  getPar("EXP",      _EXP);			// Try and get experiment name
  getPar("F1LIST",   _F1LIST);			// Try and get this
  getPar("F2LIST",   _F2LIST);			// Try and get this
  getPar("F3LIST",   _F3LIST);			// Try and get this
						// FCUCHAN# vals Parsed in Loop
  getPar("FL1",      _FL1);			// Try and get this
  getPar("FL2",      _FL2);			// Try and get this
  getPar("FL3",      _FL3);			// Try and get this
  getPar("FL4",      _FL4);			// Try and get this
  getPar("FOV",      _FOV);			// Try and get this
						// FQNLIST# vals Parsed in Loop
						// FS# values Parsed in Loop
  getPar("FTLPGN",  _FTLPGN);			// Try and get this
  getPar("FW",      _FW);			// Try and get this
  getPar("GP031",   _GP031);			// Try and get this
						// GPNAME# values Parsed in Loop
						// GPX# values Parsed in Loop
						// GPY# values Parsed in Loop
						// GPZ# values Parsed in Loop
  getPar("GRDPROG", _GRDPROG);			// Try and get this
						// HGAIN# values Parsed in Loop
  getPar("HL1",     _HL1);			// Try and get this
  getPar("HL2",     _HL2);			// Try and get this
  getPar("FL3",     _HL3);			// Try and get this
  getPar("HL4",     _HL4);			// Try and get this
  getPar("HOLDER",  _HOLDER);			// Try and get this
						// HPMOD# values Parsed in Loop
  getPar("HPPRGN",  _HPPRGN);			// Try and get this
						// IN# values Parsed in Loop
						// INP# values Parsed in Loop
  getPar("INSTRUM",  _INSTRUM);			// Try and get instrument name
						// L# values Parsed in Loop
  getPar("LFILTER",  _LFILTER);			// Try and get this
  getPar("LGAIN",    _LGAIN);			// Try and get this
  getPar("LOCKPOW",  _LOCKPOW);			// Try and get this
  getPar("LOCNUC",   _LOCNUC);			// Try and get this
  getPar("LOCPHAS",  _LOCPHAS);			// Try and get this
  getPar("LOCSHFT",  _LOCSHFT);			// Try and get this
  getPar("LTIME",    _LTIME);			// Try and get this
  getPar("MASR",     _MASR);			// Try and get this
  getPar("MASRLST",  _MASRLST);			// Try and get this
  getPar("NBL",      _NBL);			// Try and get this
  getPar("NC",       _NC);			// Try and get this
  getPar("NS",       _NS);			// Try and get this
						// NUC# values Parsed in Loop
  getPar("NUCLEI",   _NUCLEI);			// Try and get this
  getPar("NUCLEUS",  _NUCLEUS);			// Try and get this
						// O# values Parsed in Loop
						// OBSCHAN# vals Parsed in Loop
  getPar("OVERFLW",  ___OVERFLW);		// Try and get this
						// P# values Parsed in Loop
  getPar("PAPS",     _PAPS);			// Try and get this
  getPar("PARMODE",  _PARMODE);			// Try and get this
						// PCPD# values Parsed in Loop
						// PHCOR# values Parsed in Loop
  getPar("PHP",      _PHP);			// Try and get this
  getPar("PH_ref",   _PH_ref);			// Try and get this
						// PL# values Parsed in Loop
  getPar("POWMOD",   _POWMOD);			// Try and get this
  getPar("PR",       __PR);			// Try and get this
						// PRECHAN# vals Parsed in Loop
  getPar("PRGAIN",   _PRGAIN);			// Try and get this
  getPar("PROBHD",   _PROBHD);			// Try and get this
  getPar("PULPROG",  _PULPROG);			// Try and get this
  getPar("PW",       _PW);			// Try and get this
  getPar("QNP",      _QNP);			// Try and get this
  getPar("RD",       _RD);			// Try and get this
  getPar("RECPH",    _RECPH);			// Try and get this
  getPar("RG",       _RG);			// Try and get this
  getPar("RO",       _RO);			// Try and get this
						// ROUTWD1# vals Parsed in Loop
						// ROUTWD2# vals Parsed in Loop
						// RSEL# values Parsed in Loop
						// S# values Parsed in Loop
  getPar("SEOUT",    _SEOUT);			// Try and get this
						// SFO# values Parsed in Loop
  getPar("SOLVENT",  _SOLVENT);			// Try and get this
						// SP# values Parsed in Loop
  getPar("SP07",     _SP07);			// Try and get this
						// SPNAM# values Parsed in Loop
						// SPOAL# values Parsed in Loop
						// SPOFFS# values Parsed in Loop
  getPar("SW",       _SW);			// Try and get spectral width *
						// SWIBOX# values Parsed in Loop
  getPar("SW_h",     _SW_h);			// Try for spectral width (Hz)
  getPar("TD",       _TD);			// Try for number of points *
  getPar("TE",       _TE);			// Try and get the temperature
						// TL# values Parsed in Loop
						// TP# values Parsed in Loop
  getPar("TP07",     _TP07);			// Try and get this
						// TPNAME# values Parsed in Loop
						// TPOAL# values Parsed in Loop
						// TPOFFS# values Parsed in Loop
  getPar("TUNHIN",  _TUNHIN);			// Try and get this
  getPar("TUNHOUT", _TUNHOUT);			// Try and get this
  getPar("TUNXOUT", _TUNXOUT);			// Try and get this
						// USERA# values Parsed in Loop
  getPar("V9",      _V9);			// Try and get this
  getPar("VALIST",  _VALIST);			// Try and get this
  getPar("VCLIST",  _VCLIST);			// Try and get this
  getPar("VD",      _VD);			// Try and get this
  getPar("VDLIST",  _VDLIST);			// Try and get this
  getPar("VPLIST",  _VPLIST);			// Try and get this
  getPar("VTLIST",  _VTLIST);			// Try and get this
  getPar("WBST",    _WBST);			// Try and get this
  getPar("WBSW",    _WBSW);			// Try and get this
						// XGAIN# values Parsed in Loop
  getPar("YMAX_a",  _YMAX_a);			// Try and get this
  getPar("YMIN_a",  _YMIN_a);			// Try and get this
  getPar("ZL1",     _ZL1);			// Try and get this
  getPar("ZL2",     _ZL2);			// Try and get this
  getPar("ZL3",     _ZL3);			// Try and get this
  getPar("ZL4",     _ZL4);			// Try and get this

//                    Now Try And Access All Parameters In Arrays

  string pbase, pname, pidx;			// String value for # points
  for(i=0; i<4; i++)
    {
    pidx = Gdec(i+1);
    pname = "HGAIN"  + pidx; 			// Try and read the base
    getPar(pname,      _HGAIN[i]); 		// spectrometer frequencies
    pname = "XGAIN"  + pidx; 			// Try and read the base
    getPar(pname,      _XGAIN[i]); 		// spectrometer frequencies
    }

  for(i=0; i<5; i++)
    {
    pidx = Gdec(i+1);
    pname = "USERA"  + pidx; 			// Try and read the base
    getPar(pname,      _USERA[i]); 		// spectrometer frequencies
    }

  for(i=0; i<8; i++)
    {
    pidx = Gdec(i+1);
    pname = "BF"     + pidx; 			// Try and read the base
    getPar(pname,      _BF[i]); 		// spectrometer frequencies
    pname = "CPDPRG" + pidx;
    getPar(pname,      _CPDPRGN[i]);
    pname = "DBL"    + pidx;			// Try and read these
    getPar(pname,      _DBL[i]);
    pname = "DBP"    + pidx;			// Try and read these
    getPar(pname,      _DBP[i]);
    pname = "DBPNAM" + pidx;			// Try and read these
    getPar(pname,      _DBPNAM[i]);
    pname = "DBPOAL" + pidx;			// Try and read these
    getPar(pname,      _DBPOAL[i]);
    pname = "DBPOFFS"+ pidx;			// Try and read these
    getPar(pname,      _DBPOFFS[i]);
    pname = "DL"     + pidx;			// Try and read these
    getPar(pname,      _DL[i]);
    pname = "DP"     + pidx;			// Try and read these
    getPar(pname,      _DP[i]);
    pname = "DPNAME" + pidx;			// Try and read these
    getPar(pname,      _DPNAME[i]);
    pname = "DPOAL" +  pidx;			// Try and read these
    getPar(pname,      _DPOAL[i]);
    pname = "DPOFFS"+  pidx;			// Try and read these
    getPar(pname,      _DPOFFS[i]);
    pname = "FQNLIST"+ pidx;			// Try and read these
    getPar(pname,      _FQNLIST[i]);
    pname = "FS"+      pidx;			// Try and read these
    getPar(pname,      _FS[i]);
    pname = "HPMOD"+   pidx;			// Try and read these
    getPar(pname,      _HPMOD[i]);
    pname = "NUC"+     pidx;			// Try and read these
    getPar(pname,      _NUC[i]);
    pname = "O"+       pidx;			// Try and read these
    getPar(pname,      _O[i]);
    pname = "S"+       pidx;			// Try and read these
    getPar(pname,      __S[i]);
    pname = "SFO"+     pidx;			// Try and read these
    getPar(pname,      _SFO[i]);
    pname = "TL"+      pidx;			// Try and read these
    getPar(pname,      _TL[i]);
    pname = "TP"+      pidx;			// Try and read these
    getPar(pname,      _TP[i]);
    pname = "TPNAME"+  pidx;			// Try and read these
    getPar(pname,      _TPNAME[i]);
    pname = "TPOAL"+   pidx;			// Try and read these
    getPar(pname,      _TPOAL[i]);
    pname = "TPOFFS"+  pidx;			// Try and read these
    getPar(pname,      _TPOFFS[i]);
    } 	

  for(i=0; i<10; i++)
    {
    pidx = Gdec(i+1);
    pname = "FCUCHAN" + pidx; 			// Try and read the constants
    getPar(pname,       _FCUCHAN[i]);
    pname = "OBSCHAN" + pidx;			// Try and read these
    getPar(pname,       _OBSCHAN[i]);
    pname = "PCPD"    + pidx;			// Try and read these
    getPar(pname,       _PCPD[i]);
    pname = "RSEL"    + pidx;			// Try and read these
    getPar(pname,       _RSEL[i]);
    }

  for(i=0; i<16; i++)
    {
    pidx = Gdec(i+1);
    pname = "PRECHAN" + pidx;			// Try and read these
    getPar(pname,       _PRECHAN[i]);
    pname = "SP"      + pidx;			// Try and read these
    getPar(pname,       _XSP[i]);
    pname = "SPNAM"   + pidx;			// Try and read these
    getPar(pname,       _SPNAM[i]);
    pname = "SPOAL"   + pidx;			// Try and read these
    getPar(pname,       _SPOAL[i]);
    pname = "SPOFFS"  + pidx;			// Try and read these
    getPar(pname,       _SPOFFS[i]);
    pname = "SWIBOX"  + pidx;			// Try and read these
    getPar(pname,       _SWIBOX[i]);
    }

  for(i=0; i<24; i++)
    {
    pidx = Gdec(i+1);
    pname = "ROUTWD1" + pidx;			// Try and read these
    getPar(pname,       _ROUTWD1[i]);
    pname = "ROUTWD2" + pidx;			// Try and read these
    getPar(pname,       _ROUTWD2[i]);
    }

  for(i=0; i<32; i++)
    {
    pidx = Gdec(i+1);
    pname = "CNST"   + pidx; 			// Try and read the constants
    getPar(pname,      _CNST[i]);
    pname = "D"      + pidx;			// Try and read the delays
    getPar(pname,      _D_[i]);
    pname = "GPNAME" + pidx;			// Try and read these
    getPar(pname,      _GPNAME[i]);
    pname = "GPX"    + pidx;			// Try and read these
    getPar(pname,      _GPX[i]);
    pname = "GPY"    + pidx;			// Try and read these
    getPar(pname,      _GPY[i]);
    pname = "GPZ"    + pidx;			// Try and read these
    getPar(pname,      _GPZ[i]);
    pname = "IN"     + pidx;			// Try and read these
    getPar(pname,      _IN[i]);
    pname = "INP"    + pidx;			// Try and read these
    getPar(pname,      _INP[i]);
    pname = "P"      + pidx;			// Try and read these
    getPar(pname,      ___P[i]);
    pname = "PHCOR"  + pidx;			// Try and read these
    getPar(pname,      _PHCOR[i]);
    pname = "PL"     + pidx;			// Try and read these
    getPar(pname,      _PL[i]);
    }

//          Now Re-Access All REQUIRED Parameters To Insure We Know Them

  if(!getPar("BYTORDA",_BYTORDA)) _BYTORDA = 0; // Try and get byte order
  if(!getPar("TD",_TD,21,1))			// Try and get # points
    XWinAcqParfatality(2, "TD");		// We MUST know this value
  for(i=0; i<8; i++)				// Try and read in all
    { 						// spectrometer freq. offsets
    pname = "SFO"+Gdec(i+1);
    getPar(pname,_SFO[i]);
    }
  if(!getPar("SFO1",_SFO[0],23,1))		// Insure at least SFO
    XWinAcqParfatality(2, "SFO1");		// has been read in
  if(!getPar("SW",  _SW,  24,1))		// Try to get spect. width
    XWinAcqParfatality(2, "SW");		//   Problems with pname
  return true;
  }

// ____________________________________________________________________________
// D                       XWinAcqPar Output Functions
// ____________________________________________________________________________

/* These function allow for output of NMR parameters directly into a Bruker
   XWinNMR ASCII parameter file {acqu(s), acqu2(s)}. We don't write everything
   that Bruker does since users would have to set all of them if so. Note that
   the value of PARMODE, that which flags the dimension of the acquisiton, 
   will affect output somewhat. 0=1D, 1=2D, .... for PARMODE.                */  

int XWinAcqPar::writeAPar(const string& name, int warn)
  { parfile = name; return writeAPar(warn); }

int XWinAcqPar::writeAPar(int warn) const
  {  
  ofstream ofstr(parfile.c_str());		// Open ASCII file for output
  if(!ofstr.good())				// If file bad then exit
     {
     if(warn)
       {
       XWinAcqParerror(2, 1);			// Output filestream problems
       XWinAcqParerror(25,1);			// Can't write parameters out
       if(warn==1) XWinAcqParerror(21,parfile,1);// Can't write parameters out
       else        XWinAcqParfatality(21,parfile);// with out without die
       } 
     return 0;
     }  
  string nn("##");			// Bruker ## parameter line start
  string nns("##$");			// Bruker ##$ parameter line start
  string ss("$$");			// Bruker $$ parameter line start

/*   Below I've Put In Lines For Outputting A Whole Lot of Information BUT
     Most Of It Is Not Used By GAMMA On Output.  If I've Left Out Something
     Deemed Essential, Just Un-Comment Or Add The Appropriate Output Lines
     & Adjust The Class In Order To Be Able To Set the Parameter Externally  */

  ofstr << nn  << "TITLE= "    << _TITLE    << "\n";	// Starting title
  ofstr << nn  << "JCAMPDX= "  << _JCAMPDX  << "\n"; 	// JCAMP version
  ofstr << nn  << "DATATYPE= " << _DATATYPE << "\n";	// Data type in file
  if(!_PARMODE)
  {
  ofstr << nn  << "ORIGIN= "   << _ORIGIN   << "\n";
  ofstr << nn  << "OWNER= "    << _OWNER    << "\n";
  }
  ofstr << ss  << " " << _DAY;				// Time and date
  ofstr << ss  << " "        << _NAME   << "\n";	// Full file name
  ofstr << nns << "AQSEQ= "  << _AQSEQ  << "\n";	// Acquisiton sequence
  ofstr << nns << "AQ_mod= " << _AQ_mod << "\n";	// Acquisiton mode
							// qf=   single channel  
							// qseq= quadrature sequential (2)
							// qsim= quad. simultaneous (1)**
							// qdig= digital quadrature
  ofstr << nns << "AUNM= <"  << _AUNM    <<">\n"; 	// To start AU higher level
  int i=0, j=0, k=0;
  int l=8;
  if(_PARMODE) l=_PARMODE+1;
  for(i=0; i<8; i++) 					// Channel base frequencies
    ofstr << nns << "BF" << i+1
               <<        "= " << _BF[i]   << "\n";
  ofstr << nns << "BYTORDA= " << _BYTORDA << "\n";	// Set byte order flag 
  if(!_PARMODE)
  {
  ofstr << nns << "CFDGTYP= " << _CFDGTYP << "\n";	// Set 
  ofstr << nns << "CFGGTYP= " << _CFRGTYP << "\n";	// Set 
  ofstr << nns << "CNST= (0..31)\n";			// 32 flt constants
  for(i=0; i<31; i++)
     ofstr << _CNST[i] << " ";   ofstr    << "\n";
  ofstr << nns << "CPDPRG= <"  << _CPDPRG  << ">\n";	// AMX/ARC cpd on decoupler
  for(i=0; i<8; i++) 					// Pulse program commands
    ofstr << nns << "CPDPRG"  << Gdec(i+1)
	              << "= <"<<_CPDPRGN[i]<<">\n";
  ofstr << nns << "CPDPRGB= <"<< _CPDPRGB << ">\n";	// AMX/ARC cpd on 2nd decoupler
  ofstr << nns << "CPDPRGT= <"<< _CPDPRGT << ">\n";	// AMX/ARC cpd on transitter
  }
  ofstr << nns << "D= (0..31)\n";			// Delay parameters (sec)
  for(i=0; i<31; i++)
    ofstr << _D_[i] << " ";          ofstr << "\n";
  ofstr << nns << "DATE= "    << _DATE    << "\n";	// Output current time/date
  if(!_PARMODE)
  {
  ofstr << nns << "DBL= (0..7)\n";			// Low power (db) level (sq.pul.trans)
  for(i=0; i<8; i++)					// AMX/ARX channel 3 (2nd decoupler)
    ofstr << _DBL[i] << " ";        ofstr << "\n";
  ofstr << nns << "DBP= (0..7)\n";			// Shaped pulse parameter table
  for(i=0; i<8; i++)					// AMX/ARX channel 3 (2nd decoupler)
    ofstr << _DBP[i] << " ";
  ofstr << "\n";
  ofstr << nns << "DBP07= " << _DBP07 << "\n";		// Shaped pulse parameter table
  for(i=0; i<8; i++)					// AMX/ARX channel 3 (2nd decoupler)
    ofstr << nns << "DBPNAM" << Gdec(i)
	            << "= <" << _DBPNAM[i] << ">\n";
  ofstr << nns << "DBPOAL= (0..7)\n";
  ofstr << "0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5\n";
  ofstr << nns << "DBPOFFS= (0..7)\n";
  ofstr << "0 0 0 0 0 0 0 0\n";
  }
  ofstr << nns << "DE=      " << _DE      << "\n";	// Pre-scan delay (usec) 
  ofstr << nns << "DECBNUC= <" << _DECBNUC << ">\n"; 	// Channel 3 nucleus?
  if(!_PARMODE)
  ofstr << nns << "DECIM=   " << _DECIM   << "\n";	// Decimation rate (digital filter)
  ofstr << nns << "DECNUC= <" << _DECNUC  << ">\n"; 	// Channel 2 nucleus
  ofstr << nns << "DECSTAT= " << _DECSTAT << "\n"; 	// Not documented
  if(!_PARMODE)
  ofstr << nns << "DIGMOD= "  << _DIGMOD  << "\n"; 	// Digitization mode (Avance only)
  ofstr << nns << "DIGTYP=  " << _DIGTYP  << "\n";	// Digitizer type
  if(!_PARMODE)						// 0)analog 1)digital 2)homodecoupling digitial
  {
  ofstr << nns << "DL= (0..7)\n";			// Low power db levels (sq. pul. trans)
  ofstr << "0 120 120 120 120 120 120 120\n";		// AMX/ARX decoupler (2nd) channel
  }
  ofstr << nns << "DP= (0..7)\n";
  ofstr << "150 150 150 150 150 150 150 150\n";
  if(!_PARMODE)						// 0)analog 1)digital 2)homodecoupling digitial
  {
  ofstr << nns << "DP07= "    << _DP07    << "\n";	// Shaped pulse parameter table
  for(i=0; i<8; i++)					// AMX/ARX channel 3 (2nd decoupler)
    ofstr << nns << "DPNAME" << Gdec(i)
	              << "= <" << _DPNAME[i] << ">\n";
  ofstr << nns << "DPOAL= (0..7)\n";			// Undocumented
  ofstr << "0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5\n";
  ofstr << nns << "DPOFFS= (0..7)\n";			// Undocumented
  ofstr << "0 0 0 0 0 0 0 0\n";
  ofstr << nns << "DQDMODE= " << _DQDMODE << "\n";	// Undocumented
  }
  ofstr << nns << "DR= "      << _DR      << "\n";	// Digitizer resolution
  ofstr << nns << "DS= "      << _DSc      << "\n";	// Number of dummy scans 
  ofstr << nns << "DSLIST= <" << _DSLIST  << ">\n";	// Data set list file
  if(!_PARMODE)						// 0)analog 1)digital 2)homodecoupling digitial
  {
  ofstr << nns << "DSPFIRM= " << _DSPFIRM << "\n";	// Firmware dig. filtering
  ofstr << nns << "DSPFVS= "  << _DSPFVS  << "\n";	// Undocumented
  ofstr << nns << "DTYPA= "   << _DTYPA   << "\n";
  }
  ofstr << nns << "EXP= <"    << _EXP     <<">\n";	// Experiment perfomed
  ofstr << nns << "F1LIST= <" << _F1LIST  <<">\n";	// Frequency list file
  ofstr << nns << "F2LIST= <" << _F2LIST  <<">\n";	// Frequency list file
  ofstr << nns << "F3LIST= <" << _F3LIST  <<">\n";	// Frequency list file
  if(!_PARMODE)						// 0)analog 1)digital 2)homodecoupling digitial
  {
  ofstr << nns << "FCUCHAN= (0..9)\n";			// FCU for channel fx (avance)
  ofstr << "0 0 0 0 0 0 0 0 0 0\n";			// having 8 possible FCUs
  }
  ofstr << nns << "FL1=    " << _FL1      << "\n";	// F19 decoup. high pwr (db)
  ofstr << nns << "FL2=    " << _FL2      << "\n"; 	// F19 decoup. high pwr (db)
  ofstr << nns << "FL3=    " << _FL3      << "\n";	// F19 decoup. high pwr (db)
  ofstr << nns << "FL4=    " << _FL4      << "\n";	// F19 decoup. high pwr (db)
  ofstr << nns << "FOV=    " << _FOV      << "\n";
  if(!_PARMODE)						// 0)analog 1)digital 2)homodecoupling digitial
  {
  for(i=0; i<8; i++) 					// Frequency list files
    ofstr << nns << "FQ" << Gdec(i+1) 			// Avance spectrometers
	  << "LIST= <" << _FQNLIST[i] << ">\n";
  }
  ofstr << nns << "FS= (0..7)\n";			// Undocumented
  ofstr << "83 83 83 83 83 83 83 83 \n";
  if(!_PARMODE)						// 0)analog 1)digital 2)homodecoupling digitial
  ofstr << nns << "FTLPGN= " << _FTLPGN << "\n";	// Undocumented
  double swh = _SW*_SFO[0];
  ofstr << nns << "FW= " << 1.25*swh << "\n";		// Filter width in Hz (1.25*SWH)
  if(!_PARMODE)						// 0)analog 1)digital 2)homodecoupling digitial
  {
  ofstr << nns << "GP031= 0\n";				// Gradient parameter table
  for(i=0; i<10; i++)					// Output gradient parameters
    {
    ofstr << nns << "GPNAME" << Gdec(i) << "= <" << _GPNAME[i] << ">\n";
    for(j=0; j<10 && i>0 && i<4; j++)
      {
      k=10*i+j;
      if(k < 32)
        ofstr << nns << "GPNAME" << Gdec(k) << "= <" << _GPNAME[k] << ">\n";
      }
    }
  ofstr << nns << "GPX= (0..31)\n";
  ofstr << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n";
  ofstr << nns << "GPY= (0..31)\n";
  ofstr << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n";
  ofstr << nns << "GPZ= (0..31)\n";
  ofstr << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n";
  ofstr << nns << "GRDPROG= <"    << _GRDPROG << ">\n";
  ofstr << nns << "HGAIN= (0..3)\n";
  ofstr << "0 0 0 0\n";
  ofstr << nns << "HL1= "         << _HL1    << "\n";	// Decoupler high-power level
  }
  ofstr << nns << "HL2= "         << _HL2    << "\n";	// Decoupler high-power level
  ofstr << nns << "HL3= "         << _HL3    << "\n";	// Decoupler high-power level
  ofstr << nns << "HL4= "         << _HL4    << "\n";	// Decoupler high-power level
  if(!_PARMODE)						// 0)analog 1)digital 2)homodecoupling digitial
  {
  ofstr << nns << "HOLDER= "      << _HOLDER << "\n";
  ofstr << nns << "HPMOD= (0..7)"            << "\n";
  ofstr << "0 0 0 0 0 0 0 0\n";
  ofstr << nns << "HPPRGN= "      << _HPPRGN << "\n";	// Gain for HPPR preamp
  }
  ofstr << nns << "IN= (0..31)\n";			// Increments for delays
  for(i=0; i<31; i++)
    ofstr << _IN[i] << " ";
  ofstr << "\n";
  if(!_PARMODE)						// 0)analog 1)digital 2)homodecoupling digitial
  {
  ofstr << nns << "INP= (0..31)\n";			// Pulse increment values
  for(i=0; i<31; i++)
    ofstr << _INP[i] << " ";
  ofstr << "\n";
  }
  ofstr << nns << "INSTRUM= <" << _INSTRUM << ">\n";		// For autocalibration
  ofstr << nns << "L= (0..31)\n";		// Loop counter params
  ofstr << "1 28 1 1 1 28 1 1 1 1 1 1 1 1 "
        << "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 "
        << "1 1 1\n";
  if(!_PARMODE)						// 0)analog 1)digital 2)homodecoupling digitial
  {
  ofstr << nns << "LFILTER= " << _LFILTER << "\n";	// Lock filter
  ofstr << nns << "LGAIN= "   << _LGAIN   << "\n";	// Lock gain
  ofstr << nns << "LOCKPOW= " << _LOCKPOW << "\n";	// Lock power
  }
  ofstr << nns << "LOCNUC= <" << _LOCNUC  <<">\n";	// Lock nucleus
  if(!_PARMODE)						// 0)analog 1)digital 2)homodecoupling digitial
  {
  ofstr << nns << "LOCPHAS= " << _LOCPHAS << "\n";	// Undocumented
  ofstr << nns << "LOCSHFT= " << _LOCSHFT << "\n";	// Undocumented
  ofstr << nns << "LTIME= "   << _LTIME   << "\n";	// BSMS Lock Parameter
  ofstr << nns << "MASR= "    << _MASR    << "\n";	// Undocumented
  ofstr << nns << "MASRLST= <"<< _MASRLST << ">\n";	// Undocumented
  }
  ofstr << nns << "NBL= "     << _NBL     << "\n";	// Number of buffers
  ofstr << nns << "NC= "      << _NC      << "\n";	// Undocumented
  ofstr << nns << "NS= "      << _NS      << "\n";	// Number of scans
//  if(!_PARMODE)						// 0)analog 1)digital 2)homodecoupling digitial
//  {
  for(i=0; i<8; i++)					// Channel [1,8] nuclei
    ofstr << nns << "NUC" << i+1
          << "= <" << _NUC[i] << ">\n";
  ofstr << nns << "NUCLEI= "  << _NUCLEI << "\n";	// Set up nuclei
//  }
  ofstr << nns << "NUCLEUS= <" 				// Channel 1 nucleus
        << _NUCLEUS << ">\n";
  l = 8;
  if(_PARMODE) l = _PARMODE+8;
  for(i=0; i<l; i++)
    ofstr << nns << "O" << i+1 << "= " 			// Channel i+1 offset
          << _O[i] << "\n"; 				// frequency (Hz)
  if(!_PARMODE)						// 0)analog 1)digital 2)homodecoupling digitial
  {
  ofstr << nns << "OBSCHAN= (0..9)\n";			// Observation channel (avance)
  ofstr << "0 0 0 0 0 0 0 0 0 0\n";			// channel 1 is default
  ofstr << nns << "OVERFLW= 0\n";			// Undocumented
  }
  ofstr << nns << "P= (0..31)\n";			// Pulse lengths (usec)
  for(i=0; i<31; i++)
    ofstr << ___P[i] << " ";
  ofstr << "\n";
  ofstr << nns << "PAPS= "    << _PAPS    << "\n";	// Phase cycling mode
  ofstr << nns << "PARMODE= " << _PARMODE << "\n";	// Dimension of acq. data (1,2,3D)
  if(!_PARMODE)						// 0)analog 1)digital 2)homodecoupling digitial
  {
  ofstr << nns << "PCPD= (0..9)\n";
  ofstr << "100 100 100 100 100 100 100 100 100 100\n";
  ofstr << nns << "PHCOR= (0..31)\n";			// Phase correction angles
  ofstr << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
        << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n";
  ofstr << nns << "PHP= "    << _PHP << "\n";		// Preamp module select (AMX)
  ofstr << nns << "PH_ref= " << _PH_ref << "\n";	// Reference phase
  ofstr << nns << "PL= (0..31)\n";			// 32 Power levels (db)
  ofstr << "120 120 120 120 120 120 120 "
        << "120 120 120 120 120 120 120 "
        << "120 120 120 120\n";
  ofstr << "120 120 120 120 120 120 120 "
        << "120 120 120 120 120 120 120\n";
  ofstr << nns << "POWMOD= " << _POWMOD << "\n";// Power mode {high,low,linear}
  }
  ofstr << nns << "PR= "     << __PR    << "\n";	// Undocumented
  if(!_PARMODE)						// 0)analog 1)digital 2)homodecoupling digitial
  {
  ofstr << nns << "PRECHAN= (0..15)\n";		// Preamplifier channel (avance)
  ofstr << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n";
  ofstr << nns << "PRGAIN= " << _PRGAIN << "\n";		// High power preamp gain
  ofstr << nns << "PROBHD= <"<< _PROBHD << ">\n";		// Probehead
  }
  ofstr << nns << "PULPROG= <"			// Pulse program
        << _PULPROG << ">\n";
  ofstr << nns << "PW= "    << _PW    << "\n";	// Undocumented
  if(!_PARMODE)
  ofstr << nns << "QNP= "   << _QNP   << "\n";	// QNP selection [1-4]
  ofstr << nns << "RD= "    << _RD    << "\n";	// Undocumented
  if(!_PARMODE)
  ofstr << nns << "RECPH= " << _RECPH << "\n";	// Undocumented
  ofstr << nns << "RG= "    << _RG    << "\n";	// Receiver gain
  ofstr << nns << "RO= "    << _RO    << "\n";	// Spinner rotation frequency
  ofstr << nns << "ROUTWD1= (0..23)\n";		// Router control words
  ofstr << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
        << "0 0 0 0 1 1 0 0\n";
  ofstr << nns << "ROUTWD2= (0..23)\n";		// Router control workds
  ofstr << "0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 "
        << "0 0 1 0 1 1 0 0\n";
  if(!_PARMODE)
  {
  ofstr << nns << "RSEL= (0..9)\n";		// Transmitter select (avance)
  ofstr << "0 0 0 0 0 0 0 0 0 0\n";
  }
  ofstr << nns << "S= (0..7)\n";		// Decoupler power
  ofstr << "83 1 26 26 26 83 83 83\n";
  ofstr << nns << "SEOUT= " << _SEOUT << "\n";	// SE 451 receiver unit (HR, BB)
  l = 8;
  if(_PARMODE) l= _PARMODE+1;
  for(i=0; i<l; i++)
    ofstr << nns << "SFO" << i+1 << "= " 	// Channel i+1 irradiation freq.
          << _SFO[i] << "\n";
  ofstr << nns << "SOLVENT= <" 			// Solvent used
        << _SOLVENT << ">\n";
  ofstr << nns << "SP= (0..15)\n";		// Shape pulse parameters
  ofstr << "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n";
  if(!_PARMODE)
  {
  ofstr << nns << "SP07= "   << _SP07     << "\n";
  ofstr << nns << "SPNAM= <" << _SPNAM[0] << ">\n";
  ofstr << nns << "SPNAM1=<" << _SPNAM[1] << ">\n";
  for(i=10; i<16; i++)
    ofstr << nns  << "SPNAM<" << Gdec(i)
          << "= " << _SPNAM[i] << ">\n";
  for(i=2; i<10; i++)
    ofstr << nns  << "SPNAM<" << Gdec(i)
          << "= " << _SPNAM[i] << ">\n";
  ofstr << nns << "SPOAL= (0..15)\n";
  ofstr << "0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5\n";
  ofstr << nns << "SPOFFS= (0..15)\n";
  ofstr << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n";
  }
  ofstr << nns << "SW= " << _SW << "\n";	// Spectral width (PPM)
  if(!_PARMODE)
  {
  ofstr << nns << "SWIBOX= (0..15)\n";		// Cannel switching (avance)
  ofstr << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n";	// 
  }
  ofstr << nns << "SW_h= " << _SW_h << "\n";    // Spectral width (Hz)
  ofstr << nns << "TD= " << _TD << "\n";	// Set the TD value
  ofstr << nns << "TE= " << _TE << "\n";	// Temperature in K
  if(!_PARMODE)
  {
  ofstr << nns << "TL= (0..7)\n";		// Low power (db) levels (sq pul trans)
  ofstr << "0 120 120 120 120 120 120 120\n";	// for AMX/ARC transmitter (chanel 1)
  ofstr << nns << "TP= (0..7)\n";
  ofstr << "150 150 150 150 150 150 150 150\n";
  ofstr << nns << "TP07= " << _TP07 << "\n";	// Shaped pulse parameter table
  for(i=0; i<8; i++)
    ofstr << nns  << "TPNAME" << Gdec(i)
          << "= <" << _TPNAME[i] << ">\n";
  ofstr << nns << "TPOAL= (0..7)\n";
  ofstr << "0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5\n";
  ofstr << nns << "TPOFFS= (0..7)\n";
  ofstr << "0 0 0 0 0 0 0 0\n";
  ofstr << nns << "TUNHIN= "  << _TUNHIN  << "\n";
  ofstr << nns << "TUNHOUT= " << _TUNHOUT << "\n";
  ofstr << nns << "TUNXOUT= " << _TUNXOUT << "\n";
  for(i=0; i<5; i++)
    ofstr << nns  << "USERA" << Gdec(i+1) 
	  << "= <" << _USERA[i] << ">\n";
  }
  ofstr << nns << "V9= "     << _V9     << "\n";	// Variation in random delay
  if(!_PARMODE)
  ofstr << nns << "VALIST= <" << _VALIST << ">\n";
  ofstr << nns << "VCLIST= <" << _VCLIST << ">\n";	// Variable loop counter list
  ofstr << nns << "VD= "     << _VD     << "\n";	// Un-documented
  ofstr << nns << "VDLIST= <" << _VDLIST << ">\n";	// Variable delay list file
  ofstr << nns << "VPLIST= <" << _VPLIST << ">\n";	// Variable pulse list file
  ofstr << nns << "VTLIST= <" << _VTLIST << ">\n";	// Variable temp list file
  if(!_PARMODE)
  {
  ofstr << nns << "WBST= "   << _WBST   << "\n";	// # of wobble steps
  ofstr << nns << "WBSW= "   << _WBSW   << "\n";	
  ofstr << nns << "XGAIN= (0..3)\n";
  for(i=0; i<4; i++) 
    ofstr << _XGAIN[i] << " ";    ofstr << "\n"; 
  ofstr << nns << "XL= "     << _XL     << "\n";
  ofstr << nns << "YL= "     << _YL     << "\n";
  }
  ofstr << nns << "YMAX_a= " << _YMAX_a << "\n";	// Un-documented 
  ofstr << nns << "YMIN_a= " << _YMIN_a << "\n";	// Un-documented
  if(!_PARMODE)
  {
  ofstr << nns << "ZL1= "    << _ZL1    << "\n";	// 4th channel low power lev
  ofstr << nns << "ZL2= "    << _ZL2    << "\n";	// 4th channel low power lev
  ofstr << nns << "ZL3= "    << _ZL3    << "\n";	// 4th channel low power lev
  ofstr << nns << "ZL4= "    << _ZL4    << "\n";	// 4th channel low power lev
  }
  ofstr << nn  << "END= \n";				// End of acqus file
  return 1;
  }  

// ____________________________________________________________________________
// E                 XWinAcqPar Parameter Access Functions
// ____________________________________________________________________________
 
/* These functions allow access to the Bruker XWinNMR parameter set, not just
  those parameters which are directly handled within this class.
 
                       INHERITED FROM BASE CLASS XWinPSet
 
int XWinAcqPar::getPar(const string& pn,int& val,   int id,int wrn) const
int XWinAcqPar::getPar(const string& pn,double& val,int id,int wrn) const
int XWinAcqPar::getPar(const string& pn,string& val,int id,int wrn) const
ParameterSet XWinAcqPar::getPSet() const                                     */


// ____________________________________________________________________________
// F                      XWinAcqPar Auxiliary Functions
// ____________________________________________________________________________
 
 /* These are just a group of functions that return strings indicating what
    the Bruker parameter values mean.  I just add things as I learn them here
    so that GAMMA output can remind me what all of these parameters are...   */


string XWinAcqPar::AQ_modS() const
  {
  string ret;
  switch(_AQ_mod)
    {
    case 0:  ret = "qf";      break;
    case 1:  ret = "qseq";    break;
    case 2:  ret = "qsim";    break;
    case 3:  ret = "qdig";    break;
    default: ret = "Unknown";
    }
  return ret;
  }

string XWinAcqPar::BYTORDAS() const
  {
  string ret;
  switch(_BYTORDA)
    {
    case 0:  ret = "Little Endian"; break;
    case 1:  ret = "Big Endian";    break;
    default: ret = "Unknown";
    }
  return ret;
  }


string XWinAcqPar::DECSTATS() const
  {
  string ret;
  switch(_DECSTAT)
    {
    case 0:  ret = "PO"; break;
    default: ret = "Unknown";
    }
  return ret;
  }

string XWinAcqPar::DIGTYPS() const
  {
  string ret;
  switch(_DIGTYP)
    {
    case 0:  ret = "slow";            break;
    case 1:  ret = "16bit";           break;
    case 2:  ret = "fast";            break;
    case 3:  ret = "BC132-12";        break;
    case 4:  ret = "BC132-16";        break;
    case 5:  ret = "FADC (=BC133)";   break;
    case 6:  ret = "HADC (=HRD16)";   break;
    case 7:  ret = "SADC";            break;
    default: ret = "Unknown";
    }
  return ret;
  }

string XWinAcqPar::DSS() const
  { return string ("Dummy Scans"); }


string XWinAcqPar::HPPRGNS() const
  {
  string ret;
  switch(_HPPRGN)
    {
    case 0:  ret = "normal";     break;
    case 1:  ret = "plus";    break;
    default: ret = "Unknown";
    }
  return ret;
  }

string XWinAcqPar::PAPSS() const
  {
  string ret;
  switch(_PAPS)
    {
    case 0:  ret = "CP";    break;
    case 1:  ret = "AP";    break;
    case 2:  ret = "QP";    break;
    default: ret = "Unknown";
    }
  return ret;
  }

string XWinAcqPar::PARMODES() const
  {
  string ret;
  switch(_PARMODE)
    {
    case 0:  ret = "1D";          break;
    case 1:  ret = "2D";          break;
    case 2:  ret = "3D";          break;
    default: ret = "Unknown";	  break;
    }
  return ret;
  }


string XWinAcqPar::POWMODS() const
  {
  string ret;
  switch(_POWMOD)
    {
    case 0:  ret = "low";     break;
    case 1:  ret = "high";    break;
    case 2:  ret = "linear";  break;
    default: ret = "Unknown"; break;
    }
  return ret;
  }

string XWinAcqPar::PRGAINS() const
  {
  string ret;
  switch(_PRGAIN)
    {
    case 0:  ret = "high"; break;
    case 1:  ret = "low";    break;
    default: ret = "Unknown";
    }
  return ret;
  }


string XWinAcqPar::SEOUTS() const
  {
  string ret;
  switch(_SEOUT)
    {
    case 0:  ret = "HR"; break;
    case 1:  ret = "BB";    break;
    default: ret = "Unknown";
    }
  return ret;
  }


string XWinAcqPar::SWS() const
  { return string ("Sweep Width (ppm)"); }


#endif								// XWinAcqPar.cc
