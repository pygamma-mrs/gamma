/* XWinProcPar.h ************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      XWinProcPar                                   Interface		**
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
** sets. This class embodies a "set" of Bruker processing parameters,	**
** parameters that might be founde in their ASCII parameter files for	**
** any acquisiton: proc, procs, proc2, proc2s,.....  This is almost	**
** just a structure rather than a class because it is so primitive.	**
** We allow higher XWinNMR classes that deal with processing parameter	**
** files to be our friend so that they have direct access to all such	**
** parameters.								**
**									**
*************************************************************************/

#ifndef   XWinProcPar_H_ 		// Is file already included?
#  define XWinProcPar_H_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the implementation
#  endif
 
#include <string>			// Include libstdc++ strings
#include <Basics/ParamSet.h>		// Include GAMMA parameter sets
#include <GamIO/XWinPSet.h>		// Include Bruker parameter parsing

class XWinProcPar: public XWinPSet

// ----------------------------------------------------------------------------
// ------------------------------- STRUCTURE ----------------------------------
// ----------------------------------------------------------------------------
 
  {
  friend class XWinProcs;		// Allow class XWinProcs full access
  friend class XWinProc2s;		// Allow class XWinProc2s full access
  std::string       parfile;                 // Parameter file name (base)
  std::string       _TITLE;			// Parameter file title
  double       _JCAMPDX;		// JCAMP version (5.0)
  std::string       _DATATYPE;		// Type of file (Parameter Values)
  std::string       _ORIGIN;			// Where the file came from
  std::string       _OWNER;			// Who made the file
  std::string       _DATE;			// When the file was made
  std::string       _NAME;			// The full file name
  double       _ABSF1;			// BC limit
  double       _ABSF2;			// BC limit
  int          _ABSG;			// BC polynomial
  int          _ABSL;			// Integral detect
// sosi MSVC++ uses _ALPHA so changed this name to _XALPHA
  int          _XALPHA;			// Quad correct
  int          _AQORDER;		// 3D process order
  int          _ASSFAC;			// Vertical scaling
  int          _ASSFACI;		// Unknown
  int          _ASSFACX;		// Unknown
  int	       _ASSWID;			// Search interval
  std::string       _AUNMP;			// For AU
  double       _AZFE;			// Intetral ends
  double       _AZFW;			// Multiplet integral
  int          _BCFW;			// Unknown
  int          _BC_mod;			// BC modify
  int          _BYTORDP;		// Byte ordering
  int          _COROFFS;		// Unknown
  int          _DATMOD;			// Add spectra
  int 	       _DC;			// For spectra add
  std::string  _DFILT;			// Filter file
  int          _DTYPP;			// Unknown
  double       _FCOR;			// 1st pt multiply
  int          _FTSIZE;			// Unknown
  int          _FT_mod;			// FT aquire mod
  int          _GAMMA;			// Quad correct. 
  int          _GB;			// Gauss. multiply
  int          _INTBC;			// Integration param.
//  double       _INTSCL;			// Integral compare
  int          _INTSCL;			// Integral compare
  int          _ISEN;			// Small integral
  double       _LB;			// Line broadening
  int          _LEV0;			// Contouring param
  int          _LPBIN;			// Linear prediction
  int          _MAXI;			// Peak picking
  int          _MC2;			// 2D tranform type
  int          _MEAN;			// Unknown
  int          _ME_mod;			// Linear prediciton
  double       _MI;			// Peak picking
  int          _NCOEF;			// LP coefficients
  int          _NC_proc;		// 2D processing
  int          _NLEV;			// Contour levels
  double       _NOISF1;			// SN ratio
  double       _NOISF2;			// SN ratio
  int          _NSP;			// Left shift
  int          _NTH_PI;			// Unknown
  int          _NZP;			// For zp
  double       _OFFSET;			// 1st pt shift (ppm)
  double       __PC;			// Peak search
  double       _PHC0;			// Zero order phase
  double       _PHC1;			// 1st order phase
  int          _PH_mod;			// Unknown
  int          _PKNL;			// 5th order phase
  int          _PPARMOD;		// Parameter display
  int          _PSCAL;			// Peak picking
  int          _PSIGN;			// Unknown
  int          _REVERSE;		// Spectrum Reverse
  double       _SF;			// Spect. Freq.
  int          _SI;			// Size
  double       _SIGF1;			// Signal/Noise
  double       _SIGF2;			// Signal/Noise
  double       _SINO;			// Signal/Noise
  int          _SIOLD;			// Unknown
  std::string  _SREGLST;		// Unknown
  int          _SSB;			// Sin Window
  int          _STSI;			// 2D storage
  int          _STSR;			// 2D storage
  double       _SW_p;			// Sweep Width (ppm)
  int          _SYMM;			// Symmetrization
  int          _S_DEV;			// Std. Dev.
  int          _TDeff;			// Effective FFT size
  int          _TDoff;			// FFT offset
  std::string       _TI;			// Unknown
  std::string       _TILT;			// Tilt flag (2D)
  int          _TM1;			// Unknown
  int          _TM2;			// Unknown
  int          _TOPLEV;			// Contour top level
  std::string       _USERP1;			// User Permission?
  std::string       _USERP2;			// Unknown
  std::string       _USERP3;			// Unknown
  std::string       _USERP4;			// Unknown
  std::string       _USERP5;			// Unknown
  int          _WDW;			// Window function
  int          _XDIM;			// Submatrix size
  int          _YMAX_p;			// 2D Max
  int          _YMIN_p;			// 2D Min
 
// ----------------------------------------------------------------------------
// ----------------------------- PRIVATE FUNCTIONS ----------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________       
// i                      XWinNMR Procs File Error Handling
// ____________________________________________________________________________

/* These functions take care of any errors encountered when reading, writing,   
   and setting parameters in Bruker processing parameter files.
 
	Input		ProcPar : XWinNMR procs parameters (this)
			eidx    : Error index
                        pname   : Additional error message
                        noret   : Flag for linefeed (0=linefeed)
        Output          void    : An error message is output                 */

         void XWPPerror(int eidx, int noret=0) const;
         void XWPPerror(int eidx, const std::string& pname, int noret=0) const;
volatile void XWPPfatal(int eidx) const;
volatile void XWPPfatal(int eidx, const std::string& pname) const;

// ____________________________________________________________________________
// ii                 XWinNMR Procs File Default Parmameters
// ____________________________________________________________________________
 
 void SetDefaults(const  std::string& fname);
 void SetDefaults1(const std::string& fname);
 void SetDefaults2(const std::string& fname);
 void Copy(const XWinProcPar& XWPP);

// ____________________________________________________________________________ 
// iii                  XWinNMR ProcPar Inter-Related Parameters
// ____________________________________________________________________________ 
 
/* Some of the parameters in an XWinNMR processing parameter files dependent on
   on another.  Thus when one parameter is specified or altered, on or more 
   other parameters must be adjusted too.  These functions are meant to take 
   care of such details.                                                     */ 
  
void SetOffset();            // Sets OFFSET, 1st point shift (Hz) 


public:
// ____________________________________________________________________________ 
// A             XWinProcPar Parameter File Constructors, Destructor
// ____________________________________________________________________________
 
/* These are the constructors of the class handling Bruker XWinNMR output
   parameter files.  This doesn't do anything in particular, it is the write
   functions that perform the work.                                          */

        XWinProcPar();
        XWinProcPar(const std::string& name, int type=1);
        XWinProcPar(const XWinProcPar& XWP);
virtual ~XWinProcPar();
void    operator= (const XWinProcPar& XWAP);

// ____________________________________________________________________________
// B                  XWinProcPar Parameter Access Functions
// ____________________________________________________________________________
 
/* These functions allow users to set any more important parameters that
   these files contain.  The two primary values are the spectrum point size
   (SI) and the data byte order (BYTORDP).                                   */

std::string parname()   const;	// ASCII file name
int    BYTORDP()   const;	// Binary byte order
int    FT_mod()    const;	// How FFT is performed
double LB()        const;	// Line Broadening
int    MC2()       const;	// FT type on t1
double OFFSET()    const;	// Spectrum offset
double PHC0()      const;	// Zero order phase
double PHC1()      const;	// 1st order phase
int    PH_mod()    const;  // Set phase type
int    REVERSE()   const;	// Plot spectrum reverse
double SF()        const;	// Spectrometer frequency
int    SI()        const;	// Data size (re+im)
int    SSB()       const;	// Sine bell
int    STSI()      const;	// Strip size
int    STSR()      const;	// Strip start
double SW_p()      const;	// Spectral width (PPM)
double TDeff()     const;	// Effective FFT size
int    WDW()       const;	// Window function

void BYTORDP(int bo);	   // Set binary byte order
void GB(int gb);		   // Set gaussian broadening
void LB(double lb);	   // Set line broadening
void FT_mod(int ft);	   // Set transform type
void FT_mod(const std::string& ft);// Set transform type
void MC2(int mc);		   // Set acquisiiton type on t1
void MC2(const std::string& mc);   // Set acquisiiton type on t1
void OFFSET(double off);	   // Spectrum offset
void PHC0(double ph0);	   // Set 0th phase correction
void PHC1(double ph1);	   // Set 1st phase correction
void PH_mod(int phm);	   // Set phase type
void PH_mod(const std::string& p); // Set phase type
void REVERSE(int yn);	   // Set spectrum reverse
void PPARMOD(int pm);	   // Set data dimension
void SI(int si);		   // Set data size (re+im)
void SF(double SF);	   // Set spectrometer freq.
void SSB(int sb);		   // Sine offset (pi/sb)
void STSI(int sb);		   // Strip size
void STSR(int sr);		   // Strip start
void SW_p(double swp);	   // Set spectral width (ppm)
void WDW(int wd);		   // Set window function
void WDW(const std::string& wd);   // Set window function

// ____________________________________________________________________________
// C                       XWinProcPar Input Functions
// ____________________________________________________________________________

/* These functions will read all parameters from a Bruker XWinNMR processing
   parameter file. The entire Bruker parameter file is simply read into a GAMMA
   parameter set (class XWinPSet) so that ALL parameters in the file are stored
   in this class. Were the Bruker files in GAMMA format then this would be done
   exclusively with class ParamSet, but since they are not we must handle the
   file parsing explicitly. This is relegated to the base class XWinPSet. In
   turn, specific parameters of consequnce to processing (e.g.processing 
   parameter files) are parsed within this class.

      Function                               Purpose
   ____________________________________________________________________________

      readPSet          Read in parameter set (class XWinPSet).  This is done
                        special since Bruker format is NOT in GAMMA parameter
                        format so it must be parsed appropriately.
     parsePSet          Parses a Bruker parameter set (XWinPSet) for parameters
                        of consequence to processing.
       read             Read in parameter set (class XWinPset) and parse out
                        parameters of consequence to an processing.

                       INHERITED FROM BASE CLASS XWinPSet

        bool readPSet(const string& filein, int warn=1);
        bool readPSet(int warn=1);                              */
	bool readPPar(const std::string& filein, int warn=2);
	bool readPPar(int warn=2);
virtual bool parsePSet(int warn=2);

// ____________________________________________________________________________
// D                       XWinProcPar Output Functions
// ____________________________________________________________________________
 
/* These function allow for output of NMR parameters directly into a Bruker
   XWinNMR ASCII parameter file {proc(s) & proc2(s)}.                        */
 
int writePPar(const std::string& name, int warn=2);
int writePPar(int warn=2) const;
 
// ____________________________________________________________________________
// E                 XWinProcPar Parameter Access Functions
// ____________________________________________________________________________
 
/* These functions allow access to the Bruker XWinNMR parameter set, not just
  those parameters which are directly handled within this class.
 
                       INHERITED FROM BASE CLASS XWinPSet
 
int getPar(const string& pn,int& val,   int id=0,int wrn=0) const;
int getPar(const string& pn,double& val,int id=0,int wrn=0) const;
int getPar(const string& pn,string& val,int id=0,int wrn=0) const;
ParameterSet getPSet() const;                                   */

// ____________________________________________________________________________
// F                   XWinProcPar Standard Output Functions
// ____________________________________________________________________________

/* These functions allow for a quick output of the parameter file contents.
   They don't have anything to do with output while running XWinNMR, rather
   users can just glance at procs parameters or store then in a small file.

ostream& printPset(ostream& ostr) const            INHERITED    */
std::ostream& print(std::ostream& ostr, int full=0, int hdr=1) const;
friend std::ostream& operator<< (std::ostream& ostr, const XWinProcPar& P);


// ____________________________________________________________________________
// F                      XWinProcPar Auxiliary Functions
// ____________________________________________________________________________

 /* These are just a group of functions that return strings indicating what
    the Bruker parameter values mean.  I just add things as I learn them here
    so that GAMMA output can remind me what all of these parameters are...   */

std::string AQORDERS()  const;
std::string BC_modS()   const;
std::string BYTORDPS()  const;
std::string DATMODS()   const;
std::string FT_modS()   const;
std::string INTBCS()    const;
std::string MC2S()      const;
std::string ME_modS()   const;
std::string PH_modS()   const;
std::string PKNLS()     const;
std::string PPARMODS()  const;
std::string PSCALS()    const;
std::string PSIGNS()    const;
std::string REVERSES()  const;
std::string SYMMS()     const;
std::string WDWS()      const;
};	

#endif 							// XWinProcPar.h

                                                                                



