/* XWinAcqus.h **************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      XWinAcqus                                   Interface		**
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
** trivial. Note that this data file has nothing to do with GAMMA and	**
** will primarily be used for output.  It is usually found as		**
** expdir/1/acqus where expdir is a directory where an experiment is	**
** writing an acquisition.  Below is a typical Bruker directory		**
** structure in a 1D acquisiton case.					**
**									**
**                          __ acqu  (changable parameter file)		**
**			   / 						**
**                        /___ acqus (static parameter file)		**
**			 /						**
**  expname -- expnum --< ---- fid (binary data)			**
**			 \						**
**			  \___ pdata -- 1 -- proc, procs, meta		**
**									**
*************************************************************************/

#ifndef   XWinAcqus_H_ 			// Is file already included?
#  define XWinAcqus_H_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the implementation
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <string>			// Include libstdc++ strings
#include <Basics/ParamSet.h>		// Include GAMMA parameter sets
#include <GamIO/XWinPSet.h>             // Include Bruker parameter parsing
#include <GamIO/XWinAcqPar.h>           // Include Bruker acquisition pars

class XWinAcqus: public XWinAcqPar 

// ----------------------------------------------------------------------------
// ------------------------------- STRUCTURE ----------------------------------
// ----------------------------------------------------------------------------
 
  {
 
// ----------------------------------------------------------------------------
// ----------------------------- PRIVATE FUNCTIONS ----------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                      XWinNMR Acqus File Error Handling
// ____________________________________________________________________________

/* These functions take care of any errors encountered when reading, writing,   
   and setting parameters in Bruker acquisition parameter files.
 
	Input		AcqPar  : XWinNMR acqusition parameters (this)
			eidx    : Error index
                        pname   : Additional error message
                        noret   : Flag for linefeed (0=linefeed)
        Output          void    : An error message is output                 */

         void XWinAcquserror(int    eidx, int noret=0) const;
         void XWinAcquserror(int    eidx, const std::string& pname, int noret=0) const;
volatile void XWinAcqusfatality(int eidx) const;
volatile void XWinAcqusfatality(int eidx, const std::string& pname) const;

public:
// ____________________________________________________________________________ 
// A             XWinAcqus Parameter File Constructors, Destructor
// ____________________________________________________________________________
 
/* These are the constructors of the class handling Bruker XWinNMR acquisition
   parameter files.  This doesn't do anything in particular, it is the read
   and write functions that perform the work.  Thus, we only have a default
   constructor specified.  The reading and writing of the associated ASCII
   parameter file is done in one step, so we don't need anything complex.    */

        XWinAcqus();
        XWinAcqus(const std::string& name);
        XWinAcqus(const XWinAcqus& XWA);
virtual ~XWinAcqus();
void    operator= (const XWinAcqus& XWA);

// ____________________________________________________________________________
// B                  XWinAcqus Parameter Access Functions
// ____________________________________________________________________________
 
/* These functions allow direct access to some of the more important parameters
   that each acqusition file should know.  The two primary values are the fid
   point size (TD) and the data byte order (BYTORDA).

                     INHERITED FROM CLASS XWinAcqPar
 
int    acqname()  const;	// File name (short)
double field()    const;	// Field strength (T)
double AQ()	     const;	// Acquisition length (sec)
int    AQ_mod()   const;	// Acquisition mode
double BF1()      const;	// Base Spectrometer freq.
double BF2()      const;	// Base Spectrometer freq.
int    BYTORDA()  const;	// Binary byte order
int    DS()       const;	// Number of dummy scans
string EXP()      const;	// Experiment name
double IN(int i)  const;	// Dwell time
string NAME()     const;     // Full File Name
int    NS()       const;	// Number of scans
string NUC(int i) const;	// Nucleus for a channel
string NUCLEUS()  const;	// Base nucleus
double O1()       const;	// Offset freq.
double O2()       const;	// Offset freq.
int    PARMODE()  const;	// Acquisiiton dimension
string PULPROG()  const;	// Pulse program
double SFO1()     const;	// Spectrometer freq.
double SFO2()     const;	// Spectrometer freq.
double SFO3()     const;	// Spectrometer freq.
string SOLVENT()  const;	// Solvent
double SW()       const;	// Spectral width (PPM)
double SW_h()     const;	// Spectral width (Hz)
int    TD()       const;	// Total points
double TE()       const;	// Sample temperature

void AQ_mod(int aqmo);
void BF1(double bf);
void BF2(double bf);
void BYTORDA(int bo);
int  D(int idx, double tsec, int warn=2);
void DS(int ds);
void EXP(const string& exp);
void IN(int i, double in);
void O1(double of);
void O2(double of);
void NAME(const string& name);
void NS(int ns);
void NUC(int i, const string& N);
void NUCLEI(int channel, const string& I, double O, int warn=2);
void NUCLEUS(const string& I);
int  P(int idx, double tp, int warn=2);
void PARMODE(int pm);
void PULPROG(const string& P);
void SFO1(double sf);
void SFO2(double sf);
void SFO3(double sf);
void SFO(double sf, int i);
void SOLVENT(const string& S);
void SW(double sw);
void SW_h(double sw);
void TE(double te);
void TD(int npts);                                                */
 
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

                     INHERITED FROM CLASS XWinAcqPar
 
virtual bool readAPar(const string& filein, int warn=1);
virtual bool readPar(int warn=1);
vittual bool parsePSet(int warn=1);                               */

 
// ____________________________________________________________________________
// D                       XWinAcqus Output Functions
// ____________________________________________________________________________

/* These function allow for output of NMR parameters directly into a Bruker
   XWinNMR ASCII parameter file (acqus).

                       INHERITED FROM CLASS XWinAcqPar

int writeAPar(const string& name, int warn=2);
int writeAPar(int warn=2) const;                                  */

// ____________________________________________________________________________
// E                 XWinAcqus Standard Output Functions
// ____________________________________________________________________________

/* These function allow for a quick output of the parameter file contents.
   They don't have anything to do with output while running XWinNMR, rather
   users can just glance at acqus parameters or store then in a small file.

        ostream& printPset(ostream& ostr) const;	INHERITED            */
virtual std::ostream& print(std::ostream& ostr, int full=0, int hdr=1) const;
friend  std::ostream& operator<< (std::ostream& ostr, const XWinAcqus& AcqPar);

// ____________________________________________________________________________
// Z                    XWinAcqus BrukerLike Output Functions
// ____________________________________________________________________________

/* These functions perform ASCII output of acquisiton parameters.  The format
   here is similar to issuing a "dpa" command in XWinNMR or using the menu
   choice "Output/Display status pars./Acquisition only".  It isn't as complete
   as the Bruker output because I'm bored and its often a mystery how they
   get some of their values.....                                             */

std::ostream& dpa(std::ostream& ostr) const;
virtual std::ostream& printGAMMA(std::ostream& ostr) const;
};

#endif 								// XWinAcqus.h
