/* XWinProc2s.h ***************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      XWinProc2s                                   Interface		**
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
** sets. This class embodies a Bruker parameter file, proc2s, which     **
** contains parameters that control NMR processing. This		**
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
**			  \___ pdata -- 1 -- proc, procs, meta		**
**									**
*************************************************************************/

#ifndef   XWinProc2s_H_ 		// Is file already included?
#  define XWinProc2s_H_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the implementation
#  endif
 
#include <string>			// Include libstdc++ strings
#include <Basics/ParamSet.h>		// Include GAMMA parameter sets
#include <GamIO/XWinPSet.h>		// Include Bruker parameter parsing
#include <GamIO/XWinProcPar.h>		// Include Bruker prossing pars

class XWinProc2s: public XWinProcPar

// ----------------------------------------------------------------------------
// ------------------------------- STRUCTURE ----------------------------------
// ----------------------------------------------------------------------------
 
  {
 
// ----------------------------------------------------------------------------
// ----------------------------- PRIVATE FUNCTIONS ----------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________       
// i                      XWinNMR Procs File Error Handling
// ____________________________________________________________________________

/* These functions take care of any errors encountered when reading, writing,   
   and setting parameters in Bruker procs parameter files.
 
	Input		ProcPar : XWinNMR procs parameters (this)
			eidx    : Error index
                        pname   : Additional error message
                        noret   : Flag for linefeed (0=linefeed)
        Output          void    : An error message is output                 */

         void XWinPP2error(int eidx, int noret=0) const;
         void XWinPP2error(int eidx, const std::string& pname, int noret=0) const;
volatile void XWinPP2fatal(int eidx) const;
volatile void XWinPP2fatal(int eidx, const std::string& pname) const;

public:
// ____________________________________________________________________________ 
// A             XWinProc2s Parameter File Constructors, Destructor
// ____________________________________________________________________________
 
/* These are the constructors of the class handling Bruker XWinNMR output
   parameter files.  This doesn't do anything in particular, it is the write
   functions that perform the work.                                          */

        XWinProc2s();
        XWinProc2s(const std::string& name);
        XWinProc2s(const XWinProc2s& XWP);
virtual ~XWinProc2s();
void    operator= (const XWinProc2s& XWP);

// ____________________________________________________________________________
// B                  XWinProc2s Parameter Access Functions
// ____________________________________________________________________________
 
/* These functions allow direct access to some of the more important parameters
  that each processing file should know. The two primary values are the 
  spectrum point size (SI) and the data byte order (BYTORDP).

		     INHERITED FROM CLASS XWinProcPar

string parname() const;		// ASCII File name
int    BYTORDP() const;		// Byte order
int    FT_mod()  const;		// FFT mode
double LB()      const;		// Line broadening
int    MC2()     const;		// Sine bell?
double OFFSET()  const;		// Spectrum offset
double PHC0()    const;		// Zero order phase
double PHC1()    const;		// 1st order phase
string REVERSE() const;		// Spectrum reverse
double SF()      const;		// Spectrometer freq.
int    SI()      const;		// Data size (re+im)
int    SSB()     const;		// Sine bell?
int    STSI()    const;		// Sine bell?
int    STSR()    const;		// Sine bell?
double SW_p()    const;		// Spectral Width (ppm)
double TDeff()   const;		// Effective FID size
int    WDW()     const;		// Window function

void BYTORDP(bo);
void MC2(int mc);
void REVERSE(int yn);
void PPARMOD(int pm);
void SI(int si);
void SF(double SF);
void SSB(int sb);
void STSI(int sb);
void STSR(int sr);
void SW_p(double swp);
void WDW(int wd);                                               */

// ____________________________________________________________________________
// C                       XWinProc2s Input Functions
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
     parsePSet          Parses a Bruker parameter set (XWinPSet) for parameters
                        of consequence to processing.

        parse read(const string& filein, int warn=2);
virtual parse read(int warn=1);
virtual parse parsePSet(int warn=1);                             */


// ____________________________________________________________________________
// D                      XWinProc2s Output Functions
// ____________________________________________________________________________

/* These function allow for output of NMR parameters directly into a Bruker
   XWinNMR ASCII parameter file (procs).
 
virtual int      write(const string& name, int warn=2);
virtual int      write(int warn=2) const;                                    */

// ____________________________________________________________________________
// E                    XWinProc2s Standard Output Functions
// ____________________________________________________________________________

/* These function allow for a quick output of the parameter file contents.
   They don't have anything to do with output while running XWinNMR, rather
   users can just glance at procs parameters or store then in a small file.  */


//ostream& printPset(ostream& O) const;		INHERITED
//virtual ostream& print(ostream& ostr, int full=0, int hdr=1) const;
//friend ostream& operator<< (ostream& O, const XWinProc2s& P);


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

std::ostream& dpp(std::ostream& ostr, double BF0=0) const;


};

#endif 								// XWinProc2s.h

                                                                                



