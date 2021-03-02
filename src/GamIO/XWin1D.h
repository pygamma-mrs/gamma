/* XWin1D.h *****************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      XWin1D                                    Interface		**
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
** The XW* files provide an interface to Bruker XWinNMR (uxnmr) data    **
** sets. This class embodies a Bruker data set for a 1D acquisition.    **
** Each 1D acquisition generates several files (binary and ASCII) 	**
** spanning several directories.  This class is intended to handle 	**
** this structure.  Below is a typical Bruker directory structure 	**
** associated with a 1D acquisiton.					**
**									**
**                          __ acqu  (changable parameter file)		**
**			   / 						**
**                        /___ acqus (static parameter file)		**
**			 /						**
**  expname -- expnum --< ---- fid (binary data)			**
**			 \						**
**			  \___ pdata -- 1 -- proc, procs, meta		**
**									**
** This class will handle the directory hierarchy shown above as	**
** well as most of the files therein. Thus, given a base directory 	**
** for an experiment/simulation output (expname) & an experiment	**
** number (expnum) this class will generate the subdirectories shown,	**
** the acquisition parameter (acqu, acqus) and binary (fid) files	**
** as well as the processing files (proc, procs,meta).  It is NOT 	**
** intended to generated processed 1D data files (1i, 1r) - that 	**
** should be done using the XWinNMR software.  The class will also	**
** allow for import of 1D data sets.  The fid file and a processed	**
** spectrum can be read and placed into a GAMMA row_vector. Associated	**
** parameters are also provided for.					**
**									**
*************************************************************************/

#ifndef   XWin1D_H_ 				// Is file already included?
#  define XWin1D_H_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// This is the implementation
#  endif

#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <string>				// Include libstdc++ std::strings
#include <GamIO/XWinFid.h>			// Include XWinNMR fid file 
#include <GamIO/XWinAcqus.h>			// Include XWinNMR acqus file
#include <GamIO/XWinProcs.h>			// Include XWinNMR procs file
#include <GamIO/XWinSpec.h>			// Include XWinNMR spectra

class XWin1D: public XWinAcqus, public XWinFid,
              public XWinProcs, public XWinSpec

// ----------------------------------------------------------------------------
// ------------------------------- STRUCTURE ----------------------------------
// ----------------------------------------------------------------------------
 
  {
  std::string dname;					// Acquisition directory name
  int    Aidx;					// Acquisition experiment index
  int    Pidx;					// Acquisition processing index
  std::string NAIdir, NPdir, NPIdir;			// Data set directory names
  std::string Nacqu, Nacqus; 			// Acquisition ASCII file names
  std::string Nfid;					// Acquisition data file name
  std::string Nproc, Nprocs;				// Acquisition ASCII file names
  std::string Nmeta, Noutd;				// Acquisition ASCII file names
  std::string Nspec;					// Processed data base filename

 
// ----------------------------------------------------------------------------
// ----------------------------- PRIVATE FUNCTIONS ----------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________       
// i                   XWinNMR 1D Data Set Error Handling
// ____________________________________________________________________________
 
         void XWin1Derror(int eidx, int noret=1) const;
         void XWin1Derror(int eidx, const std::string& P, int nore=1) const;
volatile void XWin1Dfatality(int eidx) const;
volatile void XWin1Dfatality(int eidx, const std::string& P) const;
 
// ____________________________________________________________________________       
// ii                  XWinNMR 1D Data Set Error Handling
// ____________________________________________________________________________

int CheckDir(int TF,   int warn, const std::string& dout) const;
int CheckWrite(int TF, int warn, const std::string& dout) const;

// ____________________________________________________________________________
// iii                   XWinNMR 1D Data Set Setup Functions
// ____________________________________________________________________________
/* These functions are used to quickly set up values associated with the 1D
   data set.  In particular they handle setting up a standard Bruker directory
   structure and the file names in the data set.  Another function is used to
   insure the that the parameter sets are self-consistent. This should be
   called before the data set (or at least the parameter files) are output.  */

void SetNames();
int  MakeDirs(int warn=2);
void SetConsistent();

public:
// ____________________________________________________________________________ 
// A                XWinNMR 1D Data Set Constructors, Destructor
// ____________________________________________________________________________

/* These are the constructors of the class handling Bruker XWinNMR 1D data 
   sets. There are several files associated with a single data set FID. Those
   associated with the acquisition are the two ASCII parameter files (acqus and
   acqu) and the binary data file (fid).  In the processed data sub-directory,
   pdata will reside more files - whether or NOT the data has actually been
   processed.  These will be the parameter files (proc, procs, meta,....) and
   two binary files (1r, 1i) if the data has been processed.  GAMMA has 
   individual XWin* classes for handling these, so this class has the job of
   keeping track of directories and making sure all these are files are dealt
   with.
			bname   : Data set base directory
			mode    : I/O mode 
				  app    = append 
				  ate    = open & seek end (at end) 
				  binary = I/O in binary
				  in     = open for reading 
				  out    = open for writing 
				  trunc  = truncate file at 0 length        */
 
XWin1D();
virtual ~XWin1D();
 
// ____________________________________________________________________________
// B                    XWin1D Parameter Access Functions
// ____________________________________________________________________________

/* These functions allow direct access to some of the more important parameters
   that each 1D data set should know.  Since this class is derived from many
   classes, each of which handle different files, we simply inherit much of our
   functionality.  When parameters must be consistent across the files, we take
   the time to insure that this is indeed the case, not here, but when we write
   out the data set....

			      Acquisiton Parameters 
			      (Handled by XWinAcqus)

double acqname()    const;	// The acqus file name (short)
double AQ()         const;	// Acquisition length (sec)
int    AQ_mod()     const;	// Acquisition mode
double BF1()        const;	// Base frequency channel 1
double BF2()        const;	// Base frequency channel 2
int    BYTORDA()    const; 	// Binary byte ordering
int    DS()         const;	// Number of dummy scans
double IN(int i)    const;
std::string EXP()        const;	// Experiment name
std::string NAME()       const;	// Full File Name (long)
int    NS()         const;	// Number of scans
std::string NUC(int i)   const;	// Nucleus for a channel
std::string NUCLEUS()    const;	// Base nucleus
double O1()         const;	// Offset frequency channel 1
double O2()         const;	// Offset frequency channel 2
int    PARMODE()    const;	// Acquisiiton dimension
std::string PULPROG()    const;	// Pulse program
double SFO1()       const;	// Spectrometer freq.
double SFO2()       const;	// Spectrometer freq.
double SFO3()       const;	// Spectrometer freq.
std::string SOLVENT()    const;	// Solvent
double SW()         const;	// Spectral width (PPM)
double SW_h()       const;	// Spectral width (Hz)
int    TD()         const;	// Total points
double TE()         const;	// Sample temperature

			          Binary Access
			      (Handled by XWinFid)

std::string     XWinFid::fidname()   const;          // Fid File name (short)
int        XWinFid::size()      const;          // No. complex points
int        XWinFid::TDF()       const;          // No. total points (re+im)
bool       XWinFid::order()     const;          // Data byte order
int        XWinFid::bytes()     const;          // File size in bytes
int        XWinFid::blocks()    const;          // No. FIDs
int        XWinFid::pad()       const;          // No. padding bytes
row_vector XWinFid::data()      const;          // Data points

			          Binary Access
			      (Handled by XWinSpec)

std::string     XWinSpec::specname()  const;         // Base file name       
std::string     XWinSpec::specrname() const;         // Reals file name
std::string     XWinSpec::speciname() const;         // Imagins. file name
int        XWinSpec::size()      const;         // No. complex points
int        XWinSpec::SSI()       const;         // No. total points
bool       XWinSpec::order()     const;         // Data byte order
int        XWinSpec::bytes()     const;         // File size in bytes
int        XWinSpec::blocks()    const;         // No. FIDs
int        XWinSpec::pad()       const;         // No. padding bytes
row_vector XWinSpec::data()      const;         // Data points

			      Processing Parameters 
			      (Handled by XWinProcs)

std::string XWinProcs::parname()   const;  // ASCII File name
int    XWinProcs::BYTORDP()   const;  // Binary byte order
int    XWinProcs::FT_mod()    const;  // How FFT is performed
double XWinProcs::LB()        const;  // Line Broadening
int    XWinProcs::MC2()       const;  // FT type on t1
double XWinProcs::OFFSET()    const;  // Spectrum offset
double XWinProcs::PHC0()      const;  // Zero order phase
double XWinProcs::PHC1()      const;  // 1st order phase
std::string XWinProcs::REVERSE()   const;  // Plot spectrum reverse
double XWinProcs::SF()        const;  // Spectrometer frequency
int    XWinProcs::SI()        const;  // Data size (re+im)
int    XWinProcs::SSB()       const;  // Sine bell
int    XWinProcs::STSI()      const;  // Strip size
int    XWinProcs::STSR()      const;  // Strip start
double XWinProcs::SW_p()      const;  // Spectral width (PPM)
double XWinProcs::TDeff()     const;  // Effective FFT size
int    XWinProcs::WDW()       const;  // Window function

			      Acquisiton Parameters 
			      (Handled by XWinAcqus)

void AQ_mod(int aqmo);
int  D(int idx, double tsec, int warn=2);
void DS(int ds);
void EXP(const std::string& exp);
void NAME(const std::string& name);
void NS(int ns);
void NUC(int i, const std::string& N);
//void NUCLEI(int channel, const std::string& I, double O, int warn=2);
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
void TD(int npts);
void TE(double te);

                              Processing Parameters
                              (Handled by XWinProcs)

void MC2(int mc);
void REVERSE(int yn);
void PPARMOD(int pm);
void SI(int si);
void SF(double SF);
void SSB(int sb);
void STSI(int sb);
void STSR(int sr);
void SW_p(double swp);
void WDW(int wd); 						    */

std::string     name()     const;
row_vector FID()      const;
row_vector Spectrum() const;

// ____________________________________________________________________________
// C                   XWinNMR 1D Data Set Input Functions
// ____________________________________________________________________________

/* These functions will read in the acquisition parameters from an XWinNMR
   parameter file, typically named acqus.  Any parameters relevant to GAMMA
   will be internally stored and accessible.                                 */

int read(const std::string& dirin, int aidx=1, int pidx=1, int warn=2);
int read(int warn=2);
int readFID(const std::string& dirin, int exp=1, int warn=2);
int readFID(int warn=2);
int readSpectrum(const std::string& dirin, int aidx=1, int pidx=1, int w=2);
int readSpectrum(int warn=2);

// ____________________________________________________________________________
// C                   XWinNMR 1D Data Set Input/Output Functions
// ____________________________________________________________________________

/* These functions will write a single 1D Data Set out in Bruker XWinNMR 
   format. This can be done by directly writing a row vector after first
   specifying any acquisition parameters.  Note that in using this function
   the entire XWinNMR 1D directory structure will be output:

                           __ acqu  (changable parameter file)
			  / 
                         /___ acqus (static parameter file)
			/
   expname -- expnum --< ---- fid (binary data)
			\
			 \--- pdata -- specnum -- proc, procs, meta

   Here expname is that specified output.  The numbers expnum and specnum
   will default to 1.  The produced files are all ASCII except fid. The
   files acqu and acqus will be equivalent as will proc and procs.  All
   files in pdata, since this class will NOT do any OUTPUT of frequency
   domain 1D data, will be defaults based on values in acqus.                */
 
int write(const std::string& fout, const row_vector& data, int warn=2); 
int write(const row_vector& data, int warn=2); 

// ____________________________________________________________________________ 
// E                 XWinNMR 1D Data Set ASCII Output Functions
// ____________________________________________________________________________

/* These functions allow users to get a look at what the data set contains.
   They are set to write all sorts of information concerning both the ASCII
   parameter file (acqus) and the binary data file (fid) in a nice format.   */

virtual std::ostream& print(std::ostream& out, int full=1) const;
friend  std::ostream& operator<<(std::ostream& out, const XWin1D& XW1D);
 
        // Input                out      : output stream;
        //                      XW1D     : XWinNMR 1D Data Set
        // Output               none     : Modifies output stream


// ____________________________________________________________________________
// F              XWinNMR 1D Data Set Interactive Functions
// ____________________________________________________________________________


std::string ask_read(int argc,char* argv[],int& argn,int aidx=1,int pidx=1);

	// Input                XW1D : An XWinNMR 1D data set
	//                      argc    : Number of arguments
	//                      argv    : Vector of argc arguments
	//                      argn    : Argument index
        //                      aidx    : Flag to request experiment #
        //                                 aidx < 1: Ask for it
        //                                 aidx >=1: Use aidx itself def (1)
        //                      pidx    : Flag to request processing #
        //                                 pidx < 1: Ask for it
        //                                 pidx >=1: Use pidx itself def (1)
	// Output               dirname : The parameter argn of array argc
	//                                is used to supply a directory name
	//                                from which the 1d data set resides
	//                                If the argument argn is not in argv,
	//                                the user is asked to supply a name
	//				  and experiment number
	//                                The directory is returned
	// Note                         : The directory should contain an
	//				  XWinNMR 1D data set.
};

#endif 							// XWin1D.h
