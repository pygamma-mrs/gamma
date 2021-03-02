/* XWin2D.h *****************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      XWin2D                                    Interface		**
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
** sets. This class embodies a Bruker data set for a 2D acquisition.    **
** Each 2D acquisition generates several files (binary and ASCII) 	**
** spanning several directories.  This class is intended to handle 	**
** this structure.  Below is a typical Bruker directory structure 	**
** associated with a 2D acquisiton.					**
**									**
**                          __ acqu, acqu2 (changable parameter file)	**
**			   / 						**
**                        /___ acqus, acqu2s (static parameter file)	**
**			 /						**
**  expname -- expnum --< ---- ser (binary data)			**
**			 \						**
**			  \___ pdata -- 1 -- proc, proc2, procs, proc2s **
**									**
** This class will handle the directory hierarchy shown above as	**
** well as most of the files therein. Thus, given a base directory 	**
** for an experiment/simulation output (expname) & an experiment	**
** number (expnum) this class will generate the subdirectories shown,	**
** the acquisition parameter (acqu*) and binary (ser) files as well as	**
** the processing files (proc*, meta).  It is NOT intended to generated	**
** processed 2D/ND data files - that should be done using the XWinNMR 	**
** software.  The class will also allow for import of 2D data sets.	**
** The serial file (ser), & to a more limited extent processed spectra,	**
** can be read and placed into GAMMA row_vectors and matrices. The	**
** associated parameters are also provided for.				**
**									**
*************************************************************************/

#ifndef   XWin2D_H_ 				// Is file already included?
#  define XWin2D_H_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// This is the implementation
#  endif

#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <string>				// Include libstdc++ std::strings
#include <GamIO/XWinSer.h>			// Include XWinNMR serial file 
#include <GamIO/XWinAcqus.h>			// Include XWinNMR acqus file
#include <GamIO/XWinAcqu2s.h>			// Include XWinNMR acqu2s file
#include <GamIO/XWinProcs.h>			// Include XWinNMR procs file
#include <GamIO/XWinProc2s.h>			// Include XWinNMR proc2s file
#include <GamIO/XWinMeta.h>			// Include XWinNMR meta file
#include <GamIO/XWinOutd.h>			// Include XWinNMR outd file

class XWin2D

// ----------------------------------------------------------------------------
// ------------------------------- STRUCTURE ----------------------------------
// ----------------------------------------------------------------------------
 
  {
  int         oldMeta;				// Flag older meta format
  std::string dname;				// Acquisition directory name
  int         Aidx;				// Acquisition experiment index
  int         Pidx;				// Acquisition processing index
  std::string NAIdir, NPdir, NPIdir;		// Data set directory names
  std::string Nacqu, Nacqus, Nacqu2, Nacqu2s;	// Acquisition ASCII file names
  std::string Nser;				// Acquisition data file name
  std::string Nproc, Nprocs, Nproc2, Nproc2s;	// Processing ASCII file names
  std::string Nmeta, Noutd;			// Processing ASCII file names
  XWinAcqus   Acqs;				// Acq.  Parameter file for 1D
  XWinProcs   Procs;				// Proc. Parameter file for 1D
  XWinSer     Ser;				// Binary serial data file 
  XWinAcqu2s  Acq2s;				// Acq.  Parameter file for 2D
  XWinProc2s  Proc2s;				// Proc. Parameter file for 2D
//  double     FieldT;				// Bo Field in Teslas
 
// ----------------------------------------------------------------------------
// ----------------------------- PRIVATE FUNCTIONS ----------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________       
// i                   XWinNMR 2D Data Set Error Handling
// ____________________________________________________________________________
 
         void XWin2Derror(int eidx, int noret=1) const;
         void XWin2Derror(int eidx, const std::string& P, int noret=1) const;
volatile void XWin2Dfatality(int eidx) const;
volatile void XWin2Dfatality(int eidx, const std::string& P) const;
 
// ____________________________________________________________________________       
// ii                  XWinNMR 2D Data Set Error Handling
// ____________________________________________________________________________

int CheckDir(int TF,   int warn, const std::string& dout) const;
int CheckWrite(int TF, int warn, const std::string& dout) const;

// ____________________________________________________________________________ 
// iii                   XWinNMR 2D Data Set Setup Functions
// ____________________________________________________________________________
/* These functions are used to quickly set up values associated with the 2D
   data set.  In particular they handle setting up a standard Bruker directory
   structure and the file names in the data set.  Another function is used to
   insure the that the parameter sets are self-consistent. This should be
   called before the data set (or at least the parameter files) are output.  */

//void SetField();			// Set Spectrometer Field (T)
//void FieldReset(double BoT);		// Set Spectrometer Field (T)
void SetNames();			// Set Directory & File Names
int  MakeDirs(int warn=2);		// Construct Directories
int  ReadPars(int warn=2);		// Read All Parameter Files
void SetConsistent();
 
public:
// ____________________________________________________________________________ 
// A                XWinNMR 2D Data Set Constructors, Destructor
// ____________________________________________________________________________

/* These are the constructors of the class handling Bruker XWinNMR 2D data
   sets. There are several files associated with such a data set. Those
   associated with the acquisition are the four ASCII parameter files
  (acqu, acqus, acqu2, and acqu2s) and the binary serial data file (ser).
   In the processed data sub-directory, pdata will reside more files - whether
   or NOT the data has actually been processed.  These will be the six
   parameter files (proc, procs, proc2, proc2s, meta, and outd).  If the data
   has been processed there may also exist more binary files (e.g. 2rr) if the
   data has been processed.  GAMMA has individual XWin* classes for handling
   these individual files, so this class has the job of keeping track of
   directories and making sure all these are appropriate files are dealt with.

   Note that these constructors do NO writing! Furthermore, they will only read
   the ASCII parameter files and NOT the larger binary files.  Writing and/or
   binary reading is left to other functions.                                */

XWin2D();
XWin2D(const std::string& dname, int mode = std::ios::in, int eno=1, int pno=0);
XWin2D(const XWin2D& XW2D);
virtual ~XWin2D();
virtual XWin2D& operator=(const XWin2D& XW2D);

 
// ____________________________________________________________________________
// B                    XWin2D Parameter Access Functions
// ____________________________________________________________________________

/* These functions allow direct access to some of the more important parameters
   that each 2D data set should know.  This class contains objects of many
   types associated with the Bruker 2D data set files. Since they are quite
   similar we do not derive XWin2D from these, since this might cause conflicts
   (that may or may not be sorted out by virtual tables). Rather we rewrite the
   functions herein to match the lower classes and route to the proper
   dimension when necessary.....                                             */

/* --------------------------------------------------------------------------
                          Functions To Retrieve Parameters
   --------------------------------------------------------------------------*/
 
/*                   NOT REALLY INHERITED FROM CLASS XWinAcqPar
                           (Applies To Acqs and Acq2s)                       */
 
std::string acqname(int d=0)    const;	// Acq. Par. File name (short)
double AQ(int d=0)         const;	// Acquisition length (sec)
int    AQ_mod(int d=0)     const;	// Acquisition mode
double BF1(int d=0)        const;	// Base Spectrometer freq.
double BF2(int d=0)        const;	// Base Spectrometer freq.
int    BYTORDA(int d=0)    const;	// Binary byte order
int    DSc(int d=0)         const;	// Number of dummy scans
std::string EXP(int d=0)        const;	// Experiment name
double XW_IN(int i,int d=0)   const;	// Dwell time
std::string NAME(int d=0)       const;	// Full File Name
int    NS(int d=0)         const;	// Number of scans
std::string NUC(int i,int d=0)  const;	// Nucleus for a channel
std::string NUCLEUS(int d=0)    const;	// Base nucleus
double O1(int d=0)         const;	// Offset freq.
double O2(int d=0)         const;	// Offset freq.
int    PARMODE(int d=0)    const;	// Acquisiiton dimension
std::string PULPROG(int d=0)    const;	// Pulse program
double SFO1(int d=0)       const;	// Spectrometer freq.
double SFO2(int d=0)       const;	// Spectrometer freq.
double SFO3(int d=0)       const;	// Spectrometer freq.
std::string SOLVENT(int d=0)    const;	// Solvent
double SW(int d=0)         const;	// Spectral width (PPM)
double SW_h(int d=0)       const;	// Spectral width (Hz)
int    TD(int d=0)         const;	// Total points
double TE(int d=0)         const;	// Sample temperature

/*                   NOT REALLY INHERITED FROM CLASS XWinSer
                             (Applies to File ser)                           */
 
std::string sername() 	const;		// File name
int    TDS()		const;		// No. total points

/*                  NOT REALLY INHERITED FROM CLASS XWinProcPar
                           (Applies To Procs and Proc2s)                    */
 
std::string parname(int d=0)   const;  // ASCII File name
int    BYTORDP(int d=0)   const;  // Binary byte order
int    FT_mod(int d=0)    const;  // How FFT is performed
double LB(int d=0)        const;  // Line Broadening
int    MC2(int d=0)       const;  // FT type on t1
double OFFSET(int d=0)    const;  // Spectrum offset
double PHC0(int d=0)      const;  // Zero order phase
double PHC1(int d=0)      const;  // 1st order phase
int    PH_mod(int d=0)    const;  // Phasing type
int    REVERSE(int d=0)   const;  // Plot spectrum reverse
double SF(int d=0)        const;  // Spectrometer frequency
int    SI(int d=0)        const;  // Data size (re+im)
int    SSB(int d=0)       const;  // Sine bell
int    STSI(int d=0)      const;  // Strip size
int    STSR(int d=0)      const;  // Strip start
double SW_p(int d=0)      const;  // Spectral width (PPM)
double TDeff(int d=0)     const;  // Effective FFT size
int    WDW(int d=0)       const;  // Window function
 
/*                   ACCESS FUNCTIONS ON THIS LEVEL ONLY                     */

std::string name()  const;
double Field() const;


/* -------------------------------------------------------------------------- 
                          Functions To Specify Parameters
   --------------------------------------------------------------------------*/

/*                NOT REALLY INHERITED FROM CLASS XWinAcqPar
                         (Applies To Acqs and Acq2s)                         */
 
void AQ_mod(int aqmo, int d=0);		// Acquisiiton mode
//void BF1(double bf, int d=0);		// 1st/2nd channel Omega
//void BF2(double bf, int d=0);		// 2nd/1st channel Omega
//void BYTORDA(int bo, int d=0);	// Bin. byte order <=== arch
int  D(int idx, double tsec, int d=0, int warn=2);
void DSc(int ds, int d=0);		// Dummy scans
void EXP(const std::string& exp, int d=0);	// Experiment name
void XW_IN(int i, double in, int d=0);	// Delay increments
void O1(double of, int d=0);		// 1st/2nd channel offset
void O2(double of, int d=0);		// 2nd/1st channel offset
//void NAME(const std::string& n, int d=0);	// File name <===== on output
void NS(int ns, int d=0);		// Number of scans
//void NUC(int i, const std::string& N, int d=0);
void NUCLEI(int channel, const std::string& I, double O, int warn=2);
//void NUCLEUS(const std::string& I, int d=0);
int  P(int idx, double tp, int d=0, int warn=2);
//void PARMODE(int pm, int d=0);	// Data dimension <======== 2
void PULPROG(const std::string& P, int d=0);	// Pulse program
void SFO1(double sf, int d=0);		// 1st/2nd channel spec. freq.
void SFO2(double sf, int d=0);		// 2nd/1st channel spec. freq.
//void SFO3(double sf, int d=0);	// 3rd channel spec. freq.
//void SFO(double sf, int i, int d=0);	// Spectrometer frequencies
void SOLVENT(const std::string& S, int d=0);	// Solvent type
void SW(double sw, int d=-1);		// Spectral width in ppm 
void SW_h(double sw, int d=-1);		// Spectral width in Hz
//void TD(int npts, int d=0);		// Set size <======= on output
void TE(double te, int d=0);		// Set temperature
 
/*                  NOT REALLY INHERITED FROM CLASS XWinProcPar
                           (Applies To Procs and Proc2s)                    */

//void BYTORDP(int bo, int d=0);	// Bin. byte order <===== arch
void FT_mod(int ft, int d=-1);		// Set transform type 
void FT_mod(const std::string& ft, int d=0);	// Set transform type
void GB(int gb,  int d=-1);		// Set Gaussian broadening
void LB(int lb,  int d=-1);		// Set line broadening
void MC2(int mc, int d=0);		// Set acquisition type (t1)
void MC2(const std::string& mc, int d=0);	// Set acquisition type (t1)
void PHC0(double ph0,  int d=-1);	// Set 0th order phase correct
void PHC1(double ph1,  int d=-1);	// Set 1st order phase correct
void PH_mod(const std::string& p, int d=-1);	// Set type of phase correct
void REVERSE(int yn, int d=-1);		// Set spectrum reverse
//void PPARMOD(int pm, int d=0);	// Data dim. <============== 2
//void SI(int si, int d=0);		// Set data size (re+im)
//void SF(double SF, int d=0);		// Set spectrometer freq.
//void SR(double SR, int d=0);		// Ref. freq. <========= BF-SF
void SSB(int sb, int d=0);		// Sine offset (pi/sb)
//void STSI(int sb, int d=0);		// Strip size
//void STSR(int sr, int d=0);		// Strip start
//void SW_p(double swp, int d=0);	// SW (ppm) <========= acqus SW
void WDW(int wd, int d=-1);		// Set window function
void WDW(const std::string& wd, int d=-1);	// Set window function

/*                   ACCESS FUNCTIONS ON THIS LEVEL ONLY                     */  
 
void Field(double Bo);
void Field(double Om, const std::string& I);
void OldMeta(int om);

 

// ____________________________________________________________________________
// C                 XWinNMR 2D Data Set Input Functions
// ____________________________________________________________________________

/* These functions will read part, or the entire, Bruker XWinNMR 2D Data Set. 

                           __ acqu, acqu2 (changable parameter files)
			  / 
                         /___ acqus, acqu2s (static parameter files)
			/
   expname -- expnum --< ---- ser (binary data)
			\
			 \--- pdata -- 1 -- proc, proc2, procs, proc2s, meta

   Here expname is the base directory containing the data set. The number
   expnum is the experiment number an will default to 1. The acquisition
   parameters are taken from the file acqus and the acquisition binary
   points from fid.  The processed data, if accessed, will come from the
   parameter files procs, meta, etc. and binary data files 1r and 1i.

                       INHERITED FROM CLASS XWinSer

   Function   Arguments              Action                       Return
   --------  ------------  -----------------------------------  ----------
   readFID   fn,TD,bytord  Read 1st block in serial/FID file    row_vector
   readFID   fn,TD,bord,i  Read ith block in serial/FID file    row_vector
 + readFID   idx           Read block of index idx(-1=current)  row_vector
   readFIDs  fn,TD,bytord  Read entire serial file                matrix
 + readFIDs                Read entire serial file                matrix
		    
	     Above +=File Must Be Properly Opened For Reading

row_vector XWinSer::readFID(const std::string& fin, int TD, int bytord, int idx=-1)
row_vector XWinSer::readFID(int idx=-1)
matrix     XWinSer::readFIDs(const std::string& fin, int TD, int bytord)
matrix     XWinSer::readFIDs()                                               */

row_vector readFID(const std::string& d,int I=-1,int a=1,int p=1,int wn=2);
row_vector readFID(int idx, int warn=2);

// ____________________________________________________________________________
// D                 XWinNMR 2D Data Set Output Functions
// ____________________________________________________________________________

/* These functions will write an entire Bruker XWinNMR 2D Data Set. 

                           __ acqu, acqu2  (changable parameter files)
			  / 
                         /___ acqus, acqu2s (static parameter file)
			/
     expname -- Aidx --< ---- ser (binary data)
			\
			 \--- pdata -- Pidx -- proc, procs, proc2, proc2s

   Here expname is the base directory specified for output. The numbers Aidx
   and Pidx default to 1.  The produced files are all ASCII except ser. The
   files acqu and acqus will be equivalent as will acqu2-acqu2s, proc-procs,
   and proc2-proc2s.  All files in pdata, since this class will NOT do any
   OUTPUT of frequency domain 2D data, will be defaults based on values in
   acqus and acqu2s. 

   Function  Arguments                         Result
   ________  _________  _______________________________________________________
    write    dir,mx,i     XWinNMR 2D data set produced in dir with data mx 

   Note that GAMMA handles only some of the ASCII parameters one might want
   set in these data sets (# dimensions, # blocks, block size,...). Thus users
   should set additional parameters prior using these functions if values other
   than defaults are desired.  Suggested parameters one might set are:

   Parmeter File  Parameter                    Settings
   _____________  _________  __________________________________________________
    acqu,acqus     AQ_mod    0=qf=one channel, 1=qsim=quadrature simultaneous
			     2=qseq=quadrature sequential, 3=digital quadrature
                   BF1,BF2   Base frequencies in MHz
                   O1,O2     Offset frequencies in Hz
                  SFO{1,2,3} Spectrometer frequencies in MHz
		  SW, SW_h   Spectral Width (PPM, Hz)
									     */

int write(const std::string& Bdir, const matrix& data, int warn=2);
//int write(const row_vector& data, int warn=2);


// ____________________________________________________________________________ 
// E                 XWinNMR 2D Data Set ASCII Output Functions
// ____________________________________________________________________________

/* These functions allow for a quick output of the data set contents. They
   don't have anything to do with output while running XWinNMR, rather users
   may just glance at acquisition parameters & values or store then in a 
   small file if needed.  The two function "dpa" and "dpp" mimick the output
   that would be seen in XWinNMR.                                            */

std::ostream& dpa(std::ostream& ostr, const std::string& dirin);
std::ostream& dpa(std::ostream& ostr) const;
std::ostream& dpp(std::ostream& ostr, const std::string& dirname);
std::ostream& dpp(std::ostream& ostr) const;

       std::ostream& print(std::ostream& out, int full=1) const;
friend std::ostream& operator<<(std::ostream& out, const XWin2D& XW2D);
 
        // Input                out      : output stream;
        //                      XW2D     : XWinNMR 2D Data Set
        // Output               none     : Modifies output stream


// ____________________________________________________________________________
// F              XWinNMR 2D Data Set Interactive Functions
// ____________________________________________________________________________


//std::string ask_read(int argc, char* argv[], int& argn, int idx=1);

	// Input                XW2D : An XWinNMR 2D data set
	//                      argc    : Number of arguments
	//                      argv    : Vector of argc arguments
	//                      argn    : Argument index
	//			idx	: Flag to request expt #
	//			           idx < 1: Ask for it
	//			           idx >=1: Use idx itself def (1)
	// Output               dirname : The parameter argn of array argc
	//                                is used to supply a directory name
	//                                from which the 2D data set resides
	//                                If the argument argn is not in argv,
	//                                the user is asked to supply a name
	//				  and experiment number
	//                                The directory is returned
	// Note                         : The directory should contain an
	//				  XWinNMR 2D data set.
};

#endif 							// XWin2D.h
