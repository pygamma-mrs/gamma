/* XWinOutd.h ***************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      XWinOutd                                   Interface		**
**                                                                      **
**      Copyright (c) 1999                                              **
**      Scott Smith                                                     **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**                                                                      **
**      $Header: 
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**      Description                                                     **
**                                                                      **
** The XWin* files provide an interface to Bruker XWinNMR (uxnmr) data  **
** sets. This class embodies a Bruker parameter file, outd, which seems	**
** to control output of NMR spectra within XWinNMR. Yet another ASCII	**
** paramter file, this has nothing to do with GAMMA and we will try	**
** and have minimal dealings with it - especially since it appears	**
** associated with processed data. A typical filename will be something	**
** such as dir/1/pdata/1/outd where dir	is the directory base that	**
** contains XWinNMR data.  Note that without the existence of an outd	**
** file one cannot display status parameters under the "Output" menu	**
** in XWinNMR. I don't recall if its absence disallows proper reading	**
** of a 1D data set within XWinNMR..... below if their 1D data set	**
** directory hierarchy:							**
**                          __ acqu  (changable parameter file)		**
**			   / 						**
**                        /___ acqus (static parameter file)		**
**			 /						**
**  expname -- expnum --< ---- fid (binary data)			**
**			 \						**
**			  \___ pdata -- 1 -- proc, procs, meta, outd	**
**									**
*************************************************************************/

#ifndef   XWinOutd_H_ 			// Is file already included?
#  define XWinOutd_H_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the implementation
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <string>			// Include libstdc++ strings
#include <Basics/ParamSet.h>		// Include GAMMA parameter sets
#include <GamIO/XWinPSet.h>             // Include Bruker parameter parsing

class XWinOutd: public XWinPSet

// ----------------------------------------------------------------------------
// ------------------------------- STRUCTURE ----------------------------------
// ----------------------------------------------------------------------------
 
  {
  std::string       oname;		// Parameter file name
  std::string       _TITLE;		// File title
  std::string       _JCAMPDX;		// JCAMP version (5.0)
  std::string       _DATATYPE;		// Type of file (Parameter Values)
  std::string       _DATE;		// When the file was made
  std::string       _CURPLOT;		// The default plotter
  std::string       _CURPRIN;		// The default printer
  std::string       _DFORMAT;		// The default format
  std::string       _LFORMAT;		// The default format
  std::string       _PFORMAT;		// The default format
  int               _SURQMSG;		// Unknown

 
// ----------------------------------------------------------------------------
// ----------------------------- PRIVATE FUNCTIONS ----------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________       
// i                      XWinNMR Outd File Error Handling
// ____________________________________________________________________________

/* These functions take care of any errors encountered when reading, writing,   
   and setting parameters in Bruker outd parameter files.
 
	Input		OutPar  : XWinNMR print/ploo parameters (this)
			eidx    : Error index
                        pname   : Additional error message
                        noret   : Flag for linefeed (0=linefeed)
        Output          void    : An error message is output                 */

         void XWinOutderror(int    eidx, int noret=0) const;
         void XWinOutderror(int    eidx, const std::string& pname, int noret=0) const;
volatile void XWinOutdfatality(int eidx) const;
volatile void XWinOutdfatality(int eidx, const std::string& pname) const;

// ____________________________________________________________________________ 
// ii                XWinNMR Outd Parameter File Defaults
// ____________________________________________________________________________

void SetDefaults();

public:
// ____________________________________________________________________________ 
// A             XWinOutd Parameter File Constructors, Destructor
// ____________________________________________________________________________
 
/* These are the constructors of the class handling Bruker XWinNMR output
   parameter files.  This doesn't do anything in particular, it is the write
   functions that perform the work.                                          */

XWinOutd();
XWinOutd(const std::string& name);
virtual ~XWinOutd();                                                             

// ____________________________________________________________________________
// B                  XWinOutd Parameter Access Functions
// ____________________________________________________________________________
 
/* These functions allow users to set some of the important parameters that 
   these files contain.                                                      */
 
void Printer(std::string pname);
void Plotter(std::string pname);

// ____________________________________________________________________________
// C                       XWinOutd Input Functions
// ____________________________________________________________________________

/* These functions will read in the 1D print/plot parameters from an XWinNMR
   parameter file, typically named outd.  By design, the Bruker parameter
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
     getPar            Returns parameter found in the parameter set herein
                       This function is inherited from base class XWinPSet. */

virtual bool read(const std::string& filein, int warn=1);
virtual bool read(int warn=1);
        bool parsePSet(int warn=1);
//bool getPar(const string& n,string& v,int i=-1, int w=0) const;


// ____________________________________________________________________________
// D                       XWinOutd Output Functions
// ____________________________________________________________________________

/* These function allow for output of NMR parameters directly into a Bruker
   XWinNMR ASCII parameter file (outd).                                      */
 
virtual int      write(const std::string& name, int warn=2);
virtual int      write(int warn=2) const;


// ____________________________________________________________________________
// E                    XWinOutd Standard Output Functions
// ____________________________________________________________________________

/* These function allow for a quick output of the parameter file contents.
   They don't have anything to do with output while running XWinNMR, rather
   users can just glance at outd parameters or store then in a small file.   */

//ostream& printPset(ostream& O) const;              INHERITED
virtual std::ostream& print(std::ostream& ostr, int full=0, int hdr=1) const;
friend std::ostream& operator<< (std::ostream& O, const XWinOutd& P);

};

#endif 								// XWinOutd.h

                                                                                



