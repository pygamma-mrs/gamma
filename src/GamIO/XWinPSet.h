/* XWinPSet.h ***************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      XWinPSet                                   Interface		**
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
** sets. This class embodies a generic Bruker parameter file which will **
** contain any number of parameters detailing experimental acquisitions	**
** and processing.  Such a files are similar to GAMMA parameter sets,	**
** although in a different format.  When read, the Bruker format is	**
** converted to GAMMA format and all parameters stored in a GAMMA	**
** parameter set. Output will be generic using GAMMA parameter sets as	**
** well. It is up to derived classes to construct specificaly formatted	**
** output for particular Bruker/XWinNMR parameter files.		**
**                                                                      **
*************************************************************************/

#ifndef   XWinPSet_H_ 				// Is file already included?
#  define XWinPSet_H_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// This is the implementation
#  endif

#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <string>				// Include libstdc++ strings
#include <Basics/ParamSet.h>			// Include GAMMA parameter sets

class XWinPSet

// ----------------------------------------------------------------------------
// ------------------------------- STRUCTURE ----------------------------------
// ----------------------------------------------------------------------------
 
  {
  std::string       fname;		// Parameter set file name
  int          bigend;          // Flag for big- vs. little-endian
  ParameterSet pset;		// Parameter set matching Bruker file
 
// ----------------------------------------------------------------------------
// ----------------------------- PRIVATE FUNCTIONS ----------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                 XWinNMR Parameter Set File Error Handling
// ____________________________________________________________________________

/* These functions take care of any errors encountered when reading Bruker
   parameter files.
 
	Input		BruPSet  : XWinNMR Parameter Setition parameters (this)
			eidx    : Error index
                        pname   : Additional error message
                        noret   : Flag for linefeed (0=linefeed)
        Output          void    : An error message is output                 */

         void XWinPSeterror(int    eidx, int noret=0) const;
         void XWinPSeterror(int    eidx, const std::string& pname, int noret=0) const;
volatile void XWinPSetfatality(int eidx) const;
volatile void XWinPSetfatality(int eidx, const std::string& pname) const;


// ____________________________________________________________________________
// ii              XWinNMR Parameter File Default Parmameters
// ____________________________________________________________________________

virtual void SetDefaults();

public:
// ____________________________________________________________________________ 
// A             XWinPSet Parameter File Constructors, Destructor
// ____________________________________________________________________________
 
/* These are the constructors of the class handling Bruker XWinNMR parameter
   files.  This doesn't do anything in particular, it is the read functions
   that perform the work.  Thus, we only have a default constructor specified.
   The reading of the associated ASCII parameter file is done in one step,
   so we don't need anything complex.                                        */

XWinPSet();
XWinPSet(const std::string& name);
XWinPSet(const XWinPSet& XWP);
virtual ~XWinPSet();                                                             
void operator= (const XWinPSet& XWP);

// ____________________________________________________________________________
// B                  XWinPSet Parameter Access Functions
// ____________________________________________________________________________
 
/* These functions allow direct access to the few parameters that seem to be
   ubiquitous throughout Bruker parameter files.  

 Function         Output               Function             Output
 ---------  ----------------------    ----------   ------------------------
   name     File name as string         order       Stored data byte order   */

virtual std::string name()   const;	// File name             
//virtual bool   order()  const;	// Data byte order 

// ____________________________________________________________________________
// C                       XWinPSet Input Functions
// ____________________________________________________________________________

/* These functions will read in all parameters from a Bruker XWinNMR parameter
   file. By design, the Bruker parameter file is simply read into a GAMMA
   parameter set so that ALL parameters in the file are stored in this class.
   Were the Bruker files in GAMMA format then this would be done exclusively
   with class ParamSet, but since they are not we must handle the file parsing
   explicitly.

      Function                               Purpose
   ____________________________________________________________________________

      readPSet          Read in parameter set for class object.  This is done
                        special since Bruker format is NOT in GAMMA parameter
                        format so it must be parsed appropriately.
      getPar            Returns parameter found in the parameter set herein  */
 
bool readPSet(const std::string& filein, int warn=1);
bool readPSet(int warn=1);
ParameterSet getPSet() const;
int getPar(const std::string& pn, int& val,         int id=0, int wrn=0) const;
int getPar(const std::string& pn, double& val,      int id=0, int wrn=0) const;
int getPar(const std::string& pn, std::string& val, int id=0, int wrn=0) const;

// ____________________________________________________________________________
// D                       XWinPSet Output Functions
// ____________________________________________________________________________

/* There is nothing fancy here.  Since we don't know the particulars of any
   one Bruker parameter file, all we can do is dump out all of the parameters
   we have read.  We can keep these virtual so extended classes can make 
   output which is nicely formatted and detailing specific parameters.       */

        std::ostream& printPset(std::ostream& ostr) const;
virtual std::ostream& print(std::ostream& ostr, int full=0, int hdr=1) const;
friend  std::ostream& operator<< (std::ostream& ostr, const XWinPSet& BruPSet);

// ____________________________________________________________________________
// E                        XWinPSet Auxiliary Functions
// ____________________________________________________________________________

 /* These are just a group of functions that return strings indicating what
    the Bruker parameter values mean.  I just add things as I learn them here
    so that GAMMA output can remind me what all of these parameters are...   */

std::string FT_modS(int val) const;
std::string BYTORDS(int val) const;
std::string TDeffS(int val) const;


// ____________________________________________________________________________
// F                    XWinPSet BrukerLike Output Functions
// ____________________________________________________________________________

/* These functions perform ASCII output of processing parameters.  The format
   here is similar to issuing a "dpp" command in XWinNMR or using the menu
   choice "Output/Display status pars./(Acquisiton only OR Processing only)".
   That is, the output functions can be used to quickly generate a list in
   two columns that would be seen from "dpp" for example.                    */

int         brusize(double value) const;
std::string brustring(int befdec, int aftdec) const;
std::string bruform(double value, int aftdec) const;
void        bru8(std::ostream&  ostr, const std::string& label) const;
void        bru15(std::ostream& ostr, const std::string& label) const;

void        bru(std::ostream&   ostr, const std::string& label,
                        int value, const std::string& units, int type=0) const;
void        bru(std::ostream&   ostr, const std::string& label,
                        long value, const std::string& units, int type=0) const;
void        bru(std::ostream&   ostr, const std::string& label,
         const std::string& value, const std::string& units, int type=0) const;

void        bru(std::ostream&   ostr, const std::string& label,
 double value, const std::string& units, const std::string& parse, int type) const;
};

#endif 							// XWinPSet.h
