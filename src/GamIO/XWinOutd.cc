/* XWinOutd.cc *************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      XWinOutd                                  Implementation	**
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

#ifndef   XWinOutd_CC_			// Is file already included?
#  define XWinOutd_CC_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <GamIO/XWinOutd.h>		// Include the interface
#include <GamIO/XWinPSet.h>             // Include Bruker parameter sets
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <stdlib.h>
#include <string>			// Include libstdc++ strings
#include <time.h>			// Include time and data access
#ifndef _MSC_VER                        // If not using MSVC++ then we
 #include <sys/time.h>			// Include time and date access
 #include <unistd.h>			// Include POSIX getcwd function
#else                                   // and if using MSVC++ then I guess
 #include <time.h>			// Include time and date access
#endif

using std::string;			// Using libstdc++ strings
using std::ofstream;			// Using libstdc++ output file streams
using std::ostream;			// Using libstdc++ output streams

// ----------------------------------------------------------------------------
// ----------------------------- PRIVATE FUNCTIONS ----------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                      XWinNMR Outd File Error Handling
// ____________________________________________________________________________
 
/* These functions take care of any errors encountered when reading, writing,
   and setting parameters in Bruker acquisition parameter files.

        Input		eidx    : Error index
        		pname   : Additional error message
        		noret	: Flag for linefeed (0=linefeed)
        Output		void    : An error message is output                 */
 
void XWinOutd::XWinOutderror(int eidx, int noret) const
  {
  string hdr("XWinNMR Outd Parameter File");
  string msg;
  switch (eidx)
    {
    case 25: GAMMAerror(hdr,"Cannot Write Parameters To File",noret);break; // (25)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }  

    
void XWinOutd::XWinOutderror(int eidx, const string& pname, int noret) const
  {                                                                             
  string hdr("XWinNMR Outd Parameter File");
  string msg;
  switch(eidx)
    {
    case 21:msg = string("Cannot Write To ") + pname;
             GAMMAerror(hdr,msg+pname,noret);  break;                  // (21)
    default: GAMMAerror(hdr, eidx, pname, noret); break;// Unknown Error  (-1)
    }
  }

volatile void XWinOutd::XWinOutdfatality(int eidx) const
  {                                                                 
  XWinOutderror(eidx, 1);			// Normal non-fatal error
  if(eidx) XWinOutderror(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

volatile void XWinOutd::XWinOutdfatality(int eidx, const string& pname) const
  {                                                                 
  XWinOutderror(eidx, pname, 1);		// Normal non-fatal error
  if(eidx) XWinOutderror(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }


// ____________________________________________________________________________
// ii                 XWinNMR Outd File Default Parmameters
// ____________________________________________________________________________

void XWinOutd::SetDefaults()
  {
  oname     = "outd";
  _TITLE    = "Parameter File,  Version XWIN-NMR  GAMMA";
  _JCAMPDX  = "5.0";
  _DATATYPE = "Parameter Values";
  time_t longtime;              // Need a time structure
  longtime = 0;        //

#ifdef _MSC_VER
	struct tm newtime;
  localtime_s(&newtime,&longtime);
	char date[31];
	asctime_s(date, 30, &newtime);
  _DATE = date;
#else
  struct tm *ptr;               // For setting current date
  ptr = localtime(&longtime);
  _DATE = asctime(ptr);
#endif

  _CURPLOT = string("hpdj550");	// Default plotter
  _CURPRIN = string("$hpdj550c");// The default printer
  _DFORMAT = "normdp";		// The default format
  _LFORMAT = "normlp";		// The default format
  _PFORMAT = "normpl"; 		// The default format
  _SURQMSG = 1;			// 
  }


 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A             XWinOutd Parameter File Constructors, Destructor
// ____________________________________________________________________________
 
/* These are the constructors of the class handling Bruker XWinNMR output
   parameter files.  This doesn't do anything in particular, it is the write
   functions that perform the work.                                          */

XWinOutd::XWinOutd()  { SetDefaults(); }
XWinOutd::XWinOutd(const string& name) { SetDefaults(); oname = name; }
XWinOutd::~XWinOutd() { }                                                         


// ____________________________________________________________________________
// B                  XWinOutd Parameter Access Functions
// ____________________________________________________________________________

/* These functions allow users to set any more important parameters that
   these files contain.                                                      */

void XWinOutd::Printer(string pname) { _CURPRIN = pname; }
void XWinOutd::Plotter(string pname) { _CURPLOT = pname; }

// ____________________________________________________________________________
// C                       XWinOutd Input Functions
// ____________________________________________________________________________

/* These functions will read in the 1D plot/print parameters from an XWinNMR
   parameter file, typically named outd.  By design, the Bruker parameter
   file is initially read into a GAMMA parameter set so that ALL parameters in
   the file are stored (class XWinPSet).  Subsequently, the parameters in the
   parameter set are parsed (herein) to obtain values of consequence to GAMMA
   and these are explicitly maintained variables in this class.

     Function                               Purpose
   ____________________________________________________________________________

     read               Read in parameter set for class object.  This is done
	    	        special since Bruker format is NOT in GAMMA parameter
		        format so it must be parsed appropriately.
   parsePSet            Converts parameters in internal pset to specific values
    getPar              Returns parameter found in the parameter set herein
		        This function is inherited from base class XWinPSet. */

bool XWinOutd::read(const string& filein, int warn)
  {
  oname = filein;				// Set internal filename
  bool TF = read(warn);                          // Read in all parameters
  if(!TF && warn)                               // Output errors if trouble
    {
    XWinOutderror(1, filein, 1);               // Problems with filein
    if(warn > 1) XWinOutderror(3);             // Can't construct from pset
    else         XWinOutderror(3,1);
    }
  return TF;
  }

bool XWinOutd::read(int warn)
  {
  XWinPSet::readPSet(oname, warn);		// Use base class to read
  return parsePSet(warn);			// Parse our essential params
  }

bool XWinOutd::parsePSet(int warn)
  {
  string pvalue;					// String value for # points
  if(getPar("TITLE",   pvalue,0,0))    _TITLE = pvalue;	// Try & get the title
  if(getPar("JCAMPDX", _JCAMPDX,0,0)) 			// Try & get JCAMP version
//    _JCAMPDX = atof(pvalue.c_str());
  if(getPar("DATATYPE",pvalue,0,0)) _DATATYPE = pvalue;	// Try & get DATATYPE 
  if(getPar("DATE",    pvalue,0,0))    _DATE  = pvalue;	// Try & get the date
  if(getPar("CURPLOT", pvalue,0,0)) _CURPLOT  = pvalue;
  if(getPar("CURPRIN", pvalue,0,0)) _CURPRIN  = pvalue;
  if(getPar("DFORMAT", pvalue,0,0)) _DFORMAT  = pvalue;
  if(getPar("LFORMAT", pvalue,0,0)) _LFORMAT  = pvalue;
  if(getPar("PFORMAT", pvalue,0,0)) _PFORMAT  = pvalue;
  if(getPar("SURQMSG", pvalue,0,0))
    _SURQMSG  = atoi(pvalue.c_str());
  return true;
  }

//int XWinOutd::getPar(const string& n,string& v,int i, int w) const;

// ____________________________________________________________________________
// D                         XWinOutd Output Functions
// ____________________________________________________________________________
 
/* These function allow for output of NMR parameters directly into a Bruker
   XWinNMR ASCII parameter file (outd).                                      */

int XWinOutd::write(const string& dname, int warn)
  { oname = dname; return write(warn); }

int XWinOutd::write(int warn) const
  {
  ofstream ofstr(oname.c_str());		// Open ASCII file for output
  if(!ofstr.good())				// If file bad then exit
     {
     if(warn)
       {
       XWinOutderror(2, 1);			// Output filestream problems
       XWinOutderror(25,1);			// Can't write parameters out
       if(warn==1) XWinOutderror(21,oname,1);	// Can't write parameters out
       else        XWinOutdfatality(21,oname);	// with out without die
       } 
     return 0;
     }  
  string nn("##");			// Bruker ## parameter line start
  string nns("##$");			// Bruker ##$ parameter line start
  string ss("$$");			// Bruker $$ parameter line start
  ofstr << nn << "TITLE= "    << _TITLE << "\n";
  ofstr << nn << "JCAMPDX= "  << _JCAMPDX << "\n";
  ofstr << nn << "DATATYPE= " << _DATATYPE << "\n";
  time_t longtime;				// Need a time structure
  longtime = 0;

#ifdef _MSC_VER
	struct tm newtime;
  localtime_s(&newtime,&longtime);
	char timebuf[31];
	asctime_s(timebuf, 30, &newtime);
  ofstr << ss << " " << timebuf;
#else
  struct tm *ptr;
  ptr = localtime(&longtime);
  ofstr << ss << " " << asctime(ptr);
#endif

	ofstr << nns << "CURPLOT= <" << _CURPLOT << ">\n";
  ofstr << nns << "CURPRIN= <" << _CURPRIN << ">\n";
  ofstr << nns << "DFORMAT= <" << _DFORMAT << ">\n";
  ofstr << nns << "LFORMAT= <" << _LFORMAT << ">\n";
  ofstr << nns << "PFORMAT= <" << _PFORMAT << ">\n";
  ofstr << nns << "SURQMSG= <" << _SURQMSG << "\n";
  ofstr << nn << "END= ";
  ofstr << "\n\n";
  return 1;
  }

// ____________________________________________________________________________
// E                    XWinOutd Standard Output Functions
// ____________________________________________________________________________

/* These functions allow for a quick output of the parameter file contents.
   They don't have anything to do with output while running XWinNMR, rather
   users can just glance at procs parameters or store then in a small file.  */

// ostream& XWinOutd::printPset(ostream& O) const		INHERITED

ostream& XWinOutd::print(ostream& ostr, int full, int hdr) const
  {
  string marg(16, ' ');
  string ls = string("\n") + marg;
  if(hdr)
    ostr << ls << "Outd Parameter File:       " << oname;
  XWinPSet::print(ostr, 0, 0);
  if(!full)
    {
    ostr << ls << "Plotter:                      " << _CURPLOT;
    ostr << ls << "Printer:                      " << _CURPRIN;
    ostr << ls << "Display Format:               " << _DFORMAT;
    ostr << ls << "Print Format:                 " << _LFORMAT;
    ostr << ls << "Parameter Plot Format:        " << _PFORMAT;
    ostr << ls << "SURQMSG:   			 " << _SURQMSG;
    }
  else
    {
    ostr << ls << "Title:                        " << _TITLE;
    ostr << ls << "JCAMP Version:                " << _JCAMPDX;
    ostr << ls << "Data Type:                    " << _DATATYPE;
    ostr << ls << "File Date:                    " << _DATE;
    ostr << ls << "Plotter:                      " << _CURPLOT;
    ostr << ls << "Printer:                      " << _CURPRIN;
    ostr << ls << "Display Format:               " << _DFORMAT;
    ostr << ls << "Print Format:                 " << _LFORMAT;
    ostr << ls << "Parameter Plot Format:        " << _PFORMAT;
    ostr << ls << "SURQMSG:   			 " << _SURQMSG;
    }
  return ostr;
  }  

ostream& operator<< (ostream& O, const XWinOutd& P) { P.print(O); return O; }


#endif							// XWinOutd.cc
