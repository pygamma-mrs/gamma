/* XWinPSet.cc *************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      XWinPSet                                  Implementation	**
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
** contain any number of parameters detailing experimental acquisitions **
** and processing.  Such a files are similar to GAMMA parameter sets,   **
** although in a different format.  When read, the Bruker format is     **
** converted to GAMMA format and all parameters stored in a GAMMA       **
** parameter set. Output will be generic using GAMMA parameter sets as  **
** well. It is up to derived classes to construct specificaly formatted **
** output for particular Bruker/XWinNMR parameter files.                **
**                                                                      **
*************************************************************************/

#ifndef   XWinPSet_CC_			// Is file already included?
#  define XWinPSet_CC_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <GamIO/XWinPSet.h>		// Include the interface
#include <Basics/ParamSet.h>		// Include parameters sets
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <Basics/Gconstants.h>		// Include GAMMA constants (HZ2RAD)
#include <Basics/StringCut.h>		// Include GAMMA string parsing
#include <Basics/Isotope.h>		// Include GAMMA spin isotopes
#include <GamIO/BinIOBase.h>		// Include binary IO functions
#include <string>			// Include libstdc++ strings
#include <stdlib.h>
#ifndef _MSC_VER                        // If not using MSVC++ then we
 #include <sys/time.h>			// Include time and date access
 #include <unistd.h>			// Include POSIX getcwd function
#else                                   // and if using MSVC++ then I guess
 #include <time.h>			// Include time and date access
#endif

using std::string;			// Using libstdc++ strings
using std::ifstream;			// Using libstdc++ input file streams
using std::ostream;			// Using libstdc++ output streams


// ----------------------------------------------------------------------------
// ----------------------------- PRIVATE FUNCTIONS ----------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                      XWinNMR Parameter Set File Error Handling
// ____________________________________________________________________________

/* These functions take care of any errors encountered when reading Bruker
   parameter files.

   Input           BruPSet : XWinNMR Parameter Setition parameters (this)
		   eidx    : Error index
		   pname   : Additional error message
		   noret   : Flag for linefeed (0=linefeed)
   Output          void    : An error message is output                     */
 
void XWinPSet::XWinPSeterror(int eidx, int noret) const
  {
  string hdr("XWinNMR Parameter File");
  switch (eidx)
    {
    case 20: GAMMAerror(hdr,"Cannot Parse Bruker Parameters",noret); break; // (20)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }  

void XWinPSet::XWinPSeterror(int eidx, const string& pname, int noret) const
  {                                                                             
  string hdr("XWinNMR Parameter File");
  string msg;
  switch(eidx)
    {
    case 20:msg = string("Cannot Find ");
             GAMMAerror(hdr,msg+pname,noret);  break;                  // (20)
    case 21:msg = string("Cannot Write To ");
             GAMMAerror(hdr,msg+pname,noret);  break;                  // (21)
    case 22:msg = string("Cannot Set Parameter ");
             GAMMAerror(hdr,msg+pname,noret);  break;                  // (22)
    default: GAMMAerror(hdr, eidx, pname, noret); break;// Unknown Error  (-1)
    }
  }

volatile void XWinPSet::XWinPSetfatality(int eidx) const
  {                                                                 
  XWinPSeterror(eidx, 1);			// Normal non-fatal error
  if(eidx) XWinPSeterror(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

volatile void XWinPSet::XWinPSetfatality(int eidx, const string& pname) const
  {                                                                 
  XWinPSeterror(eidx, pname, 1);		// Normal non-fatal error
  if(eidx) XWinPSeterror(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii             XWinNMR Parameter Set Default Parmameters
// ____________________________________________________________________________
   
   void XWinPSet::SetDefaults()
     {}

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A             XWinPSet Parameter File Constructors, Destructor
// ____________________________________________________________________________

/* These are the constructors of the class handling Bruker XWinNMR parameter
   files.  This doesn't do anything in particular, it is the read functions
   that perform the work.  Thus, we only have a default constructor specified.
   The reading of the associated ASCII parameter file is done in one step, so
   we don't need anything complex.    					     */

XWinPSet::XWinPSet()
  {
  fname = string("Parameter Set");	// Set a default file name
  bigend = WeRBigEnd();			// Set output byte order
  SetDefaults();			// Set parameter defaults
  }

XWinPSet::XWinPSet(const string& name)
  {
  fname = name;				// Set a default file name
  bigend = WeRBigEnd();			// Set output byte order
  SetDefaults();			// Set parameter defaults
  }

XWinPSet::XWinPSet(const XWinPSet& XWP)
  {
  fname = XWP.fname;			// Copy file name
  bigend = XWP.bigend;          	// Copy big/little-endian
  pset = XWP.pset; 			// Copy parameter set
  }

XWinPSet::~XWinPSet() { }

void XWinPSet::operator= (const XWinPSet& XWP)
  {
  fname = XWP.fname;			// Copy file name
  bigend = XWP.bigend;          	// Copy big/little-endian
  pset = XWP.pset; 			// Copy parameter set
  }

// ____________________________________________________________________________
// B                  XWinPSet Parameter Access Functions
// ____________________________________________________________________________

/* These functions allow direct access to the few parameters that seem to be
   ubiquitous throughout Bruker parameter files.

 Function         Output               Function             Output
 ---------  ----------------------    ----------   ------------------------
   name     File name as string         size       No. complex points (int)  */

string XWinPSet::name()    const { return fname;     }
//bool   XWinPSet::order()   const { return fabs(byteordin); }

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
 
bool XWinPSet::readPSet(const string& filein, int warn)
  {
  fname = filein;				// Set internal filename
  bool TF = (readPSet(warn))?true:false;	// Read in all parameters
  if(!TF && warn)				// Output errors if trouble
    {
    XWinPSeterror(1, filein, 1);               // Problems with filein
    if(warn > 1) XWinPSeterror(3);             // Can't construct from pset
    else         XWinPSeterror(3,1);
    }
  return TF;
  }

bool XWinPSet::readPSet(int warn)
  {
  pset.clear();					// Insure empty pset
  ifstream inp(fname.c_str());			// Open up the file
  if(!inp.good())   	           	        // If file bad then exit
     {   
     if(warn)
       {
       XWinPSeterror(1, fname, 1);		//   Problems with file
       if(warn > 1) XWinPSetfatality(20); 	//   Cannot get Bruker params
       else         XWinPSeterror(20,1); 	//   either warn or fatality
       }
     return false;
     }   
  SinglePar par;				// GAMMA single parameter
  string pname, pdata;				// Strings for parameter fields
  string pcomment = "XWinNMR Parameter";	// Use this default comment 
  string input, temp="Bruker";			// Use this default name
  int cnt = 0;					// Counter for noname parameters
  while(Greadline(inp, input))			// Attempt to read a line
    {
    cutBlksXBlks(input, "##", 1);		//  Remove beginning ## if any
    cutBlksXBlks(input, "$", 1);		//  Remove beginning $ if any
    pname = cutString(input, 1);		//  Set the parameter name
    if(pname[pname.length()-1] == '=')		//  Remove any ending = from
      pname.resize(pname.length()-1);		//  the parameter name
    else if(pname[0] == '$')			//  Or perhaps it starts with $
      {						//  and then we will set a
      pname = temp + Gdec(cnt);			//  name for it
      cnt++;
      }
    pdata = input;				//  Set the parameter data
    par = SinglePar(pname,pdata,pcomment);	// Form a GAMMA parameter
    if(!pset.contains(par)) 			// Add parameter to parameter
      pset.push_back(par);			// list if not already included
    }
  return true;
  }

ParameterSet XWinPSet::getPSet() const { return pset; }

int XWinPSet::getPar(const string& pn, int& value, int eidx, int warn) const
  {
  string sval;
  int TF = pset.getString(pn,sval);		// Try and get the value
  if(!TF)					// If we can't find it then
    { 						// we issue warnings if
    if(warn) 					// desired
      { 					// 
      XWinPSeterror(eidx, 1);			//   Cant detemine size
      if(warn > 1) XWinPSetfatality(2, pn);	//   Fatal if warn > 1
      }
    }
  else value = atoi(sval.c_str());
  return TF;
  }


int XWinPSet::getPar(const string& pn, double& value, int eidx, int warn) const
  {
  string sval;
  int TF = pset.getString(pn,sval);		// Try and get the value
  if(!TF)					// If we can't find it then
    { 						// we issue warnings if
    if(warn) 					// desired
      { 					// 
      XWinPSeterror(eidx, 1);			//   Cant detemine size
      if(warn > 1) XWinPSetfatality(2, pn);	//   Fatal if warn > 1
      }
    }
  else value = atof(sval.c_str());
  return TF;
  }


int XWinPSet::getPar(const string& pn, string& value, int eidx, int warn) const
  {
  int TF = pset.getString(pn,value);		// Try and get the value
  if(!TF)					// If we can't find it then
    { 						// we issue warnings if
    if(warn) 					// desired
      { 					// 
      XWinPSeterror(eidx, 1);			//   Cant detemine size
      if(warn > 1) XWinPSetfatality(2, pn);	//   Fatal if warn > 1
      }
    }
  if(value.find('<') == 0)			// Clip any string encased
    {						// by <>!
    int l = value.length();
    if(int(value.find('>')) == l-1)
       value = value.substr(1, l-2);
    }
  return TF;
  }


// ____________________________________________________________________________
// D                         XWinPSet Output Functions
// ____________________________________________________________________________
 
 /* There is nothing fancy here.  Since we don't know the particulars of any
    one Bruker parameter file, all we can do is dump out all of the parameters
    we have read.  We can keep these virtual so extended classes can make
    output which is nicely formatted and detailing specific parameters.       */

ostream& XWinPSet::printPset(ostream& O) const { pset.print(O); return O; }

ostream& XWinPSet::print(ostream& ostr, int full, int hdr) const
  {
  string marg(16, ' ');
  string ls = string("\n") + marg;
  if(hdr)
    ostr << ls << "Bruker XWinNMR Parameter File: " << fname;
  ostr << ls << "Output Byte Ordering:       ";
  if(bigend) ostr << "big";
  else       ostr << "little";
  ostr << " endian";
  if(full) printPset(ostr);
  return ostr;
  }  

ostream& operator<< (ostream& O, const XWinPSet& A) { A.print(O); return O; }

// ____________________________________________________________________________
// E                        XWinPSet Auxiliary Functions
// ____________________________________________________________________________
 
 /* These are just a group of functions that return strings indicating what
    the Bruker parameter values mean.  I just add things as I learn them here
    so that GAMMA output can remind me what all of these parameters are...   */

string XWinPSet::FT_modS(int val) const
  { 
  string ret;
  switch(val)
    {
    case 0:  ret = "No FT Executed"; break;
    case 1:  ret = "Real FT of One Channel"; break;
    case 2:  ret = "Real FT of Quad Data"; break;
    case 3:  ret = "Complex FT of One Channel"; break;
    case 4:  ret = "Complex FT of Quad Data"; break;
    case 5:  ret = "Real Inverse FT of One Channel"; break;
    case 6:  ret = "Real Inverse FT of Quad Data"; break;
    case 7:  ret = "Complex Inverse FT of Quad Data"; break;
    default: ret = "Unknown Value"; 
    }
  return ret;
  }

string XWinPSet::BYTORDS(int val) const
  {
  string ret;
  switch(val)
    {
    case 0:  ret = "Little Endian"; break;
    case 1:  ret = "Big Endian";    break;
    default: ret = "Default";
    }
  return ret;
  }

string XWinPSet::TDeffS(int val) const
  {
  string ret("");
  if(!val) ret = string("all"); 
  return ret;
  }

// ____________________________________________________________________________
// F                    XWinPSet BrukerLike Output Functions
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
   be calculated.						             */

int XWinPSet::brusize(double dval) const
  {
  if(dval > 300) return 1;
  int ret = 0;
  while(fabs(dval) > 1)
    {
    ret++;
    dval /= 10;
    }
  if(ret == 0) ret++;
  else if(dval < 0) ret++;
  return ret;
  }

string XWinPSet::brustring(int befdec, int aftdec) const
  {
  string form = "%" + Gdec(befdec+1+aftdec)
              + "." + Gdec(aftdec) + "f";
  return form;
  }

string XWinPSet::bruform(double dval, int aftdec) const
  {
  int befdec = brusize(dval);
  return brustring(befdec, aftdec);
  }

void XWinPSet::bru8(ostream& ostr, const string& label) const
  {
  int namelen = 8;				// Length of label
  int padlen = namelen-label.length();		// Padding characters
  ostr << label << string(padlen, ' ');		// Write label, padding
  }

void XWinPSet::bru15(ostream& ostr, const string& label) const
  {
  int namelen = 15;				// Length of label
  int padlen = namelen-label.length();		// Padding characters
  int padlenst = padlen/2;			// Initial padding
  ostr << string(padlenst, ' ') << label 	// Write label, padding
       << string(padlen-padlenst, ' ');
  }

void XWinPSet::bru(ostream& ostr, const string& label, 
                     const string& value, const string& units, int type) const
  {
  string spc("  ");				// Default spacing
  if(!type) ostr << "\n    ";			// Begin new line 
  else      ostr << string(6, ' ');		// Space to new column
  bru8(ostr, label);				// Output label
  ostr << spc;					// Space after label
  if(value == "<>")       bru15(ostr, "");	// No output if <>
  else if(value == "< >") bru15(ostr, "");	// No output if < >
  else if(value == "")    bru15(ostr, "");	// This is just blanks
  else						// Output string value
    {						// but clip any starting <
    if(value.find('<') == 0)
      {
      int l = value.length();
      if(int(value.find('>')) == l-1)
           bru15(ostr, value.substr(1, l-2));
      else bru15(ostr, value);
      }
    else bru15(ostr, value);
    }
  ostr << spc;					// Space after value 
  bru8(ostr, units);				// Output units
  }

void XWinPSet::bru(ostream& ostr, const string& label, 
                   int value, const string& units, int type) const
  { 
	bru(ostr, label, Gdec(value), units, type); 
  }

void XWinPSet::bru(ostream& ostr, const string& label, 
                   long value, const string& units, int type) const
  { 
	bru(ostr, label, Gdec2(value), units, type); 
  }

void XWinPSet::bru(ostream& ostr, const string& label, 
      double value, const string& units, const string& parse, int type) const
  { 
  string sval = Gform(parse, value);
  bru(ostr, label, sval, units, type);
  }


#endif							// XWinPSet.cc
