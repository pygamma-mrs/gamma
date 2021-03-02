/* ParamSet.cc **************************************************-*-c++-*-
**									**
**                                 G A M M A				**
**									**
**      Parameter Set                             Implementation	**
**									**
**      Scott A. Smith							**
**      Copyright (c) 1997                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                      	                       	**
**      $Header: $
**                                     	                        	**
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** The class ParameterSet contains a linked list of GAMMA parameters,	**
** the latter embodied in class SinglePar.  Each entry in a Parameter	**
** Set is of type SinglePar, see class SinglePar.                       **
**                                                                      **
*************************************************************************/

#ifndef   ParameterSet_cc		// Is this file included?
#  define ParameterSet_cc 1		// If not, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler? 
#    pragma implementation		// Then this is the implementaion
#  endif

#include <GamGen.h>			// Include OS specific stuff
#include <Basics/ParamSet.h>		// Include the implementation
#include <Basics/Gutils.h>		// Include GAMMA error handling
#include <Basics/StringCut.h>		// Include Gdec function
#include <Basics/SinglePar.h>		// Include single parameters
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>			// Include libstdc++ strings
#include <iostream>			// Include io streams (cout)
#include <list>				// Include STL list class
#include <fstream>			// Include libstdc++ filestreams
#include <vector>			// Include STL vector class

using std::string;			// Using libstdc++ strings
using std::list;			// Using libstdc++ STL lists
using std::vector;			// Using libstdc++ STL vectors
using std::ostream;			// Using libstdc++ output streams
using std::endl;			// Using libstdc++ line ends

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                   PARAMETER SET ERROR FUNCTIONS
// ____________________________________________________________________________

/*      Input                   pset    : Parameter set (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
                                pname   : Added error message for output
        Output                  none    : Error message output
                                          Program execution stopped if fatal

   The following error messages use the defaults set in the Gutils package
 
                Case                          Error Message

                (0)                     Program Aborting.....
  		(1)			Problems With Input File Stream
    		(2)			Problems With Output File Stream
                default                 Unknown Error                        */

void ParameterSet::Perror(int eidx, int noret) const
  {
  string hdr("Parameter Set");
  switch(eidx)
    {
    case 3: GAMMAerror(hdr, "Cannot Read Parameters From Input", noret);// (3)
            break;
    case 4: GAMMAerror(hdr, "Cannot Write Parameters To Output", noret);// (4)
            break;
    case 5: GAMMAerror(hdr, "Cannot Fill With Any Parameters",   noret);// (5)
            break;
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }  

void ParameterSet::Perror(int eidx, const string& pname, int noret) const
  {
  string hdr("Parameter Set");
  switch(eidx)
    {
    case 1:  GAMMAerror(hdr, 1, pname, noret);	break; // File Problems  (1)
    default: GAMMAerror(hdr, -1, pname, noret);	break; // Unknown Error	 (-1)
    }
  }  

volatile void ParameterSet::Pfatality(int eidx) const
  {                                                                 
  Perror(eidx, 1);			// First output the error
  if(eidx) Perror(0);			// Now output that its fatal
  GAMMAfatal();				// Clean exit from program
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A              PARAMETER SET CONSTRUCTORS/ASSIGNMENTS/DESTRUCTOR
// ____________________________________________________________________________

/* These functions are inherited from the ANSI standard "list" class in the C++
   library libstdc++.  I am listing them here so that I & other users don't 
   have to keep looking them up all the time.

   ParameterSet()				Empty Parameter Set
   ParameterSet(int N)				Parameter Set w/ N Parameters
   ParameterSet(int N, const SinglePar& par)	Parameter Set w/ N pars
   ParameterSet(const ParameterSet& pset)	Parameter Set copy of pset
   ParameterSet assign(N)			Assign N element
   ~ParameterSet()				Destructor of Parameter Set
*/

// ____________________________________________________________________________
// B                PARAMETER SET ITERATORS & MEMBER ACCESS
// ____________________________________________________________________________

/* These functions are inherited from the ANSI standard "list" class in the C++
   library libstdc++.  I am listing them here so that I & other users don't 
   have to keep looking them up all the time.

   list<SinglePar>::iterator ParameterSet::begin()	Pointer to 1st element
   list<SinglePar>::iterator ParameterSet::end()	Pointer to last element
   SinglePar                 ParameterSet::front()	First element
   SinglePar                 ParameterSet::back()	Last element
*/

// ____________________________________________________________________________
// C                PARAMETER SET LIST & QUEUE OPERATIONS
// ____________________________________________________________________________

/* These functions are inherited from the ANSI standard "list" class in the C++
   library libstdc++.  I am listing them here so that I & other users don't 
   have to keep looking them up all the time.

   ParameterSet::push_back(const SinglePar& par)	Add par to list end
   ParameterSet::pop_back()				Remove par at list end
   ParameterSet::push_front(const SinglePar& par)	Add par to list start
   ParameterSet::pop_front(const SinglePar& par)	Remove par at list start

   ParameterSet::insert(iterator p, SinglePar& par)	Add par before p
   ParameterSet::erase(iterator p)			Remove par at p
   ParameterSet::clear()				Remove all list entries
*/

// ____________________________________________________________________________
// D                PARAMETER SET ADDITIONAL QUEUE OPERATIONS
// ____________________________________________________________________________

/* These functions are inherited from the ANSI standard "list" class in the C++
   library libstdc++.  I am listing them here so that I & other users don't 
   have to keep looking them up all the time.

   int  ParameterSet::size()				Number of entries
   bool ParameterSet::empty()				TRUE if pset empty
   bool ParameterSet::operator==(ParameterSet pset)	TRUE if psets equal
   bool ParameterSet::operator!=(ParameterSet pset)	TRUE if psets not equal
*/

// ____________________________________________________________________________
// E                     PARAMETER SET AUXILIARY FUNCTIONS
// ____________________________________________________________________________

/* These functions are "list" type of functions that have been added to make
   the list of single parameters (i.e. parameter list) do simple things needed
   for ready access of single parameters.

    Function   Arguments                     Result
   ----------  ---------  -----------------------------------------------------
    contains    string    Returns true/false if parameter with name is in pset
       "       SinglePar  Returns true/false if single parmeters is in pset
     seek       string    Returns iterator in pset for parameter with name
       "       SinglePar  Returns iterator in pset for single parameter
     strip       int      Returns parameter set with parameters named [#]name
   countpar   string,int  Counts contiguous parameters in pset of name name(#)

   Note that in the function "strip" the returned parameter set contains
   parameters whos names no longer are prefixed with [#].                   */

int ParameterSet::contains(const string& pname) const
  { return (seek(pname) != end()); }
int ParameterSet::contains(const SinglePar& par) const
  { return contains(par.name()); }

ParameterSet::const_iterator ParameterSet::seek(const string& pname) const
  {
  SinglePar par;				// Working single parameter
  list<SinglePar>::const_iterator item=begin();	// Begin at input pset start
  while(item != end())                          // and loop over all of its
    {                                           // parameters
    par = *(item);				//   This is a parameter
    if(pname == par.name()) return item;	//   Return item is name match
    item++;					//   Go to the next Pix
    }
  return item;
  }
 
ParameterSet::const_iterator ParameterSet::seek(const SinglePar& par) const
  { return seek(par.name()); }

ParameterSet ParameterSet::strip(int idx) const
  {
  ParameterSet subpset;                      // The pset to be returned
  string Nidx=string("[")+Gdec(idx)+string("]");	// Parameter name prefix
  int len = Nidx.length();                      // Prefix length
  SinglePar par;				// Temporary single parameter
  string pname;					// Temporary parameter name
  list<SinglePar>::const_iterator item=begin();	// Begin at input pset start
  while(item != end())                          // and loop over all of its
    {                                           // parameters
    par = *(item);				//   This is a parameter
    pname = par.name();                         //   & this is its name
    if(!pname.find(Nidx))			//   See if name starts with
      {						//   [#] where # is idx
      par.name(pname.substr(len)); 		//      Cut out prefix
      subpset.push_back(par);                   //      Add it to sub pset
      }
    item++;					//   Go to the next Pix
    }  
  return subpset;
  }  

int ParameterSet::countpar(const string& pnamein, int idx0)
  {
  int cnt=0, found=1, fiter;			// Set counter to zero
  string suffi("(");				// Suffix start
  string sufff(")");				// Suffix finish
  string sidx = suffi + Gdec(idx0) + sufff;	// Parameter name index 
  string pname = pnamein + sidx;		// Parameter name to look for
  list<SinglePar>::iterator item;		// Item in the parameter set
  SinglePar par;				// Single parameter
  while(found)					// Look for pname in pset
    {
    fiter = 0;					//  Havent found pname yet
    item = begin();				//  Begin at input pset start
    while((item != end()) && !fiter)		//  and loop over all of its
      {						//  parameters
      par = *(item);				//  This is a parameter
      if(pname == par.name())			//  Compare name to sought name
        fiter++;
      }
    if(fiter)					// If we found pname then
      {						// we'll update the counter
      cnt++;					// and reset for next pname
      sidx = suffi + Gdec(idx0+cnt) + sufff;
      }
    else found = 0;
    }
  return cnt;
  }

// ____________________________________________________________________________
// F                    PARAMETER SET OUTPUT FUNCTIONS
// ____________________________________________________________________________

/* These are functions that output formatted information concerning the
   Parameter set to a specified output stream of file. The latter function,
   write, will produce an ASCII file that is self-readable by this class.

                Input           pset : Parameter set (this)
                                ostr : Output ASCII file stream
                                file : ASCII output file (making pset file)
                                warn : Warning level flag for write fail
                Return          void : pset is sent to the output stream 
                                int  : T/F write to output pset file.        */

vector<string> ParameterSet::printStrings() const
  {
  vector<string> PStrings;		// Strings for parameters
  int np = int(size());				// # of parameters
  if(!np) return PStrings;			// Nothing if no parameters
  
//                    Arrays For Strings We Will Oputput

  string sts[5] = { "Integer", "Double",
                         "String",  "Coordinate",
                         "Tensor" };
  vector<string> PNs;			// For parameter names
  vector<string> PTs;			// For parameter types
  vector<string> PVs;			// For parameter values
  vector<string> PCs;			// For parameter comments
  vector<int>         MaxLens;		// For column lengths
  string pline;				// For parameter line

//           Fill String Arrays & Find Maximum Column Lengths

  list<SinglePar>::const_iterator i;
  i = begin();					// Begin at input pset start
  int k=0;					// Temp for array index
  int maxpn = 0;				// Max parameter name length
  int maxpt = 0;				// Max parameter type length
  int maxpv = 0;				// Max parameter value length
  int maxpc = 0;				// Max parameter comment length

  while(i != end())				// So long as Pix is OK
    {
    PNs.push_back((*i).name());			//   Store ith parameter name
    maxpn = gmax(maxpn, int((PNs[k]).length()));//   Set maximum name length
    PTs.push_back(sts[(*i).type()]);		//   Store ith parameter type
    maxpt = gmax(maxpt, int((PTs[k]).length()));//   Set maximum value length
    PVs.push_back((*i).data());			//   Store ith parameter value
    maxpv = gmax(maxpv, int((PVs[k]).length()));//   Set maximum value length
    PCs.push_back((*i).state());		//   Store ith parameter comment
    maxpc = gmax(maxpc, int((PCs[k]).length()));//   Set maximum comment length
    i++;					//   Get next Pix
    k++;
    }

  MaxLens.push_back(gmin(maxpn, 20));
  MaxLens.push_back(maxpt);
  MaxLens.push_back(gmin(maxpv, 30));
  MaxLens.push_back(gmin(maxpc, 50));

//                       Set Up Column Headers To Use
//        Their Widths Are Adjusted To Handle Parameters In Set

  vector<string> chdrs;
  chdrs.push_back("Name"); 			// 1st column header
  chdrs.push_back("Type"); 			// 2nd column header
  chdrs.push_back("Value"); 			// 3rd column header
  chdrs.push_back("Comment"); 			// 4th column header
  int nc = int(chdrs.size());			// Number of columns
  string tmp;
  int tl;
  for(k=0; k<nc; k++)
    {
    tl = chdrs[k].length();
    if(tl < MaxLens[k])
      {
      tmp =  string((MaxLens[k]-tl)/2, ' ');
      tmp += chdrs[k];
      tmp += string(MaxLens[k]-tmp.length(), ' ');
      chdrs[k] = tmp;
      }
    }

//                        Set Column Header Strings 

  int cspcl = 2;				// Spacing between columns
  string cspc(cspcl, ' ' );		// For spacing between cols
  int ll = -cspcl;				// Total line length
  for(k=0; k<nc; k++)				// Get total line length
    { ll+= chdrs[k].length() + cspcl; }
  string spc(40-ll/2, ' ');		// Spacer to center line

  pline = spc;					// Space to 1st column start
  for(k=0; k<nc; k++)				// Loop columns and add
    {
    pline += chdrs[k];				//   Add the column title
    if(k<nc-1) pline += cspc;			//   Space to next column
    }
  PStrings.push_back(pline);			// Store title header

  pline = spc;					// Space to 1st column start
  for(k=0; k<nc; k++)				// Loop columns and output
    {
    pline +=  string(chdrs[k].length(), '-');
    if(k<nc-1) pline += cspc;
    }
  PStrings.push_back(pline);			// Store header underlines

//                  Make And Store Strings For Each Parameter
//            (Each String Will Be A Line Containing Parameter Info)

  int ml;
  for(k=0; k<np; k++)				// Loop over all parameters
    {
    pline = spc;				// Space to 1st column start
    pline += PNs[k];				// Set parameter name
    tl = PNs[k].length();			// Length of name
    ml = MaxLens[0];				// Length of column
    if(tl<ml) pline += string(ml-tl, ' ');	// Space to end of column
    pline += cspc;				// Space to next column
    pline += PTs[k];				// Set parameter type
    tl = PTs[k].length();			// Length of name
    ml = MaxLens[1];				// Length of column
    if(tl<ml) pline += string(ml-tl, ' ');	// Space to end of column
    pline += cspc;				// Space to next column
    pline += PVs[k];				// Set parameter value
    tl = PVs[k].length();			// Length of name
    ml = MaxLens[2];				// Length of column
    if(tl<ml) pline +=  string(ml-tl, ' ');// Space to end of column
    pline += cspc;				// Space to next column
    pline += PCs[k];				// Set parameter comment
    PStrings.push_back(pline);			// Store parameter string
    }

  return PStrings;
  }

ostream& ParameterSet::print(ostream& ostr) const
  {
//              Output Header (All If Empty Parameter Set)

  string hdr("Parameter Set");		// Basic header
  int hl;					// Header length
  int np = int(size());				// # of parameters
  if(!np)					// This if no parameters
    {						// Header is centered
    hdr = "Empty " + hdr;
    hl = hdr.length();
    ostr << "\n" << string(40-hl/2, ' ')
         << hdr;
    return ostr;
    }

  hl = hdr.length();				// This if parameters
  ostr << "\n" << string(40-hl/2, ' ')	// Header is centered
       << hdr;
  ostr << endl;

  vector<string> PStrings;		// Array of strings to print
  PStrings = printStrings();			// Get our print strings
  unsigned ns = PStrings.size();		// No. of strings to print
  for(unsigned i=0; i<ns; i++)			// Loop over strings and send
    ostr << PStrings[i] << endl;		// to our string stream

  ostr << endl;				// Another line feed
  return ostr;					// Return the output stream
  }

ostream& operator<< (ostream& ostr, const ParameterSet& pset)
  { return pset.print(ostr); }

bool ParameterSet::write(const string &fileout, int warn) const
   {
   std::ofstream ofstr(fileout.c_str());	// Open fileout for output
   int w2 = 0;					// Warning level for next call
   if(warn) w2 = 1;				// Use level 1 warning
   if(!write(ofstr, w2)) 
     {
     if(warn) Perror(1,fileout,1);		// Problems with output file
     if(warn==1) Perror(4);			// Can't write parameters out
     else        Pfatality(4);			// with out without die
     return false;
     }
   return true;
   }

bool ParameterSet::write(std::ofstream& ofstr, int warn) const
   {
   if(!ofstr.good())				// If file bad then exit
     {
     if(warn)
       {
       Perror(2, 1);				// Output filestream problems
       if(warn==1) Perror(4);			// Can't write parameters out
       else        Pfatality(4);		// with out without die
       }
     return false;
     }
   list<SinglePar>::const_iterator i=begin();	// Begin at input pset start
   while(i != end())
     {
     (*i).write(ofstr);
     i++;
     }
   return true;
   }


// ____________________________________________________________________________
// G                    PARAMETER SET INPUT FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//      Direct Read of Parameter Set From An ASCII File Or Input Stream
// ----------------------------------------------------------------------------

	// Input		pset    : Parameter set (this)
	// 			filein  : Input filename
        //                 Or   inp     : Input filestream
	//			fflag   : Flag for fatal warnings.
	//				     0 : No error messages
	//				     1 : Non-fatal messages
	//				  else : Messages, stop execution
	// Output		int     : True if file filein or input stream
	//				  inp has been scanned for parameters.
	// Note				: Parameter set filled with parameters
	//				  read from the input file filein
	//				  or input file stream inp

bool ParameterSet::read(const string &filein, int fflag)
  {
  std::ifstream inp(filein.c_str());		// Open filename for input
  if(!read(inp, fflag?1:0))			// Use overloaded function
    {
    if(fflag)					//    Output any error messages
      {						//    as desired via fflag
      Perror(1,filein,1);			//      Problems with file
      (fflag>1)?Pfatality(5):Perror(5,0);	//      Can't read pset
      }
    return false;
    }
  inp.close();					// Close the input file
  return true;					// Return we succeeded
  }

bool ParameterSet::read(std::ifstream& inp, int fflag)
  {
  clear();					// Clear any existing parameters
  if(!inp)					// If file bad then exit
    {   
    if(fflag)					//    Output any error messages
      {						//    as desired via fflag
      Perror(1, 1);				//      Problems with input stream
      (fflag>1)?Pfatality(3):Perror(3,1);	//      Can't read pset
      } 
    return false;				//    Exit that we failed
    }  
  SinglePar par;				// Here's a single parameter
  int fend;
  while(!inp.eof())				// While successfully reading
    {						// look for valid parameters
    fend = par.read(inp);			// Try & read a parameter
    if(fend>0 && !contains(par)) push_back(par);// Add to list if good & unique
    }
  return true;
  }

// ____________________________________________________________________________
// H                    PARAMETER SET INTERACTIVE FUNCTIONS
// ____________________________________________________________________________
   
        // Input                pset    : Parameter set (this)
        //                      argc    : Number of arguments
        //                      argv    : Vector of argc arguments
        //                      argn    : Argument index
        // Output               filename: The parameter argn of array argc
        //                                is used to supply a filename
        //                                from which the parameter set is read
        //                                If the argument argn is not in argv,
        //                                the user is asked to supply a filename
        //                                The filename is returned
        // Note                         : The file should be an ASCII file
        //                                containing recognized sys parameters
        // Note                         : The parameter set is modifed (filled)

string ParameterSet::ask_read(int argc, char* argv[], int argn)
  {
  string fn;				// Name of spin system file
  query_parameter(argc, argv, argn,             // Get filename from command
       "\n\tASCII Parameter Set filename? ",fn);// Or ask for it
  read(fn);					// Read system from filename
  return fn;					// Give back the filename
  }


// -------------------------------------------------------------------------
//                         Parameter Value Functions
// -------------------------------------------------------------------------

        // Input                pset    : Parameter set (this)
        //                      name    : Name of parameter desired
        //                      value   : Value of parameter desired
        // Output               TF      : True if parameter named "name" exists
        //                                in pset and sucessfully parsed as
        //                                an integer, value set to param. value

bool ParameterSet::getParameter(const string& name, string& value) const
  { return getString(name,value); }

bool ParameterSet::getParameter(const string& name, int& value) const
  { return getInt(name,value); }

bool ParameterSet::getParameter(const string& name, double& value) const
  { return getDouble(name,value); }

bool ParameterSet::getString(const string& name, string& value) const
  {
  list<SinglePar>::const_iterator item;	// Pix into pset
  item = seek(name);                            // See if parameter in pset
  if(item == end()) return false;		// If not just bail
  if((*(item)).type() != 2) return false;		// Insure its a string
  value = (*(item)).data();                       // Set value to parameter data
  return value.length()?true:false;		// Fail only if no data
  }

bool ParameterSet::getInt(const string& name, int& value) const
  {
  list<SinglePar>::const_iterator item;		// Pix into pset
  item = seek(name);                            // See if parameter in pset
  if(item == end()) return false;		// If not just bail
  if((*(item)).type() != 0) return false;		// Insure its an integer
  string sval = (*(item)).data();			// Get value of parameter data
  if(!sval.length()) return false; 		// Fail only if no data
  value = atoi(sval.c_str());			// Convert to integer
  return true;
  }

bool ParameterSet::getDouble(const string& name, double& value) const
  {
  list<SinglePar>::const_iterator item;		// Pix into pset
  item = seek(name);                            // See if parameter in pset
  if(item == end()) return false;		// If not just bail
  if((*(item)).type() != 1) return false;		// Insure its a double
  string sval = (*(item)).data();			// Get value of parameter data
  if(!sval.length()) return false;		// Fail only if no data
  value = atof(sval.c_str());			// Convert to integer
  return true;
  }

#endif						// ParamSet.cc


