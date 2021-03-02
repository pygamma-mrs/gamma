/* MultiAux.cc **************************************************-*-c++-*-
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      Multiple System Auxiliary Classes	Implementations		**
**                                                                      **
**      Copyright (c) 1995                                              **
**      Nikolai Skrynnikov                                              **
**      Dr. Scott A. Smith                                              **
**      Copyright (c) 1996                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header: $
**                                                                      **
**************************************************************************

**************************************************************************
**                                                                      **
** This file contains auxiliary classes for mulitple spin systems.      **
** Specifically it handles the mapping of a spin in one system          **
** onto the spin of another system, the systems involved being          **
** those which are part of multi_sys.                                   **
**                                                                      **
*************************************************************************/

#ifndef   MultiAux_cc_			// Is the file already included?
#  define MultiAux_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif


#include <MultiSys/MultiAux.h>		// Include header file
#include <string>			// Know about libstdc++ strings
#include <cstdlib>                      // Include functions atoi, atof, etc.
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Basics/StringCut.h>		// Include string cutting 
#include <Basics/Gutils.h>		// Include GAMMA error messages
#include <fstream>

using std::string;			// Using libstdc++ strings
using std::ostream;			// Using libstdc++ output streams
using std::list;			// Using libstdc++ STL list
using std::vector;			// Using libstdc++ STL vectors

/********************************************************************************/
/********************************************************************************/
/*                             CLASS SPIN PAIR FUNCTIONS 			*/
/********************************************************************************/
/********************************************************************************/

//_________________________________________________________________________________
// A           SPIN PAIR CONSTRUCTORS, ASSIGNMENT, AND DESTRUCTOR
//_________________________________________________________________________________
 
// The spin pair constructors commented out reside in the header file.
 
//spin_pair::spin_pair()
//spin_pair::spin_pair(int, int, int, int);
//void spin_pair::operator = (const spin_pair& Sp)
//spin_pair::~spin_pair() 
 

spin_pair::spin_pair(const spin_pair& Sp)
 
        // Input                Sp      : A spin pair
        // Output               void    : New spin pair identical to Sp
 
  {
  sub1 = Sp.sub1;
  sub2 = Sp.sub2;
  sp1 = Sp.sp1;
  sp2 = Sp.sp2;
  }
 

spin_pair::spin_pair(const string& SSP)

        // Input                SSP     : A String defining a spin pair
        // Output               void    : The process is constructed
	// Note				: It is assumed the input string
	//				  has the form
	//
	//				  (C1)S1(C2)S2 --> (sub1)sp1(sub2)sp2
	//
	// Note				: No checks are preformed to insure
	//				  that any of these indices are valid
 
  {
  string sval = SSP;			// Copy of input string
  cutBlksXBlks(sval, "("); 		// Cut initial (
  sub1 = atoi(cutInt(sval).c_str());	// Cut integer (1st component)
  cutBlksXBlks(sval, ")"); 		// Cut out )
  sp1 = atoi(cutInt(sval).c_str());	// Cut out an integer (spin of C1)
  cutBlksXBlks(sval, "("); 		// Cut initial (
  sub2 = atoi(cutInt(sval).c_str());	// Cut integer (rt component)
  cutBlksXBlks(sval, ")"); 		// Cut out )
  sp2 = atoi(cutInt(sval).c_str());	// Cut out an integer (spin of C2)
  }


//_________________________________________________________________________________
// B                          SPIN PAIR ACCESS FUNCTIONS
//_________________________________________________________________________________
 
int spin_pair::Sub1() const { return sub1; }
int spin_pair::Sub2() const { return sub2; }
 
        // Input                spair   : A spin pairing (this)
        // Output               sub#    : The component index of spin #
 
int spin_pair::Spin1() const { return sp1; }
int spin_pair::Spin2() const { return sp2; }
 
        // Input                spair   : A spin pairing (this)
        // Output               sp#     : The index of spin #
 
//_________________________________________________________________________________
// C                              SPIN PAIR OUTPUT
//_________________________________________________________________________________

 
ostream& spin_pair::print(ostream& ostr) const
 
        // Input                spair   : A spin pairing (this)
        //                      ostr    : Output stream
        // Output               none    : Spin pairing writtn to output stream
        ///F_list print                 - Write system to output stream 

  {
  ostr << "Component "<< sub1 << ", Spin # " << sp1
       << " <---> "
       << "Component "<< sub2 << ", Spin # " << sp2;
  return ostr;
  }
 
 
ostream& operator<< (ostream& ostr, const spin_pair &Sp)
 
        // Input                ostr    : Output stream;
        //                      Sp	: A spin pairing
        // Output                       : Modifies output stream
        ///F_list <<                    - Standard Output
 
  { return Sp.print(ostr); }


/********************************************************************************/
/********************************************************************************/
/*                           CLASS PROCESS FUNCTIONS 				*/
/********************************************************************************/
/********************************************************************************/

//_________________________________________________________________________________
// i                       CLASS PROCESS ERROR HANDLING
//_________________________________________________________________________________


	// Input			PRO	: A process (this)
        // 				eidx	: Error flag
        //				noret   : Flag for return (0=return)
        //                              pname   : String included in error 
        // Output			none 	: Output process error message


void process::XPerror(int eidx, int noret) const
  {
  string hdr("Exchange Process");
  string msg;
  switch (eidx)
    {
//  case 0:  Program Aborting                                           // (0)
//  case 1:  Problems With Input File Stream                            // (1)
//  case 2:  Problems With Output File Stream                           // (2)
//  case 3:  Can't Construct From Parameter Set                         // (3)
//  case 4:  Cannot Construct From Input File                           // (4)
//  case 5:  Cannot Write To Parameter File                             // (5)
//  case 6:  Cannot Write To Output FileStream                          // (6)
//  case 7:  Cannot Produce Any Output                                  // (7)
//  case 8:  Problems With Parameter Set                                // (8)
//  case 9:  Problems During Construction                               // (9)
//  case 10: Cannot Access Internal Component                           // (10)
//  case 11: Cannot Read From Input FileStream                          // (11)
//  case 12: End of File Reached                                        // (12)
//  case 13: Cannot Open ASCII File For Read                            // (13)
//  case 14: Cannot Read From Input File                                // (14)
//  case 15: Cannot Write To Output File                                // (15)
//  default: Unknown Error", noret);                                    // (-1)
    default: GAMMAerror(hdr, eidx, noret); return; break;
    case 30: msg=string("Can't Read From Parameter Set");        break; // (30)
    case 31: msg=string("No Exchange Process in Parameter Set"); break; // (31)
    case 32: msg=string("No Exchange Rate in Parameter Set");    break; // (32)
    case 33: msg=string("No Left Side Components In Exchange");  break; // (33)
    case 34: msg=string("No Right Side Components In Exchange"); break; // (34)
    case 35: msg=string("Cannot Parse Right Side Exchange Cmps");break; // (35)
    case 36: msg=string("Cannot Parse LeftSide Exchange Comps"); break; // (36)
    }
  GAMMAerror(hdr, msg, noret);
  }



/*
    case 6:								// (6)
      cout << "Cannot Parse Any Spin Pair Mappings";
      break;
    case 7:								// (7)
      cout << "Cannot Re-Parse Spin Pair Mappings! This Shouldn't Happen.";
      break;
    case 10:								// (10)
      cout << "Error During Construction From Parameter Set";
      break;
    }
  if(noret)  cout << ".";
  else       cout << ".\n";
  }    
*/


void process::XPerror(int eidx, const string& pname, int noret) const
 
  {
// case 1:  Problems With File pname                                    // (1)
// case 2:  Cannot Read Parameter pname                                 // (2)
// case 3:  Invalid Use Of Function pname                               // (3)
// case 4:  Use Of Deprecated Function pname                            // (4)
// case 5:  Please Use Class Member Function pname                      // (5)
// default: Unknown Error pname                                         // (-1)
  string hdr("Exchange Process");
  string msg;
  switch (eidx)
    {
    default: GAMMAerror(hdr, eidx, pname, noret); return; break;
    case 20: msg=string("Cannot Read Parameter ") + pname; break;      	// (20)
/*
    case 2:								// (2)
      cout << "Accessed LHS Component " << pname << " Out of Range";
      break;
    case 3:								// (3)
      cout << "Accessed RHS Component " << pname << " Out of Range";
      break;
    case 4:								// (4)
      cout << "Attempted To Obtain Parameter " << pname;
      break;
*/
    }
  GAMMAerror(hdr, eidx, noret);
  }
 
 
volatile void process::XPfatal(int eidx) const
 
	// Input			PRO	: A process (this)
        // 				eidx	: Error flag
        // Output                       none    : Output process error message
        //					  Program execution stopped
 
  {
  XPerror(eidx, 1);			// Output the error message
  if(eidx) XPerror(0);			// State that its fatal
  GAMMAfatal();				// Clean exit from program
  }

//_________________________________________________________________________________
// ii                CLASS EXCHANGE PROCESS PARAMETER SET PARSING
//_________________________________________________________________________________

/* These functions allow for an exchange process to be set from parameters in a
   specified parameter set.                                                      */


//---------------------------------------------------------------------------------
//------------------------- Read in the Process Definition ------------------------
//---------------------------------------------------------------------------------

/* This will be a parameter such as    Exch(0)  (2) : (0<=>1+2) - Exchange sheme

   The string value defines the components (subsystems) involved in the exchange
   process. We only get the process definition, it doesn't bother to parse it.   */

bool process::getExch(const ParameterSet& pset, int idx,
                                                     string& exch, bool warn) const
  {
  exch = string("");                            // We don't know any exchange
  ParameterSet::const_iterator item;         // A pix into parameter list
  string Nidx = "";                             // Name addition per index
  if(idx >= 0) Nidx                             // If index exists then set up
    += string("(")+Gdec(idx)+string(")");       // the parameter name to append
  string pstate, pname = string("Exch") + Nidx;		// Parameter name for process
  item = pset.seek(pname);			// Pix in pset for Exch(i)
  if(item != pset.end())			// Retrieve the process (string)
    {
    (*item).parse(pname,exch,pstate);		//   Process as a string
    return true;				//   Return we were successful
    }
  if(warn)					// If we didnt find the exchange
    XPerror(20, pname, 1); 			// process definition, warn
  return false;
  }

//---------------------------------------------------------------------------------
//-------------- Parse The Components Involved In The Exchange Process ------------
//---------------------------------------------------------------------------------

/* Here we assume we have the string from parameter Exch, e.g. (0<=>1+2). We need
   to parse out the integer values on each side of <=>, but lying between (). 
   Also we must keep track of which are on the left and which are on the right.  */

bool process::parseExch(string& Exval, 
                               vector<int>& lhs, vector<int>& rhs, bool warn) const
  {
//		First Get Number of LHS and RHS Components in Exchange

   cutBlksXBlks(Exval, "(");			// Cut initial ( from string
   string SLHS = cutBlksXBlks(Exval, "<=>");	// Next cut past <=> part
   int nc = SLHS.find("<=>");			// Number of chars before <=>
   SLHS = SLHS.substr(0, nc);			// Get string before <=>
   string SRHS = Exval;				// This is the right side string
   if(!SLHS.length())				// If there isn't anything on
     { 						// the left of <=>, then we don't
     XPerror(33, 1);				// know left hand side components
     return false;
     }
   if(!SRHS.length())				// If there isn't anything on
     { 						// the right of <=>, then we don't
     XPerror(34, 1);				// know right hand side components
     return false;
     }
   int cl = int(SLHS.length());			// Chars in left hand side string
   int nlhs  = 1;				// Spins on left hand side
   int i;
   for(i=0; i<cl; i++)				// Look for '+' and for each there
     if(SLHS[i] == '+') nlhs++;			// should be an additional spin
   int nrhs = 1;				// Spins on right hand side
   int cr = int(SRHS.length());			// Chars on right hand side string
   for(i=0; i<cr; i++) 				// Look for '+' and for each there
     if(SRHS[i] == '+') nrhs++;			// should be an additional spin

//		   Now Set the LHS and RHS Exchanging Components

  string iX;					// String for integer parse
  for(i=0; i<nlhs; i++)				// Loop left side components
    {
    iX = cutInt(SLHS);				// Try and cut out an integer
    if(iX.length())				// If successful, we store
      { 					// component and get ready for next
      lhs.push_back(atoi(iX.c_str()));		//   Store component index
      cutBlksXBlks(SLHS, "+");			//   Cut + out for next one
      }
    else 					// If unsuccesful, we must quit
      {						// and issue warning if desired
      if(warn) XPerror(35, 1);
      return false;
      }
    }

  for(i=0; i<nrhs; i++)				// Loop right side components
    {
    iX = cutInt(SRHS);				// Try and cut out an integer
    if(iX.length())				// If successful, store 
      { 					// component and get ready for next
      rhs.push_back(atoi(iX.c_str()));		//   Store component index
      cutBlksXBlks(SRHS, "+");			//   Cut + sign out for next one
      }
    else 					// If unsuccesful, we must quit
      {						// and issue warning if desired
      if(warn) XPerror(36, 1);
      return false;
      }
    }
  return true;
  }
//XPfatal(10);				// and we can't go on.
// sosik


//---------------------------------------------------------------------------------
//------------------------ Read in the Process Exchange Rate ----------------------
//---------------------------------------------------------------------------------

// This will read a parameter such as    Kex_nm(0)  (1) : 600.0 - rate

bool process::getRate(const ParameterSet& pset, int idx,
                                                     double& rate, bool warn) const
  {
  rate = 0.0; 					// We don't know any rate
  ParameterSet::const_iterator item;         // A pix into parameter list
  string Nidx = "";                             // Name addition per index
  if(idx >= 0) Nidx                             // If index exists then set up
    += string("(")+Gdec(idx)+string(")");       // the parameter name to append
  string pstate, pname = string("Kex_nm")+Nidx;	// Parameter name for rate
  item = pset.seek(pname);			// Pix pset for Kex_nm(idx)
  if(item != pset.end())			// If we found the rate value
    {
    (*item).parse(pname,rate,pstate);           //   Retrieve the rate
    return true;				//   Return we succeeded
    }
  if(warn)					// If we didnt find the exchange
    XPerror(20, pname, 1); 			// rate definition, warn
  return false;					// Return we failed getting rate
  }

//---------------------------------------------------------------------------------
//------------------------ Read in the Process Spin Mappings ----------------------
//---------------------------------------------------------------------------------
 
// This will be a parameter such as    Smap(0,0)  (2) : (0)0(1)0 - Mapping

/*
bool process::getMappings(const ParameterSet& pset, int idx,
                                                     double& rate, bool warn) const

  {
  string pname = "Smap(";			// Basis for spin mapping
  string mapname, spair;			// Parameter name, value
  if(ip>=0) pname += Gdec(ip) + string(",");	// Adjust if process indexed
  int j=0, maps=1;				// Map index, flag for map search
  npair = 0;					// Zero the spin pair count
  while(maps)					// Search for spin pair maps
    {
    mapname = pname + Gdec(j) + string(")");	// This is parameter name
    item = pset.seek(mapname);			// Pix pset for Smap(i,j)
    if(item != pset.end()) 			// If map exists count it
       {					// and prep to read the next
       npair++;					// spin pair mapping
       j++;
       }
    else maps=0;				// Else exit search
    }
  if(!npair)					// If no spin pair mappings
    {						// have been found then we have
    XPerror(6, 1);				// some problems
    XPerror(4, mapname, 1);
    XPfatal(10);
    }
  link = new spin_pair[npair];			// Array for spin pairs exchanging
  maps = 1;					// Flag again to do map search
  for(j=0; j<npair; j++)			// Loop spin pair maps
    {
    mapname = pname + Gdec(j) + string(")");	// This is parameter name
    item = pset.seek(mapname);			// Pix pset for Smap(i,j)
    if(item != pset.end())
       {
       (*item).parse(mapname,spair,pstate); 	// 	Retrieve mapping as String
       spin_pair sp(spair);			// 	Spin pair from String
       link[j] = sp;				// 	Add mapping to process list
       }
    else					//  If one of the maps isn't found
       {					//  this is a problem
       XPerror(7,1);
       XPerror(4, mapname, 1);
       XPfatal(10);
       }
    }
  }
*/

//---------------------------------------------------------------------------------
//---------------------- Read/Set The Entire Exchange Process ---------------------
//---------------------------------------------------------------------------------

bool process::getXP(const ParameterSet& pset, int idx, bool warn) const
  {
string exch;
double rate;
  if(!getExch(pset, idx, exch, true))		// Try and get Exch(idx) string
    {						// If we have trouble with this
    if(warn) XPerror(31, 1); 			//   Issue warnings if desired
    return false;				//   Return that we failed
    }
  if(!getRate(pset, idx, rate, true))		// Try and get Kex_nm(idx) value
    {						// If we have trouble with this
    if(warn) XPerror(32, 1); 			//   Issue warnings if desired
    return false;				//   Return that we failed
    }

  return true;					// Return we read the exchange
  }

// sosik


bool process::setXP(const ParameterSet& pset, int idx, bool warn) const
  {
  return true;
  }


//_________________________________________________________________________________
//                    CLASS PROCESS CONSTRUCTORS AND DESTUCTORS
//_________________________________________________________________________________

//---------------------------------------------------------------------------------
//                              Simple Constructors
//---------------------------------------------------------------------------------

// process::process()				// Inlined in header
// process::process(int N_lhs, int N_rhs)	// Inlined in header


process::process(const process& PR)

        // Input                pro     : Process (this)
        //                      PR      : Another process
        // Output               void    : The process is constructed
        //                                that is identical to PR

  {
  npair = PR.npair;				// New number of pairs
  n_lhs = PR.n_lhs;				// New number of lhs comps
  n_rhs = PR.n_rhs;				// New number of rhs comps
  krate = PR.get_k();				// New exchange rate
  link = new spin_pair[npair];			// New array for spin maps
  comp_lhs = new int[n_lhs];			// New array for rhs comps
  comp_rhs = new int[n_rhs];			// New array for lhs comps
  for(int i=0; i<npair; i++)			// Copy all links (spin pairs)
    link[i] = PR.link[i];
  for(int j=0; j<n_lhs; j++)			// Copy all lhs comp. indices
    comp_lhs[j] = PR.comp_lhs[j];
  for(int k=0; k<n_rhs; k++)			// Copy all rhs comp. indices
    comp_rhs[k] = PR.comp_rhs[k];
  } 


process::process(string& PROC, double Kex, int maxcomp)

        // Input                pro     : Process (this) 
        //                      PROC	: String for process def.
	//			Kex     : Exchange rate (1/sec)
        //                      maxcomp	: Maximum possible components
        // Output               void    : The process is constructed

   {
//			First Parse The Input String 

   string pcopy = PROC;				// Copy of input string
   int N_lhs=0, N_rhs=0;			// Set none of either side
   int *Lcompind, *Rcompind;	// Array for component indices
   Lcompind = new int[maxcomp];
   Rcompind = new int[maxcomp];
   cutBlksXBlks(pcopy, "(");			// Cut initial (
   string X = cutInt(pcopy);			// Cut out the 1st integer
   if(X.length())Lcompind[N_lhs]=atoi(X.c_str());// This is 1st lhs component
   N_lhs++;
   while(pcopy[0] == '+')
     {
     cutBlksXBlks(pcopy, "+");			// Cut initial +
     Lcompind[N_lhs] 				// Cut an integer
          = atoi((cutInt(pcopy)).c_str());
     N_lhs++;					// Increment left side count
     }
   cutBlksXBlks(pcopy, "<");			// Cut initial <
   cutBlksXBlks(pcopy, "=");			// Cut initial =
   cutBlksXBlks(pcopy, ">");			// Cut initial >
   Rcompind[N_rhs]=atoi((cutInt(pcopy)).c_str());// Cut out the 1st integer
   N_rhs++;					// So far one on right
   while(pcopy[0] == '+')			// Loop the rest and see if more
     {
     cutBlksXBlks(pcopy, "+");			// Cut initial +
     Rcompind[N_rhs] = atoi(cutInt(pcopy).c_str());// Cut out an integer
     N_rhs++;					// Increment right side count
     }
//   if(pcopy.at(0,1) != ')')
//     multi_sys_error(2);

//			  Now Set The Process Up

   comp_lhs = new int[N_lhs];		// Construct left component array
   n_lhs = N_lhs ;			// Set number of left components
   comp_rhs = new int[N_rhs];		// Construct right component array
   n_rhs = N_rhs;			// Set number of right components
   int i;
   for(i=0; i<N_lhs; i++)
     comp_lhs[i] = Lcompind[i];		// Set left component indices
   for(i=0; i<N_rhs; i++)
     comp_rhs[i] = Rcompind[i];		// Set right component indices
   krate = Kex;				// Set the exchange rate

//			  No Specific Spins Exchanging Yet

   npair = 0;
   link = NULL;
   delete [] Lcompind;
   delete [] Rcompind;
   }
 

void process::intra_default(int ic1, int ic2, int nspins, double kr)

        // Input                pro     : Process (this) 
        //                      ic1	: A spin index 
        //                      ic2	: A second spin index 
	//			nspins  : Spins involved
	//			kr      : Exchange rate
        // Output               void    : The process is set for a 

  {
  n_lhs = 1;					// Set lhs spins to 1
  n_rhs = 1;					// Set rhs spins to 1
  comp_lhs = new int[1];			// Array for lhs spins
  comp_rhs = new int[1];			// Array for rhs spins
  comp_lhs[0] = ic1;				// Set lhs 1st spin to ic1
  comp_rhs[0] = ic2;				// Set rhs 1st spin to ic2
  krate = kr;					// Set exchange rate to 1
  for(int i=0; i<nspins; i++)
    {
    spin_pair sp(ic1, i, ic2, i);
    add_pair(sp);				
    }
  return;
  }

//----------------------------------------------------------------------------
//                    Construction From Parameter Set
//----------------------------------------------------------------------------

        // Input                pro     : Exchange process (this) 
        //                      pset    : A parameter set
        //                      ip 	: A process index 
        // Output               pro     : The process is constructed
	//				  from values in the parameter set
	// Note				: Non-indexed parameters are allowed 
	//				  and used if ip < 0 
 
process::process(const ParameterSet& pset, int ip, int warn)
  {
  ParameterSet::const_iterator item;		// A pix into parameter list
  string pstate;				// For parameter statement
  string Exval;					// For String parameter value
  double Kex;					// For exchange rate value

  string idx = ""; 				// Parameter name end will
  if(ip>=0) idx=string("(")+Gdec(ip)+string(")");// exist if indexed

//---------------------------------------------------------------------------------
//------------------------- Read in the Process Definition ------------------------
//---------------------------------------------------------------------------------

// This will be a parameter such as    Exch(0)  (2) : (0<=>1+2) - Exchange sheme

// The string value defines the components (spin systems) involved in the
// exchange process.  This only gets the process definition, it doesn't set
// it up at all!

   if(!getExch(pset, ip, Exval, true))		// Try and get Exch(ip) string
     {						// If we have trouble with this
     if(warn)					//   Issue warnings if desired
       {
       if(warn>1) XPfatal(31);			//   Cannot find exchange process
       else       XPerror(31);
       }
     return;					// Return with empty process
     }

//---------------------------------------------------------------------------------
//------------------------ Read in the Process Exchange Rate ----------------------
//---------------------------------------------------------------------------------

// This will be a parameter such as    Kex_nm(0)  (1) : 600.0 - rate
// This sets process Kex value!

   if(!getRate(pset, ip, Kex, true))		// Try and get Kex_nm(ip) value
     {						// If we have trouble with this
     if(warn)					//   Issue warnings if desired
       {
       if(warn>1) XPfatal(32);			//   Cannot find exchange rate
       else       XPerror(32);
       }
     return;					// Return with empty process
     }
   krate = Kex;

//---------------------------------------------------------------------------------
//-------------- Parse The Components Involved In The Exchange Process ------------
//---------------------------------------------------------------------------------
 
// At the moment, we have String Exval from parameter Exch, e.g. (0<=>1+2).
// We need to parse out the integer values on each side of <=> but between ().
// Also we must keep track of which are on the left and which are on the right.

//		First Get Number of LHS and RHS Components in Exchange

vector<int> RHS;
vector<int> LHS;
if(!parseExch(Exval, LHS, RHS, true))
XPfatal(10);
n_lhs = int(LHS.size());
n_rhs = int(RHS.size());
comp_lhs = new int[n_lhs];			// Allocate space for left
comp_rhs = new int[n_rhs];			// Allocate space for right
int i;
for(i=0; i<n_lhs; i++)
comp_lhs[i] = LHS[i];
for(i=0; i<n_rhs; i++)
comp_rhs[i] = RHS[i];

//---------------------------------------------------------------------------------
//------------------------ Read in the Process Spin Mappings ----------------------
//---------------------------------------------------------------------------------
 
// This will be a parameter such as    Smap(0,0)  (2) : (0)0(1)0 - Mapping

   string pname = "Smap(";                      // Basis for spin mapping
   string mapname, spair;			// Parameter name, value
   if(ip>=0) pname += Gdec(ip) + string(",");	// Adjust if process indexed
   int j=0, maps=1;				// Map index, flag for map search
   npair = 0;					// Zero the spin pair count
   while(maps)					// Search for spin pair maps
     {
     mapname = pname + Gdec(j) + string(")");	// This is parameter name
     item = pset.seek(mapname);			// Pix pset for Smap(i,j)
     if(item != pset.end()) 			// If map exists count it
       {					// and prep to read the next
       npair++;					// spin pair mapping
       j++;
       }
     else maps=0;				// Else exit search
     }
   if(!npair)					// If no spin pair mappings
     {						// have been found then we have
     XPerror(6, 1);				// some problems
     XPerror(4, mapname, 1);
     XPfatal(10);
     }
   link = new spin_pair[npair];			// Array for spin pairs exchanging
   maps = 1;					// Flag again to do map search
   for(j=0; j<npair; j++)			// Loop spin pair maps
     {
     mapname = pname + Gdec(j) + string(")");	// This is parameter name
     item = pset.seek(mapname);			// Pix pset for Smap(i,j)
     if(item != pset.end())
       {
       (*item).parse(mapname,spair,pstate); 	// 	Retrieve mapping as String
       spin_pair sp(spair);			// 	Spin pair from String
       link[j] = sp;				// 	Add mapping to process list
       }
     else					//  If one of the maps isn't found
       {					//  this is a problem
       XPerror(7,1);
       XPerror(4, mapname, 1);
       XPfatal(10);
       }
     }
   }


//----------------------------------------------------------------------------
//                         Assignment and Destruction
//----------------------------------------------------------------------------

process& process::operator=(const process& PR)

        // Input                pro     : Process (this)
        //                      PR      : Another process
        // Output               void    : The process is constructed
        //                                that is identical to PR

  {
//	    First Delete Any Existing Process Information

  if(link) delete [] link;			// Now remove existing links
  if(comp_rhs) delete [] comp_rhs;		// Remove existing rhs comps.
  if(comp_lhs) delete [] comp_lhs;		// Remove existing lhs comps.

//	        Now Copy All Process Information From Input PR

  npair = PR.npair;				// New number of pairs
  n_lhs = PR.n_lhs;				// New number of lhs comps
  n_rhs = PR.n_rhs;				// New number of rhs comps
  krate = PR.krate;				// New exchange rate
  link = new spin_pair[npair];			// New array for spin maps
  comp_lhs = new int[n_lhs];			// New array for rhs comps
  comp_rhs = new int[n_rhs];			// New array for lhs comps
  for(int i=0; i<npair; i++)			// Copy all links (spin pairs)
    link[i] = PR.link[i];
  for(int j=0; j<n_lhs; j++)			// Copy all lhs comp. indices
    comp_lhs[j] = PR.comp_lhs[j];
  for(int k=0; k<n_rhs; k++)			// Copy all rhs comp. indices
    comp_rhs[k] = PR.comp_rhs[k];
  return *this;
  }
  

process::~process() 
  {
  delete [] link;				// Delete spin mappings
  delete [] comp_lhs;				// Delete left hand indices
  delete [] comp_rhs;				// Delete right hand indices
  }

 
//_________________________________________________________________________________
// C                  Class Process Component Index Access
//_________________________________________________________________________________
 
 
int process::lhsindex(int comp)
 
        // Input                pro     : A process
        //                      comp    : A lhs component
        // Output               ic      : Index of component

  {
  if(comp<0 || comp>=n_lhs) 			// Check that the component
    {						// actually exists
    XPerror(2, Gdec(comp));
    exit(-1);
    }
  return comp_lhs[comp];			// Return component index
  }
 
 
int process::rhsindex(int comp)
 
        // Input                pro     : A process
        //                      comp    : A rhs component
        // Output               ic      : Index of rhs component
 
  {
  if(comp<0 || comp>=n_rhs) 			// Check that the component
    {						// actually exists
    XPerror(3, Gdec(comp));
    exit(-1);
    }
  return comp_rhs[comp];			// Return component index
  }

//_________________________________________________________________________________
// D                        CLASS PROCESS SPIN QUERIES
//_________________________________________________________________________________
 

int process::mixes(int ic1, int ic2)
 
        // Input                pro     : A process
	//			ic1     : Spin index 1
	//			ic2     : Spin index 2
        // Output		TF	: True if both spins are
	//				  involved in the process
 
  {
  int left=0, right=0 ;
  int i;
  for(i=0; i<n_lhs && !left; i++)		// Loop all on the lhs
    if((comp_lhs[i]==ic1) || (comp_lhs[i]==ic2) )
      left = 1;
  for(i=0; i<n_rhs && !right; i++)		// Loop all on the rhs
    if((comp_rhs[i]==ic2) || (comp_rhs[i]==ic1) )
       right = 1;
  if(left && right) return 1;
  return 0;
  }


int process::involves(int ic1, int lr)
 
        // Input                pro     : A process
	//			ic1     : Spin index 1
	//			lr      : Flag to check left and/or right
	//					0 = check left & right (def)
	//				       >0 = check right only
	//				       <0 = check left only 
        // Output		TF	: True if spin is involved in
	//				  the right and/or left process

  {
  for(int i=0; i<n_lhs && lr<=0; i++)
    if(comp_lhs[i] == ic1) return 1;
  for (int j=0; j<n_rhs && lr>=0; j++)
    if(comp_rhs[j] == ic1) return 1;
  return 0;
  }


//________________________________________________________________________________
// E                       CLASS PROCESS SPIN PAIR ACCESS
//________________________________________________________________________________


int process::pairs() const { return npair; }

        // Input                pro     : A process (this)
        // Output               npair   : Number of spin pairs defined
        //                                in the process
 

spin_pair process::get_pair(int ip) const { return link[ip]; }

        // Input                pro     : A process (this)
        //                      ip      : A spin pair index
        // Output               SP      : Spin pair (copy) is returned


void process::add_pair(spin_pair sp)

        // Input                pro     : A process (this)
        //                      SP      : A spin pair
        // Output               void    : Spin pair is set as exchanging
        //                                in the process
 
  {
  spin_pair* temp;				// Construct a new array
  temp = new spin_pair[npair+1];		// of spin pairs 
  for(int i=0; i<npair; i++)			// Copy all original spin
    temp[i] = link[i];				// pairs into the new array
  temp[npair] = sp;				// Add in the new spin pair
  if(npair) delete [] link;			// Delete the original array
  npair++;					// Update spin pair count
  link = temp;					// Now point to new array
  return;
  } 


int process::mapped(int c1, int s1, int c2, int s2) const

        // Input                pro     : A process (this)
        //                      int     : A spin pair index
        // Output               TF      : True if spins are mapped
        //                                false if not

  {
  spin_pair SP;					// An empty spin pair 
  for(int i=0; i<npair; i++)			// Loop over current spin pairs
    {
    SP = link[i];				// Get this spin pair
    if((c1==SP.sub1) && (s1==SP.sp1)		// Check if spin mapping
       && (c2==SP.sub2) && (s2==SP.sp2) )
         return 1;
    if((c1==SP.sub2) && (s1==SP.sp2)
       && (c2==SP.sub1) && (s2==SP.sp1) )
    return 1;
    }
  return 0;
  }


//__________________________________________________________________________________
// F                          CLASS PROCESS CONNECTIONS
//__________________________________________________________________________________


void process::mapping(const string& spair)

        // Input                pro     : Process (this) 
        //                      spair	: String for spin pair def.
        // Output               void    : The process is modified
	//				  to contain a spin pairing as
	//				  suggested in the input string
	// Note				: It is assumed the input string
	//				  has the form
	//
	//					(C1)S1(C2)S2
	//
	//				  where C1,S1,C2,S2 are integers
	// Note				: No checks are preformed to insure
	//				  that any of these indices are valid
 
   {
   spin_pair sp(spair);				// Generate spin pair from String
   add_pair(sp);				// Add this mapping to process list
   }


// ____________________________________________________________________________
//                       EXCHANGE PROCESS INPUT FUNCTIONS
// ____________________________________________________________________________

        // Input                pro	: An exchange process (this)
        //                      filename: Input filename
        //                      pset    : Input parameter set
	//			idx     : Process index
        //                      warn    : Warning output level
        //                                      0 = no warnings
        //                                      1 = warnings
        //                                     >1 = fatal warnings
        // Output               TF      : Process is filled with
        //                                parameters read from file
        //                                TRUE if read is successful
        // Note                         : The file should be an ASCII file
        //                                containing recognized parameters

bool process::read(const string &filename, int idx, int warn)
   {
   ParameterSet pset;                   // A GAMMA parameter set
   if(!pset.read(filename, warn?1:0))   // Try and read in pset
    {                                   // If we cannot read the file then
    if(warn)                            // we'll issue warnings as desired
      {
      XPerror(1, filename); 		//      Problems with file
      if(warn>1) XPfatal(14);		//      Can't read from input file
      else       XPerror(14);		//      Fatal if warn big enough
      }
    return false;                       // Well flag we didn't read & exit
    }
  return read(pset, idx, warn);		// Fill up spin_sys with parameters
  }

bool process::read(const ParameterSet& pset, int idx, int warn)
  {
  bool TF = setXP(pset, idx, warn?true:false);	// Use overload to read
  if(!TF)                                       // If setXP didn't handle
    {                                           // the process read from pset
    if(warn)                                    // then we'll issue warnings
      {                                         // & maybe even stop program
      XPerror(8, 1);				//   Problems with pset
      if(warn>1) XPfatal(30);			//   Can't read from pset
      else       XPerror(30, 1);		//   Fatal if warn big enough
      }
    }
  return TF;
  }


//__________________________________________________________________________________
// G                            CLASS PROCESS I/O
//__________________________________________________________________________________

 
ostream& process::print(ostream& ostr, int full) const

        // Input                pro   : A process
        //                      ostr  : Output stream
        //                      full  : Print amount flag
        //                                 0 = don't print spin mappings(def)
        //                                !0 = print individual spin mappings
        // Output               ostr  : The output stream  is returned
        //                              with the process added

    {                                                              
    int i;
    ostr << "\nRate: "  << krate << "/sec";	// Output the exchange rate
    ostr <<  " lhs = {";			// Output left hand side
    for(i=0; i<n_lhs; i++)			// components in process
      {
      ostr << comp_lhs[i];
      if(i+1 != n_lhs) ostr << ", ";
      }
    ostr << "}  rhs = {";			// Output right hand size
    for(i=0; i<n_rhs; i++)			// components in process
      {
      ostr << comp_rhs[i];
      if(i+1 != n_rhs) ostr << ", ";
      }
    ostr << "}";
    if(full)					// If full flag is set then
      {						// output all the spin mappings
      for(i=0; i<npair; i++)			// Loop over all spin pairs
        ostr << "\n\t" << link[i];		// and print spin mapping
      }
    return ostr;
    }


ostream& operator<< (ostream& ostr, const process& pro)
 
        // Input                ostr     : An output stream
        // Input                pro      : A process
        // Output               ostr     : The output stream modified by
        //                                 the defined process
 
  { return pro.print(ostr); } 			// Just use the member function


#endif
