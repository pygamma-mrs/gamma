/* ExProcess.cc *************************************************-*-c++-*-
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      Non-Mutual Exchange Process             Implementation		**
**                                                                      **
**      Copyright (c) 2001                                              **
**      Dr. Scott A. Smith                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** This class defines a single single, specific, non-mutual exchange    **
** process. This includes which components (dynamic systems) in a       **
** multi_sys system are in exchange, the exchange rate, and the         **
** spin <--> spin mappings that exist in the exchange.                  **
**                                                                      **
*************************************************************************/

#ifndef _ExProc_cc_			// Is the file already included?
#  define _ExProc_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#if defined(_MSC_VER)			// If we are using MSVC++
 #pragma warning (disable : 4786)       //   Kill STL namelength warnings
#endif

#include <MultiSys/ExProcess.h>		// Include header file
#include <cstdlib>                      // Include functions atoi, atof, etc.
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Basics/StringCut.h>		// Include string cutting 
#include <Basics/Gutils.h>		// Include GAMMA error messages
#include <string>			// Know about libstdc++ strings
#include <fstream>			// Know about file streams
#include <cmath>			// Know about fabs function
#include <list>				// Include libstdc++ STL lists
#include <vector>			// Include libstdc++ STL vectors
#include <iostream>                     // Include input output streams

using std::string;			// Using libstdc++ strings
using std::list;			// Using libstdc++ lists
using std::vector;			// Using libstdc++ vectors
using std::ostream;			// Using libstdc++ output streams

//_________________________________________________________________________________
// i                       CLASS PROCESS ERROR HANDLING
//_________________________________________________________________________________

	// Input			PRO	: A process (this)
        // 				eidx	: Error flag
        //				noret   : Flag for return (0=return)
        //                              pname   : String included in error 
        // Output			none 	: Output process error message


void ExchProc::XPerror(int eidx, int noret) const
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
    case 36: msg=string("Cannot Parse Left Side Exchange Comps");break; // (36)
    case 37: msg=string("Cant Set Up From Input Parameter Set"); break; // (37)
    case 38: msg=string("Problems Determining Exchange Comps");  break; // (38)
    case 39: msg=string("No Spin Maps Found in Parameter Set");  break; // (39)
    case 66: msg=string("Can't Find Any Spin Mappings");         break; // (66)
    }
  GAMMAerror(hdr, msg, noret);
  }



/*
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


void ExchProc::XPerror(int eidx, const string& pname, int noret) const
 
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
    default: GAMMAerror(hdr, eidx, pname, noret); return;  break;
    case 20: msg=string("Cannot Read Parameter ") + pname; break;      	// (20)
    case 22: msg=string("Accessed LHS Component ") + pname
                +string(" Out of Range");                  break;	// (22)
    case 23: msg=string("Accessed RHS Component ") + pname
                +string(" Out of Range");                  break;	// (23)
/*
    case 4:								// (4)
      cout << "Attempted To Obtain Parameter " << pname;
      break;
*/
    }
  GAMMAerror(hdr, msg, noret);
  }
 
 
volatile void ExchProc::XPfatal(int eidx) const
  {
  XPerror(eidx, 1);			// Output the error message
  if(eidx) XPerror(0,1);		// State that its fatal
  GAMMAfatal();				// Clean exit from program
  }

volatile void ExchProc::XPfatal(int eidx, const string& pname) const
  {
  XPerror(eidx, pname, 1);		// Output the error message
  if(eidx) XPerror(0,1);                // State that its fatal
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

/* This will be a parameter such as    Exch(0)  (2) : (0<=>1+2) - Exchange scheme

   The string value defines the components (subsystems) involved in the exchange
   process. The function getExch will just get the process definition, it doesn't
   bother to parse it. The fucntion parseExch will take the string from getExch
   and split it up int left hand side (L.H.S.) and right hand side (R.H.S)
   components. Function getComps just oversees the process of getting the
   components so that it looks like it is done in one step.

   The string value of Exch will be something akin to (0<=>1+2). We need to then
   to parse out the integer values on each side of <=>, but lying between ().
   Also we must keep track of which are on the left and which are on the right
   and how many of these there are (separated by + signs).                      */

bool ExchProc::getExch(const ParameterSet& pset, int idx,
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

bool ExchProc::parseExch(string& Exval, 
                               vector<int>& lhs, vector<int>& rhs, bool warn) const
  {
//	    First Get Raw Strings For LHS and RHS Components in Exchange

  int len = Exval.length();			// Length of input string
  int nc  = Exval.find("<=>");			// Number of chars before <=>
  string SLHS = Exval.substr(1, nc-1);		// Get string between ( & <=>
  string SRHS = Exval.substr(nc+3, len-6);	// Get string between <=> & )
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

//         At This Point LHS & RHS Should Be Strings Akin To #+#+#.....
//            We Need To Parse Out The Individual #'s (Components)
//                       We'll Do The L.H.S One First

  string iX;					// String for integer parse
  lhs.clear();					// Insure no components yet
  len = SLHS.length();				// Length of raw L.H.S. string
  while(len > 0)				// While the string has length
    {
    iX = cutInt(SLHS);				//   Cut out an integer
    lhs.push_back(atoi(iX.c_str()));		//   Store component index
    if(len > 1)	cutBlksXBlks(SLHS, "+");	//   Cut + out for next one
    len = SLHS.length();			//   New length
    }

//                Now Parse The Crude R.H.S Sting Into Components

  rhs.clear();					// Insure no components yet
  len = SRHS.length();				// Length of raw R.H.S. string
  while(len > 0)				// While the string has length
    {
    iX = cutInt(SRHS);				//   Cut out an integer
    rhs.push_back(atoi(iX.c_str()));		//   Store component index
    if(len > 1)	cutBlksXBlks(SRHS, "+");	//   Cut + out for next one
    len = SRHS.length();			//   New length
    }
  return true;
  }

bool ExchProc::getComps(const ParameterSet& pset, int idx,
              std::vector<int>& lhs, std::vector<int>& rhs, bool warn) const
  {
  string Exval;
  if(!getExch(pset, idx, Exval, warn))		// Try and get Exch(idx) string
    {						// If we have trouble with this
    if(warn) XPerror(31, 1); 			//   Issue warnings if desired
    return false;				//   Return that we failed
    }
  if(!parseExch(Exval, lhs, rhs, warn))		// Try for exchange components
    {						// If we have trouble with this
    if(warn) XPerror(38, 1); 			//   Issue warnings if desired
    return false;				//   Return that we failed
    }
  return true;					// Yes, we succeeded
  }

//---------------------------------------------------------------------------------
//------------------------ Read in the Process Exchange Rate ----------------------
//---------------------------------------------------------------------------------

// This will read a parameter such as    Kex_nm(0)  (1) : 600.0 - rate

bool ExchProc::getRate(const ParameterSet& pset, int idx,
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

bool ExchProc::getMappings(const ParameterSet& pset, int idx,
                                           vector<SpinMap>& smaps, bool warn) const
  {
  smaps.clear();                                // We don't know any spin maps
  SpinMap sm;					// Single spin map
  int j=0;                                      // Exchange process index
  while(sm.read(pset, idx, j, false))		// Try & read mapping (idx, j)
    {
    smaps.push_back(sm);			// If we can, store it in maps
    j++;                                        // Increment to next spin map
    }
  if(!j)                                        // If we didn't read any
    { if(warn) XPerror(66,1); return false; }	// Warn and return our failure
  return true;                                  // Return we read em OK
  }

//---------------------------------------------------------------------------------
//---------------------- Read/Set The Entire Exchange Process ---------------------
//---------------------------------------------------------------------------------

bool ExchProc::getXP(const ParameterSet& pset,
   double& rate, vector<int>& lhsc, vector<int>& rhsc, vector<SpinMap>& smaps,
                                                          int idx, bool warn) const
  {
  if(!getComps(pset, idx, lhsc, rhsc, warn))	// Try for exchange components
    {						// If we have trouble with this
    if(warn) XPerror(31, 1); 			//   Issue warnings if desired
    return false;				//   Return that we failed
    }
  if(!getRate(pset, idx, rate, true))		// Try and get Kex_nm(idx) value
    {						// If we have trouble with this
    if(warn) XPerror(32, 1); 			//   Issue warnings if desired
    return false;				//   Return that we failed
    }
  if(!getMappings(pset, idx, smaps, warn))	// Try and get Smap(idx,j)
    {						// If we have trouble with this
    if(warn) XPerror(39, 1); 			//   Issue warnings if desired
    return false;				//   Return that we failed
    }
  return true;					// Return we read the exchange
  }

bool ExchProc::setXP(const ParameterSet& pset, int idx, bool warn)
  {
  double          k;				// For exchange rate
  vector<int>     lhsc;				// For L.H.S. components
  vector<int>     rhsc;				// For R.H.S. components
  vector<SpinMap> smaps;			// For spin mappings
  if(getXP(pset,k,lhsc,rhsc,smaps,idx,warn))	// Try and read process
    {						// If we are successful
    KRate    = fabs(k);				//   Set the exchange rate
    LHSComps = lhsc;				//   Set L.H.S. components
    RHSComps = rhsc;				//   Set R.H.S. components
    SpinMaps = smaps;				//   Set spin mappings
    return true;				//   Return successful
    }
  if(warn) XPerror(37, 1); 			// We failed, issue warnings
  return false;					// Return our miserable failings
  }

//_________________________________________________________________________________
// iii             CLASS EXCHANGE PROCESS CHECKING FUNCTIONS
//_________________________________________________________________________________

/* These functions check the boundaries of the exchange components and insure 
   process integrity. Process integrity amounts to insuring that the number of
   species on the left hand side of the exchange is the same as that on the
   right hand side of the exchange.  */

bool ExchProc::CheckLHS(int comp, bool warn) const
  {
  if(comp<0 || comp>=(int)LHSComps.size()) 	// Check that the component
    {						// actually exists
    if(warn) XPerror(10, 1);
    return false;
    }
  return true;
  }

bool ExchProc::CheckRHS(int comp, bool warn) const
  {
  if(comp<0 || comp>=(int)RHSComps.size()) 	// Check that the component
    {						// actually exists
    if(warn) XPerror(10, 1);
    return false;
    }
  return true;
  }

//_________________________________________________________________________________
// A                  CLASS PROCESS CONSTRUCTORS AND DESTUCTORS
//_________________________________________________________________________________

//---------------------------------------------------------------------------------
//                              Simple Constructors
//---------------------------------------------------------------------------------

ExchProc::ExchProc() { KRate = 0.0; } 			// No exchange rate

ExchProc::ExchProc(const ExchProc& PR)
  {
  KRate    = PR.KRate;					// Copy exchange rate
  LHSComps = PR.LHSComps;				// Copy array for rhs comps
  RHSComps = PR.RHSComps;				// Copy array for lhs comps
  SpinMaps = PR.SpinMaps;				// Copy array for spin maps
  } 

//---------------------------------------------------------------------------------
//                              Specific Constructors
//---------------------------------------------------------------------------------

//          Here the arguments must define the exchange process explicitly

        // Input                pro     : Process (this) 
        //                      PROC	: String for process def.
	//			Kex     : Exchange rate (1/sec)
        //                      maxcomp	: Maximum possible components
        // Output               void    : The process is constructed

ExchProc::ExchProc(const string& PROC, double Kex, int maxcomp)
  {
LHSComps.clear();				// Zero left component array
  KRate = fabs(Kex);				// Set the exchange rate
  string pcopy = PROC;				// Copy of input string
  int N_lhs=0, N_rhs=0;				// Set none of either side
  int* Lcompind;				// Array for component indices
  int* Rcompind;
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

   LHSComps.clear();			// Zero left component array
   RHSComps.clear();			// Zero right component array
   int i;
   for(i=0; i<N_lhs; i++)
     LHSComps.push_back(Lcompind[i]);	// Set left component indices
   for(i=0; i<N_rhs; i++)
     RHSComps.push_back(Rcompind[i]);	// Set right component indices

//			  No Specific Spins Exchanging Yet

   SpinMaps.clear();
   delete [] Lcompind;
   delete [] Rcompind;
   }
 

ExchProc::ExchProc(int N_lhs, int N_rhs)

        // Input                pro     : A process (this)
        //                      N_lhs   : Number left side components
        //                      N_rhs   : Number right side components
        // Output               void    : A NULL process created
        // Note                         : No connectivities are set

  {
  KRate = 0.0;                          // Set exchange rate to zero
  LHSComps = vector<int>(N_lhs);	// Allocate space for left components
  RHSComps = vector<int>(N_rhs);	// Allocate space for left components
  }

void ExchProc::intra_default(int ic1, int ic2, int nspins, double kr)

        // Input                pro     : Process (this) 
        //                      ic1	: A spin index 
        //                      ic2	: A second spin index 
	//			nspins  : Spins involved
	//			kr      : Exchange rate
        // Output               void    : The process is set for a 

  {
  LHSComps.clear();				// Zero array for lhs spins
  RHSComps.clear();				// Array for rhs spins
  LHSComps.push_back(ic1);			// Set lhs 1st spin to ic1
  RHSComps.push_back(ic2);			// Set rhs 1st spin to ic2
  KRate = kr;					// Set exchange rate to 1
  for(int i=0; i<nspins; i++)
    {
    SpinMap sp(ic1, i, ic2, i);
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
 
ExchProc::ExchProc(const ParameterSet& pset, int ip, int warn)
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
   KRate = Kex;

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
int n_lhs = int(LHS.size());
int n_rhs = int(RHS.size());
LHSComps.clear();
RHSComps.clear();
int i;
for(i=0; i<n_lhs; i++)
LHSComps.push_back(LHS[i]);
for(i=0; i<n_rhs; i++)
RHSComps.push_back(RHS[i]);

//---------------------------------------------------------------------------------
//------------------------ Read in the Process Spin Mappings ----------------------
//---------------------------------------------------------------------------------
 
// This will be a parameter such as    Smap(0,0)  (2) : (0)0(1)0 - Mapping

   string pname = "Smap(";                      // Basis for spin mapping
   string mapname, spair;			// Parameter name, value
   if(ip>=0) pname += Gdec(ip) + string(",");	// Adjust if process indexed
   int j=0, maps=1;				// Map index, flag for map search
int npair = 0;					// Zero the spin pair count
SpinMaps.clear();
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
     XPerror(66, 1);				// some problems
     XPerror(4, mapname, 1);
     XPfatal(10);
     }
   maps = 1;					// Flag again to do map search
   for(j=0; j<npair; j++)			// Loop spin pair maps
     {
     mapname = pname + Gdec(j) + string(")");	// This is parameter name
     item = pset.seek(mapname);			// Pix pset for Smap(i,j)
     if(item != pset.end())
       {
       (*item).parse(mapname,spair,pstate); 	//  Retrieve mapping as String
       SpinMap sp(spair);			//  Spin pair from String
       SpinMaps.push_back(sp);			//  Add mapping to process list
       }
     else					//  If one of the maps isn't
       {					//  found this is a problem
       XPerror(7,1);
       XPerror(4, mapname, 1);
       XPfatal(10);
       }
     }
   }

//---------------------------------------------------------------------------------
//                         Assignment and Destruction
//---------------------------------------------------------------------------------

ExchProc& ExchProc::operator=(const ExchProc& PR)
  {
  KRate    = PR.KRate;				// New exchange rate
  SpinMaps = PR.SpinMaps;			// Copy all spin maps
  LHSComps = PR.LHSComps;			// Copy lhs components
  RHSComps = PR.RHSComps;			// New array for lhs comps
  return *this;
  } 

ExchProc::~ExchProc() { }

//________________________________________________________________________________
// B                   Class Exchange Process Access Functions
//________________________________________________________________________________

//--------------------------------------------------------------------------------
//                                 Exchange Rate
//--------------------------------------------------------------------------------

        // Input                pro     : An exchange process (this)
        //                      k       : An exchange rate (1/sec)
        // Output               void    : Exchange rate is set to k
        //                   or double  : Exchange rate is returned

double ExchProc::Kex() const   { return KRate; }
void   ExchProc::Kex(double k) { KRate = k;    }
 
//--------------------------------------------------------------------------------
//                    Class Process Component Index Access
//--------------------------------------------------------------------------------
 
        // Input                pro     : An exchange process (this)
        //                      comp    : A L.H.S. component
        // Output               ic      : Index of component
 
int ExchProc::LHSComp(int comp) const
  {
  if(!CheckLHS(comp))				// If component out of bounds
    XPfatal(22, Gdec(comp));			// quit with error messages
  return LHSComps[comp];			// Return component index
  }
 
int ExchProc::RHSComp(int comp) const
  {
  if(!CheckRHS(comp))				// If component out of bounds
    XPfatal(23, Gdec(comp));			// quit with error messages
  return RHSComps[comp];			// Return component index
  }

int ExchProc::NCompsLHS() const { return LHSComps.size(); }
int ExchProc::NCompsRHS() const { return RHSComps.size(); }

//_________________________________________________________________________________
// D                   CLASS PROCESS COMPONENT QUERIES
//_________________________________________________________________________________
 
        // Input                pro     : A process
	//			ic1     : Component index 1
	//			ic2     : Component index 2
        // Output		TF	: True if components are
	//				  involved in the exchange process

bool ExchProc::mixes(int comp1, int comp2) const
  {
  if(CompInLHS(comp1) && CompInRHS(comp2))	// If comp1 in LHS & comp2 in RHS
    { return mapped(comp1, comp2); } 		// see if any spins exchange
  if(CompInLHS(comp2) && CompInRHS(comp1))	// If comp2 in LHS & comp1 in RHS
    { return mapped(comp1, comp2); } 		// see if any spins exchange
  return false;
  }

bool ExchProc::CompInLHS(int comp) const
  {
  int nl = LHSComps.size();			// No. Components in LHS
  for(int i=0; i<nl; i++)			// Loop over LHS components
    if(LHSComps[i] == comp) return true;	// Return true if one of
  return false;					// them is comp, false else
  }

bool ExchProc::CompInRHS(int comp) const
  {
  int nr = RHSComps.size();			// No. Components in RHS
  for(int i=0; i<nr; i++)			// Loop over RHS components
    if(RHSComps[i] == comp) return true;	// Return true if one of
  return false;					// them is comp, else false
  }

bool ExchProc::involves(int comp, int lr) const
  {
  if(lr > 0) return CompInRHS(comp);		// Return if comp in RHS
  if(lr < 0) return CompInLHS(comp);		// Return if comp in LHS
  if(CompInRHS(comp) || CompInLHS(comp))	// Return if comp in exchange
    return true;
  return false;					// Comp is not in the exchange
  }


//________________________________________________________________________________
// E                       CLASS PROCESS SPIN MAP ACCESS
//________________________________________________________________________________

int  ExchProc::NSpinMaps()  const { return SpinMaps.size(); }

bool ExchProc::SMap(int comp1, int spin1, int& comp2, int& spin2) const
  {
  SpinMap SP;					// An empty spin map
  int nmaps = SpinMaps.size();			// # spin maps in this exchange
  for(int i=0; i<nmaps; i++)			// Loop over current spin maps
    {
    SP = SpinMaps[i];				//   Get spin map
    if(comp1 == SP.sub1)			//   If comp1 is 1st comp of map
      if(spin1 == SP.sp1)			//     If spin1 is in the map
        {					//     we set the component
        comp2 = SP.sub2;			//     and spin its n exchange
        spin2   = SP.sp2;			//     with & return we found it
        return true;
        }
    }
  comp2 = -1;					// Don't know partner comp
  spin2 = -1;					// Don't know partner spin
  return false;  				// No map for comp1 spin1
  }

const SpinMap& ExchProc::SMap(int ip) const { return SpinMaps[ip]; }

        // Input                pro     : A process (this)
        //                      ip      : A spin pair index
        // Output               SP      : Spin pair (copy) is returned



        // Input                pro     : A process (this)
        //                      SP      : A spin pair
        // Output               void    : Spin pair is set as exchanging
        //                                in the process
 
void ExchProc::add_pair(SpinMap sp)
  { SpinMaps.push_back(sp); }				// Add new spin pair


bool ExchProc::mapped(int comp1, int spin1, int comp2, int spin2) const
  {
  SpinMap SP;					// An empty spin map
  int nmaps = SpinMaps.size();			// # spin maps in this exchange
  for(int i=0; i<nmaps; i++)			// Loop over current spin maps
    {
    SP = SpinMaps[i];				//   Get spin map
    if(comp1 == SP.sub1)			//   If comp1 is 1st comp of map
      {						//   see if it maps spin1 
      if(spin1 == SP.sp1)			//     If spin1 is in the map
        {					//     check if it is mapped to
        if(comp2 == SP.sub2 && spin2 == SP.sp2)	//     comp2 spin2.  Return true
          return true;				//     if this is the case
        }
      }
    if(comp2 == SP.sub1)			//   If comp2 is 1st comp. of map
      {						//   see if it maps spin2 
      if(spin2 == SP.sp1)			//     If spin2 is in the map
        {					//     check if it is mapped to
        if(comp1 == SP.sub2 && spin1 == SP.sp2)	//     comp1 spin1.  Return true
          return true;				//     if this is the case
        }
      }
    }
  return false;
  }

bool ExchProc::mapped(int comp1, int comp2) const
  {
  int nmaps = SpinMaps.size();			// # spin maps in this exchange
  SpinMap SP;					// An empty spin map
  for(int i=0; i<nmaps; i++)			// Loop over current spin maps
    {
    SP = SpinMaps[i];				//   Get spin map
    if(comp1 == SP.sub1 && comp2 == SP.sub2)	//   If comp1 & comp2 in spin map
      return true;				//   some part of them exchanges
    if(comp2 == SP.sub1 && comp1 == SP.sub2)	//   If comp2 & comp1 in spin map
      return true;				//   some part of them exchanges
    }
  return false;
  }


//__________________________________________________________________________________
// F                          CLASS PROCESS CONNECTIONS
//__________________________________________________________________________________


void ExchProc::mapping(const string& spair)

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
   SpinMap sp(spair);				// Generate spin pair from String
   add_pair(sp);				// Add this mapping to process list
   }


//__________________________________________________________________________________
//                       EXCHANGE PROCESS INPUT FUNCTIONS
//__________________________________________________________________________________

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

bool ExchProc::read(const string &filename, int idx, int warn)
   {
   ParameterSet pset;                   // A GAMMA parameter set
   if(!pset.read(filename, warn?1:0))   // Try and read in pset
    {                                   // If we cannot read the file then
    if(warn)                            // we'll issue warnings as desired
      {
      XPerror(1, filename); 		//      Problems with file
      if(warn>1) XPfatal(14);		//      Can't read from input file
      else       XPerror(14,1);		//      Fatal if warn big enough
      }
    return false;                       // Well flag we didn't read & exit
    }
  return read(pset, idx, warn);		// Fill up spin_sys with parameters
  }

bool ExchProc::read(const ParameterSet& pset, int idx, int warn)
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
// G                            EXCHANGE PROCESS OUTPUT
//__________________________________________________________________________________

//----------------------------------------------------------------------------------
//              Functions To Modularize Exchange Process Output
//----------------------------------------------------------------------------------

/* These next two function return a string to indicate alphabetically which 
   components are involved in the exchange process. The string returned will look
   something like
                                A + D + E + G

   where the first component has label A, the second label B, etc.                */

char ExchProc::Label(int i)
  {
  char Caps[26] = { 'A', 'B', 'C', 'D', 'E', 'F',
                    'G', 'H', 'I', 'J', 'K', 'L',
                    'M', 'N', 'O', 'P', 'Q', 'R',
                    'S', 'T', 'U', 'V', 'W', 'X',
                    'Y', 'Z' };
  return Caps[i%26];
  }

string ExchProc::LHSStr() const
  {
  int nl = NCompsLHS();                 	// No. left components
  string lhsstr("");				// Left hand side string
  int i;
  for(i=0; i<nl-1; i++)
    { 
    lhsstr += Label(LHSComps[i]);
    lhsstr += " + ";
    }
  lhsstr += Label(LHSComps[i]);
  return lhsstr;
  }

string ExchProc::RHSStr() const
  {
  int nr = NCompsRHS();  	               	// No. left components
  string rhsstr("");				// Left hand side string
  int i;
  for(i=0; i<nr-1; i++)
    { 
    rhsstr += Label(RHSComps[i]);
    rhsstr += " + ";
    }
  rhsstr += Label(RHSComps[i]);
  return rhsstr;
  }

/* The next function returns a vector of strings that indicate the spin mappings
   that are active in the exchange process. Each string will involve one spin
   of a LHS component exchanging with one spin of the RHS component. Each string
   will be of the same length and appear something like

                                   A 0 <---> B 1

   where here A is the LHS component and 0 its spin, while B is a RHS component
   with 1 being its spin.                                                         */


vector<string> ExchProc::SpinMapStrs() const
  {
  vector<string> SpMaps;			// Strings to be returned
  int c1, c2;                                   // Component indices
  int s1, s2;                                   // Spin indices
  string middle(" <---> ");                     // Separator between LHS & RHS
  string line;					// A single line in return
  int nm = NSpinMaps();				// # of spin maps in process
  unsigned id = 1;                              // Size of spin index used
  if(nm > 9)  id++;				// Adjust this within reason
  if(nm > 99) id++;				// (assumed < 1000 spins!)
  for(int i=0; i<nm; i++)                       // Loop over spin pairs
    {
    c1 = SpinMaps[i].Sub1();			//   Get 1st component index
    c2 = SpinMaps[i].Sub2();			//   Get 2nd component index
    s1 = SpinMaps[i].Spin1();			//   Get 1st spin      index
    s2 = SpinMaps[i].Spin2();			//   Get 2nd spin      index
    line  = ExchProc::Label(c1);		//   Set line to LHS component
    line += string(" ") + Gdec(s1, id);		//   Add LHS spin to line
    line += middle;				//   Add in middle divider
    line += ExchProc::Label(c2);		//   Add RHS component to line
    line += string(" ") + Gdec(s2, id);		//   Add RHS spin to line
    SpMaps.push_back(line);			//   Store line in array
    }
  return SpMaps;				// Return the array
  }






        // Input                pro	: An exchange process (this)
        // 			ostr	: An output stream
        //                      full	: Print amount flag
        //                                 0 = don't print spin mappings(def)
        //                                !0 = print individual spin mappings
        // Output               ostr	: The output stream  is returned
        //                                with the exchange process added

ostream& ExchProc::print(ostream& ostr, int full) const
  {                                                              
  ostr << "\nRate: "  << KRate << "/sec";	// Output the exchange rate
  ostr <<  " LHS = {";				// Output left hand side
  int nl = int(LHSComps.size());		// No. LHS components
  int nr = int(RHSComps.size());		// No. RHS components
  for(int i=0; i<nl; i++)			// Loop LHS components
    {
    ostr << LHSComps[i];
    if(i+1 != nl) ostr << ", ";
    }
  ostr << "}  RHS = {";				// Output right hand size
  for(int j=0; j<nr; j++)			// components in process
    {
    ostr << RHSComps[j];
    if(j+1 != nr) ostr << ", ";
    }
  ostr << "}";
  if(full)					// If full flag is set then
    {						// output all the spin mappings
    for(unsigned k=0; k<SpinMaps.size(); k++)	// Loop over all spin pairs
      ostr << "\n\t" << SpinMaps[k];		// and print spin mapping
    }
  return ostr;
  }

ostream& operator<< (ostream& ostr, const ExchProc& pro)
  { return pro.print(ostr); } 			// Just use the member function


#endif							// ExchProc.cc
