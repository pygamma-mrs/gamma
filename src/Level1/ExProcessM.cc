/* ExProcessM.cc ************************************************-*-c++-*-
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      Mutual Exchange Process                   Implementation	**
**                                                                      **
**      Copyright (c) 2001                                              **
**      Dr. Scott A. Smith                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header: $
**                                                                      **
**************************************************************************

**************************************************************************
**                                                                      **
** This class defines a single, specific, mutual exchange process. A    **
** mutual exchange process contains an exchange rate (1/sec) as well    **
** as an array of spin indices that are in exchange. The latter array   **
** will have at least 2 spins and indicates how the spins exchange.     **
** Spin i is moving into spin i+1 where i = [0,N-2] for N spins. The    **
** lone exception is that the last spin is moving into the first.       **
**                                                                      **
**                 Initial Spin     Final Spin                          **
**                 ============     ==========                          **
**                       i      -->    i+1       (i = [0, N-2])         **
**                      N-1     -->     0                               **
**                                                                      **
** Since this class exists to be used by other classes, it has limited  **
** functionality.                                                       **
**                                                                      **
*************************************************************************/

#ifndef   ExProcM_cc_			// Is the file already included?
#  define ExProcM_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <Level1/ExProcessM.h>		// Include header file
#include <cstdlib>                      // Include functions atoi, atof, etc.
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Basics/StringCut.h>		// Include string cutting 
#include <Basics/Gutils.h>		// Include GAMMA error messages
#include <string>			// Know about libstdc++ strings
#include <fstream>			// Know about file streams
#include <cmath>			// Know about fabs function

using std::string;			// Using libstdc++ strings
using std::list;			// Using libstdc++ STL lists
using std::vector;			// Using libstdc++ STL vectors
using std::ostream;			// Using libstdc++ output streams

//_________________________________________________________________________________
// i                       MUTUAL EXCHANGE PROCESS ERROR HANDLING
//_________________________________________________________________________________


	// Input			PRO	: A process (this)
        // 				eidx	: Error flag
        //				noret   : Flag for return (0=return)
        //                              pname   : String included in error 
        // Output			none 	: Output process error message


void ExchProcM::XPerror(int eidx, int noret) const
  {
  string hdr("Mutual Exchange Process");
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
    case 33: msg=string("No Components Found In Exchange");      break; // (33)
    case 37: msg=string("Cant Set Up From Input Parameter Set"); break; // (37)
    case 38: msg=string("Problems Determining Exchange Comps");  break; // (38)
    case 39: msg=string("Repeated Spin Defined In Process!");    break; // (39)
    case 40: msg=string("Improper Process Definition!");         break; // (40)
    }
  GAMMAerror(hdr, msg, noret);
  }

void ExchProcM::XPerror(int eidx, const string& pname, int noret) const
 
  {
// case 1:  Problems With File pname                                    // (1)
// case 2:  Cannot Read Parameter pname                                 // (2)
// case 3:  Invalid Use Of Function pname                               // (3)
// case 4:  Use Of Deprecated Function pname                            // (4)
// case 5:  Please Use Class Member Function pname                      // (5)
// default: Unknown Error pname                                         // (-1)
  string hdr("Mutual Exchange Process");
  string msg;
  switch (eidx)
    {
    default: GAMMAerror(hdr, eidx, pname, noret); return;  break;
    case 20: msg=string("Cannot Read Parameter ") + pname; break;      	// (20)
    case 22: msg=string("Accessed Spin Component ") + pname
                +string(" Out of Range");                  break;	// (22)
    }
  GAMMAerror(hdr, eidx, noret);
  }
 
 
volatile void ExchProcM::XPfatal(int eidx) const
  {
  XPerror(eidx, 1);			// Output the error message
  if(eidx) XPerror(0,1);		// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

volatile void ExchProcM::XPfatal(int eidx, const string& pname) const
  {
  XPerror(eidx, pname, 1);		// Output the error message
  if(eidx) XPerror(0,1);                // State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

//_________________________________________________________________________________
// ii                MUTUAL EXCHANGE PROCESS PARAMETER SET PARSING
//_________________________________________________________________________________

/* These functions allow for an exchange process to be set from parameters in a
   specified parameter set.                                                      */

//---------------------------------------------------------------------------------
//------------------------- Read in the Process Definition ------------------------
//---------------------------------------------------------------------------------

/* This will be a parameter such as   MutExch(0)  (2) : (0,1,3,5) - Mutual exchange

   The string value defines the components (spins) involved in the exchange
   process. The function getExch will just get the process definition, it doesn't
   bother to parse it. The function parseExch will take the string from getExch
   and split it up into individual spins to be placed into the array Spins.
   Function getComps just oversees the process of getting the components so that 
   it looks like it is done in one step.

   The string value of Exch will be something akin to (0,1,3,5). We need to then
   to parse out the integer values on each side of , but lying between ().
   Also we must insure there are no repeats.                                    */

bool ExchProcM::getExch(const ParameterSet& pset, int idx,
                                                     string& exch, bool warn) const
  {
  exch = string("");                            // We don't know any exchange
  ParameterSet::const_iterator item;         // A pix into parameter list
  string Nidx = "";                             // Name addition per index
  if(idx >= 0) Nidx                             // If index exists then set up
    += string("(")+Gdec(idx)+string(")");       // the parameter name to append
  string pstate, pname=string("MutExch")+Nidx;	// Parameter name for process
  item = pset.seek(pname);			// Pix in pset for MutExch(i)
  if(item != pset.end())			// Retrieve the process (string)
    {
    (*item).parse(pname,exch,pstate);		//   Process as a string
    return true;				//   Return we were successful
    }
  if(warn)					// If we didnt find the exchange
    XPerror(20, pname, 1); 			// process definition, warn
  return false;
  }

bool ExchProcM::parseExch(string& Exval, vector<int>& sps, bool warn) const
  {
//	          First Get Raw Strings For Components in Exchange

  int len = Exval.length();			// Length of input string
  if(!len)					// If there isn't anything
    { 						// we don't know any components
    XPerror(33, 1);
    return false;
    }
  string S = Exval.substr(1, len-2);		// Get string between ( & )
  if(!S.length())				// If there isn't anything on
    { 						// we don't know any components
    XPerror(33, 1);
    return false;
    }

//         At This Point S Should Be Strings Akin To #,#,#,....
//         We Need To Parse Out The Individual #'s (Components)

  string iX;					// String for integer parse
  sps.clear();					// Insure no components yet
  len = S.length();				// Length of raw R.H.S. string
  while(len > 0)				// While the string has length
    {
    iX = cutInt(S);				//   Cut out an integer
    sps.push_back(atoi(iX.c_str()));		//   Store component index
    if(len > 1)	cutBlksXBlks(S, ",");		//   Cut , out for next one
    len = S.length();				//   New length
    }
  return true;
  }

bool ExchProcM::getComps(const ParameterSet& pset, int idx,
                                     std::vector<int>& sps, bool warn) const
  {
  string Exval;
  if(!getExch(pset, idx, Exval, warn))		// Try and get Exch(idx) string
    {						// If we have trouble with this
    if(warn) XPerror(31, 1); 			//   Issue warnings if desired
    return false;				//   Return that we failed
    }
  if(!parseExch(Exval, sps, warn))		// Try for exchange components
    {						// If we have trouble with this
    if(warn) XPerror(38, 1); 			//   Issue warnings if desired
    return false;				//   Return that we failed
    }
  return true;					// Yes, we succeeded
  }

//---------------------------------------------------------------------------------
//------------------------ Read in the Process Exchange Rate ----------------------
//---------------------------------------------------------------------------------

// This will read a parameter such as    Kex(0)  (1) : 600.0 - rate

bool ExchProcM::getRate(const ParameterSet& pset, int idx,
                                                     double& rate, bool warn) const
  {
  rate = 0.0; 					// We don't know any rate
  ParameterSet::const_iterator item;         // A pix into parameter list
  string Nidx = "";                             // Name addition per index
  if(idx >= 0) Nidx                             // If index exists then set up
    += string("(")+Gdec(idx)+string(")");       // the parameter name to append
  string pstate, pname = string("Kex")+Nidx;	// Parameter name for rate
  item = pset.seek(pname);			// Pix pset for Kex(idx)
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
//---------------------- Read/Set The Entire Exchange Process ---------------------
//---------------------------------------------------------------------------------

bool ExchProcM::getXP(const ParameterSet& pset,
                  double& rate, vector<int>& sps, int idx, bool warn) const
  {
  if(!getComps(pset, idx, sps, warn))		// Try for exchange components
    {						// If we have trouble with this
    if(warn) XPerror(31, 1); 			//   Issue warnings if desired
    return false;				//   Return that we failed
    }
  if(!getRate(pset, idx, rate, true))		// Try and get Kex_nm(idx) value
    {						// If we have trouble with this
    if(warn) XPerror(32, 1); 			//   Issue warnings if desired
    return false;				//   Return that we failed
    }
  if(!FCheck(warn))				// Make sure that the spins
    {						// defined are sensible
    if(warn) XPerror(40, 1);
    return false;				//   Return that we failed
    }
  return true;					// Return we read the exchange
  }

bool ExchProcM::setXP(const ParameterSet& pset, int idx, bool warn)
  {
  double          k;				// For exchange rate
  vector<int>     sps;				// For spin components
  if(getXP(pset,k,sps,idx,warn))		// Try and read process
    {						// If we are successful
    KRate    = fabs(k);				//   Set the exchange rate
    Spins    = sps;				//   Set the spin components
    return true;				//   Return successful
    }
  if(warn) XPerror(37, 1); 			// We failed, issue warnings
  return false;					// Return our miserable failings
  }

//__________________________________________________________________________
// iii             MUTUAL EXCHANGE PROCESS CHECKING FUNCTIONS
//__________________________________________________________________________

// The first function just checks the boundaries of the exchange components 
// The 2nd function insures the exchange process integrity

bool ExchProcM::CCheck(int comp, bool warn) const
  {
  if(comp<0 || comp>=(int)Spins.size()) 	// Check that the component
    {						// actually exists
    if(warn) XPerror(10, 1);
    return false;
    }
  return true;
  }

bool ExchProcM::FCheck(bool warn) const
  {
  int ns = Spins.size();			// # spins in process
  int i, j, k;
  for(i=0; i<ns-1; i++)
    {
    k = Spins[i];
    for(j=i+1; j<ns; j++)
      {
      if(k == Spins[j])
        {
        if(warn) XPerror(39, 1);
        return false;
        }
      }
    }
  return true;
  }

//__________________________________________________________________________

//__________________________________________________________________________
// A          MUTUAL EXCHANGE PROCESS CONSTRUCTORS AND DESTUCTORS
//__________________________________________________________________________

//---------------------------------------------------------------------------
//                            Simple Constructors
//---------------------------------------------------------------------------

ExchProcM::ExchProcM() { KRate = 0.0; } 		// No exchange rate

ExchProcM::ExchProcM(const ExchProcM& PR)
  {
  KRate    = PR.KRate;					// Copy exchange rate
  Spins    = PR.Spins;					// Copy spin index array
  } 

//--------------------------------------------------------------------------
//                     Construction From Parameter Set
//--------------------------------------------------------------------------

        // Input                pro     : Exchange process (this) 
        //                      pset    : A parameter set
        //                      ip 	: A process index 
        // Output               pro     : The process is constructed
	//				  from values in the parameter set
	// Note				: Non-indexed parameters are allowed 
	//				  and used if ip < 0 
 
ExchProcM::ExchProcM(const ParameterSet& pset, int ip, int warn)
  {
   if(!setXP(pset, ip, warn?true:false))	// Try & use MutExch(ip) 
     {						// If we have trouble with this
     if(warn)					//   Issue warnings if desired
       {
       if(warn>1) XPfatal(31);			//   Cannot find exchange process
       else       XPerror(31);
       }
     }
   }

//---------------------------------------------------------------------------
//                        Assignment and Destruction
//---------------------------------------------------------------------------

ExchProcM& ExchProcM::operator= (const ExchProcM& PR)
{
  KRate    = PR.KRate;				// New exchange rate
  Spins    = PR.Spins;				// Copy all spin indices

  return (*this);
} 

ExchProcM::~ExchProcM() { }

//___________________________________________________________________________
// B                   Mutual Exchange Process Access Functions
//___________________________________________________________________________

//---------------------------------------------------------------------------
//                             Mutual Exchange Rate
//---------------------------------------------------------------------------

        // Input                pro     : An exchange process (this)
        //                      k       : An exchange rate (1/sec)
        // Output               void    : Exchange rate is set to k
        //                   or double  : Exchange rate is returned

double ExchProcM::Kex() const   { return KRate; }
void   ExchProcM::Kex(double k) { KRate = k;    }
 
//---------------------------------------------------------------------------
//              Mutual Exchange Process Component Index Access
//---------------------------------------------------------------------------
 
        // Input                pro     : A mutual exchange process (this)
        //                      comp    : A component index
        // Output               ic      : Spin index of component
 
int ExchProcM::NComps() const { return Spins.size(); }
int ExchProcM::NSpins() const { return Spins.size(); }

int ExchProcM::Comp(int comp) const
  {
  if(!CCheck(comp))				// If component out of bounds
    XPfatal(22, Gdec(comp));			// quit with error messages
  return Spins[comp];				// Return component index
  }

//____________________________________________________________________________
// C                       MUTUAL EXCHANGE PROCESS SPIN QUERIES
//____________________________________________________________________________
 
        // Input                pro     : A mutual exchange process
	//			i	: Component index 1
	//			j       : Component index 2
        // Output		TF	: True if components are
	//				  directly exchanging
	//                   or 	: True if component i involved
	//				  in the echange at all
	// Note				: For direct exchange the two
	//				  spin must be adjacent Spins

bool ExchProcM::mixes(int i, int j) const
  {
  int ns = Spins.size();			// # of spins
  for(int k=0; k<ns; k++)			// Loop all spins
    {
    if(Spins[k]==i)				//  If found spin i, next
      {						//  must be j for exchange
      if(k<ns-1 && Spins[k+1]==j) return true;
      else if(Spins[0]==j)        return true;
      return false;
      }
    }
  return false;
  }

bool ExchProcM::involves(int i) const
  {
  int ns = Spins.size();			// # of spins
  for(int k=0; k<ns; k++)			// Loop all spins
    if(Spins[k]==i) return true;		// & return true if found
  return false;
  }

// ____________________________________________________________________________
// D                  MUTUAL EXCHANGE PROCESS INPUT FUNCTIONS
// ____________________________________________________________________________

        // Input                pro	: A mutual exchange process (this)
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

bool ExchProcM::read(const string &filename, int idx, int warn)
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

bool ExchProcM::read(const ParameterSet& pset, int idx, int warn)
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
// E                       MUTUAL EXCHANGE PROCESS OUTPUT
//__________________________________________________________________________________

string ExchProcM::ExchStr() const
  {
  string ExchS("Undefined");			// Return string
  int ns = int(Spins.size());			// No. spins involved
  if(!ns) return ExchS; 			// Done if no spins exchange
  string Sep(" --> ");				// Separator
  ExchS = string("");				// Zero return string
  int nd = 1;					// # decimals per spin
  if(ns>9)  nd++;				// to use in output string
  if(ns>99) nd++;				// (don't handle > 999 spins!)
  for(int i=0; i<ns; i++)			// Loop spins in exchange
    ExchS += Gdec(Spins[i],nd) + Sep;		//   Add exchange to string
  ExchS += Gdec(Spins[0]);			// Add last exchanging spin
  return ExchS;
  }

        // Input                pro	: A mutual exchange process (this)
        // 			ostr	: An output stream
        //                      full	: Print amount flag
        //                                 0 = don't print spin mappings(def)
        //                                !0 = print individual spin mappings
        // Output               ostr	: The output stream  is returned
        //                                with the exchange process added

ostream& ExchProcM::print(ostream& ostr, int full) const
  {                                                              
  ostr << "\nRate: "  << KRate << "/sec";	// Output the exchange rate
  ostr << "   Exchange: ";
  int ns = int(Spins.size());			// No. spins involved
  if(!ns) { ostr << "Undefined"; return ostr; }
  for(int i=0; i<ns; i++)			// Loop spins
    ostr << Spins[i] << " --> ";
  ostr << Spins[0];
  return ostr;
  }

ostream& operator<< (ostream& ostr, const ExchProcM& pro)
  { return pro.print(ostr); } 			// Just use the member function

#endif							// ExchProcM.cc
