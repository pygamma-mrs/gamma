/* MultiSys.cc **************************************************-*-c++-*-
**  									**
**                             G A M M A				**
**									**
**      Class Multiple System		            Implementation	**
**									**
**      Copyright (c) 1995						**
**      Nikolai Skrynnikov						**
**      Dr. Scott A. Smith						**
**      1800 E. Paul Dirac Drive					**
**      National High Magnetic Field Laboratory				**
**      Tallahassee FL 32306 USA					**
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
** This class embraces mulitple spin systems that have no scalar	**
** coupling between them. Each spin system is of type sys_dynamic,	**
** containing info on shifts, couplings, relaxtion, and mutual exchange	**
** (dynamics). This class allows one to account for non-mutual exchange	**
** between molecular groups that are not coupled to each other.  Note	**
** that intra-molecular exchange may be accounted for by sys_dynamic	**
** alone.								**
**									**
*************************************************************************/

#ifndef MultiSys_cc_				// Is the file already included?
#define MultiSys_cc_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)	// Using the GNU compiler?
#    pragma implementation			// Then this is the implementation
#  endif

#include <MultiSys/MultiSys.h>		// Include the interface
#include <MultiSys/ExProcess.h>		// Include exchange processes
#include <string>			// Include libstdc++ strings
#include <iostream>			// Include libstdc++ output streams
#include <Basics/ParamSet.h>		// Include GAMMA parameter sets
#include <Basics/StringCut.h>		// Include GAMMA string parsing
#include <Basics/Gutils.h>              // Include GAMMA standard errors
#include <LSLib/sys_dynamic.h>		// Include dynamic spin systems
#include <list>				// Include libstdc++ STL lists

using std::string;			// Using libstdc++ strings
using std::list;			// Using libstdc++ lists
using std::vector;			// Using libstdc++ vectors
using std::ostream;			// Using libstdc++ output streams
using std::ofstream;			// Using libstdc++ output file streams

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                     CLASS MULTI_SYS ERROR HANDLING
// ____________________________________________________________________________

void multi_sys::MSYSerror(int eidx, int noret) const
  {
  string hdr("Multiple System");
  string msg;
  switch (eidx)
    {
//  case 0:  Program Aborting						// (0)
//  case 1:  Problems With Input File Stream     			// (1)
//  case 2:  Problems With Output File Stream    			// (2)
//  case 3:  Can't Construct From Parameter Set  			// (3)
//  case 4:  Cannot Construct From Input File    			// (4)
//  case 5:  Cannot Write To Parameter File      			// (5)
//  case 6:  Cannot Write To Output FileStream   			// (6)
//  case 7:  Cannot Produce Any Output           			// (7)
//  case 8:  Problems With Parameter Set         			// (8)
//  case 9:  Problems During Construction        			// (9)
//  case 10: Cannot Access Internal Component    			// (10)
//  case 11: Cannot Read From Input FileStream   			// (11)
//  case 12: End of File Reached                 			// (12)
//  case 13: Cannot Open ASCII File For Read     			// (13)
//  case 14: Cannot Read From Input File         			// (14)
//  case 15: Cannot Write To Output File         			// (15)
//  default: Unknown Error", noret);                			// (-1)

    default: GAMMAerror(hdr, eidx, noret); return; break;
    case 30: msg=string("Can't Read From Parameter Set");        break;	// (30)
    case 31: msg=string("Can't Determine Number Of Components"); break;	// (31)
    case 32: msg=string("Can't Find Component ASCII Filename");  break; // (32)
    case 33: msg=string("Can't Find Component Population");      break; // (33)
    case 34: msg=string("Can't Get Exchange Process Definition");break; // (34)
    case 35: msg=string("Can't Find Exchange Process Rate");     break; // (35)
    case 36: msg=string("Can't Find Exchange Spin Mappings");    break; // (36)
    case 37: msg=string("Component Field Strengths Differ!");    break; // (37)
    case 38: msg=string("Can't Set Component Spin System");      break; // (38)
    case 39: msg=string("Can't Get Component Hilbert Space");    break; // (39)
    case 45: msg=string("# Spins Not Conserved In Exchange");    break;	// (45)
    case 46: msg=string("Addressed Exchange Process Undefined"); break; // (46)
    case 47: msg=string("Addressed Component Not In System");    break; // (47)
    case 48: msg=string("Mismatch Of Component Field Strength"); break; // (48)
    case 49: msg=string("Cannot Add Component Into System");     break; // (49)
    case 51: msg=string("Can't Determine System Name");          break;	// (51)
    case 52: msg=string("Can't Find All Component Filenames");   break;	// (52)
    case 53: msg=string("Unreasonable # Of Components");         break;	// (53)
    case 54: msg=string("Can't Get Component From Parameters");  break;	// (54)
    case 55: msg=string("Can't Get/Set Component Population");   break;	// (55)
    case 56: msg=string("Can't Get/Set Exchange Process");       break;	// (56)
    case 57: msg=string("Can't Get/Set Exchange Rate");          break;	// (57)
    case 60: msg=string("Can't Make Exchange Mapping Strings");  break;	// (60)
    case 61: msg=string("Can't Make Exchange Process Strings");  break;	// (61)
    case 62: msg=string("Can't Set System From Parameter Set");  break;	// (62)
    case 63: msg=string("Can't Find Any Exchange Processes");    break;	// (63)
    case 64: msg=string("Problems Setting Exchange Processes");  break;	// (64)
    case 65: msg=string("Exchange Spin Count Not Conserved!");   break;	// (65)
    case 66: msg=string("One Or More Bad Exchange Processes!");  break;	// (66)
    case 67: msg=string("Exchanging Spin Has No Spin Partner");  break;	// (67)
    }
  GAMMAerror(hdr, msg, noret);
  }  

void multi_sys::MSYSerror(int eidx, const string& pname, int noret) const
  {
// case 1:  Problems With File pname			                // (1)
// case 2:  Cannot Read Parameter pname					// (2)
// case 3:  Invalid Use Of Function pname				// (3)
// case 4:  Use Of Deprecated Function pname				// (4)
// case 5:  Please Use Class Member Function pname			// (5)
// default: Unknown Error pname						// (-1)
  string hdr("Multiple System");
  string msg;

  switch (eidx)
    {
    default: GAMMAerror(hdr, eidx, pname, noret); return; break;
    case 20: msg=string("Setting ")+pname+string(" Components"); break;	// (20)
    case 30: msg=string("Cannot Set Component ") + pname;        break;	// (30)
    case 31: msg=string("Cant Set Population Of Comp. ")+pname;  break;	// (31)
    case 32: msg=string("Cannot Get Component ") + pname;        break;	// (32)
    case 33: msg=string("Exchange Process ") + pname  
                                      + string(" Is Invalid");   break;	// (33)
    case 60: msg=string("Exchange Process ") + pname  
                                      + string(" Is Invalid");   break;	// (60)
    case 61: msg=pname + string(" Exchange Partner Is Bad");     break;	// (61)
    case 62: msg=string("Isotope Mismatch, ") + pname;           break;	// (62)
    }
  GAMMAerror(hdr, msg, noret);
  }  

volatile void multi_sys::MSYSfatal(int eidx) const
  {                                                                 
  MSYSerror(eidx, 1);				// Output error message
  if(eidx) MSYSerror(0,1);			// Write we're aborting
  GAMMAfatal();					// Clean exit from program
  }

volatile void multi_sys::MSYSfatal(int eidx, const string& pname) const
  {                                                                 
  MSYSerror(eidx, pname, 1);			// Output error message
  if(eidx) MSYSerror(0,1);			// Write we're aborting
  GAMMAfatal();					// Clean exit from program
  }
 
// ____________________________________________________________________________
// ii         MULTI SPIN SYSTEM SETUP FUNCTIONS USING PARAMETER SETS
// ____________________________________________________________________________

/* The purpose of these functions is to set up a multiple spin system from
   parameters in a specified parameter set. Many of these must be kept private
   because their misuse could lead to big problems as they totally tweak with
   the class internal  structure.  For example, the number of components can be
   set without the actual number of allocated components changing. So, their
   misuse would be a very bad thing. If they are private they won't be abused
   (except herein). In general, the "get" functions do NOT alter the spin
   system, they only glean information from the parameter set for the system.
   In contrast the "set" functions will alter the spin system using the
   parameters (gotten from "set" functions).

  Function        Type                            Purpose
 ==========  ==============  ==================================================
 getNComps   int             Get the number of components in MultiSys: NComp
 getMSName   string          Get system name:                          MSysName
 getFName    string          Get component (sub-system) file name:     Fname(#)
 getFNames   vector<string>  Get component (sub-system) file names:    Fname(#)
 getComp     sys_dynamic     Get component of MultiSys:                Various
 getComps    vector<sys_dyn> Get all components of MultiSys:           Various
 getPop      double          Get component population:		       Popul(#)
 getPops     vector<double>  Get populations of all components         Popul(#)
 getNex      int             Get the number of exchange processes:     Exch(#)
*/
 
// ----------------------------------------------------------------------------
//                 Reading Of Basics MultiSystem Parameters
// ----------------------------------------------------------------------------

/* These read very simple parameters associated with the muliple spin system.
   One reads the spin system name and this typically is set to not be of any
   consequence at all.  Another reads the number of components, i.e. the number
   of spin systems within the multiple spin system. This is usually critical
   and warning should be on if it is called.                                 */

bool multi_sys::getNComps(const ParameterSet& pset, int& ncmps, bool warn)
  {
  ncmps = 0;					// We know of no components
  ParameterSet::const_iterator item;		// A pix into parameter list
  string pstate, pname = string("NComp");	// Parameter name, comment
  item = pset.seek(pname);			// Parameter for NComp
  if(item != pset.end())			// Retrieve # of components
    {						// If successful
    (*item).parse(pname,ncmps,pstate);		//   Get number of components
    return true;				//   Return we know the number
    }
  if(warn) MSYSerror(31,1);			// Can't read number components
  return false;					// Return that we don't know it
  }

bool multi_sys::getMSName(const ParameterSet& pset, string& name, bool warn)
  {
  name = string("");				// We don't know any name
  ParameterSet::const_iterator item;         // A pix into parameter list
  string pstate, pname = "MSysName";		// Parameter name, comment
  item = pset.seek(pname);			// Empty single parameter
  if(item != pset.end()) 			// Retrieve the system name
    {						// If successful
    (*item).parse(pname,name,pstate);		//   Get the system name
    return true;   				//   Retrun we know name
    }
  if(warn)					// If not successful &
    {						// warnings are deisred
    MSYSerror(51,1);				//   Can't read system name
    MSYSerror(2, pname, 1);			//   Can't read parameter
    }
  return false;					// Return that we don't know it
  }  
 
// ----------------------------------------------------------------------------
//             Reading Of MultiSystem Components From Parameters
// ----------------------------------------------------------------------------

/* Each MultiSystem contains any number of spin systems of type sys_dynamic.
   When reading in the MultiSystem from a parameter file, all parameters for
   each component can either reside in the file or another filename can be
   specified that contains the component. If the latter is used then we must
   read in additional external files, one for each component.  If the former 
   is used then we do additional parse on the input file. Of course the files
   are immediately read into GAMMA parameter sets, so a file is equivalent
   to a parameter set.                                                        */

bool multi_sys::getFName(const ParameterSet& pset, string& name, 
                                                            int idx, bool warn)
  {
  name = string("");				// We don't know any name
  string Nidx = "";                             // Name addition per index
  if(idx >= 0) Nidx                             // If index exists then set up
    += string("(")+Gdec(idx)+string(")");	// the parameter name to append
  string pstate, pname = "Fname" + Nidx;	// Parameter name, comment
  ParameterSet::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);			// Parameter for Fname(idx)
  if(item != pset.end()) 			// Retrieve component filename
    {						// If successful
    (*item).parse(pname,name,pstate);		//   Get the file name
    return true;   				//   Retrun we know name
    }
  if(warn)					// If not successful &
    {						// warnings are deisred
    MSYSerror(32,1);				//   Can't read comp. filename
    MSYSerror(2, pname, 1);			//   Can't read parameter
    }
  return false;
  }

bool multi_sys::getFNames(const ParameterSet& pset,
                                              vector<string>& names, bool warn)
  {
  bool TF = true;				// Flag we know all names
  int nc = names.size();			// Get number of components
  string fname;					// Single component file name
  for(int i=0; i<nc; i++)			// Loop over all components
    {						// and set the ASCII filename
    TF = TF & getFName(pset, fname, i, warn);	//   Get component filename
    names[i] = fname;				//   Store component filename
    }
  if(!TF && warn)				// If we did not find all the
    MSYSerror(52,1);				// filenames & warnings wanted
  return TF;					// warn we cant find all names
  }

bool multi_sys::getComp(const ParameterSet& pset, int idx,
                                                  sys_dynamic& cmp, bool warn)
  {
  string fname;					// Component file name
  if(getFName(pset,fname,idx,false))		// If filename specified, try
    if(cmp.read(fname)) return true;		// to read compoennt from file
  if(cmp.read(pset,idx,false)) return true;	// Else try to get from pset
  if(warn) MSYSerror(54,1);			// If failed, warn if desired
  return false;
  }

bool multi_sys::getComps(const ParameterSet& pset, int ncmps, 
                                          vector<sys_dynamic>& cmps, bool warn)
  {
  cmps.clear();					// We know no components yet
  sys_dynamic cmp, cmp0;			// A single component
  bool TF =true;				// Assume we'll get them OK
  for(int i=0; i<ncmps; i++)			// Loop over all components
    { cmp = sys_dynamic();						// & try to get each one
    if(getComp(pset,i,cmp,true))		//   Look for component i
      cmps.push_back(cmp);			//   & store if found
    else					//   If not found, take some
      {						//   appropriate action
      TF = false;				//     Flag we have failed
      cmps.push_back(cmp0);			//     Set component empty 
      if(warn) MSYSerror(30, Gdec(i), 1);	//     Can't set component
      }
    }
  return TF;
  }

// ----------------------------------------------------------------------------
//           Reading Of Component Populations And Exchange Processes
// ----------------------------------------------------------------------------

bool multi_sys::getPop(const ParameterSet& pset, int idx,
                                                  double& pop, bool warn) const
  {
  ParameterSet::const_iterator item;         // A pix into parameter list
  string Nidx = "";                             // Name addition per index
  if(idx >= 0) Nidx                             // If index exists then set up
    += string("(")+Gdec(idx)+string(")");        // the parameter name to append
  string pstate, pname = "Popul" + Nidx;	// Name for population param.
  item = pset.seek(pname);			// Parameter for Popul(idx)
  if(item != pset.end()) 			// Try and get the population
    { 						// If successful
    (*item).parse(pname,pop,pstate);		//   Set population value
    return true;				//   Return we were successful
    }
  if (warn) MSYSerror(33, 1);			// Issue warning if desired
  return false;					// Flag we failed to do this
  }
 
bool multi_sys::getPops(const ParameterSet& pset, int ncmps, 
                                         vector<double>& pops, bool warn) const
  {
  pops.clear();					// We know no components yet
  double pop;					// A single population
  bool TF =true;				// Assume we'll get them OK
  for(int i=0; i<ncmps; i++)			// Loop over all components
    {						// & try to get each one
    if(getPop(pset,i,pop,true))			//   Seek comp. i population
      pops.push_back(pop);			//   & store if found
    else					//   If not found, take some
      {						//   appropriate action
      TF = false;				//     Flag we have failed
      pops.push_back(0);			//     Set component empty 
      if(warn) MSYSerror(31, Gdec(i), 1);	//     Can't set component
      }
    }
  return TF;
  }

int multi_sys::getNex(const ParameterSet& pset) const
  {
  string Pbase("Exch(");			// Base name exchange param.
  string Pend(")");				// End of exchange param. name
  int nex = 0;					// No known exchange processes
  ParameterSet::const_iterator item;         // A pix into parameter list
  string pname = Pbase + Gdec(nex) + Pend;	// Exchange process name
  item = pset.seek(pname);			// Get Pix for Exch(nex)
  while(1)					// Loop wile all is OK
    {
    if(item == pset.end()) { return nex; }	//  Stop if no Exch(nex)
    nex++;					//  If Exch(nex) count it
    pname = Pbase + Gdec(nex) + Pend;		//  Next Exch(nex) in series
    item = pset.seek(pname);			//  Pix for Exch(nex)
    }
  return nex;
  }

bool multi_sys::getProcesses(const ParameterSet& pset,
                                       vector<ExchProc>& procs, bool warn) const
  {
  procs.clear();				// We don't know any processes
  ExchProc pro;					// Single exchange processes
  int i=0;					// Exchange process index
  while(pro.read(pset, i, false))		// Try & read exchange process i
    {
    procs.push_back(pro);			// If we can, store it in procs 
    i++; 					// Increment to next process
    }
  if(!i)					// If we didn't read any 
    { if(warn) MSYSerror(63,1); return false; }	// Warn and return our failure
  return true;					// Return we read em OK
  }
 
// ----------------------------------------------------------------------------
//                    Reading Of Enire Multipule Spin System
// ----------------------------------------------------------------------------
 
bool multi_sys::getMsys(const ParameterSet& pset,
              string& name, vector<sys_dynamic>& cmps,
         vector<double>& pops, vector<ExchProc>& procs, bool warn)
  {
  getMSName(pset, name, false);			// Try for system name
  int ncmps;					// Number of components
  if(!getNComps(pset, ncmps, warn))		// Try for number of components
    return false;				// (Quit if failure, important)
  if(!CheckNComps(ncmps, warn))			// Insure # comps. reasonable
    { if(warn) MSYSerror(53,1); return false; }
  if(!getComps(pset,ncmps,cmps,warn))		// Try and read all components
    return false;				// (Quit if failure, important)
  if(!getPops(pset,ncmps,pops,warn))		// Try and read all populations
    return false;				// (Quit if failure, important)
  if(!getProcesses(pset, procs, warn))		// Try for exchange processes
    { if(warn) MSYSerror(64,1); return false; } // (Quit it failue, maybe warn)
  return true;
  }


bool multi_sys::setMsys(const ParameterSet& pset, bool warn)
  {
  string name;					// System name
  vector<sys_dynamic> cmps;			// System components
  vector<double>      pops;			// System component populations
  vector<ExchProc>    procs;			// System exchange processes
  if(!getMsys(pset,name,cmps,pops,procs,warn))	// If we cannot get all of the
    {						// multi spin system read
    if(warn) MSYSerror(62, 1);			//   Warn we cant set it
    return false;				//   Signal that we failed
    }
  _SysName = name;				// Set system name
  _Comps   = cmps;				// Set components (systems)
  _Pops    = pops;				// Set populations
  _Exs     = procs;				// Set exchange processes
  if(!CheckProcs(true))				// Insure all proceses are
    {
    if(warn) MSYSerror(66, 1);			//   Warn we cant set it
    return false;				//   Signal that we failed
    }
  return true;					// Return we succeeded
  }

// ____________________________________________________________________________
// iii                  MULTISYSTEM CHECKING FUNCTIONS
// ____________________________________________________________________________
 
/*  Function                                  Application
   ----------            ------------------------------------------------------
   CheckNComps           See if determined # of components is reasonable
   CheckRange		 See if component index is valid
   CheckProc		 See if process index is valid
   CheckProcs		 See if system exchange processes (read) are valid
   CheckProc		 See if system exchange process (read) is valid

 Every exchange process read must be checked that two rules are not violated: 

    1.) All spins that exchange cannot change their isotope type.
    2.) Exchanging components must have all spins mapped to exchange partner.

 Other checks are done implicitly by the class ExchProc and not done here.   */

bool multi_sys::CheckNComps(int nc, bool warn) const
  {
  if(nc>1 && nc < 100) return true;
  if(warn) MSYSerror(20, Gdec(nc), 1);
  return false;
  }
 
bool multi_sys::CheckRange(unsigned n, bool warn) const
  {
  if(n>=0 && n<_Comps.size()) return true;
  if(warn) MSYSerror(47, 1);
  return false;
  }

bool multi_sys::CheckProc(int ip, bool warn) const
  {
  if((ip>=0) || (ip<int(_Exs.size()))) return true;
  if(warn) MSYSerror(46, 1);
  return false;
  }

bool multi_sys::CheckProcs(bool warn) const
  {
  int np = _Exs.size();				// Number of processes
  for(int i=0; i<np; i++)			// Loop  over all processes
    if(!CheckProc(_Exs[i], warn))			// & check it is consistent
      {
      if(warn) MSYSerror(33, Gdec(i), 1);
      return false;
      }
  return true;
  }

bool multi_sys::CheckProc(const ExchProc& XP, bool warn) const
  {
  int ncl = XP.NCompsLHS();			// LHS components in XP
  int ncr = XP.NCompsRHS();			// RHS components in XP

//           Check Spin Conservation: # Spins LHS == # Spins RHS

  int i, nsl=0, nsr=0;				// Spin in (L/R)HS of exchange
  int L, R;
  for(i=0; i<ncl; i++)				// Find # of spins on LHS
    {						// of the exchange
    L = XP.LHSComp(i);
    nsl += _Comps[L].spins();
    }
  for(i=0; i<ncr; i++)				// Find # of spins on RHS
    {						// of the exchange
    R = XP.RHSComp(i);
    nsr += _Comps[R].spins();
    }
  if(nsl != nsr)				// Insure that they match!
    {
    if(warn) 
      {
      MSYSerror(65, 1);
      string pset = Gdec(nsl) + " Spin";
      if(nsl>1) pset += "s";
      pset += string(" On Left, ") 
           + Gdec(nsr) + " Spin"; 
      if(nsr>1) pset += "s";
      pset += " On Right";
      MSYSerror(60, pset, 1);
      }
    return false;
    }

//       Check Spin Type Conservation: Isotope LHS == Isotope RHS
//                          

  int s, t;
  for(i=0; i<ncl; i++)				// Loop LHS Components
    {
    L   = XP.LHSComp(i);			// LHS Component i
    nsl = _Comps[L].spins();			// Spins in this comp
    for(s=0; s<nsl; s++)			// Loop component spins
      {
      if(!XP.SMap(L,s,R,t))			//   If spin has no mapped
        {					//   partner we don't know
        if(warn)				//   how it exchanges!
          {
          MSYSerror(65, 1);
          string pname = "Spin " + Gdec(s) 
                      + " In " + _Comps[L].name();
          MSYSerror(61, pname, 1);
          MSYSerror(67, 1);
          }
        return false;
        }
      if((_Comps[L]).isotope(s) !=		//   See if same isotope
         (_Comps[R]).isotope(t))
        {					//   partner we don't know
        if(warn)				//   how it exchanges!
          {
          string pname = "XSpin " + Gdec(s) 
                      + " In " + _Comps[L].name();
          MSYSerror(61, pname, 1);
          pname = _Comps[L].symbol(s)
                + string(" Exchanging With ")
                + _Comps[R].symbol(t);
          MSYSerror(62, pname, 1);
          }
        return false;
        }
      }
    }
  return true;
  }

bool multi_sys::CheckField(const spin_system& sys, bool warn) const
  {
  if(!_Comps.size()) return true;
  if(_Comps[0].Omega() == sys.Omega( )) return true;
  if(warn) MSYSerror(48, 1);
  return false;
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A              MULTI-SPIN SYSTEM CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________
 
// ----------------------------------------------------------------------------
//                             Simple Constructors
// ----------------------------------------------------------------------------

multi_sys::multi_sys() { _SysName = ""; }

multi_sys::multi_sys(const multi_sys& msys)
  {
  _SysName = msys._SysName;		// Copy the system name
  _Comps   = msys._Comps;		// Copy the components
  _Pops    = msys._Pops;		// Copy the populations
  _Exs     = msys._Exs;			// Copy the exchange proceses
  }

// ----------------------------------------------------------------------------
//        Construction From Systems, Populations, And Exchange Rates
// ----------------------------------------------------------------------------

        // Input                msys  : Multi system
	//			pop1  : Population of component 1 
	//			sys1  : Component 1
	//			pop2  : Population of component 2 
	//			sys2  : Component 2
	//			k     : Exchange rates between two components
        // Output               none  : Multi_sys is constructed with two
	//				components exchanging
	// Note			      : There must be an equal number of
	//				spin in each compnent system and the 
	//				input spin ordering is preserved upon
	//				exchange.
     
multi_sys::multi_sys(double pop1, sys_dynamic& sys1, double pop2,
                                                   sys_dynamic& sys2, double k)
  {
  if(sys1.Omega() != sys2.Omega()) 		// Insure same Bo strength
    MSYSfatal(37);				// If not, fatal error
  _Pops  = vector<double>(2);			// Construct population array
  _Comps = vector<sys_dynamic>(2);		// Allocate components vector
  _Exs   = vector<ExchProc>(1);			// Allocate exchange proc. vect.
  _Comps[0] = sys1;				// Copy 1st system
  _Comps[1] = sys2;				// Copy 2nd system
  _Pops[0] = pop1;				// Set 1st population
  _Pops[1] = pop2; 				// Set 2nd population
  _SysName  = sys1.name()+"_"+sys2.name();	// Set the system name
  int ns = sys1.spins();			// Spins in first system
  if(ns != sys2.spins()) MSYSfatal(45);		// Insure spin conservation
  _Exs[0].intra_default(0, 1, ns, k);		// Set default exchange process
  }


// sosik - don't see the use of this one at all
/*
multi_sys::multi_sys(int Ncomp)
  {
  if(!CheckNComps(Ncomp)) MSYSfatal(9);		// Insure reasonable # comps.
  _Comps = vector<sys_dynamic>(Ncomp);		// Set components vector
  _Pops  = vector<double>(Ncomp);		// Set populations vector
  _SysName = " ";				// No know systemname
  }
*/


// ----------------------------------------------------------------------------
//                           Destruction, Assignment
// ----------------------------------------------------------------------------

multi_sys& multi_sys::operator= (const multi_sys& msys)
  {
  if(this == &msys) return *this; 	// Do nothing if already equal
  _SysName = msys._SysName;		// Copy the system name
  _Comps   = msys._Comps;		// Copy the system components
  _Pops    = msys._Pops;		// Copy the component populations
  _Exs     = msys._Exs;			// Copy the exchange processes
  return *this;				// Return duplicate operator
  }

/*
void multi_sys::operator= (const multi_sys& msys)
  {
  _SysName = msys._SysName;		// Copy the system name
  _Comps   = msys._Comps;		// Copy the system components
  _Pops    = msys._Pops;		// Copy the component populations
  _Exs     = msys._Exs;			// Copy the exchange processes
  }
*/

multi_sys::~multi_sys() { }

// ____________________________________________________________________________
// B                        MULTI-SPIN SYSTEM NAME
// ____________________________________________________________________________

        // Input                msys    : Multi-spin system
        // Output               string  : Multi-spin system name

const string& multi_sys::name() const          { return _SysName; }
      void    multi_sys::name(const string& N) { _SysName = N;    }
 
// ____________________________________________________________________________
// C                MULTI-SPIN COMPONENT POPULATION ACCESS
// ____________________________________________________________________________

        // Input                msys  : Multi system
        //                      icomp : Component index
        //                      npop  : Population of component icomp
        // Output               none  : Popullation set

void multi_sys::pop(int icomp, double npop)
  {
  if(!CheckRange(icomp))		// Insure component exists
    MSYSfatal(55, Gdec(icomp));         // If not we stop execution
  _Pops[icomp] = fabs(npop); 		// Set the component population
  }

double multi_sys::pop(int icomp) const
  {
  if(!CheckRange(icomp))		// Insure component exists
    MSYSfatal(55, Gdec(icomp));         // If not we stop execution
  return(_Pops[icomp]);                 // Return the component population
  }

double multi_sys::popmin() const
  {
  int nc = _Pops.size();		// Number of components
  if(!nc) return 0;			// Return 0 if no components
  double pmin = _Pops[0];		// Start with 1st component
  for(int i=1; i<nc; i++)		// Loop the rest of components
    if(pmin > _Pops[i]) pmin=_Pops[i];	// and find minimum population
  return pmin;
  }

double multi_sys::popmax() const
  {
  int nc = _Pops.size();		// Number of components
  if(!nc) return 0;			// Return 0 if no components
  double pmax = _Pops[0];		// Start with 1st component
  for(int i=1; i<nc; i++)		// Loop the rest of components
    if(pmax < _Pops[i]) pmax=_Pops[i];	// and find maximum population
  return pmax;
  }
 
// ____________________________________________________________________________
// D                    MULTIPLE SPIN SYSTEM COMPONENT ACCESS
// ____________________________________________________________________________


int multi_sys::NComps() const { return _Comps.size(); }

const sys_dynamic& multi_sys::Comp(int icomp) const
  {
  if(!CheckRange(icomp))                // Insure component exists
    MSYSfatal(32, Gdec(icomp));         // If not we stop execution
  return _Comps[icomp];                 // Return the component
  }

void multi_sys::Comp(int icomp, const sys_dynamic& sys)
  {
  if(!CheckRange(icomp))                // Insure component exists
    MSYSfatal(38, Gdec(icomp));         // If not we stop execution
  double big0 = _Comps[0].Omega();	// Current field strength
  if(sys.Omega() != big0) 		// Insure same Bo in new cmpt
     { 					// If not, fatal error
     MSYSerror(37,1);
     MSYSfatal(38);
     }
  _Comps[icomp] = sys;			// Set the component to sys
  }

void multi_sys::AddComp(const sys_dynamic& sys, double pval)
  {
  if(!CheckField(sys)) { MSYSfatal(49); }
  _Comps.push_back(sys);
  _Pops.push_back(fabs(pval));
  }

void multi_sys::CheckComp(unsigned n) const
  { if(!(n>=0 && n<_Comps.size())) MSYSfatal(47); }

// ____________________________________________________________________________
// E                 MULTI-SPIN EXCHANGE PROCESS ACCESS
// ____________________________________________________________________________

int  multi_sys::NExProcs() const { return _Exs.size(); }
void multi_sys::ExProc(const ExchProc& pr, int iex)
  {
  if(!CheckProc(iex)) MSYSfatal(56);	// Insure exchange process exists
  _Exs[iex] = pr;			// Set the exchange process
  } 

const ExchProc& multi_sys::ExProc(int iex) const
  {
  if(!CheckProc(iex)) MSYSfatal(56); 	// Insure exchange process exists
  return _Exs[iex];			// Return the process
  }

double multi_sys::Kex(int iex) const
  { 
  if(!CheckProc(iex)) MSYSfatal(57); 	// Insure exchange process exists
  return _Exs[iex].Kex(); 
  }

void multi_sys::Kex(double K, int iex)
  {
  if(!CheckProc(iex)) MSYSfatal(57); 	// Insure exchange process exists
  _Exs[iex].Kex(K);
  }

int multi_sys::NCompsLHS(int iex) const
  {
  if(!CheckProc(iex)) MSYSfatal(56); 	// Insure exchange process exists
  return _Exs[iex].NCompsLHS();
  }
int multi_sys::NCompsRHS(int iex) const
  {
  if(!CheckProc(iex)) MSYSfatal(56); 	// Insure exchange process exists
  return _Exs[iex].NCompsRHS();
  }

// ____________________________________________________________________________
// F            GLOBAL & INDIVIDUAL COMPONENT SPIN SYSTEM FUNCTIONS
// ____________________________________________________________________________

/* This section contains typical GAMMA spin system functions but they are
   of two types: global and specific.  The Global functions are those which
   apply to all components (spin systems). That is, use of the function on an
   individual component will produce the same result as use on the multi_sys.
   The specific functions act on a single component within the system, but not
   on all components.

   Individual components are spin systems. Since class multi_sys is NOT derived
   from any base spin system class, access functions to individual components
   are NOT automatically inherited. Thus, they must either be explicitly
   provided for here or users must obtain the system (component), alter it,
   then put it back into the multiple spin system.

   Note that for each of the functions, Fct(cmp,a,b,...), users could also just
   directly use Comp(cmp).Fct(a,b,....) to attain the same results. There are
   of course limited functions that are set here. Only those that do not ruin
   the multisys structure are allowed, & only those often used are set up.   */
 
// ----------------------------------------------------------------------------
//                      Homonuclear Versus. Heteronuclear
// ----------------------------------------------------------------------------

bool multi_sys::homonuclear(int cmp) const
  {
  if(cmp >= 0) return _Comps[cmp].homonuclear();
  for(unsigned i=0; i<_Comps.size(); i++)
    if(!_Comps[i].homonuclear()) return false;
  return true;
  }
 
bool multi_sys::heteronuclear(int cmp) const
  {            
  if(cmp >= 0) return _Comps[cmp].heteronuclear();
  for(unsigned i=0; i<_Comps.size(); i++)
    if(!_Comps[i].heteronuclear()) return false;
  return true;
  }
 
// ----------------------------------------------------------------------------
//                      Spin Hilbert And Liouville Spaces
// ----------------------------------------------------------------------------
 
int multi_sys::HS(int icomp) const
  {
  if(icomp>=0)
    {
    if(!CheckRange(icomp))              // Insure component exists
      MSYSfatal(39, Gdec(icomp));       // If not we stop execution
    return _Comps[icomp].HS();
    }
  int HS=0;
  for(unsigned i=0; i<_Comps.size(); i++)
    HS += _Comps[i].HS();
  return HS;
  }

int multi_sys::LS(int icomp) const
  {
  if(icomp>=0)
    {
    if(!CheckRange(icomp))              // Insure component exists
      MSYSfatal(39, Gdec(icomp));       // If not we stop execution
    return _Comps[icomp].HS()*_Comps[icomp].HS();
    }
  int LS = 0;
  for(unsigned i=0; i<_Comps.size(); i++)
    LS += (_Comps[i].HS())*(_Comps[i].HS());
  return LS;
  }

vector<int> multi_sys::HSs() const
  {
  int nc = _Comps.size();		// Number of components
  vector<int> hss(nc);			// Vector entry each component
  for(int i=0; i<nc; i++)
    hss[i] = _Comps[i].HS();
  return hss;
  }

vector<int> multi_sys::LSs() const
  {
  int nc = _Comps.size();		// Number of components
  vector<int> lss(nc);			// Vector entry each component
  for(int i=0; i<nc; i++)
    lss[i] = _Comps[i].HS()*_Comps[i].HS();
  return lss;
  }

// ----------------------------------------------------------------------------
//                      SPECTROMETER FREQUENCY MANIPULATIONS
// ----------------------------------------------------------------------------

///Center Multi Spin System Spectrometer Frequency Functions

//   These are set to mimic the functionality of lower spin system classes
        // Input                freq : Spectrometer frequency (MHz)
        // Output               none : System 1H Larmor set in MHz

        // Input                freq : Larmor frequency (MHz)
        //                      iso  : spin isotope type
        // Output               none : 1H Larmor stored in MHz
        // Output               double : spectrometer frequency in MHz
        // Output               iso    : isotope label


//------------------------------ Set Omega ------------------------------------

void multi_sys::Omega(double freq)
  {
  for(unsigned i=0; i<_Comps.size(); i++)	// Loop over the components
    _Comps[i].Omega(freq);			// and just use overload
  }

void multi_sys::Omega(double freq, const string& iso)
  {
  for(unsigned i=0; i<_Comps.size(); i++) 	// Loop over the components
    _Comps[i].Omega(freq, iso);			// and just use overload
  }

//------------------------------ Get Omega ------------------------------------

double multi_sys::Omega( ) const { return _Comps[0].Omega(); }
double multi_sys::Omega(const string& iso) const { return _Comps[0].Omega(iso); }

const string multi_sys::symbol(int comp, int spin) const
  { return _Comps[comp].symbol(spin); }
 
// ____________________________________________________________________________
// G                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//        Functions To Make A Parameter Set From A Multi Spin System
// ----------------------------------------------------------------------------

/* This class contains a name, a vector of components (spin systems), a vector
   of component populations, and a vector of component exchange processes. We
   will write our own name as well as our output our own component populations,
   however we will allow the components (sys_dynamic) and exchange processes
   (ExchProc) to write themselves since they are in a class structure.  For 
   them we need only loop over all components and have themselves perform this
   task with the proper indexing. Note that the multi_sys is placed into only
   one parameter set (unlike input which may occur from multiple ASCII files).
   Also note that parameter prefixing thorough use of [#] is NOT allowed here
   Only one multi_sys may be put into the same parameter set.

     Function                                 Purpose
   ------------         -------------------------------------------------------
   ParameterSet         Convert system into a parameter set
   operator +=          Adds system to existing parameter set                */

multi_sys::operator ParameterSet( ) const
  { ParameterSet pset; pset += *this; return pset; }

void operator+= (ParameterSet& pset, const multi_sys& msys)
   {
   string pname;                        // Now add unique sys_dynamic parameters
   string pdata;
   string pstate;
   SinglePar par;

   pname = string("MSysName");		// Add multi spin system name
   pstate = string("Name of the MultiSys Spin System");
   pdata = msys.name();
   par = SinglePar(pname, 2, pdata, pstate);
   pset.push_back(par);

   pname = string("NComp");		// Add Number of components
   pstate = string("Number of Component Spin Systems");
   pdata = Gdec(msys.NComps());
   par = SinglePar(pname, 0, pdata, pstate);
   pset.push_back(par);

   string pstate0 = string("Population of Component ");
   string pname0 = string("Popul(");
   int i;
   for(i=0; i< msys.NComps(); i++)
     {
     pname = pname0 + Gdec(i) + string(")");
     pstate = pstate0 + Gdec(i);
     par = SinglePar(pname,msys.pop(i),pstate);
     pset.push_back(par);
     }
// sosi - This routine is far from complete!

//   string pstate0 = string("Filename of Dynamic System ");
//   string pname0 = string("fname(");
//   for(int i=0; i< msys.NComps(); i++)
//     {
//     pname = pname0 + Gdec(i) + string(")");
//     pstate = pstat0 + Gdec(i);
//     pset.push_back(SinglePar(pname,msys.pop(i),pstate));
//     }

   string pproc = "EXCH";
   string pKex = "Kex";
   string pprocS = "Exchange Process ";
   string pKexS = "Exchange Rate (1/sec), Process ";
   for(i=0; i<msys.NExProcs(); i++)
     {
     pname = pproc + string("(") + Gdec(i) + string(")");
     pstate = pprocS + Gdec(i);
     par = SinglePar(pname,msys.Kex(i),pstate);
     pset.push_back(par);
     pname = pKex + string("(") + Gdec(i) + string(")");
     pstate = pKexS + Gdec(i);
     par = SinglePar(pname,msys.Kex(i),pstate);
//     item = pset.seek(par); 			// Pix in pset for NComp
     }
   return;
   } 
 
// ----------------------------------------------------------------------------
//    Functions To Output A Multi Spin System To ASCII From A Parameter Set
// ----------------------------------------------------------------------------

/* These functions write the multiple spin system into an ASCII file in GAMMA
   Parameter Set format.  That is, the resulting ASCII file may be read into
   any GAMMA program to create a multiple spin system identical to that written.

           Input                msys    : Multi_sys spin system (this)
                                filename: Output file name
                            or  ofstr   : Output file stream
                                warn    : Warning level
           Output               none    : System is written as a
                                          parameter set to file filename
                                          or into output file stream ofstr

        //                      basename: Component output filename base
        // Output               none    : Multi_sys spin system is written as a
        //                                parameter set to file filename & each
        //                                component written to a file with name
        //                                basename+i.dsys - basename is given
        //                                as an arguemnt and i is an integer */

void multi_sys::write(string& filename, string basename)
  {
  ofstream ofstr(filename.c_str());	// Open filename for input
  if(!ofstr.good())			// If file bad then exit
    {
    MSYSerror(2,filename,1);
    MSYSfatal(1); 
    }
  ParameterSet pset;			// Declare a parameter set
  pset += (*this);			// Add in spin system parameters
  ParameterSet::const_iterator item;	// A pix into parameter list
  item = pset.begin();			// Get Pix of first parameter
  while(item != pset.end())
     {
     (*item).write(ofstr);
     item++;
     }
  string compname;
  for(int j=0; j<int(_Comps.size()); j++)// Loop over each component
    { 
    compname = basename+Gdec(j)+".dsys";	//	Output filename
    (_Comps[j]).write(compname);		// 	Write the component
    }
  return;
  }

// ____________________________________________________________________________
// H                        SPIN SYSTEM INPUT FUNCTIONS
// ____________________________________________________________________________

        // Input                msys    : Multi_sys spin system (this)
        //                      filename: Input filename
        //                      pset    : Input parameter set
        //                      warn    : Warning output level
        //                                      0 = no warnings
        //                                      1 = warnings
        //                                     >1 = fatal warnings
        // Output               TF      : Spin system is filled with
        //                                parameters read from file
        //                                TRUE if read is successful
        // Note                         : The file should be an ASCII file
        //                                containing recognized sys parameters

bool multi_sys::read(const string& filename, int warn)
   {
   ParameterSet pset;			// A GAMMA parameter set
   if(!pset.read(filename, warn?1:0))	// Try and read in pset
    {                                   // If we cannot read the file then
    if(warn)                            // we'll issue warnings as desired
      {
      MSYSerror(1, filename);		//      Problems with file
      if(warn>1) MSYSfatal(14);		//      Can't read from input file
      else       MSYSerror(14);		//      Fatal if warn big enough
      }
    return false;			// Well flag we didn't read & exit
    }
  return read(pset, warn);		// Fill up spin_sys with parameters
  }

bool multi_sys::read(const ParameterSet& pset, int warn)
  {
  bool TF = setMsys(pset, warn?true:false);	// Use overload to read
  if(!TF)                                       // If setSsys didn't handle
    {                                           // the system read from pset
    if(warn)                                    // then we'll issue warnings
      {						// & maybe even stop program
      MSYSerror(8, 1);				//   Problems with pset
      if(warn>1) MSYSfatal(30);			//   Can't read from pset
      else       MSYSerror(30, 1);		//   Fatal if warn big enough
      }
    }
  return TF;
  }

 
        // Input                msys    : Mulit_sys spin system (this)
        //                      argc    : Number of arguments
        //                      argv    : Vector of argc arguments
        //                      argn    : Argument index
        // Output               void    : The parameter argn of array argc
        //                                is used to supply a filename
        //                                from which the spin system is read
        //                                If the argument argn is not in argv,
        //                                the user is asked to supply a filename
        // Note                         : The file should be an ASCII file
        //                                containing recognized sys parameters
        // Note                         : The spin system is modifed (filled)
 
std::string multi_sys::ask_read(int argc, char* argv[], int argn)
  {
  string filename;					// Name of spin system file
  query_parameter(argc, argv, argn,			// Get filename from command
       "\n\tSpin system filename? ", filename);		// line or ask for them
  read(filename);					// Read system from filename
  return filename;
  }
 
std::string multi_sys::ask_read(int argc, char* argv[], int argn, const std::string& def)
  {
  std::string msg = "\n\tSpin system filename ["	// Query we will ask if
                    + def + "]? ";			// it is needed
  std::string filename = def;				// Name of spin system file
  ask_set(argc,argv,argn,msg,filename);			// Or ask for it
  read(filename);					// Read system from filename
  return filename;					// Return filename
  }
 
 
// ____________________________________________________________________________
// I                  MULTIPLE SPIN SYSTEM OUTPUT FUNCTIONS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//              Functions To Modularize Exchange Process Output
//-----------------------------------------------------------------------------

/* The next function builds an array of spin mappings for a specified exchange
   process. The vector will contain a list of lines similar to
               
                            A 0 <---> B 2  13C
   
   Here the lone letter indicate components (spin systems) & the lone integers
   spins.  At the end the type of spin in the exchange is listed. Two added
   details. First, the A 0 <---> B 2 part of each line is supplied in the base
   class ExchProc so we need only ask that class to provide them. Second we
   will format the output strings so that they will precisely fit under the
   header used in the print function of this class.  Thus, the printing will
   appear akin to
                           Spin Mappings   Type
                          ---------------  ----
                           A 0 <---> B 2   13C
                           A 1 <---> B 3    1H

   Each string will be of length 22 (we use 1 column before the header too.) */

vector<string> multi_sys::SpinMapStrs(int exp) const
  {
  if(!CheckProc(exp)) MSYSfatal(60);		// Insure process exists
  vector<string> SPMaps;			// Strings to be returned
  vector<string> XMaps=_Exs[exp].SpinMapStrs();	// Strings from ExchProc
  int len    = XMaps[0].length();		// Length of spin map strings
  string sp1=string((17-len)/2, ' ');		// Spacer before & after string
  int nm     = XMaps.size();			// # of spin mappings
  string sp2, itype, line;			// For isotope type string
  int c1, s1;					// Component & spin indices
  for(int i=0; i<nm; i++)			// Loop over spin maps
    {
    c1 = (_Exs[exp].SMap(i)).Sub1();		//   Get 1st component index
    s1 = (_Exs[exp].SMap(i)).Spin1();		//   Get 1st spin      index
    itype = _Comps[c1].symbol(s1);		//   Get spin isotope type 
    len = itype.length();			//   Get isotope string length
    switch(len)					//   Buffer isotope string
      {						//   for print column alignment
      default: case 5: break;			//     No space if 5 chars
      case 4:  case 3: sp2=string(" "); break;	//     Single space here
      case 2:          sp2=string("  "); break;	//     Two space needed
      }
    line = sp1 + XMaps[i] + sp1 + sp2 + itype;	//   Set line to i spin map
    SPMaps.push_back(line);
    }
  return SPMaps;
  }

vector<string> multi_sys::LHSStrs() const
  {
  vector<string> lhss;				// Array of return strings
  int nex = _Exs.size();			// # of exchange processes
  for(int i=0; i<nex; i++)
    lhss.push_back(_Exs[i].LHSStr());
  return lhss;
  }

vector<string> multi_sys::RHSStrs() const
  {
  vector<string> rhss;				// Array of return strings
  int nex = _Exs.size();			// # of exchange processes
  for(int i=0; i<nex; i++)
    rhss.push_back(_Exs[i].RHSStr());
  return rhss;
  }

vector<string> multi_sys::EXPStrs() const
  {
  vector<string> lhss = LHSStrs();		// Array of L.H.S. strings
  vector<string> rhss = RHSStrs();		// Array of R.H.S. strings
  vector<string> exps;				// Array of return strings
  string middle(" <===> ");			// L.H.S. to R.H.S. separator
  string line;
  int nex = _Exs.size();			// # of exchange processes
  int i, ll, rl, llen=0, rlen=0;		// String lengths
  for(i=0; i<nex; i++)				// Find length of longest
    { 						// L.H.S. & R.H.S. strings
    llen = gmax(int(lhss[i].length()), llen);
    rlen = gmax(int(rhss[i].length()), rlen);
    } 						// L.H.S. & R.H.S. strings
  for(i=0; i<nex; i++)				// SFind length of longest
    {
    line = lhss[i];				//  Set L.H.S. of line
    ll = lhss[i].length();			//  Get L.H.S. length
    if(ll < llen) line += string(llen-ll, ' ');	//  Pad L.H.S. to same len
    line += middle;				//  Add <===> divider 
    line += rhss[i];				//  Add R.H.S. of line
    rl = lhss[i].length();			//  Get L.H.S. length
    if(rl < rlen)
    if(rl < rlen) line += string(rlen-rl, ' ');	//  Pad R.H.S. to same len
    exps.push_back(line);			//  Set this line
    }
  return exps;
  }


        // Input                msys  : Multi system
        //                      ostr  : Output stream
        //                      full  : Print amount flag
	//				  !0 = print individual systems (def)
	//				   0 = don't print individual systems
        // Output               ostr  : The output stream  is returned
        //                              with the spin system added 
 
ostream& multi_sys::print(ostream& ostr, int full) const
  {
//               First Insure Some Components Exist

  if(!NComps())
    {
    string hdr("Empty Multiple Spin System");
    int hl = hdr.length();
    string spacer = string(40-hl/2, ' ');
    ostr << "\n\n" << spacer << hdr << "\n";
    return ostr;
    }
  else
    {
    string hdr = string("Multiple Spin System ") + _SysName;
    int hl = hdr.length();
    string spacer = string(40-hl/2, ' ');
    ostr << "\n\n" << spacer << hdr << "\n";
    }

//	Output The Component Spin Systems & Their Populations

  int nl;
  string cname;
  string spcr(21, ' ');
  ostr << "\n" << spcr << "Component       Name        Population";
  ostr << "\n" << spcr << "---------  ---------------  ----------";
  for(int i=0; i<int(_Comps.size()); i++)
    {
    cname = _Comps[i].name();
    nl    = cname.length();
    if(nl<15) cname += string(15-nl, ' ');
    else      cname = cname.substr(0, 15); 
    ostr << "\n" << spcr << "    " << ExchProc::Label(i) << "    ";
    ostr << "  " << cname;
    ostr << "  " << Gform("%8.3f", _Pops[i]);
    }

//		Output The Defined Exchange Processes

  vector<string> EXPs = EXPStrs();			// Process strings
  string ep1, ep2, eph1, eph2, eps1, eps2;		// Col. alignment
  int pl = EXPs[0].length();				// Length of strings
  if(pl < 16)
    {
    eps1 = string((16-pl)/2, ' ');
    eps2 = string(16-pl-eps1.length(), ' ');
    }
  else if(pl > 16)
    {
    ep1  = string((pl-16)/2, ' '); 
    eph1 = string((pl-16)/2, '-');
    ep2  = string(pl-16-ep1.length(), ' ');
    eph2 = string(pl-16-ep1.length(), '-');
    }

  vector<string> SMs;					// Spin map strings
  string EXhdr = ep1 +string("Exchange Process")+ep2;	// Exchange header 
  string EXhul = eph1+string("----------------")+eph2;
  string Khdr  = " Rate(1/sec)    Spin Mappings   Type";			// Exchange header 
  string Khul  = "-------------  ---------------  ----";
  string hdr   = EXhdr + string("  ") + Khdr;
  string hdrul = EXhul + string("  ") + Khul;
  int hl = hdr.length();
  string spacer = string(40-hl/2, ' ');
  ostr << "\n";
  ostr << "\n" << spacer << hdr;
  ostr << "\n" << spacer << hdrul;

  hl = spacer.length()  + eps1.length()
     + EXPs[0].length() + eps2.length()
     + 1 + 1 + 13 + 1;
  string bigspc = string(hl, ' ');
  for(unsigned k=0; k<EXPs.size(); k++)
    {
    SMs = SpinMapStrs(k);				// Get spin mapping strings
    ostr << "\n" << spacer << eps1 << EXPs[k]		// Print process
         << eps2 << " ";
    ostr << " " << Gform("%13.4f", _Exs[k].Kex()); 	// Print the process rate
    ostr << " " << SMs[0];				// Print spin map 0, type
    for(unsigned l=1; l<SMs.size(); l++)		// Loop other spin maps &
      ostr << "\n" << bigspc << SMs[l];			// print each spin map, type
    }

// sosi - this next bit works but causes programs to later core dump.
//        probably problems with str_vector still.  9/24/97
//        there are destructors that go along with this too

//	     Store Printed Lines Of Individual Spin Systems (Components)

//  str_vector StrVec, SysStrs[ncomp];			// string vectors for systems
//  ifstream ifstr;					// Temp input file
//  ofstream ofstr;					// Temp output file
//  string cmpfile("cmp.tmp"), input;			// Temp file name
//  int lines = 0;					// File line count
//  for(i=0; i<ncomp; i++)				// Loop over components
//    {
//    ofstr.open(cmpfile);				//	Open temp file
//    (_Comps[i]).print(ofstr);				//	Write component to file
//    ofstr.close();					//	Close the file now
//    ifstr.open(cmpfile);				//	Open temp file					
//    lines = 0;						//	Zero line count
//    while( readline(ifstr, input)!=0 ) lines++; 	// 	Loop file, count lines
//    ifstr.close();					// 	Close the file
//    StrVec = str_vector(lines);				//	Make vector for strings
//    ifstr.open(cmpfile);				//	Open temp file					
//    lines = 0;						//	Zero line count
//    while( readline(ifstr, input)!=0 )	 		// 	Loop file, read lines
//      {
//      StrVec(lines) = input; 				//	Put strings into vector
//      lines++;
//      }
//    ifstr.close();					// 	Close the file
//    SysStrs[i] = StrVec;				//	Store component strings
//    }

//	       Obtain Individual Spin Systems Print Parameters

//  int maxlines[ncomp], maxwidths[ncomp];		// string vector lengths, widths
//  string Sempty[ncomp], Sdash[ncomp], Chdrs[ncomp], CHtmp;
//  int leni, lenf;
//  for(i=0; i<ncomp; i++)				// Loop over components
//    {
//    SysStrs[i].widen();					//	Set strings even length
//    maxlines[i] = (SysStrs[i]).size();			//	Store vector lengths
//    maxwidths[i] = (SysStrs[i]).maxlen();		//	Store vector widths
//    Sempty[i] = string(maxwidths[i], ' ');		//	Store blank width
//    Sdash[i] = straing(maxwidths[i], '_');		//	Store ____ width
//    CHtmp = string("Component ") + Gdec(i);		//	Base of component header
//    leni = (maxwidths[i]-CHtmp.length())/2; 
//    lenf = maxwidths[i] - CHtmp.length() - leni;
//    Chdrs[i] = string(leni, ' ') + CHtmp + string(lenf, ' ');
//    }

//				Print Individual Spin Systems

//  string sep("   ");
//  int cmpi=0, cmpf=cmpi;				// Bounds of components printed
//  int maxLine=100, maxL=0, Llen;			// Max. allowed line length
//  for(Llen=0; cmpf<ncomp && Llen<maxLine; cmpf++)	// Set cmpf to not exceed length 
//    Llen += maxwidths[j];
//  i=0;							// Index of system printing
//  while(i<ncomp)					// Loop over all systems
//    {
//    for(i=cmpi, maxL=0; i<cmpf; i++)			// Find number of lines being
//      if(maxlines[i] > maxL) maxL=maxlines[i];		// printed these components
//    ostr << "\n\n\n";					// Print a spacer
//    for(i=cmpi; i<cmpf; i++) ostr << Chdrs[i] << sep; 	// Print component headers
//    for(i=cmpi; i<cmpf; i++) ostr << Sdash[i] << sep; 	// Print component headers
//    ostr << "\n";
//    for(lines=0; lines<maxL; lines++)
//      {
//      for(i=cmpi; i<cmpf; i++)				// Loop over components
//        {
//        if(lines < maxlines[i])
//          ostr << SysStrs[i].get(lines);
//        else ostr << Sempty[i];
//        if(i < cmpf-1) ostr << sep;
//        }
//      ostr << "\n";
//      }
//    cmpi = cmpf;
//    for(Llen=0; cmpf<ncomp && Llen<maxLine; cmpf++)
//      Llen += maxwidths[j];
//    }

  if(full)					// If full output requested
    {
    ostr<<"\n";
    for(int i=0; i<int(_Comps.size()); i++)	// then print out each spin
      { 					// system also
      hdr = "Multi-System Component ";
      hdr += ExchProc::Label(i);
      hl  = hdr.length();
      spacer = string(40-hl/2, ' ');
      ostr << "\n" << spacer << hdr;
      ostr << "\n" << _Comps[i];
      }
    }

//			     Now Clean Up Before Leaving

  ostr.flush();					// Flush everything
  return ostr;					// Return the output stream
  }
 
ostream& operator<< (ostream& ostr, const multi_sys& sys)
  { return sys.print(ostr); } 			// Use member function
 
#endif								// MultiSys.cc
