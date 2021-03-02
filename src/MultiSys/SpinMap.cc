/* SpinMap.cc ***************************************************-*-c++-*-
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      Spin Mapping                            Implementation 		**
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
** This file contains a simple auxiliary class to support mulitple      **
** spin systems and non-mutual exchange. Specifically it handles the    **
** mapping of a spin in one system onto the spin of another system. The **
** class has no knowledge of of spins or spin systems, so it acutally   **
** maps an index of one component into an index of a second component.  **
**                                                                      **
** It does that by tracing 4 integers.  The first two are a component   **
** and (spin) index, as are the latter two.  The component indices      **
** may be used to index particular spin systems in a variable of type   **
** multi_sys. The (spin) indices then would index particular spins in   **
** their respective systems (or multi_sys components).                  **
**                                                                      **
*************************************************************************/

#ifndef _SpinMap_cc_			// Is the file already included?
#define _SpinMap_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// Then this is the implementation
#endif

#include <MultiSys/SpinMap.h>		// Include header file
#include <string>			// Know about libstdc++ strings
#include <cstdlib>                      // Include functions atoi, atof, etc.
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Basics/StringCut.h>		// Include string cutting 
#include <Basics/Gutils.h>		// Include GAMMA error messages
#include <fstream>

//_________________________________________________________________________________
// i                        CLASS SPIN MAP ERROR HANDLING
//_________________________________________________________________________________

        // Input			SMap	: A spin map (this)
        //                              eidx    : Error flag
        //                              noret   : Flag for return (0=return)
        //                              pname   : String included in error
        // Output                       none    : Output process error message


void SpinMap::SMerror(int eidx, int noret) const
  {
  std::string hdr("Spin Map");
  std::string msg;
  switch (eidx)
    {
//  case 0:  Program Aborting                                           // (0)
//  case 1:  Problems With Input File Stream                            // (1)
//  case 2:  Problems With Output File Stream                           // (2)
//  case 3:  Can't Construct From Parameter Set                         // (3)
//  case 4:  Cannot Construct From Input File                           // (4)
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
    case 30: msg=std::string("Can't Read From Parameter Set");        break; // (30)
    case 31: msg=std::string("Can't Set Map From Parameter Set");     break; // (31)
    case 32: msg=std::string("Unable To Set From Parameters");        break; // (32)
    case 33: msg=std::string("Cannot Map Spin To Itself!");           break; // (33)
    }
  GAMMAerror(hdr, msg, noret);
  }


void SpinMap::SMerror(int eidx, const std::string& pname, int noret) const
  {
  std::string hdr("Spin Map");
  std::string msg;
  switch (eidx)
    {
    default: GAMMAerror(hdr, eidx, pname, noret); return;  break;
//  case 1:  Problems With File pname                                   // (1)
//  case 2:  Cannot Read Parameter pname                                // (2)
//  case 3:  Invalid Use Of Function pname                              // (3)
//  case 4:  Use Of Deprecated Function pname                           // (4)
//  case 5:  Please Use Class Member Function pname                     // (5)
//  default: Unknown Error pname					// (-1)
    case 20: msg=std::string("Cannot Read Parameter ") + pname;       break; // (20)
    case 21: msg=std::string("Unreasonable Component Index") + pname; break; // (21)
    case 22: msg=std::string("Unreasonable Spin Index") + pname;      break; // (22)
    }
  GAMMAerror(hdr, eidx, noret);
  }

volatile void SpinMap::SMfatal(int eidx) const
  {
  SMerror(eidx, 1);				// Output the error message
  if(eidx) SMerror(0,1);			// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

//_________________________________________________________________________________
// ii                   CLASS SPIN MAP PARAMETER SET PARSING
//_________________________________________________________________________________

/* These functions allow for a spin map to be set from parameters in a
   specified parameter set.                                                      */

        // Input                SMap	: A spin map (this)
        //                      pset    : Input parameter set
        //                      idx     : Exchange process index
        //                      mdx     : Spin mapping index
        //                      warn    : Warning flag
        // Output               TF      : Spin map is filled with
        //                                parameters read from file

bool SpinMap::getSMStr(const ParameterSet& pset, int idx, int mdx,
                                                     std::string& sm, bool warn) const
  {
  sm = std::string("");				// We don't know any spin map
  std::string pstate, pname = "Smap(";               // Parameter name base, statement
  if(idx != -1)					// Add suffix if desired	
    {						// We allow either 1 or 2
    pname += Gdec(idx) ;			// indices here
    if(mdx !=-1) pname += std::string(",")+Gdec(mdx);// Smap, Smap(i), Smap(i,j)
    pname += std::string(")");
    }
  ParameterSet::const_iterator item;		// A pix into parameter list
  item = pset.seek(pname);                      // Pix in pset for Smap(i,j)
  if(item != pset.end())                        // Retrieve the process (string)
    {
    (*item).parse(pname,sm,pstate);		//   Process as a string
    return true;                                //   Return we were successful
    }
  if(warn)                                      // If we didnt find the spin map
    SMerror(20, pname, 1);                      // definition, warn
  return false;
  }

bool SpinMap::getSM(const ParameterSet& pset, int idx, int mdx,
                  int& comp1, int& spin1, int& comp2, int& spin2, bool warn) const
  {
  std::string sm;					// String for spin map parameter
  if(!getSMStr(pset,idx,mdx,sm,warn))		// Try to get spin map string
    { if(warn) SMerror(31, 1); return false; }	// Quit and mayb3 warn if failed
  cutBlksXBlks(sm, "(");			// Cut initial (
  comp1 = atoi(cutInt(sm).c_str());		// Cut integer (1st component)
  cutBlksXBlks(sm, ")");			// Cut out )
  spin1 = atoi(cutInt(sm).c_str());		// Cut out an integer (spin)
  cutBlksXBlks(sm, "(");			// Cut initial (
  comp2 = atoi(cutInt(sm).c_str()); 		// Cut integer (2nd component)
  cutBlksXBlks(sm, ")");			// Cut out )
  spin2 = atoi(cutInt(sm).c_str());		// Cut out an integer (spin)
  return Check(comp1,spin1,comp2,spin2,warn);	// Check these seem OK
  }

bool SpinMap::setSM(const ParameterSet& pset, int idx, int mdx, bool warn)
  {
  int s1, s2;					// Spin indices
  int c1, c2;					// Component indices
  if(!getSM(pset,idx,mdx,c1,s1,c2,s2,warn))	// Try to read spin map
    {
    if(warn) SMerror(32,1);			//   Issue warning
    return false;				//   Return we failed
    }
  sub1 = c1;					// Set 1st component index
  sp1  = s1;					// Set 1st spin index
  sub2 = c2;					// Set 2nd component index
  sp2  = s2;					// Set 2nd spin index
  return true;
  }

//_________________________________________________________________________________
// iii                   CLASS SPIN MAP CHECKING FUNCTIONS
//_________________________________________________________________________________


bool SpinMap::Check(bool warn) const
  { return Check(sub1, sp1, sub2, sp2, warn); }

bool SpinMap::Check(int c1, int s1, int c2, int s2, bool warn) const
  {
  int MaxComp=100;
  int MaxSpin=100;
  if(c1<0 || c1>MaxComp) 
    { if(warn) SMerror(21, Gdec(c1), 1); return false; }
  if(c2<0 || c2>MaxComp) 
    { if(warn) SMerror(21, Gdec(c2), 1); return false; }
  if(s1<0 || s1>MaxSpin) 
    { if(warn) SMerror(21, Gdec(s1), 1); return false; }
  if(s2<0 || s2>MaxSpin) 
    { if(warn) SMerror(21, Gdec(s2), 1); return false; }
  if(c1==c2 && s1==s2)
    { if(warn) SMerror(33); return false; }
  return true;
  }

//_________________________________________________________________________________
// A           SPIN MAP CONSTRUCTORS, ASSIGNMENT, AND DESTRUCTOR
//_________________________________________________________________________________
 
SpinMap::SpinMap()  { sub1=0; sp1=0; sub2=0; sp2=0; }

SpinMap::SpinMap(int subA, int spA, int subB, int spB)
  {
  sub1 = subA;                          // Copy the 1st component
  sp1  = spA;                           // Copy the 1st component spin
  sub2 = subB;                          // Copy the 2nd component
  sp2  = spB;                           // Copy the 2nd component spin
  }

SpinMap::SpinMap(const SpinMap& Sp)
  {
  sub1 = Sp.sub1;
  sub2 = Sp.sub2;
  sp1 = Sp.sp1;
  sp2 = Sp.sp2;
  }

SpinMap& SpinMap::operator = (const SpinMap& Sp)
  {
  sub1 = Sp.sub1;                       // Copy the 1st component
  sp1  = Sp.sp1;                        // Copy the 1st component spin
  sub2 = Sp.sub2;                       // Copy the 2nd component
  sp2  = Sp.sp2;                        // Copy the 2nd component spin
  return *this;
  }

SpinMap::SpinMap(const std::string& SSP)

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
  std::string sval = SSP;			// Copy of input string
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
// B                          SPIN MAP ACCESS FUNCTIONS
//_________________________________________________________________________________
 
int SpinMap::Sub1()  const { return sub1; }
int SpinMap::Sub2()  const { return sub2; }
int SpinMap::Spin1() const { return sp1;  }
int SpinMap::Spin2() const { return sp2;  }
 
// ____________________________________________________________________________
// C                         SPIN MAP INPUT FUNCTIONS
// ____________________________________________________________________________

        // Input                SMap	: A spin map (this)
        //                      filename: Input filename
        //                      pset    : Input parameter set
        //                      idx     : Exchange process index
        //                      mdx     : Spin mapping index
        //                      warn    : Warning output level
        //                                      0 = no warnings
        //                                      1 = warnings
        //                                     >1 = fatal warnings
        // Output               TF      : Spin map is filled with
        //                                parameters read from file
        //                                TRUE if read is successful
        // Note                         : The file should be an ASCII file
        //                                containing recognized parameters

bool SpinMap::read(const std::string &filename, int idx, int mdx, int warn)
   {
   ParameterSet pset;                   // A GAMMA parameter set
   if(!pset.read(filename, warn?1:0))   // Try and read in pset
    {                                   // If we cannot read the file then
    if(warn)                            // we'll issue warnings as desired
      {
      SMerror(1, filename);             //      Problems with file
      if(warn>1) SMfatal(14);           //      Can't read from input file
      else       SMerror(14,1);         //      Fatal if warn big enough
      }
    return false;                       // Well flag we didn't read & exit
    }
  return read(pset, idx, mdx, warn);	// Fill up spin_sys with parameters
  }

bool SpinMap::read(const ParameterSet& pset, int idx, int mdx, int warn)
  {
  bool TF = setSM(pset,idx,mdx,warn?true:false);	// Use overload to read
  if(!TF)                                       // If setSM didn't handle
    {                                           // the process read from pset
    if(warn)                                    // then we'll issue warnings
      {                                         // & maybe even stop program
      SMerror(8, 1);                            //   Problems with pset
      if(warn>1) SMfatal(30);                   //   Can't read from pset
      else       SMerror(30, 1);                //   Fatal if warn big enough
      }
    }
  return TF;
  }

//_________________________________________________________________________________
// D                              SPIN MAP OUTPUT
//_________________________________________________________________________________

        // Input                SMap	: A spin map (this)
        //                      ostr    : Output stream
        // Output               none    : Spin pairing writtn to output stream
        ///F_list print                 - Write system to output stream 

void SpinMap::print() const
  {
  std::cout << "Component " <<sub1 << " spin #" <<sp1
            <<" --- to --- component ";
  std::cout << sub2 << " spin #" << sp2 << "\n";
  }
 
std::ostream& SpinMap::print(std::ostream& ostr) const
  {
  ostr << "Component "<< sub1 << ", Spin # " << sp1
       << " <---> "
       << "Component "<< sub2 << ", Spin # " << sp2;
  return ostr;
  }
 
std::ostream& operator<< (std::ostream& ostr, const SpinMap &Sp)
  { return Sp.print(ostr); }

#endif								// SpinMap.cc
