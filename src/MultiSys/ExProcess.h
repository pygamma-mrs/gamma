/* ExProcess.h **************************************************-*-c++-*-
**									**
**                              G A M M A				**
**									**
**      NonMutual Exchange Process                    Interface		**
**									**
**      Copyright (c) 2001 						**
**      Dr. Scott A. Smith						**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**									**
**      $Header: $
**									**
**************************************************************************

**************************************************************************
**									**
** This class defines a single single, specific, non-mutual exchange	**
** process. This includes which components (dynamic systems) in a	**
** multi_sys system are in exchange, the exchange rate, and the 	**
** spin <--> spin mappings that exist in the exchange.			**
**									**
*************************************************************************/

#ifndef   GExProc_h_			// Is this file already included?
#  define GExProc_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <iostream>			// Knowledge of cout & NULL
#include <string>			// Include libstdc++ string
#include <vector>			// Include libstdc++ STL vectors
#include <Basics/ParamSet.h>		// Include parameter sets
#include <MultiSys/SpinMap.h>		// Include spin mappings

class ExchProc
  {
public:

  double               KRate;		// Exchange rate (1/sec)
  std::vector<int>     LHSComps; 	// Components on L.H.S. of reaction
  std::vector<int>     RHSComps; 	// Components on R.H.S. of reaction
  std::vector<SpinMap> SpinMaps;	// Array of connected spins

// ________________________________________________________________________________
// i                       CLASS PROCESS ERROR HANDLING
// ________________________________________________________________________________

        // Input		eidx    : Error flag
	//			pname   : String included in error
        //			noret   : Flag for return (0=return)
        // Output		none    : Output process error message
        //				  Program execution stop (fatal)
 
void XPerror(int eidx,                           int noret=0) const;
void XPerror(int eidx, const std::string& pname, int noret=0) const;
volatile void XPfatal(int eidx)                               const;
volatile void XPfatal(int eidx, const std::string& pname)     const;

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

bool getExch(const ParameterSet& pset, int idx,
                                          std::string& exch, bool warn=true) const;
bool parseExch(std::string& Exval,
               std::vector<int>& lhs, std::vector<int>& rhs, bool warn=true) const;
bool getComps(const ParameterSet& pset, int idx,
               std::vector<int>& lhs, std::vector<int>& rhs, bool warn=true) const;

//---------------------------------------------------------------------------------
//------------------------ Read in the Process Exchange Rate ----------------------
//---------------------------------------------------------------------------------

// This will read a parameter such as    Kex_nm(0)  (1) : 600.0 - rate

bool getRate(const ParameterSet& pset, int idx,
                                               double& rate, bool warn=true) const;


//---------------------------------------------------------------------------------
//------------------------ Read in the Process Spin Mappings ----------------------
//---------------------------------------------------------------------------------

// This will be a parameter such as    Smap(0,0)  (2) : (0)0(1)0 - Mapping

bool getMappings(const ParameterSet& pset, int idx,
                                std::vector<SpinMap>& smaps, bool warn=true) const;

//---------------------------------------------------------------------------------
//---------------------- Read/Set The Entire Exchange Process ---------------------
//---------------------------------------------------------------------------------

bool getXP(const ParameterSet& pset, double& rate,
  std::vector<int>& lhsc, std::vector<int>& rhsc, std::vector<SpinMap>& smaps,
                                               int idx, bool warn=true) const;
bool setXP(const ParameterSet& pset, int idx, bool warn=true);

//_________________________________________________________________________________
// iii             CLASS EXCHANGE PROCESS CHECKING FUNCTIONS
//_________________________________________________________________________________

/* These functions check the boundaries of the exchange components and insure
   process integrity.                                                            */

bool CheckLHS(int comp, bool warn=true) const;
bool CheckRHS(int comp, bool warn=true) const;

//_________________________________________________________________________________
// A		      CLASS PROCESS CONSTRUCTORS AND DESTRUCTORS
//_________________________________________________________________________________

//---------------------------------------------------------------------------------
//                           Simple Constructors
//---------------------------------------------------------------------------------

MSVCDLC ExchProc();
MSVCDLC ExchProc(const ExchProc& proc);

//---------------------------------------------------------------------------------
//              Explicitly Defined Exchange Process Constructors
//---------------------------------------------------------------------------------

MSVCDLC ExchProc(const std::string& PROC, double Kex=0, int maxcomp=20);

        // Input                pro     : Process (this)
        //                      PROC    : String for process def.
	//			Kex     : Exchange rate (1/sec)
        //                      maxcomp : Maximum possible components
        // Output               void    : The process is constructed


//---------------------------------------------------------------------------------
//                      Construction From Parameter Set
//---------------------------------------------------------------------------------

MSVCDLC ExchProc(const ParameterSet& pset, int ip=-1, int warn=2);
 
//---------------------------------------------------------------------------------
//                         Assignment and Destruction
//---------------------------------------------------------------------------------

MSVCDLL ExchProc& operator=(const ExchProc& pr);
MSVCDLC      ~ExchProc();

        // Input                pro     : Process (this) 
        // Output               void	: The process is destructed

// sosi - not sure what these are for now....
MSVCDLC ExchProc(int N_lhs, int N_rhs);

MSVCDLL void intra_default(int ic1, int ic2, int nspins, double k);

        // Input                none	: 
        // Output               pro	: A process is returned


//________________________________________________________________________________
// B	               Class Exchange Process Access Functions
//________________________________________________________________________________

//--------------------------------------------------------------------------------
//  	                           Exchange Rate
//--------------------------------------------------------------------------------

        // Input                pro     : An exchange process (this)
	//			k       : An exchange rate (1/sec)
        // Output               void	: Exchange rate is set to k
 	//		     or double  : Exchange rate is returned

MSVCDLL double Kex() const;
MSVCDLL void   Kex(double k);

//--------------------------------------------------------------------------------
//                    Class Process Component Index Access
//--------------------------------------------------------------------------------

        // Input                pro     : An exchange process (this)
        //                      comp    : A L.H.S. or R.H.S. component
        // Output               ic      : Index of component

MSVCDLL int LHSComp(int comp) const;
MSVCDLL int RHSComp(int comp) const;

MSVCDLL int NCompsLHS() const;
MSVCDLL int NCompsRHS() const;

//____________________________________________________________________________ 
// D                       CLASS PROCESS SPIN QUERIES
//____________________________________________________________________________ 
  
        // Input                pro     : A process
        //                      comp 	: A component index
        //                      comp1	: A component index
        //                      lr      : Flag to check left and/or right
        //                                      0 = check left & right (def)
        //                                     >0 = check right only
        //                                     <0 = check left only 
        // Output               TF      : True if comp & comp1 in exchange
        //                      	: True if comp in LHS of exchange
        //                      	: True if comp in RHS of exchange
        //                      	: True if comp involved in exchange

MSVCDLL bool mixes(int     comp, int comp1) const;
MSVCDLL bool CompInLHS(int comp)            const;
MSVCDLL bool CompInRHS(int comp)            const;
MSVCDLL bool involves(int  comp, int lr=0)  const;

//________________________________________________________________________________
// E		         CLASS PROCESS SPIN PAIR ACCESS
//________________________________________________________________________________

        // Input                pro     : A process (this)
	//			SP      : A spin pair
	//			int     : A spin pair index
        // Output               void	: Spin pair is set as exchanging
	//				  in the process
        // Output               TF  	: True if spins are mapped
	//				  false if not

MSVCDLL       int      NSpinMaps() const;
MSVCDLL const SpinMap& SMap(int i) const;
MSVCDLL       bool     SMap(int comp1, int sp1, int& comp2, int& sp2) const;
MSVCDLL       void     add_pair(SpinMap);
MSVCDLL       bool     mapped(int comp1, int s1, int comp2, int s2) const;
MSVCDLL       bool     mapped(int comp1,         int comp2) const;

//______________________________________________________________________________
// F                          CLASS PROCESS CONNECTIONS
//______________________________________________________________________________


MSVCDLL void mapping(const std::string& spair);

        // Input                pro     : Process (this) 
        //                      spair   : String for spin pair def.
        // Output               void    : The process is modified
        //                                to contain a spin pairing as
        //                                suggested in the input string
        // Note                         : It is assumed the input string
        //                                has the form
        //
        //                                      (C1)S1(C2)S2
        //
        //                                where C1,S1,C2,S2 are integers

// ____________________________________________________________________________
//                       EXCHANGE PROCESS INPUT FUNCTIONS
// ____________________________________________________________________________

        // Input                pro     : An exchange process (this)
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

MSVCDLL bool read(const std::string&  filename, int idx=-1, int warn=2);
MSVCDLL bool read(const ParameterSet& pset,     int idx=-1, int warn=2);

//_____________________________________________________________________________
// G                        EXCHANGE PROCESS OUTPUT
//_____________________________________________________________________________

//-----------------------------------------------------------------------------
//              Functions To Modularize Exchange Process Output
//-----------------------------------------------------------------------------

/* These next two function return a string to indicate alphabetically which
   components are involved in the exchange process. The string returned will
   look something like
                                A + D + E + G

   where the first component has label A, the second label B, etc.           */

MSVCDLL static char Label(int i);
MSVCDLL std::string LHSStr()     const;
MSVCDLL std::string RHSStr()     const;

/* The next function returns a vector of strings that indicate the spin mappings
   that are active in the exchange process. Each string will involve one spin
   of a LHS component exchanging with one spin of the RHS component. Each string
   will be of the same length and appear something like

                                   A 0 <---> B 1

   where here A is the LHS component and 0 its spin, while B is a RHS component
   with 1 being its spin.                                                         */

MSVCDLL std::vector<std::string> SpinMapStrs() const;




        // Input                pro   : A process
        //                      ostr  : Output stream
        //                      full  : Print amount flag
        //                                 0 = don't print spin mappings(def)
        //                                !0 = print individual spin mappings
        // Output               ostr  : The output stream  is returned
        //                              with the process added

MSVCDLL std::ostream& print(std::ostream& ostr, int full=0) const;
MSVCDLL friend std::ostream& operator<< (std::ostream& ostr, const ExchProc& pro);

};

#endif								// ExchProc.h

