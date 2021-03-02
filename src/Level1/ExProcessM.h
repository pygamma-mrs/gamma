/* ExProcessM.h *************************************************-*-c++-*-
**									**
**                              G A M M A				**
**									**
**      Mutual Exchange Process                       Interface		**
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
** This class defines a single, specific, mutual exchange process. A	**
** mutual exchange process contains an exchange rate (1/sec) as well	**
** as an array of spin indices that are in exchange. The latter array	**
** will have at least 2 spins and indicates how the spins exchange.	**
** Spin i is moving into spin i+1 where i = [0,N-2] for N spins. The	**
** lone exception is that the last spin is moving into the first. 	**
**									**
**                 Initial Spin     Final Spin		    		**
**                 ============     ==========				**
**                       i      -->    i+1	 (i = [0, N-2])		**
**                      N-1     -->     0				**
**									**
** Since this class exists to be used by other classes, it has limited	**
** functionality.							**
**									**
*************************************************************************/

#ifndef   GExProcM_h_			// Is this file already included?
#  define GExProcM_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <iostream>			// Knowledge of cout & NULL
#include <string>			// Include libstdc++ string
#include <vector>			// Include libstdc++ STL vectors
#include <Basics/ParamSet.h>		// Include parameter sets

class ExchProcM
  {
public:

  double           KRate;		// Exchange rate (1/sec)
  std::vector<int> Spins; 		// Spins involved in the exchange

// ________________________________________________________________________________
// i                       MUTUAL EXCHANGE PROCESS ERROR HANDLING
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
// ii                MUTUAL EXCHANGE PROCESS PARAMETER SET PARSING
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
                                      std::vector<int>& sps, bool warn=true) const;
bool getComps(const ParameterSet& pset, int idx,
                                      std::vector<int>& sps, bool warn=true) const;

//---------------------------------------------------------------------------------
//------------------------ Read in the Process Exchange Rate ----------------------
//---------------------------------------------------------------------------------

// This will read a parameter such as    Kex_nm(0)  (1) : 600.0 - rate

bool getRate(const ParameterSet& pset, int idx,
                                               double& rate, bool warn=true) const;

//---------------------------------------------------------------------------------
//---------------------- Read/Set The Entire Exchange Process ---------------------
//---------------------------------------------------------------------------------

bool getXP(const ParameterSet& pset, double& rate,
                             std::vector<int>& sps, int idx, bool warn=true) const;
bool setXP(const ParameterSet& pset, int idx, bool warn=true);

//_________________________________________________________________________________
// iii             MUTUAL EXCHANGE PROCESS CHECKING FUNCTIONS
//_________________________________________________________________________________

// The first function just checks the boundaries of the exchange components
// The 2nd function insures the exchange process integrity

bool CCheck(int comp, bool warn=true) const;
bool FCheck(          bool warn=true) const;

//_________________________________________________________________________________
// A		UTUAL EXCHANGE PROCESS CONSTRUCTORS AND DESTRUCTORS
//_________________________________________________________________________________

//---------------------------------------------------------------------------------
//                           Simple Constructors
//---------------------------------------------------------------------------------

MSVCDLC ExchProcM();
MSVCDLC ExchProcM(const ExchProcM& proc);

//---------------------------------------------------------------------------------
//                      Construction From Parameter Set
//---------------------------------------------------------------------------------

MSVCDLC ExchProcM(const ParameterSet& pset, int ip=-1, int warn=2);
 
//---------------------------------------------------------------------------------
//                         Assignment and Destruction
//---------------------------------------------------------------------------------

MSVCDLL ExchProcM& operator= (const ExchProcM& pr);
MSVCDLC      ~ExchProcM();

//________________________________________________________________________________
// B	               Class Exchange Process Access Functions
//________________________________________________________________________________

//--------------------------------------------------------------------------------
//  	                           Exchange Rate
//--------------------------------------------------------------------------------

        // Input                pro     : A mutual exchange process (this)
	//			k       : A mutual exchange rate (1/sec)
        // Output               void	: Exchange rate is set to k
 	//		     or double  : Exchange rate is returned

MSVCDLL double Kex() const;
MSVCDLL void   Kex(double k);

//--------------------------------------------------------------------------------
//                    Class Process Component Index Access
//--------------------------------------------------------------------------------

        // Input                pro     : An exchange process (this)
        //                      comp    : A component index
        // Output               ic      : Spin index of component

MSVCDLL int NComps()       const;
MSVCDLL int NSpins()       const;
MSVCDLL int Comp(int comp) const;

//____________________________________________________________________________ 
// C                       MUTUAL EXCHANGE PROCESS SPIN QUERIES
//____________________________________________________________________________ 
  
        // Input                pro     : A mutual exchange process
        //                      i       : Component index 1
        //                      j       : Component index 2
        // Output               TF      : True if components are
        //                                directly exchanging
        //                   or         : True if component i involved
        //                                in the echange at all
        // Note                         : For direct exchange the two
        //                                spin must be adjacent Spins

MSVCDLL bool mixes(int    i, int j) const;
MSVCDLL bool involves(int i)        const;

// ____________________________________________________________________________
// D                  MUTUAL EXCHANGE PROCESS INPUT FUNCTIONS
// ____________________________________________________________________________

        // Input                pro     : A mutual exchange process (this)
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

        // Input                pro   : A mutual exchange process
        //                      ostr  : Output stream
        //                      full  : Print amount flag
        //                                 0 = unused so fat
        //                                !0 = unused so far
        // Output               ostr  : The output stream  is returned
        //                              with the process added

MSVCDLL std::string ExchStr() const;
MSVCDLL std::ostream& print(std::ostream& ostr, int full=0) const;
MSVCDLL friend std::ostream& operator<< (std::ostream& ostr, const ExchProcM& pro);

};

#endif								// ExchProcM.h

