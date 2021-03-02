/* MultiAux.h ***************************************************-*-c++-*-
**									**
**                              G A M M A				**
**									**
**      Multiple System Auxiliary Classes 	    Interface		**
**									**
**      Copyright (c) 1995 						**
**      Nikolai Skrynnikov						**
**      Dr. Scott A. Smith						**
**      Copyright (c) 1996                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**									**
**      $Header: $
**									**
**************************************************************************

**************************************************************************
**									**
** This file contains auxiliary classes for mulitple spin systems.	**
** Specifically it handles the mapping of a spin in one system		**
** onto the spin of another system, the systems involved being		**
** those which are part of multi_sys.					**
**									**
** Currently, this contains the following classes:			**
**									**
** process: Tracks a single, specific, non-mutual exchange process.	**
**          This includes which components (dynamic systems) in a	**
**          multi_sys system are in exchange, the exchange rate, and	**
**          the spin <--> spin mappings that exist in the exchange	**
**									**
** spin_pair: Tracks 4 integers.  The first two are a component and	**
**            spin index as are the latter two.  The component indices	**
**            index particular spin systems in a variable of type	**
**            multi_sys.  The spin indices index particular spins in	**
**            their respective systems (or multi_sys components).	**
**									**
*************************************************************************/

#ifndef   multi_aux_h_			// Is this file already included?
#  define multi_aux_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <iostream>			// Knowledge of cout & NULL
#include <string>			// Include libstdc++ string
#include <vector>			// Include libstdc++ STL vectors
#include <Basics/ParamSet.h>		// Inlcude parameter sets

/*****************************************************************************/
/*****************************************************************************/
/*                                CLASS SPIN_PAIR			     */
/*****************************************************************************/
/*****************************************************************************/


class spin_pair
  {
  public:

    int sub1;					// First component
    int sp1;					// Spin of first component
    int sub2;					// Second component
    int sp2;					// Spin of second component

//_____________________________________________________________________________
// A	                       SPIN PAIR CONSTRUCTORS
//_____________________________________________________________________________


MSVCDLC inline spin_pair() {};

        // Input                spair   : A spin pairing (this)
        // Output               void	: A NULL spin pairing is created


MSVCDLC inline spin_pair(int, int, int, int);

        // Input                spair   : A spin pairing (this)
	//			c1,s1   : 1st component, spin involved
	//			c2,s2   : 2nd component, spin involved
        // Output               void	: A spin pairing is created


MSVCDLC spin_pair(const spin_pair& Sp);

        // Input                Sp	: A spin pair
        // Output               void	: New spin pair identical to Sp


MSVCDLC spin_pair(const std::string& SSP);   
 
        // Input                SSP     : A String defining a spin pair
        // Output               void    : The process is constructed 


MSVCDLL inline spin_pair& operator = (const spin_pair& Sp);

        // Input                spair   : A spin pairing (this)
	//			Sp      : A second spin pairing
        // Output               void	: Assigns Sp to spair


MSVCDLC inline ~spin_pair() {};

        // Input                spair   : A spin pairing (this)
        // Output               void	: Spin pair is destructed


//_____________________________________________________________________________
// B	                  SPIN PAIR ACCESS FUNCTIONS
//_____________________________________________________________________________

MSVCDLL int Sub1() const;
MSVCDLL int Sub2() const;

        // Input                spair   : A spin pairing (this)
        // Output               sub#	: The component index of spin #

MSVCDLL int Spin1() const;
MSVCDLL int Spin2() const;

        // Input                spair   : A spin pairing (this)
        // Output               sp#	: The index of spin #

//_____________________________________________________________________________
// C		                  SPIN PAIR OUTPUT
//_____________________________________________________________________________


MSVCDLL inline void print() const;

        // Input                spair   : A spin pairing (this)
 

MSVCDLL std::ostream& print(std::ostream& ostr) const;

        // Input                spair   : A spin pairing (this)
        //                      ostr	: Output stream
        // Output               none	: Spin pairing writtn to output stream
        ///F_list print			- Write system to output stream


MSVCDLL friend std::ostream& operator<< (std::ostream& ostr, const spin_pair &Sp);

        // Input                ostr	: Output stream;
        // 			spair   : A spin pairing (this)
        // Output			: Modifies output stream
        ///F_list <<			- Standard Output

};


/****************************************************************************/
/****************************************************************************/
/*                                CLASS PROCESS				    */
/****************************************************************************/
/****************************************************************************/


class process
  {
  public:

  double krate;				// Exchange rate
  int* comp_lhs;			// Components on the lhs of reaction
  int n_lhs;				// Number on the lhs of reaction
  int* comp_rhs; 			// Components on the rhs of reaction
  int n_rhs;				// Number on the rhs of the reaction
  int npair;				// Number of connected spins
  spin_pair* link;			// Array of connected spins


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

bool getExch(const ParameterSet& pset, int idx,
                                           std::string& exch, bool warn=true) const;

//---------------------------------------------------------------------------------
//-------------- Parse The Components Involved In The Exchange Process ------------
//---------------------------------------------------------------------------------

/* Here we assume we have the string from parameter Exch, e.g. (0<=>1+2). We need
   to parse out the integer values on each side of <=>, but lying between ().
   Also we must keep track of which are on the left and which are on the right.  */

bool parseExch(std::string& Exval,
                std::vector<int>& lhs, std::vector<int>& rhs, bool warn=true) const;

//---------------------------------------------------------------------------------
//------------------------ Read in the Process Exchange Rate ----------------------
//---------------------------------------------------------------------------------

// This will read a parameter such as    Kex_nm(0)  (1) : 600.0 - rate

bool getRate(const ParameterSet& pset, int idx,
                                               double& rate, bool warn=true) const;


//---------------------------------------------------------------------------------
//---------------------- Read/Set The Entire Exchange Process ---------------------
//---------------------------------------------------------------------------------

bool getXP(const ParameterSet& pset, int idx, bool warn=true) const;
bool setXP(const ParameterSet& pset, int idx, bool warn=true) const;


//_________________________________________________________________________________
// A		      CLASS PROCESS CONSTRUCTORS AND DESTRUCTORS
//_________________________________________________________________________________

//---------------------------------------------------------------------------------
//                           Simple Constructors
//---------------------------------------------------------------------------------

MSVCDLC inline process();
MSVCDLC inline process(int N_lhs, int N_rhs);
       process(const process& proc);

//---------------------------------------------------------------------------------
//                           Simple Constructors
//---------------------------------------------------------------------------------


MSVCDLC process(std::string& PROC, double Kex=0, int maxcomp=20);

        // Input                pro     : Process (this)
        //                      PROC    : String for process def.
	//			Kex     : Exchange rate (1/sec)
        //                      maxcomp : Maximum possible components
        // Output               void    : The process is constructed


//---------------------------------------------------------------------------------
//                    Construction From Parameter Set
//---------------------------------------------------------------------------------

MSVCDLC process(const ParameterSet& pset, int ip=-1, int warn=2);
 
//---------------------------------------------------------------------------------
//                         Assignment and Destruction
//---------------------------------------------------------------------------------

MSVCDLL process& operator=(const process& pr);
MSVCDLC      ~process();

        // Input                pro     : Process (this) 
        // Output               void	: The process is destructed


MSVCDLL void intra_default(int ic1, int ic2, int nspins, double k);

        // Input                none	: 
        // Output               pro	: A process is returned


//________________________________________________________________________________
// B	               Class Process Exchange Rate Access
//________________________________________________________________________________


MSVCDLL inline double get_k() const;

        // Input                pro	: A process (this)
        // Output               Kexch	: Process exchange rate


MSVCDLL inline void set_k(double k);

        // Input                pro	: A process (this)
	//			k       : An exchange rate
        // Output               Kexch	: Process exchange rate set to k

//________________________________________________________________________________
// C		       Class Process Component Index Access
//________________________________________________________________________________


MSVCDLL int lhsindex(int comp);      

        // Input                pro	: A process (this)
	//			comp    : A lhs component
        // Output               ic	: Index of component 


MSVCDLL int rhsindex(int comp);      

        // Input                pro	: A process (this)
	//			comp    : A rhs component
        // Output               ic	: Index of rhs component 

//____________________________________________________________________________ 
// D                       CLASS PROCESS SPIN QUERIES
//____________________________________________________________________________ 
  

MSVCDLL int mixes(int ic1, int ic2);

        // Input                pro     : A process
        //                      ic1     : Spin index 1
        //                      ic2     : Spin index 2
        // Output               TF      : True if both spins are
        //                                involved in the process
 

MSVCDLL int involves(int ic1, int lr=0);

        // Input                pro     : A process
        //                      ic1     : Spin index 1
        //                      lr      : Flag to check left and/or right
        //                                      0 = check left & right (def)
        //                                     >0 = check right only
        //                                     <0 = check left only 
        // Output               TF      : True if spin is involved in
        //                                the right and/or left process

//________________________________________________________________________________
// E		         CLASS PROCESS SPIN PAIR ACCESS
//________________________________________________________________________________


MSVCDLL int pairs() const;

        // Input                pro     : A process (this)
        // Output               npair   : Number of spin pairs defined
	//				  in the process


MSVCDLL void add_pair(spin_pair);

        // Input                pro     : A process (this)
	//			SP      : A spin pair
        // Output               void	: Spin pair is set as exchanging
	//				  in the process


MSVCDLL spin_pair get_pair(int) const;

        // Input                pro     : A process (this)
	//			int     : A spin pair index
        // Output               SP  	: Spin pair is returned


MSVCDLL int mapped(int c1, int s1, int c2, int s2) const;

        // Input                pro     : A process (this)
	//			int     : A spin pair index
        // Output               TF  	: True if spins are mapped
	//				  false if not


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

        // Input                pro   : A process
        //                      ostr  : Output stream
        //                      full  : Print amount flag
        //                                 0 = don't print spin mappings(def)
        //                                !0 = print individual spin mappings
        // Output               ostr  : The output stream  is returned
        //                              with the process added

MSVCDLL std::ostream& print(std::ostream& ostr, int full=0) const;
MSVCDLL friend std::ostream& operator<< (std::ostream& ostr, const process& pro);

};					// End of class process



/*********************************************************************************/
/*********************************************************************************/
/*                   MULTI_SYS AUXILIARY CLASSES INLINE FUNCTIONS		 */
/*********************************************************************************/
/*********************************************************************************/

//_________________________________________________________________________________
//_________________________________________________________________________________
//                          Class Multi_Pair Inlines
//_________________________________________________________________________________
//_________________________________________________________________________________

//_________________________________________________________________________________
// A           SPIN PAIR CONSTRUCTORS, ASSIGNMENT, AND DESTRUCTOR
//_________________________________________________________________________________
 

inline spin_pair::spin_pair(int subA, int spA, int subB, int spB)
  {
  sub1 = subA; 				// Copy the 1st component
  sp1  = spA;				// Copy the 1st component spin
  sub2 = subB;				// Copy the 2nd component
  sp2  = spB;				// Copy the 2nd component spin
  }


inline spin_pair& spin_pair::operator = (const spin_pair& Sp)
  {
  sub1 = Sp.sub1; 			// Copy the 1st component
  sp1  = Sp.sp1; 			// Copy the 1st component spin
  sub2 = Sp.sub2; 			// Copy the 2nd component
  sp2  = Sp.sp2; 			// Copy the 2nd component spin
  return *this;
  } 

//_________________________________________________________________________________
// C                              SPIN PAIR OUTPUT
//_________________________________________________________________________________
 

inline void spin_pair::print() const
  {
  std::cout << "Component " <<sub1 << " spin #" <<sp1
            <<" --- to --- component ";
  std::cout << sub2 << " spin #" << sp2 << "\n";
  }


//_________________________________________________________________________________
//_________________________________________________________________________________
//                            Class Process Inlines
//_________________________________________________________________________________
//_________________________________________________________________________________

//_________________________________________________________________________________
//                          Class Process Constructors
//_________________________________________________________________________________


inline process::process()

        // Input                pro	: A process (this)
        // Output               void	: A NULL process created
	// Note				: Connectivity is null
	// Note				: No connectivities are set

  {
  krate = 0.0;				// No exchange rate
  comp_lhs = NULL;			// No lhs components
  n_lhs = 0;				// 
  comp_rhs = NULL;			// No rhs components
  n_rhs = 0;
  npair = 0;				// No connected spin pairs
  link = NULL;				// No array of spin pairs
  }


inline process::process(int N_lhs, int N_rhs)

        // Input                pro	: A process (this)
	//			N_lhs   : Number left side components
	//			N_rhs   : Number right side components
        // Output               void	: A NULL process created
	// Note				: No connectivities are set

  {
  krate = 0.0;				// Set exchange rate to zero
  comp_lhs = new int[N_lhs];		// Allocate space for left components
  n_lhs = N_lhs ;			// Set the number of left componenets
  comp_rhs = new int[N_rhs];		// Allocate space for right components
  n_rhs = N_rhs;			// Set the number of right components 
  npair = 0;				// No connected spin pairs
  link = NULL;				// No array of spin pairs
  }

//__________________________________________________________________________
//			Class Process Exchange Rate Access
//__________________________________________________________________________


inline double process::get_k() const { return (krate); }

        // Input                pro	: A process
        // Output               Kexch	: Process exchange rate

inline void process::set_k(double k) { krate = k; }

        // Input                pro	: A process
	//			k       : An exchange rate
        // Output               Kexch	: Process exchange rate set to k



#endif								// MultiAux.h

