/* IntCSAVec.h **************************************************-*-c++-*-
**									**
** 	                         G A M M A				**
**			         				 	**
**	Shit Anisotropy Interactions Vector 		Interface	**
**								 	**
**	Scott Smith 							**
**      Copyright (c) 1996                      			**
**      National High Magnetic Field Laboratory				**
**      1800 E. Paul Dirac Drive					**
**      Tallahassee Florida, 32306-4005					**
**								 	**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**							 		**
** Description						 		**
**							 		**
** This class maintains a vector of rank2 SA interactions (from class	**
** IntCSA).  The class allows users to manipulate an array of such	**
** interactions simultaneously.						**
**						 			**
*************************************************************************/

///Chapter Class Rank 2 Shift Anisotropy Interactions Vector (IntCSAVec)
///Section Overview
///Body    The class 
///Section Available Shift Anisotropy Interactions Vector Functions

#ifndef   IntCSAVec_h_			// Is this file already included?
#  define IntCSAVec_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#if defined(_MSC_VER)			// If using MSVC++ then we
 #pragma warning (disable : 4786)       // Kill STL namelength warnings
#endif 

#include <IntRank2/IntCSA.h>		// Include SA interactions
#include <Basics/ParamSet.h>		// Include parameter sets
#include <string>			// Include libstdc++ strings
#include <vector>			// Include libstdc++ STL vectors

class IntCSAVec : public std::vector<IntCSA>
  {


// ----------------------------------------------------------------------------  
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                 CLASS SA INTERACTION VECTOR ERROR HANDLING
// ____________________________________________________________________________
 
	// Input		ICV	: CSA interaction vector (this)
        // 			eidx	: Error index
        //                      pname	: Additional string for error
        //                      noret	: Flag for return (0=linefeed)
        // Output               none	: Error message
	//				  Program execution stopped

void ICVerror(int eidx, int noret=0) const;
void ICVerror(int eidx, const std::string& pname, int noret=0) const;
volatile void ICVfatal(int eidx) const;
// ____________________________________________________________________________
// ii               CSA INTERACTIONS VECTOR SETUP FUNCTIONS
// ____________________________________________________________________________

/* These functions determine specific aspects of a shift anisotropy
   interactions vector from a specified parameter set. As the various
   interaction parameters tend to interact, the functions MUST be private
   because their misuse could produce an inconsistent interactions vector.

   The goal here is quite simple. We need to determine a series of shift
   anisotropy interactions from parameters found in an input GAMMA parameter
   set. Each interaction will be defined by unique parameters in the set that
   specify the values: { Iqn, PPM, CSA,eta, alpha, beta, gamma } Complexity
   arises because we need find a series of such interactions and because
   the interactions allow for a variety of parameters for their definition.

   In this case, since shift anisotropy interactions are singly indexed, there
   is a one to one correspondence between the interaction index and the vector
   index. We simply assume that the first interaction has an index of 0 and
   that the rest are a series spanning 0, 1, 2, .... until the first not
   found.
           Input                ICV     : CSA interaction vector (this)
                                pset    : A parameter set
                                idx     : Parameter prefix index
                                warn    : Flag for warnings
           Output               none    : System dipolar interactions are
                                          set from parameters in pset

   We Use Only Parameters With Prefix [#] (unles idx=-1). As such we off this
   prefix from all parameter names prior to looking for the interactions.
   This is done before any calls to class IntCSA which doesn't use prefixes. */
     
bool setCIV(const ParameterSet& pset,int idx=-1,bool warn=true);

// ____________________________________________________________________________
// iii         SHIFT ANISOTROPY INTERACTIONS VECTOR CHECKING FUNCTIONS
// ____________________________________________________________________________

/* These functions insure vector integrity and that access is not beyond
   the vector boundaries

           Input                ICV     : CSA interaction vector (this)
                                spin    : A CSA index
                                warn    : Warning level
           Output               TF      : Returns TRUE if spin is a valid
                                          interaction index, FALSE if not    */

bool CheckCI(int spin, int warn=0) const;

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

// ____________________________________________________________________________ 
// A              SA INTERACTION VECTOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________ 
///Center Basic Functions

// ---------------------------------------------------------------------------- 
//                  Simple Constructors That Don't Do Much
// ---------------------------------------------------------------------------- 

/* These functions are inherited from the ANSI standard "vector" class in the
   C++ library libstdc++.  I am listing them here so that I & other users
   don't have to keep looking them up all the time.                          */

MSVCDLC IntCSAVec();
//IntCSAVec(int N)                     Vector w/ N Interactions
//IntCSAVec(int N, const IntCSA& CI)   Vector w/ N pars
//IntCSAVec(const IntCSAVec& CSVec)    Vector copy of CSVec
//IntCSAVec assign(N)                  Assign N element

// ----------------------------------------------------------------------------
//                      Constructors Using Parameter Sets 
// ----------------------------------------------------------------------------

/* These functions allow users to fill up the shift anisotropy interactions
   vector from a GAMMA parameter set (external ASCII file). Variants enable
   higher classes such as spin systems to generate the vector when they have
   a knowledge of spin isotope types.                                        */
   
        // Input                ICV     : CSA interaction vector (this)
        //                      ns      : Number of spins in vector
        //                      pset    : Parameter set
        //                      idx     : Index for CSA vector
        //                                i.e [idx] is parameter name prefix
        //                      warn    : Warning level
        // Note                         : If idx is negative then no
	//				     0 - no warnings
	//				     1 - warn if e- is isotope
	//				     2 - fatal if e- is isotope
        //                                parameter prefix is used
	//				  else prefix is "[idx]"
	// Note				: If ns is not specified as an
	//				  argument the function will attempt
	//				  to read {[idx]CSA(i)} from i=0 
	//				  to i=ns-1 where ns is the number

MSVCDLC IntCSAVec(const ParameterSet& pset, int idx=-1, int warn=1);
MSVCDLC IntCSAVec(const std::vector<Isotope>& Isos,
                                         const ParameterSet& pset, int warn=2);

// ----------------------------------------------------------------------------
//               Here Be The Assignment Operator & Destructor
// ----------------------------------------------------------------------------

MSVCDLL void operator= (const IntCSAVec &ICV);
MSVCDLC      ~IntCSAVec ();

// ____________________________________________________________________________ 
// B                CSA INTERACTION ACCESS & MANIPULATIONS
// ____________________________________________________________________________ 

// ---------------------------------------------------------------------------- 
//                   Generic CSA Value Access Functions
// ---------------------------------------------------------------------------- 

	// Input		ICV	: CSA interaction vector
	// 			spin	: CSA index
	// 			val	: CSA interaction value
	//				  <=0: CSA coupling of spin (Hertz)
	//				    1: Asymmetry [0, 1]
	//				    2: Theta (degrees, down from +z)
	//				    3: Phi   (degrees, over from +z)
	// Output		none	: Get/Set CSA interaction value
	//				  for specified spin
 

MSVCDLL void   CValue(int spin, double val, int type);
MSVCDLL double CValue(int spin, int type) const;


// ---------------------------------------------------------------------------- 
//                         CSA Coupling Constants
// ---------------------------------------------------------------------------- 

	// Input		ICV	: CSA interaction vector
	// 			spin	: CSA index
	// 			dcc	: CSA coupling constant (PPM)
	// Output		none	: Get/Set spin CSA
	// Note				: Defined in class IntCSA as 1.5 times
	//			          the CSA tensor delzz value

MSVCDLL void   CSA(int  spin, double CSA);
MSVCDLL void   delz(int spin, double dcc);
MSVCDLL double CSA(int  spin) const;
MSVCDLL double delz(int spin) const;

// ---------------------------------------------------------------------------- 
//                        CSA Asymmetry Values
// ---------------------------------------------------------------------------- 

	// Input		ICV	: CSA interaction vector
	// 			spin	: CSA index
	// 			ceta	: CSA interaction asymmetry
	// Output		none	: Get/Set CSA spin asymmetry
	// Note				: Defined in class IntCSA between [0,1]

MSVCDLL void   eta(int spin, double ceta);
MSVCDLL double eta(int spin) const;

 
// ---------------------------------------------------------------------------- 
//               CSA Theta Orientation (Down From +z Axis)
// ---------------------------------------------------------------------------- 

	// Input		ICV	: CSA interaction vector
	// 			spins	: CSA index
	// 			ctheta	: CSA interaction angle (deg)
	// Output		none	: Get/Set CSA spin theta angle
	// Note				: Defined class IntCSA between [0,180]

MSVCDLL void   theta(int spin, double ctheta);
MSVCDLL double theta(int spin) const;

// ---------------------------------------------------------------------------- 
//               CSA Phi Orientation (Over From +x Axis)
// ---------------------------------------------------------------------------- 

	// Input		ICV	: CSA interaction vector
	// 			spins	: CSA index
	// 			cphi	: CSA interaction angle (deg)
	// Output		none	: Get/Set CSA spin phi angle
	// Note				: Defined in IntCSA between [0,360]

MSVCDLL void   phi(int spin, double phi);
MSVCDLL double phi(int spin) const;

// ---------------------------------------------------------------------------- 
//                        Full CSA Interaction
// ---------------------------------------------------------------------------- 

MSVCDLL IntCSA& operator() (int spins);                                          

        // Input                ICV	: CSA interaction vector (this)
        //                      spins	: Interaction index
        // Ouput                CI	: The i'th CSA interaction in ICV
        // Note				: Returns a reference to the interaction


MSVCDLL IntCSA get(int spins) const;

	// Input		ICV	: CSA interaction vector
	// 			spins	: CSA index
	// Output		CI	: Return rank 2 CSA interaction


// ----------------------------------------------------------------------------
//                 Other CSA Interaction Vector Info
// ----------------------------------------------------------------------------
 
 
// int  size() const;		INHERITED
MSVCDLL bool nonzero() const;
 
// ____________________________________________________________________________
// C                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//     Functions To Make A Parameter Set From A CSA Interaction Vector
// ----------------------------------------------------------------------------

// Note that the preferred means of specifying CSA interactions when there
// is no knowledge of spin isotope types nor spin coordinates are used for
// filling the parameter set.  For base parameters of insividual interactions
// see class IntCSA. Additional parameters for the interaction vector will 
// be defined in the IntRank2.cc file in this section.


	// Input		ICV	: CSA interaction vector
        //                      pset	: Parameter set
        // Output               pset	: Parameter set with only
        //                                CSA interaction parameters
 
	// Input		ICV	: CSA interaction vector
	//  			pset	: Parameter set
	// Output		pset	: Parameter set with all vector
	//			          CSA interaction parameters added
 
         

        // Input                ICV     : CSA interaction vector
        //                      pset    : Parameter set
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        // Output               void    : Vector parameters are
        //                                are added ot the parameter set
        //                                with interaction index idx
        // Note                         : The parameter names & types
        //                                output here MUST match those used
        //                                in setting the vector up
        //                                from parameters sets

MSVCDLL operator ParameterSet( ) const;
MSVCDLL friend void operator+= (ParameterSet& pset, const IntCSAVec &ICV);
MSVCDLL void PSetAdd(ParameterSet& pset, int idx=-1) const;



// ----------------------------------------------------------------------------
//  Functions To Output CSA Interaction Vector To ASCII From A Parameter Set
// ----------------------------------------------------------------------------
 
        // Input                ICV     : CSA interaction vector (this)
        //                      filename: Output file name
        //                      ofstr   : Output file stream
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning level
        // Output               none    : CSA interaction vector is written as
        //                                a parameter set to file filename
        //                                or to output filestream ofstr

MSVCDLL bool write(const std::string& filename,int idx=-1,int warn=2) const;
MSVCDLL bool write(std::ofstream&     ofstr,   int idx=-1,int warn=2) const;

// ____________________________________________________________________________ 
// D                             INPUT FUNCTIONS
// ____________________________________________________________________________ 

// ----------------------------------------------------------------------------
//        Direct Read of Vector From An ASCII File Or A Parameter Set
// ----------------------------------------------------------------------------

	// Input		ICV	: CSA interaction vector (this)
	// 			filename: Input filename
        //                      pset    : Parameter set
        //                      idx     : Parameter index value
	// Output		none	: CSA interaction vector filled with
	//				  parameters read from file
        //                                or with parameters read from pset
 
MSVCDLL bool read(const std::string&  filename, int idx=-1, int warn=true);
MSVCDLL bool read(const ParameterSet& pset,     int idx=-1, int warn=true);

// ----------------------------------------------------------------------------
//       Interactive Read of CSA Interaction Vector From An ASCII File
// ----------------------------------------------------------------------------

MSVCDLL std::string ask_read(int argc, char* argv[], int argn);

	// Input		ICV    : CSA interaction vector (this)
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


// ____________________________________________________________________________ 
// E                           OUTPUT FUNCTIONS
// ____________________________________________________________________________ 


MSVCDLL std::ostream& print(std::ostream& ostr, int full=0) const;

	// Input		ICV	: CSA interaction vector (this)
	// 			ostr	: Output stream
	//			full    : Flag for long vs short output
	// Output		non	: CSA interaction vector parameters
	//			          sent to the output stream


MSVCDLL friend std::ostream& operator<< (std::ostream& out, const IntCSAVec& ICV);

	// Input		ICV	: CSA interaction vector
	// 			out	: Output stream
	// Output		none	: CSA interaction vector
	//				  parameters sent to output stream
};

#endif							// IntCSAVec.h
