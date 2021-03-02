/* IntGVec.h ****************************************************-*-c++-*-
**									**
** 	                         G A M M A				**
**			         				 	**
**	Electron G Interactions Vector 			Interface	**
**								 	**
**	Scott Smith 							**
**      Copyright (c) 2001                      			**
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
** This class maintains a vector of rank 2 electron G interactions (as	**
** defined in class IntG). The class allows users, or more importantly	**
** spin systems, to manipulate an array of such	interactions either in	**
** concert or individually.						**
**						 			**
*************************************************************************/

///Chapter Class Rank 2 Electron G Interactions Vector (IntGVec)
///Section Overview
///Body    The class 
///Section Available Electron G Interactions Vector Functions

#ifndef   IntGVec_h_			// Is this file already included?
#  define IntGVec_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <IntRank2/IntG.h>		// Include electron G interactions
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Basics/Isotope.h>		// Include spin isotopes
#include <string>			// Include libstdc++ strings
#include <vector>			// Include libstdc++ STL vectors

//typedef std::vector<IntG> stdvecG;

//class IntGVec : public stdvecG
class IntGVec : public std::vector<IntG>
  {

private:

// ----------------------------------------------------------------------------  
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                 ELECTRON G INTERACTION VECTOR ERROR HANDLING
// ____________________________________________________________________________
 

/*      Input			IGV	: Electron G interaction vector (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
                                pname   : Added error message for output
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

void IGVerror(int eidx,                           int noret=0) const;
void IGVerror(int eidx, const std::string& pname, int noret=0) const;
volatile void IGVfatal(int eidx) const;

// ____________________________________________________________________________
// iii            ELECTRON G INTERACTIONS VECTOR SETUP FUNCTIONS
// ____________________________________________________________________________

/* These functions determine specific aspects of an electron g interactions
   vector from a specified parameter set. As the various interaction parameters
   tend to interact, the functions MUST be private because their misuse could
   produce an inconsistent interactions vector.

   The goal here is quite simple. We need to determine a series of electron g
   interactions from parameters found in an input GAMMA parameter set. Each
   interaction will be defined by unique parameters in the set that specify the
   values: { Iqn, g, gA, eta, alpha, beta, gamma }. Complexity arises because
   we need find a series of such interactions and because the interactions
   allow for a wide variety of parameters for their definition.

   In this case, since electron G interactions are singly indexed, there is a
   one to one correspondence between the interaction index and the vector
   index. We simply assume that the first interaction has an index of 0 and
   that the rest are a series spanning 0, 1, 2, .... until the first not
   found.
           Input                IGV     : G interactions vector (this)
                                pset    : A parameter set
                                idx     : Parameter prefix index
                                warn    : Flag for warnings
           Output               none    : Electron G interactions are
                                          set from parameters in pset

   We Use Only Parameters With Prefix [#] (unles idx=-1). As such we off this
   prefix from all parameter names prior to looking for the interactions.
   This is done before any calls to class IntG which doesn't use prefixes.  */

bool setGIV(const ParameterSet& pset, int idx=-1, bool warn=true);

// ____________________________________________________________________________
// iii            ELECTRON G INTERACTIONS VECTOR CHECKING FUNCTIONS
// ____________________________________________________________________________

/* These functions insure vector integrity and that access is not beyond
   the vector boundaries

           Input                IGV     : Electron G interaction vector (this)
                                spin    : A G index
                                warn    : Warning level
           Output               TF      : Returns TRUE if spin is a valid
                                          interaction index, FALSE if not    */

bool check_spin(int spin, int warn=0) const;

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

// ____________________________________________________________________________ 
// A           ELECTRON G INTERACTION VECTOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________ 

// ---------------------------------------------------------------------------- 
//                  Simple Constructors That Don't Do Much
// ---------------------------------------------------------------------------- 
 

/* These functions are inherited from the ANSI standard "vector" class in the
   C++ library libstdc++.  I am listing them here so that I & other users don't
   have to keep looking them up all the time.
*/

MSVCDLC IntGVec();					// Empty Interaction Vector
/*
   IntGVec(int N)                               Vector w/ N Interactions
   IntGVec(int N, const IntG& GI)               Vector w/ N pars
   IntGVec(const IntGVec& GVec)                 Vector copy of GVec
   IntGVec assign(N)                            Assign N element
   ~IntGVec()                                   Destructor of Vector         */

// ----------------------------------------------------------------------------
//                      Constructors Using Parameter Sets 
// ----------------------------------------------------------------------------
   
        // Input                IGV     : Electron G interaction vector (this)
        //                      ns      : Number of spins in vector
        //                      pset    : Parameter set
        //                      idx     : Index for Electron G vector
        //                                i.e [idx] is parameter name prefix
        //                      warn    : Warning level
        // Note                         : If idx is negative then no
        //                                parameter prefix is used
	//				  else prefix is "[idx]"
	// Note				: If ns is not specified as an
	//				  argument the function will attempt
	//				  to read {[idx]Electron G(i)} from i=0 
	//				  to i=ns-1 where ns is the number

//IntGVec(const         ParameterSet& pset, int indx=-1, int warn=1);
//IntGVec(int ns, const ParameterSet& pset, int indx=-1, int warn=1);

// ----------------------------------------------------------------------------
//          This Constructor Supports Generation From A Spin System
// ----------------------------------------------------------------------------

/* These functions allow users to fill up the electron G interactions vector
   from a GAMMA parameter set (external ASCII file). Variants enable higher
   classes such as spin systems to generate the vector when they have a
   knowledge of spin isotope types.

           Input                IGV     : Electron G interaction vector (this)
                                pset    : Parameter set
                                idx     : Index for G vector
                                          i.e [idx] is parameter name prefix
                                warn    : Warning level
                                                0 - no warnings
                                                1 - warn if spin isn't e-
                                                2 - fatal if spin isn't e-
           Note                         : If idx is negative then no
                                          parameter indicies are used
           Note                         : Since this has no idea about
                                          the number of interactions it will
                                          count through G(i) until there
                                          is none found starting with i=0     */

MSVCDLC IntGVec(const ParameterSet& pset, int indx=-1, int warn=2);
MSVCDLC IntGVec(const std::vector<Isotope>& Isos,
                                         const ParameterSet& pset, int warn=2);

// ----------------------------------------------------------------------------
//               Here Be The Assignment Operator & Destructor
// ----------------------------------------------------------------------------

//void operator= (const IntGVec &IGV);

	// Input		IGV	: Electron G interaction vector (this)
	// Output		none	: Electron G interaction vector is
	//				  constructed equivalent to sys


//~IntGVec ();

	// Input		IGV	: Electron G interaction vector (this)
	// Output		none	: System IGV is destructed


// ____________________________________________________________________________ 
// B                Electron G INTERACTION ACCESS & MANIPULATIONS
// ____________________________________________________________________________ 

// ---------------------------------------------------------------------------- 
//                   Generic Electron G Value Access Functions
// ---------------------------------------------------------------------------- 

	// Input		IGV	: Electron G interaction vector
	// 			spin	: Electron G index
	// 			val	: Electron G interaction value
	//				  <=0: Electron G coupling of spin (Hertz)
	//				    1: Asymmetry [0, 1]
	//				    2: Theta (degrees, down from +z)
	//				    3: Phi   (degrees, over from +z)
	// Output		none	: Get/Set Electron G interaction value
	//				  for specified spin
 

MSVCDLL void   CValue(int spin, double val, int type);
MSVCDLL double CValue(int spin, int type) const;


// ---------------------------------------------------------------------------- 
//                         Electron G Coupling Constants
// ---------------------------------------------------------------------------- 

	// Input		IGV	: Electron G interaction vector
	// 			spin	: Electron G index
	// 			dcc	: Electron G coupling constant (PPM)
	// Output		none	: Get/Set spin Electron G
	// Note				: Defined in class IntElectron G as 1.5 times
	//			          the Electron G tensor delzz value

//void   G(int  spin, double Electron G);
//double G(int  spin) const;
MSVCDLL void   delz(int spin, double dcc);
MSVCDLL double delz(int spin) const;

// ---------------------------------------------------------------------------- 
//                        Electron G Asymmetry Values
// ---------------------------------------------------------------------------- 

	// Input		IGV	: Electron G interaction vector
	// 			spin	: Electron G index
	// 			ceta	: Electron G interaction asymmetry
	// Output		none	: Get/Set Electron G spin asymmetry
	// Note				: Defined in class IntElectron G between [0,1]

MSVCDLL void   eta(int spin, double ceta);
MSVCDLL double eta(int spin) const;

 
// ---------------------------------------------------------------------------- 
//               Electron G Theta Orientation (Down From +z Axis)
// ---------------------------------------------------------------------------- 

	// Input		IGV	: Electron G interaction vector
	// 			spins	: Electron G index
	// 			ctheta	: Electron G interaction angle (deg)
	// Output		none	: Get/Set Electron G spin theta angle
	// Note				: Defined class IntElectron G between [0,180]

MSVCDLL void   theta(int spin, double ctheta);
MSVCDLL double theta(int spin) const;

// ---------------------------------------------------------------------------- 
//               Electron G Phi Orientation (Over From +x Axis)
// ---------------------------------------------------------------------------- 

	// Input		IGV	: Electron G interaction vector
	// 			spins	: Electron G index
	// 			cphi	: Electron G interaction angle (deg)
	// Output		none	: Get/Set Electron G spin phi angle
	// Note				: Defined in IntElectron G between [0,360]

MSVCDLL void   phi(int spin, double phi);
MSVCDLL double phi(int spin) const;

// ---------------------------------------------------------------------------- 
//                        Full Electron G Interaction
// ---------------------------------------------------------------------------- 

MSVCDLL IntG& operator() (int spins);                                          

        // Input                IGV	: Electron G interaction vector (this)
        //                      spins	: Interaction index
        // Ouput                CI	: The i'th Electron G interaction in IGV
        // Note				: Returns a reference to the interaction


MSVCDLL IntG get(int spins) const;

	// Input		IGV	: Electron G interaction vector
	// 			spins	: Electron G index
	// Output		CI	: Return rank 2 Electron G interaction


// ----------------------------------------------------------------------------
//                 Other Electron G Interaction Vector Info
// ----------------------------------------------------------------------------
 
 
//int size() const;
 
        // Input                IGV   : Electron G interaction vector
        // Output               ns    : Number of interactions in vector
 
 
MSVCDLL int nonzero() const;
 
        // Input                IGV   : Electron G interaction vector
        // Output               TF    : True if any interactions with a
        //                              finite Electron G value


// ____________________________________________________________________________
// C                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//     Functions To Make A Parameter Set From A Electron G Interaction Vector
// ----------------------------------------------------------------------------

// Note that the preferred means of specifying Electron G interactions when there
// is no knowledge of spin isotope types nor spin coordinates are used for
// filling the parameter set.  For base parameters of insividual interactions
// see class IntElectron G. Additional parameters for the interaction vector will 
// be defined in the IntRank2.cc file in this section.

MSVCDLL operator ParameterSet( ) const;

	// Input		IGV	: Electron G interaction vector
        //                      pset	: Parameter set
        // Output               pset	: Parameter set with only
        //                                Electron G interaction parameters
 

MSVCDLL friend void operator+= (ParameterSet& pset, const IntGVec &IGV);

	// Input		IGV	: Electron G interaction vector
	//  			pset	: Parameter set
	// Output		pset	: Parameter set with all vector
	//			          Electron G interaction parameters added
 
         
MSVCDLL void PSetAdd(ParameterSet& pset, int idx=-1) const;

        // Input                IGV     : Electron G interaction vector
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

// ----------------------------------------------------------------------------
//      Output Electron G Interaction Vector To ASCII From A Parameter Set
// ----------------------------------------------------------------------------
 
        // Input                IGV     : Electron G interaction vector (this)
        //               THIS   filename: Output file name
        //            OR THIS   ofstr   : Output file stream
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning level
        // Output               none    : Electron G interaction vector written
        //                                as a parameter set either to file
	//				  filename or ourput stream ofstr
 
MSVCDLL int write(const std::string &filename, int idx=-1, int warn=2) const;
MSVCDLL int write(std::ofstream& ofstr,        int idx=-1, int warn=2) const;


// ____________________________________________________________________________ 
// D                             INPUT FUNCTIONS
// ____________________________________________________________________________ 

// ----------------------------------------------------------------------------
//        Direct Read of Vector From An ASCII File Or A Parameter Set
// ----------------------------------------------------------------------------

MSVCDLL bool read(const std::string &filename, int idx=-1, int warn=2);

	// Input		IGV	: Electron G interaction vector (this)
	// 			filename: Input filename
        //                      idx     : Parameter index value
	// Output		none	: Electron G interaction vector filled with
	//				  parameters read from file

 
MSVCDLL bool read(const ParameterSet& pset, int idx=-1, int warn=2);
 
        // Input                IGV     : Electron G interaction vector (this)
        //                      pset    : Parameter set
        //                      idx     : Parameter index value used for
        //                                prefix [#] in input names
        // Output               none    : Electron G interaction vector filled
        //                                with parameters read from pset

// ----------------------------------------------------------------------------
//       Interactive Read of Electron G Interaction Vector From An ASCII File
// ----------------------------------------------------------------------------

MSVCDLL std::string ask_read(int argc, char* argv[], int argn);

	// Input		IGV	: Electron G interaction vector (this)
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

	// Input		IGV	: Electron G interaction vector (this)
	// 			ostr	: Output stream
	//			full    : Flag for long vs short output
	// Output		non	: Electron G interaction vector
	//			          parameters sent to the output stream


MSVCDLL friend std::ostream& operator<< (std::ostream& out, const IntGVec& IGV);

	// Input		IGV	: Electron G interaction vector
	// 			out	: Output stream
	// Output		none	: Electron G interaction vector
	//				  parameters sent to output stream
};

#endif							// IntGVec.h
