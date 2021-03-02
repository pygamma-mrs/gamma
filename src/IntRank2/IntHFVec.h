/* IntHFVec.h ***************************************************-*-c++-*-
**									**
** 	                         G A M M A				**
**			         				 	**
**	Hyperfine Interactions Vector 			Interface	**
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
** This class maintains a vector of rank 2 electron-nucleon hyperfine	**
** interactions (as defined in class IntHF). The class allows users, or **
** more importantly spin systems, to manipulate an array of such	**
** interactions either in concert or individually.			**
**						 			**
*************************************************************************/

///Chapter Class Rank 2 Hyperfine Interactions Vector (IntHFVec)
///Section Overview
///Body    The class 
///Section Available Hyperfine Interactions Vector Functions

#ifndef   IntHFVec_h_			// Is this file already included?
#  define IntHFVec_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <IntRank2/IntHF.h>		// Include hyperfine interactions
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Basics/Isotope.h>		// Include spin isotopes
#include <string>			// Include libstdc++ strings
#include <vector>			// Include libstdc++ STL vectors

//typedef std::vector<IntHF> stdvecHF;

//class IntHFVec : public stdvecHF
class IntHFVec : public std::vector<IntHF>
  {

private:

// ----------------------------------------------------------------------------  
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                 HYPERFINE INTERACTION VECTOR ERROR HANDLING
// ____________________________________________________________________________
 

/*      Input			IHFV	: Hyperfine interaction vector (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
                                pname   : Added error message for output
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

void IHFVerror(int eidx,                           int noret=0) const;
void IHFVerror(int eidx, const std::string& pname, int noret=0) const;
volatile void IHFVfatality(int eidx) const;

// ____________________________________________________________________________
// ii             HYPERFINE INTERACTION VECTOR SETUP FUNCTIONS
// ____________________________________________________________________________
 

bool check_spin(int spin, int warn=0) const;
 
        // Input                IHFV     : Hyperfine interaction vector (this)
        //                      spin	: A spin index
        //                      warn    : Warning level 
        // Output               TF      : Returns TRUE if spin is valid
        //                                index, FALSE if not

// ----------------------------------------------------------------------------
//     Functions To Set Hyperfine Interactions From Single Index Parameters
// ----------------------------------------------------------------------------

/* Hyperfine interactions may be specified using a single interaction index
   (as opposed to spin pair indices). The functions in this section attempt to
   read a vector of N such interactions that are index sequentially from [0,N).
   The value of N, if not set, weill be determined by the first interaction
   that i s not fully specified in the parameter set.

   The function NInts tries to determine the maximum vector length, i.e. the
   number of interactions that should be created. This is done by looking for

        // Input                IHFV    : Hyperfine interaction vector (this)
        //                      pset    : A parameter set
        //                      idx     : Parameter prefix index, [#]
                                warn    : Warning level
                                            0 - no warnings
                                            1 - warn if no spins specified
                                            2 - fatal if no spins specified
        // Output               void    : Interaction vector filled from
        //                                pset single indexed parameters    */

void setHFs(const ParameterSet& pset, int idx, int warn=2);
//int getNInts(const ParameterSet& pset, int idx=-1) const;
 
        // Input                IHFV     : Hyperfine interaction vector (this)
        //                      pset    : A parameter set
        // Output               ns      : The number of Hyperfine interactions
        //                                referenced by single insices in pset
        // Note                         : This function looks for G(#) where
        //                                # = { 0, 1, ...., ns-1 }.


//int setGs(const ParameterSet& pset, int idx=-1, int warn=0);

        // Input                IHFV     : Hyperfine interaction vector (this)
        //                      pset    : A parameter set
        //                      idx     : Parameter prefix index
        //                      warn    : Warning flag
        //                                    0 = no warnings
        //                                   !0 = warnings
        // Output               none    : System Hyperfine interactions are
        //                                set from parameters in pset
        // Note                         : Assumed that the number of spins
        //                                is ALREADY set and space allocated
  
/* This function uses the electron G interaction constructors directly.
   In turn, they use a single interaction (or spin) index.  They'll attempt
   to read, for each interaction in the vector, a set of parameters 
 
                { [Iso(i)/CI], Hyperfine, Gtheta, Gphi, Geta } 
 
   where the shielding anisotropy is both read in and maintained in PPM units.
   Hyperfine is related to the spatial tensor delzz via
 
      ^          3                                     1  
     / \ sigma = - del   = sigma  - sigma   = sigma  - - [ sigma  + sigma  ] 
     ---         2    zz        ||        |        zz  2        xx       yy
                                         --- 
   An Example Shift Anisotropy Interaction Definition Are The ASCII Lines 
 
   CI(1)        (1) : 1.5                  - Spin I quantum number         
   G(1)         (1) : 40.0                 - Shift anisotropy (PPM) 
   Gphi(1)      (1) : 45.0                 - Tensor orientation (deg) 
   Gtheta(1)    (1) : 45.0                 - Tensor orientation (deg) 
   Geta(1)      (1) : 0.0                  - Tensor asymmetry [0,1]          */ 



// ----------------------------------------------------------------------------
//   Functions To Set Hyperfine Interactions From Two Spin Indexed Parameters
// ----------------------------------------------------------------------------

/* Hyperfine interactions can be specified using two spin indices, i.e. spin
   pair indices (as opposed to using a single interaction index). These
   functions attempt to read all such interactions spanning N spins.  The
   vector vector will contain a hyperfine interaction for each unique spin
   pair, although it may be NULL if the interaction is not present in the
   parameter set.

           Input                IHFV    : Hyperfine interaction vector (this)
                                pset    : A parameter set
                                indx    : Parameter prefix index ([#])
                                warn    : Warning level
                                            0 - no warnings
                                            1 - warn if no spins specified
                                            2 - fatal if no spins specified
           Output               ns      : # of spins specified in pset       */

int  getNSpins(const ParameterSet& pset, int indx, int warn=2) const;
void setHFs(int ns, const ParameterSet& pset, int idx, int warn=2);

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

// ____________________________________________________________________________ 
// A           HYPERFINE INTERACTION VECTOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________ 

// ---------------------------------------------------------------------------- 
//                  Simple Constructors That Don't Do Much
// ---------------------------------------------------------------------------- 
 

/* These functions are inherited from the ANSI standard "vector" class in the
   C++ library libstdc++.  I am listing them here so that I & other users don't
   have to keep looking them up all the time.
*/

MSVCDLC IntHFVec();				// Empty Interaction Vector
MSVCDLC IntHFVec(const IntHFVec& HFV);	// Copy Of Interaction Vector

// ----------------------------------------------------------------------------
//                      Constructors Using Parameter Sets 
// ----------------------------------------------------------------------------
   
        // Input                IHFV     : Hyperfine interaction vector (this)
        //                      ns      : Number of spins in vector
        //                      pset    : Parameter set
        //                      idx     : Index for Hyperfine vector
        //                                i.e [idx] is parameter name prefix
        //                      warn    : Warning level
        // Note                         : If idx is negative then no
        //                                parameter prefix is used
	//				  else prefix is "[idx]"
	// Note				: If ns is not specified as an
	//				  argument the function will attempt
	//				  to read {[idx]Hyperfine(i)} from i=0 
	//				  to i=ns-1 where ns is the number

//IntHFVec(const         ParameterSet& pset, int indx=-1, int warn=1);
//IntHFVec(int ns, const ParameterSet& pset, int indx=-1, int warn=1);

// ----------------------------------------------------------------------------
//          This Constructor Supports Generation From A Spin System
// ----------------------------------------------------------------------------

/* In this case the vector of interactions is associated with a list of isotope
   types (Isos). The interaction vector length will match the vector of spin
   types. For each electron spin type a G interaction will be generated from
   parameters in the parameter set.

           Input                IHFV     : Hyperfine interaction vector (this)
                                Isos    : Array of isotope/spin types
                                pset    : Parameter set
                                warn    : Warning level
                                                0 - no warnings
                                                1 - warn if spin is nucleon
                                                2 - fatal if spin is nucleon
           Output               none    : Hyperfine interaction vector
                                          constructed with spin labels
                                          in Isos and parameters in pset
           Note                         : We don't allow users to associate
                                          a nucleon with a G interaction     */



MSVCDLC IntHFVec(const std::vector<Isotope>& Isos,
                                         const ParameterSet& pset, int warn=2);

// ----------------------------------------------------------------------------
//                    Constructors Using Parameter Sets
// ----------------------------------------------------------------------------

/* These functions allow users to fill the hyperfine interactions vector from
   a GAMMA parameter set (external ASCII file). Variants enable higher classes
   such as spin systems to generate the vector when they have a knowledge of
   spin isotope types and hyperfine parameters (parameter set)

        // Input                IHFV	: Hyperfine interaction vector (this)
        //                      pset    : Parameter set
        //                      Isos    : Array of isotope types
        //                      idx     : Index for dipolar vector
        //                                i.e [idx] is parameter name prefix
        //                      warn    : Warning level (applicable to use
        //                                of isotopes & coordinates)
        // Note                         : If idx is negative then no
        //                                parameter prefix will be used     */

//IntHFVec(const ParameterSet& pset,          int indx=-1, int warn=1);
MSVCDLC IntHFVec(const ParameterSet& pset,
                    const std::vector<Isotope>& Isos, int indx=-1, int warn=2);

// ----------------------------------------------------------------------------
//               Here Be The Assignment Operator & Destructor
// ----------------------------------------------------------------------------

MSVCDLL void operator= (const IntHFVec& IHFV);
MSVCDLC      ~IntHFVec ();

// ____________________________________________________________________________ 
// B                Hyperfine INTERACTION ACCESS & MANIPULATIONS
// ____________________________________________________________________________ 

// ---------------------------------------------------------------------------- 
//                   Generic Hyperfine Value Access Functions
// ---------------------------------------------------------------------------- 

	// Input		IHFV	: Hyperfine interaction vector
	// 			spin	: Hyperfine index
	// 			val	: Hyperfine interaction value
	//				  <=0: Hyperfine coupling of spin (Hertz)
	//				    1: Asymmetry [0, 1]
	//				    2: Theta (degrees, down from +z)
	//				    3: Phi   (degrees, over from +z)
	// Output		none	: Get/Set Hyperfine interaction value
	//				  for specified spin
 

MSVCDLL void   CValue(int spin, double val, int type);
MSVCDLL double CValue(int spin,             int type) const;


// ---------------------------------------------------------------------------- 
//                         Hyperfine Coupling Constants
// ---------------------------------------------------------------------------- 

	// Input		IHFV	: Hyperfine interaction vector
	// 			spin	: Hyperfine index
	// 			hfa	: Hyperfine anisotropy (G)
	// Output		none	: Get/Set Hyperfine anisotropy
	// Note				: Defined in class IntHF as 1.5 times
	//			          the Hyperfine tensor delzz value

//void   hfa(int  spin, double Hyperfine);
//double hfa(int  spin) const;
MSVCDLL void   delz(int spin, double dcc);
MSVCDLL double delz(int spin) const;

// ---------------------------------------------------------------------------- 
//                        Hyperfine Asymmetry Values
// ---------------------------------------------------------------------------- 

	// Input		IHFV	: Hyperfine interaction vector
	// 			spin	: Hyperfine index
	// 			ceta	: Hyperfine interaction asymmetry
	// Output		none	: Get/Set Hyperfine spin asymmetry
	// Note				: Defined in class IntHyperfine between [0,1]

MSVCDLL void   eta(int spin, double ceta);
MSVCDLL double eta(int spin) const;

 
// ---------------------------------------------------------------------------- 
//               Hyperfine Theta Orientation (Down From +z Axis)
// ---------------------------------------------------------------------------- 

	// Input		IHFV	: Hyperfine interaction vector
	// 			spins	: Hyperfine index
	// 			ctheta	: Hyperfine interaction angle (deg)
	// Output		none	: Get/Set Hyperfine spin theta angle
	// Note				: Defined class IntHyperfine between [0,180]

MSVCDLL void   theta(int spin, double ctheta);
MSVCDLL double theta(int spin) const;

// ---------------------------------------------------------------------------- 
//               Hyperfine Phi Orientation (Over From +x Axis)
// ---------------------------------------------------------------------------- 

	// Input		IHFV	: Hyperfine interaction vector
	// 			spins	: Hyperfine index
	// 			cphi	: Hyperfine interaction angle (deg)
	// Output		none	: Get/Set Hyperfine spin phi angle
	// Note				: Defined in IntHyperfine between [0,360]

MSVCDLL void   phi(int spin, double phi);
MSVCDLL double phi(int spin) const;

// ---------------------------------------------------------------------------- 
//                        Full Hyperfine Interaction
// ---------------------------------------------------------------------------- 

MSVCDLL IntHF& operator() (int spins);                                          

        // Input                IHFV	: Hyperfine interaction vector (this)
        //                      spins	: Interaction index
        // Ouput                CI	: The i'th Hyperfine interaction in IHFV
        // Note				: Returns a reference to the interaction


MSVCDLL IntHF get(int spins) const;

	// Input		IHFV	: Hyperfine interaction vector
	// 			spins	: Hyperfine index
	// Output		CI	: Return rank 2 Hyperfine interaction


// ----------------------------------------------------------------------------
//                 Other Hyperfine Interaction Vector Info
// ----------------------------------------------------------------------------
 
 
//int size() const;
 
        // Input                IHFV   : Hyperfine interaction vector
        // Output               ns    : Number of interactions in vector
 
 
MSVCDLL bool nonzero() const;
 
        // Input                IHFV   : Hyperfine interaction vector
        // Output               TF    : True if any interactions with a
        //                              finite Hyperfine value


// ____________________________________________________________________________
// C                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//     Functions To Make A Parameter Set From A Hyperfine Interaction Vector
// ----------------------------------------------------------------------------

// Note that the preferred means of specifying Hyperfine interactions when there
// is no knowledge of spin isotope types nor spin coordinates are used for
// filling the parameter set.  For base parameters of insividual interactions
// see class IntHyperfine. Additional parameters for the interaction vector will 
// be defined in the IntRank2.cc file in this section.

MSVCDLL operator ParameterSet( ) const;

	// Input		IHFV	: Hyperfine interaction vector
        //                      pset	: Parameter set
        // Output               pset	: Parameter set with only
        //                                Hyperfine interaction parameters
 

MSVCDLL friend void operator+= (ParameterSet& pset, const IntHFVec &IHFV);

	// Input		IHFV	: Hyperfine interaction vector
	//  			pset	: Parameter set
	// Output		pset	: Parameter set with all vector
	//			          Hyperfine interaction parameters added
 
         
MSVCDLL void PSetAdd(ParameterSet& pset, int idx=-1) const;

        // Input                IHFV     : Hyperfine interaction vector
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
//     Functions To Make A Hyperfine Interaction Vector From A Parameter Set
// ---------------------------------------------------------------------------- 


MSVCDLL void operator= (const ParameterSet& pset);

        // Input                IHFV     : Hyperfine interaction vector (this)
        //                      pset    : A parameter set
        // Output               none    : Hyperfine interaction vector filled with
        //                                parameters in pset
        // Note                         : We just use a member function 
        //                                which is needed to allow for 
        //                                [#] prefixed parameter names 



MSVCDLL bool setIHFVec(const ParameterSet& pset, int idx=-1, int warn=2);

	// Input		IHFV	: Hyperfine interaction vector (this)
	// 			pset	: A parameter set
        //                      idx	: Parameter index value used for
        //                                prefix [#] in input names
        //                      warn	: Warning output level
        //                                      0 = no warnings
        //                                      1 = warnings
        //                                     >1 = fatal warnings
	// Output		none	: The entire system is set
	//			  	  from parameters in pset
	// Note				: This uses the assignment from pset
	//				  It exists to support prefix insices


// ----------------------------------------------------------------------------
//      Output Hyperfine Interaction Vector To ASCII From A Parameter Set
// ----------------------------------------------------------------------------
 
        // Input                IHFV     : Hyperfine interaction vector (this)
        //               THIS   filename: Output file name
        //            OR THIS   ofstr   : Output file stream
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning level
        // Output               none    : Hyperfine interaction vector written
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

	// Input		IHFV	: Hyperfine interaction vector (this)
	// 			filename: Input filename
        //                      idx     : Parameter index value
	// Output		none	: Hyperfine interaction vector filled with
	//				  parameters read from file

 
MSVCDLL bool read(const ParameterSet& pset, int idx=-1, int warn=2);
 
        // Input                IHFV     : Hyperfine interaction vector (this)
        //                      pset    : Parameter set
        //                      idx     : Parameter index value used for
        //                                prefix [#] in input names
        // Output               none    : Hyperfine interaction vector filled
        //                                with parameters read from pset

// ----------------------------------------------------------------------------
//       Interactive Read of Hyperfine Interaction Vector From An ASCII File
// ----------------------------------------------------------------------------

MSVCDLL std::string ask_read(int argc, char* argv[], int argn);

	// Input		IHFV	: Hyperfine interaction vector (this)
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

	// Input		IHFV	: Hyperfine interaction vector (this)
	// 			ostr	: Output stream
	//			full    : Flag for long vs short output
	// Output		non	: Hyperfine interaction vector
	//			          parameters sent to the output stream


MSVCDLL friend std::ostream& operator<< (std::ostream& out, const IntHFVec& IHFV);

	// Input		IHFV	: Hyperfine interaction vector
	// 			out	: Output stream
	// Output		none	: Hyperfine interaction vector
	//				  parameters sent to output stream
};

#endif							// IntHFVec.h
