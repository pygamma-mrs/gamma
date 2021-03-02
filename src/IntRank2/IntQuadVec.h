/* IntQuadVec.h *************************************************-*-c++-*-
**									**
** 	                         G A M M A				**
**			         				 	**
**	Quadrupolar Interactions Vector 		Interface	**
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
** This class maintains a vector of rank 2 Quadrupolar interactions (as	**
** defined in class IntQuad). The class lets users, or more importantly	**
** spin systems, to manipulate an array of such	interactions either in	**
** concert or individually.						**
**						 			**
*************************************************************************/

///Chapter Class Rank 2 Quadrupolar Interactions Vector (IntQuadVec)
///Section Overview
///Body    The class 
///Section Available Quadrupolar Interactions Vector Functions

#ifndef   IntQuadVec_h_			// Is this file already included?
#  define IntQuadVec_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#if defined(_MSC_VER)			// If we are using MSVC++
 #pragma warning (disable : 4786)       // Kill STL namelength warnings
#endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <IntRank2/IntQuad.h>		// Include Quadrupolar interactions
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Basics/Isotope.h>		// Include spin isotopes
#include <string>			// Include libstdc++ strings
#include <vector>			// Include libstdc++ STL vectors

//typedef std::vector<IntQuad> stdvecG;	// These would simplify code but
//class IntQuadVec : public stdvecG	// some compilers had trouble...

class IntQuadVec : public std::vector<IntQuad>
  {

private:

// ----------------------------------------------------------------------------  
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                 QUADRUPOLAR INTERACTION VECTOR ERROR HANDLING
// ____________________________________________________________________________
 

/*      Input			IQV	: Quadrupolar interaction vector (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
                                pname   : Added error message for output
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

void IQVerror(int eidx,                           int noret=0) const;
void IQVerror(int eidx, const std::string& pname, int noret=0) const;
volatile void IQVfatality(int eidx) const;

// ____________________________________________________________________________
// ii             QUADRUPOLAR INTERACTION VECTOR SETUP FUNCTIONS
// ____________________________________________________________________________
 

bool check_spin(int spin, int warn=0) const;
 
        // Input                IQV     : Quadrupolar interaction vector (this)
        //                      spin	: A spin index
        //                      warn    : Warning level 
        // Output               TF      : Returns TRUE if spin is valid
        //                                index, FALSE if not


// ----------------------------------------------------------------------------
//      Functions To Set Quadrupolar Interactions From Single Index Parameters
// ----------------------------------------------------------------------------

int getNInts(const ParameterSet& pset, int idx=-1) const;
 
        // Input                IQV     : Quadrupolar interaction vector (this)
        //                      pset    : A parameter set
        // Output               ns      : The number of Quadrupolar interactions
        //                                referenced by single insices in pset
        // Note                         : This function looks for G(#) where
        //                                # = { 0, 1, ...., ns-1 }.


int setGs(const ParameterSet& pset, int idx=-1, int warn=0);

        // Input                IQV     : Quadrupolar interaction vector (this)
        //                      pset    : A parameter set
        //                      idx     : Parameter prefix index
        //                      warn    : Warning flag
        //                                    0 = no warnings
        //                                   !0 = warnings
        // Output               none    : System Quadrupolar interactions are
        //                                set from parameters in pset
        // Note                         : Assumed that the number of spins
        //                                is ALREADY set and space allocated
  
/* This function uses the Quadrupolar interaction constructors directly.
   In turn, they use a single interaction (or spin) index.  They'll attempt
   to read, for each interaction in the vector, a set of parameters 
 
                { [Iso(i)/CI], Quadrupolar, Gtheta, Gphi, Geta } 
 
   where the shielding anisotropy is both read in and maintained in PPM units.
   Quadrupolar is related to the spatial tensor delzz via
 
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
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

// ____________________________________________________________________________ 
// A           QUADRUPOLAR INTERACTION VECTOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________ 

// ---------------------------------------------------------------------------- 
//                  Simple Constructors That Don't Do Much
// ---------------------------------------------------------------------------- 

/* These functions are inherited from the ANSI standard "vector" class in
   the C++ library libstdc++.  I am listing them here so that I & other 
   users don't have to keep looking them up all the time.                    */

MSVCDLC    IntQuadVec();				// Empty Interaction Vector
// IntQuadVec(int N)				// Vector w/ N Interactions
// IntQuadVec(int N, const IntQuad& GI)		// Vector w/ N pars
// IntQuadVec(const IntQuadVec& QuadVec)	// Vector copy of QuadVec
// IntQuadVec assign(N)				// Assign N element
// ~IntQuadVec()				// Destructor of Vector

// ----------------------------------------------------------------------------
//                      Constructors Using Parameter Sets 
// ----------------------------------------------------------------------------
   
        // Input                IQV     : Quadrupolar interaction vector (this)
        //                      ns      : Number of spins in vector
        //                      pset    : Parameter set
        //                      idx     : Index for Quadrupolar vector
        //                                i.e [idx] is parameter name prefix
        //                      warn    : Warning level
        // Note                         : If idx is negative then no
        //                                parameter prefix is used
	//				  else prefix is "[idx]"
	// Note				: If ns is not specified as an
	//				  argument the function will attempt
	//				  to read {[idx]Quadrupolar(i)} from i=0 
	//				  to i=ns-1 where ns is the number

//IntQuadVec(const         ParameterSet& pset, int indx=-1, int warn=1);
//IntQuadVec(int ns, const ParameterSet& pset, int indx=-1, int warn=1);

// ----------------------------------------------------------------------------
//          This Constructor Supports Generation From A Spin System
// ----------------------------------------------------------------------------

/* In this case the vector of interactions is associated with a list of isotope
   types (Isos). The interaction vector length will match the vector of spin
   types. For each electron spin type a G interaction will be generated from
   parameters in the parameter set.

           Input                IQV     : Quadrupolar interaction vector (this)
                                Isos    : Array of isotope/spin types
                                pset    : Parameter set
                                warn    : Warning level
                                                0 - no warnings
                                                1 - warn if spin is nucleon
                                                2 - fatal if spin is nucleon
           Output               none    : Quadrupolar interaction vector
                                          constructed with spin labels
                                          in Isos and parameters in pset
           Note                         : We don't allow users to associate
                                          a nucleon with a G interaction     */



MSVCDLC IntQuadVec(const std::vector<Isotope>& Isos,
                                         const ParameterSet& pset, int warn=2);


// ----------------------------------------------------------------------------
//               Here Be The Assignment Operator & Destructor
// ----------------------------------------------------------------------------

//  These should be completely handled by the base class STL vector<IntQuad>

//void operator= (const IntQuadVec &IQV);
//     ~IntQuadVec ();


// ____________________________________________________________________________ 
// B                Quadrupolar INTERACTION ACCESS & MANIPULATIONS
// ____________________________________________________________________________ 

// ---------------------------------------------------------------------------- 
//                   Generic Quadrupolar Value Access Functions
// ---------------------------------------------------------------------------- 

	// Input		IQV	: Quadrupolar interaction vector
	// 			spin	: Quadrupolar index
	// 			val	: Quadrupolar interaction value
	//				  <=0: Quadrupolar coupling of spin (Hertz)
	//				    1: Asymmetry [0, 1]
	//				    2: Theta (degrees, down from +z)
	//				    3: Phi   (degrees, over from +z)
	// Output		none	: Get/Set Quadrupolar interaction value
	//				  for specified spin

MSVCDLL void   QValue(int spin, double val, int type);
MSVCDLL double QValue(int spin, int type) const;

// ---------------------------------------------------------------------------- 
//                         Quadrupolar Coupling Constants
// ---------------------------------------------------------------------------- 

	// Input		IQV	: Quadrupolar interaction vector
	// 			spin	: Quadrupolar index
	// 			qcc	: Quadrupolar coupling constant (Hz)
	// Output		none	: Get/Set spin quadrupolar coupling
	// Note				: QCC is identical to delzz for this
	//				  interaction and the anisotropy is
	//				  1.5 times the Quadrupolar coupling

MSVCDLL void   QCC(int  spin, double qcc);
MSVCDLL void   NQCC(int spin, double qcc);
MSVCDLL void   delz(int spin, double qcc);
MSVCDLL void   QA(int   spin, double qa);

MSVCDLL double QCC(int  spin) const;
MSVCDLL double NQCC(int spin) const;
MSVCDLL double delz(int spin) const;
MSVCDLL double QA(int   spin) const;

// ---------------------------------------------------------------------------- 
//                        Quadrupolar Asymmetry Values
// ---------------------------------------------------------------------------- 

	// Input		IQV	: Quadrupolar interaction vector
	// 			spin	: Quadrupolar index
        //                      qeta    : Quadrupolar interaction asymmetry
        // Output               none    : Get/Set Quadrupolar asymmetry
        // Note                         : Defined to be between [0,1]

MSVCDLL void   eta(int spin, double ceta);
MSVCDLL double eta(int spin) const;
 

//-----------------------------------------------------------------------------
//                         Orientation Angle Access
//-----------------------------------------------------------------------------

/* There are 3 angles which orient the interaction spatial tensor relative to 
   the interactions own principal axis system (PAS).  These are the set of 
   Euler angles { alpha, beta, gamma }. They corrspond to rotating the 
   coordinate axes first about the z-axis by alpha followed by a rotation about
   the new x-axis by beta and lastly about the new z-axis by gamma. For a 
   symmetric tensor the last rotation is of no consequnce and set to zero, and 
   for powder averages on the first two angles are used to sum over all 
   spatial orientations. In these two cases the angles alpha and beta are one
   and the same as the spherical coordinate angles phi and theta respectively. 
   Theta is the angle down from the PAS z-axis & the angle phi which is over 
   from the PAS x-axis. The ranges of the angles (internally imposed on 
   GAMMA's Euler angles) are alpha,gamma,phi = [0,360] and beta,theta=[0,180].
   These functions allow users to both obtain and set these angles for any
   intraction. Setting any of the angles will effectively reorient the 
   spatial tensor (interaction).                                             */

MSVCDLL double  alpha(int spin)       const;
MSVCDLL double  beta(int spin)        const;
MSVCDLL double  gamma(int spin)       const;
MSVCDLL double  phi(int spin)         const;
MSVCDLL double  theta(int spin)       const;
//EAngles orientation(int spin) const;

MSVCDLL void alpha(int spin,double A);
MSVCDLL void beta(int  spin,double B);
MSVCDLL void gamma(int spin,double G);
MSVCDLL void phi(int   spin,double P);
MSVCDLL void theta(int spin,double T);

//void orientation(int spin, const EAngles& EA);
//void orientation(int spin,
//                            double A, double B, double G, bool deg=false); 

// ---------------------------------------------------------------------------- 
//                        Full Quadrupolar Interaction
// ---------------------------------------------------------------------------- 
        // Input                IQV	: Quadrupolar interaction vector (this)
        //                      spins	: Interaction index
        // Ouput                CI	: The i'th Quadrupolar interaction in IQV
        // Note				: Returns a reference to the interaction

MSVCDLL       IntQuad& operator() (int i);
MSVCDLL const IntQuad& getcref(int     i) const;
MSVCDLL       IntQuad  get(int         i) const;

// ----------------------------------------------------------------------------
//                 Other Quadrupolar Interaction Vector Info
// ----------------------------------------------------------------------------
 
 
//int size() const;				// INHERITED
MSVCDLL int nonzero() const;

// ____________________________________________________________________________
// C                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//     Functions To Make A Parameter Set From A Quadrupolar Interaction Vector
// ----------------------------------------------------------------------------

// Note that the preferred means of specifying Quadrupolar interactions when
// there is no knowledge of spin isotope types will be used for filling the
// parameter set. For base parameters of insividual interactions see class 
// IntQuad. We simply add an index & potential prefix to all parameters here.

	// Input		IQV	: Quadrupolar interaction vector
        //                      pset	: Parameter set
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        // Output               pset	: Parameter set with (only/added)
        //                                Quad. vector interaction parameters
        // Note                         : The parameter names & types output
        //                                here MUST match those used in setting
        //                                the vector up from parameters sets
         
MSVCDLL operator ParameterSet( ) const;
MSVCDLL friend void operator+= (ParameterSet& pset, const IntQuadVec &IQV);
MSVCDLL void PSetAdd(ParameterSet& pset, int idx=-1) const;


// ---------------------------------------------------------------------------- 
//     Functions To Make A Quadrupolar Interaction Vector From A Parameter Set
// ---------------------------------------------------------------------------- 

/* These functions attempt to set an entire quadrupolar interactions vector 
   from paraemters found in a parameter set. An optional index may be specified
   for a prefix that will be used on all interaction parameters.


        // Input                IQV     : Quadrupolar interaction vector (this)
        //                      pset    : A parameter set
        //                      idx	: Parameter index value used for
        //                                prefix [#] in input names
        //                      warn	: Warning output level
        //                                      0 = no warnings
        //                                      1 = warnings
        //                                     >1 = fatal warnings
        // Output               none    : Quadrupolar interaction vector filled
        //                                with parameters in pset
        // Note                         : We just use a member function 
        //                                which is needed to allow for 
        //                                [#] prefixed parameter names 
	// Note				: This uses the assignment from pset
	//				  & exists to support prefix indices */

MSVCDLL void operator= (const ParameterSet& pset);
MSVCDLL bool setIQVec(const ParameterSet& pset, int idx=-1, int warn=2);

// ----------------------------------------------------------------------------
//      Output Quadrupolar Interaction Vector To ASCII From A Parameter Set
// ----------------------------------------------------------------------------
 
        // Input                IQV     : Quadrupolar interaction vector (this)
        //               THIS   fname   : Output file name
        //            OR THIS   ofstr   : Output file stream
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning level
        // Output               none    : Quadrupolar interaction vector written
        //                                as a parameter set either to file
	//				  fname or ourput stream ofstr
 
MSVCDLL int write(const std::string &fname, int idx=-1, int warn=2) const;
MSVCDLL int write(std::ofstream& ofstr,     int idx=-1, int warn=2) const;

// ____________________________________________________________________________ 
// D                             INPUT FUNCTIONS
// ____________________________________________________________________________ 

// ----------------------------------------------------------------------------
//        Direct Read of Vector From An ASCII File Or A Parameter Set
// ----------------------------------------------------------------------------


	// Input		IQV	: Quadrupolar interaction vector (this)
	// 			filename: Input filename
        //                      pset    : Parameter set
        //                      idx     : Parameter index value used for
        //                                prefix [#] in input names
	// Output		none	: Quadrupolar interaction vector filled
	//				  with values read from file or pset

MSVCDLL bool read(const std::string &filename, int idx=-1, int warn=2);
MSVCDLL bool read(const ParameterSet& pset,    int idx=-1, int warn=2);

// ----------------------------------------------------------------------------
//       Interactive Read of Quadrupolar Interaction Vector From An ASCII File
// ----------------------------------------------------------------------------

MSVCDLL std::string ask_read(int argc, char* argv[], int argn);

	// Input		IQV	: Quadrupolar interaction vector (this)
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

/* Functions that output the quadrupolar interaction vector to standard output
   or any output stream. These are NOT parameter formatted but formatted just
   to be presentable.

	   Input		IQV	: Quadrupolar interaction vector (this)
	   			ostr	: Output stream
	  			full    : Flag for long vs short output
	   Output		ostr    : Quadrupolar interaction vector
	  			          parameters sent to output stream   */

MSVCDLL        std::ostream& print(std::ostream& ostr, int full=0) const;
MSVCDLL friend std::ostream& operator<< (std::ostream& ostr, const IntQuadVec& IQV);

};

#endif							// IntQuadVec.h
