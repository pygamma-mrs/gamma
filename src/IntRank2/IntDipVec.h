/* IntDipVec.h **************************************************-*-c++-*-
**									**
** 	                         G A M M A				**
**			         				 	**
**	Dipole Interactions Vector 			Interface	**
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
** This class maintains a vector of rank2 dipolar interactions (from    **
** class IntRank2).  The class allows users to manipulate an array of   **
** such interactions simultaneously.                                    **
**						 			**
*************************************************************************/

///Chapter Class Rank 2 Dipolar Interactions Vector (IntDipVec)
///Section Overview
///Body    The class 
///Section Available Dipolar Interactions Vector Functions

#ifndef   IntDipVec_h_			// Is this file already included?
#  define IntDipVec_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

class ParameterSet;			// Know parameter sets
class IntDip;				// Know dipolar interactions
class coord_vec;			// Know about coordinate vectors

#if defined(_MSC_VER)			// If using MSVC++ then we
 #pragma warning (disable : 4786)       // Kill STL namelength warnings
#endif 

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <IntRank2/IntDip.h>		// Include dipolar interactions
#include <Basics/ParamSet.h>		// Inlcude parameter sets
#include <string>			// Include stdlibc++ strings

class IntDipVec : public std::vector<IntDip>
  {

private:

// ----------------------------------------------------------------------------  
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                 CLASS DIPOLE INTERACTION VECTOR ERROR HANDLING
// ____________________________________________________________________________

	// Input		IDV	: Dipole interaction vector (this)
        // 			eidx	: Error index
        //                      noret	: Flag for return (0=linefeed)
        // Output               none	: Error message

	// Input		IDV	: Dipole interaction vector (this)
        // 			eidx	: Error index
        //                      pname	: Additional string for error
        //                      noret	: Flag for return (0=linefeed)
        // Output               none	: Error message

     
void IDVerror(int eidx,                           int noret=0) const;
void IDVerror(int eidx, const std::string& pname, int noret=0) const;
volatile void IDVfatal(int eidx) const;
volatile void IDVfatal(int eidx, const std::string& pn) const;

// ____________________________________________________________________________
// ii               DIPOLE INTERACTION VECTOR SETUP FUNCTIONS
// ____________________________________________________________________________

/* These functions determine specific aspects of a dipolar interactions vector
   from a specified parameter set. As the various interaction parameters tend
   to interact, the functions MUST be private because their misuse could
   produce an inconsistent interactions vector.

   The goal here is quite simple. We need to determine a series of dipolar
   interactions from parameters found in an input GAMMA parameter set. Each
   interaction will be defined by unique parameters in the set that specify
   the values: { Iqn, Sqn, DCC ,eta, alpha, beta, gamma } Complexity arises
   because we need find a series of such interactions and because dipolar
   interactions allow for a variety of parameters for their definition.

   Consider the kth dipolar interaction in the vector. It may have been
   defined by parameters indexed with a single interaction number or it may
   have been defined by a pair of spin indices. If the former is true it is
   clear that the kth interaction is defined by the parameters using index
   k. If the latter is used what spin indices correspond to the kth dipolar
   interaction? That is a problem unless we know the total number of spins
   that are specified, and that is what we will demand. In order to use
   spin indexed parameters we MUST know the total number of spins defined.   */

// ----------------------------------------------------------------------------
//     Functions To Set Number Of Spins (For Using Spin Indexed Parameters)
// ----------------------------------------------------------------------------

/* As mentioned in the previous paragraphs, without knowledge of the total
   number of spins we cannot translate two spin indices into a the single
   interaction index the sets the position of the interaction in the vector.
   We allow several means of setting the total number of spins. These are as
   follows:

        1.) NSpins      The parameter NSpins will be assumed that number
        2.) Iso         A series of parameters Iso(i), i=0,1,2,.... sets it
	3.) Dqn		A series of parameters Dqn(i), i=0,1,2,.... sets it
        4.) DCC         A series of parameters DCC(i,j), i=0,1,2,.. sets it
        5.) Dxx         A series of parameters Dxx(i,j), i=0,1,2,.. sets it

   Note this does NOT allow for mixing up dipolar interaction designations.
   That is, one may not use Iso(i) & Iso(j) for one spin pair then switch to
   Dqn(i) & Dqn(j) for another spin pair. Stay consistent when using spin
   indices for dipolar designations.                                        */


bool getNS(const ParameterSet& pset, int indx,
                                                int& ns, bool warn=true) const;

bool getNSpins(const ParameterSet& pset, int indx,
                                                int& ns, bool warn=true) const;

bool getNIsos(const  ParameterSet& pset, int indx,
                                                int& ni, bool warn=true) const;

bool getNqns(const ParameterSet& pset,int indx,
                                                int& nq, bool warn=true) const;

bool getNdccs(const ParameterSet& pset,int indx,
                                                int& nd, bool warn=true) const;

bool getNdxxs(const ParameterSet& pset,int indx,
                                                int& nd, bool warn=true) const;

// ----------------------------------------------------------------------------
//         Functions To Set The Entire Dipolar Interactions Vector
// ----------------------------------------------------------------------------

/* We allow class IntDip to read inidividual interactions from the parameter
   set. This simplifies our task down to dealing with the interaction indexing.
   If a single index is used per interaction this is simple because then the
   interaction index is one in the same as the vector index. If two indices
   (i.e. spin indices) are used per interactions we have a bit more trouble.
   In that case we MUST know the total number of spins we are dealing with in
   order to convert indices i & j into the vector index k. So, first we look
   for the total number of spins set in the parameters.  If found we will allow
   spin indexing and interaction indexing, if not we only allow interaction
   indexing.

           Input                IDV     : Dipole interaction vector (this)
                                pset    : A parameter set
                                idx     : Parameter prefix index
                                warn 	: Flag for warnings
           Output               none    : System dipolar interactions are
                                          set from parameters in pset

   We Use Only Parameters With Prefix [#] (unles idx=-1). As such we off this
   prefix from all parameter names prior to looking for the interactions.
   This is done before any calls to class IntDip which doesn't use prefixes. */

bool setDIV(const ParameterSet& pset,int idx=-1,bool wn=true);

// ____________________________________________________________________________
// iii          DIPOLE INTERACTIONS VECTOR CHECKING FUNCTIONS
// ____________________________________________________________________________
 
/* These functions insure vector integrity and that access is not beyond
   the vector boundaries                                                    

           Input                IDV     : Dipole interaction vector (this)
                                dip     : A dipole index
                                warn    : Warning level 
           Output               TF      : Returns TRUE if dip is a valid
                                          dipole index, FALSE if not         */

bool CheckDI(int dip, int warn=0) const;

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

// ____________________________________________________________________________ 
// A              DIPOLE INTERACTION VECTOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________ 
///Center Basic Functions


// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

/* These functions should be inherited from the ANSI standard "vector" class in
   the C++ library libstdc++? I am listing them here so that I & other users
   don't have to keep looking them up all the time and because the 
   inheritance doesn't work.                                                 */

MSVCDLC IntDipVec(int N=0);
MSVCDLC IntDipVec(const IntDipVec& IDV1);

// ----------------------------------------------------------------------------
//             Constructors Using Isotopes And Spin Coordinates
// ----------------------------------------------------------------------------

/* A pair of spin isotopes and associated spin coordinates defines a set of
   dipolar interactions. Problems will occur if any of the spin pairs mix up
   electrons and nuclei, since such a pairing is disallowed. A flag is added to
   let the user ignore such pairings (leaving an empty dipolar interaction),
   ignore these pairings but just leave the iteractions empty, or fail entirely
   if such a pairing is encountered. Note that this type of constructor is very
   much tailored for use in spin systems containing mutiple dipole-dipole
   interactions                                                              */

        // Input                IDV     : Dipole interaction vector (this)
        //                      Isos    : Array of isotope types
        //                      cvec    : Array of spin coordinates
        //                      warn    : Warning level
        //                                      0 - no warnings
        //                                      1 - warn if e-/nucleus pair
        //                                      2 - fatal if e-/nucleus pair
        // Output               none    : Dipole interaction vector constructed
        //                                with isotope labels in Isos and spin
        //                                coordinates in cvec.
        // Note                         : Vector size set by cvec size,
        //                                Array Isos must contain >= this size
        // Note                         : We disallow dipolar interactions that
        //                                involve electron-nucleus pairs!  They
        //                                will trigger a NULL interaction in
        //                                the vector & issue cause warnings if
        //                                warn is set nonzero

MSVCDLC IntDipVec(const std::vector<Isotope>& Isos, 
                                            const coord_vec& cvec, int warn=2);

// ----------------------------------------------------------------------------
//                    Constructors Using Parameter Sets
// ----------------------------------------------------------------------------

/* These functions allow users to fill up the dipolar interactions vector from
   a GAMMA parameter set (external ASCII file). Variants enable higher classes
   such as spin systems to generate the vector when they have a knowledge of
   spin isotope types (and coordinates in the dipolar interaction case).

        // Input                IDV     : Dipole interaction vector (this)
        //                      pset    : Parameter set
        //                      Isos    : Array of isotope types
        //                      idx     : Index for dipolar vector
        //                                i.e [idx] is parameter name prefix
        //                      warn    : Warning level (applicable to use
        //                                of isotopes & coordinates)
        // Note                         : If idx is negative then no
        //                                parameter prefix will be used     */

MSVCDLC IntDipVec(const ParameterSet& pset,        int indx=-1, int warn=1);
MSVCDLC IntDipVec(const ParameterSet& pset, 
                    const std::vector<Isotope>& Isos, int indx=-1, int warn=2);

// ----------------------------------------------------------------------------
//               Here Be The Assignment Operator & Destructor
// ----------------------------------------------------------------------------

/* Both assignment and destruction are handled by the base vector class.     */

MSVCDLL void operator= (const IntDipVec &IDV);
MSVCDLC      ~IntDipVec ();

// ____________________________________________________________________________ 
// B                DIPOLAR INTERACTION ACCESS & MANIPULATIONS
// ____________________________________________________________________________ 

// ---------------------------------------------------------------------------- 
//                   Generic Dipolar Value Access Functions
// ---------------------------------------------------------------------------- 

	// Input		IDV	: Dipole interaction vector
	// 			dip	: Dipole index
	// 			val	: Dipolar interaction value
	//				  <=0: DCC coupling of spin (Hertz)
	//				    1: Asymmetry [0, 1]
	//				    2: Theta (degrees, down from +z)
	//				    3: Phi   (degrees, over from +z)
	// Output		none	: Get/Set dipolar interaction value
	//				  for spicified dipole dip
 

MSVCDLL void   DValue(int dip, double val, int type);
MSVCDLL double DValue(int dip, int type) const;


// ---------------------------------------------------------------------------- 
//                         Dipolar Coupling Constants
// ---------------------------------------------------------------------------- 

	// Input		IDV	: Dipole interaction vector
	// 			dip	: Dipole index
	// 			dcc	: Dipolar coupling constant (Hertz)
	// Output		none	: Get/Set dipole dip coupling
	// Note				: Defined in class IntDip as equal to
	//			          the dipolar tensor delzz value

MSVCDLL void   DCC(int dip, double dcc);
MSVCDLL double DCC(int dip) const;
MSVCDLL void   Ddelz(int dip, double dcc);
MSVCDLL double Ddelz(int dip) const;

// ---------------------------------------------------------------------------- 
//                        Dipolar Asymmetry Values
// ---------------------------------------------------------------------------- 

	// Input		IDV	: Dipole interaction vector
	// 			dip	: Dipole index
	// 			deta	: Dipolar interaction asymmetry
	// Output		none	: Get/Set dipole dip asymmetry
	// Note				: Defined in class IntDip between [0,1]
	// Note				: Very unusual if nonzero!

MSVCDLL void   Deta(int dip, double deta);
MSVCDLL double Deta(int dip) const;

 
// ---------------------------------------------------------------------------- 
//               Dipolar Theta Orientation (Down From +z Axis)
// ---------------------------------------------------------------------------- 

	// Input		IDV	: Dipole interaction vector
	// 			dip	: Dipole index
	// 			dtheta	: Dipolar interaction angle (deg)
	// Output		none	: Get/Set dipole dip theta angle
	// Note				: Defined class IntDip between [0,180]

MSVCDLL void   Dtheta(int dip, double dtheta);
MSVCDLL double Dtheta(int dip) const;

// ---------------------------------------------------------------------------- 
//               Dipolar Phi Orientation (Over From +x Axis)
// ---------------------------------------------------------------------------- 

	// Input		IDV	: Dipole interaction vector
	// 			dip	: Dipole index
	// 			dphi	: Dipolar interaction angle (deg)
	// Output		none	: Get/Set dipole dip phi angle
	// Note				: Defined in IntDip between [0,360]

MSVCDLL void   Dphi(int dip, double dphi);
MSVCDLL double Dphi(int dip) const;

// ---------------------------------------------------------------------------- 
//                        Full Dipolar Interaction
// ---------------------------------------------------------------------------- 

/*         Input                IDV     : Dipole interaction vector
                                dip     : Dipole index
           Output               DI      : Return rank 2 dipolar interaction
           Note                         : These functions/operators return
                                          either a reference, constant ref,
                                          or copy of the dipolar interaction */

MSVCDLL       IntDip& operator() (int dip);
MSVCDLL const IntDip& getcref(    int dip) const;
MSVCDLL       IntDip  get(        int dip) const;


// ----------------------------------------------------------------------------
//                 Other Dipolar Interaction Vector Info
// ----------------------------------------------------------------------------
 
MSVCDLL double Izval(int dip) const;
MSVCDLL double Szval(int dip) const;
 
        // Input                IDV   : Dipole interaction vector
        // Output               qn    : Quantum number of I or S (0.5, 1.5,..)
 
//int size() const;				INHERITED
MSVCDLL bool nonzero() const;
 
        // Input                IDV   : Dipole interaction vector
        // Output               TF    : True if any interactions with a
        //                              finite DCC value


// ____________________________________________________________________________
// C                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//     Functions To Make A Parameter Set From A Dipolar Interaction Vector
// ----------------------------------------------------------------------------

// Note that the preferred means of specifying dipolar interactions when there
// is no knowledge of spin isotope types nor spin coordinates are used for
// filling the parameter set.  For base parameters of individual interactions
// see class IntDip. Additional parameters for the interaction vector will 
// be defined in the IntRank2.cc file in this section.

	// Input		IDV	: Dipole interaction vector
        //                      pset	: Parameter set
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        // Output               pset	: Parameter set with only
        //                                dipolar interaction parameters
	//				  or with parameters added
        // Note                         : The parameter names & types
        //                                output here MUST match those used
        //                                in setting the vector up
        //                                from parameters sets

MSVCDLL             operator ParameterSet() const;
MSVCDLL friend void operator+= (ParameterSet& pset, const IntDipVec &IDV);
MSVCDLL        void PSetAdd(ParameterSet& pset, int idx=-1) const;

// ---------------------------------------------------------------------------- 
//     Functions To Make A Dipolar Interaction Vector From A Parameter Set
// ---------------------------------------------------------------------------- 

MSVCDLL void operator= (const ParameterSet& pset);

        // Input                IDV     : Dipole interaction vector (this)
        //                      pset    : A parameter set
        // Output               none    : Dipole interaction vector filled with
        //                                parameters in pset
        // Note                         : We just use a member function 
        //                                which is needed to allow for 
        //                                [#] prefixed parameter names 
// ----------------------------------------------------------------------------
// Functions To Output Dipolar Interaction Vector To ASCII From A Parameter Set
// ----------------------------------------------------------------------------
 
        // Input                IDV     : Dipole interaction vector (this)
        //                      filename: Output file name
        //                      ofstr   : Output file stream
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning level
        // Output               none    : Dipole interaction vector is written
	//				  as a parameter set to file filename
        //                                or to output filestream ofstr

MSVCDLL int write(const std::string& filename, int idx=-1, int wrn=2) const;
MSVCDLL int write(std::ofstream&        ofstr, int idx=-1, int wrn=2) const;

// ____________________________________________________________________________ 
// D                             INPUT FUNCTIONS
// ____________________________________________________________________________ 

// ----------------------------------------------------------------------------
//        Direct Read of Vector From An ASCII File Or A Parameter Set
// ----------------------------------------------------------------------------

        // Input                IDV     : Dipole interaction vector (this)
        //                      filename: Input filename
        //                      pset    : Parameter set
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning output label
        //                                 0 = no warnings
        //                                 1 = warnings
        //                                >1 = fatal warnings
        // Output               none    : Dipole interaction vector filled
        //                                with parameters read from file
        //                                or with parameters read from pset
        //                                Return is true if interaction read

MSVCDLL bool read(const std::string &filename, int idx=-1, int warn=2);
MSVCDLL bool read(const ParameterSet& pset,    int idx=-1, int warn=2);

// ----------------------------------------------------------------------------
//       Interactive Read of Dipole Interaction Vector From An ASCII File
// ----------------------------------------------------------------------------

MSVCDLL std::string ask_read(int argc, char* argv[], int argn);

	// Input		IDV    : Dipole interaction vector (this)
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
// E                             OUTPUT FUNCTIONS
// ____________________________________________________________________________ 

	// Input		IDV	: Dipole interaction vector (this)
	// 			ostr	: Output stream
	//			full    : Flag for long vs short output
	// Output		non	: Dipole interaction vector parameters
	//			          sent to the output stream
	// Input		IDV	: Dipole interaction vector
	// 			out	: Output stream
	// Output		none	: Dipole interaction vector
	//				  parameters sent to output stream

MSVCDLL std::ostream& print(std::ostream& ostr, bool full=false) const;
MSVCDLL friend std::ostream& operator<< (std::ostream& out, const IntDipVec& IDV);

};

#endif							// IntDipVec.h
