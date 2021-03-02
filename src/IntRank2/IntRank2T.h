/* IntRank2T.h **************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Rank 2 Interaction Spin Tensors   	   Interface		**
**                                                                      **
**      Scott Smith                                                     **
**      Copyright (c) 1997                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** A variable of type IntRank2T represents the spin tensor associated	**
** with a rank 2 interaction defined for a particular nuclear spin or	**
** nuclear spin pair.  The spin tensor components are stored in 	**
** irreducible spherical format. The following quantities are defined	**
** as making up objects of this class:					**
**                                                                      **
**  Ival    - Hilbert space of 1st spin involved (typically 2 for e-)   **
**  Sval    - Hilbert space of 2nd spin involved			**
**  T{0-2}  _ Spin tensor components T20, T21, and T22                  **
**  Tm{1,2} _ Spin tensor components T2-1 and T2-2                      **
**                                                                      **
** Individual spin tensors are stored in a linked list, handled by the	**
** class IntSTLList.  This prevents recalculation of the 5 arrays that	**
** are associated with such tensors when the spin(s) involved share the	**
** same Iz & Sz values.  However, the arrays are still copied into new	**
** (equivalent) spin tensors, but that is done by reference in the 	**
** GAMMA matrix class.							**
**                                                                      **
*************************************************************************/

#ifndef   IntRank2T_h_			// Is file already included?
#  define IntRank2T_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface 			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/matrix.h>		// Include class member headers
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Basics/Isotope.h>		// Include spin isotopes
#include <string>			// Include libstdc++ strings
#include <fstream>			// Include libstdc++ file streams
#include <vector>			// Include libstdc++ vectors

enum IST_type{ 				// Interaction types
               SA, 			// SA   = shift anisotropy
               DIP,			// DIP  = dipolar
               QUAD,			// QUAD = quadrupolar
               G,			// G    = electron G
               HF,			// HF   = hyperfine
               SPIN, 			// SPIN = spin with itself
               SPF,			// SPF  = spin with field
               SPSP,			// SPSP = spin-spin
               UNK=99 };		// UNK  = unknown

class IntRank2T
  {
  friend class IntRank2;		// Allow this class full access
  friend class IntDip;			// Allow this class full access
  friend class IntCSA;			// Allow this class full access
  friend class IntQuad;			// Allow this class full access
  friend class IntG;			// Allow this class full access
  friend class IntHF;			// Allow this class full access
  int Ival, Sval;			// Spin Hilbert spaces (2*{I,S}+1)
  matrix T0;				// Spherical spin component m =  0
  matrix T1;				// Spherical spin component m =  1
  matrix Tm1;				// Spherical spin component m = -1
  matrix T2;				// Spherical spin component m =  2
  matrix Tm2;				// Spherical spin component m = -2


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i             RANK 2 INTERACTION SPIN TENSOR ERROR HANDLING
// ____________________________________________________________________________

/*       Input                IST     : Rank 2 spin tensor (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pn      : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

         void ISTerror(int eidx,                        int noret=0) const;
volatile void ISTfatal(int eidx)                                     const;
         void ISTerror(int eidx, const std::string& pn, int noret=0) const;
volatile void ISTfatal(int eidx, const std::string& pn)              const;

// ____________________________________________________________________________
// ii               INTERACTION SPIN TENSOR SETUP FUNCTIONS
// ____________________________________________________________________________

/* These are PRIVATE because the value(s) of I (and S) MUST be set prior
   to the use of these functions.  These do NOT add anything to interaction
   lists (e.g. Dip_HF_list) maintained by spin systems

           Input                IST     : Interaction Spin Tensor (this)
           Output               none    : Interaction Spin Tensor spherical
                                          spin components are generated
           Cases                CSA,G   : SHIFT ANISOTROPY
           Cases                G	: ELECTRON G INTERACTION (Ival==2)
                                          (Ival = 2*0.5 + 1, where I = 1/2 e-)
                                Quad    : QUADRUPOLAR INTERACTION (I > 0.5)
                                Dip     : DIPOLAR INTERACTION
				HF      : HYPERFINE INTERACTION
        				  Always has Ival=2 and/or Sval=2
                                          (Ival = 2*0.5 + 1, where I = 1/2 e-)

                                              +
                                           m  |
                                T    = (-1)  T
                                 2,m          2,-m

                     1/2
                  [4]                      1
CSA,G:      T   = |-| * I         T    = - - I           T    = 0
             2,0  [6]    z         2,1     2  +           2,2

                 1/2
              [1]      2                     1                         1  2
Quad:   T   = |-| * [3I - I(I+1)]    T   = - -(I I + I I )      T    = - I
         2,0  [6]      z              2,1    2  + z   z +        2,2   2  +

                  1/2
               [1]                           1                       1
HF,Dip:  T   = |-| * [3I S - I.S]    T   = - -(I S + I S )    T    = - I S
          2,0  [6]      z z           2,1    2  + z   z +      2,2   2  + +  */

virtual void setSPF();		// Set for spin-field
virtual void setSPQ();		// Set for spin with itself
virtual void setSPSP();		// Set for spin-spin

// ____________________________________________________________________________
// iii           RANK 2 SPIN TENSOR FROM PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

/* These functions allow spin tensors to be set up from external ASCII files
   and/or GAMMA parameters sets. Since the parameters which define a rank 2
   spin tensor are often inter-related, we keep these private to insure that
   tensor integrity is maintained. Since objects of type IntRank2T are not
   intended to be used directly, these functions are typically used by classes
   derived from us.

           Input                IST     : Rank 2 spin tensor (this)
                                pset    : A parameter set
                                Pbase   : Parameter base name
                                idx     : Index value
                                warn    : Warning output label
                                           0 = no warnings
                                          !0 = warnings
           Output               IsoI    : Tensor spin isotope type
                                          obtained from parameters in pset
                                I       : Tensor spin quantum number
                                          obtained from parameters in pset
           Note                         : No isotope restrictions are herein
                                          considered

   Look for I value. Allowed parameters are the following:

    Iso(idx)         - Spin I type (2H, 131Xe, ...)
    Ion(idx)         - Magnetic Ion type (Ce3+, Ho3+,....)
    Pbase,  Pbase(#) - Spin I quantum number (0.5, 1.5, 2.0, 2.5,..)         */


// ----------------------------------------------------------------------------
//                          Single Spin Spin Tensors
// ----------------------------------------------------------------------------

/*     Function                                Return
   ================  ==========================================================
       getIso        Finds string for Iso(idx) if found or "" if not found
       getIqn        Returns double for Pbase(idx) if found or 0 if not found
      SpinCheck      Insures that a string signifies a valid isotope type
      SpinCheck      Insures spin for epr is electron or for nmr is nucleus
      SpinCheck      Insures that double signifies a valid quantum number
*/
bool getIso(const ParameterSet& pset, std::string& IS,
                                             int idx=-1, bool warn=true) const;
bool scanIso(const ParameterSet& pset, std::string& IS,
                                             int idx=-1, bool warn=true) const;

bool getIqn(const ParameterSet& pset, const std::string& Pbase,
                                 double& qn, int idx=-1, bool warn=true) const;
bool scanIqn(const ParameterSet& pset, const std::string& Pbase,
                                 double& qn, int idx=-1, bool warn=true) const;

bool SpinCheck(const std::string&  II,          bool warn=true) const;
bool SpinCheck(const Isotope& II,bool epr=false,bool warn=true) const;
bool SpinCheck(      double   qn,bool qud=false,bool warn=true) const;

// ----------------------------------------------------------------------------
//                          Spin Pair Spin Tensors
// ----------------------------------------------------------------------------

/* These functions try and set up a spin tensor for a spin-pair interaction
   from parameters in a parameter set. The functions are set virtual because
   they bypass the static list of spin tensor, hence any derived classes may
   wish to use their own versions in order to take advantage of the list.

       Function                                Return
   ================  ==========================================================
       getIsos       Finds strings for Iso(idxI) & Isos(IdxS), T if both found
       getIqns       Finds doubles for PbaseI(idx) & PBaseS(idx), T if found
      SpinCheck      Insures strings signify valid isotope types
      SpinCheck      Insures spin pair for epr is electron with nucleus 
                     or that spin pair for nmr is nucleus/nucleus or e-/e-
      SpinCheck      Insures doubles signify a valid quantum numbers
      SpinCheck      Insures indicies a proper spin pair

*/
bool getIsos(const ParameterSet& pset, int idxI, int idxS,
                       std::string& II, std::string& IS, bool warn=true) const;
bool scanIsos(const ParameterSet& pset, int idxI, int idxS,
                         std::string& II, std::string& IS, bool warn=true) const;
bool getIqns(const ParameterSet& pset, const std::string& Pbase,
                       double& Iqn, double& Sqn, int i=-1, bool warn=true) const;
bool getIqns(const ParameterSet& pset, const std::string& Pbase,
             double& Iqn, double& Sqn, int idxI, int idxS, bool warn=true) const;
bool scanIqns(const ParameterSet& pset, const std::string& Pbase,
                       double& Iqn, double& Sqn, int i=-1, bool warn=true) const;

bool SpinCheck(const std::string&  II, const std::string&  IS, bool warn=true) const;
bool SpinCheck(const Isotope& II, const Isotope& IS, bool epr, bool warn=true) const;
bool SpinCheck(      double   Iq,       double   Sq,           bool warn=true) const;
bool SpinCheck(      int      Ii,       int      Si,           bool warn=true) const;


//sosik




bool setTIso(const ParameterSet& pset,
                                    int ttype=0, int idx=-1, int warn=1);

bool setTqn(const ParameterSet& pset, const std::string& PBase,
                                          int ttype=0, int idx=-1, int warn=1);


bool setT(const ParameterSet& pset, const std::string& PBase,
                                          int ttype=0, int idx=-1, int warn=1);

bool setTIsos(const ParameterSet& pset,
                                    int ttype, int idxI, int idxS, int warn=1);

bool setTqns(const  ParameterSet& pset, const std::string& PBase,
                                                           int idxS, int warn);

bool setTqns(const  ParameterSet& pset, const std::string& PBase,
                                                 int idxI, int idxS, int warn);

bool setTs(const    ParameterSet& pset, const std::string& PBase,
                                                int ttype, int idxI, int warn);

bool setTs(const    ParameterSet& pset, const std::string& PBase,
                                      int ttype, int idxI, int idxS, int warn);

// ____________________________________________________________________________
// iv           RANK 2 SPIN TENSOR HILBERT SPACE TRANSFORMATIONS
// ____________________________________________________________________________

/* These functions will take a spin tensor component (matrix) defined in the
   single spin or spin pair Hilbert space and convert it into an array defined
   in any specified composite Hilbert space. Individual Hilbert spaces in the
   specified composite must match those of the spins involved in the tensor. */

static std::vector<int> HSorder(const  std::vector<int>& HSs, int i);
MSVCDLL matrix blow_up(const matrix& mx, const std::vector<int>& HSs,
                                                        int i, int j=-1) const;

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A         INTERACTION SPIN TENSOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

MSVCDLC IntRank2T();
MSVCDLC IntRank2T(const IntRank2T &IST1);

// ----------------------------------------------------------------------------
//         Direct Constructors For Interactions Involving A Single Spin
// ----------------------------------------------------------------------------

        // Input                IST     : Interaction Spin Tensor (this)
        //                      IsoI    : Spin I isotope type
        //                      IsoI    : Spin isotope
	//			Iqn	: Spin I quantum number
	//			IST_type: Interaction type (SPF, SA, G, QUAD)
	//			ttype   : Spin tensor type
	//			  	   true  = Spin Field (SA, G)
	//				   flase = Spin Self (Quad)
        //                      warn    : Warning output label
        //                                 0 = no warnings
	//				   1 = non-fatal warnings
        //                                >1 = warnings & exit
        // Output               none    : Interaction Spin Tensor constructed

MSVCDLC IntRank2T(const std::string&  IsoI, bool ttype=true, int warn=2);
MSVCDLC IntRank2T(const      Isotope& IsoI, bool ttype=true, int warn=2);
MSVCDLC IntRank2T(double Iqn,               bool ttype=true, int warn=2);

// ----------------------------------------------------------------------------
//         Direct Constructors For Interactions Involving A Spin Pair
// ----------------------------------------------------------------------------

        // Input                IST     : Interaction Spin Tensor (this)
        //                      IsoI    : Spin I isotope type
        //                      IsoS    : Spin S isotope type
	//			Iqn	: Spin I quantum number
	//			Sqn	: Spin S quantum number
	//			IST_type: Interaction type (DIP, HF)
        // Output               none    : Interaction Spin Tensor constructed


MSVCDLC IntRank2T(const std::string& IsoI, const std::string& IsoS);
MSVCDLC IntRank2T(const Isotope&     IsoI, const Isotope&     IsoS);
MSVCDLC IntRank2T(double             Iqn,  double             Sqn );


// ----------------------------------------------------------------------------
//                Now For The Assignment Operator & Destructor
// ----------------------------------------------------------------------------

MSVCDLL void    operator= (const IntRank2T &IST1);
MSVCDLC virtual ~IntRank2T();

// ____________________________________________________________________________
// B       INTERACTION SPIN TENSOR SPIN TENSOR I & S VALUE ACCESS
// ____________________________________________________________________________

/*      Input                	IST 	: Interaction Spin Tensor
	Output			I,S     : Quantum value of I/S (e.g. 0.5, 1.5)
				IV,SV   : Hilbert space of I/S (e.g. 2, 6)
				HS      : Spin Hilbert space of IST          */

MSVCDLL double Izval()           const;	// Iz: 0.5, 1.0, 1.5, ...
MSVCDLL double Szval()           const;	// Sz: 0.5, 1.0, 1.5, ...
MSVCDLL int    IV()              const;	// 2I+1: 2, 3, 4, ...
MSVCDLL int    SV()              const;	// 2S+1: 2, 3, 4, ...
MSVCDLL int    HS()              const;	// 2I+1 or (2I+1)*(2S+1)
MSVCDLL double qn(bool TF=false) const;	// Iz or Sz value

// ____________________________________________________________________________
// C       INTERACTION SPIN TENSOR SPIN TENSOR COMPONENT ACCESS
// ____________________________________________________________________________

/* These functions provide direct access to individual spin tensor spherical
   components.  Note that if the spin tensor has not been initialized each
   component will be a NULL matrix.  The index spans m = [-2, 2]. Note that the
   functions will return arrays in the single or spin pair Hilbert space
   unless a composite Hilbert space is provided.                             */

// ----------------------------------------------------------------------------
//        Single Spin Or Spin Pair Hilbert Space Spin Tensor Components
// ----------------------------------------------------------------------------

MSVCDLL matrix T20()           const;		// Return T2,0 
MSVCDLL matrix T21()           const;		// Return T2,1
MSVCDLL matrix T2m1()          const;		// Return T2,-1
MSVCDLL matrix T22()           const;		// Return T2,2
MSVCDLL matrix T2m2()          const;		// Return T2,-2
MSVCDLL matrix T2m(int m)      const;		// Return T2,m
MSVCDLL matrix Tcomp(int comp) const;		// Return T2,m

// ----------------------------------------------------------------------------
//               Composite Hilbert Space Spin Tensor Components
// ----------------------------------------------------------------------------

MSVCDLL matrix T20(       const std::vector<int>& HSs,int i,int j=-1) const;
MSVCDLL matrix T21(       const std::vector<int>& HSs,int i,int j=-1) const;
MSVCDLL matrix T2m1(      const std::vector<int>& HSs,int i,int j=-1) const;
MSVCDLL matrix T22(       const std::vector<int>& HSs,int i,int j=-1) const;
MSVCDLL matrix T2m2(      const std::vector<int>& HSs,int i,int j=-1) const;
MSVCDLL matrix T2m(int m, const std::vector<int>& HSs,int i,int j=-1) const;
MSVCDLL matrix CompSpace(const matrix& mx,
                             const std::vector<int>& HSs,int i,int j=-1) const;

// ____________________________________________________________________________
// D           INTERACTION SPIN TENSOR SPIN CHECKING FUNCTIONS
// ____________________________________________________________________________

/* These functions perform checks on the rank 2 spin tensor.
 
        Function                Purpose
   ====================
   SAspincheck             Iqn     True if I is an even multiple of 1/2
   Gspincheck              Iqn     True if I is an even multiple of 1/2
   Qspincheck              Iqn     True if I is greater than 1/2
   Dspincheck            Iqn,Sqn   True if I & S non-zero + multiples of 1/2

           Input                IST     : Interaction spin tensor (HF?)
                                Iqn,Sqn : Spin quantum values of I&S
           Output               TF      : True if both I & S have spin
                                          quantum values that are non-zero
                                          positive multipes of 1/2
                                          and at least 1 is of spin I=1/2    */

//------------------------------------------------------------------------------
//                        Isotope Designation Checks
//------------------------------------------------------------------------------
// ____________________________________________________________________________
// F                          STANDARD OUTPUT FUNCTIONS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//     Functions That Generate Strings to Simplify and Modularize Printing
//-----------------------------------------------------------------------------

MSVCDLL std::string StringI()  const;
MSVCDLL std::string StringS()  const;
MSVCDLL std::string StringIS() const;

//-----------------------------------------------------------------------------
//  Functions That Generate string Arrays to Simplify and Modularize Printing
//-----------------------------------------------------------------------------

/* These functions return spin tensor components in string format.  This is
   done to facilitate printing, in particular printing of spin tensors within
   rank 2 interactions.

           Input                IST     : Interaction Spin Tensor
                                M       : Ang. momentum component [0,4]
                                m       : Ang. momentum component [-2,2]
           Output               TSS     : Pointer to array of hs strings
                                          where hs is the spin pair
                                          Hilbert space dimension

                                [ x.x, x.x, x.x]
                        T     = [ x.x, x.x, x.x]
                         2,m    [ x.x, x.x, x.x]
                                [ x.x, x.x, x.x]

              M = { 0, 1, ..., 4 } <===> m = { 0, 1, -1, 2, -2 }             */


MSVCDLL std::vector<std::string> T2Strings(int m) const;
MSVCDLL std::vector<std::string> TStrings(int  M) const;

//-----------------------------------------------------------------------------
//     Functions That Generate Ouput Of The Rank 2 Interaction Spin Tensor
//-----------------------------------------------------------------------------

/* 	   Input                IST 	: Interaction Spin Tensor (this)
                                os	: Output stream
                                fflag   : Format flag
                                            0 - Basic Parameters
                                           !0 - Full output
           Output               none    : Spin tensor components
                                          placed into the output stream      */

MSVCDLL virtual std::ostream& print(std::ostream& os, int fflag=-1) const;
MSVCDLL friend  std::ostream& operator<< (std::ostream&      os, const IntRank2T& IST);

// ____________________________________________________________________________
// G               RANK 2 TENSOR LIST/VECTOR SUPPORT FUNCTIONS
// ____________________________________________________________________________

/*         F_list ==                    - Equality
           F_list !=                    - Inequality
   	   Input                IST     : Interaction Spin Tensor (this)
                                IST1    : Another interaction spin tensor
                                Iqn,Sqn : Spin quantum values of I&S

              Function                      Return
             ==========      ======================================
             operator==      True if I.qn() == Iqn && S.qn() =  Sqn
             operator!=      True if I.qn() != Iqn || S.qn() != Sqn
             operator<       True if I.qn() <  Iqn || S.qn() <  Sqn
             lessthan        True if I.qn() <  Iqn && S.qn() <  Sqn
             morethan	True if I.qn() >  Iqn && S.qn() >  Sqn
             equalto         True if I.qn() =  Iqn && S.qn() =  Sqn          */

MSVCDLL virtual bool operator==(const IntRank2T &IST1)   const;
MSVCDLL virtual bool operator!=(const IntRank2T &IST1)   const;
MSVCDLL virtual bool operator<(const  IntRank2T &IST1)   const;
MSVCDLL virtual bool operator>(const  IntRank2T &IST1)   const;
MSVCDLL         int  lessthan(double  Iqn, double Sqn=0) const;
MSVCDLL         int  morethan(double  Iqn, double Sqn=0) const;
MSVCDLL         int  equalto(double   Iqn, double Sqn=0) const;
  };

#endif						// IntRank2T.h
