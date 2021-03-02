/* solid_sys.h **************************************************-*-c++-*-
**									**
** 	                         G A M M A				**
**			         				 	**
**	Solid Spin System 			       Interface	**
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
** A solid spin system is a collection of spin isotopes each of	which	**
** has associated tensor properties (shift, quadrupole)	and a relative	**
** position in 3-dimensional space (dipole).  Thus the solid system 	**
** corresponds to a single molecule (or spin network, or crystallite) 	**
** with a set orientation.						**
**						 			**
*************************************************************************/

///Chapter Class Solid Spin System (solid_sys)
///Section Overview
///Body    The class 
///Section Available Solid Spin System Functions

#ifndef   Solid_sys_h_			// Is this file already included?
#  define Solid_sys_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <HSLib/SpinSystem.h>		// Base spin system class
#include <IntRank2/IntDip.h>		// Include dipolar interactions
//#include <IntRank2/IntCSA.h>		// Include shift anis. interactions
//#include <IntRank2/IntQuad.h>		// Include Quad interactions
//#include <IntRank2/IntG.h>		// Include G interactions
//#include <IntRank2/IntHF.h>		// Include Hyperfine interactions
#include <Level1/coord_vec.h>		// Include coordinate vectors
#include <Basics/ParamSet.h>		// Include GAMMA parameters
#include <Basics/Isotope.h>		// Include spin isotopes
#include <IntRank2/IntDipVec.h>         // Include dipolar interaction vectors
#include <IntRank2/IntCSAVec.h>         // Include CSA interaction vectors
#include <IntRank2/IntQuadVec.h>	// Include quad. interaction vectors
#include <IntRank2/IntGVec.h> 		// Include G   interaction vectors
#include <IntRank2/IntHFVec.h> 		// Include HF  interaction vectors
#include <string>			// Include stdlibc++ strings


class solid_sys: public spin_system 
  {
  IntQuadVec       Qvec;		// Array of quadrupolar interactions
  IntDipVec        Dvec;		// Vector of dipolar interactions
  IntCSAVec        Cvec;		// Vector of CSA interactions
  IntGVec          Gvec;		// Vector of electron G interactions
  IntHFVec         HFvec;		// Vectro of hyperfine interactions
  coord_vec        SCoords;		// Array of spin coordinates
  std::vector<int> cflags;		// Array of coordinate existence flags
  
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                  CLASS SOLID SPIN SYSTEM ERROR HANDLING
// ____________________________________________________________________________

/*       Input                ssys    : Solid spin system (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

void ssys_error(int eidx, int noret=0) const;
void ssys_error(int eidx, const std::string& pname, int noret=0) const;
volatile void ssys_fatal(int eidx) const;
volatile void ssys_fatal(int eidx, const std::string& pname) const;

// ____________________________________________________________________________
// ii                      SOLID SYSTEM SETUP FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//            Functions To Set Coordinates and Their Existence Flags
// ----------------------------------------------------------------------------
 
void zero_cindx();
 
        // Input                ssys    : Solid spin system (this)
        // Output               none    : If there are spins in the system
        //                                this zeros the coordinate existence
        //                                flags.  If the flags don't exist,
        //                                then the array is constructed.  If
        //                                no spins, the array is set to NULL.

     
int setCoords(const ParameterSet& pset);
 
        // Input                ssys    : Solid spin system (this)
        //                      pset    : A parameter set
        // Output               T/F     : True if spin coordiantes are found in
        //                                the pset & they're read into ssys.
        // Note                         : These same coordinates should be used
        //                                to construct the dipolar tensors
        // Note                         : This assumes that the array cflags
        //                                has been properly allocated


// ----------------------------------------------------------------------------
//                   Functions To Set Dipolar Interactions
// ----------------------------------------------------------------------------
 
/* For this spin system we (should) know the number of spins and their isotope
   types.  The dipolar interactions are all to be placed into the dipolar
   interaction vector, Dvec.  This must be read in using spin (NOT interaction)
   indices.  For example, { Iso(i), coord(i), Iso(j), coord(j) } for all spin
   pairs.  However, not everyone likes to use spin coordinates for designating
   a dipolar interaction, so we must allow for their setting without the use
   of coordinates, but using { Isotopes, spin indices, + other info }        */  
 

void setDs(const ParameterSet& pset, int ccount);

        // Input                ssys    : Solid spin system (this)
        //                      pset    : A parameter set
        // Output               none    : System dipolar interactions are
        //                                set from parameters in pset
 
// ----------------------------------------------------------------------------
//                Functions To Set Shift Anisotropy Interactions
// ----------------------------------------------------------------------------

void setCs(const ParameterSet& pset);

        // Input                ssys    : Solid spin system (this)
        //                      pset    : A parameter set
        // Output               none    : System shift anisotropy interactions
        //                                are set from parameters in pset
        // Note                         : It isn't mandatory that there be
        //                                any shift interactions


// ----------------------------------------------------------------------------
//                Functions To Set Quadrupolar Interactions
// ----------------------------------------------------------------------------


void setQs(const ParameterSet& pset);
 
        // Input                ssys    : Solid spin system (this)
        //                      pset    : A parameter set
        // Output               none    : Solid spin system quadrupolar
        //                                interactions are set from parameters
        //                                in pset

// ----------------------------------------------------------------------------
//                  Functions To Set Electron G Interactions
// ----------------------------------------------------------------------------

/* Each Spin That Is An Electron May Have An Electron G Tensor Interaction Set.
   We Don't Particularly Care If There Isn't A G Interaction Defined For Any
   Particular Electron Spin, The Vector Of Such Interactions Will Simply Have
   An Zero Values If No G Anisotropy Has Been Specified.

           Input                ssys    : Solid spin system (this)
                                pset    : A parameter set
           Output               none    : System electron G interactions
                                          are set from parameters in pset
           Note                         : It isn't mandatory that there be
                                          any electron G interactions
           Note                         : Isotopes MUST be set prior to here.
                                          Isotopes used to avoid confusion in
                                          IntG over possible nuclear spins   */

void setGs(const ParameterSet& pset);


// ----------------------------------------------------------------------------
//                  Functions To Set Hyperfine Interactions
// ----------------------------------------------------------------------------

/* Each Spin Pair That Involves An Electron And A Nucleon May Have A Hyperfine
   Interaction Set. We Don't Particularly Care If There Isn't A HF Interaction
   Defined For Any Particular Spin Pair, The Vector Of Such Interactions Will
   Simply Have An Zero Values If No HF Anisotropy Has Been Specified.

           Input                ssys    : Solid spin system (this)
                                pset    : A parameter set
           Output               none    : System hyperfine interactions
                                          are set from parameters in pset
           Note                         : It isn't mandatory that there be
                                          any hyperfine interactions
           Note                         : Isotopes MUST be set prior to here.
                                          Isotopes used to avoid confusion
                                          over electron vs nuclear spins     */

void setHFs(const ParameterSet& pset);


// ----------------------------------------------------------------------------
//                  Functions To Set The Full Spin System
// ----------------------------------------------------------------------------


int setSsys(const ParameterSet& pset, int indx=-1, int warn=2);

        // Input                ssys    : Solid spin system (this)
        //                      pset    : A parameter set
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning output label
        //                                 0 = no warnings
        //                                 1 = warnings
        //                                >1 = fatal warnings
        // Output               TF      : The entire system is set
        //                                from parameters in pset
        // Note                         : This uses the assignment from pset
        //                                It exists to support prefix indices


// ____________________________________________________________________________
// iii                 SOLID SYSTEM ADMINISTRATION FUNCTIONS
// ____________________________________________________________________________


void ResetSOps(int spin);
 
        // Input                ssys    : Solid spin system (this)
        // Output               none    : Resets any spin operators for
        //                                specified spin


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

// ____________________________________________________________________________
// A               SOLID SPIN SYSTEM CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

///Center Basic Functions
///F_list =		      - Assignment


MSVCDLC      solid_sys(int spins=0);
MSVCDLC      solid_sys(const solid_sys &ssys1);
MSVCDLL void operator= (const solid_sys &ssys);
MSVCDLC      ~solid_sys ();

// ____________________________________________________________________________
// B             SPIN COORDINATES VECTOR ACCESS & MANIPULATIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                          Get/Set Spin Coordinates
// ----------------------------------------------------------------------------

        // Input                ssys    : Solid spin system
	//			i	: spin index
        //                      pt      : spin coordinate
        //                      cvec    : Spin coordinates
        // Output               coord	: Spin i coordinate if requested
        //                      cvec    : Spin coordinates of all spins
        //                      none    : Spin coordinate set to pt
        //                      none    : Spin coordinates set to cvec
        // Note                         : Coordinate can only be set if there
        //                                already exists a complete set!
        // Note                         : The size of cvec must be at least as
        //                                big as the number of system spins
 
MSVCDLL coord     getCoord(int i) const;
MSVCDLL coord_vec getCoords()     const;
MSVCDLL void      setCoord(int i, coord& pt);
MSVCDLL void      setCoords(const coord_vec& cvec);
 

// ____________________________________________________________________________
// C                DIPOLAR INTERACTION ACCESS & MANIPULATIONS
// ____________________________________________________________________________

/* The functions in this section provide access to the dipolar interactions in
   the system.

           Input                ssys    : Solid spin system
                                sI      : Spin I of spin system [0, nspins-1]
                                sS      : Spin S of spin system [0, nspins-1]
                                val     : Dipolar interaction value
                                          <=0: DCC coupling of spin (Hertz)
                                            1: Asymmetry [0, 1]
                                            2: Theta (degrees, down from +z)
                                            3: Phi   (degrees, over from +z)
           Output               none    : Get/Set dipolar interaction value
                                          for specified spin pair

    Function  Arguments                        Result
    ========  =========  ======================================================
      DCC      i,j,dcc   Sets the dipolar coupling constant to dcc (Hz)
      DCC        i,j     Returns the dipolar coupling constant to dcc (Hz)
     Ddelz     i,j,dcc   Sets the dipolar coupling constant to dcc (Hz)
     Ddelz       i,j     Returns the dipolar coupling constant to dcc (Hz)
      DEta     i,j,dcc   Sets the dipolar interaction asymmetry [0, 1]
      DEta       i,j     Returns the dipolar interaction asymmetry [0,1]
     DTheta    i,j,the   Sets the dipolar interaction theta orientation (deg)
     DTheta      i,j     Returns the dipolar interaction theta orientat. (deg)
      DPhi     i,j,phi   Sets the dipolar interaction phi orientation (deg)
      DPhi       i,j     Returns the dipolar interaction phi orientat. (deg)
     IntDip      i,j     Returns the rank 2 dipolar interaction between i & j
    IntDipVec            Returns a vector of all system dipolar interactions

   Note that DCC = delzz for dipolar interactions and that setting DCC will
   kill any coordinates associated between two spins. The dipolar interaction
   asymmetry is typically zero. The theta orientation is that down from the
   lab frame to the PAS +z axis and is restricted to [0,180]. The phi 
   orientation is from the lab +x axis to the PAX +x axis and it is restricted
   to [0, 360].  Note that any changes to theta and phi will kill any spin
   coordinates. 
	*/

// ----------------------------------------------------------------------------
//                   Generic Dipolar Value Access Functions
// ----------------------------------------------------------------------------

MSVCDLL void   DValue(int sI, int sS, double val, int type);
MSVCDLL double DValue(int sI, int sS,             int type) const;

//                         Dipolar Coupling Constants

MSVCDLL void   DCC(int   spinI, int spinS, double dcc);
MSVCDLL double DCC(int   spinI, int spinS) const;
MSVCDLL double Ddelz(int spinI, int spinS) const;
MSVCDLL void   Ddelz(int spinI, int spinS, double delzz);
   
//                        Dipolar Asymmetry Values

MSVCDLL double Deta(int spinI, int spinS) const;
MSVCDLL void   Deta(int spinI, int spinS, double ETA);

//               Dipolar Theta Orientation (Down From +z Axis)

MSVCDLL double Dtheta(int spinI, int spinS) const;
MSVCDLL void   Dtheta(int spinI, int spinS, double dtheta);

//               Dipolar Phi Orientation (Over From +x Axis)
 
MSVCDLL double Dphi(int spinI, int spinS) const;
MSVCDLL void   Dphi(int spinI, int spinS, double dphi);

//                    Dipolar Spin Tensor Component
         
matrix DTcomp(int spinI, int spinS, int m) const;

//                        Full Dipolar Interaction
         
MSVCDLL const IntDip&    getDipInt(int spinI, int spinS) const;
MSVCDLL const IntDip&    getDipInt(int dip)              const;
MSVCDLL const IntDipVec& getDipVec()                     const;
 
// ____________________________________________________________________________
// D                CSA INTERACTION ACCESS & MANIPULATIONS
// ____________________________________________________________________________

/* The functions in this section provide access to the shift anisotropy 
   interactions in the system.

           Input                ssys    : Solid spin system
                                sI      : Spin I of spin system [0, nspins-1]
                                val     : Shift anisotropy interaction value
                                          <=0: CSA of spin (PPM)
                                            1: Asymmetry [0, 1]
                                            2: Theta (degrees, down from +z)
                                            3: Phi   (degrees, over from +z)
           Output               none    : Get/Set shift anisotropy interaction
                                          value for specified spin

    Function  Arguments                        Result
    ========  =========  ======================================================
      DCC       i,dcc    Sets the dipolar coupling constant to dcc (Hz)
      DCC         i      Returns the dipolar coupling constant to dcc (Hz)
     Ddelz      i,dcc    Sets the dipolar coupling constant to dcc (Hz)
     Ddelz        i      Returns the dipolar coupling constant to dcc (Hz)
      CEta      i,dcc    Sets the dipolar interaction asymmetry [0, 1]
      CEta        i      Returns the dipolar interaction asymmetry [0,1]
     CTheta     i,the    Sets the dipolar interaction theta orientation (deg)
     CTheta       i      Returns the dipolar interaction theta orientat. (deg)
      CPhi      i,phi    Sets the dipolar interaction phi orientation (deg)
      CPhi        i      Returns the dipolar interaction phi orientat. (deg)
     IntCSA       i      Returns the rank 2 dipolar interaction between i & j
    IntCSAVec            Returns a vector of all system dipolar interactions

   Note that the chemical shift anisotropy (CSA) is defined to be 1.5 times
   the interaction delzz value.  The theta orientation is that down from the
   lab frame to the PAS +z axis and is restricted to [0,180]. The phi 
   orientation is from the lab +x axis to the PAX +x axis and it is restricted
   to [0, 360]. 
*/

//                   Generic CSA Value Access Functions

MSVCDLL void   CValue(int spin, double val, int type);
MSVCDLL double CValue(int spin,             int type) const;

//                            CSA & delzz Values
 
MSVCDLL void   CSA(int   sI, double cs);
MSVCDLL double CSA(int   sI) const;
MSVCDLL void   Cdelz(int sI, double dz);
MSVCDLL double Cdelz(int sI) const;
 
//                            CSA Asymmetry Values

MSVCDLL void   Ceta(int sI, double CE);
MSVCDLL double Ceta(int sI) const;
     
//                   CSA Theta Orientation (Down From +z Axis)

MSVCDLL void   Ctheta(int sI, double ctheta);
MSVCDLL double Ctheta(int sI) const;

//                    CSA Phi Orientation (Over From +x Axis)
 
MSVCDLL void   Cphi(int sI, double cphi);
MSVCDLL double Cphi(int sI) const;

//                             Full CSA Interaction
 
MSVCDLL IntCSA getCSAInt(int sI) const;

// ____________________________________________________________________________
// E             QUADRUPOLAR INTERACTION ACCESS & MANIPULATIONS
// ____________________________________________________________________________

/* The functions in this section provide access to the quadrupolar interactions
   in the system.

           Input                ssys    : Solid spin system
                                sI      : Spin I of spin system [0, nspins-1]
                                val     : Quadrupolar interaction value
                                          <=0: QCC of spin
                                            1: Asymmetry [0, 1]
                                            2: Theta (degrees, down from +z)
                                            3: Phi   (degrees, over from +z)
           Output               none    : Get/Set quadrupolar interaction
                                          value for specified spin

    Function  Arguments                        Result
    ========  =========  ======================================================
      QCC      i,j,dcc   Sets the quadrupolar coupling constant to dcc (Hz)
      QCC        i,j     Returns the quadrupolar coupling constant to dcc (Hz)
     Qdelz     i,j,dcc   Sets the quadrupolar coupling constant to dcc (Hz)
     Qdelz       i,j     Returns the quadrupolar coupling constant to dcc (Hz)
      QEta     i,j,dcc   Sets the quadrupolar interaction asymmetry [0, 1]
      QEta       i,j     Returns the quadrupolar interaction asymmetry [0,1]
     QTheta    i,j,the   Sets the quadrupolar interaction theta orientation (deg)
     QTheta      i,j     Returns the quadrupolar interaction theta orientat. (deg)
      QPhi     i,j,phi   Sets the quadrupolar interaction phi orientation (deg)
      QPhi       i,j     Returns the quadrupolar interaction phi orientat. (deg)
     IntQuad     i,j     Returns the rank 2 quadrupolar interaction for spin i
    IntQuadVec           Returns a vector of all system quadrupolar interactions

   Note that the quadrupolar coupling (QCC) is defined to be identical to the
   interaction delzz value.  The theta orientation is that down from the
   lab frame to the PAS +z axis and is restricted to [0,180]. The phi 
   orientation is from the lab +x axis to the PAX +x axis and it is restricted
   to [0, 360]. 
*/

//                   Generic Quadrupolar Value Access Functions

MSVCDLL void   QValue(int spin, double val, int type);
MSVCDLL double QValue(int spin,             int type) const;

//                            QCC & delzz Values

MSVCDLL void   QCC(int   spin, double qcc);
MSVCDLL double QCC(int   spin) const;
MSVCDLL double Qdelz(int spin) const;
MSVCDLL void   Qdelz(int spin, double delzz);

//                        Quadrupolar Asymmetry Values

MSVCDLL double  Qeta(int   spin) const;
MSVCDLL void    Qeta(int   spin, double Qeta);

//               Quadrupolar Theta Orientation (Down From +z Axis)

MSVCDLL double  Qtheta(int spin) const;
MSVCDLL void    Qtheta(int spin, double Qtheta);

//               Quadrupolar Phi Orientation (Over From +x Axis)

MSVCDLL double  Qphi(int   spin) const;
MSVCDLL void    Qphi(int   spin, double Qphi);

//                    Quadrupolar Spin Tensor Component
         
MSVCDLL matrix QTcomp(int spin, int m) const;

//                        Full Quadrupolar Interaction
         
MSVCDLL const IntQuad&    getQuadInt(int spin) const;
MSVCDLL const IntQuadVec& getQuadVec()         const;
 
// sosi?
MSVCDLL IntQuad Qint(int   spin);
 
// ____________________________________________________________________________
// F              ELECTRON G INTERACTION ACCESS & MANIPULATIONS
// ____________________________________________________________________________

/* The functions in this section provide access to the electron G interactions
   in the system.

           Input                ssys    : Solid spin system
                                sI      : Spin I of spin system [0, nspins-1]
                                val     : Electron G interaction value
                                          <=0: DCC coupling of spin (Hertz)
                                            1: Asymmetry [0, 1]
                                            2: Theta (degrees, down from +z)
                                            3: Phi   (degrees, over from +z)
           Output               none    : Get/Set dipolar interaction value
                                          for specified spin pair

    Function  Arguments                        Result
    ========  =========  ======================================================
      GCC      i,dcc     Sets the electron G coupling constant to dcc (Hz)
      GCC        i       Returns the electron G coupling constant to dcc (Hz)
     Gdelz     i,dcc     Sets the electron G coupling constant to dcc (Hz)
     Gdelz       i       Returns the electron G coupling constant to dcc (Hz)
      GEta     i,dcc     Sets the electron G interaction asymmetry [0, 1]
      GEta       i       Returns the electron G interaction asymmetry [0,1]
     GTheta    i,the     Sets the electron G interaction theta orientation (deg)
     GTheta      i       Returns the electron G interaction theta orientat. (deg)
      GPhi     i,phi     Sets the electron G interaction phi orientation (deg)
      GPhi       i       Returns the electron G interaction phi orientat. (deg)
      IntG       i       Returns the rank 2 electron G interaction between i & j
     IntGVec             Returns a vector of all system electron G interactions

   Note that DCC = delzz for dipolar interactions and that setting DCC will
   kill any coordinates associated between two spins. The dipolar interaction
   asymmetry is typically zero. The theta orientation is that down from the
   lab frame to the PAS +z axis and is restricted to [0,180]. The phi 
   orientation is from the lab +x axis to the PAX +x axis and it is restricted
   to [0, 360].  Note that any changes to theta and phi will kill any spin
   coordinates. 
	*/

MSVCDLL void   GValue(int sI, double val, int type);
MSVCDLL double GValue(int sI,             int type) const;

//                         Electron G Coupling Constants


MSVCDLL void   DCC(int   spinI, double dcc);
MSVCDLL double DCC(int   spinI) const;
MSVCDLL double Gdelz(int spinI) const;
MSVCDLL void   Gdelz(int spinI, double delzz);

//                        Electron G Asymmetry Values

MSVCDLL double Geta(int spinI) const;
MSVCDLL void   Geta(int spinI, double ETA);

//               Electron G Theta Orientation (Down From +z Axis)

MSVCDLL double Gtheta(int spinI) const;
MSVCDLL void   Gtheta(int spinI, double gtheta);

//               Electron G Phi Orientation (Over From +x Axis)

MSVCDLL double Gphi(int spinI) const;
MSVCDLL void   Gphi(int spinI, double gphi);

//                    Electron G Spin Tensor Component

MSVCDLL matrix GTcomp(int spinI, int m) const;

//                        Full Electron G Interaction

MSVCDLL IntG    getGInt(int spinI) const;
MSVCDLL IntGVec getGVec()          const;

// ____________________________________________________________________________
// G               HYPERFINE INTERACTION ACCESS & MANIPULATIONS
// ____________________________________________________________________________


/* The functions in this section provide access to the hyperfine interactions
  in the system.

          Input                ssys    : Solid spin system
                               sI      : Spin I of spin system [0, nspins-1]
                               sS      : Spin S of spin system [0, nspins-1]
                                           1: Asymmetry [0, 1]
                                           2: Theta (degrees, down from +z)
                                           3: Phi   (degrees, over from +z)
          Output               none    : Get/Set hyperfine interaction value
                                         for specified spin pair

   Function  Arguments                        Result
   ========  =========  ======================================================
     DCC      i,j,dcc   Sets the hyperfine coupling constant to dcc (Hz)
     DCC        i,j     Returns the hyperfine coupling constant to dcc (Hz)
    HFdelz    i,j,dcc   Sets the hyperfine coupling constant to dcc (Hz)
    HFdelz      i,j     Returns the hyperfine coupling constant to dcc (Hz)
     HFEta    i,j,dcc   Sets the hyperfine interaction asymmetry [0, 1]
     HFEta      i,j     Returns the hyperfine interaction asymmetry [0,1]
    HFTheta   i,j,the   Sets the hyperfine interaction theta orientation (deg)
    HFTheta     i,j     Returns the hyperfine interaction theta orientat. (deg)
     HFPhi    i,j,phi   Sets the hyperfine interaction phi orientation (deg)
     HFPhi      i,j     Returns the hyperfine interaction phi orientat. (deg)
    IntHF       i,j     Returns the rank 2 hyperfine interaction between i & j
   IntHFVec             Returns a vector of all system hyperfine interactions

  Note that DCC = delzz for dipolar interactions and that setting DCC will
  kill any coordinates associated between two spins. The dipolar interaction
  asymmetry is typically zero. The theta orientation is that down from the
SolidSys.cc" [dos] 2159L, 91494C written
*/

//                   Generic Hyperfine Value Access Functions

MSVCDLL void   HFValue(int sI, int sS, double val, int type);
MSVCDLL double HFValue(int sI, int sS,             int type) const;

//                         Hyperfine Anisotropy Values

MSVCDLL void   HFdelz(int   spinI, int spinS, double hdelz);
MSVCDLL double HFdelz(int   spinI, int spinS) const;
   
//                        Hyperfine Asymmetry Values

MSVCDLL double HFeta(int spinI, int spinS) const;
MSVCDLL void   HFeta(int spinI, int spinS, double ETA);

//               Hyperfine Theta Orientation (Down From +z Axis)

MSVCDLL double HFtheta(int spinI, int spinS) const;
MSVCDLL void   HFtheta(int spinI, int spinS, double dtheta);

//               Hyperfine Phi Orientation (Over From +x Axis)
 
MSVCDLL double HFphi(int spinI, int spinS) const;
MSVCDLL void   HFphi(int spinI, int spinS, double dphi);

//                    Hyperfine Spin Tensor Component
         
MSVCDLL matrix HFTcomp(int spinI, int spinS, int m) const;

//                        Full Hyperfine Interaction
         
MSVCDLL IntHF    getHFInt(int spinI, int spinS) const;


// ____________________________________________________________________________
// H                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________ 

// ---------------------------------------------------------------------------- 
//        Functions To Make A Parameter Set From A Solid Spin System
// ---------------------------------------------------------------------------- 

//Center Parameter Set Functions
///F_list =		       - Conversion
///F_list +=		       - Unary Addition


	// Input		ssys   : Solid spin system
	//  			pset   : Parameter set
	// Output		pset   : Parameter set with
	//			         only system parameters


MSVCDLL             operator ParameterSet( ) const;
MSVCDLL friend void operator+= (ParameterSet& pset, const solid_sys &ssys);

	// Input		ssys   : Solid spin system
	//  			pset   : Parameter set
	// Output		pset   : Parameter set with
	//			         only solid_sys parameters



MSVCDLL void PSetAdd(ParameterSet& pset, int idx=-1) const;

	// Input		ssys	: Solid spin system
        //                      pset	: Parameter set
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        // Output               void    : System parameters are
        //                                are added ot the parameter set
        //                                with interaction index idx
        // Note                         : The parameter names & types
        //                                output here MUST match those used
        //                                in setting the spin system 
        //                                from parameters sets


// ----------------------------------------------------------------------------
//        Functions To Make A Solid Spin System From A Parameter Set
// ----------------------------------------------------------------------------

 
MSVCDLL void operator= (const ParameterSet& pset);

	// Input		ssys   : Solid spin system (this)
	// 			pset   : A parameter set
	// Output		none   : Solid spin system filled with
	//				 parameters from pset


// ----------------------------------------------------------------------------
//    Functions To Output Solid State System To ASCII From A Parameter Set
// ----------------------------------------------------------------------------

MSVCDLL int write(const std::string &filename, int idx=-1, int warn=2) const;
 
        // Input                ssys    : Solid spin system (this)
        //                      filename: Output file name
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning level
        // Output               none    : Solid state system is written
        //                                as a parameter set to file filename


MSVCDLL int write(std::ofstream& ofstr, int idx=-1, int warn=2) const;

        // Input                ssys    : Solid spin system (this)
        //                      ofstr   : Output file stream
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning level
        // Output               none    : Solid state system is written as a
        //                                parameter set to output filestream
        // Note                         : This depends on function PSetAdd!
 
         
// ____________________________________________________________________________
// F               DIPOLAR INTERACTION VECTOR INPUT FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//   Direct Read of Solid State System From An ASCII File Or A Parameter Set
// ----------------------------------------------------------------------------

/* The following read functions utilize a single spin system index, for system
   being read.  They'll will try to read the spin system parameters that are
   prefixed by [#] where #=idx.  The default value, idx=-1, flags that there
   is no prefix on the spin system paramters.

           Input                ssys    : Solid spin system (this)
                                filename: Input filename
                           or   pset    : Parameter set
                                idx     : Parameter index value used for
                                          prefix [#] in input names
                                warn    : Warning output label
                                           0 = no warnings
                                           1 = warnings
                                          >1 = fatal warnings
           Output               TF      : Solid spin system filled with
                                          parameters read from file or pset
                                          Returns true if read properly      */


MSVCDLL int read(const std::string &filename,   int idx=-1, int warn=2);
MSVCDLL int read(const ParameterSet& pset, int idx=-1, int warn=2);
 
// ----------------------------------------------------------------------------
//          Interactive Read of Solid State System From An ASCII File
// ----------------------------------------------------------------------------

/* This function will ask the user for the ASCII (GAMMA parameter set) file
   containing the spin system parameters. There are simply too many to ask for
   individually.

           Input                ssys    : Solid spin system (this)
                                argc    : Number of arguments
                                argv    : Vector of argc arguments
                                argn    : Argument index
           Output               void    : The parameter argn of array argc
                                          is used to supply a filename
                                          from which the spin system is read
                                          If the argument argn is not in argv,
                                          the user is asked to supply a filename
           Note                         : File must be ASCII containing
					  recognized sys parameters
           Note                         : Spin system is modifed (filled)    */

MSVCDLL std::string ask_read(int argc, char* argv[], int argn);
MSVCDLL std::string ask_read(int argc, char* argv[], int argn, const std::string& def);

// ____________________________________________________________________________
//                           STANDARD I/O FUNCTIONS
// ____________________________________________________________________________

///Center Solid Spin System I/O Functions
///F_list print		         - Write system to output stream
///F_list << 		         - Standard Output

/* These functions output either part or the entire solid state spin system.
   Many of these functions will print out lists of the interactions that are
   defined in the spin system as stored in the interaction vectors.  

   Function                                 Output
   --------  ------------------------------------------------------------------
   printPs   Nuclear coordinates (for dipolar interactions).  This function is
             used instead of class vector output because some spins may have no
             coordinates. Unit flag: 1=Angstroms 2=nmeters 3=meters x=none
   printDs   Dipolar interactions sent to output stream
   printCs   Shift anisotropy interactions sent to output stream
   printQs   Quadrupolar interactions sent to output stream
   printGs   Electorn G interactions sent to output stream
   printHFs  Hyperfine interactions sent to output stream
   print     Entire spin system sent to output stream
     <<      Standard output of entire spin system 
*/

MSVCDLL std::ostream& printPs(std::ostream&  ostr, int units=1) const;
MSVCDLL std::ostream& printDs(std::ostream&  ostr) const;
MSVCDLL std::ostream& printCs(std::ostream&  ostr) const;
MSVCDLL std::ostream& printQs(std::ostream&  ostr) const;
MSVCDLL std::ostream& printGs(std::ostream&  ostr, int pf=0)    const;
MSVCDLL std::ostream& printHFs(std::ostream& ostr, int pf=0)    const;
MSVCDLL std::ostream& print(std::ostream& out) const;
MSVCDLL friend std::ostream& operator<< (std::ostream& out, const solid_sys& ssys);


// ____________________________________________________________________________
// I                         BASE SPIN SYSTEM DEALINGS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                These Are All Inherited As Is From Class spin_sys
// ----------------------------------------------------------------------------

/*
inline int    spins() const;			Get number of spins
inline int    spinpairs() const;		Get number of spin pairs
inline int    HS() const;			Get system Hilbert space             
inline int    HS(int spin) const;		Get spin Hilbert space
inline const  Isotope& isotope(int) const;	Get spin isotope(eg 23Na)
inline double weight(int) const;                Get spin atomic weight
inline string symbol(int spin) const;		Get spin symbol(eg 23Na)
inline double qn(int spin) const;		Get spin I value (eg .5) 
inline double qn() const;			Get system Fz value
inline string element(int spin) const;		Get spin type (eg Carbon)
inline string momentum(int spin) const;		Get spin momentum(eg 1/2)
inline string momentum() const;			Get system momentum -1 -1
inline double gamma(int spin) const;		Get gyromag. ratio s  T
inline double gamma(const string& iso) const;	Get gyromag. ratio of iso
inline row_vector qState(int state) const;	Get vector of state Fz
inline matrix     qStates() const;		Get array of spin Iz
inline double     qnState(int state) const;
inline col_vector qnStates() const;
inline row_vector qnDist() const;		Get vector of state dist.
inline row_vector CoherDist() const;		Get vector of coher. dist.
inline int    homonuclear() const;		Test homonuclear
inline int    heteronuclear() const;		Test heteronuclear
       int    spinhalf() const;			Test all spin 1/2
       int    electrons() const;		Test for electrons
       int    nepair(int i, int j) const;	Test for nucleus/e- pair
       int    enpair(int i, int j) const;	Test for e-/nucleus pair
       int    pairidx(int i, int j) const;	Get spin pair index
inline int    isotopes() const;			Get isotope count
inline string isotopes(int idx) const;		Get isotope spin count
inline int    isotopes(const string& I) const;	Test for isotope type
inline void   flags(int TF);			Set spin flags
inline void   flag(int spin, int TF);		Set spin flag
inline void   flag(const string& isoin, int TF);Set isotope flags
inline void   flag(const Isotope& isoin,int TF);Set isotope flags
inline int    flag(int spin) const;		Set spin flag
void get_flags(int* intarray) const;		Get array of spin flags
void set_flags(int* intarray, int TF) const;	Set spin flags
void set_flag(int* intarray, int spin, int TF) const;
void set_flag(int* intarray, const string& isoin, int TF) const;
void set_flag(int* intarray, const Isotope& isoin, int TF) const;
inline void name(const string& name);		Set system name
inline const  string& name(int i=-1) const;	Get system name
void warnings(int warnf);			Set warning level
int warnings() const;				Get warning level
string IsoDefault();				Get default isotope
void IsoDefault(const string& DI);		Set default isotope
int  getSpins(const ParameterSet& pset);	Read number spins
void setName(const ParameterSet& pset);		Read system name
void setIs(const ParameterSet& pset);	        Read system isotopes         */

// ----------------------------------------------------------------------------
//    These Are All ReDefined Beyond Their Functionality In Class spin_sys
// ----------------------------------------------------------------------------

/*
virutal void     isotope(int, const Isotope&);	Set isotope type 
virtual void     isotope(int, const string&);	Set isotope type
virtual void     PSetAdd(ParameterSet& pset, int idx=-1) const;
virtual int      write(const string &filename, int ix=-1, int wn=2) const;
virtual int      write(ofstream& ofstr, int idx=-1, int warn=2) const;
virtual void     read(const string &filename);	Read spin system
virtual string   ask_read(int argc, char* argv[], int argn);
virtual ostream& print(ostream& out) const;	Print spin system           */

// ----------------------------------------------------------------------------
//          These Are The Class spin_sys Redefined Functions
// ----------------------------------------------------------------------------
 
/* Alteration of a spin isotope type will alter any interactions in which the
   spin is involved. Consequently we must not only readjust the spin Hilbert
   space but we must also adjust the appropriate interaction spin tensors that
   are associated. This function does just that.

           Input                ssys    : Solid spin system (this)
                                Iso     : Isotope type
                                symbol  : Isotope type (e.g. 1H, 13C,...)
           Output               none    : Solid spin system spin isotope
                                          type is switched to Iso            */

MSVCDLL void isotope(int spin, const std::string&  symbol);
MSVCDLL void isotope(int spin, const Isotope& Iso);
MSVCDLL const Isotope& isotope(int spin) const;

};

#endif								// solid_sys.h
