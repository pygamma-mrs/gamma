/* HSprop.h *****************************************************-*-c++-*-
**									**
**	                        G A M M A				**
**                                                                      **
**	Hilbert Space Propagator 		Interface		**
**                                                                      **
**      Copyright (c) 1998                                              **
**      Dr. Scott A. Smith                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**								 	**
** Description							 	**
**								 	**
** The class HSprop defines a propagator in spin Hilbert space.  Such	**
** propagators are simply Hilbert space operators which will evolve a	**
** density operator for a specific length of time in a particular set	**
** of rotating frames.							**
**									**
*************************************************************************/

///Chapter Class Hilbert Space Propagator 
///Section Overview
///Body    None
///Section Available Spin System Functions

#ifndef   HSprop_h_			// Is file already included?
#  define HSprop_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <HSLib/GenOp.h>		// Know operator information
#include <vector>			// Know libstdc++ STL vectors


MSVCDLL extern gen_op prop(const gen_op& ham, const double time);
MSVCDLL extern void   prop_ip(   gen_op& ham, const double time);
MSVCDLL extern gen_op evolve(const gen_op& sigma, const gen_op& ham, double time);
MSVCDLL extern void   evolve_ip(   gen_op& sigma, const gen_op& ham, double time);
MSVCDLL extern gen_op evolve(const gen_op& sigma, const gen_op& U);
MSVCDLL extern void   evolve_ip(   gen_op& sigma, const gen_op& U);

class HSprop
  {
  gen_op UOp;				// Propagator operator
  double Ut;				// Evolution time
  matrix Hceil;				// Hamiltonian ceiling array
//  std::vector<double> Frames;		// Rotating frames

 
private:
 
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i              CLASS HILBERT SPACE PROPAGATOR ERROR HANDLING
// ____________________________________________________________________________

/*       Input                U	      : HS propagator (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

void          HSPerror(int eidx, int noret=0) const;
void volatile HSPfatal(int eidx)              const;
 
// ____________________________________________________________________________
// ii            CLASS HILBERT SPACE PROPAGATOR CEILING MATRIX
// ____________________________________________________________________________

/* This sets the array Hceil which only is useful when attempting to convert
   a propagator back to an effective Hamitonian. Even worse, it is only
   useful in some instances when attempting that.                            */
                                                                
void SetCeiling(const gen_op& H, bool I=false);

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
public:
 
// ____________________________________________________________________________
// A         CLASS HILBERT SPACE PROPAGATOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

MSVCDLC HSprop();					// Null constructor
MSVCDLC HSprop(int HS);					// Identity constructor
MSVCDLC HSprop(const gen_op& H, double tevol);
MSVCDLC HSprop(const gen_op& H, double tevol, bool prop);
MSVCDLC HSprop(const HSprop& U);

MSVCDLC ~HSprop();
MSVCDLL HSprop& operator= (const HSprop& U1);

// ____________________________________________________________________________
// B               HILBERT SPACE PROPAGATOR ACCESS FUNCTIONS
// ____________________________________________________________________________

/*              Function  Output	          Purpose
                --------  ------  ------------------------------------
                  time    double  The length of U in seconds.
                 length   double  The length of U in seconds.
                  dim      int    The Hilbert space dimension of U.
                  HS       int    The Hilbert space dimension of U. 
                  LS       int    The Liouville space dimension of U.
                  H       gen_op  The effective Hamiltonian of U.            */

MSVCDLL double time()   const;
MSVCDLL double length() const;
MSVCDLL int    dim()    const;
MSVCDLL matrix Mx()     const;
MSVCDLL basis  Bs()     const;
MSVCDLL int    HS()     const;
MSVCDLL int    LS()     const;
MSVCDLL gen_op Op()     const;
MSVCDLL gen_op H()      const;

// ____________________________________________________________________________
// C                   PROPAGATOR BASIS FUNCTIONS
// ____________________________________________________________________________

        // Input                U       : A propagator (this)
        //                      Op      : An operator   
        // Output               Op1     : U put into current basis of Op

MSVCDLL void SetEBR() const;
MSVCDLL void SetBasis(const gen_op& Op);
 
// ____________________________________________________________________________
// D                   PROPAGATOR EVOLUTION FUNCTIONS
// ____________________________________________________________________________

        // Input                U       : A propagator (this)
        //                      Op      : An operator 
        // Output               Op1     : Op evolved under prop U

MSVCDLL gen_op evolve(const gen_op& Op) const;

// ____________________________________________________________________________
// E           PROPAGATOR FUNCTIONS, PROPAGATOR WITH PROPAGATOR
// ____________________________________________________________________________

/*  Operator  Output	                     Purpose
    --------  ------  ---------------------------------------------------------
       *        U     U = U1*U2.  The length of U is the sum of U1 & U2 lengths
       *=      void   U = U *U2.  The length of U is the sum of U1 & U2 lengths
       &=      void   U = U2*U.   The length of U is the sum of U1 & U2 lengths

                            Order matters - U1*U2 != U2*U1                   */

MSVCDLL HSprop operator *  (const HSprop& U) const;
MSVCDLL HSprop &   operator *= (const HSprop& U);
MSVCDLL HSprop &   operator &= (const HSprop& U);

// ____________________________________________________________________________
// F             PROPAGATOR FUNCTIONS, PROPAGATOR WITH OPERATOR
// ____________________________________________________________________________


MSVCDLL HSprop sim_trans(const gen_op& Op);
 
        // Input                U       : HS propagator (this). 
        //                      Op      : General operator 
        // Return               U1	: Applies similarity transform on U 
        //                                by Op                    t
        //                                       U1 = Op * U * [Op] 
 

MSVCDLL void sim_trans_ip(const gen_op& Op);

        // Input                U       : HS propagator (this).
        //                      Op      : General operator
        // Return               None    : Applies similarity transform on U
        //                                by Op                    t
        //                                        U = Op * U * [Op]


MSVCDLL HSprop Pow(int n) const;

        // Input                U       : HS propagator (this).
        //                      n       : power
        // Return               None    : Returns U n taken to the nth power


// ____________________________________________________________________________
// G                         PROPAGATOR FUNCTIONS
// ____________________________________________________________________________

///F_list prop  - Retrieve the propagator for the input hamiltonian
 
/* These are from before the addition of this propagator class to GAMMA. They
   will one day become deprecated.  For now the are just "friend" functions
   of this class and gleaned from what used to be nmrlib.cc                  */

	// Input 		ham   : "hamiltonian" for propagation (in Hertz)
        //                      time  : evolution time (seconds)
        // Output               U     : propagator for the input hamiltonian
        // Note                       : For the "ip" function, the evolution
	//                              is done "in place" overwriting ham

#if defined(__SUNPRO_CC) 
MSVCDLL friend gen_op prop(const gen_op& ham, const double time);
MSVCDLL friend void   prop_ip(   gen_op& ham, const double time);
#endif

// ----------------------------------------------------------------------------
//                 EVOLUTION UNDER A STATIC HAMILTONIAN
// ---------------------------------------------------------------------------
 
        // Input                sig   : Op to be propagated (dens. mx.)
        //                      ham   : Op for propagation (in Hertz)
        //                              (usually a static Hamiltonian)
        //                      time  : Evolution time (seconds)
        // Output               sigma1: Sigma evolved by ham for time
        //                              EITHER to new or in place (ip)
        // Note                       : Sigma is set to the EB of U & ham
        //                              inside sim_trans. This is also the
        //                              basis in which the computation is done
        // Note                       : As propagator U is unitless & the exp
        //                              argument is in radians, 2*PI is used


#if defined(__SUNPRO_CC)
MSVCDLL friend gen_op evolve(const gen_op& sigma, const gen_op& ham, double time);
MSVCDLL friend void   evolve_ip(   gen_op& sigma, const gen_op& ham, double time);
#endif

// ----------------------------------------------------------------------------
//                EVOLUTION UNDER A HILBERT SPACE PROPAGATOR
// ----------------------------------------------------------------------------
 
        // Input                sigma : Op to be propagated (unitless)
        //                      U     : Propagator (unitless)
        // Output               sigma1: Op propagated by U to either make a
        //                              new density operator or in place (ip)
        // Note                       : If U is in its EB, sigma set to U EB
        //                              inside sim_trans. This is then the
        //                              basis in which the computation is done
        //                              If U not in its EB it is put into the
        //                              EB of sigma and the computation done
        //                              in that eigenbase

#if defined(__SUNPRO_CC)
MSVCDLL friend gen_op evolve(const gen_op& sigma, const gen_op& U);
MSVCDLL friend void   evolve_ip(   gen_op& sigma, const gen_op& U);
#endif

// ____________________________________________________________________________
// H              HILBERT SPACE PROPAGATOR CONTAINER SUPPORT FUNCTIONS
// ____________________________________________________________________________
 
/* Aside for providing basic tests as to whether to operators are equvalent or
   not, these operators are necessary if any STL container classes are to
   be used based on propagators (e.g. list<HSprop> or vector<HSprop>)        */  

MSVCDLL bool operator== (const HSprop& U) const;
MSVCDLL bool operator!= (const HSprop& U) const;
MSVCDLL bool operator<  (const HSprop& U) const;
MSVCDLL bool operator>  (const HSprop& U) const;
 
// ____________________________________________________________________________
// Z               HILBERT SPACE PROPAGATOR OUTPUT FUNCTIONS
// ____________________________________________________________________________
 
        // Input                U       : A HS propagator
        //                      ostr    : An output stream 
	//			full	: Flag for output amount
        // Output               ostr    : Output stream that has had
        //                                U written into it 

MSVCDLL std::ostream& print(std::ostream& ostr, int full=0) const;
MSVCDLL friend std::ostream &operator << (std::ostream &ostr, const HSprop& U);

};
 
#endif							// HSprop.h
