/* HSProp.cc ****************************************************-*-c++-*-
**									**
**	                        G A M M A				**
**                                                                      **
**	Hilbert Space Propagator 		Implementation		**
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

#ifndef   HSprop_cc_			// Is file already included?
#  define HSprop_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// this is the implementation
#  endif

#include <HSLib/HSprop.h>		// Include the header
#include <HSLib/GenOp.h> 		// Must know operators
#include <Basics/Gconstants.h>		// Must know PI and PI2
#include <Basics/Gutils.h>		// Know about GAMMA errors
#include <Basics/StringCut.h>	 // Know GAMMA form Gform
#include <string> 			// Must know about strings

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------ 
// ---------------------------------------------------------------------------- 
 
// ____________________________________________________________________________ 
// i              CLASS HILBERT SPACE PROPAGATOR ERROR HANDLING
// ____________________________________________________________________________
  
/*       Input                U       : HS propagator (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

void HSprop::HSPerror(int eidx, int noret) const
  {
  std::string hdr("Hilbert Space Propagator");
  switch(eidx)
    {
    case 0:
    default: GAMMAerror(hdr, eidx, noret); break;
    }
  }

void volatile HSprop::HSPfatal(int eidx) const
  {                                                                       
  HSPerror(eidx);
  if(eidx) HSPerror(0);
  GAMMAfatal();					// Clean exit from program
  }
 
// ____________________________________________________________________________
// ii            CLASS HILBERT SPACE PROPAGATOR CEILING MATRIX
// ____________________________________________________________________________

/* This sets the array Hceil which only is useful when attempting to convert
   a propagator back to an effective Hamitonian. Even worse, it is only
   useful in some instances when attempting that.                            */

void HSprop::SetCeiling(const gen_op& H, bool I)
  {
  int hs = HS();			// Get Hilbert space dimension
  if(I)					// If ceiling is identity
    {					// set it to an I matrix
    Hceil = matrix(hs,hs,i_matrix_type);//   Set to I
    return;				//   We are finished
    }
  H.set_EBR();				// Insure H in its eigenbasis
  Hceil = matrix(hs,hs,d_matrix_type);	// Begin with empty array 
  double Hii;				// <i|H|i> element as real
  for(int i=0; i<hs; i++)		// Loop over Hilbert space
    {
    Hii = Re(H.get(i,i));		//   Get <i|H|i>
    Hceil.put(ceil(Ut*Hii), i, i);	//   Form <i|Hceil|i> from H*t
    }
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------ 
// ---------------------------------------------------------------------------- 
  
// ____________________________________________________________________________ 
// A         CLASS HILBERT SPACE PROPAGATOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________ 
  
/* A NULL propagator is evident because it will have no Hilbert space. There
   are propagators (e.g. ideal pulses) that do not evolve the system in time
   as they are infinitely short, so Ut=0 is no indicatation that the propagator
   is zero or identity. An identity propagator is determined IFF the propagator
   matrix is the identity matrix. Again, one must NOT use Ut to determine
   that. This has ramifications in other areas of this class!                */

HSprop::HSprop()       { Ut = 0; }
HSprop::HSprop(int HS) { Hceil=matrix(HS,HS,i_matrix_type); UOp=Hceil; Ut=0; }
HSprop::HSprop(const gen_op& H, double tevol)
  {
  complex z(0,-PIx2*tevol);		//   Exponential factor (Hz->rad) 
  UOp = H.exp(z);			//   Generate propagator
  Ut  = tevol;				//   Set evolution time
  SetCeiling(H);			//   Set H ceiling matrix
  }

HSprop::HSprop(const gen_op& H, double tevol, bool prop)
  {
  if(prop)				//    If H is already a propagator
    {					// ================================
    UOp = H;				//    U is the input H 
    Ut  = tevol;			//    Ut is the input time 
    int hs = H.HS();			//    This is the Hilbert space
    Hceil = matrix(hs,hs,i_matrix_type);//    Cannot regenerate H from this
    }
  else					//      If H is a Hamiltonian 
    {					// ================================
    complex z(0,-PIx2*tevol);		//   Exponential factor (Hz->rad) 
    UOp = H.exp(z);			//   Generate propagator
    Ut  = tevol;			//   Set evolution time
    SetCeiling(H);			//   Set H ceiling matrix
    }
  }
                                                                                
HSprop::HSprop(const HSprop& U) { UOp=U.UOp; Ut=U.Ut; Hceil=U.Hceil; }

HSprop::~HSprop() { }

HSprop& HSprop::operator= (const HSprop& U)
  {
  if(this == &U) return *this;
  UOp   = U.UOp;				// Copy the operator
  Ut    = U.Ut;					// Copy the evolution time
  Hceil = U.Hceil;				// Copy the H ceiling array
  return *this;
  }

// ____________________________________________________________________________
// B               HILBERT SPACE PROPAGATOR ACCESS FUNCTIONS
// ____________________________________________________________________________

/*              Function  Output                  Purpose
                --------  ------  ------------------------------------
                  time    double  The length of U in seconds.
                 length   double  The length of U in seconds.
                  dim      int    The Hilbert space dimension of U.
                  HS       int    The Hilbert space dimension of U.
                  LS       int    The Liouville space dimension of U.
                  H       gen_op  The effective Hamiltonian of U.            */

double HSprop::time()   const { return Ut; }			// Time
double HSprop::length() const { return Ut; }			// Time
int    HSprop::dim()    const { return UOp.dim(); }		// HS dimension
int    HSprop::HS()     const { return UOp.dim(); }		// HS dimension
int    HSprop::LS()     const { return UOp.dim()*UOp.dim(); }	// LS dimension
matrix HSprop::Mx()     const { return UOp.get_matrix(); }	// U matrix
basis  HSprop::Bs()     const { return UOp.get_basis(); }	// U basis
gen_op HSprop::Op()     const { return UOp; }			// Get operator 
gen_op HSprop::H() const					// Effective H
  {
  if(!HS())               return gen_op();		// Null U -> Null H
  if(Mx().test_type(i_matrix_type) == i_matrix_type) 	// U = I  -> 0 H
    return gen_op(matrix(HS(),HS(),complex0), Bs());
  if(Hceil.stored_type() == i_matrix_type)		// Arbitray U -> 
    return log(UOp)/(-PIx2*complexi*Ut);		//   H via direct ln(U)
  else if(Ut)						// Special U ->
    {							//   H by special ln(U)
    UOp.set_EBR();					//   Put U in eigenbase
    int hs = HS();					//   Hibert dimension
    matrix H(hs, hs, d_matrix_type);			//   Fresh H matrix 
    complex Ukk;
    double Xkk, Nkk;
    for(int k=0; k<hs; k++)				//   Loop diagonal
      {
      Ukk = UOp.get(k,k);				//      <k|U|k>
      Xkk = Re(log(Ukk)/(-complexi*PIx2));		//      <k|X|k>
      Nkk = Hceil.getRe(k,k);				//      <k|N|k>
      H.put((Xkk+Nkk)/Ut, k, k);			//      <k|H|k>
      }
    return gen_op(H,Bs());
    }
  return log(UOp)/(-PIx2*complexi*Ut);			// This if all else 
  }							// fails.

// ____________________________________________________________________________
// C                   PROPAGATOR BASIS FUNCTIONS
// ____________________________________________________________________________

void HSprop::SetEBR() const             { UOp.set_EBR(); }
void HSprop::SetBasis(const gen_op& Op) { UOp.Op_base(Op); }

        // Input                U	: A propagator (this)
	//			Op	: An operator
        // Output               Op1	: U put into current basis of Op

// ____________________________________________________________________________
// D                   PROPAGATOR EVOLUTION FUNCTIONS
// ____________________________________________________________________________

gen_op HSprop::evolve(const gen_op& Op) const { return UOp.sim_trans(Op); }

        // Input                U	: A propagator (this)
	//			Op	: An operator
        // Output               Op1	: Op evolved under prop U

// ____________________________________________________________________________
// E           PROPAGATOR FUNCTIONS, PROPAGATOR WITH PROPAGATOR
// ____________________________________________________________________________

/* These functions allow users to take propagator products, thus building up
   extended propagators. The time of the propagator products will be the sum
   of the individual times. One other aspect that is dealt with here is in the
   abiltiy to generate an effective Hamiltonian from the propagator. Unless the
   two propagators being multiplied are identical we loose the ability to
   regnerate a specific Hamiltonian, but must then return one that may have
   different phase factors than expected. U = exp(-i*H*t) will still hold true.

    Operator  Output                         Purpose
    --------  ------  ---------------------------------------------------------
       *        U     U = U1*U2.  The length of U is the sum of U1 & U2 lengths
       *=      void   U = U *U2.  The length of U is the sum of U1 & U2 lengths
       &=      void   U = U2*U.   The length of U is the sum of U1 & U2 lengths

                            Order matters - U1*U2 != U2*U1                   */

HSprop HSprop::operator * (const HSprop& U1) const
  {
  HSprop U(*this);			// Identical propagator to us
  U  *= U1;
  return U; 
  }      

HSprop & HSprop::operator *= (const HSprop& U1) 
{ 
    UOp *= U1.UOp; 
    Ut += U1.Ut; 
    return *this;
} 

HSprop & HSprop::operator &= (const HSprop& U1) 
{ 
    UOp &= U1.UOp; 
    Ut += U1.Ut; 
    return *this;
} 

// ____________________________________________________________________________
// F             PROPAGATOR FUNCTIONS, PROPAGATOR WITH OPERATOR
// ____________________________________________________________________________
  
        // Input                U       : HS propagator (this).  
        //                      Op      : General operator  
        // Return               U1      : Applies similarity transform on U 
        //                                by Op                    t
        //                                       U1 = Op * U * [Op]  
  
HSprop HSprop::sim_trans(const gen_op& Op)
  {
  HSprop U1(*this);
  U1.sim_trans_ip(Op);
  return U1;
  }
 
void HSprop::sim_trans_ip(const gen_op& Op) { UOp.sim_trans_ip(Op); }
 
        // Input                U       : HS propagator (this).
        //                      Op	: General operator
        // Return               None	: Applies similarity transform on U
	//				  by Op                    t
        //                                        U = Op * U * [Op]
 
        // Input                U       : HS propagator (this).
        //                      n       : power
        // Return               None    : Returns U taken to the nth power

HSprop HSprop::Pow(int n) const
  {
  HSprop Upow;					// Our return propagator
  Upow.UOp   = UOp.Pow(n);			// Set Upow to be U^n
  Upow.Ut    = n*Ut;				// Set Upow time as n*Ut
  Upow.Hceil = Hceil;				// Our ceilngs are same
  return Upow;					// Return Upow = U^n
  }

// ______________________________________________________________________
// F                      PROPAGATOR FUNCTIONS
// ______________________________________________________________________

// These are from before the invention of propagators.  The will one day
//   become deprecated.  For now the are just "friend" functions of this
//   class and gleaned from what used to be nmrlib.cc                     


        //                      ham   : "Hamiltonian" for propagation (in Hertz)
        //                      time  : Evolution time (seconds)
        // Output               U     : Propagator for the input Hamiltonian
        //                      ham   : "Hamiltonian" for propagation (in Hertz)
        //                      time  : Evolution time (seconds)
        // Output               U     : Propagator for the input Hamiltonian
        // Note                       : Done in place, overwriting ham
 

gen_op prop(const gen_op& ham, const double time)
  {
  ham.set_EBR();                        // Put ham into its eigenbase
  complex z(0,-PIx2*time);		// Exponential factor (rad->Hz)
  return ham.exp(z);			// Generate evolution propagator
  }
        
void prop_ip(gen_op& U, const double time)
  {
  U.set_EBR();                          // Put ham into its eigenbase
  complex z(0,-PIx2*time);		// Exponential factor (rad -> Hz)
  U = U.exp(z);                         // Generate evolution propagator
  }



// ____________________________________________________________________________
//                          EVOLUTION  FUNCTIONS
// ____________________________________________________________________________
 
// This set of functions were in GAMMA before the addition of this propagator
//   class.  They will one day become deprecated.  For now the are just "friend"
//   functions of this class and gleaned from what used to be nmrlib.cc        

// ----------------------------------------------------------------------------
//                   EVOLUTION UNDER A STATIC HAMILTONIAN
// ----------------------------------------------------------------------------
 
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

gen_op evolve(const gen_op& sigma, const gen_op& ham, double time)
  {
  if(!ham.exists()) return sigma; 	// No ham, no change in sigma
  gen_op U = prop(ham,time);		// Construct HS propagator
  return U.sim_trans(sigma);		// Return sigma evolved by U
  }

void evolve_ip(gen_op& sigma, const gen_op& ham, double time)
  {
  if(!ham.exists()) return;             // Do nothing if no Hamiltonian
  gen_op U = prop(ham,time);		// Construct HS propagator
  sigma.sim_trans_ip(U);                // Evolve density matrix
  return;                               // Return
  }

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

gen_op evolve(const gen_op& sigma, const gen_op& U)
  {
  if(!U.test_EBR()) 
    U.Op_base(sigma);   // Put U in sigma's basis unless EBR
  return U.sim_trans(sigma);            // Return propagated sigma
  }
 
void evolve_ip(gen_op& sigma, const gen_op& U)
  {
  if(!U.test_EBR()) U.Op_base(sigma);   // Put U in sigma's basis unless EBR
  sigma.sim_trans_ip(U);                // Evolve density matrix
  return;                               // Return sigma1
  }


// ____________________________________________________________________________
// H              HILBERT SPACE PROPAGATOR CONTAINER SUPPORT FUNCTIONS
// ____________________________________________________________________________
 
// Aside for providing basic tests as to whether to operators are equvalent or
//   not, these operators are necessary if any STL container classes are to
//   be used based on propagators (e.g. list<HSprop> or vector<HSprop>)       


bool HSprop::operator== (const HSprop& U) const
  {
  if(U.time() != Ut) return false;	// Retern FALSE if lengths differ
  return(U.Op() == UOp);		// Operators must match
  }

bool HSprop::operator!= (const HSprop& U) const { return !(*this==U); };
bool HSprop::operator<  (const HSprop& U) const { return (dim()<U.dim()); };
bool HSprop::operator>  (const HSprop& U) const { return (dim()<U.dim()); };

// ____________________________________________________________________________
// X                      COMPLEX OPERATOR FUNCTIONS
// ____________________________________________________________________________



// ____________________________________________________________________________
// Z               HILBERT SPACE PROPAGATOR OUTPUT FUNCTIONS
// ____________________________________________________________________________

        // Input                U       : A propagator (this)
	//			ostr	: An output stream
	//			full	: Flag for amount of output
        // Output               ostr	: Output stream that has had
        //              		  U written into it

std::ostream& HSprop::print(std::ostream& ostr, int full) const
  {
                       ostr << "\nEvolution Time: ";
  if(Ut > 1.e-1)       ostr << Gform("%8.3f", Ut)      << " sec\n";
  else if(Ut > 1.e-4)  ostr << Gform("%8.3f", Ut*1.e3) << " msec\n";
  else if(Ut > 1.e-7)  ostr << Gform("%8.3f", Ut*1.e6) << " usec\n";
  else                 ostr << Gform("%8.3f", Ut*1.e9) << " nsec\n";
  if(full >= 0)        UOp.print(ostr, full);
  return ostr;
  }

std::ostream &operator << (std::ostream &ostr, const HSprop &U)
  { U.print(ostr); return ostr; }

 
#endif							// HSprop.cc
