/* LSprop.cc ****************************************************-*-c++-*-
**									**
**	                        G A M M A				**
**                                                                      **
**	Liouville Space Propagator 		Implementation		**
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
** The class LSprop defines a propagator in spin Liouville space.  Such	**
** propagators are simply Liouville space superoperators which will	**
** evolve a density operator for a specific length of time in a 	**
** particular set of rotating frames.					**
**									**
*************************************************************************/

#ifndef   LSprop_cc_			// Is file already included?
#  define LSprop_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// this is the implementation
#  endif

# include <LSLib/LSprop.h>		// Include the header
# include <LSLib/SuperOp.h>		// Must know superoperators
# include <LSLib/DensOp.h>		// Know about density operators
# include <HSLib/GenOp.h>		// Know about operators
# include <Basics/Gutils.h>		// Know about GAMMA errors
# include <string>			// Must know about strings

//sosi - this needs to be fixed too 
//# include <WBR/relaxProp.h>		// Need R_prop function

using std::string;			// Using libstdc++ strings
using std::ostream;			// Using libstdc++ output streams

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------ 
// ---------------------------------------------------------------------------- 
 
// ____________________________________________________________________________ 
// i              CLASS LIOUVILLE SPACE PROPAGATOR ERROR HANDLING
// ____________________________________________________________________________
  
        // Input                G       : LS propagator (this)
        //                      eidx    : Error flag    
        //                      noret   : Return flag 
        // Output               none    : Error Message Output 
                                                               
void LSprop::LSPerror(int eidx, int noret) const
  {
  string hdr("Liouville Space Propagator");
  switch(eidx)
    {
    case 0: GAMMAerror(hdr, 0, noret); break;   // Program Aborting        (0)
    case 9: GAMMAerror(hdr, 9, noret); break;   // Construction Problems   (9)
    default:GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }
 
void volatile LSprop::LSPfatal(int eidx) const
  {                                                                       
  LSPerror(eidx, 1);
  if(eidx) LSPerror(0);
  GAMMAfatal();					// Clean exit from program
  }
 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------ 
// ---------------------------------------------------------------------------- 
  
// ____________________________________________________________________________ 
// A         CLASS LIOUVILLE SPACE PROPAGATOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________ 

/* A NULL propagator is evident because it will have no Liouville space. There
   are propagators (e.g. ideal pulses) that do not evolve the system in time
   as they are infinitely short, so Ut=0 is no indicatation that the propagator
   is zero or identity. An identity propagator is determined IFF the propagator
   matrix is the identity matrix. Again, one must NOT use Ut to determine
   that. This has ramifications in other areas of this class!                */
  
// ---------------------------------------------------------------------------- 
//                            Basic Constructors
// ---------------------------------------------------------------------------- 

LSprop::LSprop() { Gt = 0; }
LSprop::LSprop(int LS) { GOp = matrix(LS,LS,i_matrix_type); Gt = 0; }
  
// ---------------------------------------------------------------------------- 
//                          Hilbert Space Constructors
// ---------------------------------------------------------------------------- 

/* These propagators will be the Liouville space equivalents to their
   Hilbert space counterparts. Density operator evolution occurs according
   to
              -iHt            iHt                t
             e     * sigma * e    = U * sigma * U    ==> G |sigma> 

                    t
   where G = U (x) U 


                                H	: Active Hamiltonian (Hz)
	  			U	: Active propagator (1/sec)
                                tevol	: Evolution time (seconds)
           Output               LSprop	: Propagator for the input
                                          Hamiltonian or propagator         */

LSprop::LSprop(const gen_op& H, double tevol)
  {
  HSprop U(H, tevol);			// Hilbert space propagator
  gen_op UOp = U.Op();			// Get the propagator operator
  GOp =	U_transform(UOp);		// We are UOp (x) UOp*
  Gt  = tevol;				// Set our evolution time
  }

LSprop::LSprop(const gen_op& H, double tevol, bool prop)
  {
  HSprop U(H, tevol, prop);		// Hilbert space propagator
  gen_op UOp = U.Op();			// Get the propagator operator
  GOp =	U_transform(UOp);		// We are UOp (x) UOp*
  Gt  = tevol;				// Set our evolution time
  }

LSprop::LSprop(const HSprop& U)
  {
  gen_op UOp = U.Op();			// Get the propagator operator
  GOp =	U_transform(UOp);		// We are UOp (x) UOp*
  Gt  = U.time();			// Set our evolution time
  }

// ---------------------------------------------------------------------------- 
//                          Liouville Space Constructors
// ---------------------------------------------------------------------------- 

/* These propagators may or may not have equivalent Hilbert space counterparts.
   Density operator evolution occurs according either to

                            -iLt                        
                           e     |sigma>  ==> G |sigma> 

   if the Liouvillian acts directly on the density operator or as


   if the inputLouvillian acts on the difference density operator.
   Note that in these functions the Liouvillian should have the bulk of its
   frequency components as imaginary terms and the bulk of its relaxation and
   exchange components are real terms, e.g. L = -i[H, ] + R + X.  See the
   documentation (theory section) for specifics as to how the superoparator
   propagator is formed when the system evolves to a steady state.

                                Leff	: Active Liouvillian (1/sec)
	  			sigmass	: Steady state density operator
                                tevol	: Evolution time (seconds)
           Output               LSprop	: Propagator for input Liouvillian  */

LSprop::LSprop(const super_op& Leff, double tevol)
  {
  GOp =	exp(Leff, -tevol);		// Exponential Liouvillian
  Gt  = tevol;				// Set our evolution time
  }

LSprop::LSprop(const super_op& Leff, const densop& sigmass, double tevol)
  {
  super_op eLt = exp(Leff, -tevol);	// Exponential Liouvillian
GOp = R_prop(eLt, sigmass);
/*
  GOp = eLt;                            // Copy exp. to prop. superoperator
  eLt.LOp_base(sigmass);		// Put sigmass into eLt Hilbert basis
  sigmass.SetTrace(1.0);		// Set sigmass trace to 1
  int hs = sigmass.size();		// Get the Hilbert space size
  int a,aa,b,bb,g,gg,aaa=0,bbb=0,ggg=0;	// HS indicies for spanning LS
  complex Rel;				// Temp matrix element
  int nd = sigmass.in_EBR();		// See if sigmass is diagonal (nd=1)
  for(b=0; b<hs; b++)			// Fill up GOp elements
    {
    for(bb=0; bb<hs; bb++)
      {
      if(b==bb)				// When |bb>!=|b> then
        { 				// <aaa|GOp|bbb>=<aaa|eLt|bbb>
        aaa = 0; 			// so do'nt make any modifications
        for(a=0; a<hs; a++)
          {
          for(aa=0; aa<hs; aa++)
            {
            Rel = GOp.get(aaa,bbb);
            if(!nd || a==aa)		// If non-diagonal sigmass or |a>=|aa>
              Rel += sigmass.get(a,aa);	// must add <aaa|sigmass|1> to each
            ggg=0; 			// element <aaa|GOp|bbb>
            for(g=0; g<hs; g++)
              {
              for(gg=0; gg<hs; gg++)
                {
                Rel -= eLt.get(aaa,ggg)*sigmass.get(g,gg);
                ggg++;
                }
              }
            GOp.put(aaa,bbb,Rel);
            aaa++;
            }
          }
        }
      bbb++;
      }  
    }  
*/
  Gt = tevol;				// Set prop. evolution time
  }


/* The constructor below is used when it is desired to convert a superoperator
   that is already a propagator into a object of type LSprop. For example, this
   can occur if one creates a unitary transform superoperator (U_transform)
   from a general operator (gen_op) that is itself a Hilbert space propagator.
   zum beispeil: gen_op H; gen_op U = prop(H,t); super_op G = U_transform(U);
   Of course, there not many good reasons (perhaps building ideal pulse props)
   to do this because it breaks the ease of using propagators. Instead use code
   such as: gen_op H; HSprop U(H,t); LSprop(U); In the rare times where it is
   necessary to use the constructor below, make sure you specify the 
   propagator evolution time after construction or it will be left as zero.  */

LSprop::LSprop(const super_op& G)
  {
  GOp =	G;				// Set the propagator
  Gt  = 0;				// Set our evolution time
  }

// ---------------------------------------------------------------------------- 
//                Self Constructor, Assignment, Destructor
// ---------------------------------------------------------------------------- 

        LSprop::LSprop(const LSprop& G) { GOp = G.GOp; Gt  = G.Gt; }
        LSprop::~LSprop() { }
LSprop& LSprop::operator= (const LSprop& G)
  {
  if(this == &G) return *this;		// Nothing if self-assign
  GOp = G.GOp;				// Copy the propagator
  Gt  = G.Gt;				// Copy the evolution time
  return *this;				// Return ourself
  }

// ____________________________________________________________________________
// B               LIOUVILLE SPACE PROPAGATOR ACCESS FUNCTIONS
// ____________________________________________________________________________

/*              Function  Output                  Purpose
                --------  ------  ------------------------------------
                  time    double  The length of G in seconds.
                 length   double  The length of G in seconds.
                  dim      int    The Liouville space dimension of G.
                  HS       int    The Hilbert space dimension of G.
                  LS       int    The Liouville space dimension of G.        */

double   LSprop::time()   const { return Gt; }			// Evolve time
double   LSprop::length() const { return Gt; }			// Evolve time
int      LSprop::dim()    const { return GOp.LS(); }		// LS dimension
int      LSprop::HS()     const { return GOp.HS(); }		// HS dimension
int      LSprop::LS()     const { return GOp.HS(); }		// LS dimension
super_op LSprop::LOp()    const { return GOp; }			// Superoperator

void LSprop::L(const super_op& LOp) { GOp = LOp; }
void LSprop::length(double t)       { Gt  = t;   }


// ____________________________________________________________________________
// C             LIOUVILLE SPACE PROPAGATOR AUXILIARY FUNCTIONS
// ____________________________________________________________________________

void LSprop::SetEBR() const { GOp.set_EBR(); }

// ____________________________________________________________________________
// C                   PROPAGATOR BASIS FUNCTIONS
// ____________________________________________________________________________

void LSprop::SetBasis(const super_op& LOp)
  {
  GOp.LOp_Hbase(LOp);		// GOp into LOp Hilbert basis
  GOp.LOp_base(LOp);		// GOp into LOp Liouville basis
  }

// ____________________________________________________________________________
// C                   PROPAGATOR EVOLUTION FUNCTIONS
// ____________________________________________________________________________

/************** Evolution Under This Superoperator Propagator ****************/

gen_op LSprop::evolve(const gen_op& Op) { return GOp*Op; }

        // Input                G       : A LS propagator (this)
        //                      Op      : An operator
        // Output               Op1     : Op evolved under prop G

/************ Evolution Under A Static Superoperator (Liouvillian) ***********/ 

gen_op evolve(const gen_op &sigma, super_op &LOp, const double time)

	// Input		sigma : Op to be propagated (dens. mx.)
	//	 		LOp   : Propagation superoperator (rad/sec)
	//	 		time  : Evolution time (seconds)
	// Output		sigma1: Sigma evolved by LOp for time

  {
  LOp.set_EBR();			// Set LOp into its eigenbase
  complex z(-time);			// Exponential factor (sec)
  super_op U = LOp.exp(z);		// Evolve superop [exp(sec*rad/sec)]
  return U*sigma;			// Return sigma evolved by U
  }

/************** Evolution Under A Supplied Superoperator Propagator **********/

gen_op evolve(gen_op &sigma, super_op &GOp) { return (GOp*sigma); }

	// Input		sigma : Op. to be propagated (density matrix)
	//	 		GOp   : Super operator propagator
	// Output		Op    : Operator, density matrix propagated
	// Note			      : The superoperator GOp is here assumed
	//				to be a propagation superoperator!



	// Input		sigma : Operator propagated (density matrix)
	//	 		GOp   : Super operator propagator
	// Output		none  : Density matrix, sigma, is modified

void evolve_ip(gen_op &sigma, super_op &GOp)
  {
  GOp.LOp_base(sigma);			// Put sigma basis of GOp
  sigma = GOp*sigma;			// Evolve sigma under GOp
  }


// ____________________________________________________________________________
// C           PROPAGATOR FUNCTIONS, PROPAGATOR WITH PROPAGATOR
// ____________________________________________________________________________

        // Input                G1	: A propagator
        //                      G2	: Another propagator
        // Return               G	: Propagator that is the product of
	//				  the two input propagators, G = G1*G2
        // Note				: Order matters - G1*G2 != G2*G1

        // Input                G       : Super propagator (this).
        //                      G1      : Super propagator.
        // Return               Op      : Propagator which is the input
        //                                propagator multiplied into G1
        //                                        G = G * G1
        // Note                         : Order matters - G1*G2 != G2*G1
 
        // Input                G       : LS propagator (this).
        //                      G1      : LS propagator.
        // Return               Op      : Propagator which is the input
        //                                propagator multiplied by G1
        //                                        G = G1 * G
        // Note                         : Order matters - G1*G2 != G2*G1
 
LSprop LSprop::operator *  (const LSprop& G) const { LSprop G1(*this); G1 *= G; return G1; }
LSprop &  LSprop::operator *= (const LSprop& G)       { GOp *= G.GOp; Gt += G.Gt; return (*this);} 
LSprop &  LSprop::operator &= (const LSprop& G)       { GOp &= G.GOp; Gt += G.Gt; return (*this);} 

// ____________________________________________________________________________
// D           PROPAGATOR FUNCTIONS, PROPAGATOR WITH SUPEROPERATOR
// ____________________________________________________________________________

        // Input                LOp	: A superoperator
        //                      G 	: Another superoperator
        // Return               G1	: Propagator that is the product of
	//				  the superoperator and propagator,
	//					G1 = LOp*G
        // Note				: Order matters - LOp*G != G*LOp
 
LSprop operator * (super_op& LOp, LSprop& G)
  { 
  LSprop G1;
  G1.L(LOp*G.LOp());
  G1.length(G.length());
  return G1;
  }


// ____________________________________________________________________________
// Y                      TEMPORARY PROPERATOR FUNCTIONS
// ____________________________________________________________________________

/* These are taken from the relaxProp module and will either be moved here
   permanantly or placed in an appropriate spot.			     */
 
void set_trace(gen_op& sigma, double tr)
  {
  complex z;
  double acttr = Re(trace(sigma));
  double difftr = acttr - tr;
  double fact;
  int i;
  int hs = sigma.size();
  if(difftr)
    {
    fact = difftr/double(hs);
    for(i=0; i<hs; i++)
      {  
      z = sigma.get(i,i) - fact;
      sigma.put(z, i,i);
      }  
    }  
  acttr = Re(trace(sigma));
  }  


// Input		eLt	: Exponential Liouvillian to relax
//				  the density matrix
//			sigmaeq : Density matrix at equilibrium
//				  (or steady state density matrix)
//			Ho      : Operator in Hilbert space
// Output		LOp     : Superoperator which fullfills
//
//	   eLt * {sigma-sigmaeq} + sigmaeq  = LOp * sigma
//
// Note				: LOp constructed in the current basis of
//				  sigmaeq & in the default Liouville basis
// Note				: sigmaeq is copied so trace is unchanged

super_op R_prop(super_op& eLt, gen_op sigmaeq)
  {
  eLt.set_HBR();				// Put eLt in default Liouville basis
  super_op LOp(eLt);				// Copy the original superoperator
  eLt.LOp_base(sigmaeq);			// Put sigma into eLt Hilbert basis
  set_trace(sigmaeq, 1.0);			// Set the trace to 1
  int hs = sigmaeq.size();			// Get the Hilbert space size
  int a,aa,b,bb,g,gg,aaa=0,bbb=0,ggg=0;
  complex Rel;
  int nd = sigmaeq.in_EBR();			// Test if sigmaeq is diagonal
// *should add a check for i-matrix so prop is i-matrix returned
//  if(mx.stored_type() != i_matrix_type)
// *actually need to test (1-eLt)*sigmaeq and see if it is diagonal, not sigmaeq!
// should set nd=0;

  for(b=0; b<hs; b++)				// For nd=1, diagonal 
    for(bb=0; bb<hs; bb++)
      {
      if(b==bb) 				// <aaa|LOp|bbb> = <aaa|eLt|bbb> if |bb> != |b>
        {					// so do not make any modifications
        aaa = 0;		
        for(a=0; a<hs; a++)
          for(aa=0; aa<hs; aa++)
            {
            Rel = LOp.get(aaa,bbb);
            if(!nd || a==aa)			// If non-diagonal sigmaeq or |a> = |aa> must
              Rel += sigmaeq.get(a,aa);		// add <aaa|sigmaeq|1> to each <aaa|LOp|bbb>
            ggg=0;
            for(g=0; g<hs; g++)
              for(gg=0; gg<hs; gg++)
                {
                Rel -= eLt.get(aaa,ggg)*sigmaeq.get(g,gg);
                ggg++;
                }
            LOp.put(aaa,bbb,Rel);
            aaa++;
            }
        }
      bbb++;
      }
  return LOp;
  }



// Input		eLt	: Liouvillian for evolution
//			sigmaeq : Density matrix at equilibrium
//				  (or steady state density matrix)
//			t       : Evolution time
// Output		LOp     : Superoperator which fullfills
//
//	   eLt * {sigma-sigmaeq} + sigmaeq  = LOp * sigma

super_op R_prop(super_op& L, gen_op& sigmaeq, double t)
  {
  super_op LOp = exp(L, -t);		// Exponential Liouvillian
  return R_prop(LOp, sigmaeq);		// Use function overload
  }

// ____________________________________________________________________________
// Z               LIOUVILLE SPACE PROPAGATOR OUTPUT FUNCTIONS
// ____________________________________________________________________________

        // Input                G       : A propagator (this)
	//			ostr	: An output stream
	//			full	: Flag for amount of output
        // Output               ostr	: Output stream that has had
        //              		  G written into it

ostream& LSprop::print(ostream& ostr, int full) const
  {
  ostr << "\nEvolution Time: " << Gt << " s\n";
  if(full >= 0) GOp.print(ostr, full);
  return ostr;
  }

ostream &operator << (ostream &ostr, const LSprop &G)
  {
  G.print(ostr);
  return ostr;
  }

#endif						// LSprop.cc
