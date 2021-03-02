/* BlochProp.cc *************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Bloch Propagators			    Implementation 	**
**                                                                      **
**      Copyright (c) 2002                                              **
**      Dr. Scott A. Smith						**
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
** The class BlochProp defines a Bloch propagator. These are meant to   **
** evolve an Bloch magnetization vector in time as described by the	**
** phenomenological Bloch equations. A magnetization vector is embodied **
** in class MagVec in the UAMMA Bloch module. As such vectors involve   **
** N individual sub-vectors, and each of these sub-vectors has three	**
** magnetization components, Bloch matrices are similarly blocked in	**
** 3x3 sub-blocks. Cross terms, i.e. elements outside of the blocks,    **
** unless there is exchange between vectors (chemical, diffusion,...).	**
**									**
** The magnetization vector has ONLY 1 type of evolution implicitly	**
** defined, that dictated by the phenomenological Bloch equations.	**
**									**
**               -Ut							**
**     |M(t)> = e   |M(0)-M   > + |M   >  = U |M(0)-M   > + |M   >	**
**                         inf      inf	             inf      inf	**
**									**
** where the evolution matrix U is the Bloch propagator and U the Bloch	**
** propagator defined herein. The the initinite time vector M** must be	**
** specified appropriately. Bloch propagators are UAMMA matrices, but	**
** they include added functionality and added restrictions that allow	**
** them to interact in Bloch specific fashion with other UAMMA classes	**
** in support of Bloch simulations. This helps to maintain appropriate	**
** interal structures.							** 
**									**
*************************************************************************/

#ifndef   BlochU_cc_			// Is file already included?
#  define BlochU_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Bloch/BlochProp.h>		// Includes the interface 
#include <Basics/Uutils.h>              // Include UAMMA errors/queries
#include <Basics/StringCut.h>		// Include Udec function

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                      Bloch Propagator Error Handling
// ____________________________________________________________________________

        // Input                U       : Bloch matrix (this)
        //                      eidx    : Error index
        //                      noret   : Flag for linefeed (0=linefeed)
	// Output		void	: An error message is output

/* The following error messages use the defaults set in the Uutils package

                Case                          Error Message

		(0) 			Program Aborting.....
    		(3) 			Can't Construct From Parameter Set
		(4)			Cannot Construct From Input File
    		(5)			Cannot Write To Parameter File
    		(6)			Cannot Write To Output FileStream
    		(7)			Cannot Produce Any Output
    		default			Unknown Error                        */

        // Input                sys     : Spin system (this)
        //                      eidx    : Error index
	//			pname	: Additional error message
        //                      noret   : Flag for linefeed (0=linefeed)
        // Output               none    : Error message output

/* The following error messages use the defaults set in the Uutils package

                Case                          Error Message

		(1) 			Problems With File PNAME
		(2)			Cannot Read Parameter PNAME
    		default			Unknown Error - PNAME                */

        // Input                U       : Bloch matrix (this)
        //                      eidx    : Error index
	//			pname	: Additional error message
        // Output               none    : Error message output
        //                                Program execution stopped

void BlochProp::BUerror(int eidx, int noret) const
  {
  std::string hdr("Bloch Propagator");
  std::string msg;
  switch (eidx)
    {
    case 10: msg = std::string("Rectangular Array Discovered");		// (10)
             UAMMAerror(hdr,msg,noret);  break;
    case 11: msg = std::string("These Arrays Must Be Square");		// (11)
             UAMMAerror(hdr,msg,noret);  break;
    case 12: msg = std::string("Array Dimensioning Problem");		// (12)
             UAMMAerror(hdr,msg,noret);  break;
    case 13: msg = std::string("Array Not 3N x 3N with N=0,1,2,...");	// (13)
             UAMMAerror(hdr,msg,noret);  break;
    case 14: msg = std::string("Problems During Assigment From Matrix");// (14)
             UAMMAerror(hdr,msg,noret);  break;
    default: UAMMAerror(hdr, eidx, noret); break;			// (-1)
    }
  }

void BlochProp::BUerror(int eidx, const std::string& pname, int noret) const
  {
  std::string hdr("Bloch Propagator");
  std::string msg;
  switch(eidx)
    {
    case 101:                                                         // (101)
      msg = std::string("Can't Find Parameters For ")
          + pname;
     UAMMAerror(hdr, msg, noret); break;
    case 102:                                                         // (102)
      msg = std::string("Can't Find ") + pname 
          + std::string(" In Parameter Set");
     UAMMAerror(hdr, msg, noret); break;
    default: UAMMAerror(hdr, eidx, pname, noret); break;// Unknown Error   (-1)
    }
  }

volatile void BlochProp::BUfatal(int eidx) const
  {  
  BUerror(eidx, 1);				// Normal non-fatal error
  if(eidx) BUerror(0);				// Program aborting error
  UAMMAfatal();					// Clean exit from program
  }

volatile void BlochProp::BUfatal(int eidx, const std::string& pname) const
  {  
  BUerror(eidx, pname, 1);			// Normal non-fatal error
  if(eidx) BUerror(0);				// Program aborting error
  UAMMAfatal();					// Clean exit from program
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A               BLOCH PROPAGATOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

BlochProp::BlochProp()                   : matrix()   { Ut = 0;}
BlochProp::BlochProp(const BlochProp& U) : matrix(U)  { Ut = U.Ut; }
BlochProp::~BlochProp() {}

BlochProp& BlochProp::operator=(const BlochProp& U)
  {
  if(this == &U) return *this;			// Nothing if self assign
  matrix::operator=((matrix)U);			// Copy array
  Ut = U.Ut;					// Copy evolution time
  return *this;
  }

// ____________________________________________________________________________
// B               BLOCH PROPAGATOR - BLOCH PROPAGATOR INTERACTIONS
// ____________________________________________________________________________

/* These functions allow for simple mathematical operations between two Bloch
   matices.  This includes addition, subtraction, multiplication. There is
   one unary function as well, negation.
 
   Operator Arguments      Result        Operator Arguments    Result
   -------- --------- -----------------  -------- --------- -------------------   
      -        U      Returns -U            -=      U,U1     U1 subt. from U
      +      U,U1     Returns U+U1          *       U1,U2    Returns U1*U2 
     +=      U,U1     U1 added to U         *=      U,U1     U mult into U1
      -      U1,U2    Returns U1-U2                                          */

BlochProp BlochProp::operator* (const BlochProp& U1) const
  { BlochProp U(*this); U *= U1; return U;}

void BlochProp::operator*= (const BlochProp& U1)
  { matrix::operator*= (U1); Ut += U1.Ut; }

void BlochProp::operator&= (const BlochProp& U1)
  { BlochProp U(U1); U*=(*this); *this = U; }

// ____________________________________________________________________________
// B                   BLOCH PROPAGATOR - SCALAR INTERACTIONS
// ____________________________________________________________________________

/* These functions allow for two mathematical operations between a scalar &
   a Bloch matrix.

 Operator Arguments      Result        Operator Arguments    Result
 -------- --------- -----------------  -------- --------- -------------------
    *        z,U    Returns z*U          *=       U,z     U multiplied by z
    *        U,z    Returns z*U          *=       U,d     U multiplied by d
    *        d,U    Returns d*U          /        U,z     Returns (1/z)*U
    *        U,d    Returns U*d          /        U,d     Returns (1/d)*U
    /=       U,d    U mult. by (1/d)     /=       U,z     U mult. by (1/z) */


BlochProp BlochProp::operator* (const complex& z)
  { BlochProp U(*this); U *= z; return U;}

BlochProp BlochProp::operator* (double r)
  { BlochProp U(*this); U *= r; return U;}

BlochProp operator* (const complex& z, const BlochProp &U1)
  { BlochProp U(U1); U *= z; return U;}

BlochProp operator* (double r, const BlochProp &U1)
  { BlochProp U(U1); U *= r; return U;}

void BlochProp::operator*= (const complex& z)
  { matrix::operator*= (z); }

void BlochProp::operator*= (double r)  
  { matrix::operator*= (r); }

BlochProp BlochProp::operator/ (const complex& z)  { return (*this)*(1/z); }
BlochProp BlochProp::operator/ (double r)          { return (*this)*(1/r); }
void    BlochProp::operator/= (const complex& z) { (*this) *= (1/z);     }
void    BlochProp::operator/= (double r)         { (*this) *= (1/r);     }

// ____________________________________________________________________________
// C          BLOCH PROPAGATOR - MAUNETIZATION VECTOR INTERACTIONS
// ____________________________________________________________________________

/*

BlochProp exp(const BlochProp& U, double t)
  {
  BlochProp U(U1);			// Make a copy of input U
  int hs = U.HS();			// This is full spin Hilbert space
  matrix mxd, mxev;			// Used in matrix diagonalizaion
  U.blow_up();			// Uenerate full space array
  diag(U.mx, mxd, mxev);		// Diagonalize it (into mxd & mxev)
  for(int i=0; i<hs; i++)		// Exponentiate the diagonals
    mxd.put(exp(mxd.get(i,i)),i,i);
  U.mx = mxev * mxd * adjoint(mxev);	// Reform in original basis 
  U.DelSubArrays();			// Delete any sub-space arrays
  return U;				// Were all done
  }
*/

// ____________________________________________________________________________
// B                        Bloch Propagator Access
// ____________________________________________________________________________

int BlochProp::NComps() const { return rows()/3; }

#endif								// BlochProp.cc
