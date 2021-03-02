/* BlochMx.cc ****************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Bloch Matrix 				Implementation 		**
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
**                                                                      **
** Description                                                          **
**                                                                      **
** The class BlochMx defines a Bloch matrix. Such arrays are used to	**
** evolve an Bloch magnetization vector in time as described by the     **
** phenomenological Bloch equations. A magnetization vector is embodied **
** in class MagVec in the GAMMA Bloch module. As such vectors involve   **
** N individual sub-vectors, and each of these sub-vectors has three    **
** magnetization components, Bloch matrices are similarly blocked in    **
** 3x3 sub-blocks. Cross terms, i.e. elements outside of the blocks,    **
** unless there is exchange between vectors (chemical, diffusion,...).  **
**                                                                      **
** The magnetization vector has ONLY 1 type of evolution implicitly     **
** defined, that dictated by the phenomenological Bloch equations.      **
**                                                                      **
**                           -Gt                                        **
**                 |M(t)> = e   |M(0)- M   > + |M   >			**
**                                      inf      inf			**
**                                                                      **
** where the evolution matrix G is the Bloch matrix defined herein and  **
** the initinite time vector must be specified appropriately as needed.	**
** Bloch matrices are simply GAMMA matrices, but they are restricted	**
** to be square and of dimension 3N x 2N where N is the number of sub-	**
** vectors (magnetization). Bloch matrices include added functionality	**
** above class matrix that allows them to interact in a Bloch specific	**
** manner with other GAMMA classes while maintaining the appropriate	**
** interal structures.							**
**                                                                      **
*************************************************************************/

#ifndef   BlochMxx_cc_			// Is file already included?
#  define BlochMxx_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Bloch/BlochMx.h>		// Includes the interface 
#include <Basics/Gutils.h>              // Include GAMMA errors/queries
#include <Basics/StringCut.h>		// Include Gdec function

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                      Bloch Matrix Error Handling
// ____________________________________________________________________________

        // Input                G       : Bloch matrix (this)
        //                      eidx    : Error index
        //                      noret   : Flag for linefeed (0=linefeed)
	// Output		void	: An error message is output

/* The following error messages use the defaults set in the Gutils package

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

/* The following error messages use the defaults set in the Gutils package

                Case                          Error Message

		(1) 			Problems With File PNAME
		(2)			Cannot Read Parameter PNAME
    		default			Unknown Error - PNAME                */

        // Input                G       : Bloch matrix (this)
        //                      eidx    : Error index
	//			pname	: Additional error message
        // Output               none    : Error message output
        //                                Program execution stopped

void BlochMx::BMxerror(int eidx, int noret) const
  {
  std::string hdr("Bloch Matrix");
  std::string msg;
  switch (eidx)
    {
    case 10: msg = std::string("Rectangular Array Discovered");		// (10)
             GAMMAerror(hdr,msg,noret);  break;
    case 11: msg = std::string("These Arrays Must Be Square");		// (11)
             GAMMAerror(hdr,msg,noret);  break;
    case 12: msg = std::string("Array Dimensioning Problem");		// (12)
             GAMMAerror(hdr,msg,noret);  break;
    case 13: msg = std::string("Array Not 3N x 3N with N=0,1,2,...");	// (13)
             GAMMAerror(hdr,msg,noret);  break;
    case 14: msg = std::string("Problems During Assigment From Matrix");// (14)
             GAMMAerror(hdr,msg,noret);  break;
    default: GAMMAerror(hdr, eidx, noret); break;			// (-1)
    }
  }

void BlochMx::BMxerror(int eidx, const std::string& pname, int noret) const
  {
  std::string hdr("Bloch Matrix");
  std::string msg;
  switch(eidx)
    {
    case 101:                                                         // (101)
      msg = std::string("Can't Find Parameters For ")
          + pname;
     GAMMAerror(hdr, msg, noret); break;
    case 102:                                                         // (102)
      msg = std::string("Can't Find ") + pname 
          + std::string(" In Parameter Set");
     GAMMAerror(hdr, msg, noret); break;
    default: GAMMAerror(hdr, eidx, pname, noret); break;// Unknown Error   (-1)
    }
  }

volatile void BlochMx::BMxfatal(int eidx) const
  {  
  BMxerror(eidx, 1);				// Normal non-fatal error
  if(eidx) BMxerror(0);				// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

volatile void BlochMx::BMxfatal(int eidx, const std::string& pname) const
  {  
  BMxerror(eidx, pname, 1);			// Normal non-fatal error
  if(eidx) BMxerror(0);				// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A               BLOCH MATRIX CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

BlochMx::BlochMx()                 : matrix()   {}
BlochMx::BlochMx(const BlochMx& G) : matrix(G)  {}
BlochMx::BlochMx(const matrix&  G) : matrix(G)
  {
  if(!is_square())
    {
    BMxerror(10,1);				// Rectangular array
    BMxerror(11,1);				// Must be square
    BMxfatal(9);				// Error during construction
    }
  double ddim = rows();
  while (ddim > 0) ddim -= 3;
  if(ddim)
    {
    BMxerror(12,1);				// Array dimension problem
    BMxerror(12,1);				// Must be 3Nx3N, N non-neg int
    BMxfatal(9);				// Error during construction
    }
  }

BlochMx::~BlochMx() {}				// Matrix does the work

BlochMx& BlochMx::operator=(const BlochMx& G)
  {
  if(this == &G) return *this;			// Nothing if self assign
  matrix::operator=((matrix)G);			// Copy array
  return *this;
  }

BlochMx& BlochMx::operator=(const matrix& G)
  {
  matrix::operator=((matrix)G);			// Copy array
  if(!is_square())
    {
    BMxerror(10,1);				// Rectangular array
    BMxerror(11,1);				// Must be square
    BMxfatal(14);				// Error during assignment
    }
  double ddim = rows();
  while (ddim > 0) ddim -= 3;
  if(ddim)
    {
    BMxerror(12,1);				// Array dimension problem
    BMxerror(12,1);				// Must be 3Nx3N, N non-neg int
    BMxfatal(14);				// Error during assignment
    }
  return *this;
  }

// ____________________________________________________________________________
// B               BLOCH MATRIX - BLOCH MATRIX INTERACTIONS
// ____________________________________________________________________________

/* These functions allow for simple mathematical operations between two Bloch
   matices.  This includes addition, subtraction, multiplication. There is
   one unary function as well, negation.
 
   Operator Arguments      Result        Operator Arguments    Result
   -------- --------- -----------------  -------- --------- -------------------   
      -        G      Returns -G            -=      G,G1     G1 subt. from G
      +      G,G1     Returns G+G1          *       G1,G2    Returns G1*G2 
     +=      G,G1     G1 added to G         *=      G,G1     G mult into G1
      -      G1,G2    Returns G1-G2                                          */



BlochMx BlochMx::operator- () const { return BlochMx(matrix::operator- ()); }

BlochMx BlochMx::operator+ (const BlochMx& G1) const
  { BlochMx G(*this); G += G1; return G;}

void BlochMx::operator+= (const BlochMx& G1)
  { matrix::operator+= (G1); }

BlochMx BlochMx::operator- (const BlochMx &G1) const
  { BlochMx G(*this); G -= G1; return G;}

void BlochMx::operator-= (const BlochMx& G1)
  { matrix::operator-= (G1); }

BlochMx BlochMx::operator* (const BlochMx& G1) const
  { BlochMx G(*this); G *= G1; return G;}

void BlochMx::operator*= (const BlochMx& G1)
  { matrix::operator*= (G1); }

void BlochMx::operator&= (const BlochMx& G1)
  { BlochMx G(G1); G*=(*this); *this = G; }

// ____________________________________________________________________________
// B                   BLOCH MATRIX - SCALAR INTERACTIONS
// ____________________________________________________________________________

/* These functions allow for two mathematical operations between a scalar &
   a Bloch matrix.

 Operator Arguments      Result        Operator Arguments    Result
 -------- --------- -----------------  -------- --------- -------------------
    *        z,G    Returns z*G          *=       G,z     G multiplied by z
    *        G,z    Returns z*G          *=       G,d     G multiplied by d
    *        d,G    Returns d*G          /        G,z     Returns (1/z)*G
    *        G,d    Returns G*d          /        G,d     Returns (1/d)*G
    /=       G,d    G mult. by (1/d)     /=       G,z     G mult. by (1/z) */


BlochMx BlochMx::operator* (const complex& z) const
  { BlochMx G(*this); G *= z; return G;}

BlochMx BlochMx::operator* (double r) const
  { BlochMx G(*this); G *= r; return G;}

BlochMx operator* (const complex& z, const BlochMx &G1)
  { BlochMx G(G1); G *= z; return G;}

BlochMx operator* (double r, const BlochMx &G1)
  { BlochMx G(G1); G *= r; return G;}

void BlochMx::operator*= (const complex& z)
  { matrix::operator*= (z); }

void BlochMx::operator*= (double r)  
  { matrix::operator*= (r); }

BlochMx BlochMx::operator/ (const complex& z) const { return (*this)*(1/z); }
BlochMx BlochMx::operator/ (double r)         const { return (*this)*(1/r); }
void    BlochMx::operator/= (const complex& z)      { (*this) *= (1/z);     }
void    BlochMx::operator/= (double r)              { (*this) *= (1/r);     }

// ____________________________________________________________________________
// C          BLOCH MATRIX - MAGNETIZATION VECTOR INTERACTIONS
// ____________________________________________________________________________

/*

BlochMx exp(const BlochMx& G, double t)
  {
  BlochMx G(G1);			// Make a copy of input G
  int hs = G.HS();			// This is full spin Hilbert space
  matrix mxd, mxev;			// Used in matrix diagonalizaion
  G.blow_up();			// Generate full space array
  diag(G.mx, mxd, mxev);		// Diagonalize it (into mxd & mxev)
  for(int i=0; i<hs; i++)		// Exponentiate the diagonals
    mxd.put(exp(mxd.get(i,i)),i,i);
  G.mx = mxev * mxd * adjoint(mxev);	// Reform in original basis 
  G.DelSubArrays();			// Delete any sub-space arrays
  return G;				// Were all done
  }
*/

// ____________________________________________________________________________
// B                        Bloch Matrix Access
// ____________________________________________________________________________

int BlochMx::NComps() const { return rows()/3; }

#endif							// BlochMx.cc
