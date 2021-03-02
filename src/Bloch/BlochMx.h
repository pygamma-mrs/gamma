/* BlochMx.h ****************************************************-*-c++-*-
**									**
**	                        G A M M A				**
**								 	**
**	Bloch Matrix 				Interface		**
**						 			**
**      Copyright (c) 2002                                              **
**      Dr. Scott A. Smith						**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**						 			**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**								 	**
** Description							 	**
**								 	**
** The class BlochMx defines a Bloch matrix. Such arrays are meant to   **
** evolve an Bloch magnetization vector in time as described by the	**
** phenomenological Bloch equations. A magnetization vector is embodied **
** in class MagVec in the GAMMA Bloch module. As such vectors involve   **
** N individual sub-vectors, and each of these sub-vectors has three	**
** magnetization components, Bloch matrices are similarly blocked in	**
** 3x3 sub-blocks. Cross terms, i.e. elements outside of the blocks,    **
** unless there is exchange between vectors (chemical, diffusion,...).	**
**									**
** The magnetization vector has ONLY 1 type of evolution implicitly	**
** defined, that dictated by the phenomenological Bloch equations.	**
**									**
**			     -Gt					**
**                 |M(t)> = e	|M(0)-M	  > + |M   >			**
**                                     inf      inf			**	
**									**
** where the evolution matrix G is the Bloch matrix defined herein and	**
** the initinite time vector M** must be specified appropriately. Bloch **
** matrices are simply GAMMA matrices, but they include added 		**
** functionality that allows them to interact in Bloch specific fashion	**
** with other GAMMA classes & maintain appropriate interal structures.	** 
**									**
*************************************************************************/

#ifndef   BlochMxx_h_			// Is file already included?
#  define BlochMxx_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <string>			// Include libstdc++ strings
#include <Matrix/matrix.h>		// Include GAMMA matrices

class BlochMx : public matrix
  {

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                    Bloch Matrix Error Handling
// ____________________________________________________________________________

        // Input                G	: Bloch matrix (this)
        //                      ei	: Error index
        //                      nr	: Flag for linefeed (0=linefeed)
        //                      pn	: string in message
     
         void BMxerror(int ei,                        int nr=0) const;
         void BMxerror(int ei, const std::string& pn, int nr=0) const;
volatile void BMxfatal(int ei)                                  const;
volatile void BMxfatal(int ei, const std::string& pn)           const;

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

// ____________________________________________________________________________
// A                   SPIN SYSTEM CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

///Center Bloch Matrix Algebraic
///F_list BlochMx		- Constructor
///F_list ~			- Destructor
///F_list =			- Assignment	

MSVCDLC          BlochMx();
MSVCDLC          BlochMx(const BlochMx& G);
MSVCDLC          BlochMx(const matrix&  G);
MSVCDLC          ~BlochMx();
MSVCDLC BlochMx& operator=(const BlochMx& G);
MSVCDLC BlochMx& operator=(const matrix& G );

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

MSVCDLL BlochMx operator-  ()                  const;
MSVCDLL BlochMx operator+  (const BlochMx& G1) const;
MSVCDLL BlochMx operator-  (const BlochMx &G1) const;
MSVCDLL BlochMx operator*  (const BlochMx &G1) const;
MSVCDLL void    operator+= (const BlochMx& G1);
MSVCDLL void    operator-= (const BlochMx& G1);
MSVCDLL void    operator*= (const BlochMx& G1);
MSVCDLL void    operator&= (const BlochMx& G1);

// ____________________________________________________________________________
// C                   BLOCH MATRIX - SCALAR INTERACTIONS
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

MSVCDLL BlochMx operator*  (const complex& z) const;
MSVCDLL BlochMx operator*  (double r)         const;
MSVCDLL friend  BlochMx  operator*  (const complex& z,  const BlochMx &G1);
MSVCDLL friend  BlochMx  operator*  (double r,          const BlochMx &G1);
MSVCDLL void    operator*= (const complex& z);
MSVCDLL void    operator*= (double r); 

MSVCDLL BlochMx operator/  (const complex& z) const;
MSVCDLL BlochMx operator/  (double r)         const;
MSVCDLL void    operator/= (const complex& z);
MSVCDLL void    operator/= (double r);

// ____________________________________________________________________________
// I				Basics Python Code
// ____________________________________________________________________________

MSVCDLL int NComps() const;

};

#endif							// BlochMx.h
