/* BlochProp.h **************************************************-*-c++-*-
**									**
**	                        G A M M A				**
**								 	**
**	Bloch Propagators			    Interface		**
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

#ifndef _BlochU_h_			// Is file already included?
#define _BlochU_h_ 1			// If no, then remember it
#ifdef __UNUU__				// Using the UNU compiler?
#    pragma interface			// this is the interface
#endif

#include <string>			// Include libstdc++ strings
#include <Matrix/matrix.h>		// Include UAMMA matrices

class BlochProp : public matrix
  {
  double Ut;				// Evolution time

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                    Bloch Propagator Error Handling
// ____________________________________________________________________________

        // Input                U	: Bloch propagator (this)
        //                      ei	: Error index
        //                      nr	: Flag for linefeed (0=linefeed)
        //                      pn	: string in message
     
         void BlochProp::BUerror(int ei,                        int nr=0) const;
         void BlochProp::BUerror(int ei, const std::string& pn, int nr=0) const;
volatile void BlochProp::BUfatal(int ei)                                  const;
volatile void BlochProp::BUfatal(int ei, const std::string& pn)           const;

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

// ____________________________________________________________________________
// A                   SPIN SYSTEM CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

///Center Bloch Propagator Algebraic
///F_list BlochProp		- Constructor
///F_list ~			- Destructor
///F_list =			- Assignment	

           BlochProp::BlochProp();
           BlochProp::BlochProp(const BlochProp& U);
           BlochProp::~BlochProp();
BlochProp& BlochProp::operator=(const BlochProp& U);

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

BlochProp BlochProp::operator-  ()                    const;
BlochProp BlochProp::operator+  (const BlochProp& U1) const;
BlochProp BlochProp::operator-  (const BlochProp &U1) const;
BlochProp BlochProp::operator*  (const BlochProp &U1) const;
void      BlochProp::operator+= (const BlochProp& U1);
void      BlochProp::operator-= (const BlochProp& U1);
void      BlochProp::operator*= (const BlochProp& U1);
void      BlochProp::operator&= (const BlochProp& U1);

// ____________________________________________________________________________
// C                   BLOCH PROPAGATOR - SCALAR INTERACTIONS
// ____________________________________________________________________________

/* These functions allow for two mathematical operations between a scalar &
   a Bloch propagator.

 Operator Arguments      Result        Operator Arguments    Result
 -------- --------- -----------------  -------- --------- -------------------
    *        z,U    Returns z*U          *=       U,z     U multiplied by z
    *        U,z    Returns z*U          *=       U,d     U multiplied by d
    *        d,U    Returns d*U          /        U,z     Returns (1/z)*U
    *        U,d    Returns U*d          /        U,d     Returns (1/d)*U
    /=       U,d    U mult. by (1/d)     /=       U,z     U mult. by (1/z) */

BlochProp BlochProp::operator*  (const complex& z);
BlochProp BlochProp::operator*  (double r);
friend  BlochProp  operator*  (const complex& z,  const BlochProp &U1);
friend  BlochProp  operator*  (double r,          const BlochProp &U1);
void    BlochProp::operator*= (const complex& z);
void    BlochProp::operator*= (double r); 

BlochProp BlochProp::operator/  (const complex& z);
BlochProp BlochProp::operator/  (double r);
void      BlochProp::operator/= (const complex& z);
void      BlochProp::operator/= (double r);

int BlochProp::NComps() const;

};

#endif							// BlochProp.h
