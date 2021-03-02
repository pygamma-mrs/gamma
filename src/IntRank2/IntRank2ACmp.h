/* IntRank2ACmp.h ***********************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Irreducible Rank 2 Spatial Tensor Componets	Interface	**
**                                                                      **
**      Scott Smith                                                     **
**      Copyright (c) 1996                                              **
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
** Class IntRank2A represents a primitive rank 2 interaction. It is	**
** defined only in the PAS and maintains only the asymmetry value.      **
** However, it does have the ability to undergo rotations to produce    **
** both spherical and Cartesian tensor elements in frames relative to	**
** the PAS. These are only output, not stored in the class itself.      ** 
**                                                                      **
** When the tensor is rotated there are 5 spherical components and 	**
** 9 Cartesian components that are potentially useful. These should be  **
** produced simultaneously when the tensor is rotated.  		**
**                                                                      **
** Two classes are herein defined. The first handles storage and use of **
** the irreducible rank 2 spatial tensor spherical components.  The 	**
** second handles the corresponding Cartesian components.		**
**                                                                      **
*************************************************************************/
 
#ifndef   IntRank2ACmp_h_		// Is file already included?
#  define IntRank2ACmp_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface 			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/complex.h>		// Include GAMMA complex numbers
#include <iostream>			// Include libstdc++ file streams

//extern const double RT6PIO5; 		// Needed constant sqrt[6*PI/5]
//extern const double RT5O4PI;		// Needed constant sqrt[5/(4*PI)];
//extern const double RT5O24PI;		// Needed constant sqrt[5/(24*PI)];


//#include <Basics/ParamSet.h>		// Include parameter sets
//#include <Matrix/matrix.h>		// Include GAMMA matrices
//#include <Matrix/row_vector.h>		// Include GAMMA row vectors
//#include <Level1/coord.h>		// Include coordinates
//#include <string>			// Include libstdc++ strings
//#include <fstream>			// Include libstdc++ file streams
//#include <vector>			// Include libstdc++ STL vectors

class IR2ASph
  {
  complex _A20;				// Spherical A20 component
  complex _A21;				// Spherical A21 component
  complex _A22;				// Spherical A22 component

  friend class IntRank2A;		// Allow IntRank2A full access
   

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A   RANK 2 SPATIAL TENSOR SPHERICAL COMPONENTS CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ------------------------- Generic Constructors -----------------------------
 
MSVCDLC IR2ASph();
MSVCDLC IR2ASph(const IR2ASph& SC);

// ------------------------- Assignment and Destruction -----------------------

MSVCDLL void operator= (const IR2ASph &IR2Ab);
MSVCDLC      ~IR2ASph();

// ____________________________________________________________________________
// B      RANK 2 SPATIAL TENSOR SPHERICAL COMPONENTS ACCESS FUNCTIONS
// ____________________________________________________________________________
 
MSVCDLL complex A20()  const;
MSVCDLL complex A21()  const;
MSVCDLL complex A2m1() const;
MSVCDLL complex A22()  const;
MSVCDLL complex A2m2() const;

// ____________________________________________________________________________
// C   RANK 2 SPATIAL TENSOR SPHERICAL COMPONENTS CARTESIAN ACCESS FUNCTIONS
// ____________________________________________________________________________
 
MSVCDLL double Axx() const;
MSVCDLL double Axy() const;
MSVCDLL double Axz() const;
MSVCDLL double Ayx() const;
MSVCDLL double Ayy() const;
MSVCDLL double Ayz() const;
MSVCDLL double Azx() const;
MSVCDLL double Azy() const;
MSVCDLL double Azz() const;
            
/*                            1 [             ]    -1/2
                    A  (EA) = - | A   + A     | - 6     * A 
                     xx       2 [  2,2   2,-2 ]            2,0  
                  

                               1 [             ]    -1/2
                   A  (EA) = - - | A   + A     | - 6     * A 
                    yy         2 [  2,2   2,-2 ]            2,0  
                  
                                          1/2
                                       [3] 
                             A  (EA) = |-|  * A 
                              zz       [2]     2,0  
                  
                                [i]    [                      ]
                    A  (EA) = - |-|  * | A   (EA) - A    (EA) |
                     xy         [2]    [  2,2        2,-2     ] 


                                [1]    [                      ]
                    A  (EA) = - |-|  * | A   (EA) - A    (EA) |
                     xz         [2]    [  2,1        2,-1     ] 

   
                                [1]    [                      ]
                    A  (EA) = i |-|  * | A   (EA) + A    (EA) |
                     yz         [2]    [  2,1        2,-1     ]              */


// ____________________________________________________________________________
// D       RANK 2 SPATIAL TENSOR SPHERICAL COMPONENTS I/O FUNCTIONS
// ____________________________________________________________________________

MSVCDLL std::ostream& print(std::ostream& ostr, bool hdr=true) const;
MSVCDLL friend std::ostream& operator << (std::ostream& ostr, const IR2ASph& SC);

  };							// End Class IR2ASph

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

class IR2ACart
  {
  double _Axx;				// Cartesian Axx component
  double _Axy;				// Cartesian Axy component
  double _Axz;				// Cartesian Axz component
  double _Ayy;				// Cartesian Ayy component
  double _Ayz;				// Cartesian Ayz component

  friend class IntRank2A;		// Allow IntRank2A full access
   

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A   RANK 2 SPATIAL TENSOR SPHERICAL COMPONENTS CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ------------------------- Generic Constructors -----------------------------
 
MSVCDLC IR2ACart();
MSVCDLC IR2ACart(const IR2ACart& SC);

// ------------------------- Assignment and Destruction -----------------------

MSVCDLL void operator= (const IR2ACart &IR2Ab);
MSVCDLC ~IR2ACart();

// ____________________________________________________________________________
// B      RANK 2 SPATIAL TENSOR CARTESIAN COMPONENTS ACCESS FUNCTIONS
// ____________________________________________________________________________
 
MSVCDLL double Axx() const;
MSVCDLL double Axy() const;
MSVCDLL double Axz() const;
MSVCDLL double Ayx() const;
MSVCDLL double Ayy() const;
MSVCDLL double Ayz() const;
MSVCDLL double Azx() const;
MSVCDLL double Azy() const;
MSVCDLL double Azz() const;

// ____________________________________________________________________________
// C   RANK 2 SPATIAL TENSOR CARTESIAN COMPONENTS SPHERICAL ACCESS FUNCTIONS
// ____________________________________________________________________________
            
MSVCDLL complex A20()  const;
MSVCDLL complex A21()  const;
MSVCDLL complex A2m1() const;
MSVCDLL complex A22()  const;
MSVCDLL complex A2m2() const;

/*                                           1/2
                                        [ 3 ] 
                                 A    = | - |  A 
                                  2,0   [ 2 ]   zz 


                 A    = - A   + i A           A     =  A   - i A  
                  2,1      xz      yz          2,-1     xz      yz


             1 [          ]                        1 [          ]
      A    = - | A  - A   | - i A          A     = - | A  - A   | + i A       
       2,2   2 [  xx   yy ]      xy         2,-2   2 [  xx   yy ]      xy    */

// ____________________________________________________________________________
// D       RANK 2 SPATIAL TENSOR CARTESIAN COMPONENTS I/O FUNCTIONS
// ____________________________________________________________________________

MSVCDLL std::ostream& print(std::ostream& ostr, bool hdr=true) const;
MSVCDLL friend std::ostream& operator << (std::ostream& ostr, const IR2ACart& CC);



  };							// End Class IR2ACart


#endif							// IntRank2ACmp.h
