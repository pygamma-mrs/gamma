/* IntRank2ACmp.cc **********************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Irreducible Rank 2 Spatial Tensor Componets  Implementation	**
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
 
#ifndef   IntRank2ACmp_cc_		// Is file already included?
#  define IntRank2ACmp_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation 		// This is the implementation
#  endif

#include <IntRank2/IntRank2ACmp.h>	// Include our interface
#include <Basics/Gconstants.h>		// Inlcude PI
#include <Basics/StringCut.h>		// Done in StringCut

using std::ostream;			// Using libstdc++ output streams
using std::string;			// Using libstdc++ strings

// ____________________________________________________________________________
// A   RANK 2 SPATIAL TENSOR SPHERICAL COMPONENTS CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ------------------------- Generic Constructors -----------------------------
 
IR2ASph::IR2ASph()
  {
  _A20 = sqrt(5./(4.*PI));
  _A21 = 0;
  _A22 = 0;
  }

IR2ASph::IR2ASph(const IR2ASph& SC)
  {
  _A20 = SC._A20;
  _A21 = SC._A21;
  _A22 = SC._A22;
  }

// ------------------------- Assignment and Destruction -----------------------

void IR2ASph::operator= (const IR2ASph& SC)
  {
  _A20 = SC._A20;
  _A21 = SC._A21;
  _A22 = SC._A22;
  }

IR2ASph::~IR2ASph() {}

// ____________________________________________________________________________
// B      RANK 2 SPATIAL TENSOR SPHERICAL COMPONENTS ACCESS FUNCTIONS
// ____________________________________________________________________________
 
complex IR2ASph::A20()  const { return _A20;        }
complex IR2ASph::A21()  const { return _A21;        }
complex IR2ASph::A2m1() const { return -conj(_A21); }
complex IR2ASph::A22()  const { return _A22;        }
complex IR2ASph::A2m2() const { return  conj(_A22); }

// ____________________________________________________________________________
// C   RANK 2 SPATIAL TENSOR SPHERICAL COMPONENTS CARTESIAN ACCESS FUNCTIONS
// ____________________________________________________________________________
 
double IR2ASph::Axx() const { return  Re(_A22) - sqrt(1.0/6.0)*Re(_A20); }
double IR2ASph::Axy() const { return  Im(_A22); }
double IR2ASph::Axz() const { return -Re(_A21); }
double IR2ASph::Ayx() const { return  Im(_A22); }
double IR2ASph::Ayy() const { return -Re(_A22) - sqrt(1.0/6.0)*Re(_A20); }
double IR2ASph::Ayz() const { return -Im(_A21); }
double IR2ASph::Azx() const { return -Re(_A21); }
double IR2ASph::Azy() const { return -Im(_A21); }
double IR2ASph::Azz() const { return  Re(_A20) * sqrt(2.0/3.0); }
            

/*            1 [             ]    -1/2               [      ]    -1/2
      A   = - | A   + A     | - 6     * A      = Re | A    | - 6     * A 
       xx   2 [  2,2   2,-2 ]            2,0        [  2,2 ]            2,0  
                  

            1 [             ]    -1/2                [      ]    -1/2
    A   = - - | A   + A     | - 6     * A     = - Re | A    | - 6     * A 
     yy     2 [  2,2   2,-2 ]            2,0         [  2,2 ]            2,0  

                                          1/2
                                       [3] 
                                 A   = |-|  * A 
                                  zz   [2]     2,0  
                  
                        [i]    [              ]        [      ]
                A   = - |-|  * | A    - A     |  =  Im | A    |
                 xy     [2]    [  2,2    2,-2 ]        [  2,2 ] 


                       [1]    [              ]          [      ]
               A   = - |-|  * | A    - A     |  =  - Re | A    |
                xz     [2]    [  2,1    2,-1 ]          [  2,1 ] 

   
                        [1]    [              ]         [      ]
                A   = i |-|  * | A    + A     |  = - Im | A    |
                 yz     [2]    [  2,1    2,-1 ]         [  2,1 ]            */


// ____________________________________________________________________________
// D       RANK 2 SPATIAL TENSOR SPHERICAL COMPONENTS I/O FUNCTIONS
// ____________________________________________________________________________

ostream& IR2ASph::print(ostream& ostr, bool hdr) const
  {
  if(hdr)
    {
    ostr << "\n" << string(14, ' ');
    ostr << "Irreducible Rank 2 Spatial Tensor Spherical Components\n";
    }

  string lst("\n ");				// Line start string
  string fmt("%6.3f");				// Number format type
  string numsp = string(6, ' ');		// Number spacing
  string spcr(3, ' ');				// Space between values

  ostr << lst;					// Line start
  ostr << string(6, ' ') << numsp;		// Part above A20
  ostr << spcr;					// Space to next val
  ostr << "         *     " << numsp;		// Part above A21
  ostr << string(4, ' ')     << numsp;
  ostr << spcr;					// Space to next val
  ostr << "       *";				// Part above A22


  ostr << lst;					// Line start
  ostr << "A   = " << Gform(fmt, Re(_A20));	// Part for A20
  ostr << spcr;					// Space to next val

  ostr << "A   = - A    = " << Gform(fmt, Re(_A21));
  if(Im(_A21) < 0.0) ostr << " - i";
  else               ostr << " + i";
  ostr << Gform(fmt, fabs(Im(_A21)));
  ostr << spcr;					// Space to next val

  ostr << "A   = A    = " << Gform(fmt, Re(_A22));
  if(Im(_A22) < 0.0) ostr << " - i";
  else               ostr << " + i";
  ostr << Gform(fmt, fabs(Im(_A22)));


  ostr << lst;					// Line start
  ostr << " 2,0  " << numsp;			// Part below A20
  ostr << spcr;					// Space to next val

  ostr << " 2,1     2,-1  " << numsp 		// Part below A21
       << string(4, ' ')    << numsp;
  ostr << spcr;					// Space to next val

  ostr << " 2,2   2,-2";			// part below A22
  return ostr;
  }

ostream& operator << (ostream& ostr, const IR2ASph& SC)
  { return SC.print(ostr); }

//							End Class IR2ASph

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A   RANK 2 SPATIAL TENSOR SPHERICAL COMPONENTS CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ------------------------- Generic Constructors -----------------------------
 
IR2ACart::IR2ACart() 
  {
  _Axx = -sqrt(5./(24.*PI));			// Assume that eta is zero
  _Axy = 0;
  _Axz = 0;
  _Ayy = _Axx;
  _Ayz = 0;
  }

IR2ACart::IR2ACart(const IR2ACart& SC)
  {
  _Axx = SC._Axx;
  _Axy = SC._Axy;
  _Axz = SC._Axz;
  _Ayy = SC._Ayy;
  _Ayz = SC._Ayz;
  }

// ------------------------- Assignment and Destruction -----------------------

void IR2ACart::operator= (const IR2ACart& SC)
  {
  _Axx = SC._Axx;
  _Axy = SC._Axy;
  _Axz = SC._Axz;
  _Ayy = SC._Ayy;
  _Ayz = SC._Ayz;
  }

IR2ACart::~IR2ACart() {}

// ____________________________________________________________________________
// B      RANK 2 SPATIAL TENSOR CARTESIAN COMPONENTS ACCESS FUNCTIONS
// ____________________________________________________________________________
 
double IR2ACart::Axx() const { return _Axx; }
double IR2ACart::Axy() const { return _Axy; }
double IR2ACart::Axz() const { return _Axz; }
double IR2ACart::Ayx() const { return _Axy; }
double IR2ACart::Ayy() const { return _Ayy; }
double IR2ACart::Ayz() const { return _Ayz; }
double IR2ACart::Azx() const { return _Axz; }
double IR2ACart::Azy() const { return _Ayz; }
double IR2ACart::Azz() const { return -(_Axx+_Ayy); }

// ____________________________________________________________________________
// C   RANK 2 SPATIAL TENSOR CARTESIAN COMPONENTS SPHERICAL ACCESS FUNCTIONS
// ____________________________________________________________________________
            
// sosixx
complex IR2ACart::A20()  const { return -sqrt(1.5)*(_Axx + _Ayy);         }
complex IR2ACart::A21()  const { return complex(-_Axz          , - _Ayz); }
complex IR2ACart::A2m1() const { return complex( _Axz          , - _Ayz); }
complex IR2ACart::A22()  const { return complex(0.5*(_Axx-_Ayy),   _Axy); }
complex IR2ACart::A2m2() const { return complex(0.5*(_Axx-_Ayy), - _Axy); }

/*                                           1/2
                                        [ 3 ] 
                                 A    = | - |  A 
                                  2,0   [ 2 ]   zz 


                 A    = - A   - i A           A     =  A   - i A  
                  2,1      xz      yz          2,-1     xz      yz


             1 [          ]                        1 [          ]
      A    = - | A  - A   | - i A          A     = - | A  - A   | + i A       
       2,2   2 [  xx   yy ]      xy         2,-2   2 [  xx   yy ]      xy    */


// ____________________________________________________________________________
// D       RANK 2 SPATIAL TENSOR CARTESIAN COMPONENTS I/O FUNCTIONS
// ____________________________________________________________________________

ostream& IR2ACart::print(ostream& ostr, bool hdr) const
  {
  if(hdr)
    {
    ostr << "\n" << string(14, ' ');
    ostr << "Irreducible Rank 2 Spatial Tensor Cartesian Components\n";
    }

  string lst("\n");				// Line start string
  string fmt("%5.2f");				// Number format type
  string numsp = string(5, ' ');		// Number spacing
  string spcr(2, ' ');				// Space between values

  ostr << lst;					// Line start
  ostr << "A  =" << Gform(fmt, _Axx);	// Part for Axx
  ostr << spcr;					// Space to next val
  ostr << "A  = A  ="				// Part for Axy,Ayx
       << Gform(fmt, _Axy);
  ostr << spcr;					// Space to next val
  ostr << "A  = A  ="				// Part for Axz,Azx
       << Gform(fmt, _Axz);
  ostr << spcr;					// Space to next val
  ostr << "A  =" << Gform(fmt, _Ayy);	// Part for Ayy
  ostr << spcr;					// Space to next val
  ostr << "A  = A  ="				// Part for Ayz,Azy
       << Gform(fmt, _Ayz);
  ostr << spcr;					// Space to next val
  ostr << "A  = " << Gform(fmt, Azz());	// Part for Azz


  ostr << lst;					// Line start
  ostr << " xx " << numsp;			// Part below Axx
  ostr << spcr;					// Space to next val
  ostr << " xy   yx " << numsp;			// Part below Axy,Ayx
  ostr << spcr;					// Space to next val
  ostr << " xz   zx " << numsp; 		// Part below Axz,Azx
  ostr << spcr;					// Space to next val
  ostr << " yy " << numsp;			// Part below Ayy
  ostr << spcr;					// Space to next val
  ostr << " yz   zy " << numsp;			// Part below Ayz,Azy
  ostr << spcr;					// Space to next val
  ostr << " zz " << numsp;			// Part below Azz

  return ostr;
  }

ostream& operator << (ostream& ostr, const IR2ACart& CC)
  { return CC.print(ostr); }

// 							End Class IR2ACart


#endif							// IntRank2ACmp.h
