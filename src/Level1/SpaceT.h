/* SpaceT.h *****************************************************-*-c++-*-
**									**
** 	                       G A M M A				**
**									**
**	Spatial Tensors                                 Interface	**
**								 	**
**	Copyright (c) 1991, 1992				 	**
**	Scott Smith				       			**
**	Eidgenoessische Technische Hochschule			 	**
**	Labor fur physikalische Chemie				 	**
**	8092 Zurich / Switzerland				 	**
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
** Description:								**
**									**
** This file contains the definition and workings of general spatial	**
** tensors. The class tracks the tensor via its irreducible components.	**
**									**
*************************************************************************/

///Chapter Class Spatial Tensor (space_T)
///Section Overview
///Body None.
///Section Available Spatial Tensor Functions

#ifndef   Gspace_T_h_			// Is this file already included?
#  define Gspace_T_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/row_vector.h>		// Include knowledge of row vectors
#include <Level1/coord.h>		// Include knowledge of coordinates
#include <Basics/ParamSet.h>		// Include parameter sets


// Forward declaration for the benefit 
// of the following function declarations.
class space_T;

// Added the following function declarations 
// to augment the "friend" declarations inside 
// the class.

MSVCDLL space_T A1(double x, double y, double z);

MSVCDLL space_T A1(coord &pt);

// No longer supported *** Marked for deletion 
MSVCDLL space_T A1(row_vector &vx);

// l=0 no longer supported *** Marked for l=0 case removal
MSVCDLL complex A1(double x, double y, double z, int m, int l=1);

// l=0 no longer supported *** Marked for l=0 case removal
MSVCDLL complex A1(coord &pt, int m, int l=1);

// No longer supported *** Marked for deletion 
MSVCDLL complex A1(row_vector &vx, int m, int l=1);

// l=0 no longer supported *** Marked for deletion
MSVCDLL double A10(int m);

MSVCDLL complex A11(double x, double y, double z, int m);

MSVCDLL space_T SphA1(complex plus, complex zero, complex minus);

MSVCDLL space_T SphA1(coord &pt);

// No longer supported *** Marked for deletion 
MSVCDLL space_T SphA1(col_vector &vx);

MSVCDLL space_T A2(double Aiso, double delzz, double eta=0,
			double alpha=0, double beta=0, double gamma=0);

MSVCDLL space_T A2(coord &Tcomps);

MSVCDLL space_T A2(coord &Tcomps, coord &Tangles);

MSVCDLL space_T A2(const matrix &mx, double prec=1e-10);

MSVCDLL complex A2(int l, int m, double Aiso, double delzz, double eta);

MSVCDLL complex A2(int l, int m, const matrix& mx);

MSVCDLL complex A20(int m, double Aiso, double delzz, double eta);

MSVCDLL complex A20(int m, const matrix &mx);

MSVCDLL complex A21(int m, double Aiso, double delzz, double eta);


MSVCDLL complex A21(int m, const matrix &mx);

MSVCDLL complex A22(int m, double Aiso, double delzz, double eta);

MSVCDLL complex A22(int m, const matrix &mx);

MSVCDLL complex T_comp(const space_T &SphT, int L, int M);
// ?? replaced by member function component

MSVCDLL space_T T_mult(const space_T &SphT1, const space_T &SphT2);


MSVCDLL space_T T_mult(const space_T &SphT1, int l1, const space_T &SphT2, int l2);

MSVCDLL space_T T_mult(const space_T &SphT1, int l1,
					 const space_T &SphT2, int l2, int L);

MSVCDLL complex T_mult(const space_T &SphT1, int l1,
				const space_T &SphT2, int l2, int L, int M);


MSVCDLL space_T T_rot(space_T &SphT1, double alpha, double beta, double gamma);

MSVCDLL void T_rot(int num, space_T* SphT, space_T* SphTrot, double alpha, double beta, double gamma);

MSVCDLL  void T_rot(int num, space_T* SphT, space_T* SphTrot, matrix* D);

// ?? replaced by member function rotate
MSVCDLL complex T_rot(space_T &SphT1, int l, int m,
				 double alpha, double beta, double gamma);







class space_T
  {
  friend class AQuad;

// ----------------------------------------------------------------------------
// ------------------------------ STRUCTURE -----------------------------------
// ----------------------------------------------------------------------------

  int rank;			// Tensor rank
  row_vector **vx;		// Irreducible components
  coord EA;			// Current orientation (Euler angles)
  coord PAS_EA;			// PAS relative orientaton (Euler angles)
  coord PAS;			// PAS irreducible components

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
  
// ____________________________________________________________________________
// i                      SPACE TENSOR ERROR HANDLING
// ____________________________________________________________________________


void SphTerror(int eidx, int noret=0) const;

        // Input                SphT    : Spatial tensor (this)
        //                      eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message
 

void SphTerror(int eidx, const std::string& pname, int noret=0) const;

        // Input                SphT    : Spatial tensor (this)         
        //                      eidx    : Error index
        //                      pname   : Additional error message
        //                      noret   : Flag for linefeed (0=linefeed)
        // Output               none    : Error message output
 

volatile void SphTfatality(int eidx) const;
         
        // Input                SphT    : Spatial tensor (this)
        //                      eidx    : Flag for error type 
        // Output               none    : Error message output 
        //                                Program execution stopped 

volatile void SphTfatality(int eidx, const std::string& pname) const;
         
        // Input                SphT    : Spatial tensor (this)
        //                      eidx    : Flag for error type 
        //                      pname   : Additional error message
        // Output               none    : Error message output 
        //                                Program execution stopped 


friend void space_T_error(int i);

	// Input		i    : Error Flag
	// Output		none : Error Message Output


friend void volatile space_T_fatality(int error);

	// Input		none :
	// Output		none : Stops Execution & Error Message Output
  
// ____________________________________________________________________________
// ii                   SPACE TENSOR PRIVATE AUXILIARY FUNCTIONS
// ____________________________________________________________________________


void updatePAS();

	// Input		SphT : Spatial tensor (this)
	// Output		none : The interal PAS values are set
        // Note                      : Assumes that changes to A   must be
        //                             relfected in PAS         l,m 


public:

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
  
// ____________________________________________________________________________
// A                  SPACE TENSOR CONSTRUCTORS, DESTRUCTOR
// ____________________________________________________________________________


MSVCDLC space_T();

	// Input		None :
	// Return		SphT : Empty spatial tensor (this)
	///F_list space_T	     - Constructor	


MSVCDLC space_T(const space_T& SphT);

	// Input		SphT : Spatial tensor
	// Return		SphT1: Spatial tensor duplicate of input


MSVCDLC space_T(const space_T& SphT, int l);

	// Input		SphT : Spatial tensor
	//			l    : Tensor rank
	// Return		SphT1: Irreducible spatial tensor duplicate
	//			       of input spatial tensor rank l


MSVCDLC space_T(const SinglePar& par);

	// Input		par   : GAMMA parameter
	// Output		Spht  : Spatial Tensor


MSVCDLC virtual ~space_T();
	// Input		SphT : Spatial tensor (this)
	// Return		     : None, deletes current Tensor

// ______________________________________________________________________
// B                    SPACE TENSOR UNARY OPERATIONS
// ______________________________________________________________________


MSVCDLL virtual space_T& operator = (const space_T &SphT);
	// Input		SphT : spatial tensor
	// Return		SphT1: spatial tensor equivalent to the
        //	                       input spatial tensor
	// Note		             : SphT and SphT1 must be associated
	//			       with the same space system
	///F_list =		     - assignment


MSVCDLL friend space_T operator + (const space_T& SphT1, const space_T& SphT2);
	// Input		SphT1 : spatial tensor
	// 			SphT2 : spatial tensor
	// Return		SphT3 : spatial tensor which is the addition
	//			        of the two input spatial tensors
	//		                SphT3 = SphT1 + SphT2	
	// Note		              :	Tensor rank will be the higher
	//				of the input tensor ranks
	///F_list +		      - addition

MSVCDLL friend void operator += (space_T& SphT1, const space_T& SphT2);
	// Input		SphT1 : spatial tensor
	// 			SphT2 : spatial tensor
	// Return		SphT1 : spatial tensor which is the addition
	//				of the two input spatial tensors
	//		               	SphT1 = SphT1 + SphT2	
	// Note		              :	Tensor rank will be the higher
	//				of the input tensor ranks
	///F_list +=		      - unary addition


MSVCDLL friend space_T operator - (const space_T& SphT1, const space_T& SphT2);
	// Input		SphT1 : spatial tensor
	// 			SphT2 : spatial tensor
	// Return		SphT3 : spatial tensor which is the addition
	//			        of the two input spatial tensors
	//		                SphT3 = SphT1 + SphT2	
	// Note		              :	Tensor rank will be the higher
	//				of the input tensor ranks
	///F_list -		      - subtraction

MSVCDLL friend void operator -= (space_T& SphT1, const space_T& SphT2);
	// Input		SphT1 : spatial tensor
	// 			SphT2 : spatial tensor
	// Return		SphT1 : spatial tensor which is the addition
	//				of the two input spatial tensors
	//		               	SphT1 = SphT1 + SphT2	
	// Note		              :	Tensor rank will be the higher
	//				of the input tensor ranks
	///F_list -=		      - unary subtraction

MSVCDLL friend space_T operator - (const space_T &SphT);
	// Input		SphT  : spatial tensor
	// Return		SphT1 : spatial tensor which is the negated
	//			        input spatial tensor
	//		                SphT3 = - SphT2	
	///F_list -		      - negation

// ______________________________________________________________________
// C                SPACE TENSOR WITH SCALARS
// ______________________________________________________________________

///Center Tensor Functions with Scalars


MSVCDLL friend space_T operator * (const complex &z, const space_T &SphT);
	// Input		z     : complex number
	// 			SphT  : spatial tensor
	// Return		SphT1 : spatial tensor which is the input
	//			        spatial tensor multiplied by the
	//				complex scalar constant
	//		                SphT1 = z * SphT
	///F_list *		      - scalar multiplication

MSVCDLL friend space_T operator * (const space_T& SphT, const complex& z);
	// Input		SphT  : spatial tensor
	// 			z     : complex number
	// Return		SphT1 : spatial tensor which is the input
	//			        spatial tensor multiplied by the
	//				complex scalar constant
	//		                SphT1 = SphT * z


// ______________________________________________________________________
// D                  PRINCIPLE AXIS SYSTEM ACCESS
// ______________________________________________________________________

///Center Tensor Principle Axis System Functions

MSVCDLL coord PASys() const;

	// Input		SphT  : Spherical tensor (this)
	// Output		PAS   : Principal Axis System Components
	// Note			      : Only for rank 2 tensors
	// Note			      : Function name doesn't conflict with
	//				structure name
	///F_list PASys		      - PAS components (rank 2)

MSVCDLL coord PASys_EA() const;

	// Input		SphT  : Spherical tensor (this)
	// Output		pt    : Euler Angle set for orientation
	//				of the principle Axes
	// Note			      : Function name doesn't conflict with
	//				structure name
	///F_list PASys_EA	      - PAS Euler Angles (rank 2)

MSVCDLL double iso() const;

	// Input		SphT  : Spherical tensor (this)
	// Output		Aiso  : Isotropic component
	// Note			      : Only for rank 2 tensors
	///F_list iso		      - Isotropic component (rank 2)

MSVCDLL void iso(double Aiso);

	// Input		SphT  : Spherical tensor (this)
	// 			Aiso  : Isotropic tensor component
	// Output		none  : SphT isotropic component set
	//				to value Aiso
	// Note			      : Only for rank 2 tensors

MSVCDLL double delz() const;

	// Input		SphT  : Spherical tensor (this)
	// Output		delz  : Delta z component
	// Note			      : Only for rank 2 tensors
	///F_list delz		      - Delta z component (rank 2)

MSVCDLL void delz(double delzz);

	// Input		SphT  : Spherical tensor (this)
	// 			delz  : Delta z component of PAS
	// Output		      : Delta z component of PAS set to delz
	// Note			      : Only for rank 2 tensors

MSVCDLL double eta() const;

	// Input		SphT  : Spherical tensor (this)
	// Output		eta   : eta component
	// Note			      : Only for rank 2 tensors
	///F_list eta		      - Eta component (rank 2)


MSVCDLL void eta(double ETA);

	// Input		SphT  : Spherical tensor (this)
	// 			ETA   : The eta component of PAS
	// Output		      : The eta component of PAS set to ETA
	// Note			      : Only for rank 2 tensors


MSVCDLL double alpha() const;

	// Input		SphT  : Spherical tensor (this)
	// Output		alpha : alpha Euler angle of PAS
	// Note			      : Only for rank 2 tensors
	///F_list iso		      - Alpha Euler Angle (rank 2)


MSVCDLL double beta() const;

	// Input		SphT  : Spherical tensor (this)
	// Output		beta  : beta Euler angle of PAS
	// Note			      : Only for rank 2 tensors
	///F_list iso		      - Beta Euler Angle (rank 2)


MSVCDLL double gamma() const;

	// Input		SphT  : Spherical tensor (this)
	// Output		gamma : gamma Euler angle of PAS
	// Note			      : Only for rank 2 tensors
	///F_list iso		      - Gamma Euler Angle (rank 2)

// ______________________________________________________________________
// E              GENERAL RANK 1 SPHERICAL TENSOR FUNCTIONS
// ______________________________________________________________________

///Center Rank 1 Tensor Functions


MSVCDLL friend space_T A1(double x, double y, double z);

	// Input		x     : Cartesian x component
	// 			y     : Cartesian y component
	// 			z     : Cartesian z component
	// Output		SphT  : Rank 1 Spatial Tensor for
	//			        Cartesian vector specified
	///F_list A1		      - Rank 1 Tensor from Cartesian x,y,z


MSVCDLL friend space_T A1(coord &pt);

	// Input		pt    : Cartesian coordinate
	// Output		SphT  : Rank 1 Spatial Tensor for
	//			        Cartesian coordinate specified


// No longer supported *** Marked for deletion 
MSVCDLL friend space_T A1(row_vector &vx);

	// Input		vx    : Cartesian vector
	// Output		SphT  : Rank 1 Spatial Tensor for
	//			        Cartesian vector specified
 

// l=0 no longer supported *** Marked for l=0 case removal
MSVCDLL friend complex A1(double x, double y, double z, int m, int l);

	// Input		x     : Cartesian x component
	// 			y     : Cartesian y component
	// 			z     : Cartesian z component
	// 			m     : Component index
	// Output		c     : Rank 1 Spatial Tensor component
	//			        Alm for Cartesian vector specified

 
// l=0 no longer supported *** Marked for l=0 case removal
MSVCDLL friend complex A1(coord &pt, int m, int l);

	// Input		pt    : Cartesian coordinate 
	// 			m     : Component index
	// Output		c     : Rank 1 Spatial Tensor component
	//			        Alm for Cartesian coordinate specified


// No longer supported *** Marked for deletion 
MSVCDLL friend complex A1(row_vector &vx, int m, int l);

	// Input		vx    : Cartesian vector
	// 			m     : Component index
	// Output		c     : Rank 1 Spatial Tensor component
	//			        Alm for Cartesian vector specified

 
// l=0 no longer supported *** Marked for deletion
MSVCDLL friend double A10(int m);

	// Input		m     : Component index
	// 			x     : cartesian x component
	// 			y     : cartesian y component
	// 			z     : cartesian z component
	// Output		c     : Rank 1 Spatial Tensor component
	//			        A00 for Cartesian vector specified


MSVCDLL friend complex A11(double x, double y, double z, int m);

	// Input		x     : cartesian x component
	// 			y     : cartesian y component
	// 			z     : cartesian z component
	// 			m     : Component index
	// Output		c     : Rank 1 Spatial Tensor component
	//			        A1m for Cartesian vector specified


MSVCDLL friend space_T SphA1(complex plus, complex zero, complex minus);

	// 			plus  : Component A11
	// 			zero  : Component A10
	// 			minus : Component A1-1
	// Output		SphT  : Rank 1 Spatial Tensor for
	//			        spherical vector specified
	///F_list SphA1		      - Rank 1 Tensor from Spherical values


MSVCDLL friend space_T SphA1(coord &pt);

	// 			pt    : Coordinate point
	// Output		SphT  : Rank 1 Spatial Tensor for
	//			        spherical vector specified


// No longer supported *** Marked for deletion 
MSVCDLL friend space_T SphA1(col_vector &vx);

	// 			vx    : Spherical vector
	// Output		SphT  : Rank 1 Spatial Tensor for
	//			        spherical vector specified
 
// ______________________________________________________________________
// F              GENERAL RANK 2 SPHERICAL TENSOR FUNCTIONS
// ______________________________________________________________________

///Center Rank 2 Tensor Functions

// ********************* Whole Tensor Production ************************

MSVCDLL friend space_T A2(double Aiso, double delzz, double eta,
			double alpha, double beta, double gamma);

	// Input		Aiso  : Isotropic component
	//			delzz : Anisotropy parameter
	//			eta   : Asymmetry parameter 
	//			alpha : Euler angle (degrees)
	//			beta  : Euler angle (degrees)
	//			gamma : Euler angle (degrees)
	// Output		SphT  : Rank 2 Spatial Tensor
	///F_list A2		      - Rank 2 symmetric tensor


MSVCDLL friend space_T A2(coord &Tcomps);

	// Input		Tcomps : Tensor PAS components
	//				 Aiso, delzz, eta
	// Output		SphT   : Rank 2 Spatial Tensor


MSVCDLL friend space_T A2(coord &Tcomps, coord &Tangles);

	// Input		Tcomps : Tensor PAS components
	//				 Aiso, delzz, eta
	// 			Tangles: Tensor PAS Euler angles
	//				 alpha, beta, gamma (degrees)
	// Output		SphT   : Rank 2 Spatial Tensor


MSVCDLL friend space_T A2(const matrix &mx, double prec);

	// Input		mx    : Cartesian rank two tensor
	// Output		SphT  : rank 2 Spatial Tensor

 
// ************* Specific Alm Tensor Component Production ***************

MSVCDLL friend complex A2(int l, int m, double Aiso, double delzz, double eta);

	// Input		l     : component index
	//			m     : component index
	// 			Aiso  : isoropic component
	//			delzz : anisotropy parameter
	//			eta   : asymmetry parameter 
	// Output		SOp   : rank 2 Spatial Tensor component A
	//								 lm


MSVCDLL friend complex A2(int l, int m, const matrix& mx);

	// Input		l     : component index
	//			m     : component index
	// 			mx    : Cartesian rank two tensor
	// Output		SOp   : rank 2 Spatial Tensor component A
	//								 lm


MSVCDLL friend complex A20(int m, double Aiso, double delzz, double eta);

	// Input		m     : component index
	// 			Aiso  : isoropic component
	//			delzz : anisotropy parameter
	//			eta   : asymmetry parameter 
	// Output		SOp   : rank 2 Spatial Tensor component A
	//								 00


MSVCDLL friend complex A20(int m, const matrix &mx);

	// Input		m     : component index
	// 			mx    : Cartesian rank two tensor
	// Output		SOp   : rank 2 Spatial Tensor component A
	//								 00
 

MSVCDLL friend complex A21(int m, double Aiso, double delzz, double eta);

	// Input		m     : component index
	// 			Aiso  : isoropic component
	//			delzz : anisotropy parameter
	//			eta   : asymmetry parameter 
	// Output		SOp   : rank 2 Spatial Tensor component A
	//								 1m


MSVCDLL friend complex A21(int m, const matrix &mx);

	// Input		m     : component index
	// 			mx    : Cartesian rank two tensor
	// Output		SOp   : rank 2 Spatial Tensor component A
	//								 1m


MSVCDLL friend complex A22(int m, double Aiso, double delzz, double eta);

	// Input		m     : component index
	// 			Aiso  : isoropic component
	//			delzz : anisotropy parameter
	//			eta   : asymmetry parameter 
	// Output		SOp   : rank 2 Spatial Tensor component A
	//								 2m


MSVCDLL friend complex A22(int m, const matrix &mx);

	// Input		m     : component index
	// 			mx    : Cartesian rank two tensor
	// Output		SOp   : rank 2 Spatial Tensor component A
	//								 2m


// ______________________________________________________________________
// G                     SPATIAL TENSOR FUNCTIONS
// ______________________________________________________________________

///Center Tensor Auxiliary Functions

MSVCDLL int exists() const;

	// Input		SphT : Spatial tensor (this)
	// Output		TF   : Returns TRUE is SphT is not NULL
	///F_list exists	     - Existence test
 
 
MSVCDLL int exists(int l) const;
 
        // Input                SphT : Spatial tensor (this)
        //                      l    : Rank
        // Output               TF   : Returns if SphT has an l
        //                             component
 

MSVCDLL int Rank() const;

	// Input		SphT : Spatial tensor (this)
	// Output		l    : Overall tensor rank
	// Note			     : Function name needs upper case
	//			       to distinguish it from structure
	///F_list Rank		     - Tensor rank
 

MSVCDLL friend complex T_comp(const space_T &SphT, int L, int M);
// ?? replaced by member function component


MSVCDLL complex component(int L, int M) const;

	// Input		SphT : Spatial tensor (this)
	// 			L    : Momentum index
	// 			M    : Momentum index
	// Output		z    : Complex number, the L,M
	//			       component of Spatial Tensor SphT
	///F_list component	     - Tensor component
 

MSVCDLL double Ccomponent(int r, int c=0) const;

        // Input                SphT : Spatial Tensor (this)
        //                      r    : row index
        //                      c    : column index
        // Output               Arc  : Cartesian tensor component
 

MSVCDLL friend space_T T_mult(const space_T &SphT1, const space_T &SphT2);

	// Input		SphT1 : Spatial Tensor
	// 			SphT2 : Spatial Tensor
	// Output		SphT  : Spatial Tensor which is the product
	//			        of the two input tensors
	//				SphT = SphT1 x SphT2
	///F_list T_mult	      - Tensor multiplication
 

MSVCDLL friend space_T T_mult(const space_T &SphT1, int l1, const space_T &SphT2, int l2);

	// Input		SphT1 : Spatial Tensor
	//			l1    : Rank Component of Tensor 1
	// 			SphT2 : Spatial Tensor
	//			l2    : Rank Component of Tensor 2
	// Output		SphT  : Spatial Tensor which is the
	//				tensor product of the irreducible
	//			        rank l1 component of SphT1 and the
	//				irreducible l2 component of SphT2
	//				SphT = SphT1   x SphT2
	//				            l1        l2
	// Note			      : SphT will be reducible, components
	//				spanning ranks l1+l2 down to
	//				|l1-l2|
 

MSVCDLL friend space_T T_mult(const space_T &SphT1, int l1,
					 const space_T &SphT2, int l2, int L);

	// Input		SphT1 : Spatial Tensor
	//			l1    : Rank Component of Tensor 1
	// 			SphT2 : Spatial Tensor
	//			l2    : Rank Component of Tensor 2
	//			L     : Rank Component of Product Tensor
	// Output		SphT  : Irreducible Spatial Tensor which is the
	//				tensor product of the irreducible
	//			        rank l1 component of SphT1 and the
	//				irreducible l2 component of SphT2
	//				SphT  = SphT1   x SphT2
	//				    L        l1        l2
	// Note			      : L must be in the range
	//                              [l1+l2, |l1-l2|]

 

MSVCDLL friend complex T_mult(const space_T &SphT1, int l1,
				const space_T &SphT2, int l2, int L, int M);

	// Input		SphT1 : Spatial Tensor
	//			l1    : Rank Component of Tensor 1
	// 			SphT2 : Spatial Tensor
	//			l2    : Rank Component of Tensor 2
	//			L     : Rank of Irreducible Product Tensor
	//			M     : Component of Irreducible Product Tensor
	// Output		z     : complex number which is the M
	//				component of the rank L irreducible
	//				spatial tensor SphT, resulting from
	//				the tensor product of the irreducible
	//			        rank l1 component of SphT1 and the
	//				irreducible l2 component of SphT2
	//				SphT  = SphT1   x SphT2
	//				    L        l1        l2
	// Note			      : L must be in the range
	//                              [l1+l2, |l1-l2|]
 

MSVCDLL friend space_T T_rot(space_T &SphT1, double alpha, double beta, double gamma);

MSVCDLL friend void T_rot(int num, space_T* SphT, space_T* SphTrot, double alpha, double beta, double gamma);

        // Input                num   : Number of tensors
        //                      SphT  : Array of spatial tensors
        //                      SphTr : Array of rotated spatial tensors
        //                      alpha : Euler Angle (degrees)
        //                      beta  : Euler Angle (degrees)
        //                      gamma : Euler Angle (degrees)
        // Output               void  : The array SphTr is filled 
        //                              with rotated tensors in array SphT


MSVCDLL friend  void T_rot(int num, space_T* SphT, space_T* SphTrot, matrix* D);

        // Input                num   : Number of tensors
        //                      SphT  : Array of spatial tensors
        //                      SphTr : Array of rotated spatial tensors
        //                      D     : Array or Wigner rotation matrices
        // Output               void  : The array SphTr is filled
        //                              with rotated tensors in array SphT


// ?? replaced by member function rotate
MSVCDLL friend complex T_rot(space_T &SphT1, int l, int m,
				 double alpha, double beta, double gamma);
// ?? replaced by member function rotate


MSVCDLL space_T rotate(double alpha, double beta, double gamma) const;

	// Input		SphT1 : Space tensor (this)
	// 			alpha : Euler angle
	// 			beta  : Euler angle
	// 			gamma : Euler angle
	// Output		SphT  : Space tensor which is the
	//			        the input tensor in the coordinate
	//				system rotated by input Euler angles
	///F_list rotate	      - Tensor rotation
 

MSVCDLL space_T rotate(const coord &EA) const;

	// Input		SphT1 : Spatial tensor (this)
	// 			EA    : Set of Euler angles
	// Output		SphT  : Spatial Tensor which is the product
	//			        the input tensor in the coordinate
	//				system rotated by input Euler angles
 

MSVCDLL complex rotate(int l, int m, double alpha, double beta, double gamma) const;

	// Input		SphT1 : Space tensor (this)
	// 			l     : Component index
	// 			m     : Component index
	// 			alpha : Euler angle
	// 			beta  : Euler angle
	// 			gamma : Euler angle
	// Output		z     : Space tensor component of the
	//			        input tensor in the coordinate
	//				system rotated by input Euler angles
 

MSVCDLL complex rotate(int l, int m, const coord &EA) const;

	// Input		SphT1 : Spatial tensor (this)
	//			l     : Component index
	//			m     : Component index
	// 			EA    : Set of Euler angles
	// Output		Alm   : Spatial Tensor component of
	//			        the input tensor in the coordinate
	//				system rotated by input Euler angles

 
// ______________________________________________________________________
// H                      PARAMETER SET FUNCTIONS
// ______________________________________________________________________
///Center Tensor Parameter Set Functions

 
MSVCDLL SinglePar param(const std::string& pname) const;

        // Input               SphT  : A spatial tensor (this)
        //                     pname : A parameter name
        // Return              par   : A GAMMA parameter of type space_T 
        //                             with the name pname

 
MSVCDLL SinglePar param(const std::string& pname, const std::string& pstate) const;

        // Input               SphT  : A spatial tensor (this)
        //                     pname : A parameter name
        //                     pstate: A parameter statement
        // Return              par   : A GAMMA parameter of type space_T 
        //                             with the name pname and statement pstate


MSVCDLL operator ParameterSet();

	// Input		SphT  : Spatial tensor
	//  			pset  : Parameter set
	// Output		pset  : Parameter set with
	//			        only spherical tensor
	///F_list =	              - Tensor parameter conversion


MSVCDLL friend void operator+= (ParameterSet& pset, space_T &SphT);

	// Input		SphT  : Spatial tensor
	//  			pset  : Parameter set
	// Output		pset  : Parameter set with
	//			        spherical tensor added
	///F_list +=	              - Tensor parameter unary addition


MSVCDLL void operator= (const ParameterSet& pset);

	// Input		SphT	 : A spatial tensor (this)
	// 			pset     : A parameter set
	// Output		none	 : Spatial tensor filled with
	//				   data in pset
	// Note				 : Functions which place a spherical
	//				   tensor into a parameter set must
	//				   write the information read here


MSVCDLL virtual void write(const std::string &filename);

	// Input		SphT	 : Spatial tensor (this)
	//			filename : Output file name
	// Output		none 	 : Spatial tensor is written as a 
	//				   parameter set to file filename
	///F_list write	                 - Tensor write to file


MSVCDLL virtual void read(const std::string &filename);

	// Input		SphT     : Spatial tensor (this)
	// 			filename : Input filename
	// Output		none	 : Spatial tensor filled with
	//				   parameters read from file
	///F_list read	                 - Tensor read from file


// ______________________________________________________________________
// I                    SPACE TENSOR I/O FUNCTIONS
// ______________________________________________________________________


MSVCDLL friend std::ostream& operator<< (std::ostream& ostr, const space_T &SphT);

	// Input		ostr : string
	// 			SphT : spatial tensor
	// Return		     : stream, prints spatial tensor
	//			       irreducible Spherical components
	///F_list <<	             - Tensor to output stream


MSVCDLL friend void Cartesian (const space_T &SphT);

	// Input		ostr : string
	// 			SphT : spatial tensor
	// Return		     : stream, prints spatial tensor
	//			       irreducible Cartesian components
	///F_list Cartesian	     - Tensor Cartesian components to output stream

};

#endif						// SpaceT.h

