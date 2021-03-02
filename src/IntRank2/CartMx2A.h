/* CartMx2A.h ***************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Cartesian Matrix To Rank 2 Spatial Tensor           Interface	**
**                                                                      **
**      Scott Smith                                                     **
**      Copyright (c) 2001                                              **
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
** This file contains functions that allow a 3x3 matrix representing a	**
** Cartesian rank 2 tensor to be converted into a GAMMA irreducible	**
** rank 2 spherical tensor (class IntRank2A), an isotropic rank 0 	**
** component (Aiso), and a set of Euler angles that relate the rank 2	**
** PAS to the axes that were used to represent the input tensor.	**
**                                                                      **
** In short these functions are intended to do the following:		**
**                                                                      **
**      A Cartesian                       A Spherical			** 
**   -----------------        ---------------------------------------	**
**                                                                      **
**   [ A    A    A   ]        A     = Scalar, Rank 0, Isotropic Term 	**
**   |  xx   xy   xz |         iso					**
**   |               |                                                  **
**   | A    A    A   |  ===>  del   = Scaling factor			**
**   |  yx   yy   yz |           zz					**
**   |               |                                                  **
**   | A    A    A   |        eta   = Asymmetry (in Class IntRank2)	**
**   [  zz   zy   zy ]                                                  **
**			      EA    = Euler Angles {alpha,beta,gamma}	**
**                                                                      **
** The values on the right are an equivalent representation of the 	**
** Cartesian values on the left. GAMMA can easily regenerate the array 	**
** by rotating the irreducible rank 2 compoenent from the PAS to the    **
** original input orientation using the Euler angles, properly scaling  **
** it back to the original, and then adding back in the isotropic 	**
** component. This is shown below for individual elements of A 		**
**									**
**                           G						**
**                        del						**
**                           zz    G                       		**
**                  A   = ----- * A  (eta,EA) + A			**
**                   uv   del      uv            iso			**
**                           zz						**
**									**
** where the G is used to indicate spherical tensor components used	**
** internally in GAMMA.							**
**									**
** Since the conversion may not be trivial for a general Cartesian      **
** array, often an iterative fit will be performed to find the proper	**
** Euler angles. Class CartMx2A is used to control this minimization	**
** process, among other things. The Class handles the specifics of the	**
** minimization process: max. number of iteration, quality of fit,	**
** etc. Typically, the functions herein will provided both as friend	**
** and member of class CartMX2A, the latter more convenient to use.	**
**									**
*************************************************************************/

#ifndef   CartMx2A_h_			// Is file already included?
#  define CartMx2A_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface 			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/matrix.h>		// Include GAMMA matrices
#include <Matrix/col_vector.h>		// Include GAMMA column vectors
#include <Level2/EAngles.h>		// Inlcude GAMMA Euler angles
#include <vector>			// Include libstdc++ STL vectors
#include <fstream>			// Include libstdc++ streams

class CartMx2A
  {
//           Cartsian And Spherical Tensor Tensor Components

  matrix              A;		// Cartesian 3x3 rank 2 spatial tensor
  double              AisoVal;		// Spatial tensor isotropic comp.
  double              PASDelzz;		// Spatial tensor PAS delzz value
  double              EtaVal;		// Spatial tensor asymmetry [0,1]
  EAngles             EA;		// Tensor Euler Angles (from PAS)

//                     Minimization Iteration Controls 

  int                 MaxAngs;		// Maximum # angle sets       (limit)
  int                 MaxIter;		// Maximum # iteration steps  (limit)
  int                 MaxSSteps;	// Max. # repeat min steps    (limit)
  int                 Step;		// Current minimization step  (counter)
  int                 Angle;		// Current angles set         (counter)
  int                 Iter;		// Current iteration step     (counter)
  int                 JStep;		// Current Jacobian eval.     (counter)
  int                 SStep;		// Current repeat step        (counter)
  matrix              J;		// Jacobian matrix     (3x3)  (minmath)
  col_vector          F;		// Functional vector   (3x1)  (minmath)
  std::vector<double> Av;		// Vector of oriented comps.  (minmath)
  EAngles             EAo;		// Initial Euler angle guesses(minmath)
  col_vector          X;		// Minimzed angles |X> (3x1)  (minmath)
  col_vector          Y;		// Angle updates   |Y> (3x1)  (minmath)		
  double              EtaCut;		// Asymmetry zero cutoff      (cutoff)
  double              DifCut;		// Differential zero cutoff   (cutoff)
  double              NormCut;		// Minimization finish cutoff (cutoff)
  double              ZeroMxCut;	// A consistency cutoff       (cutoff)
  int                 Plevel;		// Minimization output level  (output)
  std::list<double>   Norms; 		// Iteration norm values      (track)
  bool                JOK;		// Flag if Jacobian OK        (track)
  int                 ConvMeth;		// Conversion Method	      (track)
  std::vector<int>    Flist;		// Array of function indices  (track)
  std::string         Fuse;		// Function usage in Jacobian (track)

//  Trig. Evaluations Used To Generate Jacobian Matrix & Functional Vector

  double Sa;				// sin(alpha)
  double Ca;				// cos(alpha)
  double S2a;				// sin(2*alpha)
  double C2a;				// cos(2*alpha)
  double Sb;				// sin(beta)
  double Cb;				// cos(beta)
  double S2b;				// sin(2*beta)
  double C2b;				// cos(2*beta)
  double SbCb;				// sin(beta)*cos(beta)
  double Ssqb;				// sin(beta)*sin(beta)
  double Csqb;				// cos(beta)*cos(beta)
  double Csqbp1;                        // cos(beta)*cos(beta)+1.0
  double S2g;				// sin(2*gamma)
  double C2g;				// cos(2*gamma)


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                    Class CartMx2A Error Handling
// ____________________________________________________________________________

/*       Input               C2A     : Cartesian To GAMMA A Conversion (this)
                             eidx    : Flag for error type
                             noret   : Flag for return (0=return)
                             pname   : String included in message
        Output               none    : Error message
                                       Program execution stopped if fatal   */

         void C2Aerror(int eidx, int noret=0) const;
volatile void C2Afatal(int eidx)              const;
         void C2Aerror(int eidx, const std::string& p, int n=0) const;
volatile void C2Afatal(int eidx, const std::string& p)          const;

// ____________________________________________________________________________
// ii                   Class CartMx2A SetUp Functions
// ____________________________________________________________________________

void SetDefaults();
bool SetAMx(const matrix& AC, int warn=2);
void SetAngles();

// ____________________________________________________________________________
// iii        Class CartMx2A Jacobian Matrix Singularity Test
// ____________________________________________________________________________

/* Given a possible row in a 3x3 Jacobian matrix, <i|a|0>, <i|a|1>, & <i|a|2>,
   this function does a very simple check to insure the row will not produce
   a singlular Jacobian matrix.  It does so by checking that 1.) row i is not
   zero and 2.) row j where j<i is not a constant multiple of row i. The
   latter check is done for all rows j where j<i.                            */   

bool CheckSing(double ai0, double ai1, double ai2, int i);

// ____________________________________________________________________________
// iv         Class CartMx2A Caretesian Matrix Conversion Routine
// ____________________________________________________________________________

/* This attempts to perform the conversion of the current Cartesian tensor into
   a GAMMA spherical tensor. It is separated from other routiens and private
   because it uses the current class settings which dictate how the
   conversion process will be attempted.                                     */

bool Convert(int warn=2);
 

// ____________________________________________________________________________
// v           Class CartMx2A Jacobian Function List Check
// ____________________________________________________________________________

/* This insures that the mapping of {Auv} into a Jacobian row is not improper.
   The value of findex is: {0:Axx, 1:Axy, 2:Axz, 3:Ayy, 4:Ayz, 5:Azz} and we
   make sure its range is [0,5]. There are only three rows of the Jacobian 
   so we insure that the range of frow is [0,2].                             */

bool CheckF(int findex, int frow, int warn=2) const;
bool CheckNorms() const;


// ____________________________________________________________________________
// vi          Class CartMx2A Euler Angles For Minimizaiton Search
// ____________________________________________________________________________

double               NewAlpha(int i) const;
double               NewBeta(int  i) const;
double               NewGamma(int i) const;
std::vector<EAngles> AngSeeds()      const;
std::vector<int>     FctSeeds()      const;

// ____________________________________________________________________________

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

// ____________________________________________________________________________
// A            Class CartMx2A Construction, Assignment, Destruction
// ____________________________________________________________________________

MSVCDLC      CartMx2A();
MSVCDLC      CartMx2A(const matrix& AC, int warn=2);
MSVCDLL void operator= (const CartMx2A& C2A);
MSVCDLC      ~CartMx2A();

// ____________________________________________________________________________
// B                  Class CartMx2A Access Functions
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                   Functions To Obtain Class Values
// ----------------------------------------------------------------------------

MSVCDLL matrix     ACart()         const;	// Get Cartesian Tensor
MSVCDLL double     Aiso()          const;	// Get Tensor Isotropy
MSVCDLL double     Eta()           const;	// Get Tensor Asymmetry [0,1]
MSVCDLL double     delzz()         const;	// Get Tensor PAS delzz (Anis) 
MSVCDLL EAngles    EulerAngles()   const;	// Get Tensor Euler Angles

MSVCDLL int        MaxIterations() const;	// Get Max. Allowed Iterations
MSVCDLL int        Iteration()     const;	// Get Current Iteration
MSVCDLL matrix     Jacobian()      const;	// Get Jacobian Matrix J (3x3)
MSVCDLL row_vector Functional()    const;	// Get Function Vector F (3x1) 
MSVCDLL EAngles    StartAngles()   const;	// Get Initial Euler Angles
MSVCDLL double     EtaCutoff()     const;	// Get Asymmetry Zero Cutoff
MSVCDLL double     DifCutoff()     const;	// Get Differential Zero Cutoff
MSVCDLL int        PrintLevel()    const;	// Get Current Min. Print Level
MSVCDLL int        JFunct(int i)   const;	// Get Jacobian Row i Function

// ----------------------------------------------------------------------------
//                     Functions To Set Class Values
// ----------------------------------------------------------------------------

MSVCDLL bool ACart(const matrix& AC, int w=2);// Set Cartesian Tensor
MSVCDLL void MaxIterations(int mi);		// Set Max. Allowed Iterations
MSVCDLL void EtaCutoff(double  ec);		// Set Asymmetry Zero Cutoff
MSVCDLL void StartAngles(const EAngles& EAin);// Set Initial Euler Angles
MSVCDLL void DifCutoff(double  dc);		// Set Differential Zero Cutoff
MSVCDLL void PrintLevel(int    pl);		// Set Minimization Print Level
MSVCDLL void JFuncts(int f1,int f2,int f3);	// Set 3 Jacobian row functions
MSVCDLL void JFuncts(int f1f2f3);		// Set 3 Jacobian row functions

// ____________________________________________________________________________
// C          CARTESIAN MATRIX ISOTROPY, ANISOTROPY, ASYMMETRY
// ____________________________________________________________________________

/* This function will glean the three values { Aiso, delzz, eta } from an
   input Cartesian rank 2 spatial tensor. The tensor is represented by a 3x3
   symmetric matrix. The values are determined
                                                                   A     A
          1                     2 [      1 [           ]]           xx -  yy
   A    = - Trace(A)    del   = - | A  - - | A   - A   ||    eta = ---------
    iso   3                zz   3 [  zz  2 [  xx    yy ]]              A
                                                                        zz

   where we adhere to the convention |A  | >= |A  | >= |A  | or eta = [0,1]
                                        zz       yy       xx                 */

MSVCDLL void AisoDelzEta(const matrix& A);
MSVCDLL void AisoDelzEta();

// ____________________________________________________________________________
// D                      CARTESIAN MATRIX EULER ANGLES
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//            Euler Angles For Diagonal Symmetric Cartesian Tensor
// ----------------------------------------------------------------------------

/* This function will assume that the irreducible Cartesian spatial tensor
   array A (symmetric) has no asymmetry and is diagonal. That is, eta = 0 for A
   and A contains no non-zero off-diagonal elements. It will return three Euler
   angles EA:{alpha,beta,gamma} that relate to the tensor PAS to the input
   array axes. In this case, i.e. eta=0 & A diagonal, either the tensor is
   already in its PAS or only a simple rotation is needed to switch diagonal
   element ordering so it matches the PAS.  Here are the possible Euler angles:

         EA = {0,0,0}            EA = {0,90,0}            EA = {90,90,0}

     1      [ -1  0  0 ]      1      [  2  0  0 ]      1      [ -1  0  0 ]
     - del  |  0 -1  0 |      - del  |  0 -1  0 |      - del  |  0  2  0 |
     2   zz [  0  0  2 ]      2   zz [  0  0 -1 ]      2   zz [  0  0 -1 ]  */


// ----------------------------------------------------------------------------
//                Euler Angles For Symmetric Cartesian Tensor
// ----------------------------------------------------------------------------

MSVCDLL void DiagSymCartEA();

/* This function will assume that the irreducible Cartesian spatial tensor
   array A (symmetric) has no asymmetry. That is, eta = 0 for A.  It will
   return three Euler angles EA:{alpha,beta,gamma} that relate to the tensor
   PAS to the input array axes. In this case, i.e. eta=0, the Euler angle
   gamma is ill-defined and set to zero. The two other angles alpha and beta
   may be generated from the following:

                                                        1
                                [ [   [   A         ] ] - ]
                                | | 1 |    zz       | | 2 |
                    beta = acos | | - | 2 ----- + 1 | |   |
                                | | 3 |   del       | |   |
                                [ [   [      zz     ] ]   ]

                                  [           A          ]
                                  | 4          xz        |
                     alpha = acos | - ------------------ |
                                  | 3 del  * sin(2*beta) |
                                  [      zz              ]                   */

MSVCDLL void SymCartEA();


// ----------------------------------------------------------------------------
//           Euler Angles For Diagonal Asymmetric Cartesian Tensor
// ----------------------------------------------------------------------------

/* This function will assume that the irreducible Cartesian spatial tensor
   array A (symmetric) has some asymmetry but that it is diagonal. That is,
   eta != 0 but A contains no non-zero off-diagonal elements. It will attempt
   to generate three Euler angles EA:{alpha,beta,gamma} that relate to the
   tensor PAS to the input array axes. In this case, eta=# & A diagonal, there
   are several angle combinations possible. We do not cover them all herein,
   so the function returns false if we don't find them.                     */

MSVCDLL bool DiagASymCartEA();

// ----------------------------------------------------------------------------
// Euler Angles For Asymmetric Cartesian Tensor With Some Off-Diagonals At Zero
// ----------------------------------------------------------------------------

MSVCDLL bool ASymCartODZEA(int warn=2);

// ----------------------------------------------------------------------------
//  Euler Angles For Asymmetric Cartesian Tensor With Symmetric Off-Diagonals
// ----------------------------------------------------------------------------

MSVCDLL bool ASymCartODSEA(int warn=2);

// ----------------------------------------------------------------------------
//                 Euler Angles For Generic Cartesian Tensor
// ----------------------------------------------------------------------------

/* This function will assume that the Cartesian spatial tensor array A
   (symmetric) has some asymmetry. That is, eta != 0 for A. It will attempt
   to obtain three Euler angles EA:{alpha,beta,gamma} that relate the tensor
   PAS to the input array axes. For this case, the angles must be found
   through an iterative minimization process.                                */

MSVCDLL bool ASymCartEA();
MSVCDLL int  MinAction(bool TFM, bool TFA, int NJ, int JSteps, double lastnorm) const;

// ____________________________________________________________________________
// E                 Class CartMx2A Minimization Functions
// ____________________________________________________________________________

MSVCDLL bool Minimize();


// ----------------------------------------------------------------------------
//           Jacobian Matrix & Functional Vector Generators
// ----------------------------------------------------------------------------
 
MSVCDLL void Axx(double& f, double& dfda, double& dfdb, double& dfdg) const;
MSVCDLL void Axy(double& f, double& dfda, double& dfdb, double& dfdg) const;
MSVCDLL void Axz(double& f, double& dfda, double& dfdb, double& dfdg) const;
MSVCDLL void Ayy(double& f, double& dfda, double& dfdb, double& dfdg) const;
MSVCDLL void Ayz(double& f, double& dfda, double& dfdb, double& dfdg) const;
MSVCDLL void Azz(double& f, double& dfda, double& dfdb, double& dfdg) const;

/* This function will generate a 3x3 Jacobian matrix and a 3x1 Functional
   vector for use in a Newton minimization to determine the three Euler angles
   EA:{alpha, beta, gamma} associated with an irreducible rank 2 Cartesian
   spatial tensor A. The angles relate the tensor PAS to the orientation in 
   which it exists. The Jacobian array J & Funcitonal vector F produced are


          [  @f1      @f1     @f1   ]                 [                      ]
          | ------   -----   ------ |                 | f1(alpha,beta,gamma) |
          | @alpha   @beta   @gamma |                 |                      |
          |                         |                 |                      |
          |  @f2      @f2     @f2   |                 |                      |
 J(|X>) = | ------   -----   ------ |      |F(|X>)> = | f2(alpha,beta,gamma) | 
          | @alpha   @beta   @gamma |                 |                      |
          |                         |                 |                      |
          |  @f3      @f3     @f3   |                 |                      |
          | ------   -----   ------ |                 | f3(alpha,beta,gamma) |
          [ @alpha   @beta   @gamma ]                 [                      ]


   The function takes the vector |X> containing a possible set of Euler angles.
   Were these the true angles of the problem the functional vector will be zero
   by design, i.e. |F> = 0. In addition, the user must supply the vector |A>,
   that contains 6 of the original components of the Cartesian matrix, A, being
   treated 
                       T
                    |A>  = [ Axx, Axy, Axz, Ayy, Ayz, Azz ]
  
   as well as the values for delzz and eta (which may be determined directly
   for the array A).  An option print flag allows the user to view the process
   involved in generating J and |F>.

   There are 6 possible functions that may be used for generating the rows of
   J and |F>, one for each input Cartesian element of A (Auv in vector |A>). 
   Hard-coded in are the formula for generating the elements of the Jacobian
   and functional vector for all six. However, the minimization to obtain
   the Euler angles demands we use only three of these. Hence, the routine 
   will select out three of the six that are deemed good functions.... namely
   ones that will not produce a singular Jacobian matrix. In this problem a
   singular Jacobian matrix could easily occur since the six functions are
   inter-related and may produce identical rows in J or may produce a zero
   row in J because on often finds zero angles in the vector |X>.

   For more information on the Jacobian matrix, Functional vector F, and 
   overall minimization scheme see the GAMMA documentation on class IntRank2A.


	   Input		A	: Vector - { delxx, delxy, delxz,
	  					     delyy, delyz,delzz }
	  			eta     : Asymmetry value
	   			X	: Vector - {alpha,beta,gamma} radians
	  			F       : Vector - {f1,f2,f3}
	  			J	: Jacobian Matrix (3x3)
	  			pl      : Print level for output
	   Output		TF      : Given the values in X & A the
	  				  vector F and Jacobian matrix J
	  				  are generated.  The function wil
	  				  return true if this procedure is
	  				  deemed OK, false if it appears that
	  				  one cannot build a proper Jacobian
	   Note				: FOR FITTING {alpha,beta,gamma}
	  				  WHEN WE KNOW {delzz,eta}
	   Note				: If the function returns false, one
	  				  should try a new starting vector   */

MSVCDLL bool JacobianF();

MSVCDLL void TrackMin(int info, double val=0) const;

// ____________________________________________________________________________
// F                  Class CartMx2A Auxiliary Functions
// ____________________________________________________________________________

MSVCDLL matrix      Regenerate()         const;
MSVCDLL double      Check()              const;
MSVCDLL bool        Check(double cutoff) const;
MSVCDLL std::string Type()               const;
MSVCDLL std::string Method()             const;
MSVCDLL std::string JFunctName(int i)    const;
MSVCDLL void        TrackConv(int info)  const;
MSVCDLL bool        OffDiagonals()       const;
MSVCDLL bool        SymOffDiags()        const;

// ____________________________________________________________________________
// Z               Class CartMx2A Formatted Output Functions
// ____________________________________________________________________________

/* These functions will output information concerning the conversion of the
   Cartesian spatial tensor into a GAMMA spatial tensor. The information will
   be output into any specified output stream.

		Input		C2A     : Cartesian 2 GAMMA A Conversion (this)
                                ostr    : Output stream
                                fflag   : Format flag
                                            0 - Sum rotation only
                                           !0 - Full composite rotation
           Output               none    : Conversion information placed
                                          into the output stream            */


MSVCDLL        std::ostream& print(std::ostream& ostr, int   fflag=0) const;
MSVCDLL friend std::ostream& operator<<     (std::ostream& ostr, const CartMx2A& C2A);

  };						// End Of Class CartMx2A






// ____________________________________________________________________________
// i               CARTESIAN MATRIX JACOBIAN SINGULARITY TEST
// ____________________________________________________________________________

/* Given a possible row in a 3x3 Jacobian matrix, <i|a|0>, <i|a|1>, & <i|a|2>,
   this function does a very simple check to insure the row will not produce
   a singlular Jacobian matrix.  It does so by checking that 1.) row i is not
   zero and 2.) row j where j<i is not a constant multiple of row i. The
   latter check is done for all rows j where j<i.                            */   

MSVCDLL bool CheckSing(double ai0, double ai1, double ai2, const matrix& J, int i);
 
// ____________________________________________________________________________
// A     CARTESIAN MATRIX GENERATION OF JACOBIAN MATRIX & FUNCTIONAL VECTOR
// ____________________________________________________________________________

/* This function will generate a 3x3 Jacobian matrix and a 3x1 Functional
   vector for use in a Newton minimization to determine the three Euler angles
   EA:{alpha, beta, gamma} associated with an irreducible rank 2 Cartesian
   spatial tensor A. The angles relate the tensor PAS to the orientation in 
   which it exists. The Jacobian array J & Funcitonal vector F produced are


          [  @f1      @f1     @f1   ]                 [                      ]
          | ------   -----   ------ |                 | f1(alpha,beta,gamma) |
          | @alpha   @beta   @gamma |                 |                      |
          |                         |                 |                      |
          |  @f2      @f2     @f2   |                 |                      |
 J(|X>) = | ------   -----   ------ |      |F(|X>)> = | f2(alpha,beta,gamma) | 
          | @alpha   @beta   @gamma |                 |                      |
          |                         |                 |                      |
          |  @f3      @f3     @f3   |                 |                      |
          | ------   -----   ------ |                 | f3(alpha,beta,gamma) |
          [ @alpha   @beta   @gamma ]                 [                      ]


   The function takes the vector |X> containing a possible set of Euler angles.
   Were these the true angles of the problem the functional vector will be zero
   by design, i.e. |F> = 0. In addition, the user must supply the vector |A>,
   that contains 6 of the original components of the Cartesian matrix, A, being
   treated 
                       T
                    |A>  = [ Axx, Axy, Axz, Ayy, Ayz, Azz ]
  
   as well as the values for delzz and eta (which may be determined directly
   for the array A).  An option print flag allows the user to view the process
   involved in generating J and |F>.

   There are 6 possible functions that may be used for generating the rows of
   J and |F>, one for each input Cartesian element of A (Auv in vector |A>). 
   Hard-coded in are the formula for generating the elements of the Jacobian
   and functional vector for all six. However, the minimization to obtain
   the Euler angles demands we use only three of these. Hence, the routine 
   will select out three of the six that are deemed good functions.... namely
   ones that will not produce a singular Jacobian matrix. In this problem a
   singular Jacobian matrix could easily occur since the six functions are
   inter-related and may produce identical rows in J or may produce a zero
   row in J because on often finds zero angles in the vector |X>.

   For more information on the Jacobian matrix, Functional vector F, and 
   overall minimization scheme see the GAMMA documentation on class IntRank2A.


	   Input		A	: Vector - { delxx, delxy, delxz,
	  					     delyy, delyz,delzz }
	  			eta     : Asymmetry value
	   			X	: Vector - {alpha,beta,gamma} radians
	  			F       : Vector - {f1,f2,f3}
	  			J	: Jacobian Matrix (3x3)
	  			pl      : Print level for output
	   Output		TF      : Given the values in X & A the
	  				  vector F and Jacobian matrix J
	  				  are generated.  The function wil
	  				  return true if this procedure is
	  				  deemed OK, false if it appears that
	  				  one cannot build a proper Jacobian
	   Note				: FOR FITTING {alpha,beta,gamma}
	  				  WHEN WE KNOW {delzz,eta}
	   Note				: If the function returns false, one
	  				  should try a new starting vector   */

MSVCDLL bool CMx2AJacobian(const std::vector<double>& A, double delz, double eta,
                      const col_vector& X, col_vector& F, matrix& J, int pl=0);

// ____________________________________________________________________________
// B          CARTESIAN MATRIX MINIMIZATION TO FIND EULER ANGLES
// ____________________________________________________________________________

/*
	   Input		A	: Vector - { delxx, delxy, delxz,
	  					     delyy, delyz,delzz }
				delzz   : Scaling factor for A (GAMMA)
	  			eta     : Asymmetry value
	   			X	: Vector - {alpha,beta,gamma} radians
				maxiter	: Maximum allowed iterations
				pl	: Print level for output
	   Output		TF      : Given the values in X & A the
	  				  set of Euler Angles that relate
                                          A PAS to A are found via Newton
                                          minimization. The function will
	  				  return true if this procedure is
	  				  successful, false if it appears that
	  				  one cannot find the proper angles
	   Note				: Proper angles returned in X
	   Note				: FOR FITTING {alpha,beta,gamma}
	  				  WHEN WE KNOW {delz z,eta}          */


MSVCDLL int CartMxMinimize(const std::vector<double>& A, double delz, double eta,
                                        col_vector& X, int maxiter, int pl=0);

// ____________________________________________________________________________
// C          CARTESIAN MATRIX ISOTROPY, ANISOTROPY, ASYMMETRY
// ____________________________________________________________________________

/* This function will glean the three values { Aiso, delzz, eta } from an
   input Cartesian rank 2 spatial tensor. The tensor is represented by a 3x3
   symmetric matrix. The values are determined
                                                                   A     A
          1                     2 [      1 [           ]]           xx -  yy
   A    = - Trace(A)    del   = - | A  - - | A   - A   ||    eta = ---------
    iso   3                zz   3 [  zz  2 [  xx    yy ]]              A
                                                                        zz

   where we adhere to the convention |A  | >= |A  | >= |A  | or eta = [0,1]
                                        zz       yy       xx                 */

MSVCDLL void AisoDelzEta(const matrix& A, double& AIso, 
         double& Adelzz, double& Aeta, double ecut=1.e-5, double delcut=1.e-9);

// ____________________________________________________________________________
// D              SYMMETRIC CARTESIAN MATRIX EULER ANGLES
// ____________________________________________________________________________

/* This function will assume that the irreducible Cartesian spatial tensor
   array A (symmetric, traceless) has no asymmetry. That is, eta = 0 for A.
   It will return three Euler angles EA:{alpha,beta,gamma} that related the
   tensor PAS to the input array axes. In this case, i.e. eta=0, the Euler
   angle gamma is ill-defined and set to zero. The two other angles alpha and
   beta may be generated from the following:

                                                        1
                                [ [   [   A         ] ] - ]
                                | | 1 |    zz       | | 2 |
                    beta = acos | | - | 2 ----- + 1 | |   |
                                | | 3 |   del       | |   |
                                [ [   [      zz     ] ]   ]

                                  [           A          ]
                                  | 4          xz        |
                     alpha = acos | - ------------------ |
                                  | 3 del  * sin(2*beta) |
                                  [      zz              ]                   */


MSVCDLL EAngles SymCartEA(const matrix& A, double delzz, double Aiso=0.0);

// ____________________________________________________________________________
// E              GENERAL CARTESIAN MATRIX EULER ANGLES
// ____________________________________________________________________________

/* This function will assume that the Cartesian spatial tensor array A
   (symmetric) has some asymmetry. That is, eta != 0 for A. It will return
   three Euler angles EA:{alpha,beta,gamma} that relate the tensor PAS to
   the input array axes. For this case, the angles must be found through a
   minimization process.                                                     */

MSVCDLL int ASymCartEA(const matrix& A, double Aiso, double delzz, double eta,
                                              const EAngles& EAo, EAngles& EA);

// ____________________________________________________________________________
// F           CARTESIAN MATRIX CONVERSION TO GAMMA SPATIAL TENSOR
// ____________________________________________________________________________


//void CartMx2A(const matrix& A,
//                        double& Aiso, double& delzz, double& eta, EAngles& EA);

#endif								// CartMx2A.h
