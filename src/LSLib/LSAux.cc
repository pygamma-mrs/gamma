/* LSAux.cc ******************************************************
**								**
** 	                    G A M M A				**
**								**
**	Liouville Space Auxiliary Functions    Implementation 	**
**							 	**
**	Copyright (c) 1991, 1992, 1993		 		**
**	Scott A. Smith					 	**
**								**
**	Eidgenoessische Technische Hochschule	 		**
**	Labor fuer physikalische Chemie		 		**
**	8092 Zurich / Switzerland		 		**
**								**
**	University of California, Santa Barbara			**
**	Department of Chemistry					**
**	Santa Barbara CA. 93106 USA				**
**						 		**
**      $Header: $
**							 	**
*****************************************************************/

/*****************************************************************
**							 	**
** 	Description					 	**
**							 	**
** Functions are helpful in dealing with Liouville space.	**
**							 	**
*****************************************************************/

#ifndef _LSaux_cc_		// Is this file already included?
#define _LSaux_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation          // This is the implementation
#endif

#include <LSLib/LSAux.h>		// Include the interface
#include <LSLib/SuperOp.h>		// Include superoperators
#include <HSLib/HSauxil.h>		// Include sorting routines
#include <HSLib/SpinSys.h>              // Know about spin systems
#include <HSLib/GenOp.h>                // Know about general operators
#include <Matrix/matrix.h>              // Know about matrices
#include <Basics/StringCut.h>
#include <stdlib.h>

// ____________________________________________________________________________
// *************** General Relaxation Auxiliary Functions ***************
// ____________________________________________________________________________
  
void print(const super_op& LOp, double cutoff, int nc, int ri)

	// Input		LOp   : Superoperator
	//			cotoff: Cutoff level for printing
	//			nc    : Number of columns
	//			ri    : Flag for Real(0), Complex(1)
        // Return		None  : LOp elements sent to stand. I/O
        // Note			      : Assumes the matrix LOp is real

  {
  int ls = LOp.dim();				// Get Liouville space size
  int iacc = 1;					// Set the accuracy of the
  if(ls >= 10)					// basis indicies
    iacc = 2;
  if(ls >= 100)
    iacc = 3;
  if(ls >= 1000)
    iacc =4;
  std::string Siacc = std::string("%i") + Gdec(iacc);
  complex Relc;
  double Reld;
  int k=0, np=0;
  if(ls) 					// Check for NULL LOp
    {
    for(int i=0; i< ls; i++)
      for(int j=0; j< ls; j++)
        {
        if(ri == 1)				// Deal with complex elements
          {
          Relc = LOp.get(i,j);
          if(norm(Relc) > cutoff) 
            {
            std::cout << "<" << Gdec(Siacc, i) << "|" << Gdec(Siacc, j)
                 << "> = " << Re(Relc)
                 << ", " << Im(Relc) << "\t";
            k++;
            np++;
            }
          }
        else if(ri==-1)				// Deal with imaginary elements
          {
          Reld = Im(LOp.get(i,j));
          if(fabs(Reld) > cutoff) 
            {
            std::cout << "<" << Gdec(Siacc, i) << "|" << Gdec(Siacc, j)
                 << "> = " << Reld << "\t";
            k++;
            np++;
            }
          }
        else					// Deal with real elements
          {
          Reld = Re(LOp.get(i,j));
          if(fabs(Reld) > cutoff) 
            {
            std::cout << "<" << Gdec(Siacc, i) << "|" << Gdec(Siacc, j)
                 << "> = " << Gform("%e12.5", Reld) << "\t";
//                 << "> = " << setw(12) << setprecision(5)
//                 << Reld << "\t";
// sosi - need to set the form to 'e'
            k++;
            np++;
            }
          }
        if(k == nc)				// Output nc elements per line
          {
          std::cout << "\n";
          k = 0;
          }
        }
    if(np == 0)
      std::cout << "\n\tNull Superoperator Within " << cutoff << " Magnitude Cutoff\n";
    }
  else
    std::cout << "\n\tNull Superoperator\n";
  return;
  }

  
void eigenvalues(super_op& LOp, int sort, int nc, int ri)

	// Input		LOp   : Superoperator
	//			sort  : Flag to initiate sorting
	//			nc    : Number of columns
	//			ri    : Flag Real(0), Complex(1), Im(-1)
        // Return		None  : LOp eigenvalues to standard output
        // Note			      : Sets EBR of LOp

  {
  int ls = LOp.dim();			// Get Liouville space size
  complex maxe, *ev;
  ev = new complex[ls];
  double max;
  int maxi;
  if(ls) 				// Check for NULL LOp
    {
    LOp.set_EBR();			// 1st put LOp into its EBR
    int i, j=0;
    for(i=0; i< ls; i++)		// Copy eigenvalues into ev
      ev[i] = LOp.get(i,i);

    if(sort)				// Sort the eigenvalues if needed
      {
      for(i=0; i< ls-1; i++)
        {
        if(ri == 0)			// Sort on complex eigenvalues
          max = norm(ev[i]);
        else if(ri == -1)
          max = Im(ev[i]);		// Sort on imaginary eigenvalues
        else
          max = Re(ev[i]);		// Sort on real eigenvalues
        maxe = ev[i];
        maxi = i;
        for(j=i+1; j<ls; j++)
          {
          if(ri==1 && norm(ev[j]) > max)
            {
            maxe = ev[j];
            maxi = j;
            max = norm(ev[j]);
            }
          else if(ri==-1 && Im(ev[j]) > max)
            {
            maxe = ev[j];
            maxi = j;
            max = Im(ev[j]);
            }
          else if(Re(ev[j]) > max)
            {
            maxe = ev[j];
            maxi = j;
            max = Re(ev[j]);
            }
          }
        if(maxi != i)
          {
          ev[maxi] = ev[i];
          ev[i] = maxe;
          }
        }
      }

    j=0;
    int iacc = 1;					// Set the accuracy of the
    if(ls >= 10)					// basis indicies
      iacc = 2;
    if(ls >= 100)
      iacc = 3;
    if(ls >= 1000)
      iacc =4;
    std::string Siacc = std::string("%i") + Gdec(iacc);
    for(i=0; i< ls; i++)	 	// Now print the eigenvalues
      {
      if(ri == 1)
        std::cout << Gdec(Siacc, i) << ". " << Re(ev[i])
             << ", " << Im(ev[i]) << "\t";
      else if(ri == -1)
        std::cout << Gdec(Siacc, i) << ". "
             << Gform("%e12.5", Im(ev[i])) << "\t";
// sosi - need to set the form to 'e'
//             << setw(12) << setprecision(5) << Im(ev[i]) << "\t";
      else 
// sosi - need to set the form to 'e'
        std::cout << Gdec(Siacc, i) << ". "
            << Gform("%e12.5", Re(ev[i])) << "\t";
//             << setw(12) << setprecision(5) << Re(ev[i]) << "\t";
      j++;
      if(j >= nc)		// Output nc eigenvalues per line
        {
        std::cout << "\n";
        j = 0;
        }
      }
    }
  else
    std::cout << "\nNull superoperator, all zero eigenvalues";
  delete [] ev;
  return;
  }



// ____________________________________________________________________________
// ************************** COHERENCE SORTING *******************************
// ____________________________________________________________________________

//??? Note: These are very inefficient routines for three reasons:
//	1.) The matrix which does the sorting is an integer array
//	    but is complex as used here.  Thus the adjoint is taken
//	    in one function rather than the transpose.
//	2.) Rather than a permutation matrix which requires Liouville
//	    space multiplications by 0 & 1 a function in matrix should
//	    be provided which does simple column/row exchanges instead
//	    given the sorted array (called x from sort_LOp_basis)
//	3.) Very little thought has gone into these routines in order
//	    to conserve superoperator storage.  These routines work but
//	    should be redesigned when better matrix capabilites are
//	    available and the relaxation schemes are more soundly designed.

matrix UOrderMQC(const spin_sys &sys)

	// Input		sys	: Spin system
	// Output		U	: Transformation matrix
	//				  for coherence ordering

  {
  int* x = sort_LOp_basis(sys);		// Get LOp sorted basis
  int HS=sys.HS();			// Get Hilbert Space
  int LS = HS*HS;
  matrix U(LS,LS,complex0);		// Create permutation matricx
  for(int i=0; i<LS; i++)
  U.put(complex1,i,x[i]);
  return U;				// Return permutation matrix
  }


super_op OrderMQC(const super_op &LOp, const matrix &U)

	// Input		LOp	: Superoperator
	//			U	: Transformation matrix
	// Output		R	: Superoperator LOp which has
	//				  been coherence ordered

  {
  matrix Lmx = LOp.get_mx();
  return super_op(U*times_adjoint(Lmx,U));
  }


super_op OrderMQC(super_op &LOp, spin_sys &sys)

	// Input		LOp	: Superoperator
	//			sys	: Spin system
	//			R	: Superoperator LOp which has
	//				  been coherence ordered according
	//				  to sys default basis order

  { return (OrderMQC(LOp, UOrderMQC(sys))); }


/*
super_op HsuperX(gen_op& Heff)

	// Input		Heff  : Effective Hamiltonian (Hz)
	// Output		LOp   : Hamiltonian commutation superoperator
	//
	//       <a,aa|LOp|b,bb> = del  * del     <b|H   |b> - <bb|H   |bb>
	//			      ab     aa,bb    eff           eff
	//
	// Note			      :	LOp is returned in angular frequency
	//				units

  {
  basis bs = Heff.get_basis();                    // Get Hilbert basis
  int nc = bs.sub_N();                          // Number of sub-spaces
  if(nc > 1) return Hsuper(Heff);
  const double pi2 = 6.283185307;
  int hs = Heff.dim();			// Get Hilbert space size
  int ls = hs*hs;
  Heff.set_EBR();
  matrix mx(ls,ls,0,d_matrix_type,_hermitian);	// Construct zero superoperator
  basis Hbs = Heff.get_basis();
  super_op LOp(mx, Hbs);
  int aaa=0, bbb=0;
  double wbbb;
  for(int a=0; a<hs; a++)			// Sum over transition a-aa
    for(int aa=0; aa<hs; aa++)
      {
      bbb = 0;
      for(int b=0; b<hs; b++)			// Sum over transition b-bb
        {
        wbbb = Re(Heff.get(b,b));
        for(int bb=0; bb<hs; bb++)
          {
          if(a==b && aa==bb)
            {
            wbbb -= Re(Heff.get(bb,bb));	// wbbb = wb - wb (in Hertz)
            LOp.put(aaa,bbb,complex(pi2*wbbb));// Put in matrix element (rad/sec)
            }
          bbb++;
          }
        }
      aaa++;
      }
  return LOp;
  }


super_op Hsuper(gen_op& Heff)

	// Input		Heff  : Effective Hamiltonian (Hz)
	// Output		LOp   : Hamiltonian commutation superoperator
	//
	//       <a,aa|LOp|b,bb> = del  * del     <b|H   |b> - <bb|H   |bb>
	//			      ab     aa,bb    eff           eff
	//
	//
	// Note			      :	LOp is returned in angular frequency
	//				units

  {
  Heff.set_EBR();				// Put Heff into eigenbasis
  matrix mx = Heff.get_mx();			// Get Hilbert matrix
  basis bs = Heff.get_basis();			// Get Hilbert basis
  int nc = bs.sub_N();                          // Number of sub-spaces
  int LS = bs.dim_LS();				// (Composite) Liouville space
  matrix LMx(LS,LS,0,d_matrix_type,_hermitian);	// Liouville space zero array

  int hs, hsst=0;				// Hilbert space start
  int a, aa, aaa=0, b, bb;			// Energy,transition indices
  complex wb, wbb, wbbb;
  for(int cmp=0; cmp<nc; cmp++)			// Loop subspaces (components)
    {
    hsst = bs.sub_anchor(cmp);			// Start of sub-space
    hs = bs.sub_dim(cmp);			// Dimension of sub-space
    for(a=0; a<hs; a++)				// Sum over transition a-aa
      {
      for(aa=0; aa<hs; aa++)
        {
        for(b=0; b<hs; b++)			// Sum over transition b-bb
          {
          wb = Heff.get(hsst+b,hsst+b);		//	Energy eigenstate b (Hz) 
          for(bb=0; bb<hs; bb++)
            {
            if(a==b && aa==bb)
              {
              wbb = Heff.get(hsst+bb,hsst+bb);	// 	Energy eigenstate bb (Hz)
              wbbb = PI2*(wb-wbb);		//	Transition bb-b (rad/sec)
              LMx.put(wbbb, aaa, aaa);		//      Set diagonal element (rad/sec)
              }
            }
          }
        aaa++;
        }
      }
    }
  super_op LOp(LMx, bs);
  return LOp;
  }
*/


// ____________________________________________________________________________
// ************ MATRIX FUNCTIONS THAT SHOULD BE IN CLASS MATRIX ***************
// ____________________________________________________________________________


matrix solve_it(matrix& X, matrix& Uguess, matrix& b, int lim)

	// Input		X     : Iteration matrix
	//			Uguess: Initial guess for U
	//			b     : Inhomogeneous part
	//			lim   : Maximum number of iterations allowed
        // Return		Ui    : Solution obtained by iteration to
	//
	//				     |U > = X|U   > + |b>
 	//			               i       i-1

  {
  matrix Ui, Uim1;
  Uim1 = Uguess;
std::cout << "\nIteration 0" << " matrix\n" << Uim1; 
  for(int i=1; i<=lim; i++)
    {
    Ui = X * Uim1;
    Ui += b;
std::cout << "\nIteration " << i << " matrix\n" << Ui; 
    Uim1 = Ui;
    }
  return Ui;
  }


matrix invert_it(matrix& X)

	// Input		X     : Matrix
        // Return		Y     : Inverse of X
	// Note			      : Assumes X is Hermitian

  {
std::cout << "\nIn Function invert_it\n";
  int n = X.rows();
  matrix Xeval(X);
  matrix Xevec(X);
std::cout << "\nMatrix to be Inverted\n" << X;
  diag(X,Xeval,Xevec);			// Diagonalize X
std::cout << "\nEigenvalues of Matrix to be Inverted\n" << Xeval;
std::cout << "\nEigenvectors of Matrix to be Inverted\n" << Xeval;
std::cout << "\nEigenvectors * Eigenvectors Inverse\n" << Xevec*adjoint(Xevec);
  matrix Yeval(n,n,0.0);
  for(int i=0; i<n; i++)
    Yeval.put(1.0/Xeval.get(i,i), i,i);
  matrix Y=adjoint(Xevec)*Yeval*Xevec;
std::cout << "\nEigenvalues of Inverse\n" << Yeval;
std::cout << "\nInverse\n" << Y;
std::cout << "\nProduce X * X Inverse\n" << X*Y;
  return Y;
  }


/* LU Decomposition *************************************-*-c++-*-
**							 	**
**  This routine is taken from "Numerical Recipies", W.H.Press,	**
**  B.P.Flannery, S.A.Teukolsky, & W.T.Vetterling, Cambridge	**
**  University Press, 1989.  See pages 31 through 36.  The code **
**  has been modified for GAMMA into C++.			**
**							 	**
**  The algorithm takes an input matrix A & overwrites it with  **
**  the LU decompositions of A' where A' is A with any needed	**
**  permutations. Any permutations used	are stored in the	**
**  integer array indx.  The function returns the integer d	**
**  which will be +1 if an even	number of permutations was used **
**  or -1 if an odd number of permutations was necessary.  Sub- **
**  sequently, d can be used in taking the determinant of the	**
**  input matrix A (or LU).					** 
**								**  
**  The LU decomposition formulates the equation		**
**								**  
**		     	        A = LU				**  
**								**  
**  where L is lower triangular and U is upper triangular.	**
**								**
**  [A11 A12 A13 A14]   [L11  0   0   0 ] [U11 U12 U12 U14]	**
**  |A21 A22 A23 A24]   [L21 L22  0   0 ] [ 0  U22 U23 U24|	**
**  |A31 A32 A33 A34] = [L31 L32 L33  0 ] [ 0   0  U33 U34|	**
**  [A41 A42 A43 A44]   [L41 L42 L43 L44] [ 0   0   0  U44]	**
**								**
**  Both L and U have elements on the diagonal, but those of L	**
**  are always set to be 1: <i|L|i> = 1.			**
**								**
**  [A11 A12 A13 A14]   [ 1   0   0   0 ] [U11 U12 U12 U14]	**
**  |A21 A22 A23 A24]   [L21  1   0   0 ] [ 0  U22 U23 U24|	**
**  |A31 A32 A33 A34] = [L31 L32  1   0 ] [ 0   0  U33 U34|	**
**  [A41 A42 A43 A44]   [L41 L42 L43  1 ] [ 0   0   0  U44]	**
**								**
**  Neglecting the diagonal of L, both L and U can be overlaid	**
**  for storage in a single array (here overwriting A).		**
**								**
**                      [U11 U12 U12 U14]			**
**  			[L21 U22 U23 U24|			**
**  			[L31 L32 U33 U34|			**
**  			[L41 L42 L43 U44]			**
**								**
**  Finally, the algorithm uses Crouts method to determine the	**
**  elements of L and U.  This sets the diagonal elements of L	**
**  to 1 and follows a specific order in computation of the L	**
**  and U elements which allows A to be overwritten as L and U	**
**  are determined.  For Crouts method to be stable, partial	**
**  pivoting is used (row interchanges).			**
**								**
*****************************************************************/

int LU_decomp(matrix& A, int* indx)

	// Input		A     : Input matrix
	//			indx  : Integer array for permutations
        // Return		A     : LU Decomposition of input A
        // 			indx  : Index of Permutations
	// Note			      : Assumes A is real, non-singular
	//				and square

  {
  int n = A.rows();			// Get dimension of A
  const double tiny = 1.0e-20;		// Set this to a small number
  double *vv;				// Vector for row scaling factors
  vv = new double[n];
  int d = 1;				// Permutation counter (1 even, -1 odd)

  double big, temp;
  int i;
  for(i=0; i<n; i++)			// For each row of A find the
    {					// the largest element for pivoting 
    big=0.0;				// and store corresponding scaling
    for (int j=0; j<n; j++)		// factor in array vv
      {
      temp = norm(A.get(i,j));
      if (temp > big)
        big=temp;
      }
    if(big == 0.0)			// Cannot have a zero row in A!
      {
      std::cout << "\nSingular matrix input\n";
      exit(-1);
      }
    vv[i]=1.0/big;
    }

  complex Uij, Lij;			// Element of U and L
  int imax=0;
  double dum;
  complex dumz;
  for (int j=0; j<n; j++)		// Begin going through each A column
    {

// First implement equation below to compute the values of U (except <i|U|i>)
//
//		       i-1
//		       ---
// <i|U|j> = <i|A|j> - \   <i|L|k><k|U|j>
// 		       /
//		       ---
//		       k=0
//
// This is equation 2.3.12 on page 34 of the referenced text and its use is
// restricted to terms above the diagonal.

    for (i=0; i<j; i++)
      {
      Uij = A.get(i,j);
      for(int k=0; k<i; k++)
        Uij -= A.get(i,k)*A.get(k,j);
      A.put(Uij,i,j);
      }

// Next implement the equation below to compute the values of U (i>j, i<n)
// For the case when i=j, this equation is the same as thh previous equation
// when the scaling factor is neglected (and this is done).
//
//			       j-1
//	        1	       ---
// <i|L|j> = ------- <i|A|j> - \   <i|L|k><k|U|j>
// 	     <j|U|j>	       /
//			       ---
//			       k=0

    big=0.0;				// Used to look for largest pivot
    for (i=j; i<n; i++)
      {
      Lij = A.get(i,j);
      for(int k=0; k<j; k++)
	Lij -= A.get(i,k)*A.get(k,j);
      A.put(Lij,i,j);
      dum = vv[i] * norm(Lij);		// For looking at the pivot
      if (dum >= big)			// Use it if better that big
        {
	big=dum;
	imax=i;
	}
      }

    if (j != imax)			// This part interchanges the
      {					// rows of A via permutation
      for(int k=0; k<n; k++)		// The sign of d is switched
        {				// and the permutation stored
        dumz = A.get(imax,k);		// in the array indx
        A.put(A(j,k),imax,k);
	A.put(dumz,j,k);
	}
      d = -d;
      vv[imax] = vv[j];
      }

    indx[j] = imax;
    if(norm(A.get(j,j)) == 0.0)		// If <j|A|j> is zero make it
      A.put(tiny,j,j);			// small instead!
    if (j != n-1)			// Divide by the pivot element
      {					// If <j|A|j> was zero it has 
      dumz = 1.0/A.get(j,j);		// since been set to tiny so it
      for(i=j+1; i<n; i++)		// won't blow up!
        A.put(dumz*A(i,j),i,j);
      }
    }					// Return to treat next column
  delete [] vv;
  return d;
  }


/* LU Back Substitution *********************************-*-c++-*-
**							 	**
**  This routine is taken from "Numerical Recipies", W.H.Press,	**
**  B.P.Flannery, S.A.Teukolsky, & W.T.Vetterling, Cambridge	**
**  University Press, 1989.  See pages 31 through 38.  The code **
**  has been modified for GAMMA into C++.			**
**							 	**
**  The algorithm takes an input matrix A' in its LU decomp-	**
**  osition format where A' is some original matrix A with any	**
**  needed row permutations to attain the LU form. The row	**
**  permutations used to relate A to A' are stored in the	**
**  integer array indx.  The function the proceeds to solve	**
**								**  
**			A|x> = |b>				**
**								**  
**  for |x> by considering the problem as			** 
**								**  
**              A|x> = (LU)|x> = L(U|x>) |b> 			**
**								**  
**  and first solving for the vector |y> where |y> = U|x>	**
**								**  
**			L|y> = |b>				**
**								**
**  followed by solving for |x> by				**
**								**
**			U|x> = |y>				**
**								**
**  Due to the triagular nature of the arrays L and U it is	**
**  relatively easy to solve the vector equations involving	**
**  them.  The only further complexity in this algorithm is	**
**  that the input |b> must be first permuted to |b'> so that	**
**  its element ordering match that of A'.  Following the	**
**  solution in which |x'> is obtained it is then un-permuted	**
**  to match the original A ordering.  The first equation	**
**  to be solved (actually involving L' not L) appears as	**
**								**
** 		 [ 1   0   0   0 ] [y1] = [b1]			**
** 		 [L21  1   0   0 ] [y2] = |b2|			**
** 		 [L31 L32  1   0 ] [y3] = |b3|			**
** 		 [L41 L42 L43  1 ] [y4] = [b4]			**
**								**
**  because the diagonal elements of the matrix L are always 1.	**
**  The first element of |y> (actually |y'>) is given by	**
**								**
**			y1 = b1					**
**								**
**  and then subseqent elements of this vector are found by	**
**								**
**			      [	     ---	  ]		**
**		         1    |      \	 	  |		**
**	         yi = ------- | b  - /  <i|L|j>y  |		**
**		      <i|L|i> [  i   ---	j ]		**
**								**
*****************************************************************/

void LU_backsub(matrix &ALU, int* indx, matrix& b)

	// Input		ALU   : Input matrix A in LU form
	//			indx  : Integer array for permutations
	//			b     : Input column vector b
        // Return		      :
	// Note			      :

  {
  int n = ALU.rows();			// Get dimension of array
  int ii=0, ip=0;
  complex sum=0;

//		First Solve L|y> = |b> For |y>

  int i;
  for(i=0; i<n; i++)			// Works forwards, |y(0)>
    {					// solved first because L
    ip = indx[i];			// is lower triangular
    sum = b.get(ip,0);			// A is input permuted! Get
    b.put(b.get(i,0),ip,0);		// intial b's in correct order
    if (ii)
      for (int j=0; j<i; j++)
        sum -= ALU.get(i,j)*b.get(j,0);
    else if(norm(sum))
      ii=1;
    b.put(sum,i,0);
    }

//		Now Solve U|x> = |y> For |x>

  for(i=n-1; i>=0; i--)			// Works backwards, |x(n-1)>
    {					// solved first because U
    sum = b.get(i,0);			// is upper-triangular
    for(int j=i+1; j<n; j++)
      sum -= ALU.get(i,j)*b.get(j,0);
    b.put(sum/ALU.get(i,i),i,0);	// This is now |x(i)>
    }
  return;
  }


/* Matrix Inversion *************************************-*-c++-*-
**							 	**
**  This routine is taken from "Numerical Recipies", W.H.Press,	**
**  B.P.Flannery, S.A.Teukolsky, & W.T.Vetterling, Cambridge	**
**  University Press, 1989.  See pages 31 through 38.  The code **
**  has been modified for GAMMA into C++.			**
**							 	**
**  The algorithm takes an input matrix A and intially forms	**
**  the LU decomposition of A' where A' is A with any needed	**
**  row permutations. Then, using the LU decomposition and the	**
**  identity matrix it solves for the inverse of A via back	**
**  substitution.
**								**  
**  The LU decomposition formulates the equation		**
**								**  
**		     	    A = LU				**  
**								**  
**  where L is lower triangular and U is upper triangular.	**
**								**
**  [A11 A12 A13 A14]   [L11  0   0   0 ] [U11 U12 U12 U14]	**
**  |A21 A22 A23 A24]   [L21 L22  0   0 ] [ 0  U22 U23 U24|	**
**  |A31 A32 A33 A34] = [L31 L32 L33  0 ] [ 0   0  U33 U34|	**
**  [A41 A42 A43 A44]   [L41 L42 L43 L44] [ 0   0   0  U44]	**
**								**
**  The inverse of A is then found according to the equation	**
**								**
**		         -1	 -1				**
**		       AA   = LUA   = I				**  
**								**  
**  This is done in a columnwise fashion.  Using subscript i	**
**  to indicate matrix columns, 				**
**								**
**		       -1	 -1				**
**	            A|A  > = LU|A  >  = |I >			**  
**		       i	 i	  i			**  
**								**  
**  Columns of the inverse are found by back substitution.	**
**								**
*****************************************************************/


matrix LU_invert(matrix& A)
//return Ainv(A.rows(), A.rows(), 0.0);

	// Input		A     : Input matrix
        // Return		Ainv  : Inverse of input A
	// Note			      : Assumes A is , non-singular
	//				and square

  {
  matrix ALU = A;			// Copy matrix A
  int n = A.rows();			// Get dimension of A
//  int indx[n];				// Set up indexing vector
  int *indx;
  indx = new int[n];
  LU_decomp(ALU, indx);			// LU Decomposition of A
  matrix Ii(n,1,0.0), Ainv(n,n,0.0);	// Matrix for columns of I, |Ii>       -1
  matrix Ainvi;				// Matrix for columns of A inverse, |Ai  >
  for(int i=0; i<n; i++)		// Go through the columns of I,    -1
    {					// |Ii>, solve for the columns of A  .
    Ii.put(1.0, i, 0);			//    Form column of I, |Ii>
    Ainvi = Ii; 			//			      -1 
    LU_backsub(ALU, indx, Ainvi);	//    Get column of Ainv: A|Ai  > = |Ii>
    Ainv.put_block(0, i, Ainvi);	//    Fill column of Ainv
    Ii.put(0.0, i, 0);			//    Unform column of I
    }
  delete [] indx;
  return Ainv;
  }

#endif 						// LSaux.cc

