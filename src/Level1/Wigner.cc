/* Wigner.cc ****************************************************-*-c++-*-
**									**
**	                        G A M M A				**
**								 	**
**	Wigner Rotation Matrices		    Implementation	**
**						 			**
**	Copyright (c) 1990			 			**
**	Scott Smith						 	**
**	Eidgenoessische Technische Hochschule			 	**
**	Labor fuer physikalische Chemie				 	**
**	8092 Zurich / Switzerland				 	**
**						 			**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**								 	**
** Description							 	**
**						 			**
** This module provides functions to access Wigner rotation matrices 	**
** and Wigner rotation matrix elements.					*
**								 	**
*************************************************************************/

#ifndef   GWigner_cc_				// Is this file already included?
#  define GWigner_cc_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#  endif

#include <Level1/Wigner.h>			// Include the interface
#include <Basics/Gconstants.h>			// Include DEG2RAD constant
#include <stdlib.h>

using std::cout;				// Using libstdc++ standard output
// ____________________________________________________________________________
// i                  WIGNER FUNCTION ERROR HANDLING
// ____________________________________________________________________________


void Wigner_error(int i)
  {
  cout << "\nWigner Rotation: ";
  switch (i)
    {
    case 0:
      cout <<                     "         (1/2)"
	   << "\nSpatial Function: Unknown d"
           << "\nSpatial Function:          ";
      break;
    case 1:
      cout <<                     "         (1)"
	   << "\nSpatial Function: Unknown d"
           << "\nSpatial Function:          ";
      break;
    case 2:
      cout <<                     "         (2)"
           << "\nSpatial Function: Unknown d"
           << "\nSpatial Function:          ";
    case 3:
      cout <<                     "         (0)"
           << "\nSpatial Function: Unknown d"
           << "\nSpatial Function:          ";
    case 10:
      cout << "Unable to Determine Reduced Wigner Element.\n";
      break;
    case 11:
      cout << "Unable to Determine Normalized Spherical Harmonic.\n";
      break;
    default:
      cout << "Unknown error.\n";
      break;
    }
  }

void volatile Wigner_fatality (int error)
  {
  Wigner_error (error);
  cout << "\nWigner Rotation: Program Aborting.\n";
  exit(-1);
  }

// ______________________________________________________________________
//           REDUCED WIGNER ROTATION MATRIX ELEMENTS: RANK 0 
// ______________________________________________________________________


double d0() { return 1.0; }

	// Input		none  :
	// Output		r     : rank 0 reduced Wigner rotation
	//				matrix element
 

// ______________________________________________________________________
//          REDUCED WIGNER ROTATION MATRIX ELEMENTS: RANK 1/2 
// ______________________________________________________________________


double d1half(int n, double beta )

	// Input		n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 1/2, m = 1/2 reduced Wigner
	//			        rotation matrix element
	// Note 		      : Euler angle input in degrees
 
  {
  double r = 0;
  beta = beta*PI/360.;
    switch (n)
      {
      case -1:				//  (1/2)
        r = -sin(beta);			// d        (beta) = - sin(beta/2)
        break;				//  1/2,-1/2

      case 1:				//  (1/2)
        r = cos(beta);			// d       (beta) = cos(beta/2)
        break;				//  1/2,1/2

      default:
        Wigner_error(0);
        cout << "1/2," << n << "/2";
        Wigner_fatality(10);
        break;
      }
  return r;
  } 


double dm1half(int n, double beta )

	// Input		n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 1/2, m = -1/2 reduced Wigner
	//				rotation matrix element
	// Note 		      : Euler angle input in degrees
 
  {
  double r=0;
  beta = beta*PI/360.;
    switch (n)
      {
      case -1:				//  (1/2)
        r = cos(beta);			// d         (beta) = cos(beta/2)
        break;				//  -1/2,-1/2

      case 1:				//  (1/2)
        r = sin(beta);			// d        (beta) = sin(beta/2)
        break;				//  -1/2,1/2

      default:
        Wigner_error(0);
        cout << "-1/2," << n << "/2";
        Wigner_fatality(10);
        break;
      }
  return r;
  } 


double d1half(int m, int n, double beta )

	// Input		m     : momentum index
	// 			n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 1/2 reduced Wigner
	//			        rotation matrix elements
	// Note 		      : Euler angle input in degrees
 
{
  double r = 0;
    switch (m)
      {
      case -1:
        r = dm1half(n, beta);
        break;
      case 1:
        r = d1half(n, beta);
        break;
      default:
        Wigner_error(0);
        cout << m << "/2," << n << "/2";
        Wigner_fatality(10);
        break;
      }
  return r;
} 



// ______________________________________________________________________
//           REDUCED WIGNER ROTATION MATRIX ELEMENTS: RANK 1 
// ______________________________________________________________________


double d11(int n, double beta )

	// Input		n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 1, m = 1 reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
 
{
  double r = 0;
  beta = beta*PI/180.;
    switch (n)
      {
      case -1:				//  (1)             2
        r = sin(beta/2.0);		// d    (beta) = sin (beta/2)
        r *= r;				//  1,-1
        break;

      case 0:				//  (1)		 [1]1/2
        r = -sqrt(0.5)*sin(beta);	// d   (beta) = -|-|   sin (beta)
        break;				//  1,0          [2]

      case 1:				//  (1)            2
        r = cos(beta/2.0);		// d   (beta) = cos (beta/2)
        r *= r;				//  1,1
        break;

      default:
        Wigner_error(1);
        cout << "1," << n;
        Wigner_fatality(10);
        break;
      }
  return r;
} 


double d10(int n, double beta )

	// Input		n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 1, m = 0 reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
 
  {
  double r=0;
    switch (n)
      {
      case -1:				//  (1)           (1)
        r = d11(0, beta);		// d    (beta) = d   (beta)
        break;				//  0,-1          1,0

      case 0:				//  (1)
        r = cos(beta*PI/180.);		// d   (beta) = cos(beta)
        break;				//  0,0

      case 1:				//  (1)            (1)
        r = -d11(0, beta);		// d   (beta) = - d   (beta)
        break;				//  0,1            1,0

      default:
        Wigner_error(1);
        cout << "0," << n;
        Wigner_fatality(10);
        break;
      }
  return r;
  }


double d1m1(int n, double beta )

	// Input		n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 1, m = -1 reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
 
  {
  double r=0;
    switch (n)
      {
      case -1:				//  (1)            (1)
        r = d11(1,beta);		// d     (beta) = d   (beta)
        break;				//  -1,-1          1,1

      case 0:				//  (1)		    (1)
        r = -d11(0,beta);		// d    (beta) = - d   (beta)
        break;				//  -1,0            1,0

      case 1:				//  (1)           (1)
        r = d11(-1,beta);		// d    (beta) = d    (beta)
        break;				//  -1,1          1,-1 

      default:
        Wigner_error(1);
        cout << "-1," << n;
        Wigner_fatality(10);
        break;
      }
  return r;
  }


double d1(int m, int n, double beta )

	// Input		beta  : Euler angle
	// 			m     : momentum index
	// 			n     : momentum index
	// Output		r     : rank 1 reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
 
  {
  double r = 0;
  switch (m)
    {
    case -1: r = d1m1(n, beta ); break;
    case 0:  r = d10(n, beta );  break;
    case 1:  r = d11(n, beta );  break;
    default:
      Wigner_error(1);
      cout << m << "," << n;
      Wigner_fatality(10);
      break;
    }
  return r;
  } 


// ______________________________________________________________________
//           REDUCED WIGNER ROTATION MATRIX ELEMENTS: RANK 2 
// ______________________________________________________________________


double d22(int n, double beta )

	// Input		n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 2, m = 2 reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
 
  {
  double r=0;
    switch (n)
      {
      case -2:
        r = sin(beta*PI/360.0);		//  (2)             4
        r = r*r;			// d    (beta) = sin (beta/2)
        r = r*r;			//  2,-2
        break;

      case -1:
        r = beta*DEG2RAD;		//  (2)          1
        r = 0.5*sin(r)*(cos(r)-1.0);	// d    (beta) = - sin(beta)[cos(beta)-1]
        break;				//  2,-1         2

      case 0:
        r = sin(beta*DEG2RAD);		//  (2)         [3]1/2   2
        r = sqrt(3.0/8.0)*r*r;		// d   (beta) = |-|   sin (beta)
        break;				//  2,0         [8]

      case 1:
        r = beta*DEG2RAD;		//  (2)           1
        r = -0.5*sin(r)*(1.0+cos(r));	// d   (beta) = - - sin(beta)[1+cos(beta)]
        break;				//  2,1           2

      case 2:
        r = cos(beta*PI/360.0);		//  (2)            4
        r = r*r;			// d   (beta) = cos (beta/2)
        r = r*r;			//  2,2
        break;

      default:
        Wigner_error(2);
        cout << "2," << n;
        Wigner_fatality(10);
        break;
      }
  return r;
  }


double d21(int n, double beta )

	// Input		n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 2, m = 1 reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
 
  {
  double r=0;
    switch (n)
      {
      case -2:				//  (2)           (2)
        r = d22(-1, beta);		// d    (beta) = d    (beta)
        break;				//  1,-2          2,-1

      case -1:
        r = cos(beta*PI/180.);		//  (2)                    1
        r = (r+0.5)*(1.0-r);		// d    (beta) = [cos(beta+-)][1-cos(beta)]
        break;				//  1,-1                   2

      case 0:
	r = beta*PI/180.;		//  (2)           [3]1/2
        r = -sqrt(1.5)*sin(r)*cos(r); 	// d   (beta) = - |-|   sin(beta)cos(beta)
        break;				//  1,0           [2]

      case 1:
        r = cos(beta*PI/180.);		//  (2)                   1
        r = (r-0.5)*(1+r);		// d   (beta) = [cos(beta--)][1+cos(beta)]
        break;				//  1,1                   2

      case 2:				//  (2)            (2)
        r = - d22(1, beta);		// d   (beta) = - d   (beta)
        break;				//  1,2            2,1

      default:
        Wigner_error(2);
        cout << "1," << n;
        Wigner_fatality(10);
        break;
      }
  return r;
  }


double d20(int n, double beta )

	// Input		n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 2, m = 0 reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
 
  {
  double r=0;
    switch (n)
      {
      case -2:				//  (2)           (2)
        r = d22(0, beta);		// d    (beta) = d   (beta)
        break;				//  0,-2          2,0

      case -1:				//  (2)           (2)
        r = d21(0, beta);		// d    (beta) = d   (beta)
        break;				//  0,-1          1,0

      case 0:
        r = cos(beta*PI/180.);		//  (2)         1      2
        r = (1.5*r*r - 0.5); 		// d   (beta) = - [3cos (beta) - 1]
        break;				//  0,0         2  

      case 1:				//  (2)            (2)
        r = - d21(0, beta);		// d   (beta) = - d   (beta)
        break;				//  0,1            1,0

      case 2:				//  (2)          (2)
        r = d22(0, beta);		// d   (beta) = d   (beta)
        break;				//  0,2          2,0

      default:
        Wigner_error(2);
        cout << "0," << n;
        Wigner_fatality(10);
        break;
      }
  return r;
  }


double d2m1(int n, double beta )

	// Input		n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 2, m = -1 reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
 
  {
  double r=0;
    switch (n)
      {
      case -2:				//  (2)            (2)
        r = d22(1, beta);		// d     (beta) = d   (beta)
        break;				//  -1,-2          2,1

      case -1:				//  (2)            (2)
        r = d21(1, beta);		// d     (beta) = d   (beta)
        break;				//  -1,-1          1,1

      case 0:				//  (2)             (2)
        r = - d21(0, beta);		// d    (beta) = - d   (beta)
        break;				//  -1,0            1,0

      case 1:				//  (2)           (2)
        r = d21(-1, beta);		// d    (beta) = d   (beta)
        break;				//  -1,1          1,-1

      case 2:				//  (2)             (2)
        r = - d22(-1, beta);		// d    (beta) = - d   (beta)
        break;				//  -1,2            2,-1

      default:
        Wigner_error(2);
        cout << "-1," << n;
        Wigner_fatality(10);
        break;
      }
  return r;
  }


double d2m2(int n, double beta )

	// Input		n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 2, m = -2 reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
 
  {
  double r=0;
    switch (n)
      {
      case -2:				//  (2)            (2)
        r = d22(2, beta);		// d     (beta) = d   (beta)
        break;				//  -2,-2          2,2

      case -1:				//  (2)              (2)
	r = -d22(1, beta);		// d     (beta) = - d   (beta)
        break;				//  -2,-1            2,1

      case 0:				//  (2)           (2)
        r = d22(0, beta);		// d    (beta) = d   (beta)
        break;				//  -2,0          2,0

      case 1:				//  (2)             (2)
        r = -d22(-1, beta);		// d    (beta) = - d   (beta)
        break;				//  -2,1            2,-1

      case 2:				//  (2)           (2)
        r = d22(-2, beta);		// d    (beta) = d    (beta)
        break;				//  -2,2          2,-2

      default:
        Wigner_error(2);
        cout << "-2," << n;
        Wigner_fatality(10);
        break;
      }
  return r;
  } 


double d2(int m, int n, double beta )

	// Input		beta  : Euler angle
	// 			m     : momentum index
	// 			n     : momentum index
	// Output		r     : rank 2 reduced Wigner rotation
	//				matrix element
 
{
  double r = 0;
     switch (m)
       {
       case -2:
         r = d2m2(n, beta);
         break;
       case -1:
         r = d2m1(n, beta);
         break;
       case 0:
         r = d20(n, beta);
         break;
       case 1:
         r = d21(n, beta);
         break;
       case 2:
         r = d22(n, beta);
         break;
       default:
         Wigner_error(2);
         cout << m << "," << n;
         Wigner_fatality(10);
         break;
       }
  return r;
} 


// ______________________________________________________________________
//           REDUCED WIGNER ROTATION MATRIX ELEMENTS: RANK J 
// ______________________________________________________________________


double fact(int a)

	// Input		a : an integer
	// Output		r : a double value of a!
	// Note			  : return is double precision to avoid
	//			    overflow at about 13! with int

  {
  double r = 1.0;
  while (a>0)
    {
    r *= a;
    a--;
    }
  return r;
  }


double dJ_int(int J, int m, int n, double beta )

	// Input		J     : Rank
	// 			m     : Momentum index
	// 			n     : Momentum index
	// 			beta  : Euler angle (degrees)
	// Output		r     : Reduced Wigner rot. matrix element
	// Note 		      : Euler angle input in degrees
	// Note 		      : Only INTEGRAL NON-NEGATIVE J
 
//  J                                     1/2
// d   (beta) = [(J+n)!(J-n)!(J+m)!(J-m)!]
//  m,n
//                           k                        2J+n-m-2k              m-n+2k
//        ---            (-1)            [           ]         [            ]
//     x  \   __________________________ |cos(beta/2)|         |-sin(beta/2)|
//        /   (J-m-k)!(J+n-k)!(k+m-n)!k! [           ]         [            ]
//        ---
//         k

// Where the summation index k spans all integers where the factorial arguments
// are non-negative (See Zare page 86, Equation 3.57).  Also note that these
// must adhere to several symmetry rules and that's a good check of validity.

//    J             J               (m-n) J             J
//   d   (beta) =  d     (beta) = -1     d   (-beta) = d   (-beta)
//    m,n           -n,-m                 m,n           n,m

// In this function, we'll switch the nomenclature by using the following variables:

//                                                                  1/2
// Jpn=J+n; Jmn=J-n; Jpm=J+m; Jmm=J-n; prefact=[Jpm!*Jmm!*Jpn!*Jmn!]
// d1=J-m-k=Jmm-k; d2=J+n-k=Jpn-k; d3=k+m-n; denom = d1!*d2!*d3!*k!
// COSB = cos(beta/2);  mSINB = -sin(beta/2);
// cospow = 2J+n-m-2k = Jpn+Jmm-2k; sinpow=m-n+2k

// So that the simple formula below can be used

//                               k
//  J                    --- (-1)      cospow      sinpow
// d   (beta) = prefact  \   _____ COSB       mSINB
//  m,n                  /   denom 
//                       ---        
//                        k
  {
  if((J < 0) || (abs(m) > J) || (abs(n) >J))		// 1st insure resaonable {J,m,n}  
    {
    cout << "\n                           (" << J << ")"
         << "\nSpatial Function: Unknown d"
         << "\n                           " << m << "," << n;
    Wigner_fatality(10);
    }
  double betaover2 = beta*PI/360.0;			// Get beta/2 in radians
  double COSF, COSB = cos(betaover2);			// Cosine factors needed
  double SINF, mSINB = -sin(betaover2);			// Sine factors needed
  double djval = 0.0;					// Value of Wigner element
  int Jpm=J+m, Jmm=J-m, Jpn=J+n, Jmn=J-n;		// We need these
  double prefact =					// Constant prefactor
        sqrt(fact(Jpm)*fact(Jmm)*fact(Jpn)*fact(Jmn));
  double signf=1.0, denom=0.0;
  int d1, d2, d3;
  for(int k=0; k<=2*J; k++)
    {
    d1 = Jmm-k;
    d2 = Jpn-k;
    d3 = k+m-n;
    if((d1>=0) && (d2>=0) && (d3>=0))
      {
      COSF = pow(COSB, (Jpn + Jmm - 2*k));		// Cosine to power
      SINF = pow(mSINB, (2*k + m - n));			// Sine to power
      denom = fact(d1)*fact(d2)*fact(d3)*fact(k);		// Get demoninator
      djval += signf*COSF*SINF/denom;			// Unscaled dJ contribution
      }
    signf *= -1.0;					// Switch sign 1 <-> -1
    }
  return djval*prefact;
  }


double dJ_half_int(int J, int m, int n, double beta )

	// Input		J     : rank
	// 			m     : momentum index
	// 			n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
	// Note 		      : Here J implies 1/2 units
	//				J=1->1/2; m,n=1->1/2; m,n=-1->-1/2
	//				J=2->3/2; m,n=2->3/2; m,n=1->1/2
	//				          m,n=-1->-1/2; m,n=-2->-3/2
	//				In general, J=a->J=a-1/2 and m,n then
	//				span [a-1/2, -(a-1/2)] incremented by 1
 
{
  beta = beta*PI/360.0;
  double r = 0.0;
  double COS = 0.0;
  double SIN = 0.0;
  double num = 0.0;
  double denom = 0.0;
  int t=0;   
  int actJ = 2*J - 1;
  int actm = 2*abs(m) - 1;
  int actn = 2*abs(n) - 1;
    if(m < 0)
      actm = -actm;
    if(n < 0)
      actn = -actn;
//  int addJm = J+m;
  if((J > 0) && (actm <= actJ) && (actn <= actJ) && (m != 0) && (n != 0) )
    {
    while((t <= (actJ+actm)/2) && (t <= (actJ-actn)/2))
      {
      if((t + (actn-actm)/2) >= 0)
        {
        COS = pow(cos(beta), (actJ + (actm - actn)/2 - 2*t));
        SIN = pow(sin(beta), (2*t + (actn - actm)/2));
        num = fact((actJ+actm)/2) * fact((actJ-actm)/2);
        num *= fact((actJ+actn)/2) * fact((actJ-actn)/2);
        denom = fact((actJ + actm)/2 - t) * fact((actJ - actn)/2 - t);
        denom *= fact(t) * fact(t + (actn - actm)/2);
        r += pow(-1.0,t)*sqrt(num)*COS*SIN/denom;
        }
      t++;
      }
    }
  else
    {
      cout << "\nSpatial Function:          (" << 2*J-1 << "/2)"
           << "\nSpatial Function: Unknown d"
           << "\nSpatial Function:          " << m << "/2," << n << "/2";
      Wigner_fatality(10);
    }
  return r;
}


double dJ(int J, int m, int n, double beta )

	// Input		J     : rank
	// 			m     : momentum index
	// 			n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
	// Note 		      : Negative J implies 1/2 units
	//				J=-1->1/2; m,n=1->1/2; m,n=-1->-1/2
	//				J=-2->3/2; m,n=2->3/2; m,n=1->1/2
	//				           m,n=-1->-1/2; m,n=-2->-3/2
	//				In general, J=-a->J=a-1/2 and m,n then
	//				span [a-1/2, -(a-1/2)] incremented by 1
 
{
  double r = 1.;
    switch (J)
      {
      case -1:
        r = d1half(m, n, beta);	
        break;
      case 1:
        r = d1(m, n, beta);
        break;
      case 2:
        r = d2(m, n, beta);
        break;
      default:
        if(J > 0)
          r = dJ_int(J, m, n, beta);
        else if(J < 0)
          r = dJ_half_int(-J, m, n, beta);
	else
          {
          if((m != 0) || (n != 0))
            {
            Wigner_error(3);
            cout << m << "," << n;
            Wigner_fatality(10);
            }
          }
        break;
      }
   return r;
}


// ______________________________________________________________________
//                     REDUCED WIGNER ROTATION MATRIX
// ______________________________________________________________________


matrix dJ(int J, double beta)

	// Input		J     : rank
	// 			beta  : Euler angle
	// Output		mx    : rank J reduced Wigner rotation matrix
	// Note 		      : Euler angle beta input in degrees


{
  int dim = 2*abs(J)+1;
  int m = 0;
  int n = 0;
  int i = 0;
  int j = 0;
  if(J < 0) dim--;
  matrix d(dim,dim);
    if(J < 0)
      {
      for(i=0; i<dim; i++)
        {
        m = dim/2 - i;
        if(m <= 0) m--;
        for(j=0; j<dim; j++)
          {
          n = dim/2 - j;
          if(n <= 0) n--;
          d.put(dJ(J, m, n, beta), i,j);
          }
        }
      }
    else
      {
      for(i=0; i<dim; i++)
        {
	m = J - i;
        for(j=0; j<dim; j++)
          {
          n = J - j;
          d.put(dJ(J, m, n, beta), i,j);
          }
        }
      }
  return d;
}



// ______________________________________________________________________
//                     WIGNER ROTATION MATRIX ELEMENTS
// ______________________________________________________________________


complex DJ(int J, int m, int n, double alpha, double beta, double gamma)

	// Input		J     : rank
	// 			m     : momentum index
	// 			n     : momentum index
	// 			alpha : Euler angle
	// 			beta  : Euler angle
	// 			gamma : Euler angle
	// Output		z     : rank 2 Wigner rotation matrix element
	// Note 		      : Euler angles input in degrees
	// Note 		      : Negative J implies 1/2 units
	//				J=-1->1/2; m,n=1->1/2; m,n=-1->-1/2
	//				J=-2->3/2; m,n=2->3/2; m,n=1->1/2
	//				           m,n=-1->-1/2; m,n=-2->-3/2
	//				In general, J=-a->J=a-1/2 and m,n then
	//				span [a-1/2, -(a-1/2)] incremented by 1

  {
  double r;
  complex z;
  double actm = double(m);
  double actn = double(n);
  if(J<0)
    {
    if(m) actm = double(m)-0.5;
    else actm = double(m)+0.5;
    if(n) actn = double(n)-0.5;
    else actn = double(n)+0.5;
    }
  r = (alpha*actm + gamma*actn)*PI/180.;
  z = exp(complex(0.,-r));
  z *= dJ(J, m, n, beta);		// d takes beta in degrees also
  return z;
  } 


double D0()

	// Input		none  :
	// Output		r     : rank 0 Wigner rotation matrix element

  { return d0(); }


complex D1half(int m, int n, double alpha, double beta, double gamma)

	// Input		m     : momentum index
	// 			n     : momentum index
	// 			alpha : Euler angle
	// 			beta  : Euler angle
	// 			gamma : Euler angle
	// Output		z     : rank 1/2 Wigner rotation matrix element
	// Note 		      : Euler angles input in degrees

  { return DJ(-1,m,n,alpha,beta,gamma); }


complex D1(int m, int n, double alpha, double beta, double gamma)

	// Input		m     : momentum index
	// 			n     : momentum index
	// 			alpha : Euler angle
	// 			beta  : Euler angle
	// 			gamma : Euler angle
	// Output		z     : rank 1 Wigner rotation matrix element
	// Note 		      : Euler angles input in degrees

  { return DJ(1,m,n,alpha,beta,gamma); }


complex D2(int m, int n, double alpha, double beta, double gamma)

	// Input		m     : momentum index
	// 			n     : momentum index
	// 			alpha : Euler angle
	// 			beta  : Euler angle
	// 			gamma : Euler angle
	// Output		z     : rank 2 Wigner rotation matrix element
	// Note 		      : Euler angles input in degrees

  { return DJ(2,m,n,alpha,beta,gamma); } 


// ______________________________________________________________________
//                     WIGNER ROTATION MATRICES
// ______________________________________________________________________


  matrix DJ(int J, double alpha, double beta, double gamma)

	// Input		J     : rank
	// 			alpha : Euler angle
	// 			beta  : Euler angle
	// 			gamma : Euler angle
	// Output		mx    : rank J Wigner rotation matrix
	// Note 		      : Euler angles input in degrees


  {
  int dim = 2*abs(J)+1;
  int m = 0;
  int n = 0;
  int i = 0;
  int j = 0;
  double r = 0;
  complex z;
  alpha = alpha*PI/180.;
  gamma = gamma*PI/180.;
//  int half = 1;
  if(J < 0) dim--;
  matrix D(dim,dim);
  D = dJ(J, beta);
    if(J < 0)
      {
      for(i=0; i<dim; i++)
        {
        m = dim/2 - i;
        if(m <= 0) m--;
        for(j=0; j<dim; j++)
          {
          n = dim/2 - j;
          if(n <= 0) n--;
          r = (alpha*m + gamma*n)/2.0;
          z = exp(complex(0.,-r));
          D.put(D(i,j)*z, i,j);
          }
        }
      }
    else
      {
      for(i=0; i<dim; i++)
        {
	m = J - i;
        for(j=0; j<dim; j++)
          {
          n = J - j;
          r = alpha*m + gamma*n;
          z = exp(complex(0.,-r));
          D.put(D(i,j)*z, i,j);
          }
        }
      }
  return D;
  }


matrix DJ(const matrix& dJbeta, int J, double alpha)

	// Input		dj(beta) : rank J Wigner array for beta
	//			J	 : J rank
	// 			alpha    : Euler angle (degrees)
	// Output		mx    : rank J Wigner rotation matrix

  {
  int dim = 2*abs(J)+1;
  int m = 0;
  int n = 0;
  int i = 0;
  int j = 0;
  complex z;
  alpha = alpha*PI/180.;
//  int half = 1;
  if(J < 0) dim--;
  matrix D = dJbeta;
    if(J < 0)
      {
      for(i=0; i<dim; i++)
        {
        m = dim/2 - i;
        if(m <= 0) m--;
        for(j=0; j<dim; j++)
          {
          n = dim/2 - j;
          if(n <= 0) n--;
          z = exp(complex(0.,-(alpha*m)/2));
          D.put(D(i,j)*z, i,j);
          }
        }
      }
    else
      {
      for(i=0; i<dim; i++)
        {
	m = J - i;
        for(j=0; j<dim; j++)
          {
          n = J - j;
          z = exp(complex(0.,-alpha*m));
          D.put(D(i,j)*z, i,j);
          }
        }
      }
  return D;
  }

#endif						// Wigner.cc
