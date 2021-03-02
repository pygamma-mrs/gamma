/* SpinOpSng.cc *****************************************************-*-c++-*-
**							                    **
**                                  G A M M A				    **
**									    **
**      Single Spin Operators				Implementation      **
**									    **
**      Copyright (c) 1990, 1991, 1992					    **
**      Tilo Levante, Scott A. Smith					    **
**      Eidgenoessische Technische Hochschule				    **
**      Labor fuer physikalische Chemie					    **
**      8092 Zuerich / Switzerland					    **
**									    **
**      $Header: $
**									   **
*****************************************************************************/

/*****************************************************************************
**									    **
**  Description								    **
**									    **
** This module contains functions which define the common single spin	    **
** operators for GAMMA.  These are returned as matrices, represented in the **
** common basis: alpha = (1,0,...,0) & beta = (0,1,0,...,0)		    ** 
**									    **
** Currently provided are the following:				    **
**									    **
**		        Ix, Iy, Iz, Ie, Im, Ip, Ru(phi)			    **
**									    **
*****************************************************************************/

#ifndef   single_spin_op_cc_		// Is file already included?
#  define single_spin_op_cc_ 1		// If no, then remember that it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Matrix/matrix.h>		// Include GAMMA matrices
#include <HSLib/SpinOpSng.h>		// Include the implementation
#include <list>				// Include STL lists
#include <map>				// Include STL maps 

using std::map;
using std::list;

/* All of these functions return a common spin operator that relates to a
   SINGLE spin.  The operator returned is NOT spin isotope type specific,
   rather it is spin Hilbert space specific.  The user inputs a particular
   Hilbert space and the requested operator is returned. Spin Hilbert spaces
   are related to the spin I value according to HS = 2*I+1.  Thus, for spin
   I=1/2 species (1H, 13C, 15N, ....) the Hilbert space dimension will be 2,
   whereas spin 1 species (2H) will have a Hilbert space dimension of 2, and
   so on.

   In additon, these functions also store a list of known single spin operators.
   That is, every GAMMA program that runs will store all single spin operators
   used throughout the simulation.  Since these are small arrays, the overhead
   is small.  In compensation, all such operators will be readily available for
   calculations and they are very often used to generate composite spin 
   operators, spin Hamiltonians, spin tensors, pulses, etc.  

	Input		hs  : The spin Hilbert space, (2I+1)
	Output		mx  : Matrix representation of SOp
				Ie - identity operator  (IeList)

                        I:   1/2   1   3/2   2   5/2
 			HS:   2    3    4    5    6                           */

matrix Ie(int hs)
  {
  static std::map<int, matrix> IeList;	// Map of HS vs Ie arrays
  if(hs <1) return matrix(); 		// Return NULL if I < 1/2
  std::map<int,matrix>::iterator i;	// An iterator into the list
  i = IeList.find(hs);			// See if Ie(hs) in list
  if(i == IeList.end())			// If it isn't we'll make one 
    {					// and add it to the list
    matrix IeMx(hs,hs,i_matrix_type);
    IeList.insert(std::pair<int,matrix>(hs, IeMx));
    i = IeList.find(hs);
    }
  return i->second;
  }


matrix Ix(int hs)
  {
  static std::map<int, matrix> IxList;	// Map of HS vs Ix arrays
  if(hs <1) return matrix(); 		// Return NULL if I < 1/2
  std::map<int,matrix>::iterator i;	// An iterator into the list
  i = IxList.find(hs);			// See if Ix(hs) in list
  if(i == IxList.end())			// If it isn't we'll make one 
    {					// and add it to the list
    matrix IxMx(hs,hs,0,h_matrix_type);	// Return is Hermitian
    double q  = (hs-1.0)/2;		// This is spin I value
    double q1 = q*(q+1);
    double m, tmp;
    int k;
    for(k=0,m=q-1; k<hs-1; k++,m-=1)
      {
      tmp = sqrt(q1-m*(m+1))/2;
      IxMx.put_h(tmp,k+1,k);
      }
    IxList.insert(std::pair<int,matrix>(hs, IxMx));
    i = IxList.find(hs);
    }
  return i->second;
  }


matrix Iy(int hs)
  {
  static std::map<int, matrix> IyList;	// Map of HS vs Iy arrays
  if(hs <1) return matrix(); 		// Return NULL if I < 1/2
  std::map<int,matrix>::iterator i;	// An iterator into the list
  i = IyList.find(hs);			// See if Iy(hs) in list
  if(i == IyList.end())			// If it isn't we'll make one 
    {					// and add it to the list
    matrix IyMx(hs,hs,0,h_matrix_type);	// Return is Hermitian
    double q  = (hs-1.0)/2;		// This is spin I value
    double q1 = q*(q+1);
    double m, tmp;
    int k;
    for(k=0,m=q-1; k<hs-1; k++,m-=1)
      {
      tmp = sqrt(q1-m*(m+1))/2;
      IyMx.put_h(complex(0,tmp),k+1,k);
      }
    IyList.insert(std::pair<int,matrix>(hs, IyMx));
    i = IyList.find(hs);
    }
  return i->second;
  }


matrix Iz(int hs)
  {
  static std::map<int, matrix> IzList;	// Map of HS vs Iz arrays
  if(hs <1) return matrix(); 		// Return NULL if I < 1/2
  std::map<int,matrix>::iterator i;	// An iterator into the list
  i = IzList.find(hs);			// See if Iz(hs) in list
  if(i == IzList.end())			// If it isn't we'll make one 
    {					// and add it to the list
    matrix IzMx(hs,hs,0,d_matrix_type);	// Return is diagonal
    double m, q=(hs-1.0)/2;		// This is spin I value
    int k;
    for(k=0,m=q; k<hs; k++,m-=1)
      IzMx.put(m,k,k);
    IzList.insert(std::pair<int,matrix>(hs, IzMx));
    i = IzList.find(hs);
    }
  return i->second;
  }


matrix Ip(int hs)
  {
  static std::map<int, matrix> IpList;	// Map of HS vs I+ arrays
  if(hs <1) return matrix(); 		// Return NULL if I < 1/2
  std::map<int,matrix>::iterator i;	// An iterator into the list
  i = IpList.find(hs);			// See if I+(hs) in list
  if(i == IpList.end())			// If it isn't we'll make one 
    {					// and add it to the list
    matrix IpMx(hs,hs,0,n_matrix_type);	// Return is normal
    double q  = (hs-1.0)/2;
    double q1 = q*(q+1);
    double m;
    int k;
    for(k=0,m=q-1; k<hs-1; k++,m-=1)
      IpMx.put(sqrt(q1-m*(m+1)),k,k+1);
    IpList.insert(std::pair<int,matrix>(hs, IpMx));
    i = IpList.find(hs);
    }
  return i->second;
  }


matrix Im(int hs)
  {
  static std::map<int, matrix> ImList;	// Map of HS vs I- arrays
  if(hs <1) return matrix(); 		// Return NULL if I < 1/2
  std::map<int,matrix>::iterator i;	// An iterator into the list
  i = ImList.find(hs);			// See if I-(hs) in list
  if(i == ImList.end())			// If it isn't we'll make one 
    {
    matrix ImMx(hs,hs,0,n_matrix_type);	// Return is normal
    double q  = (hs-1.0)/2;
    double q1 = q*(q+1);
    double m;
    int k;
    for(k=0,m=q-1; k<hs-1; k++,m-=1)
      ImMx.put(sqrt(q1-m*(m+1)),k+1,k);
    ImList.insert(std::pair<int,matrix>(hs, ImMx));
    i = ImList.find(hs);
    }
  return i->second;
  }


matrix Raxis(int hs, double beta, char axis)

  //  Input          hs	: The spin Hilbert space, (2I+1).
  //		   beta	: A rotation angle (degrees)
  //		   axis	: A rotation axis; x, y, or z
  //  Output    R(beta) : A rotation operator for spin
  //			  rotation about the defined axis
  //			  by angle beta.
  //  Note		: I=1/2 -> hs = 2, I=1 -> hs = 3
  //  Note		: Unlike other functions in this module,
  //			  there is no list of rotation operators
  //			  maintained herein.  That is because
  //			  these lists are ordered by hs and the
  //			  rotations must also be ordered by angle
  //			  and axis.
  // sosi		  The return matrix must somehow be
  // 			  set to unitary!!

/*		               [                                            ]
                       I=1/2   | cos(B/2)-iu sin(B/2)  (-iu -u )sin(B/2)    |
 R (B) = exp(-i*B*I ) -------> |            z              x  y             |
  u                u	       |                                            |
  			       | (-iu +u )sin(B/2)     cos(B/2)+iu sin(B/2) |
  			       [     x  y                         z         ]*/

  {
/*           We'll Store Potentially Reused Single Spin Rotation Matrices    */

  static std::map<int, matrix> RdxList;	// Map of HS vs Rx eigenvalue arrays
  static std::map<int, matrix> RdyList;	// Map of HS vs Ry eigenvalue arrays
  static std::map<int, matrix> RdzList;	// Map of HS vs Rz eigenvalue arrays
  static std::map<int, matrix> RvxList;	// Map of HS vs Rx eigenvector arrays
  static std::map<int, matrix> RvyList;	// Map of HS vs Ry eigenvector arrays
  static std::map<int, matrix> RvzList;	// Map of HS vs Rz eigenvector arrays
  double rad = beta*PI/180.0;		// Angle is beta in radians
  matrix R(hs,hs);			// Rotation operator

/*           For Idiotic Spin Hilbert Spaces We'll Handle Things Directly    */

  if(hs <1) return matrix(); 		// Return NULL if I < 1/2
  if(hs == 1)				// Return I if I = 0;
    return matrix(1,1,i_matrix_type);
  if(beta == 0)				// Return I if beta=0
    return matrix(hs,hs,i_matrix_type);

/*           For I=1/2 We'll Just Use Formulae To Fill In The Array          */

  if(hs == 2)				// Spin I=1/2 : Ru = 1*cb - 2i*Iu*sb
    {
    double cb = cos(rad/2);		// Angle is beta/2 in radians
    double sb = sin(rad/2);		// Angle is beta/2 in radians
    switch(axis)
      {
      case 'x': 			// Rx(b) = cb - 2i*Ix*sb
        R.put(cb,0,0);			// First include cosine elements
        R.put(cb,1,1);			// 1*cb with cb = cosine(beta/2)
	R.put(complex(0,-sb),0,1); 	// Now put in -2i*Ix*sb
	R.put(complex(0,-sb),1,0);
	break;
      case 'y':				// Ry(b) = cb - 2i*Iy*sb
        R.put(cb,0,0);			// First include cosine elements
        R.put(cb,1,1);			// 1*cb with cb = cosine(beta/2)
	R.put(-sb,0,1);			// Now put in -2i*Iy*sb
	R.put(sb,1,0);
	break;
      case 'z':				// Rz(b) = cb - 2i*Iz*sb
      default:
        R = matrix(hs,hs,d_matrix_type);
	R.put(complex(cb,-sb),0,0);	// These are cosine and sine terms
	R.put(complex(cb,sb),1,1);
	break;
      }
    return R;
    }

/*         For I>1/2 We'll Take The Exponential Of The Spin Operator
           & For Ru We Use The Eigensystem of Iu Which We May Store          */

  std::map<int,matrix>::iterator i;	// An iterator into the Rd list
  std::map<int,matrix>::iterator j;	// An iterator into the Rv list
  matrix Rd, Rv;			// R diagonal, R eigenvectors
  switch(axis)				// See if we have one on tap
    {
    case 'x':
      i = RdxList.find(hs);		//  Try for Rx eigenvalues
      if(i!=RdxList.end())		//  If eigenvalues exist
        {				//  We'll set both Rd & Rv
        j = RvxList.find(hs);		//    Here is Pix to eigenvectors
        Rd = i->second;			//    Get R eigenvalues
        Rv = j->second;			//    Get R eigenvectors
        }
      else				//  If eigensystem isnt stored
        {				//  we must generate & store it
        R = Ix(hs); 			//    Get Ix (in R temp) 
        diag(R, Rd, Rv);		//    Diagonalize Ix
        RdxList.insert(std::pair<int,matrix>(hs, Rd));
        RvxList.insert(std::pair<int,matrix>(hs, Rv));
        }
      break;
    case 'y':
      i = RdyList.find(hs);		//  Try for Ry eigenvalues
      if(i!=RdyList.end())		//  If eigenvalues exist
        {				//  We'll set both Rd & Rv
        j = RvyList.find(hs);		//    Here is Pix to eigenvectors
        Rd = i->second;			//    Get R eigenvalues
        Rv = j->second;			//    Get R eigenvectors
        }
      else				//  If eigensystem isnt stored
        {				//  we must generate & store it
        R = Iy(hs); 			//    Get Iy (in R temp) 
        diag(R, Rd, Rv);		//    Diagonalize Iy
        RdyList.insert(std::pair<int,matrix>(hs, Rd));
        RvyList.insert(std::pair<int,matrix>(hs, Rv));
        }
      break;
    case 'z':
      i = RdzList.find(hs);		//  Try for Rz eigenvalues
      if(i!=RdzList.end())		//  If eigenvalues exist
        {				//  We'll set both Rd & Rv
        j = RvzList.find(hs);		//    Here is Pix to eigenvectors
        Rd = i->second;			//    Get R eigenvalues
        Rv = j->second;			//    Get R eigenvectors
        }
      else				//  If eigensystem isnt stored
        {				//  we must generate & store it
        R = Iz(hs); 			//    Get Iz (in R temp) 
        diag(R, Rd, Rv);		//    Diagonalize Iz
        RdzList.insert(std::pair<int,matrix>(hs, Rd));
        RvzList.insert(std::pair<int,matrix>(hs, Rv));
        }
      break;
    }
  complex fact = -1.0*rad*complexi;	//  Finally we generate Ru by taking
  for(int k=0; k<hs; k++)		//  the exponential of Iu (diagonal)
    Rd.put(exp(fact*Rd.get(k,k)),k,k);	//  then resetting into default basis
  R = Rv * times_adjoint(Rd,Rv);
  return R;
  }

#endif 							// SpinOpSng.cc

