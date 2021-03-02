/* SphHarmic.cc *************************************************-*-c++-*-
**									**
**	                        G A M M A				**
**								 	**
**	Spherical Harmonic Functions              Implementation	**
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
** This module provides spherical harmonics to GAMMA.			**
**								 	**
*************************************************************************/

#ifndef   Gsphharm_cc_				// Is this file already included?
#  define Gsphharm_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#  endif

#include <Level1/SphHarmic.h>			// Include the interface
#include <Basics/Gconstants.h>			// Include DEG2RAD constant
#include <Basics/Gutils.h>			// Include GAMMA errors
#include <Basics/StringCut.h>

using std::string;			// Using libstdc++ strings

// ____________________________________________________________________________
// i                  SPHERICAL HARMONICS ERROR HANDLING
// ____________________________________________________________________________


void Y_error(int eidx, int noret)

        // Input                eidx    : Error index
        //                      noret   : Flag for return (0=linefeed)
	// Output		none	: Error Message Output

  {
  string hdr("Spherical Harmonic");
  string msg;
  switch(eidx)
    {
    case 0: msg = "Unknown Spherical Harmonic Y";
            GAMMAerror(hdr, msg, noret); break;	//                         (0)
    case 2: msg = "Unable to Determine Normalized Spherical Harmonic";
            GAMMAerror(hdr, msg, noret); break;	//                         (2)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }

void Y_error(int eidx, const string& pname, int noret)

        // Input                eidx    : Error index
	//			pname 	: Additional error message
        //                      noret   : Flag for return (0=linefeed)
	// Output		none	: Error Message Output

  {
  string hdr("Spherical Harmonic");
  string msg;
  switch(eidx)
    {
    case 1: msg = "                            " + pname;
            GAMMAerror(hdr, msg, noret); break;	//                         (0)
    default: GAMMAerror(hdr, -1, pname, noret); break; // Unknown Error    (-1)
    }
  }

void volatile Y_fatality(int eidx)

	// Input		none :
        // Input                eidx    : Error index

  {
  Y_error(eidx, 1);			// Output error message
  if(eidx) Y_error(0);			// Now output it fatal
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
//                 SPHERICAL HARMONICS EMPLOYING RADIANS
// ____________________________________________________________________________


double Y00rad()

	// Input		none : unnecessary
	// Output		z    : rank zero normalized spherical harmonic

 				//  0                         1/2
  { return sqrt(0.25/PI); } 	// Y (theta, phi) = [1/(4*PI)]
				//  0

double Y10rad(double theta)

	// Input		theta : spherical angle
	// Output		z     : rank one normalized spherical harmonic
	// Note			      : angle theta input in radians

					//  0              [  3 ]1/2
  { return sqrt(0.75/PI)*cos(theta); }	// Y (theta,phi) = |----|   cos(theta)
  					//  1	           [4*PI]


complex Y11rad(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 1 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians

  {
  double r = sqrt(0.375/PI);	//  1               [  3 ]1/2            i*phi
  r *= sin(theta);		// Y (theta,phi) = -|----|   sin(theta) e
  complex z(0.0, phi);		//  1		    [8*PI]
  return -r*exp(z);
  } 


complex Y1m1rad(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 1 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians

  {
  double r = sqrt(0.375/PI);	//  -1              [  3 ]1/2            -i*phi
  r *= sin(theta);		// Y  (theta,phi) = |----|   sin(theta) e
  complex z(0.0, -phi);		//  1		    [8*PI]
  return r*exp(z);
  } 


double Y20rad(double theta)

	// Input		theta : spherical angle
	// Output		z     : rank 2 normalized spherical harmonic
	// Note			      : angle theta input in radians

  {
  double r1 = sqrt(0.3125/PI);	//  0              [  5  ]1/2     2
  double r2 = cos(theta);	// Y (theta,phi) = |-----|   [3cos (theta) - 1]
  r2 = (3*r2*r2)-1;		//  2		   [16*PI]
  return r1*r2;
  } 


complex Y21rad(double  theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 2 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians

  {
  double r = sqrt(1.875/PI);	//  1               [ 15 ]1/2 		  	i*phi 
  r = r*cos(theta)*sin(theta);	// Y (theta,phi) = -|----|   cos(theta)sin(theta)e
  complex z(0.0, phi);		//  2      	  [8*PI]
  return -r*exp(z);
  } 


complex Y2m1rad(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 2 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians

  {
  double r = sqrt(1.875/PI);	//  -1              [ 15 ]1/2 		  	-i*phi 
  r = r*cos(theta)*sin(theta);	// Y  (theta,phi) = |----|   cos(theta)sin(theta)e
  complex z(0.0, -phi);		//  2      	  [8*PI]
  return r*exp(z);
  } 
 

complex Y22rad(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 2 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians

  {
  double r = sqrt(0.46875/PI);	//  2              [  15 ]1/2 	2  	 2*i*phi 
  r = r*sin(theta)*sin(theta);	// Y (theta,phi) = |-----|   sin (theta)e
  complex z(0.0, 2.0*phi);	//  2       	   [32*PI]
  return r*exp(z);
  } 


complex Y2m2rad(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 2 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians

  {
  double r = sqrt(0.46875/PI);	//  2              [  15 ]1/2 	2  	 -2*i*phi 
  r = r*sin(theta)*sin(theta);	// Y (theta,phi) = |-----|   sin (theta)e
  complex z(0.0, -2.0*phi);	//  2       	   [32*PI]
  return r*exp(z);
  } 


double Y30rad(double theta)

	// Input		theta : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angle theta input in radians

  {
  double r1 = sqrt(1.75/PI);	//  0        [  7  ]1/2     3                2
  double r2 = cos(theta);	// Y (t,p) = |-----|   [2cos (t) - 3cos(t)sin (t)
  double r3 = sin(theta);	//  3        [16*PI]
  return r1*r2*(r2*r2-1.5*r3*r3);
  } 


complex Y31rad(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians
 
  {
  double r1 = -sqrt(5.25/PI);	//  1         [ 21  ]1/2     2            3    i*phi 
  double r2 = cos(theta);	// Y (t,p) = -|-----|   [4cos(t)sin(t)-sin(t)]e
  double r3 = sin(theta);	//  3         [64*PI]
  complex z(0.0, phi);
  r2 *= r2*r3;
  r3 *= 0.25*r3*r3;
  return r1*(r2-r3)*exp(z);
  } 


complex Y3m1rad(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians

  {
  double r1 = sqrt(5.25/PI);	//  -1       [ 21  ]1/2     2            3    -i*phi 
  double r2 = cos(theta);	// Y (t,p) = |-----|   [4cos(t)sin(t)-sin(t)]e
  double r3 = sin(theta);	//  3        [64*PI]
  complex z(0.0, -phi);
  r2 *= r2*r3;
  r3 *= 0.25*r3*r3;
  return r1*(r2-r3)*exp(z);
  } 
 

complex Y32rad(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians

  {
  double r1 = sqrt(3.28125/PI);	//  2        [ 105 ]1/2         2    2*i*phi 
  double r2 = cos(theta);	// Y (t,p) = |-----|   cos(t)sin (t)e
  double r3 = sin(theta);	//  3        [32*PI]
  complex z(0.0, 2.0*phi);
  return r1*r2*r3*r3*exp(z);
  } 


complex Y3m2rad(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians

  {
  double r1 = sqrt(3.28125/PI);	//  -2       [ 105 ]1/2         2    -2*i*phi 
  double r2 = cos(theta);	// Y (t,p) = |-----|   cos(t)sin (t)e
  double r3 = sin(theta);	//  3        [32*PI]
  complex z(0.0, -2.0*phi);
  return r1*r2*r3*r3*exp(z);
  } 
 

complex Y33rad(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians

  {
  double r1 = -sqrt(.546875/PI);//  3         [ 35  ]1/2   3    3*i*phi 
  double r2 = sin(theta);	// Y (t,p) = -|-----|   sin (t)e
  complex z(0.0, 3.0*phi);	//  3         [64*PI]
  return r1*r2*r2*r2*exp(z);
  } 


complex Y3m3rad(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians

  {
  double r1 = sqrt(0.546875/PI);//  -3        [ 35  ]1/2   3    -3*i*phi 
  double r2 = sin(theta);	// Y  (t,p) = |-----|   sin (t)e
  complex z(0.0, -3.0*phi);	//  3         [64*PI]
  return r1*r2*r2*r2*exp(z);
  } 


double Y40rad(double theta)

	// Input		theta : spherical angle
	// Output		z     : rank 4 normalized spherical harmonic
	// Note			      : angle theta input in radians

  {
  double r1 = 0.1875/sqrt(PI);	//  0        [  9   ]1/2      4          2
  double r2 = cos(theta);	// Y (t,p) = |------|   [35cos(t) - 30cos(t) + 3]
  r2 = r2*r2;			//  4        [256*PI]
  return r1*(35*r2*r2-30*r2+3);
  } 


complex Y41rad(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 4 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians
 
  {
  double r1 = -0.375*sqrt(5/PI);//  1         [ 45  ]1/2           3            i*phi 
  r1 *= sin(theta);		// Y (t,p) = -|-----|   sin(t)[7cos(t)-3cos(t)]e
  double r2 = cos(theta);	//  4         [64*PI]
  complex z(0.0, phi);
  return z*(r1*(7*r2*r2*r2-3*r2));
  } 


complex Y4m1rad(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 4 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians

  {
  double r1, r2;
  r1 = 0.375*sqrt(5/PI);      //  -1       [ 45  ]1/2           3  	     -i*phi 
  r1 *= sin(theta);           // Y (t,p) = |-----|   sin(t)[7cos(t)-3cos(t)]e
  r2 = cos(theta);	      //  4        [64*PI]
  complex z(0.0, -phi);
  return z*(r1*(7*r2*r2*r2-3*r2));
  } 
 

complex Y42rad(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 4 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians

{
  double r1, r2, r3;
  complex z;
  r1 = 0.375*sqrt(2.5/PI);	//  2        [  45  ]1/2   2       2      2*i*phi 
  r2 = cos(theta);		// Y (t,p) = |------|   sin(t)[7cos(t)-1]e
  r3 = sin(theta);		//  4        [128*PI]
  z = complex(0.0, 2.0*phi);
  r2 = r2*r2;
  r3 = r3*r3;
  return r1*r3*(7*r2-1)*exp(z);
} 


complex Y4m2rad(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 4 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians

{
  double r1, r2, r3;
  complex z;
  r1 = 0.375*sqrt(2.5/PI);	//  2        [  45  ]1/2   2       2      2*i*phi 
  r2 = cos(theta);		// Y (t,p) = |------|   sin(t)[7cos(t)-1]e
  r3 = sin(theta);		//  4        [128*PI]
  z = complex(0.0, -2.0*phi);
  r2 = r2*r2;
  r3 = r3*r3;
  return r1*r3*(7*r2-1)*exp(z);
} 
 

complex Y43rad(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 4 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians

{
  double r1, r2, r3;
  complex z;
  r1 = -0.375*sqrt(35/PI);      //  3         [ 315 ]1/2   3         3*i*phi 
  r2 = sin(theta);		// Y (t,p) = -|-----|   sin(t)cos(t)e
  r3 = cos(theta);		//  4         [64*PI]
  z = complex(0.0, 3.0*phi);
  return r1*r2*r2*r2*r3*exp(z);
} 


complex Y4m3rad(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 4 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians

{
  double r1, r2, r3;
  complex z;
  r1 = 0.375*sqrt(35/PI);       //  3        [ 315 ]1/2   3         -3*i*phi 
  r2 = sin(theta);		// Y (t,p) = |-----|   sin(t)cos(t)e
  r3 = cos(theta);		//  4        [64*PI]
  z = complex(0.0, 3.0*phi);
  return r1*r2*r2*r2*r3*exp(z);
} 
 

complex Y44rad(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 4 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians

{
  double r1, r2;
  complex z;
  r1 = -0.1875*sqrt(17.5/PI);   //  4         [  315 ]1/2   4   4*i*phi 
  r2 = sin(theta);		// Y (t,p) = -|------|   sin(t)e
  r2 = r2*r2;			//  4         [512*PI]
  r2 = r2*r2;
  z = complex(0.0, 4.0*phi);
  return r1*r2*exp(z);
} 


complex Y4m4rad(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 4 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians

{
  double r1, r2;
  complex z;
  r1 = 0.1875*sqrt(17.5/PI);    //  -4        [  315 ]1/2   4   -4*i*phi 
  r2 = sin(theta);		// Y  (t,p) = |------|   sin(t)e
  r2 = r2*r2;			//  4         [512*PI]
  r2 = r2*r2;
  z = complex(0.0, -4.0*phi);
  return r1*r2*exp(z);
} 


complex Ylmrad(int l, int m, double theta, double phi)

	// Input		l     : angular momentum
	// 			m     : angular momentum component
	// 			theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank l normalized spherical harmonic
	// Note			      : angles theta and phi input in radians

{
  complex z;
     switch (l)
      {
        case 0:
          z = Y00rad();
          break;
        case 1:
          switch (m)
          {
          case 0:
            z = Y10rad(theta);
            break;
          case 1:
            z = Y11rad(theta, phi);
            break;
          case -1:
            z = Y1m1rad(theta, phi);
            break;
          default:
            Y_error(0, 1);
            Y_error(0, string("1,")+Gdec(m), 1);
            Y_fatality(2);
            break;
          }
          break;
        case 2:
          switch (m)
          {
          case 0:
            z = Y20rad(theta);
            break;
          case 1:
            z = Y21rad(theta, phi);
            break;
          case -1:
            z = Y2m1rad(theta, phi);
            break;
          case 2:
            z = Y22rad(theta, phi);
            break;
          case -2:
            z = Y2m2rad(theta, phi);
            break;
          default:
            Y_error(0, 1);
            Y_error(0, string("2,")+Gdec(m), 1);
            Y_fatality(2);
            break;
          }
          break;
        case 3:
          switch (m)
          {
          case 0:
            z = Y30rad(theta);
            break;
          case 1:
            z = Y31rad(theta, phi);
            break;
          case -1:
            z = Y3m1rad(theta, phi);
            break;
          case 2:
            z = Y32rad(theta, phi);
            break;
          case -2:
            z = Y3m2rad(theta, phi);
            break;
          case 3:
            z = Y33rad(theta, phi);
            break;
          case -3:
            z = Y3m3rad(theta, phi);
            break;
          default:
            Y_error(0, 1);
            Y_error(0, string("3,")+Gdec(m), 1);
            Y_fatality(2);
            break;
          }
          break;
        case 4:
          switch (m)
          {
          case 0:
            z = Y40rad(theta);
            break;
          case 1:
            z = Y41rad(theta, phi);
            break;
          case -1:
            z = Y4m1rad(theta, phi);
            break;
          case 2:
            z = Y42rad(theta, phi);
            break;
          case -2:
            z = Y4m2rad(theta, phi);
            break;
          case 3:
            z = Y43rad(theta, phi);
            break;
          case -3:
            z = Y4m3rad(theta, phi);
            break;
          case 4:
            z = Y44rad(theta, phi);
            break;
          case -4:
            z = Y4m4rad(theta, phi);
            break;
          default:
            Y_error(0, 1);
            Y_error(0, string("4,")+Gdec(m), 1);
            Y_fatality(2);
            break;
          }
          break;
        default:
          Y_error(0, 1);
          Y_error(0, Gdec(l)+string(",")+Gdec(m), 1);
          Y_fatality(2);
          break;
       }
  return z;
  } 

 
// ____________________________________________________________________________
//                 SPHERICAL HARMONICS EMPLOYING DEGREES
// ____________________________________________________________________________


double Y00()

	// Input		none : unnecessary
	// Output		z    : rank zero normalized spherical harmonic

{ return Y00rad(); } 


double Y10(double theta)

	// Input		theta : spherical angle
	// Output		z     : rank one normalized spherical harmonic
	// Note			      : angle theta input in degrees

{ return Y10rad(theta*DEG2RAD); } 


complex Y11(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 1 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees

{ return Y11rad(theta*DEG2RAD, phi*DEG2RAD); } 


complex Y1m1(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 1 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees

{ return Y1m1rad(theta*DEG2RAD, phi*DEG2RAD); } 


double Y20(double theta)

	// Input		theta : spherical angle
	// Output		z     : rank 2 normalized spherical harmonic
	// Note			      : angle theta input in degrees

{ return Y20rad(theta*DEG2RAD); } 


complex Y21(double theta, double  phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 2 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees

{ return Y21rad(theta*DEG2RAD, phi*DEG2RAD); } 


complex Y2m1(double theta, double  phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 2 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees

{ return Y2m1rad(theta*DEG2RAD, phi*DEG2RAD); } 
 

complex Y22(double theta, double  phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 2 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees

{ return Y22rad(theta*DEG2RAD, phi*DEG2RAD); } 


complex Y2m2(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 2 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees

{ return Y2m2rad(theta*DEG2RAD, phi*DEG2RAD); } 


double Y30(double theta)

	// Input		theta : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angle theta input in radians

{ return Y30rad(theta*DEG2RAD); } 


complex Y31(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees
 
{ return Y31rad(theta*DEG2RAD, phi*DEG2RAD); } 


complex Y3m1(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees

{ return Y3m1rad(theta*DEG2RAD, phi*DEG2RAD); } 
 

complex Y32(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees

{ return Y32rad(theta*DEG2RAD, phi*DEG2RAD); } 


complex Y3m2(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees

{ return Y3m2rad(theta*DEG2RAD, phi*DEG2RAD); } 
 

complex Y33(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees

{ return Y33rad(theta*DEG2RAD, phi*DEG2RAD); } 

 
complex Y3m3(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees

{ return Y3m3rad(theta*DEG2RAD, phi*DEG2RAD); } 


/*
double Y40(double theta)

	// Input		theta : spherical angle
	// Output		z     : rank 4 normalized spherical harmonic
	// Note			      : angle theta input in radians

	{ 
	return Y40(theta*DEG2RAD); 
	// *** Should this be return Y40rad(...) ?
	// *** currently, calling this would cause endless recursion.
	} 
*/

complex Y41(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 4 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees
 
{ return Y41rad(theta*DEG2RAD, phi*DEG2RAD); } 


complex Y4m1(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 4 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees

{ return Y4m1rad(theta*DEG2RAD, phi*DEG2RAD); } 
 

complex Y42(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 4 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees

{ return Y42rad(theta*DEG2RAD, phi*DEG2RAD); } 


complex Y4m2(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 4 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees

{ return Y4m2rad(theta*DEG2RAD, phi*DEG2RAD); } 
 

complex Y43(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 4 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees

{ return Y43rad(theta*DEG2RAD, phi*DEG2RAD); } 
 

complex Y4m3(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 4 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees

{ return Y4m3rad(theta*DEG2RAD, phi*DEG2RAD); } 
 

complex Y44(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 4 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees

{ return Y44rad(theta*DEG2RAD, phi*DEG2RAD); } 
 

complex Y4m4(double theta, double phi)

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 4 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees

{ return Y4m4rad(theta*DEG2RAD, phi*DEG2RAD); } 



complex Ylm(int l, int m, double theta, double phi)

	// Input		l     : angular momentum
	// 			m     : angular momentum component
	// 			theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank l normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees

{ return Ylmrad(l, m, theta*DEG2RAD, phi*DEG2RAD); }

#endif						// SphHarmic.cc
