/* Quaternion.cc ************************************************-*-c++-*-
**									**
** 	                         G A M M A				**
**									**
**	Quaternion		                    Implementation      **
**						 			**
**	Copyright (c) 1991, 1999				 	**
**	Scott Smith				 			**
**	Eidgenoessische Technische Hochschule	 			**
**	Labor fur physikalische Chemie		 			**
**	8092 Zurich / Switzerland		 			**
**                                                    		  	**
**      $Header: $
**								 	**
** Thanks to Jacco van Beek (ETHZ) for the suggestions on finding 	**
** proper Euler angles from Quaternions and corrections that make the	**
** definition herein consistent with the use of Euler angles in other 	**
** parts of GAMMA. These were put into GAMMA as of version 4.0.4B and 	**
** affects the functions alpha, gamma, and EA. In version 4.0.6 more	**
** modifications were added (along with testting functions) so that all	**
** seems to work fine.							**
**								 	**
*************************************************************************/

/*************************************************************************
**  									**
** This file contains class quatern, a data type which facilitates the	**
** use of quaternions.	Note that internaly all angles are used in	**
** radians. Angles in degrees occur only in two places: 1.) I/O where	**
** it is convenient to see/input angles in degrees, 2.) Using class	**
** coordinate which, unlike class EAngles which is always in radians,	**
** is assumed to be input and output in degrees. Although coordinates	**
** may be used to hold three Euler angles, THE PREFERRED GAMMA CLASS IS **
** EAngles - it is specially designed to hold three angles and will	**
** restrict those angles to be within alpha=gamma=[0,360], beta=[0,180] **
** and maintains the angles in radians.					**
**  									**
*************************************************************************/

#ifndef   Gquatern_cc_			// Is file already included?
#  define Gquatern_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// this is the implementation
#  endif

#include <Level2/Quaternion.h>		// Include the header file
#include <Basics/Gconstants.h>		// Need constant PI
#include <Basics/Gutils.h>		// Need GAMMA error messages
#include <Basics/StringCut.h>           // Additional string manipulations
#include <Basics/SinglePar.h>           // Include GAMMA single parameters
#include <Basics/ParamSet.h>		// Include GAMMA parameter sets
#include <stdlib.h>
#include <cmath>			// Include HUGE_VAL_VAL

using std::string;			// Using libstdc++ strings
using std::ostream;			// Using libstdc++ output streams
using std::ofstream;			// Using libstdc++ output file streams
using std::cout;			// Using libstdc++ standard output

int    quatern::QRange     = 1;		// Initialize static value
bool   quatern::SCTF       = false;	// Haven't set SinPos,CosPos, or TanPos
bool   quatern::SinPos     = false;	// A default for SinPos
bool   quatern::CosPos     = false;	// A default for CosPos
bool   quatern::TanPos     = false;	// A default for TanPos
double quatern::ElemCutoff = 1.e-12;	// Zero cutoff for element

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                     CLASS QUATERNION ERROR HANDLING
// ____________________________________________________________________________

/*      Input                   Qrt     : A quaternion (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
        Output                  none    : Error message output
                                          Program execution stopped if fatal */
 
void quatern::Qerror(int eidx, int noret) const
  {
  string hdr("Quaternion");
  string msg;
  switch (eidx)
    {
    case 11:GAMMAerror(hdr,"Norm Deviates From 1",         noret); break;// (11)
    case 12:GAMMAerror(hdr,"{A,B,C,D} Is Invalid",         noret); break;// (12)
    case 13:GAMMAerror(hdr,"Invalid 4x4 Rotation Matrix",  noret); break;// (13)
    case 14:GAMMAerror(hdr,"Rotation Mx Not Orthogonal",   noret); break;// (14)
    case 15:GAMMAerror(hdr,"Rotation Mx Must Be 4x4",      noret); break;// (15)
    case 16:GAMMAerror(hdr,"Cannot Rotate Via Matrix",     noret); break;// (16)
    case 17:GAMMAerror(hdr,"A,B,C & D Not In [-1,1]",      noret); break;// (17)
    case 18:GAMMAerror(hdr,"Rot. Mx Row Not Orthogonal",   noret); break;// (18)
    case 19:GAMMAerror(hdr,"Rot. Mx Column Not Orthogonal",noret); break;// (19)
    case 20:GAMMAerror(hdr,"Can't Write To Parameter File",noret); break;// (20)
    case 22:GAMMAerror(hdr,"Can't Write To Output FileStr",noret); break;// (22)
    case 23:GAMMAerror(hdr,"Cannot Output Parameters",     noret); break;// (23)
    case 33:GAMMAerror(hdr,"Cannot Set From Parameters",   noret); break;// (33)
    case 34:GAMMAerror(hdr,"Insufficient Parameters",	   noret); break;// (34)
    case 40:GAMMAerror(hdr,"Problems Finding asin() Range",noret); break;// (40)
    case 41:GAMMAerror(hdr,"Problems Finding acos() Range",noret); break;// (41)
    case 42:GAMMAerror(hdr,"Problems Finding atan() Range",noret); break;// (42)
    case 43:GAMMAerror(hdr,"Cannot Determine Euler Angles",noret); break;// (43)
    default:GAMMAerror(hdr,eidx,noret); break;
    }
  }

void quatern::Qerror(int eidx, const string& pname, int noret) const
  {
  string hdr("Quaternion");
  string msg;
  switch(eidx)
    {
    case 101:							      // (101)
      msg = string("Can't Find Parameters For ")
          + pname;
      GAMMAerror(hdr, msg, noret); break;
    default: GAMMAerror(hdr, eidx, pname, noret); break;
    }
  }

volatile void quatern::Qfatal(int eidx) const
  {
  Qerror(eidx, 1);
  if(eidx) Qerror(0);
  GAMMAfatal();					// Clean exit from program
  }
 
volatile void quatern::Qfatal(int eidx, const string& pname) const
  {
  Qerror(eidx, pname, 1);
  if(eidx) Qerror(0);
  GAMMAfatal();					// Clean exit from program
  }
 
// ____________________________________________________________________________
// ii            PRIVATE FACILITATOR FUNCTIONS FOR QUATERNIONS
// ____________________________________________________________________________

/* These are testing functions for Quaternions.  There are three tests here.
   The 2nd deals with a quaternion build from {A,B,C,D}.  In this case all
   for components must resice within [-1, 1] and the norm must always be one.
   The 3rd test deals with the 4x4 rotation matrix associated with a
   Quaternion. The array must not only be 4x4 but must also be orthonormal.  */

bool quatern::CheckNorm(int warn) const
  {
  if(fabs(norm()-1.0) > 2.*ElemCutoff)
    {
    if(warn)
      {
      Qerror(11, 1);				// Norm deviates from 1
      if(warn <=1) Qerror(12);			// {A,B,C,D} invalid
      else         Qfatal(12);			// {A,B,C,D} invalid, End
      }
    return false;
    }
  return true;
  }

bool quatern::CheckNorm(double A, double B, double C, double D, int warn) const
  {
  double NRM = A*A + B*B + C*C + D*D;
  if(   A<-1 || A>1 || B<-1 || B>1
     || C<-1 || C>1 || D<-1 || D>1 )
    {
    if(warn)
      {
      Qerror(17, 1);				// Value outside of range
      if(warn <=1) Qerror(12);			// {A,B,C,D} invalid
      else         Qfatal(12);			// {A,B,C,D} invalid, End
      }
    return false;
    }
  if(fabs(NRM-1.0) > 2.*ElemCutoff)
    {
    if(warn)
      {
      Qerror(11, 1);				// Norm deviates from 1
      if(warn <=1) Qerror(12);			// {A,B,C,D} invalid
      else         Qfatal(12);			// {A,B,C,D} invalid, End
      }
    return false;
    }
  return true;
  }

bool quatern::CheckNorm(const matrix& Rotmx, bool warn) const
  {
  if((Rotmx.rows() != 4) || (Rotmx.cols() != 4))
    {
    if(warn) Qerror(14,1);			// Rotation mx not 4x4
    return false;
    }

  quatern Qrt;
  for(int i=0; i<4; i++)			// Check Rmx rows
    {
    Qrt.AQ = Rotmx.getRe(i,0); 
    Qrt.BQ = Rotmx.getRe(i,1); 
    Qrt.CQ = Rotmx.getRe(i,2); 
    Qrt.DQ = Rotmx.getRe(i,3); 
    if(fabs(Qrt.norm() - 1.0) > 2.*ElemCutoff)
      {
      if(warn) Qerror(18,1);	
      return false;
      }
    }
  for(int j=0; j<4; j++)			// Check Rmx columns
    {
    Qrt.AQ = Rotmx.getRe(0,j); 
    Qrt.BQ = Rotmx.getRe(1,j); 
    Qrt.CQ = Rotmx.getRe(2,j); 
    Qrt.DQ = Rotmx.getRe(3,j); 
    if(fabs(Qrt.norm() - 1.0) > 2.*ElemCutoff)
      {
      if(warn) Qerror(19,1);
      return false;
      }
    }
  return true;
  }

// ____________________________________________________________________________
// iii         PRIVATE PARAMETER PARSING FUNCTIONS FOR QUATERNIONS
// ____________________________________________________________________________

/* These allow for Quaternions to be set from parameters found in GAMMA
   parameter sets.  The Quaternion is specified by Quatern(#) where # is the
   index of the Quaternion. The parameter is assumed in string format and the
   data in the form (A, B, C) which will be parsed. Once parsed, the values
   will be checked for consistency. The value of D is not specified in the
   parameter because it will be determined from the quaternion norm.         */

bool quatern::SetQuatern(const ParameterSet& pset, int idx, int warn)
  {
  string pname("Quaternion");			// Parameter name
  string app;					// Parameter name suffix
  if(idx != -1) 				// If the index is NOT -1
    app = string("(") + Gdec(idx) + ")"; 	// we need the suffix index
  pname += app;					// Now full parameter name
  ParameterSet::const_iterator item;		// Pix in parameter list
  item  = pset.seek(pname);			// Try & get parameter Pix
  if(item != pset.end())			// If parameter is found we
    {						// need to parse the data 
    string qdata = (*item).data();		//   Here is data (A,B,C,D)
    cutParBlks(qdata);				//   Remove leading ( + blanks
    string sA = cutDouble(qdata);		//   Clip double from front (A)
    cutBlksXBlks(qdata, ",");			//   Clip off blank , blank
    string sB = cutDouble(qdata);		//   Clip another double (B)
    cutBlksXBlks(qdata, ",");			//   Clip off blank , blank
    string sC = cutDouble(qdata);		//   Clip another double (C)
    AQ = atof(sA.c_str());			//   Set 1st to A value
    BQ = atof(sB.c_str());			//   Set 2nd to B value
    CQ = atof(sC.c_str());			//   Set 3rd to C value
    DQ = sqrt(1.0 - AQ*AQ - BQ*BQ - CQ*CQ);	//   Set D value so |Q|=1
    return CheckNorm(warn?1:0);
    }
  if(warn)
    {
    Qerror(33, 1);
    if(warn>1) Qfatal(34);
    else       Qerror(34,1);
    }
  return false;
  }

// ____________________________________________________________________________
// iv                  PRIVATE ANGLE DISCERNMENT FUNCTIONS
// ____________________________________________________________________________

/* These functions help converting a Quaternion into Euler angles. Theyre more
   complicated than they should be because of the multi-valued nature of the
   (inverse) trigonometric functions. Because of this I've put the separate. */

// ----------------------------------------------------------------------------
//             For Discerning Behavior Of asin, acos, & atan Functions
// ----------------------------------------------------------------------------

bool quatern::SetSinPos() const
  {
  double maxa=-10.0;
  double mina= 10.0;
  maxa = gmax(maxa,asin(-1.0));
  mina = gmax(maxa,asin(-1.0));
  maxa = gmax(maxa,asin( 1.0));
  mina = gmax(maxa,asin( 1.0));
       if(mina >= 0.0     && maxa <= PI)     SinPos = true;
  else if(mina >= -PI/2.0 && maxa <= PI/2.0) SinPos = false;
  else { quatern Q; Q.Qerror(40, 1); Q.Qfatal(43); }
  return SinPos;
  }

bool quatern::SetCosPos() const
  {
  double maxa=-10.0;
  double mina= 10.0;
  maxa = gmax(maxa,acos(-1.0));
  mina = gmax(maxa,acos(-1.0));
  maxa = gmax(maxa,acos( 1.0));
  mina = gmax(maxa,acos( 1.0));
       if(mina >= 0.0     && maxa <= PI)     CosPos = true;
  else if(mina >= -PI/2.0 && maxa <= PI/2.0) CosPos = false;
  else { quatern Q; Q.Qerror(41, 1); Q.Qfatal(43); }
  return CosPos;
  }

bool quatern::SetTanPos() const
  {
  double maxa=-10.0;
  double mina= 10.0;
  maxa = gmax(maxa,atan(-HUGE_VAL));
  mina = gmax(maxa,atan(-HUGE_VAL));
  maxa = gmax(maxa,atan( HUGE_VAL));
  mina = gmax(maxa,atan( HUGE_VAL));
       if(mina >= 0.0     && maxa <= PI)     TanPos = true;
  else if(mina >= -PI/2.0 && maxa <= PI/2.0) TanPos = false;
  else { quatern Q; Q.Qerror(42, 1); Q.Qfatal(43); }
  return TanPos;
  }

// ----------------------------------------------------------------------------
//                   For Discerning Which Beta Angle Is Best
// ----------------------------------------------------------------------------

/* The Euler angle beta is given by either of two formulae:

                                     1/2                         1/2
                      -1 {[  2   2 ]    }          -1 {[  2   2 ]    }
          beta = 2 sin   {[ A + B  ]    }  =  2 cos   {[ C + D  ]    }

  Euler angles in GAMMA restrict beta to reside within the range [0, 180].
  However, asin and acos functions in C/C++ will restrict their output values
  (principal values) to within the range [-90, 90] on some system and [0, 180]
  on other systems. So, how best to determine beta? Here we use the asin 
  function. If the range of asin matches the Euler angle range of [0, 180] we
  should have the proper angle. If the range of asin is [-90, 90] then we
  must add PI to any negative betas to map [-90, 0) to angles [90,180).     */

double quatern::FindBeta() const
  {
  double CO = 1.e-7;			// Cutoff for a zero angle
  double arg   = sqrt(AQ*AQ+BQ*BQ);	// First get sqrt(A**2+B**2)
  if(arg > 1.0) arg = 1.0;		// Insure argument is not outside
  if(arg < 0.0) arg = 0.0;		// allowed asin range of [-1,1]
  double betas = 2*asin(arg); 		// Calculate beta using asin
  if(fabs(betas) < CO) 			// If ~ 0 we have 0 or PI
    {					// check using cosine
    arg = sqrt(CQ*CQ + DQ*DQ);		//   Arg is proportional to cos(beta)
    if(arg >0) return 0;		//   If this is + ==> sin is 0
    else       return PI;		//   If this is - ==> sin is PI
    }
  if(betas < 0) betas += PI;		// Else insure beta range [0,PI]
  return betas;
  }

// ----------------------------------------------------------------------------
//                   For Discerning Which Alpha Angle Is Best
// ----------------------------------------------------------------------------

/* The Euler angle alpha is determined from the formula

                                       -1 [ A*D + B*C ]
                            alpha = tan   | --------- |
                                          [ B*D - A*C ]

   As one might imagine from the arc tangent function, there are times when the
   its argument (in the brackets above) will be undetermined. This occurs when
   then denominator is zero or when both the denominator & numerator are zero.
   It has been determined that there are two times when this may occur, when
   the angle beta is zero and when the angle beta is 180. Hence, alternative
   formulae have been derived to handle such situations.

      |        1    -1 [   2*C*D   ]           |          1    -1 [  -2*A*B   ]
 alpha|      = - sin   | --------- |      alpha|        = - sin   | --------- |
      |        2       [ C*C + D*D ]           |          2       [ A*A + B*B ]
       beta=0                                   beta=180

   So, it is evident that a clear determination of angle alpha demands that one
   know the angle beta first, then adjusts the formula used to reflect this.
   Another problem will be the range over which the angle alpha is defined.
   For GAMMA Euler angles alpha spans [0, 360). However, the inverse trig.
   functions (asin, acos, & atan) will only return principal value angles.   */

double quatern::FindAlpha() const
  { return FindAlpha(FindBeta()); }

double quatern::FindAlpha(double beta) const
  {
  double CO = 1.e-7;			// Cutoff for a zero angle
  if(fabs(beta) < CO) 			// Compute alpha when beta = 0
    { 					//	     -1          -1
    double amg = 2.0*GetAngle(CQ, DQ); 	// a-g = 2sin  (C) = 2cos  (D)
    return (amg >= PIx2)?amg-PIx2:amg;
    } 

  if(fabs(fabs(beta)-PI) < CO)		// Compute alpha when beta= PI
    { 					//	     -1          -1
    double amg =2.0*GetAngle(-AQ, BQ); 	// a-g = 2sin  (A) = 2cos  (B)
    return (amg >= PIx2)?amg-PIx2:amg;
    }
    					// 0 < beta & beta != PI 
  double X1 = - AQ / sin(beta/2.);	// sin((alpha-gamma)/2)
  double X2 =   BQ / sin(beta/2.);	// cos((alpha-gamma)/2)
  double X3 =   CQ / cos(beta/2.);	// sin((alpha+gamma)/2)
  double X4 =   DQ / cos(beta/2.);	// cos((alpha+gamma)/2)
  double aminusg = GetAngle(X1, X2);	// 1/2(alpha-gamma) in [0,360)
  double aplusg  = GetAngle(X3, X4);	// 1/2(alpha+gamma) in [0,360)
  double alp     = aplusg + aminusg;	// Here is alpha +/- 2*PI
  double gam     = aplusg - aminusg;	// Here is gamma +/- 2*PI
  if(fabs(alp) < CO) alp = 0;		// Insure no roundoff for 0
  if(fabs(gam) < CO) gam = 0;		// Insure no roundoff for 0

/* We calculate both alpha & gamma above because there is a possible dicrepancy
   due to use of GetAngle in determining alpha-gamma and alpha+gamma. Because
   GetAngle restricts output to [0, 360) we will miss times where alpha-gamma
   is in the range (-360,0) and when alpha+gamma is in the range [360, 720).
   Consequently, the calculated alpha and/or gamma may end up being outside of
   the allowed Euler angle range of [0, 360). To catch these possibilities, we
   make sure BOTH are in the proper range.  If they are not, we adjust the
   sum and difference appropriately until they both do.                      */
   
  if(gam < 0)				// If gamma < 0, adjust alpha to
    {					// compensate 
    if(alp      <     0) alp += PIx2;	//   Add 2*PI if alpha < 0
    else if(alp >= PIx2) alp -= PIx2;	//   Subtract 2*PI if alpha > 0
    }
  else if(alp <     0) alp += PIx2;	// If alpha < 0, put within
  else if(alp >= PIx2) alp -= PIx2;
  if(fabs(alp-PIx2) < 1.e-10) alp = 0; 
  return alp;
  }

// ----------------------------------------------------------------------------
//                   For Discerning Which Gamma Angle Is Best
// ----------------------------------------------------------------------------

/* If we know what the angles beta and alpha are we can use any number of
   functions to determine what the angle gamma is.
*/

double quatern::FindGamma() const
  { double beta = FindBeta(); return FindGamma(beta, FindAlpha(beta)); }

double quatern::FindGamma(double beta, double dalpha) const
  {
  double CO = 1.e-7;			// Cutoff for a zero angle
  if(fabs(beta)    < CO) return 0.0;	// Zero gamma when beta = 0
  if(fabs(beta-PI) < CO) return 0.0;	// Zero gamma when beta = PI
  double X3 =   CQ / cos(beta/2.);	// sin((alpha+gamma)/2)
  double X4 =   DQ / cos(beta/2.);	// cos((alpha+gamma)/2)
  double aplusgo2 = GetAngle(X3, X4);	// 1/2(alpha+gamma) in [0,360)
  double gam = 2.0*aplusgo2 - dalpha;	// This should be gamma
       if(gam > PIx2) gam -= PIx2;
  else if(gam < 0)    gam += PIx2;
  return gam; 
  }

// ----------------------------------------------------------------------------
//         For Discerning Which { Alpha, Beta, Gamma } Angles Are Best
// ----------------------------------------------------------------------------

/* The Euler angles { alpha, beta, gamma } are best determined at the same time
   since the formulae that relate them to { A, B, C, D } involve all three.
*/

EAngles quatern::FindEAs() const
  {
/*                                            2   2 1/2           2   2 1/2
      First Find Angle Beta:  beta/2 = sin[ (A + B )   ] = cos[ (C + D )   ]

   Using asin we will find beta in the range [-90, 90] so we need to convert to
   be in the range [0, 180]. If beta is found >= 0 we are done, except to check
   for exactly 0 & PI due to roundoff. If beta is found < 0 we add PI to get
   the angle in our range.                                                   */ 

  double CO = 1.e-7;			// Cutoff for a zero angle
  double alpha, beta, gamma;		// The Euler angles
  double arg  = sqrt(AQ*AQ+BQ*BQ);	// First get sqrt(A**2+B**2)
  if(arg > 1.0) arg = 1.0;		// Insure argument is not outside
  if(arg < 0.0) arg = 0.0;		// allowed asin range of [-1,1]
  beta = 2.*asin(arg); 			// Calculate beta using asin
  if(fabs(beta) < CO) 			// If ~ 0 we have either 0 or PI
    {					// check using cosine
    arg = sqrt(CQ*CQ + DQ*DQ);		//   Arg is proportional to cos(beta)
    if(arg >0) beta = 0;		//   If this is + ==> sin is 0
    else       beta = PI;		//   If this is - ==> sin is PI
    }
  if(beta < 0) beta += PI;		// Else insure beta range [0,PI]
  if(fabs(beta-PI) < CO) beta=PI;	// If essentially PI, set exactly PI

/* Next Find Blended Alpha & Gamma:

   We calculate both alpha & gamma above because there is a possible dicrepancy
   due to use of GetAngle in determining alpha-gamma and alpha+gamma. Because
   GetAngle restricts output to [0, 360) we will miss times where alpha-gamma
   is in the range (-360,0) and when alpha+gamma is in the range [360, 720).
   Consequently, the calculated alpha and/or gamma may end up being outside of
   the allowed Euler angle range of [0, 360). To catch these possibilities, we
*/

//                  If Beta = 0 One Can Only Determine Alpha+Gamma
//            We Then Set Gamma As Zero And Put The Rotation As All Alpha        
//
//              C = sin[(alpha+gamma)/2]      D = cos[(alpha+gamma)/2]

  if(!beta) 				//	     -1          -1
    {					// a+g = 2sin  (C) = 2cos  (D) => a+0
    gamma = 0;
    alpha = 2.0*GetAngle(CQ, DQ);       // Here, the alpha range is [0, 720)
    if(alpha >= PIx2) alpha -=PIx2;	// Insure alpha range is [0, 360)
    return EAngles(alpha, beta, gamma);	// These are the set of Euler Angles
    } 

//                 If Beta = PI One Can Only Determine Alpha-Gamma
//            We Then Set Gamma As Zero And Put The Rotation As All Alpha        
//
//             A = -sin[(alpha-gamma)/2]      B = cos[(alpha-gamma)/2]

  if(beta == PI) 			//	     -1          -1
    { 					// a-g = 2sin  (A) = 2cos  (B) => a-0
    gamma = 0;
    alpha = 2.0*GetAngle(-AQ, BQ);	// Here, the alpha range is [0, 720)
    if(alpha >= PIx2) alpha -=PIx2;	// Insure alpha range is [0, 360)
    return EAngles(alpha, beta, gamma);	// These are the set of Euler Angles
    }

//     For Beta != PI & Beta !=0 We Can Determine Alpha & Gamma Individually
//
// A = -sin(beta/2) sin[(alpha-gamma)/2]   B = sin(beta/2) cos[(alpha-gamma)/2]
// C =  cos(beta/2) sin[(alpha+gamma)/2]   D = cos(beta/2) cos[(alpha+gamma)/2]

  double X1 = - AQ / sin(beta/2.);	// This is sin((alpha-gamma)/2)
  double X2 =   BQ / sin(beta/2.);	// This is cos((alpha-gamma)/2)
  double X3 =   CQ / cos(beta/2.);	// This is sin((alpha+gamma)/2)
  double X4 =   DQ / cos(beta/2.);	// This is cos((alpha+gamma)/2)
  double aminusgo2 = GetAngle(X1, X2);	// 1/2(alpha-gamma) in [0,360)
  double aplusgo2  = GetAngle(X3, X4);	// 1/2(alpha+gamma) in [0,360)
  alpha = aplusgo2 + aminusgo2;		// Here's alpha (+/- 2PI from GetAngle)
  gamma = aplusgo2 - aminusgo2;		// Here's gamma (+/- 2PI from GetAngle)


  if(fabs(alpha) < CO) alpha = 0;	// Insure no roundoff at alpha=0
  if(fabs(gamma) < CO) gamma = 0;	// Insure no roundoff at gamma=0
   
  if(gamma < 0)				// If gamma < 0, adjust alpha to
    {					// compensate 
    if(alpha      <     0) alpha+=PIx2;	//   Add 2*PI if alpha < 0
    else if(alpha >= PIx2) alpha-=PIx2;	//   Subtract 2*PI if alpha > 0
    }
  else if(alpha <     0) alpha += PIx2;	// If alpha < 0, put within
  else if(alpha >= PIx2) alpha -= PIx2;
  if(fabs(alpha-PIx2) < 1.e-10) alpha = 0; 

       if(gamma > PIx2) gamma -= PIx2;
  else if(gamma < 0)    gamma += PIx2;

  return EAngles(alpha, beta, gamma);
  }

// ----------------------------------------------------------------------------
//             For Discerning Which Angle Is Correct In Span [0, 360]
// ----------------------------------------------------------------------------

/* This function determines the angle that suits the values of sin(angle)
   and cos(angle) within the range of [0, 360). Since sin and cos are multi-  
   valued, the proper angle cannot be determined from either exclusively.
   Therefore, both values are used to figure out the proper angles. Note that
   this function is private because it assumes that the input sine and cosine
   values stem from the same angle, that they reside between [-1, 1], and that
   the output will be within [0, 360) - the range of GAMMA Euler angles.
  
   The function divides the possible angle in to four sections based on whether
   sin and cosine are positive or negative:

            angle      | 0 | <-> | 90 | <-> | 180 | <-> | 270 |
            --------------------------|--------------------------------      
            sin(angle) | 0 |  +  | 1  |  +  |  0  |  -  | -1  |
            cos(angle) | 1 |  +  | 0  |  -  | -1  |  -  |  0  |

   As a final note, this function assumes that the inverse sine function, asin,
   returns within the range -90 to 90. The inversion trig functions provided in
   the C/C++ math library return principle values, i.e. they are always 
   appropriate for the range wherein the function is SINGLE VALUED! These
   ranges are shown below for a typical system.

			  Function	Prin.Val.Range
	  		    cos		    [0, pi]
	  		    sin		 [-pi/2, pi/2]
	  		    tan		 [-pi/2, pi/2]

   This can be problematic because the returned angle is only one of multiple
   (at least two) possibilities. So we must be careful to compensate for that
   when working in the particular angle ranges.

	   Input		sinval	: Sine value of angle
	   			cosval	: Cosine value of angle
	   Output		angle   : Angle which coincides with the
	  				  inputt angles in the range [0, 360).
	   Note			        : This routine determines an angle 
	  				  from its sine and cosine values    */

double quatern::GetAngle(double sinval, double cosval) const
  {
  double CO = 5.e-9;			// Cutoff for a zero angle
  if(fabs(sinval) < CO)			// If sin(angle) = 0, angle is
    return (cosval>0)?0:PI;		// either 0 (cos=1) or PI (cos=-1)

  if(fabs(cosval) < CO)			// If cos(angle) = 0, angle is
    return (sinval>0)?0.5*PI:1.5*PI; 	// either PI/2 or 3PI/2

  double ang = asin(sinval);		// Both sin(angle) & cos(angle) !=0
    					// non-zero. Here ang = [-pi/2, pi/2]
    
  if(sinval > 0)			// If sin(angle) > 0, angle range is
    { 					// within (0,180). Pick via cosine
    if(cosval<0) ang = PI - ang; 	//   If cos(angle) < 0 : 90<angle<180
    } 					//   If cos(angle) > 0 : 0<angle<90
  else					// If sin(angle) <0, angle range is
    {					// within (180,360). Pick via cosine
    if(cosval<0) ang = PI - ang;	//   If cos(angle) < 0 : 180<angle<270
    else	 ang += PIx2;		//   If cos(angle) > 0 : 270<angle<360
    }
  return ang;
  }



// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A               CLASS QUATERNION CONSTRUCTORS/DETRUCTOR
// ____________________________________________________________________________

/* These constructors set up a new quaternion. All demand specification of the
   four quaternion components, either explicitly or by default. These may be
   set using three Euler angles as well.
 
        Input Arguments                     Result
        ---------------       ------------------------------------------
               -              Zero Rotation Quaternion {a=0, b=0, g=0}
            EA,int            Quaternion representing EA or its inverse
            Qrt,int           Quaternion equal to Qrt or its inverse
          QA,QB,QC,QD         Quaternion having these 4 components
			      (Insures A**2 + B**2 + C**2 + D**2 = 1)
 
   Note that the relationship between the Quaternion herein and the set of
   three Euler angles is in accordance with Spiess, JMR, 61, 356 (85). In	
   particular see Equation [3] of that article. Also see the article 
   Ernst, JMR, 63, 133 (85). The norm of the 4 components is always 1.	    */

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

quatern::quatern() { AQ=0; BQ=0; CQ=0; DQ=1; }

quatern::quatern(const quatern& Qrt, bool inv)
  {
  int sign = inv?-1:1; 			// Set sign for inverse/not inverse
  AQ = sign*Qrt.AQ;			// Set {A,B,C,D} accordingly
  BQ = sign*Qrt.BQ;			// The inverse is {-A,-B,-C,D}
  CQ = sign*Qrt.CQ;
  DQ = Qrt.DQ;
  }

// ----------------------------------------------------------------------------
//              Direct Constructors Using A Rotation Designation
// ----------------------------------------------------------------------------

	// Input		Q	: Quaternion (*this)
	//			ABG     : Euler angles (angles in degrees)
	//			EA      : Euler angles (angles in radians)
	//			Qrt 	: A second quaternion
	//			A,B,C,D : Four quaternion components
	//			inv     : Flag whether inverse is desired
	// Output               none    : Dipolar interaction constructed
	// Note                         : A fatal error will result if
	//				  the quaternion norm is not 1
	// Note				: Coordinate I/O is assumed in degrees!
	//				  All other angles herein are radians.

quatern::quatern(const coord& ABC,  bool inv)
  { EAngles EA(ABC); *this = quatern(EA, inv); }
quatern::quatern(const EAngles& EA, bool inv)
  {
  double bo2   = 0.5*EA.beta();			// beta/2
  double apgo2 = 0.5*(EA.alpha()+EA.gamma());	// (alpha+gamma)/2
  double amgo2 = 0.5*(EA.alpha()-EA.gamma());	// (alpha-gamma)/2
  complex za(0,  apgo2);		// We need these two values
  complex zb(0, -amgo2);		// which are temporaries
  complex aq = cos(bo2)*exp(za); 	// cos(b/2)*exp[ i*(a+g)/2]
  complex bq = sin(bo2)*exp(zb); 	// sin(b/2)*exp[-i*(a-g)/2]
  AQ = Im(bq);				// Set 4 quaternion components
  BQ = Re(bq);				// aq = D + iC ; bq = B + iA
  CQ = Im(aq);
  DQ = Re(aq);
  if(inv) { AQ=-AQ; BQ=-BQ; CQ=-CQ; }	// If we want the inverse
  } 					// then negate A, B, and C

quatern::quatern(double QA, double QB, double QC, double QD, bool inv)
  {
  CheckNorm(QA, QB, QC, QD);		// Insure norm is proper
  int sign = inv?-1:1;			// Set sign for inverse/not inverse
  AQ = sign*QA;				// Set the quaternion A value
  BQ = sign*QB;				// Set the quaternion B value
  CQ = sign*QC;				// Set the quaternion C value
  DQ = QD;				// Set the quaternion D value
  }

quatern::quatern(const ParameterSet& pset, int idx, int warn)
  { 
  if(!SetQuatern(pset, idx, warn?1:0))
  if(warn) 
    {
    Qerror(3, 1);
    if(warn > 1) Qfatal(9);
    else         Qerror(9);
    }
  }

// ----------------------------------------------------------------------------
//                    Assignment Operators & Destructor
// ----------------------------------------------------------------------------

// Remember, we assume coord is in degrees whereas EAngles is always in radians
// Also true in class EAngles, any coordinate angles are assumed in degrees.

quatern::~quatern () { }		// No destruction needed here

quatern& quatern::operator= (const coord& ABG)
                                        { *this = quatern(ABG); return *this; }
quatern& quatern::operator= (const EAngles& EA)
                                        { *this = quatern(EA);  return *this; }
quatern& quatern::operator= (const quatern& Q)  
  { 
  if(this == &Q) return *this;		// Do nothing if already equal
  AQ=Q.AQ;BQ=Q.BQ;CQ=Q.CQ;DQ=Q.DQ;
  return *this;
  }

// ____________________________________________________________________________
// B                     QUATERNION ACCESS FUNCTIONS
// ____________________________________________________________________________

/* These allow users to directly access the components of a quaternion.     */

double quatern::A() const { return AQ; }
double quatern::B() const { return BQ; }
double quatern::C() const { return CQ; }
double quatern::D() const { return DQ; }

// ____________________________________________________________________________
//  C                 QUATERNION TO EULER ANGLE FUNCTIONS
//            Q == {A, B, C, D} ==> {alpha, beta, gamma} == EAngles
// ____________________________________________________________________________

/* These functions allow users to obtain the three Euler rotation angles 
   { alpha, beta, gamma } that are the equivalent to the Quaternion.
   The angle(s) is/are returned in units of radians.  These functions are
   implemented from the article

                          Spiess, JMR, 61, 356 (1985)
  
   The functions perform the inverse of eq. [7] in the above article. One might
   also have a look at 
                            Ernst, JMR, 63, 133 (85)

   There are a few specifics worth noting.

   1. There is an ambiguity in determining the angle beta that arises from
      use of the arc sine & arc cosine functions with respect to asin & acos.
      The function FindBeta should take care of this, but to resolve any
      discrepancy beta could be generated via two different formulae and the
      results compared. That is not currently done & should not be needed.
   2. It is far more efficinent to determine the angles alpha & beta at the
      same time since all formula for { A, B, C, D} involve sums & differences
      of these two. Furthermore, it is more efficient to determine these angles
      if the value of beta is known. The functions alpha and beta below do NOT
      take advantage of this, but the function(s) returning all Euler angles
      will.
   3. When the angle 0 or PI it is impossible to distinguish between alpha &
      beta. When beta is zero on can only determine alpha-gamma and when beta 
      is PI one can only determine alpha+gamma. The conversion will set gamma
      to zero and set the alpha angle appropriately.

	Input		Q	: Quaternion
	Output (alpha)	alpha	: Euler angles alpha [0, 360) (rad)
	Output (beta)	beta	: Euler angles beta  [0, 180] (rad)
	Output (gamma)	gamma	: Euler angles gamma [0, 360) (rad)
	Output (EA)   	EAs	: Euler angle set             (rad)
	Output (ABG)   	coord	: Euler angles in coordinate  (deg)          */

double  quatern::alpha() const { return FindAlpha(); }
double  quatern::beta()  const { return FindBeta();  }
double  quatern::gamma() const { return FindGamma(); }
EAngles quatern::EA()    const { return FindEAs(); }
coord   quatern::ABG()   const
  {
  EAngles EAs = FindEAs();
  return coord(EAs.alpha()*RAD2DEG, EAs.beta()*RAD2DEG, EAs.gamma()*RAD2DEG);
  }

// ____________________________________________________________________________
// D                  CLASS QUATERNION ROTATION FUNCTIONS
// ____________________________________________________________________________

/* These functions handle composite rotations, i.e. generation of a quarterion
   that represents two successive rotations.  Either 

  Arguments   Return                           Action 
  ---------   ------  ---------------------------------------------------------
   EA1, EA2     Qrt   Quaternion representing rotation EA1 followed by EA2
   Q1,  Q2      Qrt   Quaternion representing rotation Q1 followed by Q2
   Rmx, Q1      Qrt   Quaternion representing rotation Q1 followed Rmx
     
   Euler angles are input in radians. The composite rotation as implemented
   in found in Spiess, JMR, 61, 356 (85), eq. [9].                            */

quatern quatern::operator*  (const quatern& Q) const { return composite(Q); }
quatern&    quatern::operator*= (const quatern& Q)       { *this= composite(Q); return *this;}
quatern&    quatern::operator&= (const quatern& Q)       { *this= Q.composite(*this); return *this;}
quatern          operator*  (const matrix&  R, const quatern& Q)
                                                   { return composite(R,Q); }

quatern composite(const EAngles& EA1, const EAngles& EA2)
  {
  quatern Q1(EA1);		// Set quaternion equivalent to EA1
  quatern Q2(EA2);		// Set quaternion equivalent to EA2
  return Q1.composite(Q2);
  }

quatern composite(const coord& EA1, const coord& EA2)
  {
  quatern Q1(EA1);		// Set quaternion equivalent to EA1
  quatern Q2(EA2);		// Set quaternion equivalent to EA2
  return Q1.composite(Q2);
  }

quatern composite(const quatern& Qrt1, const quatern& Qrt2)
  { return Qrt1.composite(Qrt2); }		// Deprecated

quatern quatern::composite(const quatern& Q, bool rev) const
  {
  quatern Qcmp;			// Composite: Qcmp = this*Q
  if(rev)			// Order Reversed: Qcmp = Q * this
    {
    Qcmp.AQ  =  Q.DQ*AQ 	// Composite quaternion component A
             -  Q.CQ*BQ		// A = D2A1 - C2B1 + B2C1 + A2D1
             +  Q.BQ*CQ
             +  Q.AQ*DQ;

    Qcmp.BQ  =  Q.CQ*AQ;	// Composite quaternion component B
    Qcmp.BQ +=  Q.DQ*BQ;	// B = C2A1 + D2B1 - A2C1 + B2D1
    Qcmp.BQ -=  Q.AQ*CQ;
    Qcmp.BQ +=  Q.BQ*DQ;

    Qcmp.CQ  = -Q.BQ*AQ;	// Composite quaternion component C
    Qcmp.CQ +=  Q.AQ*BQ;	// C = -B2A1 + A2B1 + D2C1 + C2D1
    Qcmp.CQ +=  Q.DQ*CQ;
    Qcmp.CQ +=  Q.CQ*DQ;

    Qcmp.DQ  = -Q.AQ*AQ;	// Composite quaternion component D
    Qcmp.DQ -=  Q.BQ*BQ;	// A = -A2A1 - B2B1 - C2C1 + D2D1
    Qcmp.DQ -=  Q.CQ*CQ;
    Qcmp.DQ +=  Q.DQ*DQ;
    return Qcmp;
    }

  Qcmp.AQ  =  DQ*Q.AQ 		// Composite quaternion component A
           -  CQ*Q.BQ		// A = D2A1 - C2B1 + B2C1 + A2D1
           +  BQ*Q.CQ
           +  AQ*Q.DQ;

  Qcmp.BQ  =  CQ*Q.AQ;		// Composite quaternion component B
  Qcmp.BQ +=  DQ*Q.BQ;		// B = C2A1 + D2B1 - A2C1 + B2D1
  Qcmp.BQ -=  AQ*Q.CQ;
  Qcmp.BQ +=  BQ*Q.DQ;

  Qcmp.CQ  = -BQ*Q.AQ;		// Composite quaternion component C
  Qcmp.CQ +=  AQ*Q.BQ;		// C = -B2A1 + A2B1 + D2C1 + C2D1
  Qcmp.CQ +=  DQ*Q.CQ;
  Qcmp.CQ +=  CQ*Q.DQ;

  Qcmp.DQ  = -AQ*Q.AQ;		// Composite quaternion component D
  Qcmp.DQ -=  BQ*Q.BQ;		// A = -A2A1 - B2B1 - C2C1 + D2D1
  Qcmp.DQ -=  CQ*Q.CQ;
  Qcmp.DQ +=  DQ*Q.DQ;
  return Qcmp;
  }

quatern composite(const matrix& Rotmx, const quatern& Qrt1)
  {
  quatern Qrt;
  if(!(Qrt.CheckNorm(Rotmx)))
    { Qrt1.Qfatal(14); }
  Qrt.AQ  = Rotmx.getRe(0,0)*Qrt1.AQ;
  Qrt.AQ += Rotmx.getRe(0,1)*Qrt1.BQ;
  Qrt.AQ += Rotmx.getRe(0,2)*Qrt1.CQ;      // [A] = [ R11 R12 R13 R14 ][A1]
  Qrt.AQ += Rotmx.getRe(0,3)*Qrt1.DQ;      // [B] = [ R21 R22 R23 R24 ][A2]
  Qrt.BQ  = Rotmx.getRe(1,0)*Qrt1.AQ;      // [C] = [ R31 R32 R33 R34 ][A3]
  Qrt.BQ += Rotmx.getRe(1,1)*Qrt1.BQ;      // [D] = [ R41 R42 R43 R44 ][A4]
  Qrt.BQ += Rotmx.getRe(1,2)*Qrt1.CQ;
  Qrt.BQ += Rotmx.getRe(1,3)*Qrt1.DQ;
  Qrt.CQ  = Rotmx.getRe(2,0)*Qrt1.AQ;
  Qrt.CQ += Rotmx.getRe(2,1)*Qrt1.BQ;
  Qrt.CQ += Rotmx.getRe(2,2)*Qrt1.CQ;
  Qrt.CQ += Rotmx.getRe(2,3)*Qrt1.DQ;
  Qrt.DQ  = Rotmx.getRe(3,0)*Qrt1.AQ;
  Qrt.DQ += Rotmx.getRe(3,1)*Qrt1.BQ;
  Qrt.DQ += Rotmx.getRe(3,2)*Qrt1.CQ;
  Qrt.DQ += Rotmx.getRe(3,3)*Qrt1.DQ;
  return Qrt;
  }

	// Input		Qrt	: Quaternion
	// Output		Rmx	: 4x4 quaternion rotation matrix
	// Note			        : See Spiess, JMR, 61, 356 (85)	
	//				  their eq. [9] is implemented
     
matrix quatern::RotMx() const { return RMx(); }	// DEPRECATED
matrix quatern::RMx() const
  {
  matrix Rmx(4,4);
  Rmx.put( DQ,0,0); Rmx.put(-CQ,0,1); Rmx.put( BQ,0,2); Rmx.put(AQ,0,3);
  Rmx.put( CQ,1,0); Rmx.put( DQ,1,1); Rmx.put(-AQ,1,2); Rmx.put(BQ,1,3);
  Rmx.put(-BQ,2,0); Rmx.put( AQ,2,1); Rmx.put( DQ,2,2); Rmx.put(CQ,2,3);
  Rmx.put(-AQ,3,0); Rmx.put(-BQ,3,1); Rmx.put(-CQ,3,2); Rmx.put(DQ,3,3);
  return Rmx;
  }
 
// ____________________________________________________________________________
// E                  CLASS QUATERNION AUXILIARY FUNCTIONS
// ____________________________________________________________________________

	// Input		Qrt	: Quaternion (this)
	// Output		Qnorm	: Computes a "norm" for the input
	//				  Quaternion by the formula
	//				      A**2+B**2+C**2+D**2
	// Note				: All valid quaternions should
	//				  have a norm of 1
 
double  quatern::norm()    const { return AQ*AQ + BQ*BQ + CQ*CQ + DQ*DQ; }
quatern quatern::inverse() const { return quatern(*this, true); }

// ____________________________________________________________________________
// F                     CLASS QUATERNION I/O FUNCTIONS
// ____________________________________________________________________________

// ------------------------ ASCII Output Functions ----------------------------

/*              Input           Quar : Quaternion (this)
                                ostr : Output ASCII file stream
                                full : Flag for amount of output
                Return          void : Quart is sent to the output stream   */

ostream& quatern::print(ostream& ostr, bool nf, bool hdr) const
  {
  if(hdr) ostr << "Quaternion: ";
  ostr << "A= " << Gform("%7.4f", AQ) << "  " 
       << "B= " << Gform("%7.4f", BQ) << "  "
       << "C= " << Gform("%7.4f", CQ) << "  "
       << "D= " << Gform("%7.4f", DQ);
  if(nf) ostr << "  " << "|Q|= " << Gform("%7.4f", norm());
  return ostr;
  }

ostream& operator << (ostream& ostr, const quatern& Quar)
  { return Quar.print(ostr); }

// ____________________________________________________________________________
// G              Class Quaternion Container Support Functions
// ____________________________________________________________________________

/* These allow the user to make STL lists and vectors from quaternions. If
   not present some compilers complain about code such as list<quatern>      */    

bool quatern::operator== (const quatern& Quar) const
  {
  if(fabs(AQ-Quar.AQ) > ElemCutoff) return false;
  if(fabs(BQ-Quar.BQ) > ElemCutoff) return false;
  if(fabs(CQ-Quar.CQ) > ElemCutoff) return false;
  if(fabs(DQ-Quar.DQ) > ElemCutoff) return false;
  return true;
  }

bool quatern::operator!= (const quatern& Quar) const
  {
  if(fabs(AQ-Quar.AQ) > ElemCutoff) return true;
  if(fabs(BQ-Quar.BQ) > ElemCutoff) return true;
  if(fabs(CQ-Quar.CQ) > ElemCutoff) return true;
  if(fabs(DQ-Quar.DQ) > ElemCutoff) return true;
  return false;
  }

bool quatern::operator<  (const quatern& Quar) const
  {
  if(AQ<Quar.AQ) return true;
  if(AQ>Quar.AQ) return false;
  if(BQ<Quar.BQ) return true;
  if(BQ>Quar.BQ) return false;
  if(CQ<Quar.CQ) return true;
  if(CQ>Quar.CQ) return false;
  if(DQ<Quar.DQ) return true;
  return false;
  }

bool quatern::operator>  (const quatern& Quar) const
  {
  if(AQ>Quar.AQ) return true;
  if(AQ<Quar.AQ) return false;
  if(BQ>Quar.BQ) return true;
  if(BQ<Quar.BQ) return false;
  if(CQ>Quar.CQ) return true;
  if(CQ<Quar.CQ) return false;
  if(DQ>Quar.DQ) return true;
  return false;
  }

// ____________________________________________________________________________
// H                  Class Quaternion Range Functions
// ____________________________________________________________________________
 
/* The Euler angles associated with a Quaternion can in principle have any 
   value. However, most users prefer to limit the range over which the three
   angles are defined.  These functions provide a means by which the Euler
   angles input to and returned from this class can be limited to reside
   within a specified range.  The range depends upon the flag "_range". Users
   can associate values of range with implied Euler angle limits by adjusting
   the functions below.  Zero is reserved for no imposed limits.  As long as
   all 3D space can be sampled by the Euler angle range, this results in no
   loss in gerenality (multiple Euler angles sets produce the same rotation).
   They only affect of range will be in the Euler angle values output.  The
   quarternions and the associated rotation(s) will remain the same.        */

//int  quatern::range() const { return quatern::_Qrange; } 
//void quatern::range(int r)  { quatern::_Qrange = r; }

// ____________________________________________________________________________
// I            Class Quaternion Parameter & Parameter Set Functions
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//                      Single Parameter From Quaternion
//-----------------------------------------------------------------------------

	// Input		Qrt	: Quaternion (this)
        //			pname	: A parameter name
        //			pstate	: A parameter statement
        // Return		par	: A GAMMA parameter of type Quaternion
        //				  with default name and comment unless
	//				  specified

SinglePar quatern::param() const
  {
  string pname  = "Quaternion";		// Default statement
  string pstate = "A Quaternion";	// Default statement
  return param(pname, pstate);          // Use overload function
  }

SinglePar quatern::param(const string& pname) const
  {
  string pstate = "A Quaternion";	// Default statement
  return param(pname, pstate);          // Use overload function
  }

SinglePar quatern::param(const string& pname, const string& pstate) const
  {
  string ot("%g");
  string pdata = "(";			// Forming coordinate data string for
  pdata += Gform(ot, AQ);		// parameter output in ASCII. This is
  pdata += string(", ");
  pdata += Gform(ot, BQ);		//      ( A.aaa, B.bbb, C.ccc)
  pdata += string(", "); 
  pdata += Gform(ot, CQ);		// The value of D is not output to
  pdata += string(") ");		// avoid roundoff errors since |Q|=1
  SinglePar par(pname,3,pdata,pstate);	// Now construct parameter (coord)
  return par;				// Return the parameter
  }

//-----------------------------------------------------------------------------
//                       Parameter Set From Quaternion
//-----------------------------------------------------------------------------

quatern::operator ParameterSet( ) const
  { ParameterSet pset; pset += *this; return pset; }

void operator+= (ParameterSet& pset, const quatern& Qrt)
  { Qrt.PSetAdd(pset); }

bool quatern::PSetAdd(ParameterSet& pset, int idx, int pfx) const
  {
  string suffx;                                 // Parameter suffix
  if(idx != -1)                                 // Only use suffix if idx
    suffx = string("(")+Gdec(idx)+string(")");  // is NOT -1
  string prefx;                                 // Parameter prefix
  if(pfx != -1)                                 // Only use prefix if pfx
    prefx = string("[")+Gdec(pfx)+string("]");  // is NOT -1
  string pname = prefx+string("Quatern")+suffx;	// Quaternion parameter name
  SinglePar par = param(pname);			// Quaternion parameter
  pset.push_back(par);				// Add parameter to pset
  return true;
  }

//-----------------------------------------------------------------------------
//                       Parameter Set File From Quaternion
//-----------------------------------------------------------------------------

	// Input		Qrt	: Quaternion (this)
        //                      filename: Output file name
        //                      ofstr   : Output file stream
        //                      sfx     : Quaternion index  (default -1)
        //                      pfx     : Quaternion prefix (default -1)
        //                      warn    : Warning level
        // Output               none    : Quaterion is written as a parameter
	//				  to file or output file stream

bool quatern::write(const string& filename, int idx, int pfx, int warn) const
  {
  ofstream ofstr(filename.c_str());     // Open filename for input
  if(!write(ofstr, idx, pfx, warn?1:0))	// If file bad then exit
    {
    Qerror(40, filename, 1);		// Filename problems
    if(warn>1) Qfatal(20);		// Fatal error
    return false;
    }
  ofstr.close();                        // Close it now
  return true;
  }

bool quatern::write(ofstream& ofstr, int idx, int pfx, int warn) const
  {
  ParameterSet pset;                    // Declare a parameter set
  PSetAdd(pset, idx, pfx);              // Add in interaction parameters
  if(!pset.write(ofstr, warn?1:0)) 	// Use parameter set to write
    {
    if(warn)
      {
      Qerror(22, 1);			// Problems writing to filestream
      if (warn>1) Qfatal(23); 		// Fatal error
      }
    return false;
    }
  return true;
  }

//-----------------------------------------------------------------------------
//                       Quaternion From Parameter Set
//-----------------------------------------------------------------------------

	// Input		Qrt	: Quaternion (this)
        //                      filein  : An input (ASCII) file
        //                      pset    : A parameter set
        //                      indx    : Quaterion index
        //                      warn    : Warning level
        //                                  0 = no warnings
        //                                  1 = non-fatal warnings
        //                                  2 = fatal warnings
        // Output               none    : Quaterion filled with
        //                                values specified in filein
	//				  or in pset

bool quatern::read(const string& filein, int indx, int warn)
  {
  ParameterSet pset;			// Declare a parameter set
  if(!pset.read(filein, warn?1:0))      // Read in pset from file
    {
    if(warn)
      {
      Qerror(40, filein, 1);		// Problems with file filein
      if(warn>1) Qfatal(101,filein);	// Can't read from file filein
      else       Qerror(101,1);
      }
    return 0;
    }
  return read(pset, indx, warn);
  }

bool quatern::read(const ParameterSet& pset, int indx, int warn)
  {
  if(SetQuatern(pset,indx,warn?1:0))
    return true;
  if(warn)                                      // Looks like we can't read
    {                                           // the quaternion at all
    string sl("");
    if(indx != -1) sl = string(Gdec(indx));
    Qerror(77, sl, 1);
    if(warn > 1)  Qfatal(5);
    else          Qerror(78, sl);
    }
  return false;
  }

// ____________________________________________________________________________
// J             INVERSE TRIGONOMETRIC ANGLE DISCERNMENT FUNCTIONS
// ____________________________________________________________________________

/* These function help converting a Quaternion into Euler angles. Such
   conversions are more complicated than one would like because of the 
   multi-valued nature of the (inverse) trigonometric functions. Because of
   this I've tried to simplyfy the procedure by 1st determininng how the
   computer sees the inverse functions asin & acos.

   The inverse sine and cosine functions in C/C++, asin and acos, will restrict
   their output values to be "principal values", always residing within a range
   over which sine and cosine are single valued. On some systems this is set
   to be the range of angles [-90, 90] whereas on some systems it is set to be
   the range [0, 180]. For example, on my Sun SPARC system  asin() returns in 
   the range [-pi/2, pi/2] radians and acos() returns in the range [0,pi] 
   radians.  The same is true for my SGI Irix systems. However, on my PC/Cygwin
   systems the ranges on both are [-pi/2, pi/2]. These functions just test a 
   few values to discern what the standard is on the computer running this.

   Since acos and asin must span [-1,1], the range of asin will likely be
   [-pi/2, pi/2] and the range of acos can be either.                        */

bool quatern::ASinPos()
  {
  if(SCTF) return SinPos;
  quatern Q;
  Q.SetSinPos();
  Q.SetCosPos();
  Q.SetTanPos();
  return SinPos;
  }

bool quatern::ACosPos()
  {
  if(SCTF) return CosPos;
  quatern Q;
  Q.SetSinPos();
  Q.SetCosPos();
  Q.SetTanPos();
  return CosPos;
  }

bool quatern::ATanPos()
  {
  if(SCTF) return TanPos;
  quatern Q;
  Q.SetSinPos();
  Q.SetCosPos();
  Q.SetTanPos();
  return TanPos;
  }

// ----------------------------------------------------------------------------
//                 For Outputting Tests During Conversions
// ----------------------------------------------------------------------------

void quatern::ShowConversion() const
  {
  cout << "\nConversion Of Quaternion To Euler Angles\n";
  double alphas=0, betas=0, gammas=0;
  double arg;
  double CO    = 1.e-7;			// Cutoff for a zero angle

//                     Show Generation of Angle Beta
//             This Should Mimick The Get

  double arg1  = sqrt(AQ*AQ+BQ*BQ);	// First get sqrt(A**2+B**2)
  if(arg1>1.0) arg1=1.0;
  if(arg1<0.0) arg1=0.0;
  betas = 2*asin(arg1); 		// Calculate beta using asin
  cout << "\n\n";
  cout << "\n\t" << *this;
  cout << "\n\tDeterminining Angle Beta";
  cout << "\n\t\tsqrt(AQ^2 + BQ^2) = " << arg1;
  cout << "\n\t\tasin(sqrt(AQ^2 + BQ^2)) = " << asin(arg1);
  cout << "\n\t\t2*asin(sqrt(AQ^2 + BQ^2)) = beta = " << betas*RAD2DEG;
  if(fabs(betas) < CO)
    {
    cout << "\n\t\tNear To Zero -> Exactly Zero Or PI";
    arg = sqrt(CQ*CQ + DQ*DQ);          //   Arg is proportional to cos(beta)
    if(arg >0) betas = 0;                //   If this is + ==> sin is 0
    else       betas = PI;               //   If this is - ==> sin is PI
    }
  else if(betas < 0)
    {
    cout << "\n\t\tNegative Beta, Adding PI";
    betas += PI;            // Else insure beta range [0,PI]
    }
  cout << "\n\t\tFinal Beta Value: " << betas*RAD2DEG;

// FindAl
  cout << "\n\n\tDeterminining Angle Alpha";
  if(fabs(betas) < CO) 			// Compute alpha when beta = 0
    { 					//	                 C
    alphas = 2*atan(CQ/DQ);		//    alpha = 2 * arctan -
          				//                       D
    cout << "\n\t\tFound That Beta Is Zero";
    cout << "\n\t\tAlpha = 2atan(C/D): " << alphas*RAD2DEG;
    }

  else if(fabs(fabs(betas)-PI) < CO)		// Compute alpha when beta= PI
    {
    double X1 = - AQ / sin(betas/2.);	// sin((alpha-gamma)/2)
    double X2 =   BQ / sin(betas/2.);	// cos((alpha-gamma)/2)
    cout << "\n\t\tFound That Beta Is PI";
    cout << "\n\t\tAlpha = -AQ/sin(beta/2): " << X1*RAD2DEG;
    cout << "\n\t\tAlpha =  BQ/sin(beta/2): " << X2*RAD2DEG;
    alphas = 2.*GetAngle(X1, X2);		// Determine angle
    cout << "\n\t\tChoosen Alpha Is: " << alphas;
    }
    					// 0 < beta & beta != PI 
  else
    {
    double X1 = - AQ / sin(betas/2.);		// X1 = sin((alpha-gamma)/2)
    double X2 =   BQ / sin(betas/2.);		// X2 = cos((alpha-gamma)/2)
    double X3 =   CQ / cos(betas/2.);		// X3 = sin((alpha+gamma)/2)
    double X4 =   DQ / cos(betas/2.);		// X4 = cos((alpha+gamma)/2)
    double aminusgo2 = GetAngle(X1, X2);	// 1/2(alpha-gamma) in [0,360)
    double aplusgo2  = GetAngle(X3, X4);	// 1/2(alpha+gamma) in [0,360)
    cout << "\n";
    cout << "\n\t\t-AQ/sin(beta/2): " << X1;
    cout << "\n\t\t BQ/sin(beta/2): " << X2;
    cout << "\n\t\t CQ/cos(beta/2): " << X3;
    cout << "\n\t\t DQ/cos(beta/2): " << X4;
    cout << "\n";
    cout << "\n\t\tPossible 1/2(Alpha - Gamma) Is: " << asin(X1)*RAD2DEG  << " (asin)";
    cout << "\n\t\tPossible 1/2(Alpha - Gamma) Is: " << acos(X2)*RAD2DEG  << " (acos)";
    cout << "\n\t\tChoosen  1/2(Alpha - Gamma) Is: " << aminusgo2*RAD2DEG << " (comp)";
    cout << "\n\t\t-AQ/sin(beta/2) Is Then:        " << sin(aminusgo2);
    cout << "\n\t\t BQ/sin(beta/2) Is Then:        " << cos(aminusgo2);
    cout << "\n\t\tPossible 1/2(Alpha + Gamma) Is: " << asin(X3)*RAD2DEG << " (asin)";
    cout << "\n\t\tPossible 1/2(Alpha + Gamma) Is: " << acos(X4)*RAD2DEG << " (acos)";
    cout << "\n\t\tChoosen  1/2(Alpha + Gamma) Is: " << aplusgo2*RAD2DEG << " (comp)";
    cout << "\n\t\t CQ/cos(beta/2) Is Then:        " << sin(aplusgo2);
    cout << "\n\t\t DQ/cos(beta/2) Is Then:        " << cos(aplusgo2);

    cout << "\n";
    cout << "\n\t\tPossible Alpha - Gamma Is: " << 2.0*asin(X1)*RAD2DEG << " (asin)";
    cout << "\n\t\tPossible Alpha - Gamma Is: " << 2.0*acos(X2)*RAD2DEG << " (acos)";
    cout << "\n\t\tChoosen  Alpha - Gamma Is: " << 2.0*aminusgo2*RAD2DEG  << " (comp)";
    cout << "\n\t\tPossible Alpha + Gamma Is: " << 2.0*asin(X3)*RAD2DEG << " (asin)";
    cout << "\n\t\tPossible Alpha + Gamma Is: " << 2.0*acos(X4)*RAD2DEG << " (acos)";
    cout << "\n\t\tChoosen  Alpha + Gamma Is: " << 2.0*aplusgo2*RAD2DEG   << " (comp)";

    cout << "\n";
    cout << "\n\t\tPossible Alpha - Gamma Value 1:" << 2.0*aminusgo2*RAD2DEG;
    cout << "\n\t\tPossible Alpha - Gamma Value 2:" << (2.0*aminusgo2+PIx2)*RAD2DEG;
    cout << "\n\t\tPossible Alpha - Gamma Value 3:" << (2.0*aminusgo2-PIx2)*RAD2DEG;
    cout << "\n\t\tPossible Alpha + Gamma Value 1:" << 2.0*aplusgo2*RAD2DEG;
    cout << "\n\t\tPossible Alpha + Gamma Value 2:" << (2.0*aplusgo2+PIx2)*RAD2DEG;
    cout << "\n\t\tPossible Alpha + Gamma Value 3:" << (2.0*aplusgo2-PIx2)*RAD2DEG;
//    double alp     = aplusgo2 + aminusgo2;	// Here is alpha +/- 2*PI
//    double gam     = aplusgo2 - aminusgo2;	// Here is gamma +/- 2*PI

   
    }
  cout << "\n\tFinal Result: " << EAngles(alphas, betas, gammas); 
  cout << "\nEnd Conversion Of Quaternion To Euler Angles\n\n";
  }

// ----------------------------------------------------------------------------
//               Check For Valid Quaternion Rotation Matrix
// ----------------------------------------------------------------------------

bool quatern::ValidRMx(const matrix& R, bool msgs)
  {

//                        Insure Array Is 4 X 4

  if(R.rows() != 4 || R.cols() !=4)
    {
    if(msgs)
      {
      cout << "\n\tInvalid Quaternion Rotation Matrix!";
      cout << "\n\tQuaternion Rotation Matrices Must Be (4 x 4)";
      cout << "\n\tTested Array Is Of Dimension (" << R.rows()
                                          << " x " << R.cols() << ")";
      }
    return false;
    }

//                        Insure Array Is Real

  if(!R.is_real())
    {
    if(msgs)
      {
      cout << "\n\tInvalid Quaternion Rotation Matrix!";
      cout << "\n\tQuaternion Rotation Matrices Must Be (4 x 4)";
      cout << "\n\tTested Array Is Not Real";
      }
    return false;
    }

//                            Insure Rows Are Normal

  quatern Qrt;					// Container for 4 elems
  int i;
  for(i=0; i<4; i++)                            // Check Rmx rows i
    {						// Store row as a quaternion
    Qrt.AQ = R.getRe(i,0);	
    Qrt.BQ = R.getRe(i,1);
    Qrt.CQ = R.getRe(i,2);
    Qrt.DQ = R.getRe(i,3);
    if(fabs(Qrt.norm() - 1) > 2.*ElemCutoff)
      {
      if(msgs)
        {
        cout << "\n\tInvalid Quaternion Rotation Matrix!";
        cout << "\n\tQuaternion Rotation Matrix Rows Are Normal";
        cout << "\n\tTested Array Row " << i << " Norm Is " << Qrt.norm();
        }
      return false;
      }
    }

//                            Insure Columns Are Normal

  int j;
  for(j=0; j<4; j++)                            // Check Rmx rows i
    {						// Store row as a quaternion
    Qrt.AQ = R.getRe(0,j);	
    Qrt.BQ = R.getRe(1,j);
    Qrt.CQ = R.getRe(2,j);
    Qrt.DQ = R.getRe(3,j);
    if(fabs(Qrt.norm() - 1) > 2.*ElemCutoff)
      {
      if(msgs)
        {
        cout << "\n\tInvalid Quaternion Rotation Matrix!";
        cout << "\n\tQuaternion Rotation Matrix Columns Are Normal";
        cout << "\n\tTested Array Column " << j << " Norm Is " << Qrt.norm();
        }
      return false;
      }
    }

  return true;
  }

#endif						// Quaternion.cc
