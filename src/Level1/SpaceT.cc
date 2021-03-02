/* SpaceT.cc ****************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**	Spatial Tensors	                         Implementation		**
**                                                                      **
**	Copyright (c) 1990, 1991, 1992					**
**	Scott Smith							**
**	Eidgenoessische Technische Hochschule				**
**	Labor fur physikalische Chemie					**
**	8092 Zurich / Switzerland					**
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description:	                                                        **
**                                                                      **
** This file contains the implementation of general spatial tensors.    **
** Note that for Magnetic Resonance simulations, interactions are by	**
** and large of rank 2.  Thus, it is rank 2 tensors which are almost	**
** always used.  With this in mind, GAMMA also supplies typical rank	**
** 2 interactions explicitly and they are usually more conveient to use **
** in simulations.  These should be used when higher ranked tensors are	**
** desired (rare) and/or one wishes to work with reducible rank 2       **
** tensors, i.e. ranks 0,1,&2 simultaneously. The latter is done for	**
** example in treating anti-symmetric shift anisotropy (rank 1) along	**
** with the isotropic (rank 0, typical shift) and symmetric (rank 2) 	**
** parts of the SA tensor.                                              **
**                                                                      **
*************************************************************************/

#ifndef   Gspace_T_cc_			// Is this file already included?
#  define Gspace_T_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Level1/SpaceT.h>		// Include the inferface
#include <Level1/SpinT.h>		// Know about spin tensors
#include <Level1/Wigner.h>		// Know these spatial functions
#include <Level1/coord.h>		// Know coordinates (convenience)
#include <Basics/StringCut.h>		// Include string parsing
#include <Basics/Gutils.h>		// Include GAMMA std errors
#include <cstdlib>			// Include functions atoi, atof

using std::string;			// Using libstdc++ strings
using std::cout;			// Using libstdc++ standard output
using std::list;			// Using libstdc++ STL lists
using std::ostream;			// Using libstdc++ output streams
using std::ofstream;			// Using libstdc++ output file streams
using std::ifstream;			// Using libstdc++ input file streams

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// i                     SPACE TENSOR ERROR HANDLING
// ____________________________________________________________________________

                                                                               
void space_T::SphTerror(int eidx, int noret) const
 
        // Input                SphT	: Spatial tensor (this)
        //                      eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message

 
  {
  string hdr("Spatial Tensor");
  string msg;
  switch (eidx)
    {
    case 6:  GAMMAerror(hdr,"Bad Internal Component Access",noret);break;//(6)
    case 7:  GAMMAerror(hdr,"Accessing Empty Operator",noret);    break;// (7)
    case 20: msg = string("Current Tensor Rank is ") + Gdec(rank);
             GAMMAerror(hdr, msg, noret);  break;			// (20)
    case 21: msg = string("Component Rank Must Be [0,")
                 + Gdec(rank) + string("]");
             GAMMAerror(hdr, msg, noret);  break;			// (21)
    case 22: msg = string("Component m Must Be [")
                 + Gdec(-rank) + string(",");
//                 + Gdec(rank) + string("]");
             GAMMAerror(hdr, msg, noret);  break;			// (22)
    case 23: GAMMAerror(hdr,"Returning Zero Value", noret);	  break;// (23)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  if(!noret) cout << "\n";
  }  

void space_T::SphTerror(int eidx, const string& pname, int noret) const
 
        // Input                SphT	: Spatial tensor (this)
        //                      eidx    : Error index
        //                      pname   : Additional error message
        //                      noret   : Flag for linefeed (0=linefeed)
        // Output               none    : Error message output

/* The following error messages use the defaults set in the Gutils package

                Case                          Error Message                

                (1)                     Problems With File PNAME
                (2)                     Cannot Read Parameter PNAME
                (3)                     Invalid Use Of Function PNAME
                default                 Unknown Error - PNAME                */  

  {                                                                             
  string hdr("Spatial Tensor");
  string msg;
  switch(eidx)
    {
    case 5: msg = string("Bad Use Of ") + pname + string(" Function ");
             GAMMAerror(hdr, msg, noret);  break;			// (5)
    case 10:msg = string("Use of Deprecated Function ") + pname;
             GAMMAerror(hdr, msg, noret);  break;			// (10)
    case 11:msg = string("Please Use Function ") + pname;
             GAMMAerror(hdr, msg, noret);  break;			// (11)
    case 20:msg = string("Requested Component With Rank ") + pname;
             GAMMAerror(hdr, msg, noret);  break;			// (20)
    case 21:msg = string("Requested Component With m = ") + pname;
             GAMMAerror(hdr, msg, noret);  break;			// (21)
    default: GAMMAerror(hdr, eidx, pname, noret); break;// Unknown Error  (-1)
    }
  }


volatile void space_T::SphTfatality(int eidx) const
       
        // Input                SphT	: Spatial tensor (this)
        //                      eidx    : Flag for error type
        // Output               none    : Error message output
        //                                Program execution stopped

  {                                                                 
  SphTerror(eidx, 1);				// Output error message
  if(eidx) SphTerror(0);			// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }


volatile void space_T::SphTfatality(int eidx, const string& pname) const
       
        // Input                SphT	: Spatial tensor (this)
        //                      eidx    : Flag for error type
        //                      pname   : Additional error message
        // Output               none    : Error message output
        //                                Program execution stopped

  {                                                                 
  SphTerror(eidx, pname, 1);			// Output error message
  if(eidx) SphTerror(0);			// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

void space_T_error(int i)

	// Input		i    : Error Flag
	// Output		none : Error Message Output

  {
  cout << "\nSpace_T: ";
  switch (i)
    {
    case 0:							// (0)
      cout << "\nSpace_T: Unknown A"
           << "\nSpace_T:          ";
      break;
    case 4:							// (4)
      cout << "Construction From Improper GAMMA Parameter";
      break;
    case 5:							// (5)
      cout << "Non Rank 2 Construction From GAMMA Parameter";
      break;
    case 8:							// (8)
      cout << "Negative Rank Designation";
      break;
    case 9:							// (9)
      cout << "Multiplication of Irreducible Spherical Tensors";
      break;
    case 10:
      cout << "Unable to Determine Spherical Tensor Component.\n";
      break;
    case 11:
      cout << "Wrong Size Vector for Construction of Rank 1 A Tensor.\n";
      break;
    case 12:
      cout << "Wrong Size Matrix for Construction of Rank 2 A Tensor.\n";
      break;
    case 20:							// (20)
      cout << "Cannot Set Isotropic Component, Tensor Rank Is 0.\n";
      break;
    case 21:							// (21)
      cout << "Principal Axes Are Defined Only For Tensors Of Rank 2.\n";
      break;
    default:
      cout << "Unknown error.\n";
      break;
    }
  return;
  }


void volatile space_T_fatality (int error)

	// Input		none :
	// Output		none : Stops Execution & Error Message Output

  {
  space_T_error(error);
  cout << "\nSpace_T: Program Aborting.\n";
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii                   SPACE TENSOR PRIVATE AUXILIARY FUNCTIONS
// ____________________________________________________________________________

 
void space_T::updatePAS()
 
        // Input                SphT : Spatial tensor (this)
        // Output               none : The interal PAS values are set 
	// Note			     : Assumes that changes to A   must be
	//			       relfected in PAS         l,m 

  {
  double Aiso=0;
  double delzz=0;
  double eta=0;
  if(vx[0])				// See if A   value exists
    { 					//         00
    complex A00 = (*(vx[0]))(0);
    Aiso = -zRe(A00)/sqrt(3.0);
    }
  if(vx[2])				// See if A   values exists
    { 					//         2m
    complex A20 = (*(vx[2]))(2);
    complex A22 = (*(vx[2]))(4);
    delzz = zRe(A20)/sqrt(1.5);
    eta   = zRe(A22)/(0.5*delzz);
    }
  PAS = coord(Aiso, delzz, eta);
  return;	
  }


// ____________________________________________________________________________
// A                  SPACE TENSOR CONSTRUCTORS, DESTRUCTOR
// ____________________________________________________________________________


space_T::space_T() { rank = 0; vx = NULL; }

	// Input		none :
	// Return		SphT : Empty spatial tensor
	// Note 		     : Euler angles and PAS are zeroed
	//			       in class coordinate



space_T::space_T(const space_T& SphT)

	// Input		SphT	: Spatial tensor
	// Return		SphT1	: Spatial tensor (this) which is
	//				  a duplicate of input spatial tensor
	// Note				: Remember, vx is a pointer to a an
	//				  array of row vector pointers, so
	//				  vx[n] is a pointer to the nth row
	//				  vector and *vx[n] is the nth vector.

  {
  rank = SphT.rank;			// Copy Tensor rank
  vx = NULL;				// Delete any exiting vector
  EA = SphT.EA;				// Copy current Euler angles
  PAS_EA = SphT.PAS_EA;			// Copy current PAS Euler angles
  PAS = SphT.PAS;			// Copy current PAS
  if(!SphT.vx) return;			// Done if SphT is empty 
  vx = new row_vector*[rank+1];		// Set array of [A0, A1, ...., Al]
  for(int l=0; l<=rank; l++)		// Loop over l components and copy
    {					// the { Alm } per each l present 
    if(SphT.vx[l])			// (only if l terms are in SphT)
      vx[l]=new row_vector(*(SphT.vx[l]));
    else vx[l] = NULL;
    }
  }


space_T::space_T(const space_T& SphT, int l)

	// Input		SphT : Spatial tensor
	//			l    : tensor rank
	// Return		SphT1: Irreducible spatial tensor from input
        //	                       spatial tensor rank l

  {
  vx = NULL;					// No components at all
  rank = 0;					// No rank at all
  if(!SphT.vx) return;				// Return if SphT empty	
  if((SphT.rank <= l) && SphT.vx[l])
    {
    rank = l;
    vx = new row_vector *[rank+1];
    for (l=0; l<rank; l++)
    vx[l] = NULL;
    l = rank;
    vx[l] = new row_vector(*(SphT.vx[l]));
    }
  }


space_T::space_T(const SinglePar& par)

	// Input		par   : A single GAMMA parameter
	//  			pset  : Parameter set
	// Output		pset  : Parameter set with
	//			        only spherical tensor
        // Note                      : This function assumes that a
        //                             GAMMA parameter for a spatial tensor
        //                             is of type 4 and stored as a string
        //                             as "rank, ( #, #, # ), (#, #, # )".
	//			       All I/O with GAMMA parameters must
	//			       follow this format!      
// sosi - this is only good for rank 2 & symmetric so far!

  {
  int trank=0;
  double tiso=0, tdelz=0, teta=0;
  double talpha=0, tbeta=0, tgamma=0;

  if(par.type() != 4)			// Insure parameter type is tensor
    space_T_fatality(4);

  string val1, val2;
  val1 = par.data();			// Get the parameter data
  cutWhite(val1);			// Remove any initial blanks
  val2 = cutInt(val1, 0);		// Cut out an integer - rank
  trank = atoi(val2.c_str());		// Get the tensor rank

  cutBlksXBlks(val1, ",");		// Cut out any blanks-,-blanks
  cutBlksXBlks(val1, "(");		// Cut out any blanks-(-blanks
  val2 = cutDouble(val1);		// Cut out a double value - Aiso
  tiso = atof(val2.c_str());
  cutBlksXBlks(val1, ",");		// Cut out any blanks-,-blanks
  val2 = cutDouble(val1);		// Cut out a double value - delz
  tdelz = atof(val2.c_str());
  cutBlksXBlks(val1, ",");		// Cut out any blanks-,-blanks
  val2 = cutDouble(val1);		// Cut out a double value - eta
  teta = atof(val2.c_str());

  cutBlksXBlks(val1, ")");		// Cut out any blanks-)-blanks
  cutBlksXBlks(val1, ",");		// Cut out any blanks-,-blanks
  cutBlksXBlks(val1, "(");		// Cut out any blanks-)-blanks
  val2 = cutDouble(val1);		// Cut out a double value - alpha
  talpha = atof(val2.c_str());
  cutBlksXBlks(val1, ",");		// Cut out any blanks-,-blanks
  val2 = cutDouble(val1);		// Cut out a double value - beta
  tbeta = atof(val2.c_str());			
  cutBlksXBlks(val1, ",");		// Cut out any blanks-,-blanks
  val2 = cutDouble(val1);		// Cut out a double value - gamma
  tgamma = atof(val2.c_str());

  rank = trank;				// Tensor rank
  vx = new row_vector *[rank+1];	// Create vectors for tensor components
  if(rank == 2)				// Do this for a rank 2 tensor
    {
    PAS.xyz(tiso, tdelz, teta);		// Set up Principle Axis System
    PAS_EA.xyz(talpha, tbeta, tgamma);
    if(tiso != 0)			// A  exists only if trace exists
      {					//  0
      vx[0] = new row_vector(1);
      (*(vx[0])).put(A2(0,0,tiso,tdelz,teta), 0);
      }
    else
      vx[0] = NULL;
    vx[1] = NULL;			// A  not present for symmetric tensor
    					//  1

    if(tdelz)				// A  not present for spherically
      {					//  2	         symmetric tensor
      vx[2] = new row_vector(5);
      for(int m=0; m<5; m++)
        (*(vx[2])).put(A2(2,2-m,tiso,tdelz,teta), m);
      }
    else
      vx[2] = NULL;
    }
  else
    space_T_fatality(5);
  return;				// Return SphT
  }


space_T::~space_T()

	// Input		SphT : Sperical Tensor
	// Return		     : none, deletes current Tensor

  {
  if(!vx) return;			// If no vectors, do nothing
  for(int l=0; l<=rank; l++)		// Loop over all l components
    if(vx[l]) delete vx[l];		// and delete any Al vector
  delete [] vx;				// Now delete array of vector
  vx = NULL;				// pointers and insure its NULL
  }


// ______________________________________________________________________
//                      SPACE TENSOR UNARY OPERATIONS
// ______________________________________________________________________

// ********************** Tensors with Tensors **************************

space_T& space_T::operator = (const space_T &SphT)

	// Input		SphT : spatial tensor
	// Return		SphT1: spatial tensor equivalent to the
        //	                       input spatial tensor

  {
  int l = 0;
  if(this == &SphT) return *this;		// Check if already equal
  if(SphT.vx)				// Input Tensor exists
    {
    if(vx)				// Initially zero output Tensor
      {
      for (l=0; l<=rank; l++)
        if(vx[l])
          delete vx[l];
      delete [] vx;
      }
    rank = SphT.rank;			// Copy Tensor rank
    vx = new row_vector *[SphT.rank+1];	
    for (l=0; l<=SphT.rank; l++)	// Copy all ranks
      {
      if(SphT.vx[l])
        vx[l] = new row_vector(*(SphT.vx[l]));
      else
        vx[l] = NULL;
      }
    if(SphT.rank == 2)			// Copy PAS if rank 2
      {
      PAS = SphT.PAS;
      PAS_EA = SphT.PAS_EA;
      EA = SphT.EA;
      }
    }
  else					// Input Tensor NULL
    {					// Insure output Tensor NULL
    if(vx)
      {
      for (l=0; l<=rank; l++)
        if(vx[l])
          delete vx[l];
      delete [] vx;
      vx = NULL;
      }
    }
  return *this;
  }


space_T operator + (const space_T& SphT1, const space_T& SphT2)
//return SphT3;

	// Input		SphT1 : spatial tensor
	// 			SphT2 : spatial tensor
	// Return		SphT3 : spatial tensor which is the addition
	//			        of the two input spatial tensors
	//		                SphT3 = SphT1 + SphT2	

{
space_T SphT3;
if(SphT1.rank > SphT2.rank) 
  {
  SphT3 = SphT1;
  SphT3 += SphT2;
  }
else 
  {
  SphT3 = SphT2;
  SphT3 += SphT1;
  }
return SphT3;					// Return SphT3
}


void operator += (space_T& SphT1, const space_T& SphT2)

	// Input		SphT1 : spatial tensor
	// 			SphT2 : spatial tensor
	// Return		SphT1 : spatial tensor which is the addition
	//				of the two input spatial tensors
	//		               	SphT1 = SphT1 + SphT2	
	// Note		              :	The rank of SphT1 will equal the
	//				higher of SphT1 or SphT2 input	

{
if(SphT2.vx)				// Add only if SphT2 is not NULL
  {
  if(SphT1.vx)				// Add if SphT1 is not NULL
    {					// SphT1 must be expanded
    if(SphT2.rank > SphT1.rank)		// If rank of SphT2 is larger the
      {					// SphT1 must be expanded
      row_vector **vxtemp;
      vxtemp = new row_vector *[SphT2.rank+1];
      int i=0;
      for(i=0; i<=SphT2.rank; i++)
        vxtemp[i] = NULL;
      for(i=0; i<=SphT1.rank; i++)
        vxtemp[i] = SphT1.vx[i];
      delete [] SphT1.vx;
      SphT1.vx = vxtemp;
      SphT1.rank = SphT2.rank;
      }
    for(int l=0; l<=SphT2.rank; l++)	// Add all irreducible tensors
      if(SphT2.vx[l])			// Add only if rank exists
        {
        if(SphT1.vx[l])
          *(SphT1.vx[l]) += *(SphT2.vx[l]);
        else
          SphT1.vx[l] = new row_vector(*SphT2.vx[l]);
// sosi - switched to the above in alpha version
//          {
//          SphT1.vx[l] = new row_vector(2*l+1);
//          *(SphT1.vx[l]) = *(SphT2.vx[l]);
//          }
        }
    }
  else					// Set to SphT2 if SphT1 is NULL
    SphT1 = SphT2;
  }
return;
}


space_T operator - (const space_T& SphT1, const space_T& SphT2)

	// Input		SphT1 : spatial tensor
	// 			SphT2 : spatial tensor
	// Return		SphT3 : spatial tensor which is the subtraction
	//			        of the two input spatial tensors
	//		                SphT3 = SphT1 - SphT2	

  { space_T SphT3(SphT1); SphT3 -= SphT2; return SphT3; }



void operator -= (space_T& SphT1, const space_T& SphT2)

	// Input		SphT1 : spatial tensor
	// 			SphT2 : spatial tensor
	// Return		SphT1 : spatial tensor which is the subtraction
	//				of the two input spatial tensors
	//		               	SphT1 = SphT1 - SphT2	
	// Note		              :	The rank of SphT1 will equal the
	//				higher of SphT1 or SphT2 input	

{
if(SphT2.vx)				// Add only if SphT2 is not NULL
  {
  if(SphT1.vx)				// Add if SphT1 is not NULL
    {					// SphT1 must be expanded
    if(SphT2.rank > SphT1.rank)		// If rank of SphT2 is larger the
      {					// SphT1 must be expanded
      row_vector **vxtemp;
      vxtemp = new row_vector *[SphT2.rank+1];
      int i=0;
      for(i=0; i<=SphT2.rank; i++)
        vxtemp[i] = NULL;
      for(i=0; i<=SphT1.rank; i++)
        vxtemp[i] = SphT1.vx[i];
      delete [] SphT1.vx;
      SphT1.vx = vxtemp;
      SphT1.rank = SphT2.rank;
      }
    for(int l=0; l<=SphT2.rank; l++)	// Add all irreducible tensors
      if(SphT2.vx[l])			// Add only if rank exists
        {
        if(SphT1.vx[l])
          *(SphT1.vx[l]) -= *(SphT2.vx[l]);
        else
          SphT1.vx[l] = new row_vector((*SphT2.vx[l]).matrix::operator-());
// sosi - switched to the above in alpha version
//          {
//          SphT1.vx[l] = new row_vector(2*l+1);
//          *(SphT1.vx[l]) = -(*(SphT2.vx[l]));
//          }
        }
    }
  else SphT1 = -SphT2; 			// Set to SphT2 if SphT1 is NULL
  }
return;
}

space_T operator - (const space_T &SphT)

	// Input		SphT : spatial tensor
	// Return		SphT1 : spatial tensor which is the negated
	//			        input spatial tensor
	//		                SphT3 = - SphT2	

  { return (complex(-1.0)*SphT); }


// ********************** Tensors with Scalars **************************


space_T operator * (const complex &z, const space_T& SphT)
//return SphT1(SphT);

	// Input		z     : complex number
	// 			SphT  : spatial tensor
	// Return		SphT1 : spatial tensor which is the input
	//			        spatial tensor multiplied by the
	//				complex scalar constant
	//		                SphT1 = z * SphT

{
// ???? should add a check for zero z
space_T SphT1(SphT);
if(SphT1.vx)			// Multiply only if SphT2 is not NULL
  for(int i=0; i<=SphT1.rank; i++)
    if(SphT1.vx[i])		// Multiply only if rank is not NULL
      *(SphT1.vx[i]) *= z;
return SphT1;				// Return SphT1
}


space_T operator * (const space_T& SphT, const complex& z)

	// Input		SphT  : spatial tensor
	// 			z     : complex number
	// Return		SphT1 : spatial tensor which is the input
	//			        spatial tensor multiplied by the
	//				complex scalar constant
	//		                SphT1 = SphT * z

{ return (z*SphT); }


// ______________________________________________________________________
//                GENERAL RANK 1 SPHERICAL TENSOR FUNCTIONS
// ______________________________________________________________________


space_T A1(double x, double y, double z)
//return SphT;

	// Input		x     : Cartesian x component
	// 			y     : Cartesian y component
	// 			z     : Cartesian z component
	// Output		SphT  : Rank 1 Spatial Tensor for
	//			        Cartesian vector specified
 
{
  space_T SphT;
  SphT.rank = 1;
  SphT.vx = new row_vector *[2];
  SphT.vx[0] = NULL;
  SphT.vx[1] = new row_vector(3);
  for (int m=0; m<3; m++)
    (*(SphT.vx[1])).put(A1(x, y, z, 1-m), m);
  return SphT;				// Return SphT
}


space_T A1(coord &pt)

	// Input		pt   : Cartesian coordinate 
	// Output		SphT : Rank 1 Spatial Tensor for
	//			       Cartesian coordinate specified
 
{ return A1(pt.x(), pt.y(), pt.z()); }


// No longer supported *** Delete eventually
space_T A1(row_vector &vec)

	// Input		vec  : Cartesian vector
	// Output		SphT : Rank 1 Spatial Tensor for
	//			       Cartesian vector specified
 
{
  if(vec.elements() != 3)
    space_T_fatality(11);
  return A1(Re(vec(0)), Re(vec(1)), Re(vec(2)));
}


// Case l=0 no longer supported *** remove eventually
complex A1(double x, double y, double z, int m, int l)
//return c;

	// Input		x : Cartesian x component
	// 			y : Cartesian y component
	// 			z : Cartesian z component
	// 			m : Tensor component
	// Output		c : Spatial Tensor component for
	//			    Cartesian vector specified
 
{
  complex c;
  switch (l)
    {
    case 0:
      c = A10(m);			// A  from A , A , and A
      break;				//  00      x   y       z
    case 1:
      c = A11(x, y, z, m);		// A  , A  , A    from A , A , and A
      break;				//  11   10   1-1       x   y       z
    default:
      space_T_error(0);
      cout << " 0," << m << "\n";
      break;
    }
  return c;				// Return c
}


// Case l=0 no longer supported *** remove eventually
complex A1(coord &pt, int m, int l)

	// Input		pt  : Cartesian vector
	// 			m   : Tensor component
	// Output		c   : Spatial Tensor component for
	//			      Cartesian coordinate specified

{ return A1(pt.x(), pt.y(), pt.z(), m, l); }


// No longer supported *** Delete eventually
complex A1(row_vector &vec, int m, int l)

	// Input		l   : Tensor component
	// 			m   : Tensor component
	// 			vec : Cartesian vector
	// Output		c   : Spatial Tensor component for
	//			     Cartesian vector specified

{ return A1(Re(vec(0)), Re(vec(1)), Re(vec(2)), m, l); }


// Case l=0 no longer supported *** remove eventually
double A10(int m)

{
  switch (m)
    {
    case 0:
      break;
    default:
      space_T_error(0);
      cout << " 0," << m << "\n";
      break;
    }
  return 1.0;
}


complex A11(double x, double y, double z, int m)
//return c;

{
  complex c;
  switch (m)
    {
					//       [1} 
    case -1:				// A   = |-|{A  - iA }
      c = sqrt(0.5)*complex(x,-y);	//  1-1  [2]  x     y
      break;

    case 0:				// A   = A
      c = z;				//  10    z
      break;
					//         [1} 
    case 1:				// A   = - |-|{A  + iA }
      c = -sqrt(0.5)*complex(x,y);	//  11     [2]  x     y
      break;

    default:
      space_T_error(0);
      cout << " 1," << m << "\n";
      break;
    }
  return c;				// Return c
}


space_T SphA1(complex plus, complex zero, complex minus)
//return SphT;

	// 			plus  : Component A11
	// 			zero  : Component A10
	// 			minus : Component A1-1
	// Output		SphT  : Rank 1 Spatial Tensor for
	//			        spherical vector specified
 
  {
  space_T SphT;
  SphT.rank = 1;
  SphT.vx = new row_vector *[2];
  SphT.vx[0] = NULL;
  SphT.vx[1] = new row_vector(3);
  (*(SphT.vx[1])).put(plus,0);
  (*(SphT.vx[1])).put(zero,1);
  (*(SphT.vx[1])).put(minus,2);
  return SphT;				// Return SphT
  }


space_T SphA1(coord &pt)

	// Input		pt   : Coordinate point (spherical)
	// Output		SphT : Rank 1 Spatial Tensor for
	//			       spherical coordinate specified
 
  { return SphA1(pt.x(), pt.y(), pt.z()); }


// No longer supported *** Delete eventually
space_T SphA1(row_vector &vec)

	// Input		vec  : Spherical vector
	// Output		SphT : Rank 1 Spatial Tensor for
	//			       spherical vector specified
 
  {
  if(vec.elements() != 3)
    space_T_fatality(11);
  return SphA1(vec(0), vec(1), vec(2));
  }


// ____________________________________________________________________________ 
//                  GENERAL RANK 2 SPHERICAL TENSOR FUNCTIONS
// ____________________________________________________________________________ 

// ************************ Whole Tensor Production ***************************

space_T A2(double Aiso, double delzz, double eta,
			              double alpha, double beta, double gamma)
//return SphT;

	// Input		Aiso  : Isotropic component
	//			delzz : Anisotropy parameter
	//			eta   : Asymmetry parameter 
	//			alpha : Euler angle (degrees)
	//			beta  : Euler angle (degrees)
	//			gamma : Euler angle (degrees)
	// Output		SphT  : Rank 2 Spatial Tensor
 
  {
  space_T SphT;
  SphT.rank = 2;			// Set rank to 2
  SphT.vx = new row_vector *[3];	// Set up components 0, 1, & 2
  SphT.PAS.xyz(Aiso, delzz, eta);	// Set up Principle Axis System
  SphT.PAS_EA.xyz(alpha,beta,gamma);

  if(Aiso != 0.)			// A  exists only if trace exists
    {					//  0
    SphT.vx[0] = new row_vector(1);
    (*(SphT.vx[0])).put(A2(0,0, Aiso, delzz, eta), 0);
    }
  else
    SphT.vx[0] = NULL;

  SphT.vx[1] = NULL;			// A  not present for symmetric tensor
    					//  1

  if(delzz)				// A  not present for spherically
    {					//  2			symmetric tensor
    SphT.vx[2] = new row_vector(5);
    for (int m=0; m<5; m++)
      (*(SphT.vx[2])).put(A2(2,2-m, Aiso, delzz, eta),m);
    }
  else
    SphT.vx[2] = NULL;
  return SphT;				// Return SphT
  }


space_T A2(coord &Tcomps) { return A2(Tcomps.x(), Tcomps.y(), Tcomps.z()); }

	// Input		Tcomps : Three tensor PAS components
	// Output		SphT   : Rank 2 Spatial Tensor
 


space_T A2(coord &Tcomps, coord &Tangles)

	// Input		Tcomps : Three tensor PAS components
	//                               Aiso, delzz, eta
	//			Tangles: Three tensor Euler angles
	//                               alpha, beta, gamma (degrees)
	// Output		SphT   : Rank 2 Spatial Tensor
 
  {
  return A2(Tcomps.x(), Tcomps.y(), Tcomps.z(),
			Tangles.x(), Tangles.y(), Tangles.z());
  }



space_T A2(const matrix& mx, double prec)
//return SphT;	

	// Input		mx	: Rank 2 Cartesian tensor
	//			prec    : Matrix precision
	// Output		SphT	: Rank 2 Spatial Tensor
	// Note				: The orientation of mx will be		
	//				  retained in the returned SphT 
	// Note				: SphT is maintained in irreducible
	//				  spherical format

// sosi - why doesn't this set the PAS components?
//       when you set the tensor iso value it sets both A00 & Aiso!
//       I set this to make it so 10/28/97, but should be consistent!
 
  {
  space_T SphT;
  if(mx.cols() != 3 || mx.rows() != 3) 		// Die if array not 3x3
    space_T_fatality(12);
  SphT.rank = 2;				// Set the rank to 2
  int span = 0;					// The m Span per l value
  SphT.vx = new row_vector *[3];		// Construct over 0, 1, 2
  for(int l=0; l<3; l++)			// Loop the ranks
    {
    span = 2*l+1;				//	The m span this l
    SphT.vx[l] = new row_vector(span);		//	Construct array
    for(int m=0; m<span; m++)			//	Fill in components
      (*(SphT.vx[l])).put(A2(l,l-m, mx), m);
    }
  if(norm(trace(mx)) < prec)			// Zero A    if no mx trace
    {						//       0,0
    delete SphT.vx[0];
    SphT.vx[0] = NULL;
    } 
  if(mx.test_type(h_matrix_type,1.e-10)
                             == h_matrix_type)	// Zero A   if symmetric mx
    {						//       1,m
    delete SphT.vx[1];
    SphT.vx[1] = NULL;
    } 
  SphT.updatePAS();				// Fill up PAS values
  return SphT;					// Return SphT
  }

// ************* Specific Alm Tensor Component Production ***************

complex A2(int l, int m, double Aiso, double delzz, double eta)

	// Input		l     : component index
	//			m     : component index
	// 			Aiso  : isoropic component
	//			delzz : anisotropy parameter
	//			eta   : asymmetry parameter 
	// Output		SOp   : rank 2 Spatial Tensor component
	//			        Alm

  {
  complex z;
  switch (l)
    {
    case 0: z = A20(m, Aiso, delzz, eta); break;
    case 1: z = A21(m, Aiso, delzz, eta); break;
    case 2: z = A22(m, Aiso, delzz, eta); break;
    default:
      space_T_error(0);
      cout << " " << l << "0," << m << "\n";
      space_T_fatality(10);
      break;
    }
  return z;				// Return z
  }


complex A2(int l, int m, const matrix &mx)

	// Input		l     : component index
	//			m     : component index
	// 			mx    : Cartesian rank 2 tensor
	// Output		SOp   : rank 2 Spatial Tensor component
	//			        Alm

  {
  complex z;
  switch (l)
    {
    case 0:
      z = A20(m, mx);
      break;
    case 1:
      z = A21(m, mx);
      break;
    case 2:
      z = A22(m, mx);
      break;
    default:
      space_T_error(0);
      cout << " " << l << "0," << m << "\n";
      space_T_fatality(10);
      break;
    }
  return z;				// Return z
  }


complex A20(int m, double Aiso, double delzz, double eta)

	// Input		m     : component index
	// 			Aiso  : isoropic component
	//			delzz : anisotropy parameter
	//			eta   : asymmetry parameter 
	// Output		SOp   : rank 2 Spatial Tensor component A
	//								 00
 
  {
  complex z;
  switch (m)
    {
    case 0:
      z = -1.0*sqrt(3.0)*Aiso;
      break;
    default:
      space_T_error(0);
      cout << " 0," << m << "\n";
      break;
    }
  return z;				// Return z
  delzz = 0;				// Compiler likes this used
  eta = 0;				// Compiler likes this used
  }


complex A20(int m, const matrix &mx)

	// Input		m     : component index
	// 			mx    : Cartesian rank 2 tensor
	// Output		SOp   : rank 2 Spatial Tensor component A
	//								 00
 
  {
  complex z;
  switch (m)
    {
    case 0:
      z = -1.0*trace(mx)/sqrt(3.0); break;
    default:
      space_T_error(0);
      cout << " 0," << m << "\n";
      break;
    }
  return z;				// Return z
  }


complex A21(int m, double Aiso, double delzz, double eta)
// sosix - what the heck is this supposed to be?!

	// Input		m     : component index
	// 			Aiso  : isoropic component
	//			delzz : anisotropy parameter
	//			eta   : asymmetry parameter 
	// Output		SOp   : rank 2 Spatial Tensor component A
	//								 1m

  {
  switch (m)
    {
    case -1: return complex0; break;
    case 0:  return complex0; break;
    case 1:  return complex0; break;
    default:
      space_T_error(0);
      cout << " 1," << m << "\n";
      break;
    }
  return complex0;			// Return z
  Aiso = 0;				// Compiler likes this used
  delzz = 0;				// Compiler likes this used
  eta = 0;				// Compiler likes this used
  }


complex A21(int m, const matrix &mx)

	// Input		m     : component index
	// 			mx    : Cartesian rank 2 tensor
	// Output		SOp   : rank 2 Spatial Tensor component A
	//								 1m

  {
  complex z;
  switch (m)
    {
    case -1:
      z = -0.5*((mx.get(2,0)-mx.get(0,2)) - complexi*(mx.get(2,1)-mx.get(1,2)));
      break;
    case 0:
      z = (-sqrt(0.5))*complexi*(mx.get(0,1)-mx.get(1,0));
      break;
    case 1:
      z = -0.5*((mx.get(2,0)-mx.get(0,2)) + complexi*(mx.get(2,1)-mx.get(1,2)));
      break;
    default:
      space_T_error(0);
      cout << " 1," << m << "\n";
      break;
    }
  return z;
  }


complex A22(int m, double Aiso, double delzz, double eta)

	// Input		m     : component index
	// 			Aiso  : isoropic component
	//			delzz : anisotropy parameter
	//			eta   : asymmetry parameter 
	// Output		SOp   : rank 2 Spatial Tensor component A
	//								 2m
  {
  complex z;
  switch (m)
    {
					//        [1.0]
    case -2:				// A    = |---| * del  * eta
      z = 0.5*delzz*eta;		//  2-2   [2.0]      zz
      break;

    case -1:				// A    = 0
      z = 0.0;				//  2-1
      break;
					//       [3.0]1/2
    case 0:				// A   = |---|   * del
      z = sqrt(1.5)*delzz;		//  20   [2.0]        zz
      break;

    case 1:				// A   = 0
      z = 0.0;				//  21
      break;
					//       [1.0]
    case 2:				// A   = |---| * del  * eta
      z = 0.5*delzz*eta;		//  22   [2.0]      zz
      break;

    default:
      space_T_error(0);
      cout << " 2," << m << "\n";
      break;
    }
  return z;				// Return z
  Aiso = 0;				// Compiler likes this used
  }


complex A22(int m, const matrix &mx)
//return z;

	// Input		m     : Component index
	// 			mx    : Cartesian rank 2 tensor
	// Output		SOp   : Rank 2 Spatial Tensor component A
	//								 2m

  {
  complex z;
  switch (m)
    {
    case -2:
      z = 0.5*(mx.get(0,0) - mx.get(1,1) - complexi*(mx.get(0,1) + mx.get(1,0)));
      break;
    case -1:
      z = 0.5*(mx.get(0,2) + mx.get(2,0) - complexi*(mx.get(1,2) + mx.get(2,1)));
      break;
    case 0:
      z = (2.0*mx.get(2,2) - mx.get(1,1) - mx.get(0,0))/sqrt(6.0);
      break;
    case 1:
      z = -0.5*(mx.get(0,2) + mx.get(2,0) + complexi*(mx.get(1,2) + mx.get(2,1)));
      break;
    case 2:
      z = 0.5*(mx.get(0,0) - mx.get(1,1) + complexi*(mx.get(0,1) + mx.get(1,0)));
      break;
    default:
      space_T_error(0);
      cout << " 2," << m << "\n";
      break;
    }
  return z;			// Return z
  }

// ****************** Principle Axis System Access ********************

coord space_T::PASys() const

	// Input		SphT  : Spherical tensor (this)
	// Output		PAS   : Principal Axis System Components
	// Note			      : Only for rank 2 tensors
	// Note			      : Function named soas not to conflict
	//			   	with structure name
  {
  if(rank !=2) space_T_fatality(21);	//Can't get PAS if rank != 2
  return PAS;				// Return Principle Axes Values
  }


coord space_T::PASys_EA() const

	// Input		SphT  : Spherical tensor (this)
	// Output		pt    : Euler Angle set for orientation
	//				of the principle Axes
	// Note			      : Function named soas not to conflict
	//			   	with structure name
  {
  if(rank !=2) space_T_fatality(21);	//Can't get PAS if rank != 2
  return PAS_EA;		// Return Principle Axes Values
  }


double space_T::iso() const

	// Input		SphT  : Spherical tensor (this)
	// Output		Aiso  : Isotropic component of PAS
	// Note			      : Only for rank 2 tensors
{
// ??? check rank here first
  return PAS.x();		// Return Aiso
}


void space_T::iso(double Aiso)

	// Input		SphT  : Spherical tensor (this)
	// 			Aiso  : Isotropic tensor component
	// Output		none  : SphT isotropic component set
	//				to value Aiso

  {
  if(!rank) space_T_fatality(20);	// Can't set rank if no tensor
  PAS.x(Aiso);				// Set the PAS Aiso value
  if(!Aiso)				// If setting Aiso to zero
    {					// then we must make sure that
    if(vx[0])				// the A   component is gone
      {					//      00
      delete vx[0];
      vx[0] = NULL;
      }
    }
  else
    {
    if(!vx[0])				// Set the tensor A   value
      vx[0] = new row_vector(1);	//                 00
//  (*(SphT.vx[0]))(0) =
    (*(vx[0])).put(-1.0*sqrt(3.0)*Aiso, 0); 
    }
  return;	
  }


double space_T::delz() const

	// Input		SphT  : Spherical tensor (this)
	// Output		delz  : Delta z component of PAS
	// Note			      : Only for rank 2 tensors
{
// ??? check rank here first
  return PAS.y();		// Return delzz
}


void space_T::delz(double delzz)

	// Input		SphT  : Spherical tensor (this)
	// 			delz  : Delta z component of PAS
	// Output		      : Delta z component of PAS set to delz
	// Note			      : Only for rank 2 tensors
  {
// ??? check rank here first
// 7/23/92 - Not a complete function, must reset tensor - sosi
  PAS.y(delzz);			// Set the PAS delz value
  return;
  }


double space_T::eta() const

	// Input		SphT  : Spherical tensor (this)
	// Output		eta   : eta component of PAS
	// Note			      : Only for rank 2 tensors
{
// ??? check rank here first
  return PAS.z();		// Return eta
}


void space_T::eta(double ETA)

	// Input		SphT  : Spherical tensor (this)
	// 			ETA   : The eta component of PAS
	// Output		      : The eta component of PAS set to ETA
	// Note			      : Only for rank 2 tensors
  {
// ??? check rank here first
// 7/23/92 - Not a complete function, must reset tensor - sosi
  PAS.z(ETA);			// Set the PAS delz value
  return;
  }


double space_T::alpha() const

	// Input		SphT  : Spherical tensor (this)
	// Output		alpha : alpha Euler angle of PAS
	// Note			      : Only for rank 2 tensors
{
// ??? check rank here first
  return PAS_EA.x();	// Return alpha
}


double space_T::beta() const

	// Input		SphT  : Spherical tensor (this)
	// Output		beta  : beta Euler angle of PAS
	// Note			      : Only for rank 2 tensors
{
// ??? check rank here first
  return PAS_EA.y();	// Return beta
}


double space_T::gamma() const

	// Input		SphT  : Spherical tensor (this)
	// Output		gamma : gamma Euler angle of PAS
	// Note			      : Only for rank 2 tensors
{
// ??? check rank here first
  return PAS_EA.z();	// Return gamma
}


// ____________________________________________________________________________ 
//                          SPACE TENSOR FUNCTIONS
// ____________________________________________________________________________ 


int space_T::exists() const

	// Input		SphT : Spatial tensor (this)
	// Output		TF   : Returns TRUE if SphT is not NULL
 
  {
  for(int i=0; i<rank; i++)
    if(vx[i]) return(1);
  if(iso())   return(1);
  if(delz())  return(1);
  if(eta())   return(1);
  return 0;
  }


int space_T::exists(int l) const

	// Input		SphT : Spatial tensor (this)
	//			l    : Rank
	// Output		TF   : Returns if SphT has an l
	//			       component
 
  {
  if(l>=0 && l<=rank) if(vx[l]) return(1);
  return 0;
  }

int space_T::Rank() const { return rank; }

	// Input		SphT : Spatial tensor (this)
	// Output		l    : Overall tensor rank
	// Note			     : Function name needs uppercase R
	//			       to distinguish it from structure
 
// **************************************************************************** 
//                       Spatial Tensor Component Functions
// **************************************************************************** 

complex T_comp(const space_T &SphT, int l, int m)
  { 
  SphT.SphTerror(10, "T_comp", 1);		// Deprecated function
  SphT.SphTerror(11, "component");		// Use component instead
  return SphT.component(l,m);			// Now use member function
  }

complex space_T::component(int l, int m) const

	// Input		SphT : Spatial Tensor (this)
	// 			l    : momentum index (e.g. 0, 1, 2...)
	// 			m    : momentum index (e.g. -2, -1, 0, 1, 2)
	// Output		z    : complex number, the l,m
	//			       component of Spatial Tensor SphT
 
  {
  if(l>rank || l<0)
    {
    SphTerror(20, Gdec(l), 1);		// Bad rank for requested component
    SphTerror(20, 1);			// Tell current tensor rank
    SphTerror(21, 1);			// Tell current component rank range
    SphTerror(23, 0);			// Tell we're returning zero
    return complex0;
    }
  else if(abs(m) > l)
    {
    SphTerror(21, Gdec(m), 1);		// Bad m for requested component
    SphTerror(20, 1);			// Tell current tensor rank
    SphTerror(22, 1);			// Tell current component m range
    SphTerror(23, 0);			// Tell we're returning zero
    return complex0;
    }
  else if(vx[l]) return (*(vx[l]))(l-m);
  return complex0;;
  }



double space_T::Ccomponent(int r, int c) const

	// Input		SphT : Spatial Tensor (this)
	// 			r    : row index
	// 			c    : column index
	// Output		Arc  : Cartesian tensor component

  {
  double val = 0.0;
  complex z;
  if(!vx) return val;
  switch(rank)
    {
    case 0:
      z = component(0,0);
      return zRe(z);
      break;
    case 1:
      {
      if(!vx[1] && !vx[0]) return val;
      switch(r)
        {
        case 0:	
          z = sqrt(0.5)*(-component(1,1) + component(1,-1));
          z += component(0, 0);
          break;
        case 1:
          z = sqrt(0.5)*(component(1,1) + component(1,-1));
          z = complexi*z + component(0,0);
          break;
        case 2:
          z = component(1,0) + component(0,0);
          break;
        default:
          z = complex0;
          break;
        }
      val = Re(z);
      if(fabs(val) < 1.e-10) val = 0.0;
      return val;
      }
      break;
    case 2:
      {
      if(!vx[2] && vx[1] && vx[0]) return val;
      switch(r)
        {
        case 0:					// Axx, Axy, or Axz
          {
          if(c==0)
            {
            z = 0.5 * (component(2,2) + component(2,-2));
            z = z - component(2, 0)/sqrt(6.0) - component(0, 0)/sqrt(3.0);
            }
          else if(c==1)
            {
            z = -0.5 * (component(2, 2) - component(2, -2));
            z += component(1, 0)/sqrt(2.0);
            z *= complexi;
            }
          else if(c==2)
            {
            z = 0.5 * (component(1, 1) + component(1, -1));
            z = z - 0.5 * (component(2, 1) - component(2, -1));
            }
          else z = complex0;
          }
          break;
        case 1:					// Ayx, Ayy, or Ayz
          {
          if(c==0)
            {
            z = -0.5 * (component(2, 2) - component(2, -2));
            z = z - component(1, 0)/sqrt(2.0);
            z = z*complex(0,1.0);
            }
          else if(c==1)
            {
            z = -0.5 * (component(2, 2) + component(2, -2));
            z = z - component(2, 0)/sqrt(6.0) - component(0, 0)/sqrt(3.0);
            }
          else if(c==2)
            {
            z = 0.5 * (component(2, 1) + component(2, -1));
            z = z - 0.5 * (component(1, 1) - component(1, -1));
            z *= complexi;
            }
          else z = complex0;
          }
          break;
        case 2:					// Azx, Azy, or Azz
          {
          if(c==0)
            {
            z = -0.5 * (component(1, 1) + component(1, -1));
            z = z - 0.5 * (component(2, 1) - component(2, -1));
            }
          else if(c==1)
            {
            z = 0.5 * (component(2, 1) + component(2, -1));
            z = z + 0.5 * (component(1, 1) - component(1, -1));
            z = z * complex(0,1.0);
            }
          else if(c==2)
            { z = sqrt(2.0/3.0)*component(2,0) - sqrt(1.0/3.0)* component(0,0); }
          else z = complex0;
          }
          break;
        default:
          z = complex0;
          break;
        }
      val = Re(z);
      if(fabs(val) < 1.e-10) val = 0.0;
      return val;
      break;
      }
    default:
      return val;
      break;
    }
  return val;
  }


space_T T_mult(const space_T &SphT1, const space_T &SphT2)
//return SphT;

	// Input		SphT1 : Spatial Tensor
	// 			SphT2 : Spatial Tensor
	// Output		SphT  : Spatial Tensor which is the product
	//			        of the two input tensors
	//				SphT = SphT1 x SphT2
 
{
space_T SphT;
if(SphT1.vx && SphT2.vx)		// Return NULL if either is NULL
  {
  for(int l1=0; l1<=SphT1.rank; l1++)   // Sum over all ranks of both
    for(int l2=0; l2<=SphT2.rank; l2++) // input tensors
      if(SphT1.vx[l1] && SphT2.vx[l2])
        SphT += T_mult(SphT1, l1, SphT2, l2);
  }
return SphT;					// Return SphT
}


space_T T_mult(const space_T &SphT1, int l1, const space_T &SphT2, int l2)
//return SphT;

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
 
{
space_T SphT;
if ((l1 > SphT1.rank) || (l2 > SphT2.rank))	// Check components exist
  {
  cout << "\nSpace_T: A  x A  Irreducible Product on Tensors A  and A ."
       << "\nSpace_T:  " << l1 << "    " << l2 << "\t\t"  << SphT1.rank
       << "      " << SphT2.rank;
  space_T_fatality(9);
  }
if ((l1 < 0) || (l2 < 0))			// Check for negative rank
  {
  space_T_error(8);
  space_T_fatality(9);
  }
if(SphT1.vx && SphT2.vx)			// Return NULL if either is NULL
  if(SphT1.vx[l1] && SphT2.vx[l2])
    {
    SphT.rank = l1 + l2;
    SphT.vx = new row_vector *[SphT.rank+1];
    int L=0;
    for(L=0; L<=SphT.rank; L++)			// Initialize all irreducible comps 
      SphT.vx[L] = NULL;
    for(L=0; L<=SphT.rank; L++)			// Compute all irreducible comps 
      {
      if(L < abs(l1 - l2)) 
        SphT.vx[L] = NULL;
      else 
        SphT += T_mult(SphT1, l1, SphT2, l2, L);
      }
    }
return SphT;						// Return SphT
}


space_T T_mult(const space_T &SphT1, int l1, const space_T &SphT2, int l2, int L)
//return SphT;

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
 
{
space_T SphT;
if ((l1 > SphT1.rank) || (l2 > SphT2.rank))	// Check components exist
  {
  cout << "\nSpace_T: A  x A  Irreducible Product on Tensors A  and A ."
       << "\nSpace_T:  " << l1 << "    " << l2 << "\t\t"  << SphT1.rank
       << "      " << SphT2.rank;
  space_T_fatality(9);
  }
if ((l1 < 0) || (l2 < 0))			// Check for negative rank
  {
  space_T_error(8);
  space_T_fatality(9);
  }
if ((L <= l1+l2) && (L >= abs(l1-l2)))		// Return NULL if L outside range
  {
  SphT.rank = L;
  int span = 2*L+1;
  SphT.vx = new row_vector *[SphT.rank+1];
  for(int l=0; l<L; l++)			// Set all ranks to NULL
    SphT.vx[l] = NULL;
  SphT.vx[L] = new row_vector(span);		// Irreducible rank L output
  for (int M=0; M<span; M++)
    (*(SphT.vx[L])).put(T_mult(SphT1, l1, SphT2, l2, L, L-M), M);
  }
return SphT;						// Return SphT
}


complex T_mult(const space_T &SphT1, int l1,
				 const space_T &SphT2, int l2, int L, int M)

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
 
  {
  complex z = complex0;
  double norm = 0.0;
  double prefact = 0.0;
  int span1 = 2*l1+1;
  int span2 = 2*l2+1;
  int m1 = 0;
  int m2 = 0;
  norm = sqrt(double(2*L+1));
  prefact = pow(double(-1.0), double(abs(2*l2+L-M)));
  for(int M1=0; M1<span1; M1++)
    for(int M2=0; M2<span2; M2++)
      if((l1-M1+l2-M2) == M)
        {
        m1 = l1-M1;
        m2 = l2-M2;
        z += SphT1.component(l1, m1) * SphT2.component(l2, m2)
			 	* Wigner_3j(L, l1, l2, M, -m1, -m2);
        }
  z *= prefact*norm;
  return z;
  }


// ____________________________________________________________________________ 
//                        SPACE TENSOR ROTATION FUNCTIONS
// ____________________________________________________________________________

/* These functions rotate a spatial tensor or spatial tensor component by user
   specified input Euler angles. 

	   Input		SphT  : Spatial Tensor (this)
	   			alpha : Euler angle alpha 
	   			beta  : Euler angle beta
	   			gamma : Euler angle gamma
				l,m   : Component indices
	   Output		SphT  : Spatial tensor (component) which is the
					input tensor (component) rotated from 
				        its position by specified Euler angles

                                 +/- l
                                  ---    l
       A   (alpha, beta, gamma) = \     D    (alpha,beta,gamma) A    
        l,m                       /      m',m                    l,m'
                                  ---
                                   m'
           l
   Above, D    (alpha,beta,gamma) is the Wigner rotation matrix element for the
           m',m                l
   rotation (also written <m'|D (alpha,beta,gamma)|m>. The rotated elements are
   linear combinations of the un-rotated elements within the same rank (same l).

   Note that the non-member "friend" functions T_rot are now deprecated
   and destined to be removed in further GAMMA versions.                     */

complex T_rot(space_T &SphT1,int l,int m,double alpha,double beta,double gamma)
  { 
  SphT1.SphTerror(10, "T_rot", 1);			// Deprecated function
  SphT1.SphTerror(11, "rotate");			// Use rotate instead
  return SphT1.rotate(l,m,alpha,beta,gamma);
  }

space_T T_rot(space_T &SphT1, double alpha, double beta, double gamma)
  {
  SphT1.SphTerror(10, "T_rot", 1);			// Deprecated function
  SphT1.SphTerror(11, "rotate", 0);			// Use rotate instead
  return SphT1.rotate(alpha,beta,gamma); 
  }

void T_rot(int num, space_T* SphT, space_T* SphTrot,
                                        double alpha, double beta, double gamma)

	// Input		num   : Number of tensors
	//			SphT  : Array of spatial tensors
	//			SphTr : Array of rotated spatial tensors
	// 			alpha : Euler Angle (degrees)
	// 			beta  : Euler Angle (degrees)
	// 			gamma : Euler Angle (degrees)
	// Output		void  : The array SphTr is filled 
	//			        with rotated tensors in array SphT

  {
  int lmax = 0;
  for(int i=0; i<num; i++)			// Loop over input tensors
    if(SphT[i].rank > lmax) lmax=SphT[i].rank;	// Get largest rank
  matrix* D;
  D = new matrix[lmax+1];
//  matrix D[lmax+1];				// Array of Wigner matrices
  for(int l=0; l<=lmax; l++)			// Loop over all ranks
    D[l] = DJ(l,alpha,beta,gamma);		// Get Wigner rotation matrices
  T_rot(num, SphT, SphTrot, D);			// Use function overload
  delete [] D;
  return;
  }


void T_rot(int num, space_T* SphT, space_T* SphTrot, matrix* D)

	// Input		num   : Number of tensors
	//			SphT  : Array of spatial tensors
	//			SphTr : Array of rotated spatial tensors
	// 			D     : Array or Wigner rotation matrices
	// Output		void  : The array SphTr is filled 
	//			        with rotated tensors in array SphT
	// Note			      : Assumes SphTrot tensors have been
	//				initialized to be the unrotated
	//				SphT tensors

  {
  complex Alm = 0;
  int span=0, l=0, n=0;
  for(int i=0; i<num; i++)			// Loop over input tensors
    {
    if(SphTrot[i].vx)				// Rotate only if tensor exists
      {
//      for(l=0; l<=SphT[i].rank; l++)		// Loop over all tensor ranks
//      Take l=0 as rotationally invarient, & preset to SphT value
      for(l=1; l<=SphT[i].rank; l++)		// Loop over all tensor ranks
        {
        if(SphTrot[i].vx[l])			// Rotate this rank if it exists
          {
          span = 2*l+1;				// Number components this rank
          for(int m=0; m<span; m++)		// Loop over all components
            {
            Alm=0;
            for(n=0; n<span; n++)
              Alm += (*(SphT[i].vx[l]))(n)
                   * D[l].get(span-1-n, span-1-m);
            (*(SphTrot[i].vx[l])).put(Alm, m);
            }
          }
        }					// Rotate next rank
      }						// Rotate next tensor if exists
    }						// Get next tensor
  return;
  }


space_T space_T::rotate(double alpha, double beta, double gamma) const
  {
  space_T SphT(*this);			// Copy the input tensor
  if(!SphT.vx) return SphT;		// Bail if empty tensor
  for(int l=0; l<=SphT.rank; l++)	// Loop over all ranks
    {
    if(SphT.vx[l])			// If rank l exists then we must
      {					// loop over all m components and
      int span = 2*l+1;			// and rotate them (linear combos)
      for(int m=0; m<span; m++)		// Recall *(SphT.vx[l]) is the lth
        (*(SphT.vx[l])).put(rotate(l, l-m, alpha, beta, gamma), m);
      }
    }
  return SphT;				// Return SphT
  }


space_T space_T::rotate(const coord& EA) const
   { return rotate(EA.x(), EA.y(), EA.z()); }


complex space_T::rotate(int l, int m, double alp, double bet, double gam) const
  {
  complex z(0);					// Begin with zero
  int span = 2*l+1;				// Number of m components
  for(int n=0; n<span; n++)			// Loop current components
    z += (*(vx[l]))(n) * DJ(l, l-n, m, alp, bet, gam);
  return z;
  }


complex space_T::rotate(int l, int m, const coord &EA) const

	// Input		SphT1 : Spatial tensor (this)
	//			l     : Component index
	//			m     : Component index
	// 			EA    : Set of Euler angles
	// Output		Alm   : Spatial Tensor component of
	//			        the input tensor in the coordinate
	//				system rotated by input Euler angles
 
{ return rotate(l, m, EA.x(), EA.y(), EA.z()); }


// ____________________________________________________________________________ 
//                        PARAMETER SET FUNCTIONS
// ____________________________________________________________________________ 


SinglePar space_T::param(const string& pname) const
 
        // Input               SphT  : A spatial tensor (this)
        //                     pname : A parameter name
        // Return              par   : A GAMMA parameter of type space_T 
        //                             with the name pname

  {
  string pstate = "Coordinate Point";
  return param(pname, pstate);
  }
 
 
SinglePar space_T::param(const string& pname, const string& pstate) const
 
        // Input               SphT  : A spatial tensor (this)
        //                     pname : A parameter name
        //                     pstate: A parameter statement
        // Return              par   : A GAMMA parameter of type space_T 
        //                             with the name pname and statement pstate
 
   {
   string pdata;
   pdata = string(Gdec(rank));
   pdata += string(", ( ");
   pdata += Gform("%g", iso());
   pdata += string(", ");
   pdata += Gform("%g", delz());
   pdata += string(", ");
   pdata += Gform("%g", eta());
   pdata += string(" ), ( ");
   pdata += Gform("%g", alpha());
   pdata += string(", ");
   pdata += Gform("%g", beta());
   pdata += string(", ");
   pdata += Gform("%g", gamma());
   pdata += string(" )");
   // pdata += string(", ( ") + string(dtoa(iso()));
   //pdata += string(", ") + string(dtoa(delz()));
   //pdata += string(", ") + string(dtoa(eta()));
   //pdata += string(" ), ( ") + string(dtoa(alpha()));
   //pdata += string(", ") + string(dtoa(beta()));
   //pdata += string(", ") + string(dtoa(gamma()));
   //pdata += string(" )");
   SinglePar par(pname, 4, pdata, pstate);
   return par;
   }


// sosi - this should be deleted
   space_T::operator ParameterSet( )
//   return pset();

	// Input		SphT  : Spatial tensor
	//  			pset  : Parameter set
	// Output		pset  : Parameter set with
	//			        only spherical tensor

 { ParameterSet pset; pset += *this;  return pset; }		// Add in space tensor values


// sosi - this should be deleted
 void operator+= (ParameterSet& pset, space_T &SphT)

	// Input		SphT  : Spatial tensor
	//  			pset  : Parameter set
	// Output		pset  : Parameter set with
	//			        spherical tensor added

   {
   string pname = string("SphT");
   string pstate = string("Spatial Tensor");
   SinglePar par = SphT.param(pname,pstate);
   pset.push_back(par);
   return;
   } 


void space_T::operator= (const ParameterSet& pset)

	// Input		SphT	 : A spatial tensor (this)
	// 			pset     : A parameter set
	// Output		none	 : Spatial tensor filled with
	//				   data in pset
	// Note				 : Functions which place a spherical
	//				   tensor into a parameter set must
	//				   write the information read here

   {
   string pstate;
   int r;
   double is=0,dz=0,et=0;
   double al=0,be=0,ga=0;
   string pname = "SphT";
   SinglePar par(pname);
   ParameterSet::const_iterator item;
   item = pset.seek(par);
   if(item != pset.end())
     {
     par.parse(pname,r,is,dz,et,al,be,ga,pstate);
     (*this) = A2(is, dz, et, al, be, ga);
// sosi - this assumes rank 2 still, so update someday
     }
   else
     {
     cout << "\nClass space_T: Cannot Read Parameter "
          << pname << " from Parameter Set\n";
     space_T_error(3);
     }
   return;
   }


void space_T::write(const string &filename)

	// Input		SphT	 : Spatial tensor (this)
	//			filename : Output file name
	// Output		none 	 : Spatial tensor is written as a 
	//				   parameter set to file filename

  {
  ofstream ofstr(filename.c_str());	// Open filename for output
  if(!ofstr.good())			// If file bad then exit
    {
    cout << "\nClass space_T: "
         << "Problems with File " << filename;
    space_T_error(5);
    }

  ParameterSet pset;                // Declare a parameter set
  pset += (*this);                     // Add in spin system parameters
  pset.write(ofstr);
  return;
  }


void space_T::read(const string &filename)

	// Input		SphT     : Spatial tensor (this)
	// 			filename : Input filename
	// Output		none	 : Spatial tensor filled with
	//				   parameters read from file

  {
  ParameterSet pset;			// Declare a parameter set
  ifstream inp(filename.c_str());	// Open filename for input
  if(!inp.good())                      // If file bad then exit
    {
    cout << "\nClass space_T: "
         << "Problems with File " << filename;
    space_T_error(4);
    }
  SinglePar par;
  while(par.read(inp))                 // Read all file parameters
    {
    if(!pset.contains(par))            // Add them to the parameter list if
      pset.push_back(par); 		// not already included
    }
   (*this) = pset;                      // Fill up space_T with a parameter
   }


// ____________________________________________________________________________ 
//                       SPACE TENSOR I/O FUNCTIONS
// ____________________________________________________________________________ 


ostream& operator<< (ostream& ostr, const space_T &SphT)

	// Input		ostr : string
	// 			SphT : spatial tensor
	// Return		     : stream, prints spatial tensor
	//			       irreducible spherical components

  {
  int span = 0;
  complex plus, minus, zero;
  if(SphT.vx)
    {
    for(int l=0; l<=SphT.rank; l++)
      {
      if(SphT.vx[l])      
        {
        switch (l)
          {
          case 0:
            zero = (*(SphT.vx[0]))(0);
            ostr << "\nA  = " << zero << "\n 00\n";
            break;
          case 1:
            plus = (*(SphT.vx[1]))(0);
            zero = (*(SphT.vx[1]))(1);
            minus = (*(SphT.vx[1]))(2);
            ostr << "\n(A  , A  , A   ) = ("
                 << plus << ", " << zero << ", " << minus
                 << ")\n  11   10   1-1\n";
            break;
          case 2:
            plus = (*(SphT.vx[2]))(0);
            minus = (*(SphT.vx[2]))(1);
            zero = (*(SphT.vx[2]))(2);
            ostr << "\n                             " << plus
                 << "\n                             " << minus
                 << "\n(A  , A  , A  , A  , A   ) = " << zero;
            minus = (*(SphT.vx[2]))(3);
            plus = (*(SphT.vx[2]))(4);
            ostr << "\n  22   21   20   2-1  2-2    " << minus
                 << "\n                             " << plus << "\n"; 
            break;
          case 3:
            plus = (*(SphT.vx[3]))(0);
            zero = (*(SphT.vx[3]))(1);
            minus = (*(SphT.vx[3]))(2);
            ostr << "\n                                         " << plus
                 << "\n                                         " << zero
                 << "\n                                         " << minus;
            plus = (*(SphT.vx[3]))(3);
            zero = (*(SphT.vx[3]))(4);
            minus = (*(SphT.vx[3]))(5);
            ostr << "\n(A  , A  , A  , A  , A   , A   , A   ) = " << plus 
                 << "\n  33   32   31   30   3-1   3-2   3-3    " << zero
                 << "\n                                         " << minus;
            plus = (*(SphT.vx[3]))(6);
            ostr << "\n                                         " << plus << "\n"; 
            break;
          default:
            span = 2*l+1;
            for (int m=0; m<span; m++)
              {
              zero = (*(SphT.vx[l]))(m);
              ostr << "\n\tA    = " << zero
	           << "\n\t " << l << "," << l-m << "\n";
              }
            break;
          }
        }
      }
    }
  else ostr << "\n\tSpatial Tensor is Currently Empty\n";
  return ostr; 
  }


void Cartesian(const space_T &SphT)

	// Input	 	SphT : spatial tensor
	// Return		     : Sends spatial tensor SphT
	//			       Cartesian components to
	//			       standard output

  {
  double x,y,z;
  double Axx, Axy, Axz;
  double Ayx, Ayy, Ayz;
  double Azx, Azy, Azz;
  complex c; 
  if(SphT.vx)
    switch(SphT.rank)
      {
      case 0:
        if (SphT.vx[0])
          cout << "\nA = " << SphT.component(0, 0) << "\n";
        else
          cout << "\nSpatial Tensor Currently NULL\n";
        break;
      case 1:
        if (SphT.vx[1] || SphT.vx[0])
          {
          c = sqrt(0.5)*(-SphT.component(1, 1) + SphT.component(1,-1));
          c += SphT.component(0, 0);
          x = Re(c);
          if (x < 1.e-10) x = 0.0;

          c = sqrt(0.5)*(SphT.component(1, 1) + SphT.component(1, -1));
          c = complex(0.0,1.0)*c + SphT.component(0, 0);
          y = Re(c);
          if (y < 1.e-10) y = 0.0;

          c = SphT.component(1, 0) + SphT.component(0, 0);
          z = Re(c);
          if (z < 1.e-10) z = 0.0;

          cout << "\n(A , A , A ) = (" << x << ", " << y << ", " << z
               << ")\n  x   y   z\n";
          }
        else
          cout << "\nSpatial Tensor Currently NULL\n";
        break;
      case 2:
        if (SphT.vx[2] || SphT.vx[1] || SphT.vx[0])
          {
          c = 0.5 * (SphT.component(2, 2) + SphT.component(2, -2));
          c = c - SphT.component(2, 0)/sqrt(6.0) - SphT.component(0, 0)/sqrt(3.0);
          Axx = Re(c);
          c = -0.5 * (SphT.component(2, 2) - SphT.component(2, -2));
          c = c + SphT.component(1, 0)/sqrt(2.0);
          c = c*complex(0,1.0);
          Axy = Re(c);
          c = 0.5 * (SphT.component(1, 1) + SphT.component(1, -1));
          c = c - 0.5 * (SphT.component(2, 1) - SphT.component(2, -1));
          Axz = Re(c);
          c = -0.5 * (SphT.component(2, 2) - SphT.component(2, -2));
          c = c - SphT.component(1, 0)/sqrt(2.0);
          c = c*complex(0,1.0);
          Ayx = Re(c);
          c = -0.5 * (SphT.component(2, 2) + SphT.component(2, -2));
          c = c - SphT.component(2, 0)/sqrt(6.0) - SphT.component(0, 0)/sqrt(3.0);
          Ayy = Re(c);
          c = 0.5 * (SphT.component(2, 1) + SphT.component(2, -1));
          c = c - 0.5 * (SphT.component(1, 1) - SphT.component(1, -1));
          c = c * complex(0,1.0);
          Ayz = Re(c);
          c = -0.5 * (SphT.component(1, 1) + SphT.component(1, -1));
          c = c - 0.5 * (SphT.component(2, 1) - SphT.component(2, -1));
          Azx = Re(c);
          c = 0.5 * (SphT.component(2, 1) + SphT.component(2, -1));
          c = c + 0.5 * (SphT.component(1, 1) - SphT.component(1, -1));
          c = c * complex(0,1.0);
          Azy = Re(c);
          c = sqrt(2.0/3.0)*SphT.component(2,0) - sqrt(1.0/3.0)* SphT.component(0,0);
          Azz = Re(c);
          cout << "\n[A  , A  , A  ]"
               << "\n[ xx   xy   xz]   ["
	       << Axx << ", " << Axy << ", " << Axz << "]"
               << "\n[A  , A  , A  ] = ["
	       << Ayx << ", " << Ayy << ", " << Ayz << "]"
               << "\n[ yx   yy   yz]   ["
	       << Azx << ", " << Azy << ", " << Azz << "]"
               << "\n[A  , A  , A  ]"
               << "\n[ zx   zy   zz]\n";
          }
        else
          cout << "\nSpatial Tensor Currently NULL\n";
        break;
      default:
        cout << "\nCurrently Unable to Output Cartesian Components"
             << "\nfor Spatial Tensor Rank " << SphT.rank;
        break;
      }
  else
    cout << "\n\tSpatial Tensor is Currently Empty\n";
  return; 
  }

#endif 						// Space_T.cc
