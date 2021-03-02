/* SpinT.cc *********************************************-*-c++-*-
**								**
**								**
** 	                       G A M M A			**
**								**
**	Spin Tensors			   Implementation   	**
**						 		**
**	Copyright (c) 1990, 1991, 1992			 	**
**	Scott Smith					        **
**	Eidgenoessische Technische Hochschule	 		**
**	Labor fur physikalische Chemie		 		**
**	8092 Zurich / Switzerland			 	**
**						 		**
**      $Header: $
**								**
*****************************************************************/

/*****************************************************************
**								**
**  Description:						**
**								**
**  This file contains the implementation of Spin Tensors.	**
**								**
*****************************************************************/

#ifndef _spin_T_cc_			// Is file already included?
#define _spin_T_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

#include <Level1/SpinT.h>		// Include our own interface
#include <HSLib/SpinOpCmp.h>		// Include composite spin operators
#include <Level1/Wigner.h>		// Include Wigner rotations
#include <stdlib.h>

typedef spin_op **pp_spin_op;	        // pointer to a pointer to spin_op

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
//                         SPIN TENSOR ERROR HANDLING
// ____________________________________________________________________________


	// Input		i    : Error Flag
	// Output		none : Error Message Output

void spin_T_error(const int i)
  {
  std::cout << "\nSpin_T: ";
  switch (i)
    {
    case 0:					// (0)
      std::cout << "Program Aborting.\n";
      break;
    case 1:					// (1)
      std::cout <<           "         (0)"
	   << "\nSpin_T: Unknown T"
           << "\nSpin_T:          ";
      break;
    case 2:					// (2)
      std::cout <<           "         (1)"
	   << "\nSpin_T: Unknown T"
           << "\nSpin_T:          ";
      break;
    case 3:					// (3)
      std::cout <<           "         (2)"
           << "\nSpin_T: Unknown T"
           << "\nSpin_T:          ";
      break;
    case 10:
      std::cout << "Unable to Determine Spherical Tensor Component.\n";
      break;
    default:
      std::cout << "Unknown error.\n";
      break;
    }
  return;
  }



	// Input		none :
	// Output		none : Stops Execution & Error Message Output

void spin_T_fatality(int error)
  {
  spin_T_error (error);
  if(error) spin_T_error(0);
  exit(-1);
  return;
  }


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                  SPIN TENSOR CONSTRUCTORS, DESTRUCTOR
// ____________________________________________________________________________


	// Input		none : 
	// Return		SphT : null spin tensor

spin_T::spin_T()
  {
  sys = NULL;			// Mark spin system
  rank = 0;			// Tensor rank
  pr = NULL;			// Currently empty tensor
  }




	// Input		sys1 : spin system
	// Return		SphT : spin tensor associated with
        //	                       the input spin system

spin_T::spin_T(const spin_sys& sys1)
  {
// sosi - this is a totally lame kludge
spin_sys sys2(sys1);
sys = &sys2;			// Copy spin system
//  sys = sys1;			// Copy spin system
  rank = 0;			// Tensor rank
  pr = NULL;			// Currently empty tensor
  }



	// Input		SphT : spin tensor
	// Return		SphT1: spin tensor duplicate of input
        //	                       spin tensor

spin_T::spin_T(const spin_T& SphT)
  {
  sys = SphT.sys;			// copy spin system pointer
  rank = SphT.rank;			// copy Tensor rank
  pr = NULL;
  int span = 0;
  if (SphT.pr)				// copy all ranks (l+1) included
    {
    //pr = new spin_op **[rank+1];	// work around bug in g++
    pr = new pp_spin_op [rank+1];		// work around bug in g++
    for (int l=0; l<=rank; l++)
      {
      if (SphT.pr[l])			// copy all rank l (2l+1) components
        {
        span = 2*l+1;
        pr[l] = new spin_op* [span];
        for (int m=0; m<span; m++)
          if (SphT.pr[l][m])		// copy only if SOp exists
            pr[l][m] = new spin_op(*(SphT.pr[l][m]));
          else
            pr[l][m] = NULL;
        }
      else pr[l] = NULL;
      }
    }
  }



	// Input		SphT : spin tensor
	//                      l    : rank
	// Return		SphT1: spin tensor duplicate of input
        //	                       spin tensor

spin_T::spin_T(const spin_T& SphT, int l)
  {
  sys = SphT.sys;			// copy spin system pointer
  rank = 0;				// Temporary Tensor rank
  pr = NULL;				// currently empty tensor
  if (SphT.pr)				// copy if SphT not NULL
    if ((SphT.rank <= l) && SphT.pr[l])	// copy if SphT  not NULL
      {					//	       l	
      int span = 2*l+1;
      rank = l;				// set Tensor rank
      //pr = new spin_op **[l+1];	// work around bug in g++
      pr = new pp_spin_op[l+1];		// work around bug in g++
      for(l=0; l< rank; l++)
        pr[l] = NULL;
      l = rank;
      pr[l] = new spin_op* [span];
        for (int m=0; m<span; m++)
          if (SphT.pr[l][m])		// copy only if SOp exists
            pr[l][m] = new spin_op(*(SphT.pr[l][m]));
          else
            pr[l][m] = NULL;
      }
  }



	// Input		SphT : Sperical Tensor
	// Return		     : none, deletes current Tensor

spin_T::~spin_T()
  {
  int span = 0;
  if (pr)				// delete if Tensor not NULL
    {
    for (int l=0; l<=rank; l++)
      if (pr[l])			// delete entire rank l if non NULL
        {
        span = 2*l+1;
        for (int m=0; m<span; m++)
          if (pr[l][m])			// delete all rank l non-NULL SOps
            delete pr[l][m];
        delete [ ] pr[l];	
        }
    delete [] pr;	
    pr = NULL;
    }
  }


// ____________________________________________________________________________
//                      SPIN TENSOR UNARY OPERATIONS
// ____________________________________________________________________________


	// Input		SphT : spin tensor
	// Return		SphT1: spin tensor equivalent to the
        //	                       input spin tensor
	// Note		             : SphT and SphT1 must be associated
	//			       with the same spin system

spin_T & spin_T::operator = (const spin_T& SphT)
  {
  if(this == &SphT) return (*this); 		// Exit if already equal
  int span = 0;
  int l = 0;
  if(SphT.pr)				// Copy if non-NULL SphT
    {
    if (pr)				// delete current spin tensor
      {
      for (l=0; l<=rank; l++)
        if (pr[l])			// delete entire rank l if non NULL
          {
          span = 2*l+1;
          for (int m=0; m<span; m++)
            if (pr[l][m])		// delete all rank l non-NULL SOps
              delete pr[l][m];	
          delete [] pr[l];	
          }
      delete [] pr;	
      }
    //pr = new spin_op **[rank+1];	// work around bug in g++
    pr = new pp_spin_op [SphT.rank+1];	// work around bug in g++
    for(l=0; l<=SphT.rank; l++)
      {
      if(SphT.pr[l])			// copy all rank l (2l+1) components
        {
        span = 2*l+1;
        pr[l] = new spin_op* [span];
        for (int m=0; m<span; m++)
          {
          if (SphT.pr[l][m])		// copy only if SOp exists
            pr[l][m] = new spin_op(*(SphT.pr[l][m]));
          else pr[l][m] = NULL;
          }
        }
      else pr[l] = NULL;
      }
    }
  else
    {
    if (pr)
      {
      for(l=0; l<=rank; l++)
        if (pr[l])			// delete entire rank l if non NULL
          {
          span = 2*l+1;
          for (int m=0; m<span; m++)
            if (pr[l][m])		// delete all rank l non-NULL SOps
              delete pr[l][m];	
          delete [] pr[l];	
          }
      delete [] pr;	
      pr = NULL;
      }
    }
  sys = SphT.sys;			// Copy spin system pointer
  rank = SphT.rank;			// Copy Tensor rank
  
  return (*this);
  }


// ____________________________________________________________________________
//               GENERAL RANK 1 SPHERICAL SPIN TENSOR FUNCTIONS
// ____________________________________________________________________________

// ***********************  Full Rank One Treatment ***************************

	// Input		sys   : spin system
	// 			spin  : spin index
	// Output		SphT  : rank 1 Spin Tensor for spins specified
 
spin_T T1(const spin_sys& sys, int spin)
  {
  spin_T SphT(sys);	
  SphT.rank = 1;
  spin_op SOp;
  //SphT.pr = new spin_op **[rank+1];	// work around bug in g++
  typedef spin_op **tmp;	        // work around bug in g++
  SphT.pr = new tmp [2];		// work around bug in g++
  for (int l=0; l<2; l++)
    {
    int span = 2*l+1;
    SphT.pr[l] = new spin_op *[span];
// ???? but l=0 case is always null for rank 1 !!!!
    for (int m=0; m<span; m++)
      {
      SOp = T1(sys, spin, l, l-m);
      SphT.pr[l][m] = new spin_op(SOp);
      }
    }
  return SphT;
  }

// ********************** Irreducible Rank One Tensor *************************


	// Input		sys   : Spin system
	// 			spin  : Spin index	    1
	// Output		SOp   : Rank 1 spin tensor T (i)

spin_T T11(const spin_sys& sys, int spin)
  {
  spin_T SphT;
  SphT.rank = 1;
  spin_op SOp;
  //SphT.pr = new spin_op **[rank+1];	// work around bug in g++
  typedef spin_op **tmp;	        // work around bug in g++
  SphT.pr = new tmp [2];		// work around bug in g++
  SphT.pr[0] = NULL;
  SphT.pr[1] = new spin_op *[3];
    for (int m=0; m<3; m++)
      {
      SOp = T1(sys, spin, 1, 1-m);
      SphT.pr[1][m] = new spin_op(SOp);
      }
  return SphT;				// Return SphT
  }


// *************  Specific Irreducible Rank One Tensor Components *************

	// Input		sys   : spin system
	// 			spin  : spin index
	//			l     : component index
	//			m     : Component index		      1
	// Output		SOp   : Rank 1 spin tensor component T  (i)
	//							      lm
 
spin_op T1(const spin_sys& sys, int spin, int l, int m)
  {
  spin_op SOp;
  switch (l)
    {
    case 0:			//  (1)
      SOp = T10(sys, spin, m);	// T   (i) 
      break;			//  0,m

    case 1:			//  (1)
      SOp = T11(sys, spin, m);  // T   (i)
      break;			//  1,m
    default:
      spin_T_error(2);
      std::cout << " " << l << "0," << m << "\n";
      spin_T_fatality(10);
      break;
    }
  return SOp;
  }



	// Input		sys   : spin system
	// 			spin  : spin index
	//			m     : Component index		      1
	// Output		SOp   : Rank 1 spin tensor component T  (i)
	//							      0m

spin_op T10(const spin_sys& sys, int spin, int m)
  {
  spin_op SOp;
  switch (m)
    {
    case 0:			//  (1)
      SOp = 0.0*Ie(sys,spin);	// T   (i) = 0
      break;			//  00
    default:
      spin_T_error(2);
      std::cout << " 0," << m << "\n";
      break;
    }
  return SOp;
  }

	// Input		Ie    : Identity matrix
	//			m     : Component index		      1
	// Output		SOp   : Rank 1 spin tensor component T
	//							      0m

spin_op T10(spin_op& Ie, int m)
  {
  spin_op SOp;
  switch (m)
    {
    case 0:			//  (1)
      SOp = 0.0*Ie;		// T    = 0
      break;			//  00
    default:
      spin_T_error(2);
      std::cout << " 0," << m << "\n";
      break;
    }
  return SOp;			// Return SOp
  }


	// Input		sys   : spin system
	// 			spin  : spin index
	//			m     : Component index		      1
	// Output		SOp   : Rank 1 spin tensor component T  (i)
	//							      1m

spin_op T11(const spin_sys& sys, int spin, int m)
  {
  spin_op SOp;
  switch (m)
    {
    case -1:				//  (1)	         [1]
      SOp = sqrt(0.5)*Im(sys,spin);	// T   (i) = sqrt|-| I (i)	
      break;				//  1-1          [2]  -

    case 0:				//  (1)
      SOp = Iz(sys,spin);		// T   (i) = I (i)
      break;				//  10        z

    case 1:				//  (1)           [1]
      SOp = -sqrt(0.5)*Ip(sys,spin);	// T   (i) = -sqrt|-| I (i)
      break;				//  11            [2]  +

    default:
      spin_T_error(2);
      std::cout << " 1," << m << "\n";
      break;
    }
  return SOp;				// Return SOp
  }

	// Input		Im    : Spin operator I-
	// 			Iz    : Spin operator Iz
	// 			Ip    : Spin operator I+
	//			m     : Component index		      1
	// Output		SOp   : Rank 1 spin tensor component T
	//							      1m

spin_op T11(spin_op& Im, spin_op& Iz, spin_op& Ip, int m)
  {
  spin_op SOp;
  switch (m)
    {
    case -1:			//  (1)       [1]
      SOp = sqrt(0.5)*Im;	// T    = sqrt|-| I
      break;			//  1-1       [2]  -

    case 0:			//  (1)
      SOp = Iz;			// T    = I
      break;			//  10     z

    case 1:			//  (1)        [1]
      SOp = -sqrt(0.5)*Ip;	// T    = -sqrt|-| I
      break;			//  11         [2]  +

    default:
      spin_T_error(2);
      std::cout << " 1," << m << "\n";
      break;
    }
  return SOp;			// Return SOp
  }


// ____________________________________________________________________________
//                GENERAL RANK 2 SPHERICAL SPIN TENSOR FUNCTIONS
// ____________________________________________________________________________


// ***********************  Full Rank Two Treatment ***************************

	// Input		sys   : Spin system
	// 			spin1 : Spin index
	// 			spin2 : Spin index 		        2
	// Output		SOp   : Irreducible rank 2 spin tensor T (ij)
 
spin_T T2(const spin_sys& sys, int spin1, int spin2)
  {
  spin_T SphT;
  SphT.rank = 2;
  spin_op SOp;
  //SphT.pr = new spin_op **[rank+1];	// work around bug in g++
  typedef spin_op **tmp;	        // work around bug in g++
  SphT.pr = new tmp [3];		// work around bug in g++
  for (int l=0; l<3; l++)
    {
    int span = 2*l+1;
    SphT.pr[l] = new spin_op *[span];
    for (int m=0; m<span; m++)
      {
      SOp = T2(sys, spin1, spin2, l, l-m);
      SphT.pr[l][m] = new spin_op(SOp);
      }
    }
  return SphT;
  }


// ********************* Irreducible Rank Two Spin Tensor *********************


spin_T T22(const spin_sys& sys, int spin1, int spin2)

	// Input		sys   : Spin system
	// 			spin1 : Spin index
	// 			spin2 : Spin index 		        2
	// Output		SOp   : Irreducible rank 2 spin tensor T (ij)
	//							        2

  {
  spin_T SphT;
  SphT.rank = 2;
  spin_op SOp;
  //SphT.pr = new spin_op **[rank+1];	// work around bug in g++
  typedef spin_op **tmp;	        // work around bug in g++
  SphT.pr = new tmp [3];		// work around bug in g++
  for (int l=0; l<2; l++)		// set ranks 0 & 1 to NULL
    {
    SphT.pr[l] = new spin_op* [2*l+1];
    for(int k=0; k<2*l+1; k++)
      SphT.pr[l][k] = NULL;
    }
// sosi 11/11/92 replaced SphT.pr[l] = NULL; with {} above 
  SphT.pr[2] = new spin_op *[5];
    for (int m=0; m<5; m++)
      {
      SOp = T2(sys, spin1, spin2, 2, 2-m);
      SphT.pr[2][m] = new spin_op(SOp);
      }
  return SphT;				// Return SphT
  }



	// Input		sys   : Spin system
	// 			spin1 : Spin index
	// 			spin2 : Spin index 		        2
	// Output		SOp   : Irreducible rank 2 spin tensor T (ij)
	//							        2
	// Note			      : This routine will leave all but
	//				the m=0 component empty if spin1 and
	//				spin 2 are not of the same spin type

spin_T T22wh(const spin_sys& sys, int spin1, int spin2)
  {
  spin_T SphT;
  SphT.rank = 2;			// Set the tensor rank
  spin_op SOp;			// Begin with an empty spin op
  //SphT.pr = new spin_op **[rank+1];	// work around bug in g++
  typedef spin_op **tmp;	        // work around bug in g++
  SphT.pr = new tmp [3];		// work around bug in g++
  for(int l=0; l<=2; l++)		// Initialize lower ranks to NULL
    {
    SphT.pr[l] = new spin_op* [2*l+1];
    for(int k=0; k<2*l+1; k++)
      SphT.pr[l][k] = NULL;
    }
  if(sys.symbol(spin1)			// For heteronuclear case, just 
                 != sys.symbol(spin2))	// use the 2,0 zz component
    {    				//  (2)		  [1] [
    SOp = (2.0/sqrt(6.0))* 		// T   (ij) = sqrt|-| |2 I (i)I (j)
           Iz(sys,spin1)*Iz(sys,spin2);	//  20		  [6] [   z    z
    SphT.pr[2][2] = new spin_op(SOp);
//    SOp = T2(sys, spin1, spin2, 2, 0);
    }
  else					// For a homonuclear case, must fill
    { 					// up all five l=2 components
    for(int m=0; m<5; m++)
      {
      SOp = T2(sys, spin1, spin2, 2, 2-m);
      SphT.pr[2][m] = new spin_op(SOp);
      }
    }
  return SphT;				// Return SphT
  }



	// Input		sys   : Spin system
	// 			Im1   : I- spin 1
	// 			Iz1   : Iz spin 1
	// 			Ip1   : I+ spin 1
	// 			Im2   : I- spin 2
	// 			Iz2   : Iz spin 2
	// 			Ip2   : I+ spin 2			2
	// Output		SOp   : Irreducible rank 2 spin tensor T (12)
	//							        2

spin_T T22(const spin_sys& sys, spin_op& Im1, spin_op& Iz1, spin_op& Ip1,
			       spin_op& Im2, spin_op& Iz2, spin_op& Ip2)
  {
  spin_T SphT;	
  SphT.rank = 2;
  spin_op SOp;
  //SphT.pr = new spin_op **[rank+1];	// work around bug in g++
  typedef spin_op **tmp;	        // work around bug in g++
  SphT.pr = new tmp [3];		// work around bug in g++
  for (int l=0; l<2; l++)
    SphT.pr[l] = NULL;
  SphT.pr[2] = new spin_op *[5];
    for (int m=0; m<5; m++)
      {
      SOp = T2(Im1, Iz1, Ip1, Im2, Iz2, Ip2, 2, 2-m);
      SphT.pr[2][m] = new spin_op(SOp);
      }
  return SphT;				// Return SphT
  }


// **********  Specific Irreducible Rank Two Tensor Components **********



	// Input		sys   : spin system
	// 			spin1 : spin index
	// 			spin2 : spin index
	//			l     : component index
	//			m     : Component index               2
	// Output		SOp   : Rank 2 Spin Tensor component T  (ij)
	//							      lm

spin_op T2(const spin_sys& sys, int spin1, int spin2, int l, int m)
  {
  spin_op SOp;
  switch(l)
    {
    case 0:				//	  (2)
      SOp = T20(sys, spin1, spin2, m);	// S0p = T   (i,j)
      break;				//        0,m

    case 1:				//        (2)
      SOp = T21(sys, spin1, spin2, m);	// SOp = T   (i,j)
      break;				//        1,m

    case 2:				//        (2)
      SOp = T22(sys, spin1, spin2, m);	// SOp = T   (i,j)
      break;				//        2,m

    default:
      spin_T_error(3);
      std::cout << " " << l << "," << m << "\n";
      spin_T_fatality(10);
      break;
    }
  return SOp;
  }



	// Input		sys   : spin system
	// 			Im1   : I- spin 1
	// 			Iz1   : Iz spin 1
	// 			Ip1   : I+ spin 1
	// 			Im2   : I- spin 2
	// 			Iz2   : Iz spin 2
	// 			Ip2   : I+ spin 2
	//			l     : component index
	//			m     : Component index               2
	// Output		SOp   : Rank 2 Spin Tensor component T  (12)
	//							      lm

spin_op T2(spin_op& Im1, spin_op& Iz1, spin_op& Ip1,
		spin_op& Im2, spin_op& Iz2, spin_op& Ip2, int l, int m)
  {
  spin_op SOp;
  switch (l)
    {
    case 0:
      SOp = T20(Im1, Iz1, Ip1, Im2, Iz2, Ip2, m);
      break;
    case 1:
      SOp = T21(Im1, Iz1, Ip1, Im2, Iz2, Ip2, m);
      break;
    case 2:
      SOp = T22(Im1, Iz1, Ip1, Im2, Iz2, Ip2, m);
      break;
    default:
      spin_T_error(3);
      std::cout << " " << l << "," << m << "\n";
      spin_T_fatality(10);
      break;
    }
  return SOp;
  }



	// Input		sys   : spin system
	// 			spin1 : spin index
	// 			spin2 : spin index
	//			m     : Component index               2
	// Output		SOp   : Rank 2 Spin Tensor component T  (ij)
	//							      0m
 
spin_op T20(const spin_sys& sys, int spin1, int spin2, int m)
  {
  spin_op SOp;
  switch (m)
    {
    case 0:				//  (2)		    [1] [
      SOp = -1.0/sqrt(3.0)* 		// T   (ij) = - sqrt|-| | I (i)I (j)
            (Iz(sys,spin1)*Iz(sys,spin2)//  00		    [3] [  z    z
            +0.5*(Ip(sys,spin1)*	//
            Im(sys,spin2)+Im(sys,spin1) //	  [1] {			      } ]
	    *Ip(sys,spin2)));		//	+ |-| {I (i)I (j) + I (i)I (j)} |    
      break;				//   	  [2] { +    -       -    +   } ]

    default:
      spin_T_error(3);
      std::cout << " 0," << m << "\n";
      break;
    }
  return SOp;
  }



	// Input		Im1   : I- spin 1
	// 			Iz1   : Iz spin 1
	// 			Ip1   : I+ spin 1
	// 			Im2   : I- spin 2
	// 			Iz2   : Iz spin 2
	// 			Ip2   : I+ spin 2
	//			m     : Component index               2
	// Output		SOp   : Rank 2 Spin Tensor component T  (12)
	//							      0m
 
spin_op T20(spin_op& Im1, spin_op& Iz1, spin_op& Ip1,
		spin_op& Im2, spin_op& Iz2, spin_op& Ip2, int m)
  {
  spin_op SOp;
  switch(m)
    {
    case 0:				//  (2)		    [1] [
      SOp = -1.0/sqrt(3.0)*		// T   (ij) = - sqrt|-| | I (1)I (2)
            (Iz1*Iz2			//  00		    [3] [  z    z
              +0.5*(Ip1*Im2+Im1*Ip2)); 	//	  [1] {			      } ]
      break;				//	+ |-| {I (1)I (2) + I (1)I (2)} |    
      					//   	  [2] { +    -       -    +   } ]

    default:
      spin_T_error(3);
      std::cout << " 0," << m << "\n";
      break;
    }
  return SOp;
  }

	// Input		sys   : spin system
	// 			spin1 : spin index
	// 			spin2 : spin index
	//			m     : Component index               2
	// Output		SOp   : Rank 2 Spin Tensor component T  (ij)
	//							      1m
 
spin_op T21(const spin_sys& sys, int spin1, int spin2, int m)
  {
  spin_op SOp;
  switch (m)
    {
    case -1:				//  (2)		1 [
      SOp = -0.5*	 		// T   (ij) = - - |I (i)I (j)
            (Im(sys,spin1)*Iz(sys,spin2)//  1-1		2 [ -    z
          -Iz(sys,spin1)*Im(sys,spin2));//			              ]
      break;				//			  - I (i)I (j)|
					//			     z    -   ]

    case 0:				//  (2)	          -1     [
      SOp = -0.5/sqrt(2.0)* 		// T   (ij) =  --------- |I (i)I (j)
            (Ip(sys,spin1)*Im(sys,spin2)//  10	       2*sqrt(2) [ +    -
          -Im(sys,spin1)*Ip(sys,spin2));//			              ]
      break;				//			  - I (i)I (j)|
					//			     -    +   ]

    case 1:				//  (2)		1 [
      SOp = -0.5*	 		// T   (ij) = - - |I (i)I (j)
            (Ip(sys,spin1)*Iz(sys,spin2)//  11		2 [ +    z
          -Iz(sys,spin1)*Ip(sys,spin2));//			              ]
      break;				//			  - I (i)I (j)|
					//			     z    +   ]
    default:
      spin_T_error(3);
      std::cout << " 1," << m << "\n";
      break;
    }
  return SOp;				// Return SOp
  }



	// Input		Im1   : I- spin 1
	// 			Iz1   : Iz spin 1
	// 			Ip1   : I+ spin 1
	// 			Im2   : I- spin 2
	// 			Iz2   : Iz spin 2
	// 			Ip2   : I+ spin 2
	//			m     : Component index               2
	// Output		SOp   : Rank 2 Spin Tensor component T  (12)
	//							      1m
 
spin_op T21(spin_op& Im1, spin_op& Iz1, spin_op& Ip1,
		spin_op& Im2, spin_op& Iz2, spin_op& Ip2, int m)
  {
  spin_op SOp;
  switch (m)
    {
    case -1:				//  (2)		1 [                       ]
      SOp = -0.5 * (Im1*Iz2-Iz1*Im2);	// T   (12) = - - |I (1)I (2) - I (1)I (2)|
      break; 				//  1-1		2 [ -    z       z    -   ]

    case 0:				//  (2)	      [   -1    ] [
      SOp = 0.5/sqrt(2.0)* 		// T   (12) = |---------| |I (1)I (2)
                 (Ip1*Im2-Im1*Ip2);	//  10	      [2*sqrt(2)] [ +    -
      break;				//			                  ]
      					//			      - I (1)I (2)|
					//			         -    +   ]

    case 1:				//  (2)		1 [                       ]
      SOp = -0.5 * (Ip1*Iz2-Iz1*Ip2);	// T   (12) = - - |I (1)I (2) - I (1)I (2)|
      break;				//  11		2 [ +    z       z    +   ]

    default:
      spin_T_error(3);
      std::cout << " 1," << m << "\n";
      break;
    }
  return SOp;				// Return SOp
  }



	// Input		sys   : spin system
	// 			spin1 : spin index
	// 			spin2 : spin index
	//			m     : Component index               2
	// Output		SOp   : Rank 2 Spin Tensor component T  (ij)
	//							      2m

spin_op T22(const spin_sys& sys, int spin1, int spin2, int m)
  {
  spin_op SOp;
  switch (m)
    {
    case -2:				//  (2)	      [1] 
      SOp = 0.5* 			// T   (ij) = |-| I (i) I (j)
            Im(sys,spin1)*Im(sys,spin2);//  2-2	      [2]  -     -
      break;

    case -1:				//  (2)	      [1] [
      SOp = 0.5* 			// T   (ij) = |-| |I (i)I (j)
            (Im(sys,spin1)*Iz(sys,spin2)//  2-1	      [2] [ -    z
          +Iz(sys,spin1)*Im(sys,spin2));//			              ]
      break;				//			  + I (i)I (j)|
					//			     z    -   ]

    case 0:				//  (2)		  [1] [
      SOp = (1.0/sqrt(6.0))* 		// T   (ij) = sqrt|-| |2 I (i)I (j)
     (2.0*Iz(sys,spin1)*Iz(sys,spin2)	//  20		  [6] [   z    z
       - (Ix(sys,spin1)*Ix(sys,spin2))	//		  			  ]
       - (Iy(sys,spin1)*Iy(sys,spin2)));//		- I (i)I (j) - I (i)I (j) |    
      break;				//   		   x    x       y    y    ]

    case 1:				//  (2)		[1] [
      SOp = -0.5* 			// T   (ij) = - |-| |I (i)I (j)
            (Ip(sys,spin1)*Iz(sys,spin2)//  21		[2] [ +    z
          +Iz(sys,spin1)*Ip(sys,spin2));//			              ]
      break;				//			  + I (i)I (j)|
					//			     z    +   ]

    case 2:				//  (2)	      [1] 
      SOp = 0.5* 			// T   (ij) = |-| I (i) I (j)
            Ip(sys,spin1)*Ip(sys,spin2);//  22	      [2]  +     +
      break;

    default:
      spin_T_error(3);
      std::cout << " 2," << m << "\n";
      break;
    }
  return SOp;				// Return spin operator
  }



	// Input		Im1   : I- spin 1
	// 			Iz1   : Iz spin 1
	// 			Ip1   : I+ spin 1
	// 			Im2   : I- spin 2
	// 			Iz2   : Iz spin 2
	// 			Ip2   : I+ spin 2
	//			m     : Component index               2
	// Output		SOp   : Rank 2 Spin Tensor component T  (12)
	//							      2m

spin_op T22(spin_op& Im1, spin_op& Iz1, spin_op& Ip1,
		spin_op& Im2, spin_op& Iz2, spin_op& Ip2, int m)
  {
  spin_op SOp;
  switch (m)
    {
    case -2:				//  (2)	      [1] 
      SOp = 0.5*Im1*Im2;		// T   (12) = |-| I (1) I (2)
      break;				//  2-2	      [2]  -     -


    case -1:				//  (2)	      [1] [
      SOp = 0.5*(Im1*Iz2+Iz1*Im2); 	// T   (12) = |-| |I (1)I (2)
      break;   				//  2-1	      [2] [ -    z
          				//			              ]
 					//			  + I (1)I (2)|
					//			     z    -   ]

    case 0:				//  (2)		  [1] [
      SOp = 1.0/sqrt(6.0)* 		// T   (12) = sqrt|-| |2 I (1)I (2)
        (2.0*Iz1*Iz2			//  20		  [6] [   z    z
        - 0.5*(Ip1*Im2+Im1*Ip2));   	//	      1	  			  ]
      break;				//	    - - {I (1)I (2) + I (1)I (2)} |    
      					//   	      2   +    -       -    +     ]

    case 1:				//  (2)		[1] [
      SOp = -0.5*(Ip1*Iz2+Iz1*Ip2);	// T   (12) = - |-| |I (1)I (2)
      break;				//  21		[2] [ +    z
          				//			              ]
      					//			  + I (1)I (2)|
					//			     z    +   ]

    case 2:				//  (2)	      [1] 
      SOp = 0.5*Ip1*Ip2;		// T   (12) = |-| I (1) I (2)
      break;				//  22	      [2]  +     +

    default:
      spin_T_error(3);
      std::cout << " 2," << m << "\n";
      break;
    }
  return SOp;				// Return SOp
  }


// ______________________________________________________________________
//           GENERAL RANK 2 SPHERICAL SPIN-SPACE TENSOR FUNCTIONS
// ______________________________________________________________________


// ********************  Full Rank Two Treatment ************************

spin_T T2(const spin_sys& sys, int spin, const coord& vect)
  {return T2SS(sys,spin,vect);}
spin_T T2(const spin_sys& sys, spin_op& Im, spin_op& Iz, spin_op& Ip, const coord& vect)
  {return T2SS(sys,Im,Iz,Ip,vect);}


	// Input		sys   : Spin system
	// 			spin  : Spin index
	// 			vect  : Cartesian vector
	//			rev   : Reverse tensor order
	// Output		SphT  : Rank 2 "Spin" Tensor for spin specified
 
spin_T T2SS(const spin_sys& sys, int spin, const coord& vect, int rev)
  {
  spin_T SphT(sys);
  SphT.rank = 2;
  spin_op SOp;
  //SphT.pr = new spin_op **[rank+1];	// work around bug in g++
  typedef spin_op **tmp;	        // work around bug in g++
  SphT.pr = new tmp [3];		// work around bug in g++
  for (int l=0; l<3; l++)
    {
    int span = 2*l+1;
    SphT.pr[l] = new spin_op *[span];
    for (int m=0; m<span; m++)
      {
      SOp = T2SS(sys, spin, vect, l, l-m, rev);
      SphT.pr[l][m] = new spin_op(SOp);
      }
    }
  return SphT;				// Return SphT
  }



	// Input		sys   : Spin system
	// 			Im    : Spin operator I-
	// 			Iz    : spin operator Iz
	// 			Ip    : Spin operator I+
	// 			vect  : Cartesian vector
	//			rev   : Reverse tensor order
	// Output		SphT  : rank 2 "Spin" Tensor for spin specified
 
spin_T T2SS(const spin_sys& sys, spin_op& Im, spin_op& Iz,
				 spin_op& Ip, const coord& vect, int rev)
  {
  spin_T SphT(sys);
  SphT.rank = 2;
  spin_op SOp;
  //SphT.pr = new spin_op **[rank+1];	// work around bug in g++
  typedef spin_op **tmp;	        // work around bug in g++
  SphT.pr = new tmp [3];		// work around bug in g++
  for (int l=0; l<3; l++)
    {
    int span = 2*l+1;
    SphT.pr[l] = new spin_op *[span];
    for (int m=0; m<span; m++)
      {
      SOp = T2SS(Im, Iz, Ip, vect, l, l-m, rev);
      SphT.pr[l][m] = new spin_op(SOp);
      }
    }
  return SphT;				// Return SphT
  }


// ****************** Irreducible Rank Two Spin Tensor ******************


spin_T T22(const spin_sys& sys, int spin, const coord& vect)
  {return T22SSirr(sys,spin,vect);}
spin_T T22(const spin_sys& sys, spin_op& Im, spin_op& Iz, spin_op& Ip, const coord& vect)
  {return T22SSirr(sys,Im,Iz,Ip,vect);}


	// Input		sys   : spin system
	// 			spin  : spin index
	// 			vect  : Cartesian row vector
	//			rev   : Reverse tensor order
	// Output		SphT  : irreducible rank 2 "Spin" Tensor
	//				for spin specified

spin_T T22SSirr(const spin_sys& sys, int spin, const coord& vect, int rev)
  {
  spin_T SphT(sys);
  SphT.rank = 2;
  spin_op SOp;
  //SphT.pr = new spin_op **[rank+1];	// work around bug in g++
  typedef spin_op **tmp;	        // work around bug in g++
  SphT.pr = new tmp [3];		// work around bug in g++
  for (int l=0; l<2; l++)
    {
    SphT.pr[l] = new spin_op* [2*l+1];
    for(int k=0; k<2*l+1; k++)
      SphT.pr[l][k] = NULL;
    }
// sosi 11/11/92 replaced SphT.pr[l] = NULL; with {} above 
  SphT.pr[2] = new spin_op *[5];
    for (int m=0; m<5; m++)
      {
      SOp = T2SS(sys, spin, vect, 2, 2-m, rev);
      SphT.pr[2][m] = new spin_op(SOp);
      }
  return SphT;				// Return SphT
  }



	// Input		sys   : Spin system
	// 			Im    : Spin operator I-
	// 			Iz    : Spin operator Iz
	// 			Ip    : Spin operator I+
	// 			vect  : Cartesian vector
	//			rev   : Reverse tensor order
	// Output		SphT  : Irreducible rank 2 "Spin" Tensor
	//				for spin specified

spin_T T22SSirr(const spin_sys& sys, spin_op& Im, spin_op& Iz,
				 spin_op& Ip, const coord& vect, int rev)
  {
  spin_T SphT(sys);
  SphT.rank = 2;
  spin_op SOp;
  //SphT.pr = new spin_op **[rank+1];	// work around bug in g++
  typedef spin_op **tmp;	        // work around bug in g++
  SphT.pr = new tmp [3];		// work around bug in g++
  for (int l=0; l<2; l++)
    SphT.pr[l] = NULL;
  SphT.pr[2] = new spin_op *[5];
    for (int m=0; m<5; m++)
      {
      SOp = T2SS(Im, Iz, Ip, vect, 2, 2-m, rev);
      SphT.pr[2][m] = new spin_op(SOp);
      }
  return SphT;				// Return SphT
  }


// **********  Specific Irreducible Rank Two Tensor Components **********


spin_op T2(const spin_sys& sys, int spin, const coord& vect, int l, int m)
  {return T2SS(sys,spin,vect,l,m);}
spin_op T2(spin_op& Im, spin_op& Iz, spin_op& Ip, const coord& vect, int l, int m)
  {return T2SS(Im,Iz,Ip,vect,l,m);}



	// Input		sys   : Spin system
	// 			spin  : Spin index
	// 			vect  : Cartesian vector
	//			l     : Rank index
	//			m     : Component index
	//			rev   : Reverse tensor order
	// Output		SOp   : Spin Operator for SphT Component T
	//								  lm 

spin_op T2SS(const spin_sys& sys, int spin, const coord& vect, int l, int m, int rev)
  {
  spin_op SOp;
  switch (l)
    {
    case 0:				//  (2)
      SOp = T20SS(sys, spin, vect, m);	// T   (i)
      break;				//  0,m

    case 1:				//  (2)
      SOp = T21SS(sys,spin,vect,m,rev);	// T   (i)
      break;				//  1,m

    case 2:				//  (2)
      SOp = T22SS(sys, spin, vect, m);	// T   (i)
      break;				//  2,m
    default:
      spin_T_error(3);
      std::cout << " " << l << "," << m << "\n";
      spin_T_fatality(10);
      break;
    }
  return SOp;				// Return SOp
  }



	// Input		Im    : Spin operator I-
	// 			Iz    : Spin operator Iz
	// 			Ip    : Spin operator I+
	// 			vect  : Cartesian vector
	//			l     : Rank index
	//			m     : Component index
	//			rev   : Reverse tensor order
	// Output		SOp   : Spin Operator for SphT Component T
	//								  lm 

spin_op T2SS(spin_op& Im, spin_op& Iz, spin_op& Ip,
				 const coord& vect, int l, int m, int rev)
  {
  spin_op SOp;
  switch (l)
    {
    case 0:				//  (2)
      SOp = T20SS(Im,Iz,Ip,vect,m);	// T   (i)
      break;				//  0,m

    case 1:				//  (2)
      SOp = T21SS(Im,Iz,Ip,vect,m,rev);	// T   (i)
      break;				//  1,m

    case 2:				//  (2)
      SOp = T22SS(Im,Iz,Ip,vect,m);	// T   (i)
      break;				//  2,m
    default:
      spin_T_error(3);
      std::cout << " " << l << "," << m << "\n";
      spin_T_fatality(10);
      break;
    }
  return SOp;				// Return SOp
  }

spin_op T20(const spin_sys& sys, int spin, const coord& vect, int m)
  {return T20SS(sys,spin,vect,m);}
spin_op T20(spin_op& Im, spin_op& Iz, spin_op& Ip, const coord& vect, int m)
  {return T20SS(Im,Iz,Ip,vect,m);}
spin_op T21(const spin_sys& sys, int spin, const coord& vect, int m)
  {return T21SS(sys,spin,vect,m);}
spin_op T21(spin_op& Im, spin_op& Iz, spin_op& Ip, const coord& vect, int m)
  {return T21SS(Im,Iz,Ip,vect,m);}
spin_op T22(const spin_sys& sys, int spin, const coord& vect, int m)
  {return T22SS(sys,spin,vect,m);}
spin_op T22(spin_op& Im, spin_op& Iz, spin_op& Ip, const coord& vect, int m)
  {return T22SS(Im,Iz,Ip,vect,m);}


	// Input		sys   : Spin system
	// 			spin  : Spin index
	// 			vect  : Cartesian vector
	//			m     : Component index
	//                      rev   : Reverse tensor order (unused)

	// Output		SOp   : Spin Operator for SphT Component T
	//								  0m 

spin_op T20SS(const spin_sys& sys, int spin, const coord& vect, int m, int rev)
  {
  spin_op SOp;
  complex vp,vm;
  double vz;
  vp = complex(vect.x(),vect.y());	// u+ = ux + iuy
  vm = complex(vect.x(),-vect.y());	// u- = ux - iuy
  vz = vect.z();			// uz
  switch (m)
    {
    case 0:				//  (2)            [1]       1 [          ]
      SOp = -1.0/sqrt(3.0)		// T   (i) = - sqrt|-| I u + - |I u + I u |
	     *(Iz(sys,spin)*vect.z() 	//  00             [3]  z z  2 [ - -   + +]
              + 0.5*(Im(sys,spin)*vm
               + Ip(sys,spin)*vp));
      break;
    default:
      spin_T_error(3);
      std::cout << " 0," << m << "\n";
      break;
    }
  return SOp;				// Return SOp
  rev = 0;				// Compiler likes this to be used
  }

	// Input		Im    : Spin operator I-
	// 			Iz    : Spin operator Iz
	// 			Ip    : Spin operator I+
	// 			vect  : Cartesian vector
	//			m     : Component index
	//                      rev   : Reverse tensor order (unused)
	// Output		SOp   : Spin Operator for SphT Component T
	//								  0m 

spin_op T20SS(spin_op& Im, spin_op& Iz, spin_op& Ip, const coord& vect, int m, int rev)
  {
  spin_op SOp;
  complex vp,vm;
  double vz;
  vp = complex(vect.x(),vect.y());	// u+ = ux + iuy
  vm = complex(vect.x(),-vect.y());     // u- = ux - iuy
  vz = vect.z();			// uz
  switch (m)
    {
    case 0:				//  (2)            [1]       1 [          ]
      SOp = -1.0/sqrt(3.0)*(Iz*vect.z() // T   (i) = - sqrt|-| I u + - |I u + I u |
                   + 0.5*(Im*vm+Ip*vp));//  00             [3]  z z  2 [ - -   + +]
      break;
    default:
      spin_T_error(3);
      std::cout << " 0," << m << "\n";
      break;
    }
  return SOp;				// Return SOp
  rev = 0;				// Compiler likes this to be used
  }

	// Input		sys   : spin system
	// 			spin  : spin index
	// 			vect  : Cartesian vector
	//			m     : component index
	//			rev   : Reverse tensor order
	// Output		SOp   : Spin Operator for SphT Component T
	//								  1m 
	// Note			      : Assumes the first component is
	//				the vector vect

spin_op T21SS(const spin_sys& sys, int spin, const coord& vect, int m, int rev)
  {
  spin_op SOp;
  complex vp,vm;
  double vz;
  double fact = -0.5;
  if(rev) fact = 0.5;
  vp = complex(vect.x(),vect.y());	// u+ = ux + iuy
  vm = complex(vect.x(),-vect.y());     // u- = ux - iuy
  vz = vect.z();			// uz
  switch (m)
    {
    case -1:				//  (2)         1 [           ]
      SOp = fact*(Im(sys,spin)*vz	// T    (i) = - - |I u  - I u |
                     - Iz(sys,spin)*vm);//  1,-1        2 [ - z    z -]
      break;
    case 0:				//  (2)            [1][          ]
      SOp = -1.0/sqrt(8.0)		// T   (i) = - sqrt|-||I u - I u |
                 *(Ip(sys, spin)*vm	//  1,0            [8][ + -   - +]
                     - Im(sys,spin)*vp);
      break;
    case 1:				//  (2)        1 [           ]
      SOp = fact*(Ip(sys,spin)*vz	// T   (i) = - - |I u  - I u |
                     - Iz(sys,spin)*vp);//  1,1        2 [ + z    z +]
      break;
    default:
      spin_T_error(3);
      std::cout << " 1," << m << "\n";
      break;
    }
  return SOp;				// Return SOp
  }



	// Input		Im    : Spin operator I-
	// 			Iz    : Spin operator Iz
	// 			Ip    : Spin operator I+
	// 			vect  : Cartesian vector
	//			m     : Component index
	//			rev   : Reverse tensor order
	// Output		SOp   : Spin Operator for SphT Component T
	//								  1m 

spin_op T21SS(spin_op& Im, spin_op& Iz, spin_op& Ip,
					 const coord& vect, int m, int rev)
  {
  spin_op SOp;
  double fact = -0.5;
  if(rev) fact = 0.5;
  complex vp,vm;
  double vz;
  vp = complex(vect.x(),vect.y());	// u+ = ux + iuy
  vm = complex(vect.x(),-vect.y());     // u- = ux - iuy
  vz = vect.z();			// uz
  switch (m)
    {
    case -1:				//  (2)         1 [           ]
      SOp = fact*(Im*vz - Iz*vm);	// T    (i) = - - |I u  - I u |
      break;				//  1,-1        2 [ - z    z -]

    case 0:				//  (2)            [1][          ]
      SOp = -1.0/sqrt(8.0)		// T   (i) = - sqrt|-||I u - I u |
                 *(Ip*vm - Im*vp);	//  1,0            [8][ + -   - +]
      break;
    case 1:				//  (2)        1 [           ]
      SOp = fact*(Ip*vz - Iz*vp);	// T   (i) = - - |I u  - I u |
      break; 				//  1,1        2 [ + z    z +]
    default:
      spin_T_error(3);
      std::cout << " 1," << m << "\n";
      break;
    }
  return SOp;				// Return SOp
  }

	// Input		sys   : Spin system
	// 			spin  : Spin index
	// 			coord : Cartesian vector
	//			m     : Component index
	//                      rev   : Reverse tensor order (unused)
	// Output		SOp   : Spin Operator for SphT Component T
	//								  2m 

spin_op T22SS(const spin_sys& sys, int spin, const coord& vect, int m, int rev)
  {
  spin_op SOp;
  complex vp,vm;
  double vz, vy, vx;
  vp = complex(vect.x(),vect.y());	// u+ = ux + iuy
  vm = complex(vect.x(),-vect.y());     // u- = ux - iuy
  vz = vect.z();			// uz
  vx = vect.x();
  vy = vect.y();
  switch (m)
    {
    case -2:				//  (2)       1
      SOp = 0.5*Im(sys,spin)*vm;	// T    (i) = - I * u
      break;				//  2,-2      2  -   -

    case -1:				//  (2)       1 [           ]
      SOp = 0.5*(Iz(sys,spin)*vm	// T    (i) = - | I u + I u |
                    + Im(sys,spin)*vz);	//  2,-1      2 [  z -   - z]
      break;
    case 0:				//  (2)          [1][          ]
      SOp = 1.0/sqrt(6.0)*		// T   (i) = sqrt|-||3I u - I.u|
              (3.0*Iz(sys,spin)*vz	//  2,0          [6][  z z     ]
                -(Ix(sys,spin)*vx
                   + Iy(sys,spin)*vy
                    + Iz(sys,spin)*vz));
      break;
    case 1:				//  (2)        1 [           ]
      SOp = -0.5*(Iz(sys,spin)*vp	// T   (i) = - - | I u + I u |
                    + Ip(sys,spin)*vz);	//  2,1        2 [  z +   + z]
      break;
    case 2:				//  (2)      1
      SOp = 0.5*Ip(sys,spin)*vp;	// T   (i) = - I * u
      break;				//  2,2      2  +   +
    default:
      spin_T_error(3);
      std::cout << " 2," << m << "\n";
      break;
    }
  return SOp;				// Return SOp
  rev = 0;				// Compiler likes this to be used
  }

	// Input		Im    : Spin operator I-
	// 			Iz    : spin operator Iz
	// 			Ip    : Spin operator I+
	// 			vect  : Cartesian vector
	//			m     : component index
	//                      rev   : Reverse tensor order (unused)
	// Output		SOp   : Spin Operator for SphT Component T
	//								  2m 

spin_op T22SS(spin_op& Im, spin_op& Iz, spin_op& Ip,
				             const coord& vect, int m, int rev)
  {
  spin_op SOp;
  complex vp,vm;
  double vz, vy, vx;
  vp = complex(vect.x(),vect.y());	// u+ = ux + iuy
  vm = complex(vect.x(),-vect.y());     // u- = ux - iuy
  vz = vect.z();			// uz
  vx = vect.x();
  vy = vect.y();
  switch (m)
    {
    case -2:				//  (2)       1
      SOp = 0.5*Im*vm;			// T    (i) = - I * u
      break;				//  2,-2      2  -   -

    case -1:				//  (2)       1 [           ]
      SOp = 0.5*(Iz*vm + Im*vz);	// T    (i) = - | I u + I u |
      break; 				//  2,-1      2 [  z -   - z]

    case 0:				//  (2)          [1][       1             ]
      SOp = 1.0/sqrt(6.0)*		// T   (i) = sqrt|-||2I u - - (I v + I v )|
              (2.0*Iz*vz		//  2,0          [6][  z z  2   m m   p p ]
            -0.5*(Im*vm + Ip*vp));
      break;
    case 1:				//  (2)        1 [           ]
      SOp = -0.5*(Iz*vp	+ Ip*vz);	// T   (i) = - - | I u + I u |
      break;				//  2,1        2 [  z +   + z]

    case 2:				//  (2)      1
      SOp = 0.5*Ip*vp;			// T   (i) = - I * u
      break;				//  2,2      2  +   +
    default:
      spin_T_error(3);
      std::cout << " 2," << m << "\n";
      break;
    }
  return SOp;				// Return SOp
  rev = 0;				// Compiler likes this to be used
  }
 

// ____________________________________________________________________________
//                           SPIN TENSOR FUNCTIONS
// ____________________________________________________________________________

	// Input		SphT : Spin Tensor(this)
	// 			l    : momentum index
	// 			m    : momentum index
	// Output		SOp  : Spin Operator, the l,m
	//			     : component of Spin Tensor SphT

spin_op T_comp(spin_T& SphT, int l, int m) { return SphT.component(l,m); }
 
spin_op spin_T::component(int l, int m)
  {
  spin_op SOp;
//  int comp = 0;
  if(l > rank || l < 0)
    std::cout << "\n\tRequested Spin Tensor component l=" << l << " on rank "
	 << rank << " tensor";
  else if(abs(m) > l)
    std::cout << "\n\tRequested Spin Tensor component l=" << l << " and m="
	 << m << " ?";
// sosi 11/11/92 replaced - else - with the else if line below
  else if(pr[l]) SOp = *(pr[l][l-m]);
  return SOp;				// Return SOp
  }


/*                       l1  l2
                        --- ---
                  (R)   \   \         
                 T    = /   /   T1      T2      <l1l2m1m2|LM>
                  L,M   --- ---   l1,m1   l2,m2
                         m1  m2
   where
  	 L = [l1+l2, |l1-l2|],  M = m1 + m2, & <l1l2m1m2|LM> = C-G coefficient

	   Input		SphT1 : Irreducible Spin Tensor
	   			SphT2 : Irreducible Spin Tensor
	  			L     : Spin Tensor Component Index
	  			M     : Spin Tensor Component Index
	   Output		SOp   : Spin operator which is the L,M
	  				component of tensor SphT, the
	  			        product of the two input irreducible
	  				tensors	SphT = SphT1 x SphT2         */
 
spin_op T_mult(spin_T& SphT1, spin_T& SphT2, int L, int M)
  {
  spin_op SOp;					// Empty spin operator
  int l1 = SphT1.rank;				// Get l1 value
  int l2 = SphT2.rank;				// Get l2 value
  int span1 = 2*l1+1;				// Range of m1
  int span2 = 2*l2+1;				// Range of m2
  int m1 = 0;					// Initial m1
  int m2 = 0;					// Initial m2

	//Changed long(2*L+1) to double(2*L+1)
  double norm = sqrt(double(2*L+1));			// For Wigner 3J use over CG

	// Changed long(-1.0) to double(-1.0)
	double prefact = 				// For Wigner 3J use over CG
   pow(double(-1.0), abs(2*l2+L-M));
  
	for(int M1=0; M1<span1; M1++)			// Loop over m1 values
    for(int M2=0; M2<span2; M2++)		// Loop over m2 values
      {
      m1 = l1-M1;				// Calculate actual m1
      m2 = l2-M2;				// Calculate actual m2
      if((m1 + m2) == M)			// Restrict sum M = m1 + m2
        SOp += T_comp(SphT1, l1, m1)
              * T_comp(SphT2, l2, m2)
		* Wigner_3j(L,l1,l2,M,-m1,-m2);
      }
    SOp *= prefact*norm;			// Scale the result
  return SOp;					// Return SOp
  }


/*                       l1  l2
                        --- ---
                  (R)   \   \         
                 T    = /   /   T1      T2      <l1l2m1m2|LM>
                  L,M   --- ---   l1,m1   l2,m2
                         m1  m2
   where
  	 L = [l1+l2, |l1-l2|],  M = m1 + m2, & <l1l2m1m2|LM> = C-G coefficient

	   Input		SphT1 : Irreducible Spin Tensor
	   			SphT2 : Irreducible Spin Tensor
	   Output		SphT  : Spin Tensor which is the product
	  			        of the two irreducible input tensors
	  				SphT = SphT1 x SphT2                 */

spin_T T_mult(spin_T& SphT1, spin_T& SphT2)
  {
  spin_T SphT(*SphT1.sys);
  spin_op SOp;					// Working spin operator
  if(SphT1.pr && SphT2.pr)			// Check that both tensors exist		
    if(SphT1.pr[SphT1.rank]
                       && SphT2.pr[SphT2.rank])		
      {
      SphT.rank = SphT1.rank + SphT2.rank;	// Rank = l1 + l2
      //SphT.pr = new spin_op **[SphT1.rank+1];	// work around bug in g++
      typedef spin_op **tmp;	   		//   "     "    "   "  "
      SphT.pr = new tmp [SphT.rank+1]; 		//   "     "    "   "  "
      for(int l=0; l<=SphT.rank; l++)		// Loop through all L 
        {
        int span = 2*l+1;
        SphT.pr[l] = new spin_op *[span];	// Allocate rank l storage
        for (int m=0; m<span; m++)		// Compute all l components
          {
          SOp = T_mult(SphT1, SphT2, l, l-m);
          SphT.pr[l][m] = new spin_op(SOp);
          }
        }
      }
  return SphT;					// Return SphT
  }
 
// sosi - need the above two functions for blending spatial and spin tensors!

// ____________________________________________________________________________
//                           SPIN TENSOR ROTATIONS
// ____________________________________________________________________________

// sosi Next two functions should be removed/replaced by the member function!

spin_T T_rot(spin_T& SphT1, double alpha, double beta, double gamma)
{ return SphT1.rotate(alpha, beta, gamma); }

spin_op T_rot(spin_T& SphT1, int l, int m, double alpha, double beta, double gamma)
  { return SphT1.rotate(l, m, alpha, beta, gamma); }



	// Input		SphT1 : Spin Tensor(this)
	// 			alpha : Euler Angle
	// 			beta  : Euler Angle
	// 			gamma : Euler Angle
	// Output		SphT  : Spin Tensor which is the input
	//			        tensor in the coordinate system
	//				rotated by input Euler angles
 
spin_T spin_T::rotate(double alpha, double beta, double gamma)
  {
  spin_T SphT(*((*this).sys));
  spin_op SOp;
  if(pr)					// rotate only non-NULL SphT1
    {
    //SphT.pr = new spin_op **[SphT1.rank+1];	// work around bug in g++
    typedef spin_op **tmp;	   		// work around bug in g++
    SphT.pr = new tmp [rank+1];			// work around bug in g++
    for(int l=0; l<=rank; l++)
      {
      if(pr[l])					// rotate only non-Null SphT1
        {					//                           l
        int span = 2*l+1;
        SphT.pr[l] = new spin_op *[span];
        for (int m=0; m<span; m++)
          {
          SOp = (*this).rotate(l, l-m, alpha, beta, gamma);
          SphT.pr[l][m] = new spin_op(SOp);
          }
        }
      else
        SphT.pr[l] = NULL;
      }
    }
  SphT.rank = rank;
  return SphT;					// Return SphT
  }


spin_T spin_T::rotate(const coord& EA) { return rotate(EA.x(), EA.y(), EA.z()); }

	// Input		SphT1 : Spin Tensor(this)
	// 			EA    : Set of Euler angles
	// Output		SphT  : Spin Tensor which is the input
	//			        tensor in the coordinate system
	//				rotated by input Euler angles



	// Input		SphT1 : Spin tensor(this)
	//			l     : Component index
	//			m     : Component index
	// 			alpha : Euler Angle
	// 			beta  : Euler Angle
	// 			gamma : Euler Angle
	// Output		Tlm   : Spatial Tensor component of
	//			        the input tensor in the coordinate
	//				system rotated by input Euler angles
 
spin_op spin_T::rotate(int l, int m, double alpha, double beta, double gamma)
  {
  spin_op SOp;
  int span = 2*l+1;
  for (int n=0; n<span; n++)
    SOp += *(pr[l][m]) * DJ(l, l-n, m, alpha, beta, gamma);
  return SOp;
  }

	// Input		SphT1 : Spin Tensor(this)
	//			l     : Component index
	//			m     : Component index
	// 			EA    : Set of Euler angles
	// Output		Tlm   : Spatial Tensor component of
	//			        the input tensor in the coordinate
	//				system rotated by input Euler angles

spin_op spin_T::rotate(int l, int m, const coord& EA)
  { return rotate(l, m, EA.x(), EA.y(), EA.z()); }

/*
                                            m                  m            
        A(l,m).T(l,m) = T(l,m).A(l,m) = (-1) * T  * A    = (-1) * A  * T    
                                                lm   l-m           lm   l-m 

	   Input		SphT  : Irreducible Spin Tensor
	   			SphA  : Irreducible Spatial Tensor
	  			l     : Tensor Rank Index
	  			m     : Tensor Component Index
	   Output		SOp   : Spin operator which is the result
	  				of the scalar product of  the l,m 
	  			        component ot SphT with the l-m of SphA
	   Note			      : Both function forms should produce 
	  				equivalent results                   */

spin_op T_prod(spin_T& SphT, space_T& SphA, int l, int m)
	{ 
	// changed long(-1.0) to double(-1.0)
	return pow(double(-1.0),abs(m)) * T_comp(SphT,l, m) * SphA.component(l, -m); 
	}

spin_op T_prod(space_T& SphA, spin_T& SphT, int l, int m)
	{ 
	return pow(double(-1.0), abs(m)) * SphA.component(l, m) * T_comp(SphT, l, -m); 
	}


/*
                         Sum  [     m            ]   Sum  [     m            ]
 A(l).T(l) = T(l).A(l) = over | (-1) * T  * A    | = over | (-1) * A  * T    |
                          m   [         lm   l-m ]    m   [         lm   l-m ]

	   Input		SphT  : Irreducible Spin Tensor
	   			SphA  : Irreducible Spatial Tensor
	  			l     : Tensor Rank Index
	   Output		SOp   : Spin operator which is the result
	  				of the scalar product of  the rank l 
	  			        components of SphT with SphA
	   Note			      : Both function forms should produce 
	  				equivalent results                   */
 
spin_op T_prod(spin_T& SphT, space_T& SphA, int l)
  {
  spin_op SOp;					// Empty spin operator
  int span = 2*l+1;				// Span of m this rank
  for(int m=0; m<span; m++)			// Loop over m values
    SOp += T_prod(SphT, SphA, l, l-m); 		// Sum product
  return SOp;					// Return spin operator
  }
 
spin_op T_prod(space_T& SphA, spin_T& SphT, int l)
  {
  spin_op SOp;					// Empty spin operator
  int span = 2*l+1;				// Span of m this rank
  for(int m=0; m<span; m++)			// Loop over m values
    SOp += T_prod(SphA, SphT, l, l-m); 		// Sum product 
  return SOp;					// Return spin operator
  }

/*
             Sum  Sum  [     m             ]   Sum  Sum  [     m             ]
 A.T = T.A = over over | (-1)  * T  * A    | = over over | (-1)  * A  * T    |
              l    m   [          lm   l-m ]    l    m   [          lm   l-m ]

	   Input		SphT  : Irreducible Spin Tensor
	   			SphA  : Irreducible Spatial Tensor
	   Output		SOp   : Spin operator which is the result
	  				of the scalar product of  the rank l 
	  			        components of SphT with SphA
	   Note			      : The functions should produce 
	  				equivalent results                   */
 
spin_op T_prod(spin_T& SphT, space_T& SphA)
  {
  spin_op SOp;					// Empty spin operator
  for(int l=0; l<=SphT.rank; l++)		// Loop over rank of A & T
    SOp += T_prod(SphT, SphA, l); 		// Sum tensor product this l
  return SOp;					// Return spin operaotr
  }

spin_op T_prod(space_T&SphA, spin_T& SphT)
  {
  spin_op SOp;					// Empty spin operator
  for(int l=0; l<=SphT.rank; l++)		// Loop over rank of A & T
    SOp += T_prod(SphA, SphT, l); 		// Sum tensor product this l	
  return SOp;					// Return spin operator
  }

// ____________________________________________________________________________
//                       AUXILIARY SPIN TENSOR FUNCTIONS
// ____________________________________________________________________________


int spin_T::Rank() { return rank; }

	// Input		SphT  : Spin Tensor (this)
	// Output		r     : Tensor rank
 



	// Input		a     : momentum index
	// 			b     : momentum index
	// 			alpha : z-momentum index
	// 			beta  : z-momentum index
	// 			c     : momentum index
	// 			gamma : z-momentum index
	// Output		r     : Clebsch-Gordan coefficients
	// Note 		      : Only Integral (Non-Negative) J
 
double Clebsch_Gordan(int a, int b, int alpha, int beta, int c, int gamma)
  {
  double r = 0.0;
  double delabc = 0.0;
  double term1 = 0.0;
  double term2 = 0.0;
  double term3 = 0.0;
  int nu = 0;
  if((alpha + beta - gamma) == 0)
    {
    delabc = fact(a+b-c) * fact(a+c-b) * fact(b+c-a);
    delabc = sqrt(delabc/fact(a+b+c+1));
    term1 = (2*c+1)*fact(a+alpha)*fact(a-alpha)*fact(b+beta);
    term1 *= fact(b-beta)*fact(c+gamma)*fact(c-gamma);
    term1 = sqrt(term1);
    while((nu <= (a-alpha)) && (nu <= (b+beta)) && (nu <= (a+b-c)))
      {
      if((c-b+alpha+nu >= 0) && (c-a-beta+nu >= 0))
        {
        term3 = fact(a-alpha-nu) * fact(c-b+alpha+nu) * fact(b+beta-nu);
        term3 *= fact(c-a-beta+nu) * fact(nu) * fact(a+b-c-nu);
				// changed long(-1.0) to double(-1.0)
        term2 += pow(double(-1.0), abs(nu))/term3;
        }
      nu++;
      }
    r += delabc*term1*term2;
    }
  return r;
  }


//double Wigner3J_half_int(const int J, const int m, const int n, double beta )

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
 
//{
// beta = beta*PI/360.0;
// double r = 0.0;
// double COS = 0.0;
// double SIN = 0.0;
// int num = 0;
// int denom = 0;
// int t=0;   
// int actJ = 2*J - 1;
// int actm = 2*abs(m) - 1;
// int actn = 2*abs(n) - 1;
//   if(m < 0)
//     actm = -actm;
//   if(n < 0)
//     actn = -actn;
// int addJm = J+m;
// if((J > 0) && (actm <= actJ) && (actn <= actJ) && (m != 0) && (n != 0) )
//   {
//   while((t <= (actJ+actm)/2) && (t <= (actJ-actn)/2))
//     {
//     if((t + (actn-actm)/2) >= 0)
//       {
//       COS = pow(cos(beta), (actJ + (actm - actn)/2 - 2*t));
//       SIN = pow(sin(beta), (2*t + (actn - actm)/2));
//       num = fact((actJ+actm)/2) * fact((actJ-actm)/2);
//       num *= fact((actJ+actn)/2) * fact((actJ-actn)/2);
//       denom = fact((actJ + actm)/2 - t) * fact((actJ - actn)/2 - t);
//       denom *= fact(t) * fact(t + (actn - actm)/2);
//       r += pow(-1.0,t)*sqrt(num)*COS*SIN/denom;
//       }
//     t++;
//     }
//   }
// else
//   {
//     cout << "\nSpatial Function:          (" << 2*J-1 << "/2)"
//          << "\nSpatial Function: Unknown d"
//          << "\nSpatial Function:          " << m << "/2," << n << "/2";
//     spatial_fatality(10);
//   }
// return r;
//}



	// Input		a     : momentum index
	// 			b     : momentum index
	// 			c     : momentum index
	// 			alpha : z-momentum index
	// 			beta  : z-momentum index
	// 			gamma : z-momentum index
	// Output		r     : Wigner 3-j coefficients
 
double Wigner_3j(int a, int b, int c, int alpha, int beta, int gamma)
  {
	// changed long(-1.0) to double(-1.0)
	// also changed sqrt(2*c+1) to sqrt(double(2*c+1)
  double r = pow(double(-1.0),abs(a-b-gamma))/sqrt(double(2*c+1));
  return r*Clebsch_Gordan(a, b, alpha, beta, c, -gamma);
  }


// ____________________________________________________________________________
//                       SPIN TENSOR I/O FUNCTIONS
// ____________________________________________________________________________

	// Input		ostr : string
	// 			SphT : spin tensor
	// Return		     : stream, prints spin tensor components

std::ostream& operator<< (std::ostream& ostr, const spin_T& SphT)
  {
  if(!SphT.pr)						// If tensor is NULL
    {							// quick print/return
    ostr << "\n\tSpin Tensor is Currently NULL\n";	
    return ostr; 
    }
  int span = 0;						// Span of m per l
  for(int l=0; l<=SphT.rank; l++)			// Loop over all ranks
    {
    if(SphT.pr[l])					// Print T  if rank !0
    span = 2*l+1;					//        l
    for(int m=0; m<span; m++)
      {
      if(SphT.pr[l][m])					// Print T  if SOp !NULL
        {						//        2
        ostr << "\n\tT"
             << "\n\t " << l << "," << l-m
             << "\n" << *(SphT.pr[l][m]);
        }
      else
        ostr << "\n\tT   = 0"
             << "\n\t " << l << "," << l-m << "\n";
      }
    }
  return ostr; 
}

#endif 						// SpinT.cc
