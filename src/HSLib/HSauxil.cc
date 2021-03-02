/* HSauxil.cc ******************************************-*-c++-*-*
**								**
** 	                     G A M M A				**
**								**
**	NMR Library	                Implementation  	**
**							 	**
**	Copyright (c) 1991, 1992			 	**
**	Scott Smith				 		**
**	Eidgenoessische Technische Hochschule	 		**
**	Labor fuer physikalische Chemie		 		**
**	8092 Zurich / Switzerland		 		**
**						 		**
**      $Header: $
**							 	**
*****************************************************************/

/*****************************************************************
**							 	**
** Description						 	**
**							 	**
** The NMR Library provides functions for the simulation of	**
** magnetic resonance experiments and associated mathematical	**
** capabilities. Most of the library exists in modules called	**
** nmr_"name", where "name" indicates the type of functions 	**
** contained therein.  Included here are functions of a	more	**
** general nature not yet associated with a group of routines.	**
**							 	**
*****************************************************************/

#ifndef _nmrlib_cc_			// Is this file already included?
#define _nmrlib_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#endif

#include <HSLib/HSauxil.h>		// Include the header 
#include <Basics/Gutils.h>		// Include query_parameter 
#include <HSLib/SpinSystem.h>		// Include spin systems
#include <HSLib/SpinOp.h>		// Include spin operators
#include <HSLib/SpinOpCmp.h>		// Include composite spin ops
#include <Basics/Gutils.h>		// Include error handling
#include <Basics/StringCut.h>		// Include string parsing
#include <stdlib.h>
#include <iostream>
#include <cmath>			// Inlcude HUGE_VAL_VAL

// ____________________________________________________________________________
//                           DENSITY MATRIX FUNCTIONS
// ____________________________________________________________________________


	// Input	       sys : A spin system
	// Output		Op : The density matrix describing
	//			     the system under a high-temp
	//			     approximation.
	// Note			   : The identity matrix is subtracted
	//			     out, leaving the operator traceless
	//			      (normally Tr{density matrix} = 1)
	// Note			   : The density matrix is scaled by a
	//			     constant scaling factor.  This does
	//			     not change the relative intensities
	//			     of any observed transitions

/* By definition

               1    [ -H ]   1 [          2    3        ]           -H
     sigmaeq = - exp|----] = - | 1 + X + X  + X  + .... | where X = --
               Z    [ kT ]   Z [                        ]           kT

   When kT >> -H, the identity matrix is neglected, & the factor ZkT removed,

             ~ 1 [     -H ]     -H
     sigmaeq = - | 1 + -- | --> --- --> -H 
               Z [     kT ]     ZkT

   For the isotropic static Hamiltonian in NMR, neglecting chemical shifts and
   coupling constants (e.g. energy level populations are insignificantly 
   affected by these contributions, so H is the Zeeman Hamiltonian)

                 ---                           ---  gamma
             ~   \                             \         i     homonuclear
     sigmaeq = - /   - hbar * gamma  * Iz  --> /    ------ Iz  -----------> Fz
                 ---               i     i     ---  gamma    i
                  i                             i        0

   where the operator has be again rescaled by hbar*gamma
  							 0                   */

gen_op sigma_eq(const spin_sys& sys)
  {
  if(!sys.spins()) return gen_op();		// No system? Return null Op
  if(sys.homonuclear()) 			// If the system is Homonuclear
    { 						// just return Fz
    gen_op Op(Fz(sys), sys.get_basis()); 
    Op.name("Equilibrium Density");
    return Op;
    }
  gen_op Op = Iz(sys,0);			// For Heteronuclear systems,
  double gamma0 = sys.gamma(0);			// must return summed Iz scaled
  for(int i=1; i<sys.spins(); i++)		// by spin gyromagnetic ratio
    Op += (sys.gamma(i)/gamma0)*Iz(sys,i);	// for proper spin polarization
  Op.name("Equilibrium Density");
  return Op;
  } 

gen_op sigma_eq(const spin_sys& sys, const Isotope& I)
  {
  if(!sys.spins()) return gen_op();		// No system? Return null Op
  double gamma0 = I.gamma();			// Gyromagnetic ratio of I
  double sf;					// Scaling factor
  if(sys.homonuclear()) 			// If system is homonuclear
    { 						// just return a scaled Fz,
    sf = sys.gamma(0)/gamma0;			// scaled by gyromagnetic
    return gen_op(Fz(sys), sys.get_basis());	// ratios
    }
  gen_op Op;					// For Heteronuclear systems,
  for(int i=0; i<sys.spins(); i++)		// by spin gyromagnetic ratio
    Op += (sys.gamma(i)/gamma0)*Iz(sys,i);	// for proper spin polarization
  return Op;
  } 




//gen_op sigma_eq(const spin_sys& sys, basis& bs)
//{
//  spin_op Z = Fz(sys);
//  Op = gen_op(Z.matrix(), bs);
//  return;
//} 




// ____________________________________________________________________________
//                        PROPAGATOR FUNCTIONS
// ____________________________________________________________________________

//              These now reside in HSprop I believe 

/*

gen_op prop(gen_op &ham, const double time)
  return U

	//	 		ham   : "Hamiltonian" for propagation (in Hertz)
	//	 		time  : Evolution time (seconds)
	// Output		U     : Propagator for the input Hamiltonian

  {
  ham.set_EBR();			// Put ham into its eigenbase
  complex z(0,-2.0*PI*time);		// Exponential factor (rad/Hz)
  U = ham.exp(z);			// Generate evolution propagator
  return;
  }


  void prop_ip(gen_op &U, const double time)

	//	 		ham   : "Hamiltonian" for propagation (in Hertz)
	//	 		time  : Evolution time (seconds)
	// Output		U     : Propagator for the input Hamiltonian
	// Note			      : Done in place, overwriting ham

  {
  U.set_EBR();				// Put ham into its eigenbase
  complex z(0,-2.0*PI*time);		// Exponential factor (rad/Hz)
  U = U.exp(z);				// Generate evolution propagator
  return;
  }

*/

// ____________________________________________________________________________
//                           COHERENCE SELECTION
// ____________________________________________________________________________


void zero_mqc(const spin_sys &sys, gen_op &Op, int order=-1, int type=-1)

	// Input	sys	: Spin system
	// Input	Op	: General operator (associated to sys)
	//		order	: Coherence order selected
	//			    <0 - populations (default)
	//			     0 - zero quantum coherence
	//		             1 - single quantum coherence
	//			     n - n quantum coherence		
	//		type	: Type of zeroing
	//			    <0 - zero coherences except "order" (DEF)
	//			     0 - zero coherences of "order"
	//			     1 - zero coherences >= "order"
	//			    >1 - zero coherences <= "order"
	// Output	none	: Function is void.  Op has had
	//			  elements associated with one or
	//			  more coherence types zeroed

  {
  Op.set_DBR();				// Set default basis for easy selection
  if(order > 2*sys.qn())                // Insure order is reasonable, it can't
     order = int(2*sys.qn());           // exceed the maximum coherence order!
  if(type < 0) type = -1;               // Insure that the type is proper, it
  if(type > 2) type = 2;                // must be in the range [-1, 1]
  int size = sys.HS(), i, j;            // System Hilbert space dimension
  if(size != Op.dim())
    GAMMAerror("zero_mqc", "Hilbert Space Mismatch");
  col_vector qnstates;			// Total Fz values of basis functions 
  if(order >= 0)			// Deal with "order" quantum coherence
    {
    qnstates = sys.qnStates(); 		// Get a list of elements vs. coherence
    switch(type)
      {
      default:
      case(-1):				// All coherences zeroed except "order"
	for(i=0; i<size; i++)
	  for(j=0; j<size; j++)
	    if(abs(int(Re(qnstates(i)-qnstates(j))))!=order)
	      Op.put(complex0,i,j);
	break;
      case(0):				// All coherences of "order" to zero
	for(i=0; i<size; i++)
	  for(j=0; j<size; j++)
	    if((i != j)&(abs(int(Re(qnstates(i)-qnstates(j))))==order))
	      Op.put(complex0,i,j);
	break;
      case(1):				// All coherences >= "order" to zero
	for(i=0; i<size; i++)
	  for(j=0; j<size; j++)
	    if((i != j) & (abs(int(Re(qnstates(i)-qnstates(j))))>=order))
	      Op.put(complex0,i,j);
	break;
      case(2):				// All coherenceces <= "order" and
	for(i=0; i<size; i++)		// populations to zero
	  for(j=0; j<size; j++)
	    if(abs(int(Re(qnstates(i)-qnstates(j))))<=order)
	      Op.put(complex0,i,j);
	break;
      }
    }
  else                                  // Here if dealing with populations
    {
    if(type == 0)                               // Zero ONLY populations
      for(i=0; i<size; i++)                     // Just set all diagonal
        Op.put(complex0,i,i);                   // elements to zero
    else					// Zero ALL BUT populations
      {
      matrix mx(size, size, d_matrix_type);     // Op matrix will be diagonal
      for(i=0; i<size; i++)                     // Copy the Op diagonal
        mx.put(Op.get(i,i),i,i);                // elements into its new array
      Op = gen_op(mx, Op.get_basis());          // Reform the operator
      }
    }
  return;
  }

  
// ____________________________________________________________________________
// ************************ Single Transition Operators ***********************
// ____________________________________________________________________________


gen_op st_Op(gen_op &Ham, int lev1, int lev2, char axis)

	// Input		Ham  : General operator (Hamiltonian).
	// 			lev1 : Energy level 1
	// 			lev2 : Energy level 2
	// 			type : Transition operator type
	//				x = Ix 		p = I+
	//				y = Iy 		m = I-
	//				z = Iz
	// Return		Op   : General operator which is a
	//			       single transition operator
	//			       for the transtion between levels
	//			       lev1 and lev2.
	// Note			     : Result EXCLUSIVELY in EBR of Ham
	// Note			     : See EBW, page 34.

  {
  int hs = Ham.dim();
  if(!hs) 				// Check for NULL Op
    {
    std::string hdr("Single Transition Operator");
    GAMMAerror(hdr,"Input Hamiltonian NULL",1);
    GAMMAerror(hdr,"Cannot Generate Operator",1);
    std::cout << "\n\n";
    exit(-1);
    }
  matrix mx(hs,hs,0);			// For return op rep
  Ham.set_EBR();			// Put Hamiltonian in EBR
  switch (axis)
    {
    case 'x':				// Ix operator
      mx.put(0.5,lev1,lev2); 
      mx.put(0.5,lev2,lev1); 
      break;
    case 'y':				// Iy operator
      mx.put(complex(0,-0.5),lev1,lev2); 
      mx.put(complex(0,0.5),lev2,lev1); 
      break;
    case 'z':				// Iz operator
      mx.put(0.5,lev1,lev1); 
      mx.put(-0.5,lev2,lev2); 
      break;
    case 'p':				// I+ operator
      mx.put(1.0,lev1,lev2); 
      break;
    case 'm':				// I- operator
      mx.put(0.0,lev1,lev2); 
      mx.put(1.0,lev2,lev1); 
      break;
    }
  return gen_op(mx, Ham.get_basis());
  }


void sqt_v(gen_op &Ham)

	// Input		Ham  : General operator (Hamiltonian).
	// Return		none : All transitions associated with
	//			       the input Hamiltonian are sent to
	//			       standard output
	// Note			     : Sets Ham to its EBR

  {
  int size = Ham.dim( );
  if(size) 				// Check for NULL Op
    {
    Ham.set_EBR();			// Put Hamiltonian in EBR
    for(int i=0; i<size; i++)
      for(int j=i+1; j<size; j++)
         std::cout << "\n" << i << " --> " << j << " : " << Ham(i,i) - Ham(j,j);
    }
  std::cout << "\n";
  return;
  }


// ____________________________________________________________________________
//                         VARIOUS USEFUL FUNCTIONS
// ____________________________________________________________________________

/*
int read_units (string xstr, double &x)

	// Input		xstr	: string containing input value
        //                                of the format
        //                                  double unit
        //                                with unit in {GMKmunpf}
	//			x       : actual value
	// Output		int	: TRUE if time interpreted O.K.
	//				  FALSE if not
	// Note				: x is set (if TRUE) as written in xstr
	// Note				: error message if x not determined
	// Note				: default units are unscaled
  
  {
  cutWhite(xstr);				// Remove beginning white space
  string val = cutDouble(xstr);			// Cut double from string
  cut(xstr,RXdouble);
  string unit = cut(xstr,Regex("[GMKmunpf]"));
  if ((xstr=="")&(val!=""))
    {
      double t = atof(val);
      if (unit!="") 
	{
	  switch (unit[0])
	    {
	    case 'G':
	      x = t*1.e9;	// Giga units
	      break;
	    case 'M':
	      x = t*1.e6;	// Mega units
	      break;
	    case 'K':
	      x = t*1.e3;	// Giga units
	      break;
	    case 'm':
	      x = t*1.e-3;	// Mega units
	      break;
	    case 'u':
	      x = t*1.e-6;	// Mega units
	      break;
	    case 'n':
	      x = t*1.e-9;	// Mega units
	      break;
	    case 'p':
	      x = t*1.e-12;	// Mega units
	      break;
	    case 'f':
	      x = t*1.e-15;	// Mega units
	      break;
	    }
	}
      else 
	x = t*1.0e-6;
      return 1;
    }
  else  return 0;
}
*/

int* sort_super_op_basis (const spin_sys& sys)
  {
  return sort_LOp_basis(sys);
  }


int* sort_LOp_basis (const spin_sys& sys)
//   return index;

//	Input			sys	: Spin system
//	Output			index	: An integer vector containing
//					  the basis function order which
//					  has coherence sorting.
//	Note				: Output order has populations 1st,
//					  followed by ZQC, then SQC, ...
//					  up to the highest coherenece order
//					  of the spin system
//	Note				: The default basis order is 
//					  0,1,2,...,LS-1 and this is
//					  usually not coherence ordered due
//					  to the tensor products used in
//					  formulating super operators

{
  int HS = sys.HS();
  int* index;
  index = new int[HS*HS];
  int mqc = int(2*sys.qn());	// maximum coherence order
  int* pos;
  pos = new int[2*mqc+1];
  int i,j,p,c;

  row_vector coher=sys.CoherDist();
  col_vector states=sys.qnStates();

  pos[0]=HS;			// start of ZQC 
  pos[1]=int(Re(coher(mqc)));		// start of (+1)QC

  for (i=1; i<mqc; i++)
    {
      pos[2*i]   = pos[2*i-1]+int(Re(coher(i+mqc)));
      pos[2*i+1] = pos[2*i]+int(Re(coher(mqc-i)));
    }

  pos[2*mqc]   = pos[2*mqc-1]+int(Re(coher(2*mqc)));

  for (i=0; i<HS; i++)		// first do the populations
    index[i] = i*(HS+1);

  for (p=0,i=0; i<HS; i++)		// do the coherences
    for (j=0; j<HS; j++,p++)
      if (i!=j)
	{
	  c=int(fabs( 2*(Re(states(i))-Re(states(j)))-0.5) -0.5 );
	  index[pos[c]] = p;
	  pos[c]++;
	}
  delete [] pos;
  return index;
}


int* sort_Op_basis (const spin_sys& sys)


//	Input			sys	: Spin system
//	Output			index	: An integer vector containing
//					  the basis function order from
//					  from highest Fz to lowest
//	Note				: The default basis order is 
//					  0,1,2,...,HS-1 and this is
//					  usually not Fz ordered due to
//					  single spin tensor products used
//					  in formulating spin operators

{
  int HS = sys.HS();			// Get Hilbert space
  col_vector Fztot = sys.qnStates();	// Array of Fz each state
  int *index;
  index = new int[HS];			// Array for indices, fill with
  for (int i=0; i<HS; i++)		// tensor product basis order
    index[i] = i;
  complex max;
  int pos = 1;
  int imax,j;
  while(pos < HS-1)			// Loop over the Hilbert space
    {					// The first and last are always
    max = Fztot(pos);			// ordered properly in all bases
    imax = pos;
    for(j=pos+1; j<HS-1; j++)
      {
      if(Re(max) < Re(Fztot(j)))	// Find largest Fz
        {
        max = Fztot(j);
        imax = j;
        }
      }
    if(imax != pos)			// Reorder elements if needed
      {
      Fztot.put(Fztot(pos), imax);
      Fztot.put(max, pos);
      j = index[imax];
      index[imax] = index[pos];
      index[pos] = j;
      }
    pos++;
    }
  return index;
}


void mqt_v(const spin_sys& sys, gen_op &H, int qn=1, int type=0, int ncols=3)

	// Input	sys	: A spin system
	// 		H	: General operator (associated to sys)
	//		qn	: Transition quantum number
	//			     0 - zero quantum transitions
	//		             1 - single quantum transitions (DEFAULT)
	//			     n - n quantum transitions		
	//		type	: Type of transition output
	//			    <0 - all transitions except qn quantum transitions
	//			     0 - only "qn" quantum transitions (DEFAULT)
	//			     1 - all transitions >= "qn" quantum transitions
	//			    >1 - all transitions <= "qn" quantum transitions
	//		ncols   : Number of columns to output
	// Return	none 	: All indicated transitions associated with
	//			  the input Hamiltonian are sent to standard output
	// Note		     	: Sets H to its EBR

    {
    int size = H.dim();				// Get the operator dimension
    if(!size) return; 				// Exit if NULL operator
    if(type < 0) type = -1;			// Just use type = {-1, 0, 1, 2}
    else if(type > 1) type = 2;
    if(qn < 0) qn *= -1;			// Insure the quantum # is non-negative
    H.set_EBR();				// Set operator in its eigenbasis
    spin_op FZ = Fz(sys);			// Get the Fz spin operator
    gen_op FzOp(FZ);			// Put Fz as a general operator
//    gen_op FzOp(FZ.matrix());			// Put Fz as a general operator
//    gen_op FzOp(FZ.matrix(), sys.get_basis());	// Put Fz as a general operator
    FzOp.Op_base(H);				// Put Fz in eigebasis of H

//			     Begin Printing Setup

// Each printed row of transitions will appear in the following format:

//  [--rst--]([# --> #][sep1][energy]{  Fz})[sep2]([# --> #][sep1][energy]{  Fz})[sep2]... 

//  [--rst--]          (column 1)      [sep2]            (column2)     [sep2]...

// Each column will have length: collen = idig + 5 + idig + sep1 + edig + {2 + idig}
// Each row will have length: rowlen =  rst + ncols*collen + (ncols-1)*sep
// To center each row (for 80 columns) we need rst = 80 - {ncols*collen + (ncols-1)*sep}/2

    int idig = 1;				// Transition index printed size
    while(pow(double(10.0),idig) < size)		// Get the transition index printed size
      idig++;
    std::string Sidig = std::string("%d") + Gdec(idig);
    int tlen = 2*idig+5;			// This length for # --> #
    int sep1 = 1;				// This length for the separator
    if(tlen < 10) sep1 += (10-tlen)/2;		// Insure room for column title
    int edig = 10;				// This length for energies
    int sep2 = 3;				// This length for the separator
    int fzlen = 0;				// This length for Fz
    if(type) fzlen = 2 + idig;
    int collen = tlen + sep1 + edig + fzlen;	// This length for column
    int rst = 80-(ncols*collen+(ncols-1)*sep2);	// This length for the start;
    rst /= 2;

// Above the rows will appear column titles in the following format:

//  [--rstt--](Transition[T2E]Energy[E2F]{Fz}[--sept-])(Transition[T2E]Energy[E2F]{Fz}[--sept-])...
//  [--rst--] ([# --> #][energy]{, Fz})[--sep--]([# --> #] [energy]{, Fz})[--sep--]... 

    std::string T = "Transition";			// Column label for transition, length 10 
    int rstt = rst - (10-tlen)/2;		// Start at rst, adjusted for T
    std::string E = "Energy";			// Column label for energy, length 10;
    int sep1t = (edig - 6)/2 + 1;

//    int fzprn = type;				// Print del Fz in type != 0
    std::string numb[9] = {"Zero","Single","Double",	// Labels for lower quantum orders
                      "Triple", "Four", "Five",
                      "Six", "Seven", "Eight"};
    std::string mqlabel = Gdec(Sidig,qn);		// Label for the used quantum orders
    if((qn) < 9) mqlabel = numb[qn];
    int mqlen = mqlabel.length();
    int spaces = 0;
    std::string blanks = "                                       ";
    switch(type)				// Output an Initial Header
      {
      case(-1):	
        spaces = (80 - (30 + mqlen + 20))/2;
        std::cout << "\n" << std::string(spaces, ' ')
             << "All System Transitions Except "
             << mqlabel << " Quantum Transitions";
        break;
      default:
      case(0):					// Output only qn transitions
        spaces = (80 - (7 + mqlen + 20))/2;
        std::cout << "\n" << std::string(spaces, ' ')
             << "System "
             << mqlabel << " Quantum Transitions";
             break;
      case(1):					// Output transitions >= qn transitions
        if(qn == 0)
          {
          spaces = (80 - 22)/2;
          std::cout << "\n" << std::string(spaces, ' ')
               << "All System Transitions";
          }
        else
          {
          spaces = (80 - (29 + mqlen + 20))/2;
          std::cout << "\n" << std::string(spaces, ' ')
               << "All System Transitions Above "
               << mqlabel << " Quantum Transitions";
          }
        break;
      case(2):	 				// Output transitions <= qn transitions
        if(qn == size)
          {
          spaces = (80 - 22)/2;
          std::cout << "\n" << std::string(spaces, ' ')
               << "All System Transitions";
          }
        else
          {
          spaces = (80 - (29 + mqlen + 20))/2;
          std::cout << "\n" << std::string(spaces, ' ')
               << "All System Transitions Below "
               << mqlabel << " Quantum Transitions";
          }
         break;
       }
    std::string ssep(sep2, ' ');			// Separator between columns
    std::string sst(rst, ' ');			// Spaces before 1st column start
    std::string F = "Fz";
    std::string T2E = std::string(sep1t, ' ');
    std::string E2F = std::string(1+(edig-5)/2, ' ');
    std::string F2T = std::string(1+(10-tlen)/2, ' ');
    std::cout << "\n\n" << std::string(rstt, ' ');	// Begin Column Labels
    for(int n=0; n<ncols; n++)			// Output Each Column Label
      {
      std::cout << T << T2E << E;
      if(type)
        std::cout << E2F << F << F2T;
      }
    std::cout << "\n" << sst;
    double delFz, v;
    std::string levi, levf;				// strings for transition levels 
    int prnt = 1;				// Flag if levition is printed
    int cols = 0;				// Keep track of transitions printed
    for(int i=0; i<size; i++)			// Loop through the initial levels
      {
      levi = std::string(Gdec(Sidig,i)) + " --> ";	// Label for the transition initial level
      for(int j=i+1; j<size; j++)		// Loop through all the final levels
        {
        delFz = Re(FzOp.get(i,i)-FzOp.get(j,j));// This is del Fz of the transition
        if(delFz > -0.1)			// Insure del Fz is positive
          {
          prnt = 0;				// Assume this transition not printed
          switch(type)				// Determine if transition is printed
	    {
	    case(-1):				// Output all but qn transitions
	      if(delFz != qn) prnt = 1; break;
	    default:
            case(0):				// Output only qn transitions
	      if(delFz == qn) prnt = 1; break;
            case(1):				// Output transitions >= qn transitions
	      if(delFz >= qn) prnt = 1; break;
            case(2): 				// Output transitions <= qn transitions
	      if(delFz <= qn) prnt = 1; break;
            }
          if(prnt)				// See if transition is to be printed
            {
            levf = std::string(Gdec(Sidig,j));		//	string for transition final level
            v = Re(H.get(i,i))-Re(H.get(j,j));	// 	This is the transition value
            if(fabs(v) < 1.e-10) v=0;
            std::cout << levi << levf		// 	Print the transition label
                 << std::string(sep1, ' ');
            std::string vform = std::string("%f") + Gdec(edig)
                         + std::string(".2");
            std::cout << Gform(vform,v);
//            std::cout << setw(edig)
//                 << setprecision(2) << v;	//	Print the transition value
            if(fabs(delFz) < 1.e-10) delFz=0;
            if(type) 				//	Print delFz value (if not just 1)
              std::cout << "  " << Gdec(Sidig,int(delFz));
            cols++;				// 	Keep Track of # columns printed
            if(cols >= ncols)			//      If # printed columns of
              { 				// 	transitions is ncol, then
	      std::cout << "\n" << sst;		// 		Begin a new line
              cols = 0;				//		Restart column counter
              }
            else				// 	Else, for a just another value
              std::cout << ssep;			//		output the separator
            }
          }
        }
      }
  std::cout << "\n";
  return;
  }


void wavefunction(const spin_sys& sys, gen_op &Op, int wf, int pbf=0)

	// Input	sys	: A spin system
	// 		Op	: General operator (associated to sys)
	//		wf	: Wavefunction number
	//		pbf	: Flag to specify how to print product basis function
	//			  -1 - designate pdct basis functions by total Fz
	//			   0 - designate pdct basis functions by number (DEFAULT)
	//			   1 - designate pdct basis functions by alpha(a) & beta(b)
// *** pbg=1 only works for systems with all I=1/2.
	// Return	none 	: Indicated wavefunction wf for the working basis
	//			  of operator Op is printed to standard output

  {
  int size = Op.dim();				// Operator dimension size
  if(!size) return; 				// Exit if NULL operator
  matrix B = Op.get_basis().U();		// Retrieve the basis array
  matrix PBF = sys.qStates();			// Retrieve pdct basis functions
  double rcoeff, icoeff;			// These are real & imag coeffs
  int bf = 0;
  int bout = 0;
  int idig = 1;					// Basis fct index printed size
  while(pow(double(10.0),idig) < size)		// Get the index printed size
    idig++;
  std::string ind = Gdec(wf);
  std::string blanks = "                       ";
  std::cout << ind << "." << std::string(idig-ind.length(), ' ');
  for(int i=0; i<size; i++)		// Go through each basis fct for the state
    {
    rcoeff = Re(B.get(i,wf));
    icoeff = Im(B.get(i,wf));
    if(fabs(rcoeff) > 1.e-10)
      {
      if((bout >= 5) && (i < size-1))
        {
        std::cout << "\n     ";
        bout = 0;
        }
      bf++;
      if((rcoeff > 0) && (bf > 1)) std::cout << " + ";
      else if((rcoeff > 0) && (bf == 1)) std::cout << "   ";
      else if((rcoeff < 0) && (bf > 1)) std::cout << " - ";
      else if((rcoeff < 0) && (bf == 1)) std::cout << "  -";
//      cout << setw(5) << setprecision(3) << fabs(rcoeff) << "|";
      std::cout << Gform("%f5.3", fabs(rcoeff)) << "|";
      if(pbf < -1) pbf = -1;
      if(pbf > 1) pbf = 1;
      switch(pbf)
	{
	case(-1):			// Label Basis Functions by Total Fz
          ind = Gdec(i);
          std::cout << ind << ","
               << std::string(idig-ind.length(), ' ');
          if(sys.qnState(i) < 0)
            std::cout << "-";
          else
            std::cout << " ";
          std::cout << fabs(sys.qnState(i));
	  break;
	case(0):			// Label Basis Functions by number
	default:
          std::cout << i;
	  break;
	case(1):			// Label Basis Functions by alpha (a) & beta (b) products
          for(int j=0; j<sys.spins(); j++)
            {
            if(Re(PBF(i,j)) > 0)
              std::cout << "a";
            else if(Re(PBF(i,j)) < 0)
              std::cout << "b";
            }
	  break;
	}
      std::cout << ">";
      bout++;
      }
    }
  std::cout << "\n";
  return;
  }


void wavefunctions(const spin_sys& sys, gen_op &Op, int pbf=0)

	// Input	sys	: A spin system
	// 		Op	: General operator (associated to sys)
	//			  -1 - designate product basis functions by total Fz
	//			   0 - designate product basis functions by number (DEFAULT)
	//			   1 - designate product basis functions by alpha(a) & beta(b)
	// Return	none 	: All wavefunctions of the working basis
	//			  of operator Op are printed to standard output

  {
  int size = Op.dim();			// Operator dimension size
  if(!size) return; 			// Exit if NULL operator
  if(Op.in_EBR())
    std::cout << "\n\tCurrent System Eigenbasis Functions\n\n";
  else
    std::cout << "\n\tCurrent System Basis Functions\n\n";
  for(int i=0; i<size; i++)
    wavefunction(sys, Op, i, pbf);
  std::cout << "\n";
  return;
  }

 void eigensystem(std::ostream& ostr, gen_op Op)

	// Input 	Op	: General operator
	//			  -1 - designate product basis functions by total Fz
	//			   0 - designate product basis functions by number (DEFAULT)
	//			   1 - designate product basis functions by alpha(a) & beta(b)
	// Return	none 	: All wavefunctions of the working basis
	//			  of operator Op are printed to standard output

  {
  Op.set_EBR();				// Place Op into its eigenbasis
  int size = Op.dim();			// Operator dimension size
  if(!size) return; 			// Exit if NULL operator
  matrix eval = Op.get_mx();		// Retrieve matrix: diagonals are eigenvalues
  matrix evec = Op.get_basis().U();	// Retrieve basis: cols are eigenvectors
  int* ri;
  ri = new int[size+1];			// Flag columns as -1=imaginary, 0=complex, 1=real
  int tf = 0;
  int tf1 = 0; 
  int etf = 0;
  int etf1 = 0; 
  complex zco, zen;
  double zre,zim;
  for(int j=0; j<size; j++)		// Loop over all the columns
    {
    if(!etf)
      {
      zen = eval.get(j,j);
      zre = Re(zen);
      zim = Im(zen);
      if(fabs(zre) > 1.e-10)		// Energy has a real component
        {
        if(fabs(zim) > 1.e-10)		// This energy is complex
          {				// so the column is complex
          etf = 1;
          ri[size] = 0;
          }
        else if(etf1 < 0)		// This energy real, previous energy imaginay
          {				// so the column is complex
          etf = 1;
          ri[size] = 0;
          }
        else				// This energy real, so set that it is
          {
          etf1 = 1;
          ri[size] = 1;
          }
        }
      else if(fabs(zim) > 1.e-10)	// Energy is imaginary
        {
        if(etf1 > 1)			// A previous real energy was found
          {				// so the energy is complex
          etf = 1;
          ri[size] = 0;
          }
        else				// Nothing else known
          {				// so flag that an imaginary found
          etf1 = -1;
          ri[size] = -1;
          }
        }
      }
    ri[j] = 1; 				// Set column j print flag to real
    tf = 0;
    tf1 = 0;
    for(int i=0; i<size && !tf; i++)	// Loop over all the rows (coefficients)
      {
      zco = evec.get(i,j);
      zre = Re(zco);
      zim = Im(zco);
      if(fabs(zre) > 1.e-10)		// Element has a real component
        {
        if(fabs(zim) > 1.e-10)		// This element complex
          {				// so the column is complex
          tf = 1;
          ri[j] = 0;
          }
        else if(tf1 < 0)		// This is real, but previous imaginay
          {				// so the column is complex
          tf = 1;
          ri[j] = 0;
          }
        else				// This is real, so set that it is
          {
          tf1 = 1;
          ri[j] = 1;
          }
        }
      else if(fabs(zim) > 1.e-10)	// Element is imaginary
        {
        if(tf1 > 1)			// A previous real was found
          {				// so the column is complex
          tf = 1;
          ri[j] = 0;
          }
        else				// Nothing else known
          {				// so flag that an imaginary found
          tf1 = -1;
          ri[j] = -1;
          }
        }
      }
    }

  complex ev, coeff;
  double rcoeff, icoeff;
  int rf;
  double x;
  std::cout << "\n\t\tCurrent Eigensystem: Eigenvalues and Eigenfunctions\n\n";
  for(int ef=0; ef<size; ef++)		// Loop through all the eigenfunctions
    {
    rf = ri[size];			// Get print flag for eigenvalue column
    ev = eval.get(ef,ef);		// Get eigenvalue for eigenfunction ef 
    if(!rf)
      ostr << ev << " ";		// Output the complex eigenvalue
    else if(rf>0)
      {
      if(Re(ev) > 0)
        ostr << "  ";
      else
        ostr << "- ";
      ostr << fabs(Re(ev)) << "  ";	// Output the real eigenvalue
      }
    else
      {
      if(Im(ev) > 0)
       ostr << "  ";
      else
        ostr << "- ";
      ostr << fabs(Im(ev)) << "i ";	// Output the real eigenvalue
      }
    for(int bf=0; bf<size; bf++)	// Loop through all the basis functions 
      {
      rf = ri[bf];			// Get print flag for basis column
      coeff = evec.get(bf,ef);		// Get bf component of eigenfunction ef
      rcoeff = Re(coeff);		// This is the real component
      icoeff = Im(coeff);		// This is the imaginary component
      if(!rf)
        ostr << coeff;			// Output the complex component
      else if(rf>0) 			// Output the real component
        {
        x = fabs(rcoeff);
        if(x < 1.e-10)
          ostr << " 0    ";
        else if(rcoeff < 0)
          {
          ostr << "-";
//         ostr << setw(5) << setprecision(3) << x;
         ostr << Gform("%f5.3", x);
          }
        else 
          {
          ostr << " ";
 //         ostr << setw(5) << setprecision(3) << x;
          ostr << Gform("%f5.3", x);
         }
        ostr << "  ";
        }
      else				// Output the imaginary component
        {
        x = fabs(icoeff);
        if(x < 1.e-10)
          std::cout << " 0   ";
        else if(icoeff < 0)
          {
          ostr << "-";
//          :ostr << setw(5) << setprecision(3) << x;
          ostr << Gform("%f5.3", x);
         }
        else 
          {
          ostr << " ";
 //         ostr << setw(5) << setprecision(3) << x;
          ostr << Gform("%f5.3", x);
         }
        ostr << "i ";
        }
      }
    ostr << "\n";
    }
  delete [] ri;
  return;
  }


 double vecmax(row_vector &vx)

	// Input	vx       : Data vector
	//		i        : Empty integer
	//              max      : Empty double
	// Return		 : i & max are filled with
	//			   the vector values at it maximum

  {
  int np = vx.elements();
  double max = -HUGE_VAL;		// Find maximum value
  double t;
  for(int i=0; i<np; i++)
    {
    t = Re(vx(i));
    if(t > max)
      max = t;
    }
  return max;
  }


complex integral(const row_vector &vx)

	// Input		vx  : Data vector
	// Return		z   : Integral, sum of all vector elements
	//			      between the limits specified.
	// Note			    : Range includes high & low [low, high]
	// Note			    : It is allowed that low = high
	// Note			    : Move to row_vector!

  {
  complex z = 0.0;
  int nel = vx.elements();
  double re,im; 
  for(int i=0; i<nel; i++)
    {
    re = Re(vx.get(i));
    im = Im(vx.get(i));
    z += complex(fabs(re),fabs(im));
    }
  return z;
  }


 void lwhh(row_vector &vx, int& i1, int& i2)

	// Input	vx       : Data vector
        //		ri	 : Flag for real versus imaginary
	//			   0=reals, non-zero=imaginaries
	// Return		 : Min and Max are returned

  {
  int np = vx.elements();
  double max = -HUGE_VAL;		// First find maximum value
  int maxi = 0;			// and maximum coordinate
  double t;
  int i = 0;
  for(i=0; i<np; i++)
    {
    t = Re(vx(i));
    if(t > max)
      {
      max = t;
      maxi = i;
      }
    }
  double maxhh = max/2.0;	// Set half-height value expected
  double del, delmin = HUGE_VAL;
  i = maxi;
  while(i >= 0) 		// Look for half height to left
    {
    t = Re(vx(i));
    del = fabs(t-maxhh);
    if(del < delmin)
      {
      delmin = del;
      i1 = i;
      }
    i--;
    }     
  i = maxi;
  delmin = HUGE_VAL;
  while(i < np)			// Look for half height to right
    {
    t = Re(vx(i));
    del = fabs(t-maxhh);
    if(del < delmin)
      {
      delmin = del;
      i2 = i;
      }
    i++;
    }     
  }


// ______________________________________________________________________
//                  PULSE SEQUENCE PARAMETER QUERY FUNCTIONS
// ______________________________________________________________________


  int query_isotope(const spin_sys& sys, std::string& Isotype)

	// Input		sys     : Spin system
	//       		Isotype : Isotope type
	// Output		isoset  : Index of a spin having the
	//				  choosen isotope type
	// Note       		Isotype : Isotope type

    {
    int isoset = 0;			// Set default to spin 0 type
    Isotype = sys.symbol(0);		// Set default to spin 0 isotope type
    if(sys.heteronuclear())		// Check if system is heteronuclear
      {
      isoset--;
      while(isoset < 0)
        {
        std::cout << "\n\tWhich Isotope Type? ";
        std::cin >> Isotype;
        for(int k=0; k<sys.spins(); k++)
          if(Isotype == sys.symbol(k))
          {
          isoset = k;
          break;
          }
        if(isoset < 0)
          {
          std::cout << "\n\tSystem Contains No Spin of That Type!\n";
          std::cout << "\n\tChoices are " << sys.symbol(0);
          for(int l=1; l<sys.spins()-1; l++)
            std::cout << ", " << sys.symbol(l);
          std::cout << " and " << sys.symbol(sys.spins()-1) << "\n"; 
          }  
        }
      }
    return isoset;
    }


  int query_isotope(const spin_sys& sys, std::string& Isotype, const std::string& Query)

	// Input		sys     : Spin system
	//       		Isotype : Isotope type
	//       		Query	: Question to be asked
	// Output		isoset  : Index of a spin having the
	//				  choosen isotope type
	// Note       		Isotype : Isotope type

    {
    int isoset = 0;			// Set default to spin 0 type
    Isotype = sys.symbol(0);		// Set default to spin 0 isotope type
    if(sys.heteronuclear())		// Check if system is heteronuclear
      {
      isoset--;
      while(isoset < 0)
        {
        std::cout << Query;
        std::cin >> Isotype;
        for(int k=0; k<sys.spins(); k++)
          if(Isotype == sys.symbol(k))
          {
          isoset = k;
          break;
          }
        if(isoset < 0)
          {
          std::cout << "\n\tSystem Contains No Spin of That Type!\n";
          std::cout << "\n\tChoices are " << sys.symbol(0);
          for(int l=1; l<sys.spins()-1; l++)
            std::cout << ", " << sys.symbol(l);
          std::cout << " and " << sys.symbol(sys.spins()-1) << "\n"; 
          }  
        }
      }
    return isoset;
    }


  int query_isotope(int argc, char* argv[], int argn, const spin_sys& sys, std::string& Isotype)

	// Input		argc    : Number of arguments available
	//       		argv    : Vector of argument values
	// 			argn    : Argument index
	// 			sys     : Spin system
	//       		Isotype : Isotope type
	// Output		isoset  : Index of a spin having the
	//				  choosen isotope type

    {
    int isoset = 0;			// Set default to spin 0 type
    Isotype = sys.symbol(0);		// Set default to spin 0 isotope type
    if(sys.heteronuclear())		// Check if system is heteronuclear
      {
      query_parameter(argc, argv, argn,	// Get filename from command
        "\n\tDesired Isotope Type? ",	// line or ask for them
    	     		      Isotype);
      isoset--;
      while(isoset < 0)
        {
        for(int k=0; k<sys.spins(); k++)
          if(Isotype == sys.symbol(k))
          {
          isoset = k;
          break;
          }
        if(isoset < 0)
          {
          std::cout << "\n\tSystem Contains No Spin of the Specified Type!\n";
          std::cout << "\n\tChoices are " << sys.symbol(0);
          for(int l=1; l<sys.spins()-1; l++)
            std::cout << ", " << sys.symbol(l);
          std::cout << " and " << sys.symbol(sys.spins()-1) << "\n"; 
          std::cout << "\n\tDesired Isotopye Type? ";
          std::cin >> Isotype;
          }  
        }
      }
    return isoset;
    }



double query_offset(int argc, char* argv[], int& argn, spin_system& sys, int isoset)

	// Input		argc    : Number of arguments available
	//       		argv    : Vector of argument values
	// 			argn    : Argument index
	// 			sys     : Spin system
	//       		Isotype : Isotope type
	// Output		isoset  : Index of a spin having the
	//				  choosen isotope type

    {
    double offset = sys.center(isoset);		// Find the center of the system
    std::string typ = "n";				// Flag for type of change
    std::string q = std::string("\n\tTo Center the ")	// string for the question
             + sys.symbol(isoset) 
             + std::string(" Spectrum, An Offset of ")
             + std::string(Gform("%12.4f", offset))
             + std::string(" Hertz is Needed.\n")
             + std::string("\n\tDesired Offset ")
             + std::string("[u=use this value,")
             + std::string(" n=no offset,")
             + std::string(" c=choose an offset]? ");
    query_parameter(argc,argv,argn++,q,typ);	// See what change desired
    if(typ == "u")				// Do this if we're to use
      {						// the calculated offset value 
      sys.offsetShifts(offset, isoset);		// 	Offset the rotating frame
      argn--;					//	Reset query index back
      }
    else if(typ == "n")				// Do this if we don't do any
      {						// offset at all
      offset = 0;				//	No offset
      argn--;					//	Reset query index back
      }
    else if(typ == "c")				// Do this if we want to set a
      {						// specific offset value
      query_parameter(argc,argv,argn,		//	Ask for the value
       "\n\tPlease Input an Offset Value: ",
                                       offset);
      }
    else					// Do this if the first response
      {						// wasn't {u, c, n}
      std::cout << "\n\n\tYour input response to "	//	Indicate a bad response
           << " a query on how to handle the "
	   << " offset was unacceptable.\n\t"
           << "We'll assume the offset is 0!\n";
      offset = 0;				//	Set the offset to 0
      argn--;					//	Reset query index back
      }
    return offset;
    }


double query_offset(spin_system& sys, int isoset, int ask)

	// Input		sys     : Spin system
	//       		isoset  : Index of a spin having the
	//				  choosen isotope type
	//			ask     : Flag to force a query
	//				   0 - ask if non-zero center (Def)
	//				  >0 - always ask
	//				  <0 - never ask
	// Output		offset  : Carrier offset value (Hz) for
	//				  the choosen isotope type

    {
    int repeat = 1;
    double offset = sys.center(isoset);			// Find the center of the system
    std::string typ = "n";
    if(ask>0 || (offset && ask==0))
      {
      std::cout << "\n\tTo Center the " << sys.symbol(isoset) 
           << " Spectrum, An Offset of " << offset
           << " Hertz is Needed.\n";
      while(repeat)
        {
        std::cout << "\n\tDesired Offset "
             << "[u=use this value,"
             << " n=no offset,"
             << " c=choose an offset]? ";
        std::cin >> typ;
        if(typ == std::string("u"))
          {
          sys.offsetShifts(offset, isoset);		// Offset the rotating frame
          repeat = 0;
          }
        else if(typ == std::string("n"))
          {
          offset = 0;
          repeat = 0;
          }
        else if(typ == std::string("c"))
          {
          std::cout << "\n\tPlease Input an Offset Value: ";
          std::cin >> offset;
          std::cout << "\n";
          sys.offsetShifts(offset, isoset);		// Offset the rotating frame
          repeat = 0;
          }
        else
          repeat = 1;
        }
      }
    else if(offset && ask<0)
      sys.offsetShifts(offset, isoset);			// Offset the rotating frame
    return offset;
    }


double query_offset(spin_system& sys, std::string& Isotype, int ask)


	// Input		sys     : Spin system
	//       		Isotype : Isotope type
	//			ask     : Flag to force a query
	//				   0 - ask if non-zero center
	//				  >0 - always ask
	//				  <0 - never ask
	// Output		offset  : Carrier offset value (Hz) for
	//				  the choosen isotope type

    {
    int isoset = -1;
    for(int i=0; i<sys.spins() && isoset<0; i++)
      if(Isotype == sys.symbol(i))
        isoset = i;
    return query_offset(sys, isoset, ask);
    }


double query_Nyquist(const spin_system& sys, int isoset, double lw, double fact)

	// Input		sys     : Spin system
	//       		isoset  : Spin index (spin of choosen isotope type)
	//			lw	: An expected linewidth
	//			fact    : A scaling factor (1.2 = add 20%)
	// Output		Nyqf    : An estimated Nyquist frequency

    {
    int repeat = 1;
    std::string typ("u");
    double NyqF = sys.Nyquist(isoset, 1.2, lw);		// Approximate Nyquist frequency
    std::cout << "\n\tAn Approximate " << sys.symbol(isoset)
         << " Nyquist Frequency is " << NyqF << " Hertz";
    if(lw)
      std::cout << "\n\tAssuming Linewidths of Approximately " << lw << " Hertz";
    while(repeat)
      {
      std::cout << "\n\n\tDesired " << sys.symbol(isoset) << " Nyquist Frequency"
           << "[u=use this value,"
           << " c=choose a frequency]? ";
      std::cin >> typ;
      if(typ == std::string("c"))
        {
        std::cout << "\n\tPlease Input a " << sys.symbol(isoset)
             << " Nyquist Frequency: ";
        std::cin >> NyqF;
        repeat = 0;
        }
      else if(typ == std::string("u"))
        repeat = 0;
      }
    return NyqF;
    //fact = 0.0; // this is not needed, nor is it ever executed.
    }


  double query_Nyquist(const spin_system& sys, std::string& Isotype, double lw, double fact)

	// Input		sys     : Spin system
	//       		Isotype : Isotope type
	//			lw	: An expected linewidth
	//			fact    : A scaling factor (1.2 = add 20%)
	// Output		Nyqf    : An estimated Nyquist frequency

    {
    int isoset = -1;
    for(int i=0; i<sys.spins() && isoset<0; i++)
      if(Isotype == sys.symbol(i))
        isoset = i;
    return query_Nyquist(sys, isoset, lw, fact);
    }


void query_file1D(std::string& filename, int& type)

	// Input	filename: I/O filename
	//		type    : I/O file type
	// Output	none	: Function is void.  The string filename
	//		          is filled with an input name.  The type of
	//			  file is also set to one of the following:
	//			     1 - FrameMaker .mif format
	//			     2 - Felix .dat format
	//			     3 - NMRi format
	//			     4 - MatLab format

  {
  while((type <=0) || (type > 4))
    {
    std::cout << "\n\tPlease Choose a File Format";
    std::cout << "\n\n\t\t1. FrameMaker";
    std::cout << "\n\t\t2. Felix";
    std::cout << "\n\t\t3. NMRi";
    std::cout << "\n\t\t4. MatLab";
    std::cout << "\n\n\tFile Format? ";
    std::cin >> type;
    std::cout << "\n";
    }
  if(type == 2)
    query_FelixFile1D(std::cout, filename);
  else
    {
    std::cout << "\n\tFile Name - Please Include any Extension? ";
    std::cin >> filename;
    }
  return;
  }


void query_FelixFile1D(std::ostream& ostr, std::string& filename)

	// Input	filename: I/O filename
	// Output	none	: Function is void.  It asks for a filename
	//		          assuming the file will be a 1D Felix file
  {
  ostr << "\n\tFelix 1D File Name - Please Include any Extension? ";
  ostr << "\n\tFelix 1D Files Normally Use a .dat Extension";
  ostr << "\n\t(Remember Felix reads only lower case names!)    ";
  std::cin >> filename;
  return;
  }

//void query_output1D(const spin_system& sys, std::string& filename, int& type=-1, int FFT=0)

	//// Input	filename: I/O filename
	//// Input	filename: I/O filename
	//		type//    : I/O file type
	//			//     1 - FrameMaker .mif format
	//			//     2 - Felix .dat format
	//			//     3 - NMRi format
	//			//     4 - MatLab format
	//		sw	:// spectral width (in Hz)
	//		offset//  : Spectral offset (in Hz)
	//		Omega	:
	//		FFT	:// Flag for application of a Fourier Transform
	//// Output	none	: Function is void.  The string filename
	//		//          is filled with an input name.  The type of
	//			//  file is also set to one of the following:

//  {
//  int inPPM = 2;
//  double lf=0, rt=0;
//  switch(type)
//    {
//    case 1:					// FrameMaker output
//    default:
//      if(!FFT)					// If no FFT, output FID
//        {
//        lf = 0.0;
//        rt = data.size() - 1;
//        }
//      else
//        {
//        lf = offset-(sw/2.0);			// Set plot limits in Hertz
//        rt = offset+(sw/2.0);
//        while(inPPM <0 || inPPM >1)		// Querie for Hertz or PPM axis
//          {
//          cout << "\n\n\tAxes in Hz(0) or PPM(1)? ";
//          cin >> inPPM;
//          cout << "\n";
//          }
//        if(inPPM)				// For PPM adjust plot limits
//          {
//          lf = lf/Omega;
//          rt = rt/Omega;
//          }
//        exponential_multiply(data,-8);		// Apodize the FID
//        data = FFT(data);           		// Fourier transform the FID
//        }
//      FM_1D(filename,data,14,5,lf,rt);	 	// Write spectrum to a FrameMaker file
//      break;
//    case 2:					// Felix output
//      Felix(filename,data);	 		// Output FID in Felix format
//      Felix1D_params(cout, Omega, sw, data.size(), offset);
//      break;
//    case 3:
//      NMRi(filename,data);		 	// Output FID in NMRi format
//      break;
//    case 4:
//      MATLAB(filename,"spectrum",data);		// Output FID in MatLab format
//      cout << "\n\t\tMATLAB Spectrum is Internally Named: spectrum";
//      break;
//    cout << "\n";
//    }
//  return;
//  }


void Felix1D_params(std::ostream& ostr, double Omega, double sw, int npts, double offset)

	// Input	ostr    : Output stream
	//			: Spectrometer frequency (MHz)
	//			: Spectral width (Hz)
	//			: Number of points
	//			: Offset frequency (Hz)
	// Output	none	: Function is void.  Felix parameters for 1D
	//		          spectral workup are sent into the output stream
  {
  ostr << "\n\n\t\tParameters Needed for Felix 1D Data Workup\n";
  ostr << "\n\tTo set OMEGA         : sf " << Omega;
  ostr << "\n\tTo set spectral width: sw " << sw;
  ostr << "\n\tTo set offset (no zf): ref " << npts/2 << " " << offset;
  ostr << "\n\tTo set axis in Hertz : ax 2" << "\n\n";
  return;
  }


void Felix2D_params(std::ostream& ostr, double Omega, double sw, int npts, double offset=0, int PPM=3)

	// Input	ostr    : Output stream
	//			: Spectrometer frequency (MHz)
	//			: Spectral width (Hz)
	//			: Number of points
	//			: Offset frequency (Hz)
	// Output	none	: Function is void.  Felix parameters for 2D
	//		          2D homonuclear spectral workup are sent
        //                        into the output stream

  {
  ostr << "\n\n\tParameters Needed for Felix 2D Data Workup\n";
  ostr << "\n\tTo set t2/f2 dimension : rmx 1 " << Omega << " " << sw << " " << PPM << " "
                                                << npts/2 << " " << offset << " F2";
  ostr << "\n\tTo set t1/f1 dimension : rmx 2 " << Omega << " " << sw << " " << PPM << " "
                                                << npts/2 << " " << offset << " F1";
  return;
  }


void Felix2D_params(std::ostream& ostr, double O2, double sw2, int npts2, double off2,
                                   double O1, double sw1, int npts1, double off1, int PPM=3)

	// Input	ostr    : Output stream
	//		O1	: Spectrometer frequency (MHz)
	//			: Spectral width (Hz)
	//			: Number of points
	//			: Offset frequency (Hz)
	// Output	none	: Function is void.  Felix parameters for 2D
	//		          2D homonuclear spectral workup are sent
        //                        into the output stream

  {
  ostr << "\n\n\tParameters Needed for Felix 2D Data Workup\n";
  ostr << "\n\tTo set t2/f2 dimension : rmx 1 " << O2 << " " << sw2 << " " << PPM << " "
                                                << npts2/2 << " " << off2 << " F2";
  ostr << "\n\tTo set t1/f1 dimension : rmx 2 " << O1 << " " << sw1 << " " << PPM << " "
                                                << npts1/2 << " " << off1 << " F1";
  return;
  }

#endif 							// HSauxil.cc
