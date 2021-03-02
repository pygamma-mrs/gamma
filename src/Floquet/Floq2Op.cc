/* Floq2Op.cc *******************************************-*-c++-*-
**                                                              **
**                            G A M M A                         **
**								**
**	Floquet Operator		   Implementation	**
**							 	**
**	Copyright (c) 1992          			 	**
**	Eidgenoessische Technische Hochschule	 		**
**	Labor fuer physikalische Chemie		 		**
**	8092 Zurich / Switzerland		 		**
**						 		**
**      $Header: $
**							 	**
*****************************************************************/

/*****************************************************************
**							 	**
**  Description						 	**
**							 	**
**  The Class FLOQUET OPERATOR Defines the Properties and	**
**  Allowed Operations of a Floquet Operator within the	       	**
**  GAMMA library.						**
**							 	**
*****************************************************************/

#ifndef   Floq2Op_cc_			// Is file already included?
#  define Floq2Op_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Floquet/Floq2Op.h>		// Include the header file
#include <stdlib.h>

// ----------------------------------------------------------------------
// ------------------------ PRIVATE FUNCTIONS ---------------------------
// ----------------------------------------------------------------------

// ______________________________________________________________________
//              CLASS FLOQUET OPERATOR ERROR HANDLING
// ______________________________________________________________________

           

void floq2_op_error (int error)

	// Input		error: Error Flag
	// Output		none : Error Message Output

{
  std::cout << "\nClass Floquet2 Operator: ";
  switch (error)
    {
     case 0:
      std::cout << "Program aborting";
      break;
     case 1:
      std::cout << "Error during construction";
      break;
     case 2:
      std::cout << "Error during floq_op-floq_op Operation";
      break;
     case 3:
      std::cout << "Error during floq_op-matrix Operation";
      break;
     case 4:
      std::cout << "Error accessing Floquet Submatrix";
      break;
     case 5:
      std::cout << "Error accessing internal Floquet component";
      break;
     case 10:
      std::cout << "Inconsistent Floquet-space dimensions";
      break;
     case 20:
      std::cout << "Inconsistent Dimensions during addition";
      break;
     case 21:
      std::cout << "Inconsistent Dimensions during subraction";
      break;
     case 22:
      std::cout << "Inconsistent Dimensions during multiplication";
      break;
     case 23:
      std::cout << "Inconsistent Dimensions during copying";
      break;
     case 30:
      std::cout << "Unable to perform matrix addition";
      break;
     case 31:
      std::cout << "Unable to perform matrix subtraction";
      break;
     case 32:
      std::cout << "Unable to perform matrix multiplication";
      break; 
     case 40:
      std::cout << "Element access beyond floq_op range";
      break;
     case 41:
      std::cout << "Specified position exceeds Floquet dimension";
      break;
     case 42:
      std::cout << "Inconsistent Dimension of gen_op";
      break;
     case 43:
      std::cout << "Sidediagonal number inconsistent with Floquet Dimension";
      break;
     case 50:
      std::cout << "Element access beyond floq_op range";
      break;
     default:
      std::cout << "Unknown error (Number "<< error << ")";
     
  }
  std::cout << ".\n";
}

void volatile floq2_op_fatality (int error)
  {
  floq2_op_error (error);
  if(error) floq2_op_error(0);
  exit(-1);
  } 
                    


// ______________________________________________________________________       
//               INTERNAL FLOQUET OPERATOR MANIPULATIONS
// ______________________________________________________________________        


void floq2_op::add_omegas()

        // Input                     : Floquet Operator (this)
        // Output            F_op    : Floquet operator with omegas added
        //                           : on main diagonal
            
{ 
  int fsize=(2*N1+1)*(2*N2+1)*hs;
  basis bs (matrix(fsize,fsize,i_matrix_type));
  gen_op::Op_base(bs);
  
  for (int n1=0;n1<=2*N1;n1++)

    {
       for (int n2=0;n2<=2*N2;n2++)
         {
           int actpos = (n1*(2*N2+1)+n2)*hs;
           double omegaeff = (n1-N1)*_omega1+(n2-N2)*_omega2;

            for (int h=0;h<hs;h++)
      gen_op::put( gen_op::get(actpos+h,actpos+h) + omegaeff, actpos+h, actpos+h);
         }
    }
} 



void floq2_op::sub_omegas ()

        // Input                     : Floquet Operator (this)
        // Output            F_op    : Floquet operator with omegas subtracted
        //                           : on main diagonal   

{ 
  int fsize=(2*N1+1)*(2*N2+1)*hs;
  basis bs (matrix(fsize,fsize,i_matrix_type));
  gen_op::Op_base(bs);
  
  for (int n1=0;n1<=2*N1;n1++)

    {
       for (int n2=0;n2<=2*N2;n2++)
         {
           int actpos = (n1*(2*N2+1)+n2)*hs;
           double omegaeff = (n1-N1)*_omega1+(n2-N2)*_omega2;

            for (int h=0;h<hs;h++)
      gen_op::put( gen_op::get(actpos+h,actpos+h) - omegaeff, actpos+h, actpos+h);
         }
    }
}


// ----------------------------------------------------------------------
// ------------------------- PUBLIC FUNCTIONS ---------------------------
// ----------------------------------------------------------------------

// ______________________________________________________________________
//               CLASS FLOQUET OPERATOR CONSTRUCTORS/DETRUCTOR
// ______________________________________________________________________


MSVCDLC floq2_op::floq2_op( )
                : gen_op(),_omega1(0),_omega2(0),N1(0),N2(0),hs(0) {};

floq2_op::floq2_op(const floq2_op& FOp)
  {
  N1      = FOp.N1;
  N2      = FOp.N2;
  hs      = FOp.hs;
  _omega1 = FOp._omega1;
  _omega2 = FOp._omega2;
  gen_op::operator=((const gen_op&) FOp);
  }  

floq2_op::floq2_op (int N1_,int N2_,int hs_,double omega1_,double omega2_,matrix &mx) 
         :gen_op(mx),_omega1(omega1_),_omega2(omega2_),N1(N1_),N2(N2_),hs(hs_) 
       
        // Input                N   : 2*N+1: Dimension of Floquet space
        //                      hs  : 2*hs+1: Dimension of spin space
        //                   omega  : Rotation frequency
        // Output                   : returns a new floq_op with Dimension
        //                            2*N+1*hs and k*wr(k=-N...N) on main
        //                            diagonal.

{
 int floqsize=(2*N1_+1)*(2*N2_+1)*hs_;            // make sure matrix has correct size
 if ((mx.cols()!=floqsize) || (mx.rows()!=floqsize))
   {
     floq2_op_error(1);               // Error during construction
     floq2_op_fatality(10);           // Inconsistent Dimensions
   }
}

	// Input                N_   : Floquet space dimension
        //                      hs_  : Hilbert space dimension
        //                   omega_  : Rotation frequency
        // 	                mx   : Matrix
        //                      bs   : Basis
        // Output		Op   : Floquet operator (this) constructed
	//			       with matrix mx and basis bs
        // Note			     : Assumes mx DBR of Op

floq2_op::floq2_op(int N1_,int N2_,int hs_,double omega1_,double omega2_, matrix& mx, basis& bs)
         :gen_op(mx,bs),_omega1(omega1_),_omega2(omega2_),N1(N1_),N2(N2_),hs(hs_)
  {
  int floqsize=(2*N1_+1)*(2*N2_+1)*hs;	// This is the Floquet dimension
  if((mx.cols()!=floqsize) 		// Insure matrix has correct size
  || (mx.rows()!=floqsize))
     {
     floq2_op_error(1);			// Error during construction
     floq2_op_fatality(10);		// Inconsistent dimensions
     }
  }

floq2_op& floq2_op::operator = (const floq2_op& FOp)
  {
  if(this == &FOp) return *this;		// Do nothing if already equal
  N1      = FOp.N1;				// Equate first  photon space
  N2      = FOp.N2;				// Equate second photon space
  hs      = FOp.hs;				// Equate Floquet spaces
  _omega1 = FOp._omega1;			// Equate first  Fourier frequency
  _omega2 = FOp._omega2;			// Equate second Fourier frequency
  gen_op::operator=((const gen_op&)FOp);	// Equate operators
  return *this;
  }

floq2_op::~floq2_op () { }

// ______________________________________________________________________
//  FLOQUET OPERATOR FUNCTIONS, FLOQUET OPERATOR WITH FLOQUET OPERATOR 
// ______________________________________________________________________


floq2_op operator + (floq2_op &Op1, floq2_op &Op2)

	
        // Input		Op1  : Floquet operator.
	// 			Op2  : Floquet operator 
	// Return		Op   : Operator which is the sum
        //	                       of Op1 and Op2,
	//			       Op = Op1 + Op2.
	

  { 
  floq2_op s(Op1.N1, Op1.N2,Op1.hs, Op1._omega1,Op1._omega2);
 if ((Op1.N1 == Op2.N1) && (Op1.hs == Op2.hs) && (Op1._omega1 == Op2._omega1) &&
     (Op1.N2 == Op2.N2) && (Op1._omega1==Op2._omega2))
    s.gen_op::operator=( (gen_op&)Op1 + (gen_op&) Op2);
  else
    { 
       floq2_op_error(2);       // Error during F.Op.-F.Op. -Operation
       floq2_op_fatality(20);   // F.Op.-F.Op. addition error
     }
  return s;
  }

	// Input		Op1  : Floquet operator.
	// 			Op   : Floquet operator (this).
	// Return		Op   : Operator which is the input
        //	                       operator with Op1 added,
	//			       Op = Op + Op1.
	
void floq2_op::operator += (floq2_op &Op1)
{
 if ((Op1.N1 == N1) && (Op1.hs == hs) && (Op1._omega1 == _omega1) &&
     (Op1.N2 == N2) && (Op1._omega2 ==_omega2))
  
    gen_op::operator+=((gen_op&) Op1);
  else
    { 
       floq2_op_error(2);       // Error during F.Op.-F.Op. -Operation
       floq2_op_fatality(20);   // F.Op.-F.Op. addition error
    }
  return;
}


floq2_op operator - (floq2_op &Op1, floq2_op &Op2)

         
         // Input                Op1  : Floquet operator.
         //                      Op2  : Floquet operator.
         // Return               Op   : Operator difference of two input
         //                             operators, Op = Op1 - Op2

  {
  floq2_op d(Op1.N1, Op1.N2, Op1.hs, Op1._omega1,Op1._omega2);
  if ((Op1.N1 == Op2.N1) && (Op1.hs == Op2.hs) && (Op1._omega1 == Op2._omega1) &&
     (Op1.N2 == Op2.N2) && (Op1._omega2 ==Op2._omega2))
  
 
    d.gen_op::operator=((gen_op&) Op1 - (gen_op&) Op2);
  else
    { 
       floq2_op_error(2);       // Error during F.Op.-F.Op. -Operation
       floq2_op_fatality(21);   // F.Op.-F.Op. subtraction error
    }
  return d;
}


floq2_op operator - (floq2_op &Op1)

        // Input                Op1  : Floquet operator
        // Return               Op   : Floquet operator which is the
        //                             negated input operator, Op = - Op1
  {
  floq2_op d(Op1.N1, Op1.N2, Op1.hs, Op1._omega1,Op1._omega2);
  d.gen_op::operator=(-(gen_op&) Op1);
  return d;
  }


void floq2_op::operator -= (floq2_op &Op1)

	// Input		Op1  : Floquet operator.
	// 			Op   : Floquet operator (this).
	// Return		Op   : Operator which is the input
        //	                       operator with Op1 added,
	//			       Op = Op + Op1.
	
{
  if ((Op1.N1 == N1) && (Op1.hs == hs) && (Op1._omega1 == _omega1) &&
     (Op1.N2 == N2) && (Op1._omega2 == _omega2))
  
    gen_op::operator-=((gen_op&) Op1);
  else
    { 
       floq2_op_error(2);       // Error during F.Op.-F.Op. -Operation
       floq2_op_fatality(21);   // F.Op.-F.Op. subtraction error
     }
  return;

}

floq2_op operator * (floq2_op &Op1, floq2_op &Op2)



	// Input		Op1  : Floquet operator.
	// 			Op2  : Floquet operator.
	// Return		Op   : Operator product of the two input
	//			       operators, Op =  Op1 * Op2.

  {
  floq2_op p(Op1.N1, Op1.N2, Op1.hs, Op1._omega1,Op1._omega2);
  if ((Op1.N1 == Op2.N1) && (Op1.hs == Op2.hs) && (Op1._omega1 == Op2._omega1) &&
     (Op1.N2 == Op2.N2) && (Op1._omega2 == Op2._omega2))

    p.gen_op::operator=((gen_op&) Op1 * (gen_op&) Op2);
  else
    { 
       floq2_op_error(2);       // Error during F.Op.-F.Op. -Operation
       floq2_op_fatality(22);   // F.Op.-F.Op. product error
    }
  return p;
  }



void floq2_op::operator *= (floq2_op &Op1)

	// Input		Op1  : Floquet operator.
	// 			Op   : Floquet operator (this).
	// Return		Op   : Operator which is the input
        //	                       operator multiplied into Op1
	//			       Op = Op * Op1.
	// Note		             : Result EXCLUSIVELY in WBR of Op
	// Note			     : Order matters - Op*=Op1 != Op1*=Op

{
 if ((Op1.N1 == N1) && (Op1.hs == hs) && (Op1._omega1 == _omega1) &&
     (Op1.N2 == N2) && (Op1._omega2 == _omega2))
 
    gen_op::operator*=((gen_op&) Op1);
  else
    { 
       floq2_op_error(2);       // Error during F.Op.-F.Op. -Operation
       floq2_op_fatality(22);   // F.Op.-F.Op. product error
    }
 
}



void floq2_op::operator &= (floq2_op &Op1)

	// Input		Op1  : Floquet operator.
	// 			Op   : Floquet operator (this).
	// Return		Op   : Operator which is the input
        //	                       operator multiplied by Op1
	//			       Op = Op1 * Op.
	// Note		             : Result EXCLUSIVELY in WBR of Op
	// Note			     : Order matters - Op&=Op1 != Op1&=Op
{
  if ((Op1.N1 == N1) && (Op1.hs == hs) && (Op1._omega1 == _omega1) &&
     (Op1.N2 == N2) && (Op1._omega2 == _omega2))
   
    gen_op::operator&=((gen_op&) Op1);
  else
    { 
      floq2_op_error(2);       // Error during F.Op.-F.Op. -Operation
      floq2_op_fatality(22);   // F.Op.-F.Op. product error
    }
}





// ______________________________________________________________________
//       FLOQUET OPERATOR FUNCTIONS, FLOQUET OPERATOR WITH MATRIX
// ______________________________________________________________________

	// Input		Op1  : Floquet operator.
	//			mx   : A matrix
	// Return		Op   : Operator sum of two input
	//			       and matrix, Op =  Op1 + mx.

floq2_op operator + (floq2_op &Op1, matrix &mx)

 {
  floq2_op s(Op1.N1, Op1.N2, Op1.hs, Op1._omega1,Op1._omega2);
  int fsize=(2*Op1.N1+1)*(2*Op1.N2+1)*Op1.hs;
  if ((fsize == mx.rows()) && (fsize == mx.cols()))    
                                                            // check dimension

  s.gen_op::operator=((gen_op&)Op1 + (matrix&) mx);
  else
   { 
       floq2_op_error(3);       // Error during F.Op.-Matrix -Operation
       floq2_op_fatality(30);   // F.Op.-Matrix sum  error
    }
  return s;
  }


floq2_op operator + (matrix &mx, floq2_op &Op1)

	
        // Input		mx   : A matrix
	//			Op1  : Floquet operator.
	// Return		Op   : Operator sum of two input
	//			       and matrix, Op =  mx + Op1.

  { 
  floq2_op s(Op1.N1, Op1.N2, Op1.hs, Op1._omega1,Op1._omega2);
  int fsize=(2*Op1.N1+1)*(2*Op1.N2+1)*Op1.hs;
  if ((fsize == mx.rows()) && (fsize == mx.cols()))    
   
  s.gen_op::operator=((matrix&) mx + (gen_op&) Op1);
  else
   { 
       floq2_op_error(3);       // Error during F.Op.-Matrix -Operation
       floq2_op_fatality(30);   // F.Op.-Matrix sum  error
    }
  return s;
  }

	// Input		mx   : A matrix
	// 			Op   : floq operator (this).
	// Return		Op   : Operator which is the input
        //	                       operator with mx added,
	//			       Op = Op + mx.
	// Note		             : Result EXCLUSIVELY in DBR

void floq2_op::operator += (const matrix& mx)
  { 
  int fsize = (2*N1+1)*(2*N2+1)*hs;
  if ((fsize == mx.rows()) && (fsize == mx.cols()))    
   gen_op::operator+=((const matrix&) mx);
   else
    { 
    floq2_op_error(3);       // Error during F.Op.-Matrix -Operation
    floq2_op_fatality(30);   // F.Op.-Matrix sum  error
    }
  return;
  }


floq2_op operator - (floq2_op &Op1, matrix &mx)

	

	// Input		Op1  : Floquet operator.
	//			mx   : A matrix
	// Return		Op   : Operator difference of input
	//			       operator and matrix,
	//			       Op = Op1 - mx
	// Note		             : Result EXCLUSIVELY in DBR
	
  { 
  floq2_op d(Op1.N1, Op1.N2, Op1.hs, Op1._omega1,Op1._omega2);
  int fsize=(2*Op1.N1+1)*(2*Op1.N2+1)*Op1.hs;
  if ((fsize == mx.rows()) && (fsize == mx.cols()))    
  
 
  d.gen_op::operator=((gen_op&) Op1 - (matrix&) mx);
  else
   { 
       floq2_op_error(3);       // Error during F.Op.-Matrix -Operation
       floq2_op_fatality(31);   // F.Op.-Matrix difference error
    }
  return d;
}


floq2_op operator - (matrix &mx, floq2_op &Op1)
	
	// Input		mx   : A matrix
	// 			Op1  : Floquet operator.
	// Return		Op   : Operator difference of input
	//			       matrix and operator
	//			       Op = mx - Op1
	// Note		             : Result EXCLUSIVELY in DBR
	
  { 
  floq2_op d(Op1.N1, Op1.N2, Op1.hs, Op1._omega1,Op1._omega2);
  int fsize=(2*Op1.N1+1)*(2*Op1.N2+1)*Op1.hs;
  if ((fsize == mx.rows()) && (fsize == mx.cols()))    
 
  
  d.gen_op::operator=((matrix&) mx - (gen_op&) Op1);
  else
   { 
       floq2_op_error(3);       // Error during F.Op.-Matrix -Operation
       floq2_op_fatality(31);   // F.Op.-Matrix difference error
   }
  return d; 
  }
	

	// Input		mx   : A matrix
	// 			Op   : Floquet operator (this).
	// Return		Op   : Operator input minus input matrix
	//			       Op = Op - mx
	// Note		             : Result EXCLUSIVELY in DBR

void floq2_op::operator -= (const matrix& mx)
  {
  int fsize=(2*N1+1)*(2*N2+1)*hs;
  if ((fsize == mx.rows()) && (fsize == mx.cols()))    
  
  gen_op::operator-=((const matrix&) mx);
   else
    { 
    floq2_op_error(3);       // Error during F.Op.-Matrix -Operation
    floq2_op_fatality(31);   // F.Op.-Matrix difference error
    } 

  return;
  }


floq2_op operator * (floq2_op &Op1, matrix &mx)


	// Input		Op1  : Floquet operator.
	//			mx   : A matrix
	// Return		Op   : Operator product of the input
	//			       operator and matrix
	//			       Op =  Op1 * mx.
	// Note			     : Result EXCLUSIVELY in DBR
	// Note		             : Order matters - Op1*mx != mx*Op1

{ 
floq2_op p(Op1.N1, Op1.N2, Op1.hs, Op1._omega1,Op1._omega2);
int fsize=(2*Op1.N1+1)*(2*Op1.N2+1)*Op1.hs;
  if ((fsize == mx.rows()) && (fsize == mx.cols()))    
  
  p.gen_op::operator=((gen_op&) Op1*(matrix&) mx);
  else
    { 
       floq2_op_error(3);       // Error during F.Op.-Matrix -Operation
       floq2_op_fatality(32);   // F.Op.-Matrix product error
    } 
  return p;
}
  


floq2_op operator * (matrix &mx, floq2_op &Op1)



	// Input		mx   : A matrix
	//			Op1  : Floquet operator.
	// Return		Op   : Operator product of the input
	//			       matrix and operator
	//			       Op =  mx * Op1
	// Note			     : Result EXCLUSIVELY in DBR
	// Note		             : Order matters - Op1*mx != mx*Op1

{
floq2_op p(Op1.N1, Op1.N2, Op1.hs, Op1._omega1,Op1._omega2);
  int fsize=(2*Op1.N1+1)*(2*Op1.N2+1)*Op1.hs;
  if ((fsize == mx.rows()) && (fsize == mx.cols()))    
  
  p.gen_op::operator=((matrix&) mx* (gen_op&) Op1);
  else
    {  
       floq2_op_error(3);       // Error during F.Op.-Matrix -Operation
       floq2_op_fatality(32);   // F.Op.-Matrix product error
    } 
  return p;
}


void floq2_op::operator *= (matrix &mx)

	// Input		mx   : A matrix
	// 			Op   : Floquet operator (this).
	// Return		Op   : Operator input times input matrix
	//			       Op = Op * mx
	// Note		             : Result EXCLUSIVELY in DBR

{ int fsize=(2*N1+1)*(2*N2+1)*hs;
  if ((fsize == mx.rows()) && (fsize == mx.cols()))    
  
  gen_op::operator*=((matrix&) mx);
   else
     {  
       floq2_op_error(3);       // Error during F.Op.-Matrix -Operation
       floq2_op_fatality(32);   // F.Op.-Matrix product error
     } 
 return;
 
}
	

void floq2_op::operator &= (matrix &mx)

	// Input		mx   : A matrix
	// 			Op   : Floquet operator (this).
	// Return		Op   : Operator multiplied by input matrix
	//			       Op = mx * Op
	// Note		             : Result EXCLUSIVELY in DBR

{ int fsize=(2*N1+1)*(2*N2+1)*hs;
  if ((fsize == mx.rows()) && (fsize == mx.cols()))    
  
  gen_op::operator&=((matrix&) mx);
   else
     {  
       floq2_op_error(3);       // Error during F.Op.-Matrix -Operation
       floq2_op_fatality(32);   // F.Op.-Matrix product error
     } 
 return;
}



// ______________________________________________________________________
//       FLOQUET  OPERATOR FUNCTIONS, FLOQUET OPERATOR WITH SCALAR
// ______________________________________________________________________

/* These functions handle all dealings between Floquet operators and scalar
   values (int, double, complex).

 Function   Arguments      Return       Function   Arguments      Return
 --------   ---------   ------------    --------   ---------  ----------------
    *        FOp1,z     FOp = z*FOp1       *=       this, z   void, z*this
    *        z,FOp1     FOp = z*FOp1       *=       this, d   void, d*this
    *        FOp1,d     FOp = d*FOp1        /       FOp1, z   FOp = (1/d)*this
    *        d,FOp1     FOp = d*FOp1       /=       this, z   void, (1/d)*this
                                                                               */

floq2_op operator * (floq2_op &Op1, complex& z)
{
floq2_op p(Op1.N1, Op1.N2, Op1.hs, Op1._omega1, Op1._omega2);
   p.gen_op::operator=((gen_op&)Op1 * z);
  return p;
}


floq2_op operator * (complex& z, floq2_op &Op1)
{
floq2_op p(Op1.N1, Op1.N2, Op1.hs, Op1._omega1, Op1._omega2);
  p.gen_op::operator=(z*(gen_op&)Op1);
return p;
}

floq2_op operator * (floq2_op& Op1, double d)
{
floq2_op p(Op1.N1, Op1.N2, Op1.hs, Op1._omega1, Op1._omega2);
p.gen_op::operator=(d*(gen_op&)Op1); 
return p;
}

floq2_op operator * (double d, floq2_op &Op1)

{
floq2_op p(Op1.N1, Op1.N2, Op1.hs, Op1._omega1, Op1._omega2);
p.gen_op::operator=(d*(gen_op&)Op1); 
return p;
}
  
void floq2_op::operator *= (const complex& z) { gen_op::operator*=(z); }
void floq2_op::operator *= (      double   d) { gen_op::operator*=(d); }

floq2_op floq2_op::operator / (const complex& z) const
  {
  floq2_op s(N1, N2, hs, _omega1, _omega2);
  s.gen_op::operator=((const gen_op&) (*this) / z);
  return s;
  }

floq2_op floq2_op::operator / (double d) const
  {
  floq2_op s(N1, N2, hs, _omega1, _omega2);
  s.gen_op::operator=((const gen_op&) (*this) / d);
  return s;
  }

void floq2_op::operator /= (const complex& z) { gen_op::operator/=(z); }
void floq2_op::operator /= (      double   d) { gen_op::operator/=(d); }

// ______________________________________________________________________
//               COMPLEX  FLOQUET OPERATOR FUNCTIONS
// ______________________________________________________________________



  
/* to be done ...
gen_op pho_trace (floq2_op &Op)
return t(Op.hs);

        // Input                Op   : Floquet operator
        // Output          gen_op    : Operator containing trace(s) over    
        //                             Photon space

{
  complex sum;
  for (int s1=0; s1<Op.hs; s1++)
    for (int s2=0; s2<Op.hs; s2++)
      {
        sum = 0;
        for (int pho=-Op.N; pho<= Op.N; pho++)
          sum += Op.get(pho,pho,s1,s2);
        t.put(sum,s1,s2);
      }
}        
*/        

	// Input		Op   : Floquet operator (this)
        // Return		int  : Op floquet space dimension
	// Input		Op   : Floquet operator
        // Return		int  : Op Hilbert space dimension 
        // Input		Op   : Floquet operator
        // Return		int  : Op Photon space dimension
        // Input		Op   : Floquet operator
        // Return		int  : Op Photon space dimension

int floq2_op::dim()       const { return gen_op::dim(); }
int floq2_op::size()      const { return gen_op::dim(); }
int floq2_op::hsdim ()    const { return hs; }
int floq2_op::phodim1()   const { return N1; }
int floq2_op::phodim2()   const { return N2; }
double floq2_op::omega1() const { return _omega1; }
double floq2_op::omega2() const { return _omega2; }


	// Input		Op1   : Floquet operator
        // Return		Op    : exponential of Op1
	//				Op = exp(Op1)
        // Note			      : Computed in EBR of Op1

floq2_op exp(floq2_op &Op1)
  {  
  floq2_op EXpOp1(Op1.N1, Op1.N2, Op1.hs, Op1._omega1, Op1._omega2);
   EXpOp1.gen_op::operator= (exp((gen_op&) Op1));  // calculate exp(Op1)
  return EXpOp1;
  }

/*

floq2_op prop (floq2_op &ham , double &time)

	// Input		Op1   : Floquet Hamiltonian
        //                      time  : Evolution time
        // Return		FP    : Floquet Propagator in DBR
	
{     
  ham.set_EBR();
  return exp( complex(0,-2*PI)*ham*time );
}


floq_op fprop (floq_op &FLOQHAM , double &time)
 return FP(FLOQHAM.N, FLOQHAM.hs, FLOQHAM._omega);

	// Input		Op1   : Floquet Hamiltonian
        //                      time  : Evolution time
        // Return		FP    : Floquet Propagator in DBR
	
{     FP=FLOQHAM;

      FP.add_omega ();           // calculate 'full' Floquet Hamiltonian
      //std::cout << "nach omega add.="<<FP;
     
      FP.set_EBR();
      //std::cout <<FP;

      FP=complex(0,-2*PI)*FP*time;  // i * Hf * time
      //std::cout << " nach i*Hf*time=" << FP;

      FP=exp(FP);        
      //std::cout << " nach exp... " <<FP;      
      
      FP.set_DBR();
      //std::cout << " nach set_dbr " << FP;
}
*/
// ______________________________________________________________________
//             FLOQUET OPERATOR COMPONENT MANIPULATIONS
// ______________________________________________________________________

   
// ------------------ Floquet Submatrix Manipulations -----------------

/*gen_op  floq_op::operator () ( int N1, int N2)
return Cop(hs);
        
        // Input             floq_op : Floquet operator (this)
        //                    N1     : Row index of photon dimension
        //                    N2     : Column index of photon dimension
        // Output            gen_op  : at position (N1,N2) of Floquet Matrix

{
  if ( (abs(N1)<=N)&&(abs(N2)<=N) )
    {
      
      for (int h1=(N1+N)*hs; h1<(N1+N+1)*hs; h1++)
	{
	  for (int h2=(N2+N)*hs; h2<(N2+N+1)*hs; h2++)
	    Cop.put(gen_op::get(h1,h2), h1-(N1+N)*hs, h2-(N2+N)*hs);
	}
    }
  else
    {  
      floq_op_error(4);       // Error accessing Floquet submatrix           
      floq_op_fatality(40);   // Element access beyond Floq_operator range
    } 
}

gen_op floq_op::get_block (int N1, int N2)
return Cop(hs);

        // Input             floq_op : Floquet operator (this)
        //                    N1     : Row index of photon dimension
        //                    N2     : Column index of photon dimension
        // Output            gen_op  : at position (N1,N2) of Floquet Matrix

{
   if ((abs(N1)<=N)&&(abs(N2)<=N))
     {
       
       for (int h1=(N1+N)*hs; h1<(N1+N+1)*hs; h1++)
	 {
	   for (int h2=(N2+N)*hs; h2<(N2+N+1)*hs; h2++)
	     
	     Cop.put(gen_op::get(h1,h2), h1-(N1+N)*hs, h2-(N2+N)*hs);
	 }
     }
   else
     {  
       floq_op_error(4);       // Error accessing Floquet submatrix   
       floq_op_fatality(40);      // Element access beyond Floq_operator range
     } 
}

*/

void floq2_op::put_block   (gen_op &Op1, int N1x,int N1y, int N2x, int N2y)
       
        // Input                 Op1 : General operator
        //                    N1,N2  : Position <N2|Op|N1> in Floquet Matrix 
        // Output            floq_op : Floquet operator with gen_op Op1
        //                           : set at position <N2|Op1|N1>

{
  if ((abs(N1x*N1y) <= N1*N1) && (abs(N2x*N2y)<=N2*N2) && (Op1.dim()==hs))
                               
                                // Input Dimensions consistent with Dimensions
                                // of Floq_op (this) ?       
       
    {
	// set basis of this to tensorproduct e * Op1.get_bs
        { int psize=(2*N1+1)*(2*N2+1);
          matrix e   = matrix(psize,psize,i_matrix_type);
          basis  obs = Op1.get_basis();
          matrix m   = tensor_product(e,obs.U());
          basis  fbs(m);
          Op_base(fbs);
        }

	matrix mx  = get_mx();         // Retrieve Floq.- operator matrix
	matrix omx = Op1.get_mx();     // Retrieve gen_op matrix
        
	int xpos=((N1x+N1)*(2*N2+1)+(N2x+N2))*hs;
        int ypos=((N1y+N1)*(2*N2+1)+(N2y+N2))*hs;       
        mx.put_block(xpos,ypos,omx); // put block on desired position

	put_mx(mx);
      }
	
    
   else
      {  
       if (Op1.dim()==hs)
        {
       floq2_op_error(4);           // Error accessing Floquet submatrix   
       floq2_op_fatality(41);       // desired position exceeds Floquet dimension
        }
       else
        {
       floq2_op_error(4);            // Error accessing Floquet submatrix 
       floq2_op_fatality(42);        // Incorrect Dimension of gen_op Op1
        }
      } 
 }

void floq2_op::put_sdiag (gen_op &Op1, int sdn1, int sdn2)
 
        // Input                 Op1 : General operator
        //                       sdn : Number of side-diagonal to be set
        // Output            floq_op : Floquet operator with sidediagonal set

 {     
    if ((abs(sdn1) <= 2*N1) && (abs(sdn2) <=N2) && (Op1.dim()==hs))
            
                        // number of sidediag in dimension of F.-Op.?
      {
	// set basis of this to tensorproduct e * Op1.get_bs
	{ 
          int psize=(2*N1+1)*(2*N2+1);
	  matrix e   = matrix(psize,psize,i_matrix_type);
	  basis  obs = Op1.get_basis();
	  matrix m   = tensor_product(e,obs.U());
	  basis  fbs(m);
	  Op_base(fbs);
	}

	matrix mx  = get_mx();         // Retrieve Floq.- operator matrix
        matrix omx = Op1.get_mx();     // Retrieve gen_op matrix
	int actposx,actposy;
        int lowposb,highposb;
        int lowposa  = sdn1;
        int highposa = 2*N1;
        int d1,d2,d3;
        if (sdn1<0)                     // Diag. on left side of main diag.
          {
          lowposa  = 0;
          highposa = 2*N1+sdn1;
          }
       for (d1=lowposa; d1 <= highposa; d1++)
          {
             actposx=d1*(2*N2+1);
             actposy=(d1-sdn1)*(2*N2+1);

             lowposb=sdn2;
             highposb=2*N2;
             if (sdn2<0)                     // Diag. on left side of main diag.
              {
               lowposb  = 0;
               highposb = 2*N2+sdn2;
              }
             for (d2=actposx+lowposb; d2 <= actposx+highposb; d2++)
               {
                  d3=d2-sdn2-d1*(2*N2+1)+actposy;
                  mx.put_block(d3*hs,d2*hs,omx);
               }
	   }
      put_mx(mx);
      }
	
    
   else
     { 
       floq2_op_error(4);             // Error accessing Floquet submatrix
       floq2_op_fatality(43);         // Sidediagonal number inconsistent with
                                     // Floquet Dimension
     }
 
}



  


// ----------------- Individual Element Manipulations -------------------



/*void floq_op::put (complex &z, int N1, int N2, int H1, int H2)
     
       // Input                      : Floq_operator (this)
       //                   N1,N2    : Indices of photon space:-N...N
       //                   H1,H2    : Indices of Hilbert space:0...hs
       // Output            none     : Value of <N2,H2|Op|H1,N1> will be
       //                            : set to z
       // Note                       : floq_op will be changed



  {
   if ((abs(N1)<=N)&&(abs(N2)<=N)&&(abs(H1)<=hs)&&(abs(H2)<=(H2)))
   
     {
       int actposa=(N+N1)*hs+H1;  // matrix starts with (0,0) ...
       int actposb=(N+N2)*hs+H2;
       gen_op::put(z,actposa,actposb);
     }
  else
     {  
       floq_op_error(5);        // Error accessing internal component
       floq_op_fatality(50);    // Element access beyond Floq_operator range
     } 
  }


complex floq_op::get (int N1, int N2, int H1, int H2 )

       // Input                      : Floq_operator (this)
       //                   N1,N2    : Indices of photon space:-N...N
       //                   H1,H2    : Indices of Hilbert space:0...hs
       // Output                z    : complex value of <N2,H2|Op|H1,N1>
       // Note                       : floq_op remains unchanged


{
  if ((abs(N1)<=N)&&(abs(N2)<=N)&&(abs(H1)<=hs)&&(abs(H2)<=(H2)))
   
     {
       int actposa=(N+N1)*hs+H1;  // matrix starts with (0,0) ...
       int actposb=(N+N2)*hs+H2;
       //std::cout << "actposa="<<actposa<<"\n";
       //std::cout << "actposb="<<actposb<<"\n";
       return gen_op::get(actposa,actposb);
     }
  else
     {  
       floq_op_error(5);        // Error accessing internal component
       floq_op_fatality(50);    // Element access beyond Floq_operator range
     } 
 }


void floq2_op::put (complex &z, int row, int col )
       
       // Input                      : Floq_operator (this)
       //                 row,col    : Indices 
       // Output            none     : Value of <row|Op|col> in WBR
       //                            : set to z
       // Note                       : floq_op will be changed
{
  gen_op::put(z,row,col);
}
 
complex floq2_op::get ( int row, int col)


       // Input                      : Floq_operator (this)
       //                 row,col    : Indices 
       // Output            none     : Value of <row|Op|col> in WBR
       // Note                       : floq_op remains unchanged
{
  return gen_op::get(row,col);
}
*/

// ______________________________________________________________________
//                  OPERATOR REPRESENTATION MANIPULATIONS
// ______________________________________________________________________

void floq2_op::set_DBR () { gen_op::set_DBR(); }
void floq2_op::set_EBR () { gen_op::set_EBR(); }

// ______________________________________________________________________
//                   CLASS FLOQUET OPERATOR I/O FUNCTION
// ______________________________________________________________________

	// Input		out  : Output stream
	// 			Op   : A Floquet operator
	// Return		     : Output stream into which Op is put

std::ostream& operator<< (std::ostream& out, const floq2_op &Op)
  {
  out  << "Photon space dimension1  : " << Op.N1 <<"\n"
       << "Photon space dimension2  : " << Op.N2 <<"\n"
       << "Hilbert space dimension : "  << Op.hs <<"\n"
       << "Omega1                   : " << Op._omega1<< "\n"
       << "Omega2                   : " << Op._omega2<< "\n";
  return out<< (const gen_op&) Op;
  }



#endif							// Floq2Op.h
 



