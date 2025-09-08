/* SuperOp.cc ***************************************************-*-c++-*-
**									**
** 	                          G A M M A				**
**									**
**	Superoperators                            Implementation	**
**								 	**
**	Copyright (c) 1991, 1992, 1993				 	**
**	Scott Smith				 			**
**	Eidgenoessische Technische Hochschule	 			**
**	Labor fuer physikalische Chemie		 			**
**	8092 Zurich / Switzerland		 			**
**						 			**
**      $Header: $
**						 			**
*************************************************************************/

/*************************************************************************
**							 		**
** Description							 	**
**								 	**
** The Class SUPER OPERATOR (super_op) Defines the Properties & Allowed	**
** Operations of a Quantum Mechanical Superoperator in Liouville Space 	**
** for the GAMMA platform.						**
**								 	**
*************************************************************************/

#ifndef _super_op_cc_			// Is file already included?
#define _super_op_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// this is the implementation
#endif

#include <LSLib/SuperOp.h>			// Include the interaface
#include <Basics/Gconstants.h>			// Include PI definition
#include <Basics/Gutils.h>			// Include GAMMA errors
#include <Basics/StringCut.h>


//*** :FIXME: Used to round double precision square roots to integers.
// e.g. static_cast<int>( sqrt(double(9)) + SUPEROP_EPSILON ) = 3
// With more understanding of how this is used, 
// this could be changed or removed.
const double SUPEROP_EPSILON = 0.000001;

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
	
// ____________________________________________________________________________
// i                   CLASS SUPER OPERATOR ERROR HANDLING
// ____________________________________________________________________________

        // Input                LOp     : Superoperator (this)
        //                      eidx    : Error index
        //                      noret   : Flag for linefeed (0=linefeed)
        // Output               void    : An error message is output

/* The following error messages use the defaults set in the Gutils package

                Case                          Error Message                

                (0)                     Program Aborting.....
                (9)                     Problems During Construction
	        (10)			Cannot Access Internal Component
                default                 Unknown Error                        */
 
void super_op::LOperror(int eidx, int noret) const
  {
  std::string hdr("SuperOperator");
  std::string msg;
  switch (eidx)
    {
    case 7:  GAMMAerror(hdr,"Accessing Null Superoperator",noret);break;// (7)
    case 20: GAMMAerror(hdr,"Unable To Perform Addition",noret);  break;// (20)
    case 21: GAMMAerror(hdr,"Cannot Perform Subtraction",noret);  break;// (21)
    case 22: GAMMAerror(hdr,"Unable To Do Multiplication",noret); break;// (22)
    case 39: GAMMAerror(hdr,"Bad Mix With Superoperator",noret);  break;// (39)
    case 40: GAMMAerror(hdr,"Bad Mix With Operator",noret);       break;// (40)
    case 41: GAMMAerror(hdr,"Trouble Mixing With Matrix",noret);  break;// (41)
    case 49: GAMMAerror(hdr,"No Associated Hilbert Space",noret); break;// (49)
    case 50: GAMMAerror(hdr,"Rectangular Array Use",noret);       break;// (50)
    case 51: GAMMAerror(hdr,"Matrix-Basis Mismatch",noret);       break;// (51)
    case 52: GAMMAerror(hdr,"LOp-LOp Dimensions Mismatch",noret); break;// (52)
    case 53: GAMMAerror(hdr,"LOp-Op Dimensions Mismatch",noret);  break;// (53)
    case 56: GAMMAerror(hdr,"Element Access Out Of Range",noret); break;// (56)
    case 67: GAMMAerror(hdr,"Cannot Access Element",noret);       break;// (67)
    case 77: GAMMAerror(hdr,"Problems Setting Basis",noret);      break;// (77)
    case 90: GAMMAerror(hdr,"Mixing LOp-LOp Hilbert Bases",noret);break;// (90)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }    
  }


/* The following error messages use the defaults set in the Gutils package

                Case                          Error Message                

                (1)                     Problems With File PNAME
                (2)                     Cannot Read Parameter PNAME
                (3)                     Invalid Use Of Function PNAME
                default                 Unknown Error - PNAME                */  
    
void super_op::LOperror(int eidx, const std::string& pname, int noret) const
  {                                                                             
  std::string hdr("SuperOperator");
  std::string msg;
  switch(eidx)
    {
    case 5: msg = std::string("Bad Use Of ") + pname + std::string(" Function");
             GAMMAerror(hdr,msg,noret);  break;   	//                (5)
    case 6: msg = std::string("Cannot Access Element ") + pname;
             GAMMAerror(hdr,msg,noret);  break;   	//                (6)
    default: GAMMAerror(hdr, eidx, pname, noret); break;// Unknown Error  (-1)
    }
  }


void volatile super_op::LOpfatal(int eidx) const
  {
  LOperror(eidx, 1);
  if(eidx) LOperror(0);
  GAMMAfatal();					// Clean exit from program
  }
 

volatile void super_op::LOpfatal(int eidx, const std::string& pname) const
  {                                                                 
  LOperror(eidx, pname, 1);			// Normal non-fatal error
  if(eidx) LOperror(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

	
// ____________________________________________________________________________
// ii             USEFUL INTERNAL FUNCTIONS FOR SUPER OPERATORS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A            CLASS SUPER OPERATOR CONSTRUCTORS/DESTRUCTOR
// ____________________________________________________________________________

	
super_op::super_op ( ) { HSp=0; LSp=0; }

        // Input		none :
	// Output		LOp  : A NULL super operator (this)


super_op::super_op(const matrix& mx1)

	// Input		mx1  : A matrix
	// Output		LOp  : A super operator (this)
	//			       constructed with matrix mx1
  //			       as default basis representation

  {
 
	if(!checkLOp(mx1, 1))		// Insure matrix acceptable
	{
    LOpfatal(9);		// Error in LOp construction
	}

	// Set internal class variable "mx".
	mx = mx1;			      // Copy input matrix
  LSp = mx.rows();		// Set Liouville space

  // HSp = int(sqrt(LSp));		// Set Hilbert space
	// replace with line below. int(sqrt(int)) fix.
	// **** add a little to make sure we get the right int back.
	// for sqrt(9), int(2.999999999... + .000001) = 3
	HSp = static_cast<int>(sqrt(static_cast<double>(LSp)) + SUPEROP_EPSILON); 
  
	Hbs = basis(HSp);		// Set Hilbert space basis to default
  Lbs = basis(LSp);		// Set Liouville space basis to default

  }


super_op::super_op(const matrix& mx1, const matrix& mx2)

	// Input		mx1  : A matrix in Liouville space
        //                      mx2  : A matrix in Hilbert space
	// Output		LOp  : A superoperator (this)
        //                             with rep mx1 & basis mx2

  {
  if(!mx1.cols()) {HSp=0; LSp=0; return;}	// Check for NULL matrix
  if(!checkLOp(mx1, mx2, 1))			// Insure acceptable pair
    LOpfatal(9);				// Error during construction
  mx = mx1;					// Copy input matrix
  Hbs = basis(mx2);				// Create HS basis from matrix
  LSp = mx1.rows();				// Set Liouville space
  HSp = mx2.rows();				// Set Hilbert space
  Lbs = basis(LSp);				// Set LS basis to default
  }
        

super_op::super_op(const matrix& mx1, const basis& bs1)

	// Input		mx1  : A matrix (in L. space)
        //                      bs1  : A basis (in H. space)
	// Output		LOp  : A super operator (this)
        //                             with rep mx1 & basis bs1

  {
  if(!checkLOp(mx1, bs1, 1))		// Insure acceptable pair
    LOpfatal(9);			// Error during construction
  if(!mx1.cols()){ HSp=0; LSp=0; return;}// NULL LOp if NULL matrix
  mx = mx1;				// Copy input Liouville matrix
  Hbs = bs1;				// Copy input Hilbert space basis
  LSp = mx.rows();			// Set Liouville space
  HSp = bs1.dim();			// Set Hilbert space
  Lbs = defLSbasis(bs1);		// Set Liouville basis to default
  }


super_op::super_op(const std::vector<matrix>& mxc, const std::vector<matrix>& bsc)
  {
  LSp = 0;				// Start with no Liouville space
  HSp = 0;				// Start with no Hilbert space
  int nc = mxc.size();			// # of components in composite
  int* ncd;				// space.  Array ncd will hold
  ncd = new int[nc];			// sub-space Liouville dimensions
  int i=0;
  for(i=0; i<nc; i++)			// Loop over components and 
    {					// figure out what the total
    LSp += mxc[i].rows();		// Liouville & Hilbert dimensions
    if(!bsc.size())			// Default basis if no bsc set
      {
      ncd[i] = bsc[i].rows();		//   Hilbert dimension this comp.
      HSp += ncd[i];			//   Add to total Hilbert dimenson
      basis hbs(bsc[i]);		//   Construct Hilbert space bassi
      if(!checkLOp(mxc[i],hbs,1))	//   Insure {mx,bs} pairing is OK
        {
        LOpfatal(9);			// Construction trouble
        }
      }
    else				// If no vector of basis arrays
      {					// were input, set each basis to

      //ncd[i] = int(sqrt(mxc[i].rows()));// a default Hilbert space basis
			// *** int(sqrt(int)) fix.
			ncd[i] = static_cast<int>(sqrt(static_cast<double>(mxc[i].rows())) + SUPEROP_EPSILON);
      
			HSp += ncd[i];			// and store it in basis vector
      }  
    }
  matrix MX(LSp,LSp,i_matrix_type);	// Construct LOp Liouville array
  matrix bs(HSp,HSp,i_matrix_type);	// Construct LOp Hilbert basis
  int pinblockL = 0, pinblockH = 0;	// Now we must fill up the {mx,bs}
  for(i=0; i<nc; i++)			// Loop over the components and
    {
    MX.put_block(pinblockL,pinblockL,mxc[i]);
    pinblockL += mxc[i].rows();		// Next Liouville component position
    if(bsc.size())			// If basis isn't NULL then
      {					// explicity set up HS basis array
      bs.put_block(pinblockH,pinblockH,bsc[i]);
      pinblockH += bsc[i].rows();	// Next Hilbert component position
      }
    }
  basis BS(bs, nc, ncd);		// Construct Hilbert basis
  mx = MX;				// Set LOp Liouville array
  Hbs = BS;				// Set LOp Hilbert basis
  Lbs = basis(LSp);			// LS basis (no component info)
  }


super_op::super_op(matrix* mxc, int nc, matrix* bsc)

        // Input                nc  : Number of components in multi-sys
        //                      mxc : SOp submatrices associated with comps
        //                      bsc : HS basis matrices associated with comps
        // Output               SOp : Super_op (this) containing
        //                            mxc as diagonal blocks of mx, and bsc
        //                            as diagonal blocks of Hbs

  {
  LSp = 0;
  HSp = 0;
  int* ncd;
  ncd = new int[nc];
  int i=0;
  for(i=0; i<nc; i++)
    {
    LSp += mxc[i].rows();
    if(bsc != NULL)                        // bsc!=NULL explicit basis 
      {
      ncd[i] = bsc[i].rows();
      HSp += ncd[i];
      basis hbs(bsc[i]);
      if( !checkLOp(mxc[i], hbs, 1)  )        // bsc and mxc at odds
        {
        LOpfatal(9);			// Construction trouble
        }
      }
    else
      {					  // this is default HS basis
      
			//ncd[i] = int( sqrt(mxc[i].rows()) );
			// *** int(sqrt(int)) fix.
			ncd[i] = static_cast<int>(sqrt(static_cast<double>(mxc[i].rows())) + SUPEROP_EPSILON);

      HSp += ncd[i];
      }  
    }
  matrix MX(LSp, LSp, i_matrix_type), bs(HSp, HSp, i_matrix_type);
  int pinblockL = 0, pinblockH = 0;
  for(i=0; i<nc; i++)
    {
    MX.put_block(pinblockL, pinblockL, mxc[i]);
    pinblockL += mxc[i].rows();
    if( bsc != NULL )
      {
      bs.put_block(pinblockH, pinblockH, bsc[i]);
      pinblockH += bsc[i].rows();
      }
    }
  basis BS(bs, nc, ncd);
  mx = MX;
  Hbs = BS;
  Lbs = basis(LSp);	     	      // LS basis doesn't contain info on comps
  return;
  }

super_op::super_op(const super_op& LOp1)
  {
  LSp = LOp1.LSp;		// Set Liouville space of LOp
  HSp = LOp1.HSp;		// Set Hilbert space of LOp
  if(LSp)			// Check for NULL LOp1
    {
    mx  = LOp1.mx;		// Copy input super operator matrix
    Hbs = LOp1.Hbs;		// Copy input super operator basis
    Lbs = LOp1.Lbs;		// Copy input super operator Liouville basis
    }
  }
        

super_op::super_op(const gen_op& Op1, const gen_op& Op2)

	// Input		Op1  : A general operator
	//			Op2  : A general operator
	// Output		LOp  : A superoperator (this) as given by
        //				     LOp = |Op1><Op2|
	// Note			     : Performed in the basis of Op1

  {
  HSp = Op1.dim();			// Set Hilbert space of LOp
  LSp = HSp*HSp;			// Set Liouville space of LOp
  matrix mx1, mx2;
  if(HSp && Op2.dim())			// Check for NULL Operators
    {
    Op2.Op_base(Op1);			// Put Op2 into the basis of Op1
    Hbs = Op1.get_basis();		// LOp HS basis equal to Op basis
    Lbs = basis(LSp);			// Set LS basis to default
    mx1=(Op1.get_mx()).resize(LSp, 1);	// Op1 matrix as a column vector
    mx2=(Op2.get_mx()).resize(1, LSp);	// Op2 matrix as a row
    mx = mx1*mx2;			// Compute super operator matrix
    }
  }

super_op::~super_op () { LSp = 0; HSp = 0; }

// ____________________________________________________________________________
// B      SUPER OPERATOR FUNCTIONS: SUPER OPERATOR WITH SUPER OPERATOR 
// ____________________________________________________________________________
                                                                                
/* These functions allow for simple mathematical operations between two super-
   operators.  This includes addition, subtraction, multiplication. There is
   one unary function as well, negation.

   Operator Arguments      Result        Operator Arguments    Result
   -------- --------- -----------------  -------- --------- -------------------
      -        LOp    Returns -LOp          -=    LOp,LOp1  LOp1 subt. from LOp
      +     LOp,LOp1  Returns LOp+LOp1       *    LOp1,LOp2 Returns LOp1*LOp2
     +=     LOp,LOp1  LOp1 added to LOp     *=    LOp,LOp1  LOp mult into LOp1
      -     LOp1,LOp2 Returns LOp1-LOp2                                      */

super_op super_op::operator + (const super_op& LOp1) const
  { super_op LOp(*this); LOp += LOp1; return LOp; }

//void super_op::operator += (const super_op& LOp1)
super_op& super_op::operator += (const super_op& LOp1)
  {
  if(!LOp1.LSp) return (*this);			// Just exit if NULL LOp1
  if(!LSp) { *this = LOp1; return (*this); } 	// Use assign if NULL LOp
  if(!checkLOp(LOp1, 1))		// Insure LOp-LOp1 compatibility
    {
    LOperror(39, 1);			// LOp, LOp function error
    LOpfatal(20);			// LOp, LOp addition problems
    }
  LOp_Hbase(LOp1);			// LOp into LOp1 Hilbert basis
  LOp_base(LOp1);			// LOp into LOp1 Liouville basis
  mx += LOp1.mx;			// Add LOp1 matrix to LOp
  return (*this);
  }

super_op super_op::operator - (const super_op& LOp1) const
  { super_op LOp(*this); LOp -= LOp1; return LOp; }

//void super_op::operator -= (const super_op& LOp1)
super_op& super_op::operator -= (const super_op& LOp1)
  {
  if(!LOp1.LSp) return (*this);			// No change if NULL LOp1
  if(!LSp) { *this = -LOp1; return (*this); }	// If NULL LOp, return -LOp1
  if(!checkLOp(LOp1, 1))		// Check LOp-LOp1 compatibility
    {
    LOperror(39, 1);			// LOp, LOp function error
    LOpfatal(21);			// LOp, LOp subtraction problems
    }
  LOp_Hbase(LOp1);			// LOp into LOp1 Hilbert basis
  LOp_base(LOp1);			// LOp into LOp1 Liouville basis
  mx -= LOp1.mx;			// Add LOp1 matrix to LOp
  return (*this);
  }

super_op super_op::operator - () const
  { 
  super_op LOp(mx.matrix::operator-(), Hbs);
  LOp.Lbs = Lbs; 
  return LOp;
  }

super_op super_op::operator * (const super_op& LOp1) const
  { super_op LOp(*this); LOp *= LOp1; return LOp; }

//void super_op::operator *= (const super_op& LOp1)
super_op& super_op::operator *= (const super_op& LOp1)
  {
  if(!LSp)      { *this=super_op(); return (*this); } 	// LOp NULL,  product is NULL
  if(!LOp1.LSp) { *this=super_op(); return (*this); }	// LOp1 NULL, product is NULL
  if(!checkLOp(LOp1, 1))			// Insure HS compatibility
    {						// If incompatible, then quit
    LOperror(39, 1);			// LOp, LOp function error
    LOpfatal(22);			// LOp*LOp problems
    }
  LOp_Hbase(LOp1);				// Put LOp in LOp1 HS basis
  LOp_base(LOp1);				// Put LOp in LOp1 LS basis
  mx *= LOp1.mx;				// Multiply into LOp1 matrix
  return (*this);
  }

super_op& super_op::operator &= (const super_op& LOp1) 
{                                                            
  if(!LSp)           
    *this=super_op();  // LOp is NULL, product NULL
  else if(!LOp1.LSp) 
    *this=super_op();  // LOp1 is NULL, product NULL
  else                                  // Both LOp & LOp1 exist so take
  {                                   // their product
    if(!checkLOp(LOp1, 1)) 		// Insure Hilbert space compatibility
    {                                 // If incompatible, then quit
      LOperror(39, 1);		// LOp, LOp function error
      LOpfatal(22);		// LOp, LOp multiplication problems
    }  
    LOp_Hbase(LOp1);                    // Put LOp in LOp1 Hilbert basis
    LOp_base(LOp1);                     // Put LOp in LOp1 Liouville basis
    mx = LOp1.mx*mx;                    // Multiply by LOp1 matrix
  }
  
  return *this;
}  

void super_op::operator = (const super_op& LOp1)

	// Input		LOp  : A super operator (this)
	// 			LOp1 : A super operator
	// Return		LOp  : Super Operator (this) which
	//			       has been set equivalent to LOp1
        // Note			     : result in the basis of LOp1

  {
  LSp = LOp1.LSp;			// Set Liouville space of LOp
  HSp = LOp1.HSp;			// Set Hilbert space of LOp
  mx  = LOp1.mx;			// Copy LOp1 matrix
  Hbs = LOp1.Hbs;			// Copy LOp1 Hilbert space basis
  Lbs = LOp1.Lbs;			// Copy LOp1 Liouville space basis
  }


// ____________________________________________________________________________
// C       SUPER OPERATOR FUNCTIONS: SUPER OPERATOR WITH OPERATOR 
// ____________________________________________________________________________

/* These functions allow for simple mathematical operations between a super-
   operator in Liouville space and an operator in Hilbert space.  This only
   includes multiplication in a bra/ket fashion currently, 

   Operator Arguments      Result        Operator Arguments       Result
   -------- --------- -----------------  -------- --------- -------------------
      *      LOp, Op      LOp|Op>           *      Op,LOp         <Op|LOp
     
   The result of test multiplcations will be an operator in Hilbert space
   in the (HS) basis of LOp.  Not that this operation must occur when LOp is
   in its default LS basis. If LOp is not in its LS default basis, it will be
   converted according to                      -1
                            LOp   = U * LOp * U
                               DB
          -l
   where U   will be the same as the adjoint of U if U is unitary.           */


gen_op operator * (const super_op& LOp, const gen_op& Op1)
  {
  if(!Op1.dim() || !LOp.LSp) return gen_op();	// NULL Op if Op1 or LOp Null
  if(!LOp.checkLOp(Op1, 1))			// Check LOp - Op compatibility
    {
    LOp.LOperror(40,1);			// LOp, Op function error
    LOp.LOpfatal(22);			// Cannot do the multiplication
    }
  LOp.set_HBR();				// Put LOp in default LS basis
  LOp.LOp_base(Op1);				// Put Op into HS basis of LOp
  col_vector mk = Op1.superket();		// Set Op as a superket
  mk = LOp.mx * mk;				// Do multiplcation with arrays
  gen_op Op2 = Op1;				// This is the Op for return
  Op2.desuperket(mk);				// Make it from the superket
  return Op2;
  }


gen_op operator * (const gen_op& Op1, const super_op& LOp)
  {
  if(!Op1.dim() || !LOp.LSp) return gen_op();	// NULL Op if Op1 or LOp Null
  if(!LOp.checkLOp(Op1, 1))			// Check LOp - Op compatibility
    {
    LOp.LOperror(40,1);			// LOp, Op function error
    LOp.LOpfatal(22);			// Cannot do the multiplication
    }
  LOp.set_HBR();				// Put LOp into default LS basis
  LOp.LOp_base(Op1);				// Put Op into HS basis of LOp
  col_vector mk = Op1.superket();		// Set Op as a superket
  mk = mk * LOp.mx;				// Do multiplication with arrays
  gen_op Op2 = Op1;				// This is the Op for return
  Op2.desuperket(mk);				// Make from the superket
  return Op2;
  }


// ____________________________________________________________________________
// D            SUPER OPERATOR FUNCTIONS, SUPER OPERATOR WITH SCALAR
// ___________________________________________________________________________

/* These functions allow for simple mathematical operations between a super-
   operators in Liouville space and an constant.  Functions exist that provide
   users with the ability to multiply and divide superoperators with constants.

   Operator Arguments      Result        Operator Arguments       Result
   -------- --------- -----------------  -------- --------- -------------------
      *      LOp, z        z*LOp             /      LOp,z       (1/z)*LOp
      *      z, LOp        z*LOp            /=      this,z      (1/z)* this
      *=     this, z       z*this                                            */

super_op operator * (const super_op& LOp1, const complex& z)
  { super_op LOp(LOp1); LOp.mx *= z; return LOp; }

super_op operator * (const complex& z, const super_op& LOp1)
  { super_op LOp(LOp1); LOp.mx *= z; return LOp;}

super_op operator * (const super_op& LOp1, double d)
  { super_op LOp(LOp1); LOp.mx *= d; return LOp; }

super_op operator * (double d, const super_op& LOp1)
  { super_op LOp(LOp1); LOp.mx *= d; return LOp;}

super_op& super_op::operator *= (const complex& z) 
  { 
    if(LSp) mx *= z; 
    return (*this);
  }

super_op& super_op::operator *= (double d)         
  { 
    if(LSp) mx *= d; 
    return (*this);
  }

super_op operator / (const super_op& LOp1, const complex& z)
  { super_op LOp(LOp1); LOp /= z; return LOp; }

super_op operator / (const super_op& LOp1, double d)
  { super_op LOp(LOp1); LOp /= d; return LOp; }

super_op& super_op::operator /= (const complex& z) 
  { 
    if(LSp) mx /= z; 
    return (*this);
  }

super_op& super_op::operator /= (double d)         
  { 
    if(LSp) mx /= d; 
    return (*this);
  }

// ____________________________________________________________________________
// E                     SUPER OPERATOR COMPLEX FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//          Left And Right Transformation (Translation) Superoperators
// ----------------------------------------------------------------------------

/* Left:   LOp*Op1 = Op*Op1           where             LOp = Op (X) E

  					                              t
   Right:  LOp*Op1 = Op1*Op           where             LOp = E (X) Op       */

super_op left(const gen_op& Op)
  {
  super_op LOp;
  matrix mx = Op.get_mx();			// Get Hilbert matrix
  basis bs  = Op.get_basis();			// Get Hilbert basis
  int nc    = bs.sub_N();			// Number of sub-spaces

// 		This Is Normal Non-Composite Superoperator

  if(nc == 1) return left(mx,bs);		// If only one sub-space

// 		This Is Composite Liouville Space Superoperator
//	      (Now We Must Deal With All Sub-Spaces Individually)

  LOp.HSp = Op.dim();				// Set Hilbert dimension
  LOp.LSp = bs.dim_LS();			// Set Liouville dimension 
  if(!LOp.LSp) return LOp;			// Exit if empty Op
  matrix *mxc;					// Array of Liouville sub-matrices
  mxc = new matrix[nc];				// Array of Liouville sub-matrices
  int *ncd;					// Array of Liouville sub-space dims
  ncd = new int[nc];
  int cmp, ist, idim, ls=0;			// Indicies for sub-spaces
  matrix I, mxs;				// Arrays for sub-space treatment
  for(cmp=0; cmp<nc; cmp++)			// Loop subspaces (components)
    {
    ist = bs.sub_anchor(cmp);			//	Start of sub-space
    idim = bs.sub_dim(cmp);			//	Dimension of sub-space
    mxs = mx.get_block(ist,ist,idim,idim);	//	Op matrix, this component
    I = matrix(idim,idim,i_matrix_type);	// 	Required Identity matrix
    mxc[cmp] = tensor_product(mxs, I);		//	LOp matrix, this component
    ncd[cmp] = mxc[cmp].rows();			//	Store subspace dimension
    ls += ncd[cmp];				//	Track total Liouville space
    }
  matrix LMx(ls,ls, complex0);			// Array for superop matrix
  ist = 0;					// Starting index of subspace
  for(cmp=0; cmp<nc; cmp++)			// Loop subspaces (components)
    {
    idim = (mxc[cmp]).rows();			//	Dimension of sub-space
    LMx.put_block(ist, ist, mxc[cmp]);		//	Put sub-space array in
    ist += idim;				// 	Update subspace index
    }
  LOp.mx = LMx; 				// Set Liouville matrix
  I = matrix(ls,ls,i_matrix_type);		// Identity  matrix for basis
  LOp.Lbs = basis(I, nc, ncd);			// Set Liouville basis (Default)
  LOp.Hbs = bs;					// LOp Hilbert basis = Op basis
  delete [] mxc;
  delete [] ncd;
  return LOp;
  }

super_op left(const matrix& mx, const basis& bs)
  {
  int hs = mx.rows();				// Hilbert space dimension
  if(!hs) return super_op(); 			// Exit if empty array
  matrix I = matrix(hs, hs, i_matrix_type);
  return super_op(tensor_product(mx,I), bs);
  }

super_op left(const matrix& mx)
  {
  int hs = mx.rows();				// Hilbert space dimension
  if(!hs) return super_op(); 			// Exit if empty array
  matrix I = matrix(hs, hs, i_matrix_type);
  return super_op(tensor_product(mx,I));
  }

super_op right(const gen_op& Op)
  {
  super_op LOp;
  matrix mx = Op.get_mx();			// Get Hilbert matrix
  basis bs = Op.get_basis();			// Get Hilbert basis
  int nc = bs.sub_N();				// Number of sub-spaces

// 		This Is Normal Non-Composite Superoperator

  if(nc == 1) return right(mx,bs);		// If only one sub-space

// 		This Is Composite Liouville Space Superoperator
//	      (Now We Must Deal With All Sub-Spaces Individually)

  LOp.HSp = Op.dim();				// Set Hilbert dimension
  LOp.LSp = bs.dim_LS();			// Set Liouville dimension 
  if(!LOp.LSp) return LOp;				// Exit if empty Op
  matrix *mxc;				// Array of Liouville sub-matrices
  int *ncd;					// Array of Liouville sub-space dims
  mxc = new matrix[nc];
  ncd = new int[nc];
  int cmp, ist, idim, ls=0;			// Indicies for sub-spaces
  matrix I, mxs;				// Arrays for sub-space treatment
  for(cmp=0; cmp<nc; cmp++)			// Loop subspaces (components)
    {
    ist = bs.sub_anchor(cmp);			//	Start of sub-space
    idim = bs.sub_dim(cmp);			//	Dimension of sub-space
    mxs = mx.get_block(ist,ist,idim,idim);	//	Op matrix, this component
    I = matrix(idim,idim,i_matrix_type);	// 	Required Identity matrix
    mxc[cmp]=tensor_product(I, transpose(mxs));	//	LOp matrix, this component
    ncd[cmp] = mxc[cmp].rows();			//	Store subspace dimension
    ls += ncd[cmp];				//	Track total Liouville space
    }
  matrix LMx(ls,ls, complex0);			// Array for superop matrix
  ist = 0;					// Starting index of subspace
  for(cmp=0; cmp<nc; cmp++)			// Loop subspaces (components)
    {
    idim = (mxc[cmp]).rows();			//	Dimension of sub-space
    LMx.put_block(ist, ist, mxc[cmp]);		//	Put sub-space array in
    ist += idim;				// 	Update subspace index
    }
  LOp.mx = LMx; 				// Set Liouville matrix
  I = matrix(ls,ls,i_matrix_type);		// Identity  matrix for basis
  LOp.Lbs = basis(I, nc, ncd);			// Set Liouville basis (Default)
  LOp.Hbs = bs;					// LOp Hilbert basis = Op basis
  delete [] mxc;
  delete [] ncd;
  return LOp;
  }

super_op right(const matrix& mx, const basis& bs)
  {
  int hs = mx.rows();				// Hilbert space dimension
  if(!hs) return super_op(); 			// Exit if empty array
  matrix I = matrix(hs, hs, i_matrix_type);
  return super_op(tensor_product(I,mx), bs);
  }

super_op right(const matrix& mx)
  {
  int hs = mx.rows();				// Hilbert space dimension
  if(!hs) return super_op(); 			// Exit if empty array
  matrix I = matrix(hs, hs, i_matrix_type);
  return super_op(tensor_product(I,mx));
  }


// ----------------------------------------------------------------------------
//                          Commutation Superoperators
// ----------------------------------------------------------------------------

//                           t
//  LOp = Op (X) E - E (x) Op          where       LOp*Op1 = [Op, Op1]

	// Input		Op    : General operator
	// Output		LOp   : Superoperator equal to
	//			        the commutation operator

super_op commutator(const gen_op& Op)
  {
  super_op LOp;
  matrix mx = Op.get_mx();			// Get Hilbert matrix
  basis bs = Op.get_basis();			// Get Hilbert basis
  int nc = bs.sub_N();				// Number of sub-spaces

// 		This Is Normal Non-Composite Superoperator

  if(nc == 1)					// If only one sub-space
    {
    LOp = commutator(mx);			// Use function overload
    LOp.Hbs = Op.get_basis();			// Hilbert basis = Op basis
    return LOp;
    }

// 		This Is Composite Liouville Space Superoperator
//	      (Now We Must Deal With All Sub-Spaces Individually)

  LOp.HSp = Op.dim();				// Set Hilbert dimension
  LOp.LSp = bs.dim_LS();			// Set Liouville dimension 
  if(!LOp.LSp) return LOp;			// Exit if empty Op
  matrix *mxc;					// Array Liouv. sub-matrices
  int *ncd;					// Array Liouv. sub-space dims
  mxc = new matrix[nc];
  ncd = new int[nc];
  int cmp, ist, idim, ls=0;			// Indicies for sub-spaces
  matrix I, mxs;				// Arrays sub-space treatment
  for(cmp=0; cmp<nc; cmp++)			// Loop subspaces (components)
    {
    ist = bs.sub_anchor(cmp);			//	Start of sub-space
    idim = bs.sub_dim(cmp);			//	Dimension of sub-space
    mxs = mx.get_block(ist,ist,idim,idim);	//	Op matrix, this compon.
    I = matrix(idim,idim,i_matrix_type);	// 	Required Identity matrix
    mxc[cmp] = tensor_product(mxs, I)		//	LOp matrix, this compon.
             - tensor_product(I,transpose(mxs));
    ncd[cmp] = mxc[cmp].rows();			//	Store subspace dimension
    ls += ncd[cmp];				//	Track total Liouv. space
    }
  matrix LMx(ls,ls, complex0, d_matrix_type);	// Array for superop matrix
  ist = 0;					// Starting index of subspace
  for(cmp=0; cmp<nc; cmp++)			// Loop subspaces (components)
    {
    idim = (mxc[cmp]).rows();			//	Dimension of sub-space
    LMx.put_block(ist, ist, mxc[cmp]);		//	Put sub-space array in
    ist += idim;				// 	Update subspace index
    }
  LOp.mx = LMx; 				// Set Liouville matrix
  I = matrix(ls,ls,i_matrix_type);		// Identity  matrix for basis
  LOp.Lbs = basis(I, nc, ncd);			// Set Liouville basis (Default)
  LOp.Hbs = bs;					// LOp Hilbert basis = Op basis
  delete [] mxc;
  delete [] ncd;
  return LOp;
  }


super_op commutator(const matrix& mx)

	// Input		mx    : Matrix
	// Output		LOp   : Superoperator equal to
	//			        the commutation operator

  {
  super_op LOp;
  int hs = mx.rows();				// Hilbert space dimension 
  if(hs!=mx.cols())				// Insure square input matrix
    {
    LOp.LOperror(5, "commutator", 1);	// Trouble in commutator
    LOp.LOpfatal(41);			// Can't mix LOp & mx
    }
  LOp.HSp = hs;					// Set Hilbert space
  LOp.LSp = hs*hs;				// Set Liouville space
  if(!hs) return LOp;				// Exit if NULL Op
  matrix I = matrix(hs, hs, i_matrix_type);	// Make required I matrix
  LOp.mx = tensor_product(mx, I);		// Set mx (X) I
  LOp.mx -= tensor_product(I, transpose(mx));	// Subtract I (X) mxt
  LOp.Hbs = basis(LOp.HSp);			// Default Hilbert space basis 
  LOp.Lbs = basis(LOp.LSp);			// Default Liouville space basis
  return LOp;
  }


// ----------------------------------------------------------------------------
//                      Double Commutation Superoperators
// ----------------------------------------------------------------------------

// sosi - these still need to be made valid for composite Liouville space
//                                       T           T           T     T 
//        LOp = Op1*Op2 (X) E - Op1 X Op2 - Op2 X Op1 + E (X) Op1 * Op2
// where
//                         LOp*Op3 = [Op1,[Op2, Op3]]


super_op d_commutator(const gen_op& Op, const complex& z)

	// Input		Op    : General operator
	// Output		LOp   : Superoperator equal to
	//			        the double commutation operator
	//				LOp*Op1 = z*[Op,[Op, Op1]]

  {
  super_op LOp;
  //  int hs = SuperOp.dim(Op);				// Hilbert space dimension
  int hs = Op.dim();				// Hilbert space dimension
  if(!hs) return LOp;				// Exit if NULL Op
  LOp.HSp = hs;
  LOp.LSp = hs*hs;
  matrix I = matrix(hs, hs, i_matrix_type);
  matrix Opmx = Op.get_mx();
  matrix z2Opmx = (z*2.0)*Opmx;
  matrix zOp_sq = z * (Opmx*Opmx);

  basis bs = Op.get_basis();			// Get Hilbert basis
  int nc = bs.sub_N();				// Number of sub-spaces

// 		This Is Normal Non-Composite Superoperator

  if(nc == 1)					// If only one sub-space
    {
    LOp.mx = tensor_product(zOp_sq, I);
    LOp.mx -= tensor_product(z2Opmx, transpose(Opmx));
    LOp.mx += tensor_product(I, transpose(zOp_sq));
    LOp.Hbs = Op.get_basis();			// Hilbert space basis = Op basis
    LOp.Lbs = basis(LOp.LSp);			// Default Liouville space basis
    return LOp;
    }


// 		This Is Composite Liouville Space Superoperator
//	      (Now We Must Deal With All Sub-Spaces Individually)
//
  gen_op GOp;
  int ls = 0;
  int cmp;
  matrix *mxc, *bsc;					// Array Liouv. sub-matrices
  mxc = new matrix[nc];
  bsc = new matrix[nc];
  for(cmp=0; cmp<nc; cmp++)			// Loop subspaces (components)
    {
    GOp = Op.project_sub(cmp);
    GOp.set_DBR();
    hs = GOp.dim();
    I = matrix(hs, hs, i_matrix_type);
    ls += hs*hs;
    
    Opmx = GOp.get_mx();
    z2Opmx = (z*2.0)*Opmx;
    zOp_sq = z * (Opmx*Opmx);
    
    LOp.mx = tensor_product(zOp_sq, I);
    LOp.mx -= tensor_product(z2Opmx, transpose(Opmx));
    LOp.mx += tensor_product(I, transpose(zOp_sq));

    mxc[cmp] = LOp.get_mx();			//store mutual exchange matrix, this compon.
    bsc[cmp] = GOp.get_basis().get_mx();	//store Hilbert space basis matrix, this compon.
    }

  LOp = super_op(mxc, nc, bsc);
  delete [] mxc;
  delete [] bsc;

  return LOp;
  }

	// Input		mx    : Matrix
	// Output		LOp   : Superoperator equal to
	//			        the double commutation operator
	//				LOp*mx1 = [mx,[mx, mx1]]

super_op d_commutator(const matrix& mx)
  {
  super_op LOp;
  int hs = mx.rows();				// Hilbert space dimension
  if(hs != mx.cols())				// Insure the input martix
    {						// is square.
    LOp.LOperror(5, "d_commutator", 1);	// Problems in d_commutator
    LOp.LOpfatal(41);			// Can't mix matrix with LOp
    }
  LOp.HSp = hs;					// Set Hilbert space dimension
  LOp.LSp = hs*hs;				// Set Liouville space dimension
  if(!hs) return LOp;				// Exit if for NULL matrix
  complex z(2);					// Needed cmplx 2 for *
  matrix I = matrix(hs, hs, i_matrix_type);	// Needed I matrix 
  matrix mx_sq = mx * mx;			// Needed mx*mx
  LOp.mx = tensor_product(mx_sq, I);		// Set mx*mx (X) I
  LOp.mx -= tensor_product(z*mx, transpose(mx));// Subtract 2*mx (X) I
  LOp.mx += tensor_product(I, transpose(mx_sq));// Add I (X) mxt*mxt
  LOp.Hbs = basis(LOp.HSp);			// Default Hilbert space basis 
  LOp.Lbs = basis(LOp.LSp);			// Default Liouville space basis
  return LOp;
  }




	// Input		Op1   : General operator
	// 			Op2   : General operator
	// Output		LOp   : Double commutation superoperator
	//				LOp = [Op1, [Op2, ]]
	// Note			      : No operator basis checking
	// Note			      : A NULL Op returns NULL LOp
	// Note			      : No check is made to see if
	//				Op1=Op2 for faster function

super_op d_commutator(const gen_op& Op1, const gen_op& Op2)
  {
  super_op LOp;
  if(Op1.exists() && Op2.exists())	// Check for NULL Op1 or Op2
    {
    LOp.HSp = Op1.dim();
    LOp.LSp = LOp.HSp*LOp.HSp;
    matrix I = matrix(LOp.HSp, LOp.HSp, i_matrix_type);
    matrix Op1Op2 = Op1.get_mx() * Op2.get_mx();
    matrix Op1Op2T = transpose(Op1.get_mx()) * transpose(Op2.get_mx());
    LOp.mx = tensor_product(Op1Op2, I);
    LOp.mx -= tensor_product(Op1.get_mx(), transpose(Op2.get_mx()));
    LOp.mx -= tensor_product(Op2.get_mx(), transpose(Op1.get_mx()));
    LOp.mx += tensor_product(I, Op1Op2T);
    LOp.Hbs = Op1.get_basis();		// LOp Hilbert basis = Op1 basis
    LOp.Lbs = basis(LOp.LSp);		// LOp Liouville basis = default
    }
  return LOp;
  }


super_op d_commutator(const gen_op& Op1, const gen_op& Op2, const complex& z)
//  return LOp;

	// Input		Op1   : General operator
	// 			Op2   : General operator
	//			z     : Scaling factor
	// Output		LOp   : Double commutation superoperator
	//				LOp = [Op1, [Op2, ]]
	// Note			      : No operator basis checking
	// Note			      : A NULL Op returns NULL LOp
	// Note			      : No check is made to see if
	//				Op1=Op2 for faster function

  {
  super_op LOp;
  if(Op1.exists() && Op2.exists())	// Check for NULL Op1 or Op2
    {
    LOp.HSp = Op1.dim();
    LOp.LSp = LOp.HSp*LOp.HSp;
    matrix I = matrix(LOp.HSp, LOp.HSp, i_matrix_type);
    matrix Op2mx = Op2.get_mx();
    matrix zOp2mx = z * Op2mx;
    matrix Op1mx = Op1.get_mx();
    matrix zOp1Op2 = Op1mx*zOp2mx;
    matrix zOp1Op2T = transpose(Op1mx) * transpose(zOp2mx);
    LOp.mx = tensor_product(zOp1Op2, I);
    LOp.mx -= tensor_product(Op1mx, transpose(zOp2mx));
    LOp.mx -= tensor_product(zOp2mx, transpose(Op1mx));
    LOp.mx += tensor_product(I, zOp1Op2T);
    LOp.Hbs = Op1.get_basis();		// LOp Hilbert basis = Op1 basis
    LOp.Lbs = basis(LOp.LSp);		// LOp Liouville basis = default
    }
  return LOp;
  }


super_op d_commutator(const matrix& mx1, const matrix& mx2)

	// Input		mx1   : Matrix
	// 			mx2   : Matrix
	// Output		LOp   : Superoperator equal to
	//			        the double commutation operator
	//				LOp*Op = [mx1, [mx2, Op]]

//                                     T           T           T     T 
//      LOp = mx1*mx2 (X) E - mx1 X mx2 - mx2 X mx1 + E (X) mx1 * mx2

  {
  super_op LOp;
  if(mx1 == mx2) return d_commutator(mx1);	// Faster to do this if mx1 = mx2
  int hs = mx1.rows();				// This is Hilbert space dimension 
  if(hs!=mx1.cols() || mx2.rows()!=mx2.cols())	// Insure the input matrices
    {						// are square.
    LOp.LOperror(5, "d_commutator", 1);	// Problems in d_commutator
    LOp.LOpfatal(41);			// Can't mix matrix with LOp
    }
  if(hs != mx2.rows())				// Insure both arrays reside in
    {						// the same Hilbert space
    LOp.LOperror(5, "d_commutator", 1);	// Problems in d_commutator
    LOp.LOpfatal(41);			// Can't mix LOp & mx
    }
  LOp.HSp = hs;					// Set Hilbert space dimension
  LOp.LSp = hs*hs;				// Set Liouville space dimension
  if(!hs) return LOp;				// Quit if formed from NULL mx's
  matrix I = matrix(hs, hs, i_matrix_type);	// Begin with I matrix
  matrix m1m2 = mx1 * mx2;			// Product mx1*mx2 in H.S.
  matrix m1m2T = transpose(mx1)*transpose(mx2);	// Product mx1t*mx2t in HS
  LOp.mx = tensor_product(m1m2, I);		// Make (mx1t*mx2t) X I in LS
  LOp.mx -= tensor_product(mx1,transpose(mx2));	// Subtract mx1 X mx2t in LS
  LOp.mx -= tensor_product(mx2,transpose(mx1));	// Subtract mx2t X mx1 in LS
  LOp.mx += tensor_product(I, m1m2T);		// Add I x (mx1t*mx2t) 
  LOp.Hbs = basis(LOp.HSp);			// Set default Hilbert basis 
  LOp.Lbs = basis(LOp.LSp);			// Set default Liouville basis
  return LOp;
  }

// ----------------------------------------------------------------------------
//                    Unitary Transformation Superoperators
// ----------------------------------------------------------------------------

//                     *                                                    -1
//	LOp = Op (X) Op                 where        LOp*Op1 = Op * Op1 * Op 


super_op U_transform(const gen_op& Op)

	// Input		Op    : General operator
	// Output		LOp   : Unitary transformation superoperator
	//			        formed from the input operator Op

  {
  super_op LOp;
  matrix mx = Op.get_mx();			// Get Hilbert matrix
  basis bs = Op.get_basis();			// Get Hilbert basis
  int nc = bs.sub_N();				// Number of sub-spaces

// 		This Is Normal Non-Composite Superoperator

  if(nc == 1)					// If only one sub-space
    {						// then quick and dirty
    LOp = U_transform(mx);			// by using overload and
    LOp.Hbs = bs;				// assigning Hilbert basis
    return LOp;
    }

// 		This Is Composite Liouville Space Superoperator
//	      (Now We Must Deal With All Sub-Spaces Individually)

  LOp.HSp = Op.dim();				// Set Hilbert dimension
  LOp.LSp = bs.dim_LS();			// Set Liouville dimension 
  if(!LOp.LSp) return LOp;			// Exit if empty Op
  matrix *mxc;				// Array of Liouville sub-matrices
  int *ncd;					// Array of Liouville sub-space dims
  mxc = new matrix[nc];
  ncd = new int[nc];
  int cmp, ist, idim, ls=0;			// Indicies for sub-spaces
  matrix mxs;					// Array for sub-space
  for(cmp=0; cmp<nc; cmp++)			// Loop subspaces (components)
    {
    ist = bs.sub_anchor(cmp);			//	Start of sub-space
    idim = bs.sub_dim(cmp);			//	Dimension of sub-space
    mxs = mx.get_block(ist,ist,idim,idim);	//	Op matrix, this component
    mxc[cmp] = tensor_product(mxs,conj(mxs));	//	LOp matrix, this component
    ncd[cmp] = mxc[cmp].rows();			//	Store subspace dimension
    ls += ncd[cmp];				//	Track total Liouville space
    }
  matrix LMx(ls,ls, complex0);			// Array for superop matrix
  ist = 0;					// Starting index of subspace
  for(cmp=0; cmp<nc; cmp++)			// Loop subspaces (components)
    {
    idim = (mxc[cmp]).rows();			//	Dimension of sub-space
    LMx.put_block(ist, ist, mxc[cmp]);		//	Put sub-space array in
    ist += idim;				// 	Update subspace index
    }
  LOp.mx = LMx; 				// Set Liouville matrix
  matrix Imx(ls,ls,i_matrix_type);		// Identity  matrix for basis
  LOp.Lbs = basis(Imx, nc, ncd);		// Set Liouville basis (Default)
  LOp.Hbs = bs;					// LOp Hilbert basis = Op basis
  delete [] ncd;
  delete [] mxc;
  return LOp;
  }


super_op U_transform(const matrix& mx)

	// Input		mx    : matrix
	// Output		LOp   : Superoperator (this) equal to
	//			        the unitary transformation operation

  {
  super_op LOp;
  int hs = mx.rows();			// Hilbert space dimension
  LOp.HSp = hs;				// Set Hilbert space
  LOp.LSp = hs*hs;			// Set Liouville space
  if(!hs) return LOp;			// Exit if NULL matrix
  LOp.mx = tensor_product(mx, conj(mx));// Form mx (X) mx*
  LOp.Hbs = basis(LOp.HSp);		// Set default Hilbert basis 
  LOp.Lbs = basis(LOp.LSp);		// Set default Liouville basis
  return LOp;
  }

// ----------------------------------------------------------------------------
//                         Projection Superoperators
// ----------------------------------------------------------------------------

	// Input		Op    : General operator
	// Output		LOp   : Projection superoperator
	//				LOp*Op1 = z*Op, w/ Op1 = z*Op + .....
	//					     t
	//				LOp = Op(x)Op
super_op project(const gen_op& Op)
// sosi this function does not work properly
  {
  super_op LOp;
  LOp.HSp = dim(Op);
  LOp.LSp = LOp.HSp*LOp.HSp;
  if (LOp.LSp)			// Check for NULL Op
    {
    matrix mxt = transpose(Op.get_mx()); 
    complex z = 1.0/(trace(Op.get_mx()*mxt));
//    LOp.mx = tensor_product(z*Op.get_mx(), mxt);
    LOp.mx = tensor_product(mxt,z*Op.get_mx());
    LOp.Hbs = Op.get_basis();	// LOp Hilbert space basis equal to Op basis
    LOp.Lbs = basis(LOp.LSp);	// LOp Liouville space basis to default basis
    }
  return LOp;
  }



	// Input		mx    : matrix
	// Output		LOp   : Projection superoperator
	//				LOp*Op1 = z*mx, w/ Op1 = z*mx + .....
	//					     t
	//				LOp = mx(x)mx
// sosi this function does not work properly

super_op project(const matrix& mx)
{
  super_op LOp;
  LOp.HSp = mx.rows();
  LOp.LSp = LOp.HSp*LOp.HSp;
  if (LOp.LSp)			// Check for NULL Op
    {
    matrix mxt = transpose(mx);
    complex z = 1.0/(trace(mx*mxt));
//    LOp.mx = tensor_product(z*mx, mxt);
    LOp.mx = tensor_product(z*mxt, mx);
    LOp.Hbs = basis(LOp.HSp);	// LOp Hilbert space basis equal to default basis
    LOp.Lbs = basis(LOp.LSp);	// LOp Liouville space basis to default basis
    }
  return LOp;
}

// ----------------------------------------------------------------------------
//                          Exponential Superoperators
// ----------------------------------------------------------------------------

  


	// Input		LOp   : Superoperator (this)
        // Return		ExpLOp: Exponential of LOp
	//				ExpLOp = exp(LOp)
        // Note			      : Exponential output in EBR of LOp
        // Note			      : L0p's EBR is generated herein

super_op super_op::exp() const
  {
  if(!HSp)				// Check for NULL LOp
    {
    LOperror(5, "exp", 1);		// Error during LOp exp function
    LOpfatal(7);			// Accessing Null LOp
    }
  set_EBR();				// Put LOp in its EBR
  super_op ExpLOp(*this);		// Copy LOp in its EBR
  complex z;
  for(int i=0; i<ExpLOp.LSp; i++)	// Exponentiate the diagonal
    {
    z = ExpLOp.get(i,i);
    ExpLOp.put(i,i,z.Zexp());
    }
  return ExpLOp;
  }


super_op super_op::exp(const complex& t, double cutoff) const

        // Input                LOp   : Superoperator (this)
        //                      t     : Exponential factor
        //                      cutoff: Exponential factor roundoff
        // Return               ExpLOp: Exponential of LOp
        //                              ExpLOp = exp(t*LOp)
        // Note                       : Exponential output in EBR of LOp
        // Note                       : L0p's EBR is generated herein
        // Note                       : Value of t is considered 0 if
        //                              it's magnituded is less than cutoff

  {
  if(!HSp)				// Check for NULL LOp
    {
    LOperror(5, "exp", 1);		// Error during LOp exp function
    LOpfatal(7);			// Accessing NULL Lop
    }
  set_EBR();				// Put LOp in its EBR
  super_op ExpLOp(*this);		// Copy LOp in its EBR
  if(norm(t) < fabs(cutoff))
    {
    matrix mx(LSp, LSp, i_matrix_type);	// ExpLOp is the identity superop
    ExpLOp.mx = mx;
    }
  else
    {
    complex z;
    for(int i=0; i<LSp; i++)		// Exponentiate the diagonal
      {
      z = t*ExpLOp.get(i,i);
      ExpLOp.put(i,i,z.Zexp());
      }
    }
  return ExpLOp;
  }

  


	// Input		LOp1  : Superoperator
        // Return		LOp   : Exponential of LOp1
	//				LOp = exp(LOp1)
        // Note			      : Computed in EBR of LOp1
        // Note			      : Superoperator output in EBR of LOp1

super_op exp(const super_op& LOp1)
  {
  super_op LOp(matrix(LOp1.LSp, LOp1.LSp, d_matrix_type));
  if(LOp1.HSp)				// Check for NULL LOp
    {
    LOp1.set_EBR();			// Put LOp1 in EBR
    for(int i=0; i<LOp.LSp; i++)
      LOp.put(i,i,exp(LOp1(i,i)));
    LOp.Hbs = LOp1.Hbs;			// LOp Hilbert space basis same as LOp1's
    LOp.Lbs = LOp1.Lbs;			// LOp Liouville space basis same as LOp1's
    return LOp;
    }
  else
    {
    LOp1.LOperror(5, "exp", 1);	// Error during LOp exp function
    LOp1.LOpfatal(7);		// Accessing Null LOp
    }
  return LOp;
  }

  
super_op exp(const super_op& LOp1, const complex& t)
//return LOp(matrix(LOp1.LSp, LOp1.LSp, d_matrix_type));

	// Input		LOp1  : Superoperator
	//			t     : A time
        // Return		LOp   : Exponential of LOp1
	//				LOp = exp(LOp1*t)
        // Note			      : Computed in EBR of LOp1

  {
  super_op LOp(matrix(LOp1.LSp, LOp1.LSp, d_matrix_type));
  if(LOp.HSp) 				// Check for NULL LOp
    {
    if(t != complex0)
      {
      LOp1.set_EBR();			// First put LOp1 in its eigenbasis
      for(int i=0; i<LOp.LSp; i++)	// Exponentiate diagonal elements
        LOp.put(i,i,exp(LOp1(i,i)*t));
      }
    else
      {
      matrix mx(LOp1.LSp, LOp1.LSp, i_matrix_type);	// LOp is the identity superoperator
      LOp.mx = mx;
      }
    LOp.Hbs = LOp1.Hbs;			// LOp Hilbert space basis same as LOp1's
    LOp.Lbs = LOp1.Lbs;			// LOp Liouville space basis same as LOp1's
    }
  else
    {
    LOp1.LOperror(5, "exp", 1);	// Error during LOp exp function
    LOp1.LOpfatal(7);		// Accessing Null LOp
    }
  return LOp;
  }

// ----------------------------------------------------------------------------
//                 Exponential Superoperators using Pade Method
// ----------------------------------------------------------------------------


	// Input		LOp   : Superoperator (this)
        // Return		ExpLOp: Exponential of LOp
	//				ExpLOp = exp(LOp) (Pade method)
        // Note			      : Exponential output in same base as LOp

super_op super_op::expm() const
  {
  if(!HSp)				// Check for NULL LOp
    {
    LOperror(5, "exp", 1);		// Error during LOp exp function
    LOpfatal(7);			// Accessing Null LOp
    }
  super_op ExpLOp(*this);		// Copy LOp in its current base
  matrix mx1 = this->mx;
  
  ExpLOp.mx = mx1.expm();
  
  return ExpLOp;
  }


super_op super_op::expm(const complex& t, double cutoff) const

        // Input                LOp   : Superoperator (this)
        //                      t     : Exponential factor
        //                      cutoff: Exponential factor roundoff
        // Return               ExpLOp: Exponential of LOp
        //                              ExpLOp = exp(t*LOp)
        // Note                       : Exponential output in EBR of LOp
        // Note                       : L0p's EBR is generated herein
        // Note                       : Value of t is considered 0 if
        //                              it's magnituded is less than cutoff

  {
  if(!HSp)				// Check for NULL LOp
    {
    LOperror(5, "exp", 1);		// Error during LOp exp function
    LOpfatal(7);			// Accessing NULL Lop
    }
  super_op ExpLOp(*this);		// Copy LOp in its current base
  if(norm(t) < fabs(cutoff))
    {
    matrix mx(LSp, LSp, i_matrix_type);	// ExpLOp is the identity superop
    ExpLOp.mx = mx;
    }
  else
    {
    matrix mx1 = this->mx*t;
    ExpLOp.mx = mx1.expm();
    }
  return ExpLOp;
  }

  


	// Input		LOp1  : Superoperator
        // Return		LOp   : Exponential of LOp1
	//				LOp = exp(LOp1)
        // Note			      : Computed in EBR of LOp1
        // Note			      : Superoperator output in EBR of LOp1

super_op expm(const super_op& LOp1)
  {
  super_op LOp = LOp1;
  if(LOp1.HSp)				// Check for NULL LOp
    {
    matrix m1 = LOp.get_mx();
    LOp.put_mx(m1.expm());
    return LOp;
    }
  else
    {
    LOp1.LOperror(5, "exp", 1);	// Error during LOp exp function
    LOp1.LOpfatal(7);		// Accessing Null LOp
    }
  return LOp;
  }

  
super_op expm(const super_op& LOp1, const complex& t)
//return LOp(matrix(LOp1.LSp, LOp1.LSp, d_matrix_type));

	// Input		LOp1  : Superoperator
	//			t     : A time
        // Return		LOp   : Exponential of LOp1
	//				LOp = exp(LOp1*t)
        // Note			      : Computed in EBR of LOp1

  {
  super_op LOp = LOp1;
  if(LOp.HSp) 				// Check for NULL LOp
    {
    if(t != complex0)
      {
      matrix m1 = LOp.get_mx()*t;
      LOp.put_mx(m1.expm());
      return LOp;
      }
    else
      {
      matrix mx(LOp1.LSp, LOp1.LSp, i_matrix_type);	// LOp is the identity superoperator
      LOp.mx = mx;
      }
    }
  else
    {
    LOp1.LOperror(5, "exp", 1);	// Error during LOp exp function
    LOp1.LOpfatal(7);		// Accessing Null LOp
    }
  return LOp;
  }


// ----------------------------------------------------------------------------
//                      Superoperators Taken To A Power
// ----------------------------------------------------------------------------

	// Input		LOp   : Superoperator (this)
	//			power : A power factor
        // Return		LOp1  : Initial superoperator taken
	//				to a particular power
        // Note			      : Computed in EBR of LOp1

super_op pow(const super_op& LOp, int power)
  {
  super_op LOp1(matrix(LOp.LSp, LOp.LSp, d_matrix_type));
  complex z, zpow(power);
  if(LOp.HSp) 			// Check for NULL LOp
    {
    LOp.set_EBR();		// First put LOp in its eigenbasis
    for(int i=0; i<LOp.LSp; i++)
      {
      z = LOp.get(i,i);
      LOp1.put(i,i, pow(z, zpow));
      }
    LOp1.Hbs = LOp.Hbs;	// LOp Hilbert space basis same as LOp1's
    LOp1.Lbs = LOp.Lbs;	// LOp Liouville space basis same as LOp1's
    }
  else
    {
    LOp1.LOperror(5, "pow", 1);	// Error during LOp pow function
    LOp1.LOpfatal(7);		// Accessing Null LOp
    }
  return LOp1;
  }

// ____________________________________________________________________________
// F                  SUPER OPERATOR BASIS MANIPULATIONS
// ____________________________________________________________________________


void super_op::set_EBR() const

	// Input		LOp  : Superoperator (this)
	// Output		none : LOp set to eigenbasis
	// Note			     : Hilbert space basis unaltered

  {
  matrix mxd,bsd;
  if(HSp) 			// Check for NULL LOp
    if(d_matrix_type		// Check for EBR existence 
	   != mx.stored_type())
      {
      set_HBR();		// First set to default basis
      diag(mx, mxd, bsd);	// Diagonalize LOp
      mx = mxd;			// Set Lop matrix to diagonal form
      Lbs = basis(bsd);		// Set Liouville space matrix
      }
  return;
  }

// sosi 9/9/97 - Letting basis handle unitary vs. non-unitary issue now.
//               Thus, checking is done along rather than forced invert
	
void super_op::set_HBR() const

	// Input		LOp  : Superoperator (this)
	// Output		none : LOp set to default Liouville space basis
	// Note			     : LOp Hilbert space basis is unaltered

  {
  if(!HSp) return; 			// Do nothing on a NULL LOp
  if(Lbs.isDefaultBasis()) return;	// LOp LS basis already default!
  if(mx.stored_type() != i_matrix_type)	// If there are bases to deal with
    mx = Lbs.convert_back(mx);		// kick back out of them
//    mx = Lbs.U() * mx * inv(Lbs.U());	// kick back out of them
  Lbs = basis(LSp);			// Liouville basis to default basis
  return;
  }

	
void super_op::set_DBR() const

	// Input		LOp  : Superoperator (this)
	// Output		none : LOp set into a default Hilbert
	//			       space basis.  When this is done
	// Note			     : LOp Hilbert space basis is unaltered

  {
  if(!HSp) return; 			// Do nothing on a NULL LOp
  if(Hbs.isDefaultBasis()) return;	// LOp HS basis already default!
  set_HBR();				// Insure default Liouville space
  basis S = get_basis();		// Get LOp Hilbert space basis
  basis DS = defbasis(S);		// Default basis with S structure
  gen_op OpX(S.U(), DS);		// Dummy operator to enable U_trans
  super_op LOpX = U_transform(OpX);	// Dummy superop to get U_trans of S
  basis LS(LOpX.mx);			// Here is basis transformation matrix
  mx = LS.convert_back(mx);		// Kick LOp matrix out of HS basis
  Hbs = defbasis(S);			// Set default Hilbert space basis
  }



	// Input		LOp  : Superoperator (this).
	//			LOp1 : Superoperator.
	// Return		none : Superoperator LOp is put
	//			       into the Liouville basis of LOp1
	// Note			     : This function performs similarity
	//			       transformations in Liouville space.
	//			       It is likely to be computationally
	//			       intensive, avoid it if possible.

// sosi 9/9/97 - Letting basis handle unitary vs. non-unitary issue now.
//               Thus, checking is done along rather than forced invert
//    mx = inv(LOp1.Lbs.U()) * mx * LOp1.Lbs.U();

void super_op::LOp_base(const super_op& LOp1) const
  {
  if(Lbs != LOp1.Lbs)			// Only if LOp & LOp1 bases differ!
    {
    set_HBR();				// Put Lop in default LS basis
    mx = (LOp1.Lbs).convert(mx);	// Put LOp matrix into basis of LOp1
    Lbs = LOp1.Lbs;	 		// LOp LS basis to LOp1 LS basis
    }
  return;
  }


void super_op::LOp_Hbase(const super_op& LOp1, int warn) const

	// Input		LOp  : Superoperator (this).
	//			LOp1 : Superoperator.
        //                      warn : Flag if warnings desired (off)
	// Return		none : Superoperator LOp is put
	//			       into the Hilbert basis of LOp1
	// Note			     : This function performs similarity
	//			       transformations in Liouville space.
	//			       It is likely to be computationally
	//			       intensive, avoid it if possible.

//                            -1                      *              t
//    LOp(hbs1) = U*LOp(hbs)*U       where U = u (X) u   and u = hbs1 * hbs    

// sosi: 9/9/97 -As in function set_HBR, there's  difficulty in non-unitary basis
// transformations, but for now assume the adjoint will suffice.  I've rewritten
// the code to hopefully handle Liouville composite spaces, but there are still
// a few questions regarding the basis handling.  Testing will insure it is proper.

  {
  if(Hbs == LOp1.Hbs) return;		// Nothing if LOp & LOp1 share Hbs
  if(warn) LOperror(90);		// Warn of LOp basis mixing if desired
  set_DBR();				// Put in default LS and HS bases
  basis S = LOp1.get_basis();		// Get LOp1 Hilbert space basis
  basis DS = defbasis(S);		// Default basis with S structure
  gen_op OpX(adjoint(S.U()), DS);	// Dummy operator to enable U_trans
  super_op LOpX = U_transform(OpX);	// Dummy superop to get U_trans of S
  basis LS(LOpX.mx);			// Here is basis transformation matrix
  mx = LS.convert_back(mx);		// Kick LOp matrix out of HS basis
  Hbs = S;				// Set new Hilbert space basis

//  set_HBR();
//  basis S = get_basis();		// Get LOp Hilbert space basis
//  basis S1 = LOp1.get_basis();		// Get LOp1 Hilbert space basis
//  matrix U = adjoint_times(S1.U(), S.U());
//  super_op Ubt = U_transform(U);	// Superoperator for basis change(s)
//  U = Ubt.get_mx();
//  matrix S2 = U*times_adjoint(get_mx(),U); 
//  put_mx(S2);				// Set new superoperator matrix
//  put_basis(LOp1.get_basis());		// Set new Hilbert space basis
  return;
  }

	// Input		LOp  : Superoperator (this).
	// 			Op   : General operator.
	// Return		none : General Operator Op put into
	//			       the Hilbert space basis of LOp

void super_op::LOp_base(const    gen_op& Op) const { Op.Op_base(Hbs); }
void super_op::SetHSBaseOf(const gen_op& Op) const { Op.Op_base(Hbs); }

// ____________________________________________________________________________
// G                   SUPER OPERATOR AUXILIARY FUNCTIONS
// ____________________________________________________________________________


int super_op::HS()   const { return HSp; }
int super_op::dim()  const { return LSp; }
int super_op::size() const { return LSp; }
int super_op::LS()   const { return LSp; }

void super_op::eigenvalues(int nc, int ri) const

	// Input		LOp   : Superoperator (this)
	//			nc    : Nubmer of columns
	//			ri    : Real/Imag/Complex flag
        // Return		None  : LOp eigenvalues sent to
	//				standard output
        // Note			      : Sets EBR of LOp

  {
  if(HSp) 				// Check for NULL LOp
    {
    set_EBR();				// Put LOp into its EBR
    int j=0;
    std::cout << "\n\tSuperoperator Eigenvalues:\n";
    for(int i=0; i< LSp; i++)
      {
      if(ri == 1)			// Print Complex
        std::cout << i << ". " << get(i,i) << "\t";
      else if(ri == -1)			// Print Imaginaries
        std::cout << i << ". " << Im(get(i,i)) << "\t";
      else				// Print Reals
        std::cout << i << ". " << Re(get(i,i)) << "\t";
      j++;
      if(j == nc)			// Output nc eigenvalues per line
        {
        std::cout << "\n";
        j = 0;
        }
      }
    return;
    }
  else
    std::cout << "\nNull superoperator, all zero eigenvalues";
  return;
  }

// ____________________________________________________________________________
//                  SUPER OPERATOR COMPONENT MANIPULATIONS
// ____________________________________________________________________________

// -------------------------- Matrix Manipulations ----------------------------

matrix super_op::Mx()     const { return (LSp)?mx:matrix(0,0); }
matrix super_op::get_mx() const { return (LSp)?mx:matrix(0,0); }

void super_op::put_mx(const matrix& mx1)

	// Input		LOp   : Superoperator (this)
        // 			mx1   : Matrix (in L. space) 
        // Output		none  : LOp altered to have mx as
	//				its representation

  {
  if(!LSp) {*this=super_op(mx); return;}// Null LOp, set LOp in DBR
  if(!checkLOp(mx1, Hbs, 1))		// Insure proper dimensions 
    {
    LOperror(5, "put_matrix", 1);	//   Problems in put_matrix
    LOpfatal(41);		// Can't mix LOp and matrix
    }
  mx = mx1;				// Replace LOp matrix
  }


// ----------------------- Basis Manipulations --------------------------

basis super_op::Bs()        const { return (LSp)?Hbs:basis(1); }
basis super_op::get_basis() const { return (LSp)?Hbs:basis(1); }

void super_op::put_basis(const basis& bs1)

	// Input		LOp   : Superoperator (this)
        // 			bs1   : Basis
        // Output		none  : LOp altered to have bs as
	//				its current basis without a
	//			        change in the matrix

  {
  if(LSp)			// Check for NULL LOp
    Hbs = bs1;			// Replace LOp basis
  else
    {
    *this = super_op();		// Return NULL LOp
    LOperror(7,1);		// Accessing Null LOp
    LOperror(77);		// Error accessing LOp basis
    }
  }


basis super_op::LBs()        const { return (LSp)?Lbs:basis(1); }
basis super_op::get_Lbasis() const { return (LSp)?Lbs:basis(1); }


void super_op::put_Lbasis(const basis& bs1)

	// Input		LOp   : Superoperator (this)
        // 			bs1   : Basis
        // Output		none  : LOp altered to have bs as
	//				its current Liouville basis
	//				without a change in the matrix

  {
  if(LSp)			// Check for NULL LOp
    Lbs = bs1;			// Replace LOp basis
  else
    {
    *this = super_op();		// Return NULL LOp
    LOperror(7,1);		// Accessing Null LOp
    LOperror(77);		// Error accessing LOp basis
    }
  }

// -------------------- Individual Element Manipulations ----------------------


complex super_op::operator() (int row, int col) const { return get(row,col); }

	// Input		LOp   : Superoperator (this)
        //                      row   : Row index
        //                      col   : Column index
        // Output               z     : Value of <row|LOp|col>
	// Note			      : LOp remains static, so LOp(i,j) = z
	//				will not work, unlike the matrix analog


void super_op::put(int row, int col, const complex& z)

	// Input		LOp   : Superoperator (this)
        //                      row   : Row index
        //                      col   : Column index
	//			z     : Complex number
        // Output               none  : Value of <row|LOp|col> set to z

  {
  if(!LSp)				// Check for NULL LOp
    {
    LOperror(7, 1);			// Dealing with Null operator
    LOperror(5, "put", 1);		// Error in the put funciton
    LOpfatal(67);		// Error accessing element
    }
  if(!checkLOp(row,col,1))		// Check element indices
    {
    LOperror(5, "put", 1);		// Error in the put function
    std::string ele = "<" + std::string(Gdec(row))
               + std::string("|LOp|")
               + std::string(Gdec(col))
               + std::string(">");
    LOperror(6, ele, 1);		// Error in element access
    LOpfatal(67);		// Error accessing element
    }
  mx.put(z,row,col);			// Set the element value
  }


complex super_op::get(int row, int col) const

	// Input		LOp   : Superoperator (this)
        //                      row   : Row index
        //                      col   : Column index
        // Output               z     : Value of <row|LOp|col>
	// Note			      : LOp remains unchanged, element copy 	

  {
  if(!LSp) return complex0;		// 0 for elements of NULL superop
  if(!checkLOp(row,col,1))		// Check element indices
    {
    LOperror(5, "get", 1);		// Error in the put function
    LOpfatal(67);		// Error accessing element
    }
  return (mx.get(row,col));		// Output the element value
  }


// ____________________________________________________________________________
// H                       CLASS SUPER OPERATOR CHECKS
// ____________________________________________________________________________


bool super_op::checkLOp(const matrix& mx1, const matrix& mx2, int warn) const
  {
  bool TF = checkLOp(mx1, 1) 		// Insure square matrix mx1
          & checkLOp(mx2, 1);		// Insure square matrix mx2
  if(TF) TF = (mx1.cols()==mx2.cols());	// Dimension equivalence
  if(TF) return true;
  if(warn)
    {
    if(warn > 1) LOpfatal(51);	// Dimensioning mismatch!
    else         LOperror(51,1);
    }
  return false;
  }


bool super_op::checkLOp(const matrix& mx, const basis& Hbs, int warn) const
  {
  int TF = checkLOp(mx, 1);		// Insure square matrix mx
  TF *= (mx.cols() == Hbs.dim_LS());	// Insure dimension equivalence
  if(TF) return true;
  if(warn)
    {
    if(warn > 1) LOpfatal(51);	// Dimensioning mismatch!
    else         LOperror(51,1);
    }
  return false;
  }

 
/* This set of functions insures that the LOp/Op/mx/bs is properly square &/or 
   of compatible dimension to mix with another LOp/Op/mx/bs.  All return TRUE
   if dimensions look OK.  The functions take a warning level flag which will
   be used to decide what happens under a FAIL.  If the flag is 0 then nothing 
   will result.  If the flags is 1 a non-fatal error is output.  If the flag
   is >1 a fatal error is output and program execution stopped.              */   

bool super_op::checkLOp(const super_op& LOp1, int warn) const
  {
  if(mx.cols() == LOp1.mx.cols()) return true;
  if(warn)
    {
    if(warn > 1) LOpfatal(52);	// LOp-LOp dimension mismatch
    else         LOperror(52,1);
    }
  return false;
  }


bool super_op::checkLOp(const gen_op& Op, int warn) const

	// Input		LOp  : super operator (this)
	// 			Op   : general operator
	// Output		T_F  : TRUE if LOp and Op are
	//			       compatible, FALSE if not
  {
  basis bsOp = Op.get_basis();		// Get Hilbert basis
//  int ncOp = bsOp.sub_N();		// Number of sub-spaces in Op
//  int ncLOp = Hbs.sub_N();		// Number of sub-spaces in LOp
  int ls = Op.dim_LS();			// Get Op Hilbert space dimension 
//  if ((mx.cols() == ls) && (ncOp == ncLOp)) return true;	// Compare Liouville spaces
  if (mx.cols() == ls) return true;	// Compare Liouville spaces
  if(warn)
    {
    if(warn > 1) LOpfatal(53);	// LOp-Op Dimensioning mismatch!
    else         LOperror(53,1);
    }
  return false;
  }


bool super_op::checkLOp(const matrix& mx, int warn) const
  {
  if(mx.cols() != mx.rows()) 			// Insure square matrix
    {						// If it isn't square
    if(warn)					// we might want to issue
      {						// a warning and even quit...
      if(warn > 1) LOpfatal(50);		// Rectangular array!
      else         LOperror(50,1);		// Rectangular array
      }
    return false;				// Return that test failed
    }
// The next insures ls is square of some hs. But this messes up
// on composite space stuff   
/*
  int hs = int(sqrt(mx.cols()));		// This is the Hilbert space
  int ls = hs*hs;				// Recalculate Liouville dim
  if(ls != mx.cols())				// The two should agree unless
    {						// the array has a "non-square"
    if(warn)					// Liouville space, its bad.
      {						// Issue warning or even quit...
      if(warn > 1) LOpfatal(49);		// No Hilbert space
      else         LOperror(49,1);		// Note: This will fail for some
      }						// composite Liouville spaces
    return false;				// Return that test failed
    }
*/
  return true;					// Test pass
  }


bool super_op::checkLOp(int row, int col, int warn) const

	// Input		LOp   : Superoperator (this)
        //                      row   : Row index
        //                      col   : Column index
        // Output               T_F   : TRUE if <row|LOp|col> exists

  {
  int TF = 1;					// Assume all is O.K.
  if(row<0 || row>=LSp) TF=0;			// Check row index
  if(col<0 || col>=LSp) TF=0;			// Check column index
  if(TF) return true;				// Return TRUE if element OK
  if(warn)					// Else issue warnings as
    {						// desired, perhaps fatal ones
    if(warn > 1) LOpfatal(56);		//   Element out of range
    else         LOperror(56,1);
    }
  return false;					// Problems with this element
  }


// ____________________________________________________________________________
// I                       CLASS SUPER OPERATOR TESTS
// ____________________________________________________________________________

///Center Test Functions

void super_op::status() const

	// Input		LOp   : Superoperator (this)
	// Output		void  : Outputs superoperator status

  {
  std::cout << "\n\tLOp matrix references = " << mx.refs();
  std::cout << "\n\t           type       = " << mx.mxtype();
  std::cout << "\n\tLOp basis  references = " << Lbs.refs();
  std::cout << "\n\t           type       = " << (Lbs.U()).mxtype();
  std::cout << "\n\tHS  basis  references = " << Hbs.refs();
  std::cout << "\n\t           type       = " << (Hbs.U()).mxtype();
//  if(mx.rows()) mx.status(1);
  std::cout << "\n";
  }


int super_op::operator == (const super_op& LOp1)

	// Input		LOp   : Superoperator (this)
        //                      LOp   : Superoperator
        // Output               T_F   : TRUE if LOp = LOp1

  {
  if(mx  != LOp1.mx)   return 0;		// Check that LOp reps equal
  if(Lbs != LOp1.Lbs) return 0;			// Check that bases are equal 
  if(Hbs != LOp1.Hbs) return 0;			// Check HS bases are equal
  return 1;
  }


int super_op::below(double d) const

	// Input		LOp   : Superoperator (this)
        //                      d     : Number
        // Output               T_F   : TRUE if all norm(<i|LOp|j>) <= d

  {
  int T_F = 1;
  for(int i=0; i<LSp && T_F; i++)
    for(int j=0; j<LSp && T_F; j++)
      if(norm(get(i,j)) > d)
        {
        T_F = 0;
        break;
        }
  return T_F;
  }


// ____________________________________________________________________________
// J                    CLASS OPERATOR I/O FUNCTIONS
// ____________________________________________________________________________

	// Input		LOp   : Superoperator (this)
	// 			ostr : Output stream
	//			full : Amount of output
	// Return		     : Superoperator output to stream

std::ostream& super_op::print(std::ostream& ostr, int full) const
  {
  if(!LSp) return ostr << "NULL Super Operator\n";
  ostr << "Matrix:\n" << mx
       << "Hilbert Space Basis:\n" << Hbs
       << "Liouville Space Basis:\n" << Lbs;
  if(full) ostr << "\n";
  return ostr;
  }

 
std::ostream& operator << (std::ostream& ostr, const super_op& LOp)

	// Input		ostr : Output stream
	// 			LOp  : Superoperator
	// Return		     : Superoperator output to stream

  { return LOp.print(ostr); }

// -------------------------- Binary Output Functions -------------------------
 
/*         Input                LOp  : Superoperator (this)
                                fn   : Output binary filename
				fp   : File stream (pointing at LOp spot)
           Return               void : LOp is written in binary format to
				       the file or filestream.               */

void super_op::write(const std::string& fn) const
  {
  std::ofstream fp;					// Construct a file
  fp.open(fn.c_str(), std::ios::out|std::ios::binary);	// Open file
  write(fp);						// Write LOp
  fp.close();						// Close file
  }

std::ofstream& super_op::write(std::ofstream& fp) const
  {
  int dim = HSp;			// Get the Hilbert space dimension
  fp.write((char*)&dim,sizeof(int));	// Write the Hilbert space dimension
  dim = LSp;				// Get the Liouville space dimension
  fp.write((char*)&dim,sizeof(int));	// Write the Liouville space dimension
  mx.write(fp,0);			// Write the Liouville space matrix
  Hbs.write(fp);			// Write the Hilbert space basis
  Lbs.write(fp);			// Write the Liouville space basis
  return fp;
  }

// --------------------------- Binary Input Functions -------------------------
 
/*         Input                LOp  : Superoperator (this)
                                fn   : Input binary filename
				fp   : File stream (pointing at LOp spot)
				Op   : Operator used to set working HS basis
				LOp  : Superoperator used to set working bases
           Return               void : LOp is read in binary format to
				       the file or filestream.
	   Note                      : Superoperators read in this fashion will
				       NOT share the same bases even if they
				       did so when written.  This must be done
                                       explicitly if desired.                */


void super_op::read(const std::string& fn)

	// Input		LOp  : Superoperator (this)
	// 			fn   : Filename
	// Return		void : Superoperator LOp is read from
	//			       a file called filename
	// Note			     : No basis sharing with other
	//			       Ops/LOps will be assumed here!

  {
  std::ifstream fp;					// Construct a file
  fp.open(fn.c_str(),std::ios::in|std::ios::binary);	// Open file
  read(fp);						// Read LOp
  fp.close();						// Close file
  }


void super_op::read(const std::string& fn, const gen_op& Op)

	// Input		LOp  : Superoperator (this)
	// 			fn   : Filename
	//			Op   : General operator
	// Return		void : Superoperator LOp is read from
	//			       a file called filename.  LOp
	//			       Hilbert space basis will be set
	//			       to the working basis of Op if
	//			       they match.

  {
  std::ifstream fp;					// Construct a file
  fp.open(fn.c_str(), std::ios::in|std::ios::binary);	// Open file
  read(fp);					// Read LOp
  if(Hbs == Op.get_basis())			// Compare Hilbert bases
    Hbs = Op.get_basis();
  fp.close();					// Close file
  }


void super_op::read(const std::string& fn, const super_op& LOp1)

	// Input		LOp  : Superoperator (this)
	// 			fn   : Filename
	//			LOp1 : 2nd superoperator
	// Return		void : Superoperator LOp is read from
	//			       a file called filename.  LOp
	//			       bases will be set to the bases
	//			       of LOp1 if they match.

  {
  std::ifstream fp;					// Construct a file
  fp.open(fn.c_str(),std::ios::in);			// Open file
  read(fp);					// Read LOp
  if(Hbs == LOp1.Hbs) Hbs = LOp1.Hbs;
  if(Lbs == LOp1.Lbs) Lbs = LOp1.Lbs;
  fp.close();					// Close file
  }


std::ifstream& super_op::read(std::ifstream& fp)

	// Input		LOp  : Superoperator (this)
	// 			fp   : File (pointing at LOp spot)
	// Return		void : Superoperator LOp is read from
	//			       file fp at current location
	// Note			     : No basis sharing with other
	//			       Ops/LOps will be assumed here!

  {
  int dim;
  matrix mxtmp;
  fp.read((char*)&dim,sizeof(int));	// Read the Hilbert space dimension
  HSp = dim;			// Set the Hilbert space dimension
  fp.read((char*)&dim,sizeof(int));	// Read the Liouville space dimension
  LSp = dim;			// Set the Liouville space dimension
  mx.read(fp);			// Read the Liouville space matrix
  mxtmp.read(fp);		// Read the Hilbert space basis
  Hbs = mxtmp;			// Set the Hilbert space basis
  mxtmp.read(fp);		// Read the Liouville space basis
  Lbs = mxtmp;			// Set the Liouville space basis
  return fp;
  }


// ____________________________________________________________________________
// L                 CLASS SUPEROPERATOR LEFTOVER FUNCTIONS
// ____________________________________________________________________________
 

super_op HsuperX(const gen_op& Heff)

	// Input		Heff  : Effective Hamiltonian (Hz)
	// Output		LOp   : Hamiltonian commutation superoperator
	//
	//       <a,aa|LOp|b,bb> = del  * del     <b|H   |b> - <bb|H   |bb>
	//			      ab     aa,bb    eff           eff
	//
	// Note			      :	LOp is returned in angular frequency
	//				units

  {
  basis bs = Heff.get_basis();			// Get Hilbert space basis
  int nc = bs.sub_N();                          // Number of sub-spaces
  if(nc > 1) return Hsuper(Heff);
  const double pi2 = 6.283185307;
  int hs = Heff.size();				// Get Hilbert space size
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
            LOp.put(aaa,bbb,complex(pi2*wbbb));	// Set matrix element (rad/sec)
            }
          bbb++;
          }
        }
      aaa++;
      }
  return LOp;
  }


super_op Hsuper(const gen_op& Heff)

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

#endif 						// SuperOp.cc
