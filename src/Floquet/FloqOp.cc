/* FloqOp.cc ****************************************************-*-c++-*-
**									**
**	                     G A M M A					**
**						 			**
**	Floquet Operators                  Implementation 		**
**								 	**
**	Copyright (c) 1992  	        			 	**
**	                  					 	**
**	Eidgenoessische Technische Hochschule			 	**
**	Labor fuer physikalische Chemie		 			**
**	8092 Zurich / Switzerland		 			**
**						 			**
**      Scott A. Smith                                                  **
**      Copyright (c) 1994, 1997                                        **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** Class floq_op embodies Floquet operators in GAMMA.  The class sets	**
** defines the properties and allowed operations such operators such	**
** that they may be freely used in GAMMA based programs.		**
**                                                                      **
*************************************************************************/

#ifndef   Floq_op_cc_			// Is file already included?
#  define Floq_op_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Floquet/FloqOp.h>		// Include the header file
#include <Basics/Gconstants.h>		// Include Pi and other values
#include <Basics/Gutils.h>		// Include GAMMA error handling
#include <stdlib.h>

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                 CLASS FLOQUET OPERATOR ERROR HANDLING
// ____________________________________________________________________________


 
        // Input		FOp	: Floquet operator (this)
        //                      eidx    : Error index
        //                      noret   : Flag for return (0=linefeed)
        // Output               none    : Error message

void floq_op::FOperror(int eidx, int noret) const
  {
  std::string hdr("Floquet Operator");
  std::string msg;
  switch(eidx)
    {
    case 0: GAMMAerror(hdr, 0, noret); break;   // Program Aborting        (0)
    case 1: GAMMAerror(hdr, 1, noret); break;   // Input File Stream Bad   (1)
    case 2: GAMMAerror(hdr, 2, noret); break;   // Output File Stream Bad  (2)
    case 6: msg = "Fourier Expansion Frequency Mismatch!";	
            GAMMAerror(hdr, msg, noret); break;	//                         (6)
    case 7: msg = "Floquet Space Mismatch Between Operators!";	
            GAMMAerror(hdr, msg, noret); break;	//                         (7)
    case 8: msg = "Photon Space Mismatch Between Operators!";	
            GAMMAerror(hdr, msg, noret); break;	//                         (8)
    case 9: GAMMAerror(hdr, 9, noret); break;   // Construction Error      (9)
    case 10:msg = "Inconsistent Floquet Space Dimensions";	
            GAMMAerror(hdr, msg, noret); break;	//                         (10)
    case 11: msg = "Floquet Space Mismatch Between Operator & Matrix!";	
            GAMMAerror(hdr, msg, noret); break;	//                         (11)
    case 12: msg = "Error Accessing Operator Sub-Matrix!";	
            GAMMAerror(hdr, msg, noret); break;	//                         (12)
    case 20:msg = "Dimensioning Problems During Addition!";	
            GAMMAerror(hdr, msg, noret); break;	//                         (20)
    case 21:msg = "Dimensioning Problems During Subtraction!";	
            GAMMAerror(hdr, msg, noret); break;	//                         (21)
    case 22:msg = "Dimensioning Problems During Multiplication!";	
            GAMMAerror(hdr, msg, noret); break;	//                         (22)
    case 41:msg = "Specified Position Exceeds Operator Dimension!";	
            GAMMAerror(hdr, msg, noret); break;	//                         (41)
    case 42:msg = "Inconsistent Dimension of Input Hilbert Space Operator!";	
            GAMMAerror(hdr, msg, noret); break;	//                         (42)
    case 43:msg = "Side-Diagonal # Inconsistent With Floquet Dimension!";	
            GAMMAerror(hdr, msg, noret); break;	//                         (42)
    case 49:msg = "Sub-Block Access Out Of Range!";	
            GAMMAerror(hdr, msg, noret); break;	//                         (49)
    case 50:msg = "Element Access Out Of Range!";	
            GAMMAerror(hdr, msg, noret); break;	//                         (50)
    case 51:msg = "Error Accessing Internal Floquet Component";	
            GAMMAerror(hdr, msg, noret); break;	//                         (51)
    default:GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error   (-1)
    }
  }
     

volatile void floq_op::FOpfatality(int eidx) const

        // Input		FOp	: Floquet operator (this)
        //                      eidx    : Error index
        // Output               none    : Error message output
        //                                Program execution stopped
 
  {
  FOperror(eidx, 1);			// First output the error
  if(eidx) FOperror(0);			// Now output it fatal
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii            FLOQUET OPERATOR COMPATIBILITY CHECKING
// ____________________________________________________________________________


int floq_op::FOpCheck(const floq_op& FOp1, int warn) const
 
        // Input		FOp	: Floquet operator (this)
	//			FOp1	: A second floquet operator
	//			warn    : Warning level
	//				      0 = No warnings (default)
	//				      1 = Non-fatal warnings
	//				     >1 = Fatal warnings
        // Output               TF      : True if the two operators match in
	//				  three areas -
        //                                     1. Phonon space sizes (N)
        //                                     2. Hilbert space sizes (hs)
        //                                     3. Fourier frequencies (Om)

  {
  if(Om != FOp1.Om)				// Check Fourier frequencies
    {
    if(!warn)       return 0; 
    else if(warn>1) FOperror(6, 1);
    else            FOpfatality(6);
    }
  if(hs != FOp1.hs)				// Check Floquet dimensions
    {
    if(!warn)       return 0; 
    else if(warn>1) FOperror(7, 1);
    else            FOpfatality(7);
    }
  if(N != FOp1.N)				// Check photon space size
    {
    if(!warn)       return 0;
    else if(warn>1) FOperror(8, 1);
    else            FOpfatality(8);
    }
  return 1;
  }


int floq_op::FOpCheck(const matrix& Fmx, int warn) const
 
        // Input		FOp	: Floquet operator (this)
	//			Fmx	: An array in Floquet space
	//			warn    : Warning level
	//				      0 = No warnings (default)
	//				      1 = Non-fatal warnings
	//				     >1 = Fatal warnings
        // Output               TF      : True if the operator & matrix match
        //                                Floquet space dimensions

  {
  int Fdim = (2*N+1)*hs;			// Floquet space dimension
  if(Fdim!=Fmx.rows() || Fdim==Fmx.cols()) 	// Check Floquet space size
    {
    if(!warn)       return 0;
    else if(warn>1) FOperror(11, 1);
    else            FOpfatality(11);
    }
  return 1;
  }


int floq_op::FOpCheck(int N1, int N2, int H1, int H2, int warn) const
 
        // Input		FOp	: Floquet operator (this)
	//			N1,N2   : Photon indices
	//			H1,H2	: Operator indices (Hilbert space)
	//			warn    : Warning level
	//				      0 = No warnings (default)
	//				      1 = Non-fatal warnings
	//				     >1 = Fatal warnings
        // Output               TF      : True if the element <I|FOp|J> is
        //                                within Floquet space dimensions where
	//				  I=(N+N1)*hs+H1 & J=(N+N2)*hs+H2

  {
  int TF = 1;					// Assume element is fine
  if(abs(N1)>N  || abs(N2)>N) TF=0;		// Bad if outside photon range
  if(abs(H1)>hs || abs(H2)>hs) TF*=0;		// Bad if outside Hilbert range
  if(!TF) 					// If element bad, here's  how
    {						// we handle it
    if(!warn)       return 0;
    else if(warn>1) FOperror(51, 1);		// Element access out of range
    else            FOpfatality(51);
    }
  return 1;
  }

// ____________________________________________________________________________
// iii            FLOQUET OPERATOR - HILBERT OPERATOR BASES
// ____________________________________________________________________________

void floq_op::SetBasis(const gen_op& Op)

        // Input		FOp	: Floquet operator (this)
	//			Op	: A Hilbert space opertor
	// Output		void	: FOp has its basis set to the
	//				  tensor product E X OpBs
	// Note				: PRIVATE - no FOp/Op compatibility
	//				  checks are made herein

  {
  int Fdim = 2*N+1;				// Full Floquet space dimension
  matrix E(Fdim, Fdim, i_matrix_type);		// I matrix, full Floquet space
  basis OpBs = Op.get_basis();			// Get Op basis (Hilbert space)
  matrix FBsMx = tensor_product(E,OpBs.U());	// Form tensor product E X OpBs
  basis Fbs(FBsMx);				// Set this now as a basis
  Op_base(Fbs);					// Put FOp into basis Fbs
  }


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A             CLASS FLOQUET OPERATOR CONSTRUCTORS/DETRUCTOR
// ____________________________________________________________________________ 


floq_op::floq_op() : gen_op(), Om(0), N(0), hs(0) {}

        // Input			: none
        //			FOp	: NULL Floquet operator (this)


floq_op::floq_op(const floq_op &FOp1)

        // Input                   FOp	: Floquet operator (this)
        //                         FOp1 : Floquet operator
        // Output			: Returns FOp equivalent to FOp1

  {
  N  = FOp1.N;					// Equate photon spaces 
  hs = FOp1.hs;					// Equate Floquet spaces
  Om = FOp1.Om;					// Equate Fourier frequencies
  gen_op::operator=((const gen_op&) FOp1);	// Equate the operators
  }  


 
        // Input                   NP   : Integer for photon space
        //                         HS   : Dimension of spin Hilbert space
        //                         Omega: Fourier (rotation) frequency (Hz)
        // Output                  FOp  : Returns a new FlOp of dimension
        //                                2*N+1*hs and k*wr(k=-N...N) on main
        //                                diagonal.

floq_op::floq_op(int NP, int HS, double Omega)
  {
  int Fdim = (2*NP+1)*HS;			// Floquet space dimension
  N  = NP;					// Set photon space 
  hs = HS;					// Set spin Hilbert space
  Om = Omega;					// Set Fourier expansion freq.
  matrix Fmx(Fdim,Fdim,complex0,n_matrix_type);	// Op. Matrix (Floquet space)
  gen_op::operator=(Fmx);			// Set the operator represent.
  }


floq_op::floq_op(int NP, int HS, double Omega, const matrix& mx)
        :gen_op(mx),Om(Omega),N(NP),hs(HS) 
       
        // Input                   FOp	: Floquet operator (this)
        //                         NP	: Integer for photon space
        //                         HS   : Dimension of spin Hilbert space
        //                         Omega: Fourier (rotation) frequency (Hz)
	//			   mx   : Matrix for FOp representation 
        // Output                   	: Returns a new FOp with dimension
        //                                2*N+1*hs and k*wr(k=-N...N) on main
        //                                diagonal.
	// Note				: The input matrix MUST properly
	//				  relate to the specified photon and 
	//				  spin Hilbert space dimensions!
	//				  mx.rows()=mx.cols()=(2NP+1)HS=Fdim

  {
  int Fdim = (2*NP+1)*HS;		// Floquet space dimension
  if(mx.cols() != Fdim)			// Insure proper matrix size!
    {
    FOperror(9,1); 			// Error during construction
    FOpfatality(10); 			// Inconsistent Dimensions
    }					// --> mx is forced to be square when
  }					// we set the operator gen_op(mx)!
 


floq_op::floq_op(int NP,int HS,double Omega, const matrix& mx, const basis& bs)
        :gen_op(mx,bs),Om(Omega),N(NP),hs(HS)
 
	// Input                NP	: Floquet space dimension
        //			HS	: Hilbert space dimension
        //			Omega	: Fourier (rotation) frequency (Hz)
	//			mx	: Matrix for FOp representation 
	//			bs	: Basis for FOp representation 
        // Output                   	: Returns a new FOp with dimension
        //                                2*N+1*hs and k*wr(k=-N...N) on main
        //                                diagonal.
	// Note				: The input {mx,bs} MUST properly
	//				  relate to the specified photon and 
	//				  spin Hilbert space dimensions!
	//				  mx.rows()=mx.cols()=(2NP+1)HS=Fdim
        // Note				: Assumes mx DBR of Op

  {
  int Fdim = (2*NP+1)*HS;		// Floquet space dimension
  if(mx.cols() != Fdim)			// Insure proper {mx,bs} size!
    {
    FOperror(9,1); 			// Error during construction
    FOpfatality(10); 			// Inconsistent Dimensions
    }					// --> mx is forced to be square when
  }					// --> & {mx,bs} compatibility also

	// Input		Op1  : Floquet operator.
	// 			Op   : Floquet operator (this).
	// Return		Op   : Operator which is a copy of
	//			       the input operator, Op = Op1.
	// Note		             : Result EXCLUSIVELY in WBR of Op1

floq_op& floq_op::operator = (const floq_op& FOp1)
  {
  if(this == &FOp1) return *this;		// Do nothing if already equal
  N  = FOp1.N;					// Equate photon spaces
  hs = FOp1.hs;					// Equate Floquet spaces
  Om = FOp1.Om;					// Equate Fourier frequencies
  gen_op::operator=((const gen_op&)FOp1);	// Equate operators
  return *this;
  }

floq_op::~floq_op () {}

// ____________________________________________________________________________
// B                   FLOQUET OPERATOR ACCESS FUNCTIONS
// ____________________________________________________________________________
 
 
int floq_op::hsdim() const { return hs; }
 
        // Input                Op   : Floquet operator (this)
        // Return               int  : Op Hilbert space dimension
	// Note			     : This is the base dimension from
	//			       which Floquet space is constructed
	//					FS = (2*NP+1)*HS
 

int floq_op::phodim() const { return N; }
 
        // Input                Op   : Floquet operator (this)
        // Return               int  : Op Photon space dimension
	// Note			     : Used to build up Floquet space as
	//					FS = (2*NP+1)*HS
 

double floq_op::omega() const { return Om; }
 
        // Input                Op   : Floquet operator (this)
        // Return             double : Fourier expansion frequency
	// Note			     : Usually set in Hz


int floq_op::dim()  const { return (gen_op::dim()); }
int floq_op::size() const { return (gen_op::dim()); }
 
        // Input                Op   : Floquet operator (this)
        // Return               int  : Op Hilbert space dimension
	// Note			     : This is FULL Floquet space dimension
	//					FS = (2*NP+1)*HS
	//			       which is larger than the Hilbert space
	//			       dimension (unless NP = 0)

// ____________________________________________________________________________
// C   FLOQUET OPERATOR FUNCTIONS, FLOQUET OPERATOR WITH FLOQUET OPERATOR 
// ____________________________________________________________________________

/* These functions handle all dealings between Floquet operators and Floquet
   operators.

 Function   Arguments      Return       Function   Arguments      Return
 --------   ---------   -------------   --------   ---------   -------------
    +       FOp1,FOp2   FOp=FOp1+FOp2      *       FOp1,FOp2   FOp=FOp1*FOp2
    +=      FOp1        FOp=FOp+FOp1       *=      FOp1        FOp=FOp*FOp1
    -       FOp1,FOp2   FOp=FOp1-FOp2      &=      FOp1        FOp=FOp1*FOp 
    -=      FOp1        FOp=FOp-FOp1
    -       FOp1        FOp=-FOp1                                            */

floq_op floq_op::operator + (const floq_op& FOp2) const
  {
  floq_op FOp(*this);			// Set FOp to be FOp1
  FOp += FOp2;				// Use operator overload
  return FOp;				// Return FOp1+FOp2
  }

void floq_op::operator += (const floq_op& FOp1)
  {
  if(!FOpCheck(FOp1,1)) FOpfatality(20);// Insure dimensions match
  gen_op::operator+=((gen_op&)FOp1);	// Just add the operators
  }

floq_op floq_op::operator - (const floq_op &FOp2) const
  { 
  floq_op FOp(*this);			// Set FOp to be FOp1
  FOp -= FOp2;				// Use operator overload
  return FOp;				// Return FOp1-FOp2
  }

floq_op floq_op::operator - () const
  {
  floq_op FOp(N,hs,Om);
  FOp.gen_op::operator=(-(gen_op&)(*this));
  return FOp; 
  }

void floq_op::operator -= (const floq_op &FOp1)
  {
  if(!FOpCheck(FOp1,1)) FOpfatality(21);// Insure dimensions match
  gen_op::operator-=((gen_op&) FOp1);	// Just subtract the operators
  }

floq_op floq_op::operator * (const floq_op &FOp2) const
  {
  floq_op FOp(*this);			// Set FOp to be FOp1
  FOp *= FOp2;				// Use operator overload
  return FOp;				// Return FOp1*FOp2
  }

void floq_op::operator *= (const floq_op& FOp1)
  {
  if(!FOpCheck(FOp1,1)) FOpfatality(22);// Insure dimensions match
  gen_op::operator*=((gen_op&) FOp1);	// Just multiply the operators
  }

void floq_op::operator &= (const floq_op& FOp1)
  {
  if(!FOpCheck(FOp1,1)) FOpfatality(22);// Insure dimensions match
  gen_op::operator&=((gen_op&) FOp1);	// Use genop for this
  }

// ____________________________________________________________________________
// D        FLOQUET OPERATOR FUNCTIONS, FLOQUET OPERATOR WITH MATRIX
// ____________________________________________________________________________

/* These functions handle all dealings between Floquet operators and Floquet
   dimension matrices.

 Function   Arguments      Return       Function   Arguments      Return
 --------   ---------   -------------   --------   ---------   -------------
    +       FOp1,Fmx    FOp=FOp1+Fmx       *       FOp1,FMx    FOp=FOp1*FMx
    +       Fmx,FOp1    FOp=FOp1+Fmx       *       FMx,FOp1    FOp=FMx*FOp1
    +=      FOp1        FOp=FOp+FOp1       *=      FMx         FOp=FOp*FMx
    -       FOp1,Fmx    FOp=FOp1-Fmx       &=      FMx         FOp=FMx*FOp 
    -=      Fmx         FOp=FOp-Fmx        -       Fmx,FOp1    FOp=Fmx-FOp1  */

/*****************************************************************************
**                                                                          **
** Note: All of these functions except for operator + (FOp, Fmx) checked    **
**       that the dimension of Fmx would match N*hs of FOp rather than      **
**       the Floquet space of (2*N+1)*hs.  The latter was enforced for all  **
**       FOp-FOp operations as well as the FOp+Fmx function below.  I don't **
**       see how this is possible!  So either these functions DEMAND that   **
**       the array dimensions match the Floquet dimensions or Marc had      **
**       something in mind here that I am not grasping...... in which case  **
**       it looks like the functions will ALWAYS FAIL the way they were	    **
**       written.  ----> I am switching all of these to use (2*N+1)*hs as   **
**       the enforced dimension.  Things will become clearer when I write   **
**       some test programs.  I'll eventually talk to Marc about it too.    **
**                                                                          **
**       sosi circa Feb 1999  GAMMA version 4.0                             **
**                                                                          **
** P.S.  I've made all arguments contant here because all these functions   **
**       must work in the default basis which every FOp must have. However  **
**       this can only work if, when the input FOp (FOp1) is copied it must **
**       copy the DBR and be able to access it.  Is that true? Gotta check. **
**                                                                          **
*****************************************************************************/

floq_op floq_op::operator + (const matrix& Fmx) const
  {
  floq_op FOp(*this);			// Set FOp to be FOp1
  FOp += Fmx;				// Use operator overload
  return FOp;				// Return FOp1+Fmx
  }

floq_op operator + (const matrix& Fmx, const floq_op &FOp1)
  {
  floq_op FOp(FOp1);			// Set FOp to be FOp1
  FOp += Fmx;				// Use operator overload
  return FOp;				// Return FOp1+Fmx
  }

void floq_op::operator += (const matrix& Fmx)
  {
  if(!FOpCheck(Fmx,1)) FOpfatality(20);	// Insure dimensions match
  gen_op::operator+=((matrix&)Fmx);	// Just add Fmx to FOp
  }

floq_op floq_op::operator - (const matrix& Fmx) const
  {
  floq_op FOp(*this);			// Set FOp to be FOp1
  FOp -= Fmx;				// Use operator overload
  return FOp;				// Return FOp1-Fmx
  }

floq_op operator - (const matrix& Fmx, const floq_op& FOp1)
  {
  if(!FOp1.FOpCheck(Fmx,1)) 			// Insure dimensions match
    FOp1.FOpfatality(21);			// (needed to get proper error)
  floq_op FOp(FOp1.N,FOp1.hs,FOp1.Om,Fmx);	// Set FOp to be Fmx
  FOp -= FOp1.get_mx();				// Use operator overload
  return FOp;					// Return Fmx-FOp1
  }

void floq_op::operator -= (const matrix& Fmx)
  {
  if(!FOpCheck(Fmx,1)) FOpfatality(21);	// Insure dimensions match
  gen_op::operator+=(Fmx);		// Just subtract Fmx from FOp
  }

floq_op floq_op::operator * (const matrix& Fmx) const
  {
  floq_op FOp(*this);			// Set FOp to be FOp1
  FOp *= Fmx;				// Use operator overload
  return FOp;				// Return FOp1*Fmx
  }

floq_op operator * (const matrix& Fmx, const floq_op& FOp1)
  {
  floq_op FOp(FOp1);			// Set FOp to be FOp1
  FOp &= Fmx;				// Use operator overload
  return FOp;				// Return FOp1*Fmx
  }

void floq_op::operator *= (const matrix& Fmx)
  {
  if(!FOpCheck(Fmx,1)) FOpfatality(22);	// Insure dimensions match
  gen_op::operator*=(Fmx);		// Just multiply FOp into Fmx
  }
	
void floq_op::operator &= (const matrix& Fmx)
  {
  if(!FOpCheck(Fmx,1)) FOpfatality(22);	// Insure dimensions match
  gen_op::operator&=(Fmx);		// Just multiply Fmx into FOp
  }


// ____________________________________________________________________________
// E     FLOQUET  OPERATOR FUNCTIONS, FLOQUET OPERATOR WITH SCALAR
// ____________________________________________________________________________

/* These functions handle all dealings between Floquet operators and scalar
   values (int, double, complex).

 Function   Arguments      Return       Function   Arguments      Return
 --------   ---------   ------------    --------   ---------  ----------------
    *        FOp1,z     FOp = z*FOp1       *=       this, z   void, z*this
    *        z,FOp1     FOp = z*FOp1       *=       this, d   void, d*this
    *        FOp1,d     FOp = d*FOp1        /       FOp1, z   FOp = (1/d)*this
    *        d,FOp1     FOp = d*FOp1       /=       this, z   void, (1/d)*this
                                                                             */

floq_op floq_op::operator * (const complex& z) const
  { floq_op FOp(*this); FOp *= z; return FOp; }

floq_op operator * (const complex &z, const floq_op &FOp1)
  { floq_op FOp(FOp1); FOp *= z; return FOp; }

floq_op floq_op::operator * (double d) const
  { floq_op FOp(*this); FOp *= d; return FOp; }

floq_op operator * (double d, const floq_op &FOp1)
  { floq_op FOp(FOp1); FOp *= d; return FOp; }

void floq_op::operator *= (const complex &z) { gen_op::operator*=(z); }
void floq_op::operator *= (double d)         { gen_op::operator*=(d); }

floq_op floq_op::operator / (const complex &z) const
  { floq_op FOp(*this); FOp /= z; return FOp; }

floq_op floq_op::operator / (double d) const
  { floq_op FOp(*this); FOp /= d; return FOp; }

void floq_op::operator /= (const complex &z) { gen_op::operator/=(z); }
void floq_op::operator /= (double d)         { gen_op::operator/=(d); }

// ___________________________________________________________________________
// F                  COMPLEX  FLOQUET OPERATOR FUNCTIONS
// ___________________________________________________________________________

  

        // Input               FOp	: Floquet operator
        // Output          	Op	: General operator containing
	//			          trace(s) over photon space

gen_op pho_trace(const floq_op &FOp)
  {
  gen_op Op(matrix(FOp.hs, FOp.hs));		// The return operator
  complex sum;					// Element (for trace sum)
  for(int s1=0; s1<FOp.hs; s1++)		// Loop Floquet dimension
    for(int s2=0; s2<FOp.hs; s2++)		// i.e. over each element
      {
      sum = complex0;				//   Set element to 0
      for(int pho=-FOp.N; pho<= FOp.N; pho++)	//   Loop photon dimension
         sum += FOp.get(pho,pho,s1,s2);		//   Sum from all N operators
      Op.put(sum,s1,s2);			//   Store the Op element 
      }
  return Op;
  }        

        // Input		FOp	: Floquet operator
        // Return		ExpFOp	: Exponential of FOp
	//				   ExpFOp = exp(FOp)
        // Note				: Computed in EBR of FOp

floq_op exp(const floq_op& FOp)
  {
  floq_op ExpFOp(FOp.N, FOp.hs, FOp.Om);	// Allocate return operator
  ExpFOp.gen_op::operator= (FOp.gen_op::exp());	// Set exponential array
  return ExpFOp;				// Return exponentiated FOp
  }


floq_op floq_op::exp()

        // Input                FlOp    : Floquet operator (*this)
        // Return               ExpFlOp : Exponential of Op1
        //                                   Op = exp(Op1)
        // Note                         : Computed in EBR of Op1
	// sosi				: Be nice to set this up as in
	//				  class gen_op.....

  {
  floq_op ExpFOp(N, hs, Om);	// Allocate return operator
  ExpFOp=*this;
  ExpFOp.gen_op::operator= (ExpFOp.gen_op::exp());	// Set exponential array
  return ExpFOp;
  }


floq_op prop(floq_op& FH, double time)

        // Input                FH	: Floquet Hamiltonian operator
        //                      time	: Evolution time
        // Return		UFOp	: Floquet Propagator in EBR
	
  {     
  floq_op UFOp(FH.N, FH.hs, FH.Om);		// Allocate return operator
  FH.set_EBR();					// Set FH into an eigenbasis
  complex z(0,-2*PI*time);			// Exponential factor
  UFOp.gen_op::operator= (FH.gen_op::exp(z));	// Set exponential array
  return UFOp;					// Return propagator
  }


floq_op fprop(floq_op& FH , double time)

        // Input                FH	: Floquet Hamiltonian operator
        //                      time	: Evolution time (sec)
        // Return		UF	: Floquet Propagator in DBR
	
  {
  floq_op UFOp(FH);			// Begin with Hamiltonian copy
  UFOp.add_omega(); 			// Set to 'full' Floquet Hamiltonian
  UFOp.set_EBR();			// Set this into its eigenbasis
  complex z(0,-2*PI*time);		// Exponential factor
  UFOp = UFOp*z; 			// This is now i*FH*t (unitless)
  UFOp = exp(UFOp);        		// Take the exponential 
  UFOp.set_DBR();			// Put into default basis
  return UFOp;
  }

// ___________________________________________________________________________
// G               FLOQUET OPERATOR COMPONENT MANIPULATIONS
// ___________________________________________________________________________
   
// --------------------- Floquet Submatrix Manipulations ---------------------
        
        // Input		FOp	: Floquet operator (this)
        //			N1	: Row index of photon dimension
        //			N2	: Column index of photon dimension
        // Output		Op	: Operator at <N1|FOp|N2>

gen_op floq_op::get_block  (int N1, int N2) const { return (*this)(N1,N2); }
gen_op floq_op::operator() (int N1, int N2) const
  {
  if((abs(N1)>N) || (abs(N2)>N))
    {  
    FOperror(12, 1); 			// Error accessing Floquet submatrix
    FOpfatality(49);			// Element access out of range
    } 
  gen_op Op(matrix(hs,hs));		// Here is the return operator
  int Io = (N1+N)*hs;			// Full space row starting index
  int Jo = (N2+N)*hs;			// Full space column starting index
  for(int I=Io,i=0; i<hs; I++,i++)	// Copy the requested elements into
    for(int J=Jo,j=0; j<hs; J++,j++)	// the subspace operator elements
      Op.put(gen_op::get(I,J), i, j);
  return Op;
  }


       
        // Input		FOp	: Floquet operator (this)
        // 			Op	: General operator
        //			N1,N2	: Position <N2|Op|N1> in Floquet Matrix 
        // Output		void	: Floquet operator with gen_op Op
        //				  set at photon position <N2|Op|N1>

void floq_op::put_block(const gen_op &Op, int N1, int N2)
  {
  if(Op.dim() != hs)			// Insure Op dimension proper 
      {
      FOperror(12, 1); 			// Error accessing Floquet submatrix
      FOpfatality(42);			// Improper Op dimension
      }
  if((abs(N1)>N) || (abs(N2)>N))	// Insure put within FOp limits
      {
      FOperror(11, 1); 			// Error accessing Floquet submatrix
      FOpfatality(49);			// Floquet dimension exeeded
      }
  SetBasis(Op);				// Put FOp into Op basis
  matrix FMx  = get_mx();		// Get the FOp matrix
  matrix OpMx = Op.get_mx();		// Get the Op matrix
  int I =(N+N1)*hs;			// Block row start, full Floquet space
  int J =(N+N2)*hs;			// Block col start, full Floquet space
  FMx.put_block(I, J, OpMx);		// Put block on desired position
  put_mx(FMx);
  }


void floq_op::put_sdiag(const gen_op& Op, int sdn)
         
        // Input                FOp     : Floquet operator (this)
        //                       Op     : General operator
        //                      sdn_    : Number of side diagonals to be set
        // Output                FOp1   : Floquet operator with side diagonals
        //                              : set to Op.  The number of side
        //                                diagonals set will be sdn_

  {     
  if(Op.dim() != hs)			// Insure Op dimension proper 
    {
    FOperror(12, 1); 			// Error accessing Floquet submatrix
    FOpfatality(42);			// Improper Op dimension
    }
  if(abs(sdn) > 2*N)			// Insure side-diagonal number is NOT
    {					// beyond photon dimension
    FOperror(12, 1); 			// Error accessing Floquet submatrix
    FOpfatality(43);			// Side-diagonal > photon dimension
    }
  SetBasis(Op);				// Put FOp into Op basis
  matrix FMx  = get_mx();		// Get the FOp matrix
  matrix OpMx = Op.get_mx();		// Get the Op matrix
  int lowpos  = sdn;
  int highpos = 2*N;			// Last photon space position
  if(sdn<0)				// Force the diagonal specified
    { 					// on left side of FOp diag.
    lowpos=0;
    highpos=2*N+sdn;
    }
  int d1,d2;
  for(d1=lowpos; d1 <= highpos; d1++)	// Loop FOp block posiitons
    {					// and insert Op into FOp
    d2=d1-sdn;				// at each spot
    FMx.put_block(d2*hs,d1*hs,OpMx);
    }
  put_mx(FMx);				// Reset the FOp matrix
  }


// -------------------- Individual Element Manipulations ----------------------


void floq_op::put(const complex &z, int N1, int N2, int H1, int H2)
     
        // Input                FOp     : Floquet operator (this)
	//			N1,N2	: Indices of photon space:-N...N
	//			H1,H2	: Indices of Hilbert space:0...hs
	// Output		none	: Value of <N2,H2|Op|H1,N1> will be
	//				  set to z

  {
  if(!FOpCheck(N1,N2,H1,H2,1)) 		// Insure element in FOp range
    FOpfatality(50);			// Quit if not
  int I=(N+N1)*hs+H1;			// Floquet space row
  int J=(N+N2)*hs+H2;			// Floquet space column
  gen_op::put(z,I,J);			// Set the element
  }


complex floq_op::get(int N1, int N2, int H1, int H2 ) const

        // Input                FOp     : Floquet operator (this)
	//                   N1,N2    : Indices of photon space:-N...N
	//                   H1,H2    : Indices of Hilbert space:0...hs
	// Output                z    : complex value of <N2,H2|Op|H1,N1>
	// Note                       : floq_op remains unchanged

  {
  if(!FOpCheck(N1,N2,H1,H2,1)) 		// Insure element in FOp range
    FOpfatality(50);			// Quit if not
  int I=(N+N1)*hs+H1;			// Floquet space row
  int J=(N+N2)*hs+H2;			// Floquet space column
  return gen_op::get(I,J);		// Return the element
  }


void floq_op::put(const complex &z,int row,int col) {gen_op::put(z,row,col);}
       
       // Input                      : Floq_operator (this)
       //                 row,col    : Indices 
       // Output            none     : Value of <row|Op|col> in WBR
       //                            : set to z
       // Note                       : floq_op will be changed

complex floq_op::get(int row, int col) const { return gen_op::get(row,col); }

       // Input                      : Floq_operator (this)
       //                 row,col    : Indices 
       // Output            none     : Value of <row|Op|col> in WBR
       // Note                       : floq_op remains unchanged

// ___________________________________________________________________________
// H                OPERATOR REPRESENTATION MANIPULATIONS
// ___________________________________________________________________________

	// Input		Op    : Floquet operator (this)
	// Output set_DBR	none  : Op WBR set to DBR 
	// Output set_EBR	none : Op WBR set to EBR 
	
void floq_op::set_DBR() const { gen_op::set_DBR(); }
void floq_op::set_EBR() const { gen_op::set_EBR(); }

// ____________________________________________________________________________
// I                 CLASS FLOQUET OPERATOR I/O FUNCTION
// ____________________________________________________________________________

	// Input		FOp	: Floquet operator (this)
	//			ostr    : An output stream
	//			full	: Output level
	//			spacer	: Initial spacing level
	// Output		void    : Output stream altered by FOp 
	//				  details being placed into it

void floq_op::print(std::ostream& ostr, int full, std::string spacer) const
  {
  int Fhs = (2*N+1)*hs;
  ostr << "\n" << spacer << "Photon Space Dimension:      " << N;
  ostr << "\n" << spacer << "Hilbert Space Dimension:     " << hs;
  ostr << "\n" << spacer << "Fourier Expansion Frequency: " << Om << " Hz";
  ostr << "\n" << spacer << "Full Space Dimension:        " << Fhs;
  if(!hs)
    ostr << "\n" << spacer << "Operator Representation:     NULL";
  else if(full)
    ostr << "\n" << spacer << "Operator Representation:\n"
         << (const gen_op&)(*this);
  }

std::ostream& operator<< (std::ostream& ostr, const floq_op& FOp)
  { FOp.print(ostr); return ostr;}


// ____________________________________________________________________________
// J              INTERNAL FLOQUET OPERATOR MANIPULATIONS
// ____________________________________________________________________________


void floq_op::add_omega()

        // Input             FOp     : Floquet operator (this)
        // Output            FOp     : Floquet operator with omegas added
        //                             on main diagonal
            
  {
  int fs = (2*N+1)*hs;			// Floquet space dimension
  matrix mx(fs, fs, i_matrix_type);	// I matrix is Floquet space
  basis bs(mx);				// Default basis in Floquet space
  gen_op::Op_base(bs);			// Put FOp in default basis
  for(int n=0; n<=2*N; n++)		// Add Omega's to the operator blocks
    for(int h=0; h<hs; h++)
      gen_op::put( gen_op::get(n*hs+h,n*hs+h) + (n-N)*Om, n*hs+h, n*hs+h);
  }

        // Input             FOp     : Floquet operator (this)
        // Output            FOp     : Floquet operator with omegas subtracted
        //                           : on main diagonal   

void floq_op::sub_omega()
  {
  int fs = (2*N+1)*hs;			// Floquet space dimension
  matrix mx(fs, fs, i_matrix_type);	// I matrix is Floquet space
  basis bs(mx);				// Default basis in Floquet space
  gen_op::Op_base(bs);			// Put FOp in default basis
  for(int n=0; n<=2*N; n++)		// Subtr. Omegas from op. blocks
    for (int h=0;h<hs;h++)
      gen_op::put( gen_op::get(n*hs+h,n*hs+h) - (n-N)*Om, n*hs+h, n*hs+h);
  }

#endif						// FloqOp.cc
 
