/* basis.cc *****************************************************-*-c++-*-
**									**
**                                 G A M M A				**
**									**
**	BASIS					   Implementation	**
**									**
**	Copyright (c) 1991, 1992, 1997					**
**									**
**      Scott Smith                                                     **
**      Copyright (c) 1996                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**									**
**	Tilo Levante							**
**	Eidgenoessische Technische Hochschule				**
**	Labor fuer physikalische Chemie					**
**	8092 Zuerich / Switzerland					**
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
** Description								**
**									**
** The class basis provides arrays that defines a basis in GAMMA for	**
** operators and superoperators.  These are static private structures	**
** that are difficult to alter.  The class currently supports composite	**
** spaces, i.e. bases that span multiple sub-spaces.			**
**									**
*************************************************************************/

#ifndef   Gbasis_cc_			// Is file already included?
#  define Gbasis_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif


#include <HSLib/Basis.h>		// Include the header info
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <Basics/StringCut.h>		// Include Gdec function
 
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                       CLASS BASIS ERROR HANDLING
// ____________________________________________________________________________

/*      Input		       bs      : Basis (this)
                               eidx    : Error index
                               noret   : Flag for linefeed (0=linefeed)
                               pname   : Added error message for output
       Output                  none    : Error message output
                                         Program execution stopped if fatal

  The following error messages use the defaults set in the Gutils package

               Case                          Error Message

               (0)                     Program Aborting.....
               (1)                     Problems With Input File Stream
               (2)                     Problems With Output File Stream
               default                 Unknown Error                        */

void basis::BSerror(int eidx, int noret) const
  {
  std::string hdr("Basis");
  std::string msg;
  switch(eidx)
    {    
    case 11: msg = std::string("Basis must be square!");         // (11)
             GAMMAerror(hdr, msg, noret); break;
    case 12: msg = std::string("Improper # of sub- components"); // (12)
             GAMMAerror(hdr, msg, noret); break;
    case 13: msg = std::string("Dimensioning mismatch");         // (13)
             GAMMAerror(hdr, msg, noret); break;
    case 14: msg = std::string("Basis component out of range");  // (14)
             GAMMAerror(hdr, eidx, noret); break;
    case 15: msg = std::string("Need at least 1 sub-component"); // (15)
             GAMMAerror(hdr, msg, noret); break;
    case 16: msg = std::string("Bad assigment from matrix");     // (16)
             GAMMAerror(hdr, msg, noret); break;
    default: GAMMAerror(hdr, eidx, noret); break;	// Unknown Error       (-1)
    }
  }

void basis::BSerror(int eidx, const std::string& pname, int noret) const
  {
  std::string hdr("Basis");
  std::string msg;
  switch(eidx)
    {    
    case 16: msg = std::string("Summed Subspace Dimension Is "); // (16)
             GAMMAerror(hdr, msg + pname, noret); break;
    case 17: msg = std::string("Full Dimension Is ");            // (17)
             GAMMAerror(hdr, msg + pname, noret); break;
    case 18: msg = std::string("Number Of Components Is ");      // (18)
             GAMMAerror(hdr, msg + pname, noret); break;
    case 19: msg = std::string("Index Of Accessed Component ");  // (19)
             GAMMAerror(hdr, msg + pname, noret); break;
    default: GAMMAerror(hdr, eidx, pname, noret); break;// Unknown Error       (-1)
    }
  }

volatile void basis::BSfatal(int eidx) const
  {                                                                 
  BSerror(eidx, 1);				// Output error message
  if(eidx) BSerror(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

volatile void basis::BSfatal(int eidx, const std::string& pname) const
  {                                                                 
  BSerror(eidx, pname, 1);			// Output error message
  if(eidx) BSerror(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }
 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
  
// ____________________________________________________________________________
// A                     BASIS CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ---------------------------------------------------------------------------
//                           Simple Constructors
// ---------------------------------------------------------------------------

/* The produce default bases, essentially identity matrices either in an
   single space or in a composite space.                                    */

basis::basis() : matrix()
  {
  nc     = 1;					// Only 1 component
  ncd    = new int[nc];				// Array of dimensions
  ncd[0] = 0;					// No dimension
  bname  = "";					// No name
  }

basis::basis(int dim) : matrix (dim, dim, i_matrix_type)
  {
  nc     = 1;					// Only 1 component
  ncd    = new int[nc];				// Array of dimensions
  ncd[0] = dim;					// Set dimension
  bname  = "Default Basis";			// Set basis name
  }

basis::basis(const std::vector<int> dims)
  {
  nc = dims.size();				// This many componets
  int dim = 0;
  ncd    = new int[nc];				// Array of dimensions
  for(int i=0; i<nc; i++) 			// Set dimensions
    { ncd[i] = dims[i]; dim += dims[i]; }
  matrix imx(dim,dim,i_matrix_type);		// Default basis array
  matrix::operator=(imx);			// Set basis array
  bname  = "";					// No name
  }

        // Input             mx : Conversion matrix
        // Output               : Returns a basis of matrix mx

basis::basis(const matrix& mx, int NC, int* NCD) : matrix(mx) 
  { 
  if(mx.rows() != mx.cols())			// Insure squre array
    { 
    BSerror(11,1);				//   Must use square array
    BSfatal(9);					//   Error in construction
    }
  if(!NC)					// Must be >= 1 component
    {
    BSerror(12,1);				//   Improper subcomponents
    BSerror(15,1);				//   Need at least 1
    BSfatal(9);					//   Error in construction
    } 
  ncd = new int[NC];				// Array for component dims
  nc = NC;					// Store number components
  if(nc==1) ncd[0] = mx.rows();			// If 1 component, store size
  else 						// For multi-components
    {
    int dimS = 0;				// 	Track sub-spaces
    for(int i=0; i<nc; i++)			//	Loop all components
      {
      dimS += NCD[i]; 				//	Total space size
      ncd[i] = NCD[i];				//	Store sub-space size
      }
    if(dimS != mx.rows())			// 	Insure spaces match!
      {
      BSerror(13,1);				//      Dimension mismatch
      BSerror(16, Gdec(dimS), 1);
      BSerror(17, Gdec(mx.rows()), 1);
      BSfatal(9);				//   Error in construction
      }
    }
  bname  = "";					// No name
  }


basis::basis(const basis& bs) : matrix((const matrix&) bs)
  {
  nc = bs.nc;					// Copy number of components
  ncd = new int[nc];				// Storage for subspaces
  for(int i=0; i<nc; i++)			// Copy sub space dimensions
    ncd[i] = bs.ncd[i];
  bname  = bs.bname;				// Copy name
  }

        // Input             bs : An existing basis
        //                   mx : A matrix for a basis
        // Output               : New basis, which is the old basis
	//			  but with a new matrix
	// Note			: This constructor allows a new basis
	//			  to be formed (say from a similarity
	//			  transformation) which retains the
	//			  subspace structure of a related basis

basis::basis(const basis& bs, const matrix& mx) : matrix(mx)

  {
  nc = bs.nc;					// Copy number of components
  ncd = new int[nc];				// Array for component sizes
  for(int i=0; i<nc; i++)			// Copy the subspace sizes
    ncd[i] = bs.ncd[i];
  bname  = "";					// No name
  }


basis::~basis() { delete [] ncd; }		// Just kill components array

basis& basis::operator=(const basis& bs)
  {
  if(this == &bs) return *this;			// Nothing if self assign
  matrix::operator=((matrix)bs);		// Copy basis array
  nc = bs.nc;					// Copy number of components
  delete [] ncd;				// Remove component array
  ncd = new int[nc];				// New component array
  for(int i=0; i<nc; i++)			// Copy all subspace sizes
    ncd[i] = bs.ncd[i];
  bname  = bs.bname;				// Copy name
  return *this;
  }

basis& basis::operator=(const matrix& mx)
  {
  if(mx.rows() != mx.cols())
  if(mx.rows() != mx.cols())			// Insure squre array
    { 
    BSerror(11,1);				//   Must use square array
    BSfatal(15);				//   Bad assignment
    }
  matrix::operator=(mx);			// Copy basis array
  nc = 1;					// Copy number of components
  delete [] ncd;				// Remove component array
  ncd = new int[nc];				// New component array
  ncd[0] = mx.rows();				// 1 component, store size
  bname  = "";					// No name
  return *this;
  }

  
// ____________________________________________________________________________
// B                       BASIC BASIS FUNCTIONS
// ____________________________________________________________________________
 
int         basis::dim()  const { return rows(); }
int         basis::size() const { return rows(); }
std::string basis::name() const { return bname;  }
void        basis::name(const std::string& nm) { bname = nm; }

// ____________________________________________________________________________
// C            BASIC BASIS FUNCTIONS SUPPORTING SUBSPACES
// ____________________________________________________________________________

/* These function deal with the dimension of the basis space.  Note that when
   a basis contains sub-spaces its Hilbert space dimension will simply sum all
   sub-spaces together.  In contrast, the total Liouville space will be the sum
   of individual Liouvillie spaces, the sum of the squares of each Hilbert sub
   space (NOT the square of the total Hilbert space)

    Function    Arguments                   Output
   ----------   ---------  ----------------------------------------------------
   dim_LS        ------    Dimension of Basis Liouville space (sum of squares)
   sub_N         ------    Number of Basis components (sub-spaces)
   sub_dim        icmp     Dimension of component icmp's Hilbert space 
   sub_anchor     icmp     Start index of component icmp in full Hilbert space
   sub_anchor_LS  icmp     Start index of component icmp in full Liouville sp.
           
*/

int basis::dim_LS() const
  {
  int dimLS = 0;			// Start with zero
  for(int i=0; i<nc; i++)		// Loop over components
    dimLS += ncd[i]*ncd[i];		// Add square of component subspace
  return dimLS;
  }

int basis::sub_N() const { return nc; }
int basis::sub_dim(int ic) const
  { 
  if((ic>=0) && (ic<nc)) 		// Check if component exists
    return ncd[ic]; 			// Return subspace dimension
  else return 0;
  }

int basis::sub_anchor(int ic) const
  {
  int anchor = 0;
  int i=0;
  for(; i<ic; i++) anchor += ncd[i];
  return anchor;
  } 

int basis::sub_anchor_LS(int ic) const
  {
  if(ic > nc)
     {
     BSerror(14,1);			//   Component out of range
     BSerror(18, Gdec(nc), 1);		//   Number of components
     BSfatal(19, Gdec(ic));
     }
  int anchor = 0;			// Start at <0|bs|0>
  for(int i=0; i<ic; i++)		// Loop over the components
    anchor += ncd[i]*ncd[i];		// Add the block's Liouvilles space
  return anchor;
  }



	// Input		 i    : Number of basis state 
	// Output		 int  : Number of sub-component
	// 			        to which this state belongs 
// sosi - should clean this up someday, give compiler warnings?

int basis::which_sub_LS(int i) const
  {
  int j=0;
  if(i<0 || i>=dim_LS()) 
     {
     BSerror(14,1);			//   Component out of range
     BSerror(18, Gdec(i), 1);		//   Number of components
     BSfatal(19, Gdec(dim_LS()));
     }
  for(j=0; j<nc; j++)			// Loop over all components
    if(i< sub_anchor_LS(j))		// If state resides below this
      return(j-1);			// block then return previous one
  return(nc-1);
  }
  

basis defLSbasis(const basis& bs)

        // Input            bs  : A (Hilbert space) basis
        // Output           Lbs : A default Liouville space basis
        //                        which has a sub-space Liouville structure
	//			  that coincides with bs Hilbert sub-space(s)

  {
  int ls = bs.dim_LS();				// Get Liouville dimension
  basis bs1(matrix(ls,ls,i_matrix_type)); 	// Begin new basis
  bs1.nc = bs.nc;				// Copy number of components
  bs1.ncd = new int[bs1.nc];			// Storage for subspaces
  for(int i=0; i<bs1.nc; i++)			// Copy sub space dimensions
    bs1.ncd[i] = (bs.ncd[i])*(bs.ncd[i]);
  return bs1;
  }
 

basis defbasis(const basis& bs)
 
        // Input            bs  : A (Hilbert space) basis
        // Output           hbs : A default basis which has a sub-space
        //                        structure that coincides with bs

  {
  int hs = bs.dim();				// Get basis dimension
  basis bs1(hs); 				// Begin new basis
  bs1.nc = bs.nc;				// Copy number of components
  bs1.ncd = new int[bs1.nc];			// Storage for subspaces
  for(int i=0; i<bs1.nc; i++)			// Copy sub space dimensions
    bs1.ncd[i] = bs.ncd[i];
  return bs1;
  }
 
// ____________________________________________________________________________
// D                 BASIS ARRAY & TRANSFORMATION ACCESS
// ____________________________________________________________________________

/* These functions allow access to the basis transformation array and apply
   that transformation to an input matrix.

   Function   Arguments                           Result 
 ------------ --------- ------------------------------------------------------- 
      U         ----    Returns the basis array as a matrix
   convert       mx     Applies basis transformation (put mx in bs)
 convert_back    mx     Applies inverse basis transform (take mx out of bs)

                           -1             bs          *t
        convert:   mx1 = bs  * mx * bs =======> = (bs)  * mx * bs
                                       unitary

	                             -1   bs                   *t
   convert_back:   mx1 = bs * mx * bs   =======> bs * mx * (bs)
	                                unitary

The two "convert" functions below are at the moment very poor. This is because
GAMMA quickly looses track of whether an array (basis) is unitary... forcing
us to to a time-consuming check to see if it is.                             */
 
matrix basis::U()          const { return (matrix)*this; }
matrix basis::get_matrix() const { return (matrix)*this; }
matrix basis::get_mx()     const { return (matrix)*this; }

matrix basis::convert(const matrix& mx) const
  {
  if(((matrix)*this).is_unitary())
    return adjoint_times((matrix)*this,mx) * (matrix)*this;
  else
    return inv((matrix)*this) * mx * (matrix)*this; 
  }

matrix basis::convert_back(const matrix& mx) const
  {
  if(((matrix)*this).is_unitary())
    return ((matrix)*this)*times_adjoint(mx,(matrix)*this);
  else
    return ((matrix)*this) * mx * inv((matrix)*this); 
  }

// ____________________________________________________________________________
// E                    MORE COMPLEX BASIS FUNCTIONS
// ____________________________________________________________________________
 
// These functions are also "kludges" in the sense that better matrix handling
// would make them totally unnecessary.  A scan of the GAMMA sources should
// reveal that these functions are rarely used (I Hope)

// ***************************** tensor product *******************************

basis tensor_product(const basis& bs1, const basis& bs2)

        // Input            bs1 : A basis
        //                  bs2 : A second basis
        // Output           pdt : A basis which is the tensor product
        //                        of the two input bases
	// Note			: This mirrors the handling of tensor
	//			  products by class matrix, but because
	//			  of private inheritance in class basis
	//			  it needs to be in this class too.

        //                           pdt        =   bs1 (x) bs2

        //                       (m*o x n*p)       (mxn)   (oxp)

        //                    <i*o+k|pdt|j*p+l> = <i|bs1|j><k|bs2|l>
	// Note			: This could have been used to produce
	//			  mult_sys bases, but we choose not to.
	//			  Now tensor_product cannot be used for
	// 			  multi_sys bases.

  {
  matrix mx1 = bs1.U();				// Copy basis 1 matrix
  matrix mx2 = bs2.U();				// Copy basis 2 matrix
  return basis(tensor_product(mx1, mx2));	// Return product
  }


//basis U_transform(const basis& bs1);

        // Input            bs1 : A basis
        // Output           pdt : A unitary transform basis which is
	//			  a Liouvile space equivalent of basis
        //                        transform.
	// Note			: This mirrors the handling of unitary
	//			  transformation superoperators done
	//			  by class super_op.


 
// ____________________________________________________________________________
// F                    BASIS COMPARISON/TEST FUNCTIONS
// ____________________________________________________________________________


bool basis::operator==(const basis& bs2) const
  {
  if(!check(bs2)) return false;		// Check compatibility 1st
  if(isDefaultBasis())			// Next quick check if both
    return(bs2.isDefaultBasis());	// default bases
  else					// Now check arrays
    {
    if(bs2.isDefaultBasis()) return false;
    else return((const matrix&)*this == (const matrix&)bs2);
    }
  }


bool basis::operator!=(const basis& bs2) const { return !(*this == bs2); }

bool basis::isDefaultBasis() const
  {
  if(stored_type() == i_matrix_type) return true;
  else
    return((const matrix&)*this==matrix(rows(),cols(),i_matrix_type));
  }


int basis::refs() const { return ((const matrix&)*this).refs(); }

bool basis::check(const basis& bs1) const
  {
  if(nc != bs1.sub_N()) return false;		// Not so if different components
  for(int i=0; i<nc; i++)			// Loop components and die if
    if(ncd[i] != bs1.sub_dim(i)) return false;	// any mismatches
  return true;
  }   
    
// ____________________________________________________________________________
// G                         BASIS I/O FUNCTIONS
// ____________________________________________________________________________

// ------------------------ ASCII Output Functions ----------------------------

/*              Input		bs   : Basis (this)
                                ostr : Output ASCII file stream
                                full : Flag for amount of output
                Return          void : bs is sent to the output stream       */


std::ostream& basis::print(std::ostream& ostr, int full) const
  {
  int N = sub_N();
  if(isDefaultBasis())
    ostr << "Default Basis ("
         << dim() << " x "
         << dim() << ") Identity Matrix";
  else ostr << (matrix)*this;
  if(N>1) 
    {
    ostr << "\n" << N << " Components: ";
    for(int i=0; i<N; i++)
      {
      ostr << "Subspace " << i << " = " << sub_dim(i);
      if(i != N-1) ostr << ", ";
      }
    }
  if(full) ostr << "\n";
  return ostr;
  }

std::ostream& operator << (std::ostream& ostr, const basis& bs) {return bs.print(ostr);}


// ------------------------ Binary Output Functions ---------------------------

/*              Input		bs   : Basis (this)
                                fn   : Input binary filename
                                fp   : File stream (pointing at bs
                Return          void : Op is written to either the
                                       specified file or filestream.
                Note                 : Output format is partially set by
                                       class matrix (matrix typing)          */
 
void basis::write(const std::string& fn) const
  {
  std::ofstream fp;				  	// Construct a file stream
  fp.open(fn.c_str(), std::ios::out|std::ios::binary);	// Open the file for output
  write(fp);						// Write Op to file stream
  fp.close();						// Close file stream
  }


std::ofstream& basis::write(std::ofstream& fp) const
  {
  fp.write((char*)&nc,sizeof(int));		// Output the number of components
  for(int i=0; i<nc; i++)			// Output the component dimension
    fp.write((char*)&(ncd[i]),sizeof(int));
  ((matrix)*this).write(fp, 0);			// Now output bs matrix, GAMMA format
  return fp;
  }

// ------------------------ Binary Input Functions ----------------------------

/*              Input		bs   : Basis (this)
                                fn   : Input binary filename
                                fp   : File stream (pointing at bs)
                Return          void : bs is read in from either the
                                       specified file or filestream.         */  
 
void basis::read(const std::string& fn)
  {
   std::ifstream fp;					// Construct a file stream
  fp.open(fn.c_str(),std::ios::in|std::ios::binary);	// Open file for reading
  read(fp);						// Read bs, use funct. overload
  fp.close();						// Close the file stream
  }

 std::ifstream& basis::read( std::ifstream &fp)
  {
  fp.read((char*)&nc,sizeof(int));		// Read the number of components
  if(ncd) delete [] ncd;			// Remove existing dimensions
  ncd = new int[nc];				// Allocate array for dimensions
  for(int i=0; i<nc; i++)			// Input the component dimension
    fp.read((char*)&(ncd[i]),sizeof(int));
  matrix::read(fp);				// Read in basis array
  return fp;
  }
    
// ____________________________________________________________________________
// H                         BASIS TESTING FUNCTIONS
// ____________________________________________________________________________

/* Since a Basis is fundamental to both operators and superoperators in GAMMA
   I have added these tests to provide a means to insure they are working
   correctly. It is paramount that we know how to properly convert an operator
   (or superoperator) between bases. In order for that to occur we MUST be able
   to obtain the inverse of the basis transformation array. In certains cases
   this is just the adjoint (if the basis is unitary) but often times we will
   not be so fortunate, such as when dealing with relaxation superoperator.
   A good test would then be U * inv(U) = 1. That will be proivded here.     */

double basis::TestBasis(int pf) const
  {
  matrix S   = (matrix)(*this);			// The U array
  matrix Sm1 = inv(S);				// The inverse of U
  matrix I   = S*Sm1;				// Hope to generate I
  int nr     = rows();
  matrix Z   = I - matrix(nr,nr,i_matrix_type);	// Now this is zero?
  double dev = norm(Z.maxZ()); 			// This is max. deviation
  if(pf)
    {
    if(pf>1)
      {
      std::cout << "\n\n" << std::string(30, ' ') << "Basis Being Tested" << S;
      if(pf > 2)
        {
        std::cout << "\n\n" << std::string(30, ' ') << "Basis Inverse"      << Sm1;
        std::cout << "\n\n" << std::string(23, ' ') << "S * inverse(S) = I?" << I;
        }
      }
    std::cout << "\n\t\tLargest Deviation From S*Sm1-I: " << dev;
    }
  return dev;
  }

#endif 							// Basis.cc
