/* decomp.cc ********************************************-*-c++-*-
**                                                              **
**                            G A M M A                         **
**                                                              **
**      Decomp 					Implementation	**
**                                                              **
**      Copyright (c) 1996                                      **
**      Scott Smith                                             **
**      National High Magnetic Field Laboratory                 **
**      1800 E. Paul Dirac Drive                                **
**      Tallahassee Florida, 32306-4005                         **
**                                                              **
**      $Header: $
**                                                              **
*****************************************************************/

/*****************************************************************
**                                                              **
**  Description                                                 **
**                                                              **
**  Class decomp provides the means with which an operator can  **
**  be readily decomposed into a set of basis operators.  Each  **
**  variable of class decomp may contains a set of such basis   **
**  operators.  Functions are provided which will construct     **
**  the basis set and which will break down a provided operator **
**  into the basis components.                                  **
**                                                              **
*****************************************************************/

#ifndef   Gdecomp_cc_			// Is file already included?
#  define Gdecomp_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <GamGen.h>			// Include OS specific stuff
#include <Level2/BaseDecomp.h>		// Include the header
#include <string>			// Knowledge of strings
#include <Matrix/complex.h>		// Knowledge of complex numbers
#include <HSLib/GenOp.h>		// Knowledge of operators
#include <HSLib/SpinSys.h>		// Knowledge of spin systems
#include <HSLib/SpinOp.h>		// Knowledge of spin operators
#include <HSLib/SpinOpCmp.h>		// Knowledge of spinoperators
#include <Matrix/row_vector.h>		// Knowledge of row vectors
#include <Basics/StringCut.h>
#include <Basics/Gutils.h>              // Need GAMMA error messages

// sosi: Still think about the following things 
// 1.) Test the product base decomposition
// 2.) Completely replace product base function currently in GAMMA
// 3.) Add a sorting by coherence order, name
// 4.) Perhaps must set default values for the arrays?
// 5.) Add different decompositions such as single transition...
// 6.) Add the ability to set a basis
// 7.) Insure all base operators share a basis
// 8.) Generalize the tensor base expansion
// 9.) Generalize the constructors so the user has more control

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                   CLASS DECOMPOSITION ERROR HANDLING
// ____________________________________________________________________________

/*      Input                   OD      : Operator decomposition (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

void decomp::ODerror(int eidx, int noret) const
  {
  std::string hdr("Operator Basis");
  std::string msg;
  switch (eidx)
    {
    case  3:GAMMAerror(hdr,"Unknown Basis Expansion",      noret); break;// (3)
    case  4:GAMMAerror(hdr,"Access Is Out Of Range",       noret); break;// (4)
    case  5:GAMMAerror(hdr,"Cant Get Operator Name",       noret); break;// (5)
    case  6:GAMMAerror(hdr,"Cant Get Requested Operator",  noret); break;// (6)
    case  7:GAMMAerror(hdr,"Cant Get Requested Intensity", noret); break;// (7)
    case  8:GAMMAerror(hdr,"Cant Get Requested Coherence", noret); break;// (8)
    case  9:GAMMAerror(hdr,"Base Operator Does Not Exist", noret); break;// (9)
    case 10:GAMMAerror(hdr,"Liouville Dimension Mismatch", noret); break;// (10)
    default:GAMMAerror(hdr,eidx,noret); break;
    }
  }

volatile void decomp::ODfatal(int eidx) const
  {
  ODerror(eidx, 1);
  if(eidx) ODerror(0);
  GAMMAfatal();					// Clean exit from program
  }

// ______________________________________________________________________
// II                     SPECIFIC DECOMPOSITIONS
// ______________________________________________________________________


std::vector<int> decomp::sub_indices(int ils, int ssdim, int nss)
  {
  std::vector<int> idxs(nss, 0);		// Array of spin Op indices
  int iss = 0;				// Begin with 1st component
  int nssm1 = nss-1;			// Done in reverse order
  while(ils != 0)			// Split up the full index
    {					// from last component to
    idxs[nssm1-iss] = ils % ssdim;	// first component.  But store 
    ils /= ssdim;			// the indices from first to last
    iss++;
    }
  return idxs;
  }

void decomp::sub_indices(int indices[], int ils, int ssdim, int nss)

	// Input		indices	: Basis indices per subspace
	//			ils     : Full space index
	//			ssdim   : Dimension of a single subspace
	//			nss	: Number components in full space
	// Output		void    : Fills array "indices" with
	//				  sub-space indices corresponding
	//				  to the full space index ils

  {
  int iss = 0;				// Begin with 1st component
  int nssm1 = nss-1;			// Done in reverse order
  while(ils != 0)			// Split up the full index
    {					// from last component to
    indices[nssm1-iss] = ils % ssdim;	// first component.  But store 
    ils /= ssdim;			// the indices from first to last
    iss++;
    }
  }

void decomp::spin3halves(const spin_sys& sys) 

	// Input		sys	: Spin system
	// Output		void    : Fills dec with irreducible
	//				  spherical tensor components
	// Note				: This pertains ONLY to a system
	//				  containing a single spin I=3/2
	//				  (so it's private & taken care of)

// The basis operators for this treatment were taken out of the article
// "Relaxation-Induced Violations of Coherence Transfer Selection Rules
// in Nuclear Magnetic Resonance", N. Mueller, G. Bodenhausen, and
// R.R. Ernst, J. Magn. Reson., 75, 297, 334 (1987).

  {
//                 Generate Required Base Operators

  gen_op FE = Fe(sys);                          // Here are the Cartesian Ops.
  gen_op FX = Fx(sys);                          // which we will use to
  gen_op FY = Fy(sys);                          // build up the spin operators
  gen_op FZ = Fz(sys);                          // for the spherical tensors
  gen_op FP = Fp(sys);
  gen_op FM = Fm(sys);

//                      Required Scaling Factors

  BaseSpins = std::vector<int>(16, 0);		// All Ops involve a single spin

//                 Generate Rank 0 Basis Operators
//                           (1 Operator)

  double l0fact = 0.5; 				// Scaling factor, l=0 terms
  BaseOps[0] = l0fact*FE;                      	// This is T00
  BaseNames[0] = "T00 ";			// This is T00
  BaseAltNames[0] = " [1/2]*E";			// This is T00
  BaseCoherences[0] = 0;			// This is T00
  BaseCoefficients[0] = l0fact;			// This is T00
  BaseCoeffNames[0] = "2";			// This is T00
 
//                 Generate Rank 1 Basis Operators
//                          (3 Operators)
 
  double l1fact = 1.0/sqrt(5.0);                // Scaling factor, l=1 terms
  BaseOps[1] = (l1fact/sqrt(2.0))*FM;		// This is T1-1
  BaseOps[2] = l1fact*FZ;                      	// This is T10
  BaseOps[3] = (-l1fact/sqrt(2.0))*FP; 		// This is T11
  BaseNames[1] = "T1-1";			// This is T1-1
  BaseNames[2] = "T10 ";			// This is T10
  BaseNames[3] = "T11 ";			// This is T11
  BaseAltNames[1] = " [1/sqrt(10.0)]*I-";	// This is T1-1
  BaseAltNames[2] = " [1/sqrt(5.0)] *Iz";	// This is T10
  BaseAltNames[3] = "-[1/sqrt(10.9)]*I+";	// This is T11
  BaseCoherences[1] = -1;			// This is T1-1
  BaseCoherences[2] = 0;			// This is T10
  BaseCoherences[3] = 1;			// This is T11
  BaseCoefficients[1] = 1/sqrt(10.0);		// This is T1-1
  BaseCoefficients[2] = 1/sqrt(20.0);		// This is T10
  BaseCoefficients[3] = 1/sqrt(10.0);		// This is T11
  BaseCoeffNames[1] = "sqrt(10.0)";		// This is T1-1
  BaseCoeffNames[2] = "sqrt(20.0)";		// This is T10
  BaseCoeffNames[3] = "sqrt(10.0)";		// This is T11
 
//                 Generate Rank 2 Basis Operators
//                        (5 Operators)
 
  double l2fact = 1.0/sqrt(6.0);			// Scaling factor
  BaseOps[4] = 0.5*l2fact*FM*FM;                	// This is T2-2
  BaseOps[5] = 0.5*l2fact*(FZ*FM+FM*FZ);        	// This is T2-1
  BaseOps[6] = (l2fact/sqrt(6.0))			// This is T20
             * (3.0*FZ*FZ-(15./4.)*FE);
  BaseOps[7] =-0.5*l2fact*(FZ*FP+FP*FZ);        	// This is T21
  BaseOps[8] = 0.5*l2fact*FP*FP;                	// This is T22
  BaseNames[4] = "T2-2";				// This is T2-2
  BaseNames[5] = "T2-1";				// This is T2-1
  BaseNames[6] = "T20 ";				// This is T20
  BaseNames[7] = "T21 ";				// This is T21
  BaseNames[8] = "T22 ";				// This is T22
  BaseAltNames[4] = " [1/sqrt(24.)]*I-I-";		// This is T2-2
  BaseAltNames[5] = " [1/sqrt(24.)]*(IzI- + I-Iz)";	// This is T2-1
  BaseAltNames[6] = " [1/sqrt(36.)]*(3IzIz - I(I+1))";	// This is T20
  BaseAltNames[7] = "-[1/sqrt(24.)]*(IzI+ + I+Iz)";	// This is T21
  BaseAltNames[8] = " [1/sqrt(24.)]*I+I+";		// This is T22
  BaseCoherences[4] = -2;				// This is T2-2
  BaseCoherences[5] = -1;				// This is T2-1
  BaseCoherences[6] = 0;				// This is T20
  BaseCoherences[7] = 1;				// This is T21
  BaseCoherences[8] = 2;				// This is T22
  BaseCoefficients[4] = 1/sqrt(2.);			// This is T2-2
  BaseCoefficients[5] = 1/sqrt(2.);			// This is T2-1
  BaseCoefficients[6] = 0.5;				// This is T20
  BaseCoefficients[7] = 1/sqrt(2.);			// This is T21
  BaseCoefficients[8] = 1/sqrt(2.);			// This is T22
  BaseCoeffNames[4] = "sqrt(2.)";			// This is T2-2
  BaseCoeffNames[5] = "sqrt(2.)";			// This is T2-1
  BaseCoeffNames[6] = "2";				// This is T20
  BaseCoeffNames[7] = "sqrt(2.)";			// This is T21
  BaseCoeffNames[8] = "sqrt(2.)";			// This is T22
 
//                 Generate Rank 3 Basis Operators
//                         (7 Operators)
   
  double l3fact = sqrt(2.0/9.0);			// Scaling factor
  gen_op X = 5.0*FZ*FZ-(17.0/4.0)*FE;
  BaseOps[9] = (l3fact/sqrt(8.))*FM*FM*FM;		// This is T3-3
  BaseOps[10]= (l3fact*sqrt(0.75))*FM*FZ*FM;     	// This is T3-2
  BaseOps[11]= (l3fact*sqrt(3./160.)) * (FM*X + X*FM); 	// This is T3-1
  BaseOps[12]= (l3fact/sqrt(10.))			// This is T30
             * (5.*FZ*FZ - (41./4.)*FE)*FZ;
  BaseOps[13]=-(l3fact*sqrt(3./160.)) * (FP*X + X*FP);	// This is T31
  BaseOps[14]= (l3fact*sqrt(0.75))*FP*FZ*FP;     	// This is T32
  BaseOps[15]=-l3fact/sqrt(8.)*FP*FP*FP;	        // This is T33
  BaseNames[9]  = "T3-3";				// This is T3-3
  BaseNames[10] = "T3-2";				// This is T3-2
  BaseNames[11] = "T3-1";				// This is T3-1
  BaseNames[12] = "T30 ";				// This is T30
  BaseNames[13] = "T31 ";				// This is T31
  BaseNames[14] = "T32 ";				// This is T32
  BaseNames[15] = "T33 ";				// This is T33
  BaseAltNames[9]  = " [1/sqrt(36.)]*I-I-I-";		// This is T3-3
  BaseAltNames[10] = " [1/sqrt(6.)] *I-IzI-";		// This is T3-2
  BaseAltNames[11] = std::string(" [1/sqrt(15)]*") 		// This is T3-1
                   + std::string("[I-(5IzIz-I(I+1)-E/2)")
                   + std::string(" + (5IzIz-I(I+1)-E/2)I-]"); 
  BaseAltNames[12] =" [1/sqrt(45.)]*[5IzIz-3I(I+1)+E]Iz";// This is T30
  BaseAltNames[13] = std::string("-[1/sqrt(15)]*") 		// This is T31
                   + std::string("[I+(5IzIz-I(I+1)-E/2)")
                   + std::string(" + (5IzIz-I(I+1)-E/2)I+]"); 
  BaseAltNames[14] = "-[1/sqrt(6.)] *I+IzI+";		// This is T32
  BaseAltNames[15] = "-[1/sqrt(36.)]*I+I+I+";		// This is T33
  BaseCoherences[9]  = -3;				// This is T3-3
  BaseCoherences[10] = -2;				// This is T3-2
  BaseCoherences[11] = -1;				// This is T3-1
  BaseCoherences[12] = 0;				// This is T30
  BaseCoherences[13] = 1;				// This is T31
  BaseCoherences[14] = 2;				// This is T32
  BaseCoherences[15] = 3;				// This is T33
  BaseCoefficients[9]  = 1;				// This is T3-3
  BaseCoefficients[10] = 1/sqrt(2.);			// This is T3-2
  BaseCoefficients[11] = 1/sqrt(5.);			// This is T3-1
  BaseCoefficients[12] = 1/sqrt(20.);			// This is T30
  BaseCoefficients[13] = 1/sqrt(5.);			// This is T31
  BaseCoefficients[14] = 1/sqrt(2.);			// This is T32
  BaseCoefficients[15] = 1;				// This is T33
  BaseCoeffNames[9]  = "1";				// This is T3-3
  BaseCoeffNames[10] = "sqrt(2.)";			// This is T3-2
  BaseCoeffNames[11] = "sqrt(5.)"; 			// This is T3-1
  BaseCoeffNames[12] = "sqrt(20.)";			// This is T30
  BaseCoeffNames[13] = "sqrt(5.)"; 			// This is T31
  BaseCoeffNames[14] = "sqrt(2.)";			// This is T32
  BaseCoeffNames[15] = "1";				// This is T33

  dname = "Single Spin 3/2";
  }


 
        // Input                sys     : Spin system
        // Output               void    : Fills dec with product operator
        //                                components
        // Note                         : This pertains ONLY to a system
        //                                containing a all spin I=1/2
        //                                (so it's private & taken care of) 

void decomp::product_operators(const spin_sys& sys)
  {
//                 Generate Required Base Operators

  int i, j, ns=sys.spins();		// Number of spins
  int ssdim = 4;			// 1 spin dimension {E, Ix, Iy, Iz}
  gen_op** Ixyz;			 // Array {Ix,Iy,Iz} operators
  Ixyz = new gen_op*[ns];
  for(i=0; i<ns; i++)			// Fill 'em up for each spin
    {					// in the system
    Ixyz[i] = new gen_op[ssdim-1];
    Ixyz[i][0] = Ix(sys,i);
    Ixyz[i][1] = Iy(sys,i);
    Ixyz[i][2] = Iz(sys,i);
    }

//                 Generate Required Base Operators

  std::string name;			// For individual operator name
  gen_op POp;				// For individual product operator
  gen_op E = Ie(sys,0);			// An identity operator
  long int q;				// Needed for scaling
  long int scaling = 1;			// Needed for scaling
  std::vector<int> ssindxs;
  int delm = 0;				// Coherence order of Op
  for(i=0; i<_LS; i++)
    {
    name = "";				// Start with no name
    POp = E;				// Start with "no" operator
    q =0;				// Start with "no" scaling factor
    delm =0;				// Start with no coherence order
    ssindxs = sub_indices(i,ssdim,ns);	// Generate base function code

// Given Individual Spin Components To Use, Here We Take Their Product
//      To Produce The Product Operator And A String Defining It

    for(j=0; j<ns; j++)			// Loop over spins
      {					// and form both the product
      switch(ssindxs[j])		// operator and its name
        {
        case 0:				//	Spin's contribution is E
          if(!i && !j) { name += "E"; } //	so don't do anything
          break;
        case 1:				// 	Spin's contribution is Ix
          POp &= Ixyz[j][0];		//	so add this one in
          name += "Ix" + Gdec(j);
          q++;
          delm += 1;
          break;
        case 2:				// 	Spin's contribution is Iy
          POp &= Ixyz[j][1];		//	so add this one in
          name += "Iy" + Gdec(j);
          q++;
          delm += 1;
          break;
        case 3:				// 	Spin's contribution is Iz
          POp &= Ixyz[j][2];		//	so add this one in
          name += "Iz" + Gdec(j);
          q++;
          break;
        default:			// 	Shouldn't ever get here
          break;
        }
      }
    if(i) scaling = long(pow(2.0,double(q-1)));	// POp normalization via EBW (2.1.87)
    BaseOps[i] = complex(scaling)*POp;  	// This is the product operator
    BaseSpins[i] = q;				// Spins involved in the operator
    if(scaling == 1) BaseNames[i] = name;  	// Add the scaling factor to the name
    else BaseNames[i] = Gdec(scaling)+name;  	// Add the scaling factor to the name
    BaseAltNames[i] = "";			// No alternate name for this
    BaseCoherences[i] = delm;			// This is the coherence order
    BaseCoefficients[i] = scaling;	
    if(scaling == 1)
      BaseCoeffNames[i] = std::string("1");  	// Add the scaling factor to the name
    else
      BaseCoeffNames[i] = std::string("1/")
                        + Gdec(scaling);			// This is T2-2
    }
  dname = "Product Operator";
  delete Ixyz;
  }


// ____________________________________________________________________________
//                  DECOMPOSITION COMMON INTERNAL FUNCTIONS
// ____________________________________________________________________________

bool decomp::ChkIndex(int i, bool warn) const
  {
  if(i>=0 && i<_LS) return true;
  if(warn) ODerror(4, 1);
  return false;
  }

bool decomp::ChkSize(int ls, bool warn) const
  {
  if(ls == _LS) return true;
  if(warn) ODerror(10, 1);
  return false;
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A              CLASS DECOMPOSITION CONSTRUCTORS/DETRUCTOR
// ____________________________________________________________________________

	
decomp::decomp ( )
  {
  _NS     = 0;				// No spins
  _LS     = 0;				// No Liouville space
  thresh  = 1.e-3;			// Set output threshold
  prif    = 0;				// Set print for reals
  }

decomp::decomp(const decomp &dec1)
  {
  _NS              = dec1._NS;			// Copy the number of spins
  _LS              = dec1._LS;			// Copy the Liouville space
  thresh           = dec1.thresh;		// Copy the output threshold
  prif             = dec1.prif;			// Copy the print flag
  BaseVals         = dec1.BaseVals;		// Copy the intensities
  BaseOps          = dec1.BaseOps;		// Copy basis operators
  BaseNames        = dec1.BaseNames;		// Copy base names
  BaseAltNames     = dec1.BaseAltNames;		// Copy base alternate names
  BaseCoeffNames   = dec1.BaseCoeffNames;	// Copy base coefficient names
  BaseCoherences   = dec1.BaseCoherences;	// Copy base coherence values
  BaseCoefficients = dec1.BaseCoefficients;	// Copy base coefficients
  BaseSpins        = dec1.BaseSpins;		// Copy # of spins
  }

decomp::decomp(const spin_sys& sys)
  {
  int hs = sys.HS();			// Get the system Hilbert space
  _LS    = hs*hs;			// Spin system Liouville space
  _NS    = sys.spins();			// Get the number of spins
  thresh = 1.e-3;			// Set output threshold
  prif = 0;				// Set print for reals
  if(_LS)
    {
    BaseOps          = std::vector<gen_op> (_LS);
    BaseNames        = std::vector<std::string> (_LS);
    BaseAltNames     = std::vector<std::string> (_LS);
    BaseCoeffNames   = std::vector<std::string> (_LS);
    BaseCoherences   = std::vector<int>    (_LS);
    BaseCoefficients = std::vector<double> (_LS);	
    BaseSpins        = std::vector<int>    (_LS);	// Copy # of spins
    BaseVals         = row_vector(_LS,complex0);
    }
// sosi - have to generalize this! Currently only two bases provided!

  if(sys.spins()==1 && hs==4)
    spin3halves(sys);
  else if(sys.spinhalf())
    product_operators(sys);
  else ODfatal(3);			// Unknown basis expansion
  return;
  }

decomp::~decomp () { }

decomp& decomp::operator= (const decomp &dec1)
  {
  _NS              = dec1._NS;			// Copy the number of spins
  _LS              = dec1._LS;			// Copy the Liouville space
  thresh           = dec1.thresh;		// Copy the output threshold
  prif             = dec1.prif;			// Copy the print flag
  BaseVals         = dec1.BaseVals;		// Copy the intensities
  BaseOps          = dec1.BaseOps;		// Copy basis operators
  BaseNames        = dec1.BaseNames;		// Copy basis names
  BaseAltNames     = dec1.BaseAltNames;		// Copy alternate basis names
  BaseCoeffNames   = dec1.BaseCoeffNames;	// Copy basis coefficient names
  BaseCoherences   = dec1.BaseCoherences;	// Copy base coherences
  BaseCoefficients = dec1.BaseCoefficients;	// Copy base coefficients
  BaseSpins        = dec1.BaseSpins;		// Copy # of spins
  return *this;
  }

// ______________________________________________________________________
// B                      DECOMPOSITION FUNCTIONS
// ______________________________________________________________________

        // Input                dec     : Decomposition (this)
	//			Op      : Operator to decompose
        // Output               void    : The operator is decomposed,
        //                                projected intensities are put
	//				  into the internal array

void decomp::decompose(const gen_op& Op) 
  {
  gen_op X(Op);				// Use a copy to do this
  for(int i=0; i<_LS; i++)		// Loop over all basis operators
    BaseVals.put(proj(X,BaseOps[i]),i);	// Project Op onto base Op i
  }

// ____________________________________________________________________________
// C                           ACCESS FUNCTIONS
// ____________________________________________________________________________

// ---------------------------- Dimension Access ------------------------------

int decomp::size() const { return _LS; }
int decomp::LS()   const { return _LS; }
int decomp::HS()   const { return _LS?BaseOps[0].dim():0; }

// ------------------------------- Name Access --------------------------------

        // Input                dec     : Decomposition (this)
	//			onames  : String for names
        // Output               void	: Names of the base operators
	//				  are copied into the String
	//				  array onames 

std::vector<std::string> decomp::Names() const { return BaseNames; }
std::vector<std::string> decomp::Names(int m) const
  {
  std::vector<std::string> onames;
  for(int i=0; i<_LS; i++) 
    {
    if(BaseCoherences[i] == m)
      onames.push_back(BaseNames[i]); 
    }
  return onames;
  }

void   decomp::Name(const std::string& name) { dname = name; }
std::string decomp::Name() const             { return dname; }

std::string decomp::OpName(int i) const
  {
  if(!ChkIndex(i)) ODfatal(5);
  return BaseNames[i];
  }

std::string decomp::AltOpName(int i) const
  {
  if(!ChkIndex(i)) ODfatal(5);
  return BaseAltNames[i];
  }

int decomp::MaxOpNameLen() const
  {
  int len = 0;
  for(int i=0; i<_LS; i++)
    len = gmax(len, int(BaseNames[i].length()));
  return len;
  }

int decomp::MaxOpAltNameLen() const
  {
  int len = 0;
  for(int i=0; i<_LS; i++)
    len = gmax(len, int(BaseAltNames[i].length()));
  return len;
  }

// ---------------------------- Coherence Access ------------------------------

int decomp::Coherence(int i) const
  {
  if(!ChkIndex(i)) ODfatal(8);
  return BaseCoherences[i];
  }

int decomp::MaxCoherence() const
  {
  int mc = 0;
  for(int i=0; i<_LS; i++)
    mc = gmax(mc, BaseCoherences[i]);
  return mc;
  }

// ----------------------------- Operator Access ------------------------------

        // Input                dec     : Decomposition (this)
        //                      Opname  : Name of a basis operator
        // Output               Op      : The base operator with
	//				  the name Opname 

gen_op decomp::Op(const std::string& Opname) const
  {
  int i = index(Opname);		// Get base index
  if(i != -1) return BaseOps[i];	// If good, return Op
  ODfatal(9);				// If no good, quit
  return complex0*BaseOps[0];
  }

gen_op decomp::Op(int i) const
  {
  if(!ChkIndex(i)) ODfatal(6);
  return BaseOps[i];
  }


// --------------------------- Value Access -----------------------------



        // Input                dec     : Decomposition (this)
        // Output               vx	: Row vector of the current
	//				  decomposition values

        // Input                dec     : Decomposition (this)
        //                      m       : Coherence order
        // Output               vx      : Row vector of the current
        //                                decomposition values
        //                                for base operators associated
        //                                with the coherence order m 

row_vector decomp::values() const { return BaseVals; }
row_vector decomp::values(int m) const
  {
  int i, j, num=0;
  for(i=0; i<_LS; i++) 
    if(BaseCoherences[i] == m) num++;
  row_vector vx(num);
  for(i=0, j=0; i<_LS; i++) 
    if(BaseCoherences[i] == m)
      vx.put(BaseVals.get(j++), i);
  return vx;
  }



        // Input                dec     : Decomposition (this)
        //                      i       : An index
        // Output               z	: Current decomposition value
	//				  of the ith basis operator

complex decomp::value(int i) const
  { 
  if(!ChkIndex(i)) ODfatal(7);
  return BaseVals.get(i);
  }


// -------------------- Operator Coefficient Access ---------------------


        // Input                dec     : Decomposition (this)
        //                      i       : An index
        // Output               z	: A coefficient which, when
	//				  taken out of the ith operator
	//				  matrix representation, makes
	//				  the matrrix elemnts nice

double decomp::bcoefficient(int i) const
  { 
  if(!ChkIndex(i)) ODfatal(7);
  return BaseCoefficients[i];
  }
 

// --------------------------- Index Access -----------------------------

 
int decomp::index(const std::string& Opname) const
 
        // Input                dec     : Decomposition (this)
        //                      Opname  : Name of a basis operator
        // Output               i       : Index of specified basis
        //                                operator.
        // Note                         : If requested operator does
        //                                not exist, return is -1
// sosi - this is probably broken now that Regex is gone....
// I don't remember why I parsed the operator name before the
// comparison with Opname & don't have time for it now... 
 
  {
  int j, found=0, i=-1;
//  Regex RXbaseOp("-?[0-9A-Za-z]+");
  std::string ntmp;
  for(j=0; j<_LS && !found; j++)
    {
    ntmp = BaseNames[j];
//    if(Opname == cut(ntmp, RXbaseOp))
    if(Opname == ntmp)

       { 
       found++;
       i = j; 
       } 
    }
  return i;
  }

// ____________________________________________________________________________
// U                          SORTING FUNCTIONS
// ____________________________________________________________________________

/* Allows for the basis operators to be sorted based on three possibilities.

         1.) The number of spin components in the operator
         2.) The name   of the  operator
         3.) The coherenece order of the operator

   Product operators are best displayed by first sorting by the number of
   spins involved, followed by the coherence order, followed by the operator
   name.                                                                     */

/*
void decomp::Sort(int A, int B, int C)
  {
  }
*/

std::vector<int> decomp::SortBySpins() const
  {
  std::vector<int> idx;				// Index array
  for(int n=0; n<=_NS; n++)			// Loop over max spins
    for(int i=0; i<_LS; i++)			// Loop over all basis Ops
      if(BaseSpins[i] == n)
        idx.push_back(i);
  return idx;
  }

// ____________________________________________________________________________
// X                         FRIEND FUNCTIONS
// ____________________________________________________________________________


void PB_dec(const spin_sys &sys, const gen_op &Op)
  { decomp D(sys); D.print(std::cout, Op); }


// ____________________________________________________________________________
// Z                 CLASS DECOMPOSITION I/O FUNCTIONS
// ____________________________________________________________________________

        // Input                dec     : Decomposition (this)
        //                      ostr    : An output stream
        // Output               ostr    : The output stream, modified to
        //                                contain the information about
        //                                the decomposition

std::ostream& decomp::print(std::ostream& ostr, int nc) const
  {
int an = 0;
  if(nc < 1) nc = 3;				// Insure at least 1 column
  std::string hdr = dname + " Basis";		// Header to output
  int len = hdr.length(); 			// Length of the header
  ostr << std::string(40 - len/2,' ') << hdr << "\n";// Output the header

  std::string BO("Base Op");				// Odd column header
  std::string BU("=======");				// Odd column underline
  std::string MQ("MQC");				// Even column header
  std::string MU("===");				// Even column underline
  int blen  = BO.length();			// Length of odd columns
  int cslen = 3;				// Separation between cols
  std::string csep = std::string(cslen, ' ');		// Column separation string
  int mslen = 4;				// Midcolumn separation
  std::string msep(mslen, ' ');


  int nlmax = MaxOpNameLen();			// Maximum name length
  int anmax = MaxOpAltNameLen();		// Maximum alt name length
  if(an)           nlmax += cslen + anmax;	// Add alternate names
  if(blen > nlmax) nlmax = blen;		// Names at least this long
  else if(nlmax > blen)				// If names longer than odd
    {						// column header then adjust
    if(blen > nlmax)
      {
      BO =  std::string((blen-nlmax)/2, ' ') + BO;	// the column header so it
      BU =  std::string((blen-nlmax)/2, '=') + BU;	// is this width also
      }
    if(nlmax > int(BO.length()))
      {
      BO += std::string(nlmax-BO.length(), ' ');
      BU += std::string(nlmax-BU.length(), '=');
      }
    }
  len = nlmax + cslen + 3;
  len *= nc;
  len += (nc-1) * 4;
  std::string spc(40-len/2, ' ');
  int i;


  std::vector<int> SI = SortBySpins();		// Sorted indices
  ostr << "\n\n" << spc;
  for(i=0; i<nc; i++) 
    { 
    ostr << BO << csep << MQ;
    if(i<nc-1) ostr << msep;
    }
  ostr <<   "\n" << spc;
  for(i=0; i<nc; i++)
    {
    ostr << BU << csep << MU;
    if(i<nc-1) ostr << msep;
    }
  ostr << "\n" << spc;
  std::string oc;
  int ccnt = 0;
  int q, lastq = -1;
  for(int j=0; j<_LS; j++)			// Loop over base operators
    {
    i = SI[j];					// Sorted index to use
    q = BaseSpins[i];
    if(q != lastq && ccnt)
      {
      ostr << "\n" << spc;
      ccnt = 0;
      }
    oc = BaseNames[i];				// Odd column string (name)
    if(an) oc += csep + BaseAltNames[i];	// Add in alt name if wanted
    ostr << oc
         << std::string(nlmax - oc.length(), ' ');	// Output odd column
    ostr << csep				// Output even column
         << Gform("%3i", BaseCoherences[i]);
    ccnt++;					// Increment column count
    if(ccnt == nc)
      {
      ostr << "\n" << spc;
      ccnt = 0;
      }
    else ostr << msep;
    lastq = q;
    }
  ostr << "\n";
  return ostr;
  }


 
        // Input                dec     : Decomposition (this)
        //                      ostr    : An output stream
        // Output               ostr    : The output stream, modified to
        //                                contain the information about
        //                                the decomposition
	// Note				: To simplify the output operators
	//				  they are scaled by a coefficient.
	//				  This is reflected in the printed
	//				  operator names as well

std::ostream& decomp::printOps(std::ostream& ostr, bool bs) const
  {
  std::string hdr = dname + " Base Operators";
  int len = hdr.length(); 
  ostr << std::string(40 - len/2, ' ') << hdr << "\n\n";
  complex z;
  for(int i=0; i<_LS; i++)		// Loop over all basis operators
    {
    z = 1.0/BaseCoefficients[i];			// Op coefficient
    hdr = BaseCoeffNames[i] + " * " + BaseNames[i];	// Op name
    len = hdr.length();
    ostr << std::string(40 - len/2, ' ') << hdr;
    if(!bs) ostr << "\n" << z*(BaseOps[i]).get_matrix() << "\n";
    else    ostr << "\n" << z*BaseOps[i] << "\n";
    }
  ostr << "\n";
  return ostr;
  }


 
        // Input                dec     : Decomposition (this)
        //                      ostr    : An output stream
        // Output               ostr    : The output stream, modified to
        //                                contain the information about
        //                                the decomposition

std::ostream& decomp::print(std::ostream& ostr, const gen_op& Op)
  {
  decompose(Op);				// Decompose (intensities)
  std::string hdr = dname + " Decomposition";	// Header for output
  int len = hdr.length(); 			// Length of the header
  ostr << std::string(40 - len/2,' ') << hdr 	// Output the header
       << std::endl << std::endl;
  int nlmax = MaxOpNameLen();			// Maximum name length
  int anmax = MaxOpAltNameLen();		// Maximum alt name length
  len = 7 + 1 + nlmax + 3 + 1 + anmax;
  std::string spc = std::string(40 - len/2,' ');
  std::string csp = " ";
  for(int i=0; i<_LS; i++)			// Loop over all basis ops
    {						// Output only those with norms
    if(norm(BaseVals.get(i)) > thresh)		// above the specified threshold
      {
      if(!prif) ostr << spc << Gform("%7.3f", BaseVals.getRe(i));
      else      ostr << "  " << BaseVals.get(i);
      ostr << csp << BaseNames[i] << "\t"
           << csp << Gform("%3i", BaseCoherences[i]) ;
      if(anmax) ostr << "csp" << BaseAltNames[i];
      ostr << "\n";
      }
    }
  ostr << "\n";
  return ostr;
  }

std::ostream &operator << (std::ostream& ostr, const decomp& dec)
  { return dec.print(ostr); }


#endif					// Class decomp implementation
