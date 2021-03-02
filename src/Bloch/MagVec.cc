/* MagVec.cc ****************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Magnetization Vector 				Implementation 	**
**                                                                      **
**      Copyright (c) 2002                                              **
**      Dr. Scott A. Smith						**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** The class MagVec defines a Bloch magnetization vector. Such vectors  **
** are meant to evolve in time as described by the phenomenological     **
** Bloch equations. A magnetization vector as specified herein involves **
** N individual sub-vectors, each of these sub-vectors has three        **
** magnetization components. The vector with N components takes the     **
** blocked form                                                         **
**                   |M> = ||M > |M > |M > ..... |M   >                 **
**                            0    1    2          N-1                  **
**                                                                      **
** & each sub-vector is                                                 **
**                        |M > = |M  , M   , M  >                       **
**                          i      ix   iy    iz                        **
**                                                                      **
** The components of |M> are then given by                              **
**                                                                      **
**           <3*i|M> = M     <3*i+1|M> = M     <3*i+2|M> = M            **
**                      ix                iy                iz          **
**                                                                      **
** a static Bo field, an applied rf-field (B1), simplistic relaxation   **
** (T1 & T2), and possibly exchange between magnetizaiton vectors.      **
**                                                                      **
** The magnetization vector has ONLY 1 type of evolution implicitly     **
** defined, that dictated by the phenomenological Bloch equations.      **
**                                                                      **
**                           -Gt                                        **
**                 |M(t)> = e   |M(0)-M   > + |M   >                    **
**                                     inf      inf                     **
**                                                                      **
** where the evolution matrix G and the initinite time vector M         **
** must be specified appropriately.                            inf      **
**                                                                      **
*************************************************************************/

#ifndef   MagVec_cc_			// Is file already included?
#  define MagVec_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#if defined(_MSC_VER)			// If we are using MSVC++
 #pragma warning (disable : 4786)       // Kill STL namelength warnings
#endif

#include <Bloch/MagVec.h>		// Includes the interface 
#include <Basics/Gutils.h>              // Include GAMMA errors/queries
#include <Basics/StringCut.h>		// Include Gdec function
#include <HSLib/HSLibIF.h>		// Include Hilbert space module
#include <Level2/acquire1D.h>		// Inlcude 1D rapid acquisitions

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                   Magnetization Vector Error Handling
// ____________________________________________________________________________


        // Input                sys     : Spin system (this)
        //                      eidx    : Error index
        //                      noret   : Flag for linefeed (0=linefeed)
	// Output		void	: An error message is output

/* The following error messages use the defaults set in the Gutils package

                Case                          Error Message

		(0) 			Program Aborting.....
    		(3) 			Can't Construct From Parameter Set
		(4)			Cannot Construct From Input File
    		(5)			Cannot Write To Parameter File
    		(6)			Cannot Write To Output FileStream
    		(7)			Cannot Produce Any Output
    		default			Unknown Error                        */

        // Input                sys     : Spin system (this)
        //                      eidx    : Error index
	//			pname	: Additional error message
        //                      noret   : Flag for linefeed (0=linefeed)
        // Output               none    : Error message output

/* The following error messages use the defaults set in the Gutils package

                Case                          Error Message

		(1) 			Problems With File PNAME
		(2)			Cannot Read Parameter PNAME
    		default			Unknown Error - PNAME                */
        // Input                sys     : Spin system (this)
        //                      eidx    : Error index
	//			pname	: Additional error message
        // Output               none    : Error message output
        //                                Program execution stopped


void MagVec::MVerror(int eidx, int noret) const
  {
  std::string hdr("Magnetization Vector");
  std::string msg;
  switch (eidx)
    {
    case 10: msg = std::string("Dimensioning Mismatch"); 		// (10)
             GAMMAerror(hdr,msg,noret);  break;
    case 15: msg = std::string("Improper Number Of SubVector Norms");// (15)
             GAMMAerror(hdr,msg,noret);  break;
    case 16: msg = std::string("Accessed Subvector Out Of Range");	// (16)
             GAMMAerror(hdr,msg,noret);  break;
    case 20: msg = std::string("Can't Write To Parameter File");	// (20)
             GAMMAerror(hdr,msg,noret);  break;
    case 22: msg = std::string("Can't Write To Output FileStr"); 	// (22)
             GAMMAerror(hdr,msg,noret);  break;
    case 23: msg = std::string("Cannot Output Parameters");		// (23)
             GAMMAerror(hdr,msg,noret);  break;
    case 27: msg = std::string("Negative SubVector Norm Specified");	// (27)
             GAMMAerror(hdr,msg,noret);  break;
    case 35: msg = std::string("Cant Set SubVector Norms");		// (35)
             GAMMAerror(hdr,msg,noret);  break;
    case 40: msg = std::string("Parameters Insufficient");		// (40)
             GAMMAerror(hdr,msg,noret);  break;
    case 54: msg = std::string("Cant Determine Sub-Vector Cmpnts."); // (54)
             GAMMAerror(hdr,msg,noret);  break;
    case 56: msg = std::string("Cant Determine # Sub-Vectors");	// (56)
             GAMMAerror(hdr,msg,noret);  break;
    case 60: msg = std::string("Cant Determine X Magnetization");	// (60)
             GAMMAerror(hdr,msg,noret);  break;
    case 61: msg = std::string("Cant Determine Y Magnetization");	// (61)
             GAMMAerror(hdr,msg,noret);  break;
    case 62: msg = std::string("Cant Determine Z Magnetization");	// (62)
             GAMMAerror(hdr,msg,noret);  break;
    case 63: msg = std::string("Cant Determine All Sub-Vectors");	// (63)
             GAMMAerror(hdr,msg,noret);  break;
    default: GAMMAerror(hdr, eidx, noret); break;		// (-1)
    }
  }

void MagVec::MVerror(int eidx, const std::string& pname, int noret) const
  {
  std::string hdr("Magnetization Vector");
  std::string msg;
  switch(eidx)
    {
    case 101:                                                         // (101)
      msg = std::string("Can't Find Parameters For ")
          + pname;
     GAMMAerror(hdr, msg, noret); break;
    case 102:                                                         // (102)
      msg = std::string("Can't Find ") + pname 
          + std::string(" In Parameter Set");
     GAMMAerror(hdr, msg, noret); break;
    case 103:                                                         // (103)
      msg = std::string("Can't Find Components For Sub-Vector ") + pname;
     GAMMAerror(hdr, msg, noret); break;
    case 104:                                                         // (104)
      msg = std::string("Can't Get Sub-Vector ") + pname + " Norm";
     GAMMAerror(hdr, msg, noret); break;
    case 105:                                                         // (105)
      msg = std::string("Can't Set Sub-Vector ") + pname + " Norm";
     GAMMAerror(hdr, msg, noret); break;
    default: GAMMAerror(hdr, eidx, pname, noret); break;// Unknown Error   (-1)
    }
  }

volatile void MagVec::MVfatal(int eidx) const
  {  
  MVerror(eidx, 1);				// Normal non-fatal error
  if(eidx) MVerror(0);				// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

volatile void MagVec::MVfatal(int eidx, const std::string& pname) const
  {  
  MVerror(eidx, pname, 1);			// Normal non-fatal error
  if(eidx) MVerror(0);				// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii                   Magnetization Vector Checking Functions
// ____________________________________________________________________________


bool MagVec::CheckNorms(const std::vector<double>& norms, bool warn) const
  {
  int nn = norms.size();			// Number of norms
  int nc = NComps();				// Number of sub-vectors
  if(nn != nc)					// Insure number of norms
    {						// is proper for this vector
    if(warn)
      {
      MVerror(10, 1);				// Dimensioning mismatch
      MVerror(15, 1);				// Improper # of norms
      }
    return false;
    }
  for(int i=0; i<nn; i++)			// Insure no norms are
    if(norms[i] < 0)				// negative
      {
      if(warn) MVerror(27, 1);			// Negative norm set
      return false;
      }
  return true;
  }

bool MagVec::CheckRange(int cmp, bool warn) const
  {
  int nc = NComps();				// Number of sub-vectors
  if(cmp >= nc)
    {
    if(warn)
      {
      MVerror(10, 1);				// Dimensioning problem
      MVerror(16, 1);				// Subvector out of range
      }
    return false;
    }
  return true;
  }

// ____________________________________________________________________________
// iii                Magnetization Vector Setup Functions
// ____________________________________________________________________________

/* These are protected functions because they allow specific aspects of the
   vector to be set up without worrying about consistency! Basically what we
   must determine are

                          1.) Number of sub-vectors
                          2.) Mx, My, and Mz for each sub-vector

   Parameter NMagVecs is used indicate how many sub-vectors are defined.
   If that is NOT found, it is assumed that there is only 1 sub-vector.     */

bool MagVec::SetVector(const ParameterSet& pset, int pfx, bool warn)
  {
  ParameterSet  subpset;                        // Working parameter set
  if(pfx != -1) subpset = pset.strip(pfx);      // Get only params with [#]
  else          subpset = pset;                 // Or use full pset

  int nsv;					// For # of sub-vectors
  if(!GetNVects(subpset, nsv, false))		// If # sub-vectors not found
    {						// assume that theres only 1
    double Mx, My, Mz;				//   Magnetization components
    if(GetMxMyMz(subpset,Mx,My,Mz,-1,warn))	//   Look for Mx, My, Mz
      {						//     If found, we know M
      *this = MagVec(1);			//     Set for 1 sub-vector
      put(Mx, 0);				//     Set Mx at <0|M>
      put(My, 1);				//     Set My at <1|M>
      put(Mz, 2);				//     Set Mz at <2|M>
      return true;				//     All finished
      }
    else					//   No MagVects, No Mx,My,Mz
      {						//   so don't know how we work
      if(warn) MVerror(56, 1);			//     Warn no # mag vects
      if(warn) MVerror(40, 1);			//     Warn not enough params
      return false;				//     Return we failed
      }
    }
  if(SetSubVects(subpset,nsv,warn)) return true;// Set up nsv sub-vectors

  if(warn)					// If here, we cannot find
    {						// all of the sub-vectors
    MVerror(63, 1);				//   Cannot get all sub-vectors
    MVerror(40, 1);				//   Warn insufficient params
    }
  return false;					//   Return we failed
  }

bool MagVec::GetNVects(const ParameterSet& pset, int& nm, bool warn) const
  {
  std::string pstate, pname("NMagVecs");		// Parameter name
  ParameterSet::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);                      // Pix in param. list for NSpins
  if(item != pset.end())                        // If it was found in the list
    {
    (*item).parse(pname,nm,pstate);		//   Retrieve # of spins
    return true;				//   Return we found it
    }
  if(warn) MVerror(102, pname, 1);
  return false;
  }

bool MagVec::SetSubVects(const ParameterSet& pset, int N, bool warn)
  {
  MagVec MV(N);					// New vector
  coord   pt;					// For sub-vect components
  double  Mx, My, Mz;				// For sub-vect components
  for(int i=0; i<N; i++)			// Loop over sub-vectors
    {						// and try & find Mx,My,Mz
    if(GetCoord(pset, pt, i, false))		//   If in coordinate form
      {						//   then parse out the
      Mx = pt.x();				//   three values from the
      My = pt.y();				//   coordinate
      Mz = pt.z();
      }						//   If not coordinate form
    else if(!GetMxMyMz(pset,Mx,My,Mz,i,false)) 	//   try for Mx,My,Mz
      {
      if(warn)					//   If not coordinate or
        {					//   Mx, My, Mz we failed
        MVerror(103, Gdec(i), 1);		//     No components found 
        }
      return false;
      }
    MV.put(Mx, 3*i);				// Set Mx ith subvector
    MV.put(My, 3*i+1);				// Set My ith subvector
    MV.put(Mz, 3*i+2);				// Set Mz ith subvector
    }
  *this = MV;
  return true;					// Set all N sub-vecotrs
  }

bool MagVec::GetCoord(const ParameterSet& pset, coord& Pt, int i, bool warn) const
  {
  std::string sufx;					// Name suffix
  if(i != -1) sufx = "(" + Gdec(i) + ")";	// Use suffix if idx != -1
  std::string pname = "MCoord" + sufx; 		// Parameter name
  ParameterSet::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);                      // Pix in list for MCoord(i)
  std::string pstate, pname1;                        // Temp string
  if(item != pset.end())                        // Retrieve the offset value
    {
    Pt = coord(*item);
    return true;
    }
  pname1 = "Coord(" + Gdec(i) + ")";            // Parameter name
  item = pset.seek(pname1);                     // Pix in list for MCoord(i)
  if(item != pset.end())                        // Retrieve the coordinate
    {
    Pt = coord(*item);
    return true;
    }
  if(warn)
    {
    MVerror(54, 1);
    MVerror(102, pname, 1);
    }
  Pt = UnitZ;
  return false;
  }

bool MagVec::GetMxMyMz(const ParameterSet& pset,
                   double& Mx, double& My, double& Mz, int idx, bool warn) const
  {
  std::string sufx;					// Name suffix
  if(idx != -1) sufx = "(" + Gdec(idx) + ")";	// Use suffix if idx != -1
  std::string pname = "Mx" + sufx;			// Parameter name
  ParameterSet::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);                      // Pix in list for Mx(i)
  if(item == pset.end())                        // If Mx(i) not is pset, we
    {						// have failed miserably
    if(warn) 					//   Issue warnings if desired
      { MVerror(60,1); MVerror(102,pname,1); }
    return false;				//   Return our failure
    }
  std::string pstate; 				// Temp string
  (*item).parse(pname,Mx,pstate);		// Retrieve Mx value
  pname = "My" + sufx;				// Parameter name
  item = pset.seek(pname);                      // Pix in list for My(i)
  if(item == pset.end())                        // If My(i) not is pset, we
    {						// have failed miserably
    if(warn) 					//   Issue warnings if desired
      { MVerror(61,1); MVerror(102,pname,1); }
    return false;				//   Return our failure
    }
  (*item).parse(pname,My,pstate);		// Retrieve My value
  pname = "Mz" + sufx;				// Parameter name
  item = pset.seek(pname);                      // Pix in list for Mz(i)
  if(item == pset.end())                        // If Mz(i) not is pset, we
    {						// have failed miserably
    if(warn) 					//   Issue warnings if desired
      { MVerror(62,1); MVerror(102,pname,1); }
    return false;				//   Return our failure
    }
  (*item).parse(pname,Mz,pstate);		// Retrieve My value
  return true;					// We found Mx, My and Mz
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A               MAGNETIZATION VECTOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                            Simple Constructors
// ----------------------------------------------------------------------------

MagVec::MagVec(int N) : col_vector(3*N)
  {							// For N sub-vectors
  for(int i=0; i<N; i++)				// Fill each sub-vector
    {							// with {0,0,1} by
    put(0,3*i);						// Default
    put(0,3*i+1);
    put(1,3*i+2);
    }
  }

MagVec::MagVec(const MagVec& M)      : col_vector(M)  {}
MagVec::MagVec(const col_vector& CV) : col_vector(CV) {}

// ----------------------------------------------------------------------------
//              Constructors Using Individual Magnetization Vectors
// ----------------------------------------------------------------------------

MagVec::MagVec(double Mx, double My, double Mz) : col_vector(3)
  { put(Mx, 0); put(My, 1); put(Mz, 2); }

MagVec::MagVec(const coord& M) : col_vector(3)
  { put(M.x(), 0); put(M.y(), 1); put(M.z(), 2); }

MagVec::MagVec(double Mx1,double My1,double Mz1,double Mx2,double My2,double Mz2)
       : col_vector(6)
  { put(Mx1,0); put(My1,1); put(Mz1,2); put(Mx2,3); put(My2,4); put(Mz2,5); }

MagVec::MagVec(const coord& M1, const coord& M2) : col_vector(6)
  { 
  put(M1.x(),0); put(M1.y(),1); put(M1.z(),2);
  put(M2.x(),3); put(M2.y(),4); put(M2.z(),5); 
  }

MagVec::MagVec(const std::vector<coord>& Ms) : col_vector(3*Ms.size())
  { 
  int N = Ms.size();
  for(int i=0; i<N; i++)
    {	
    put(Ms[i].x(), 3*i);
    put(Ms[i].y(), 3*i+1);
    put(Ms[i].z(), 3*i+2);
    }
  }

// ----------------------------------------------------------------------------
//                        Assignment and Destruction
// ----------------------------------------------------------------------------

// MagVec& MagVec::operator= (const MagVec& MV)		// Inherited
// MagVec::~MagVec () 					// Inherited

// ____________________________________________________________________________
// B           Magnetization Vector - Magnetization Vector Interactions
// ____________________________________________________________________________

/* These functions allow for simple mathematical operations between two Bloch
   magnetization vectors.  This includes addition & subtraction. There is
   one unary function as well, negation.
 
   Operator Arguments      Result        Operator Arguments    Result
   -------- --------- -----------------  -------- --------- -------------------   
      -        M      Returns -M            -=      M,M1     M1 subt. from M
      +      M,M1     Returns M+M1          +=      M,M1     M1 added to M
      -      M1,M2    Returns M1-M2                                          */



//MagVec MagVec::operator- () const { return MagVec(col_vector::operator- ()); }

MagVec MagVec::operator+ (const MagVec& M1) const
  { MagVec M(*this); M += M1; return M;}

MagVec& MagVec::operator+= (const MagVec& M1)
  { col_vector::operator+= (M1); return (*this); }

MagVec MagVec::operator- (const MagVec &M1) const
  { MagVec M(*this); M -= M1; return M;}

MagVec& MagVec::operator-= (const MagVec& M1)
  { col_vector::operator-= (M1); return (*this); }

// ____________________________________________________________________________
// B                     Magnetization Vector Access
// ____________________________________________________________________________

//int MagVec::size()					// Inherited
int MagVec::NComps() const { return size()/3; }

// ----------------------------------------------------------------------------
//                  Magnetization (Sub)Vector Components
// ----------------------------------------------------------------------------

double MagVec::Mx(int cmp) const { return getRe(3*cmp);   }
double MagVec::My(int cmp) const { return getRe(3*cmp+1); }
double MagVec::Mz(int cmp) const { return getRe(3*cmp+2); }

void MagVec::Mx(int cmp, double mx) { put(mx, 3*cmp);   }
void MagVec::My(int cmp, double my) { put(my, 3*cmp+1); }
void MagVec::Mz(int cmp, double mz) { put(mz, 3*cmp+2); }

double MagVec::x(int    cmp) const { return getRe(3*cmp);   }
double MagVec::y(int    cmp) const { return getRe(3*cmp+1); }
double MagVec::z(int    cmp) const { return getRe(3*cmp+2); }

double MagVec::norm(int cmp) const
  {
  double MX = x(cmp);
  double MY = y(cmp);
  double MZ = z(cmp);
  return sqrt(MX*MX + MY*MY + MZ*MZ);
  }

double MagVec::theta(int cmp) const
  { coord pt(x(cmp), y(cmp), z(cmp)); return pt.theta(); }

double MagVec::phi(int cmp) const
  { coord pt(x(cmp), y(cmp), z(cmp)); return pt.phi(); }


// ----------------------------------------------------------------------------
//                    Magnetization (Sub)Vector Norms
// ----------------------------------------------------------------------------

std::vector<double> MagVec::Norms() const
  {
  std::vector<double> Nrms; 
  int nm = NComps();
  for(int i=0; i<nm; i++)
    Nrms.push_back(Norm(i));
  return Nrms;
  }

void MagVec::Norms(const std::vector<double>& Ns)
  {
  if(!CheckNorms(Ns))  MVfatal(35);		// Check norm count matches
  int nm = NComps();
  for(int i=0; i<nm; i++)
    Norm(Ns[i], i);
  }

double MagVec::Norm(int cmp) const
  {
  if(!CheckRange(cmp)) MVfatal(104, Gdec(cmp));
  int     I = 3*cmp;				// Index of component
  double Mx = getRe(I);				// Get Mx of component cmp
  double My = getRe(I+1);			// Get My of component cmp
  double Mz = getRe(I+2);			// Get Mz of component cmp
  return sqrt(Mx*Mx + My*My + Mz*Mz);		// Return norm of component
  }

void MagVec::Norm(double nv, int cmp) 
  {
  if(!CheckRange(cmp)) MVfatal(105, Gdec(cmp));
  double sf =  nv/Norm(cmp);			// Scaling to adjust norm
  int I = 3*cmp;				// Index of component
  put(get(I)*sf, I);				// Set Mx of component cmp
  put(get(I+1)*sf, I+1);			// Set My of component cmp
  put(get(I+2)*sf, I+2);			// Set Mz of component cmp
  }

// ____________________________________________________________________________
// C              Magnetization Vector Parameter Set Functions
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//       Functions To Make A Magnetization Vector From A Parameter Set
// ----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                   Parameter Set From Magnetization Vector
//-----------------------------------------------------------------------------

MagVec::operator ParameterSet( ) const
  { ParameterSet pset; pset += *this; return pset; }

void operator+= (ParameterSet& pset, const MagVec& MV)
  { MV.PSetAdd(pset); }

bool MagVec::PSetAdd(ParameterSet& pset, int pfx) const
  {
  std::string prefx;						// Parameter prefix
  if(pfx != -1)						// Use prefix if pfx
    prefx = std::string("[")+Gdec(pfx)+std::string("]");		// is NOT -1
  int NC = NComps();					// Get # subvectors
  std::string pname = prefx+std::string("NMagVecs");		// # of subvectors
  std::string pstate("No. of Magnetization SubVectors");	// Parameter statement
  SinglePar par(pname, NC, pstate);			// # subvects parameter
  pset.push_back(par);					// Add param to pset
  coord MSV;						// Working coordinate
  for(int i=0; i<NC; i++)				// Loop sub-vectors
    {
    pname = prefx + "MCoord(" + Gdec(i) + ")";		//   Coordinate name
    pstate = std::string("Magnetization SubVector Comps");	//   Coordinate state
    MSV = coord(getRe(3*i), getRe(3*i+1), getRe(3*i+2));//   Coordiante
    par = MSV.param(pname, pstate);			//   Coordinate param
    pset.push_back(par);				//   Add parameter
    }
  return true;
  }
 
//-----------------------------------------------------------------------------
//            Parameter Set File From Magnetization Vector
//-----------------------------------------------------------------------------

        // Input                M       : Magnetization vector (this)
        //                      filename: Output file name
        //                      ofstr   : Output file stream
        //                      pfx     : Vector prefix (default -1)
        //                      warn    : Warning level
        // Output               none    : Vector is written as a parameter
        //                                to file or output file stream

bool MagVec::write(const std::string& filename, int pfx, int warn) const
  {
  std::ofstream ofstr(filename.c_str());     // Open filename for input
  if(!write(ofstr, pfx, warn?1:0))      // If file bad then exit
    {
    MVerror(40, filename, 1);           // Filename problems
    if(warn>1) MVfatal(20);              // Fatal error
    return false;
    }
  ofstr.close();                        // Close it now
  return true;
  }

bool MagVec::write(std::ofstream& ofstr, int pfx, int warn) const
  {
  ParameterSet pset;                    // Declare a parameter set
  PSetAdd(pset, pfx);                   // Add in interaction parameters
  if(!pset.write(ofstr, warn?1:0))      // Use parameter set to write
    {
    if(warn)
      {
      MVerror(22, 1);                    // Problems writing to filestream
      if (warn>1) MVfatal(23);           // Fatal error
      }
    return false;
    }
  return true;
  }

// ____________________________________________________________________________
// D                  Magnetization Vector Input Functions
// ____________________________________________________________________________

bool MagVec::read(const std::string& filein, int idx, int warn)
  {
  ParameterSet pset;                    // Declare a parameter set
  if(!pset.read(filein, warn?1:0))      // Read in pset from file
    {
    if(warn)
      {
      MVerror(40, filein, 1);		// Problems with file filein
      if(warn>1) MVfatal(101,filein);	// Can't read from file filein
      else       MVerror(101,1);
      }
    return false;
    }
  return read(pset, idx, warn);
  }

bool MagVec::read(const ParameterSet& pset,int idx, int warn)
  {
  if(SetVector(pset, idx, warn?true:false))	// Try and set the vector
    return true;
  if(warn)                                      // Looks like we can't read
    {                                           // the system at all
    std::string sl("");				//   String for index
    if(idx != -1) sl = std::string(Gdec(idx));	//   Set it if one exists
                  MVerror(77, sl, 1);
    if(warn > 1)  MVfatal(5);
    else          MVerror(78, sl);
    }
  return false;
  }

std::string MagVec::ask_read(int argc, char* argv[], int argn)
  {
  std::string filename;                             	 // Name of spin system file
  std::string M("\n\tMagnetization vector filename? ");// Query we will ask if
  query_parameter(argc, argv, argn,M,filename);	// Get/ask for filename
  read(filename);                               // Read system from filename
  return filename;				// Return filename used
  }

std::string MagVec::ask_read(int argc, char* argv[], int argn, const std::string& def)
  {
  std::string M = "\n\tMagnetization vector filename"// Query we will ask if
           + std::string(" [") + def + std::string("]? ");// it is needed
  std::string filename = def;                        // Name of spin system file
  ask_set(argc,argv,argn,M,filename);		// or ask for it
  read(filename);                               // Read system from filename
  return filename;                              // Return filename used
  }
 
// ____________________________________________________________________________
// E                      SYSTEM STANDARD I/O FUNCTIONS
// ____________________________________________________________________________

        // Input                MV      : A magnetization vector (this)
        //                      ostr    : Output stream
        // Output               none    : Vector info is sent
        //                                into the output stream

std::vector<std::string> MagVec::printStrings() const
  {
  std::vector<std::string> PStrings;		// Strings we will return
  int NS = NComps();				// Get number components
  std::string hdr("Magnetization Vector");	// Output header
  int hl;					// Header length
  if(!NS)					// If there are no vectors
    {						// just write system is empty
    hdr = "Empty " + hdr;
    hl = hdr.length();
    hdr = std::string(40-hl/2, ' ') + hdr;
    PStrings.push_back(hdr);
    return PStrings;
    }

/*          Output Some Printed Lines Which Will Look Like The Following

                                 Magnetization Vector

                     |Mi>     Mx       My       Mz    ||Mi>| 
                    ------ -------- -------- ------- -------- 
                     ...     ...      ...      ...      ...
                     ...     ...      ...      ...      ...
                     ...     ...      ...      ...      ...
                     ...     ...      ...      ...      ...
                     ...     ...      ...      ...      ...                  */

//                          This Outputs The Main Header

  hl = hdr.length();				// Main header length
  PStrings.push_back("\n\n");			// Start with new lines
  hdr = std::string(40-hl/2,' ') + hdr;		// Output main header
  PStrings.push_back(hdr);			// Store main header
  PStrings.push_back("");

//                            Set Up Output Columns

  std::string MXC("    Mx    ");			// 2nd column header
  std::string MYC("    My    ");			// 3rd column header
  std::string MZC("    Mz    ");			// 4th column header
  std::string MNC("   |Mi>|  ");			// 6th column header

  int    cgap = 1;					// Gap between columns
  std::string cspc = std::string(cgap, ' ');		// String for gap

  std::vector<std::string> chdrs;			// Column headers
  chdrs.push_back(std::string("|Mi>"));			// 1st column header
  chdrs.push_back(MXC);					// 2nd column header
  chdrs.push_back(MYC);					// 3rd column header
  chdrs.push_back(MZC);					// 4th column header
  chdrs.push_back(MNC);					// 5th column header
  int nc = chdrs.size();				// Number of columns

  int vlen = Gdec(NS).length();				// Length # of components
  int c1l = chdrs[0].length();				// Current col width
  std::string vst, vfi;					// To center row index
  int x;
  if(c1l > vlen)					// Need to center if
    {							// column header wider
    x = (c1l-vlen)/2; 					// than largest index
    if(x) vst = std::string(x, ' ');
    x = c1l-vlen-vst.length();
    if(x) vfi = std::string(x, ' ');
    }
  else if(c1l < vlen)					// Need to widen col
    {							// header if #s wider
    x = (vlen-c1l)/2; 					// than largest index
    if(x) chdrs[0] = std::string(x, ' ') + chdrs[0];
    x = vlen-chdrs[0].length();
    if(x) chdrs[0] = chdrs[0] + std::string(x, ' ');
    }
    
  int i, llen=0;					// Line length
  for(i=0; i<nc-1; i++)					// Determine line length
    llen += chdrs[i].size() + cgap;
  llen += chdrs[i].size();

  std::string spc = std::string(40-llen/2, ' ');	// Spacer to center line

//                         Output Column Names
  
  std::string pline = spc;				// Start of header line
  for (i=0; i<nc; i++)
    {
    pline += chdrs[i];
    if(i<nc-1) pline += cspc;
    }
  PStrings.push_back(pline);

//                      Output Column Name Underlines

  pline = spc;					// Start header underline
  for (i=0; i<nc; i++)
    {
    pline += std::string(chdrs[i].length(), '-');
    if(i<nc-1) pline += cspc;
    }
  PStrings.push_back(pline);


  for(i=0; i<NS; i++)				// Loop all sub-vectors
    {
    pline = spc;
    pline += vst + Gdec(i, vlen) + vfi + cspc;
    pline += Gform("%10.4f", getRe(3*i))   + cspc;
    pline += Gform("%10.4f", getRe(3*i+1)) + cspc;
    pline += Gform("%10.4f", getRe(3*i+2)) + cspc;
    pline += Gform("%10.4f", Norm(i))      + cspc;
    PStrings.push_back(pline);
    }

  return PStrings;
  }

std::ostream& MagVec::print(std::ostream& ostr, int np) const
  {
  int NS = NComps();				// Get number components
  std::string hdr("Magnetization Vector");	// Output header
  int hl;					// Header length
  if(!NS)					// If there are no vectors
    {						// just write system is empty
    hdr = "Empty " + hdr;
    hl = hdr.length();
    ostr << "\n" << std::string(40-hl/2, ' ') << hdr;
    return ostr;
    }

/*          Output Some Printed Lines Which Will Look Like The Following

                                 Magnetization Vector

                     |Mi>     Mx       My       Mz    ||Mi>| 
                    ------ -------- -------- ------- -------- 
                     ...     ...      ...      ...      ...
                     ...     ...      ...      ...      ...
                     ...     ...      ...      ...      ...
                     ...     ...      ...      ...      ...
                     ...     ...      ...      ...      ...                  */

//                          This Outputs The Main Header

  hl = hdr.length();				// Main header length
  ostr << "\n\n";				// Begin with new lines
  ostr << std::string(40-hl/2,' ') << hdr;		// Output main header
  if(NS > np)
    {
    hdr = "(" + Gdec(np) + " Components of " + Gdec(NS) + " Output)";
    hl = hdr.length();				// Sub header length
    ostr << "\n" << std::string(40-hl/2,' ') << hdr;	// Output sub header
    NS = np;					// Reset to min comps
    }

//                            Set Up Output Columns

  std::string MXC("    Mx    ");			// 2nd column header
  std::string MYC("    My    ");			// 3rd column header
  std::string MZC("    Mz    ");			// 4th column header
  std::string MNC("   |Mi>|  ");			// 6th column header

  int    cgap = 1;                              // Gap between columns
  std::string cspc = std::string(cgap, ' ');              // String for gap

  std::vector<std::string> chdrs;                         // Column headers
  chdrs.push_back(std::string("|Mi>"));		// 1st column header
  chdrs.push_back(MXC);				// 2nd column header
  chdrs.push_back(MYC);				// 3rd column header
  chdrs.push_back(MZC);				// 4th column header
  chdrs.push_back(MNC);				// 5th column header
  int nc = chdrs.size();			// Number of columns

  int x = NS;					// Copy # of components
  int vlen = 1;					// Spacing on 1st column
  while(x > 10) { vlen++; x/=10; }		// Space depends on #
  int c1l = chdrs[0].length();			// Current col width
  std::string vst, vfi;				// To center row index
  if(c1l > vlen)				// Need to center if
    {						// column header wider
    x = (c1l-vlen)/2; 				// than largest index
    if(x) vst = std::string(x, ' ');
    x = c1l-vlen-vst.length();
    if(x) vfi = std::string(x, ' ');
    }
  else if(c1l < vlen)				// Need to widen col
    {						// header if #s wider
    x = (vlen-c1l)/2; 				// than largest index
    if(x) chdrs[0] = std::string(x, ' ') + chdrs[0];
    x = vlen-chdrs[0].length();
    if(x) chdrs[0] = chdrs[0] + std::string(x, ' ');
    }
    
  int i, llen=0;				// Line length
  for(i=0; i<nc-1; i++)				// Determine line length
    llen += chdrs[i].size() + cgap;
  llen += chdrs[i].size();

  std::string spc = std::string(40-llen/2, ' ');			// Spacer to center line

//                         Output Column Names

  ostr << "\n\n" << spc;			// Start of header line
  for (i=0; i<nc; i++)
    {
    ostr << chdrs[i];
    if(i<nc-1) ostr << cspc;
    }

//                      Output Column Name Underlines

  ostr << "\n" << spc;				// Start header underline
  for (i=0; i<nc; i++)
    {
    ostr << std::string(chdrs[i].length(), '-');
    if(i<nc-1) ostr << cspc;
    }

  for(i=0; i<NS; i++)				// Loop all sub-vectors
    {
    ostr << "\n" << spc;
    ostr << vst << Gdec(i, vlen) << vfi   << cspc;
    ostr << Gform("%10.4f", getRe(3*i))   << cspc;
    ostr << Gform("%10.4f", getRe(3*i+1)) << cspc;
    ostr << Gform("%10.4f", getRe(3*i+2)) << cspc;
    ostr << Gform("%10.4f", Norm(i))      << cspc;
    }
  return ostr;
  }

std::ostream& operator<< (std::ostream& ostr, const MagVec& MV)
  { return MV.print(ostr); }

// ____________________________________________________________________________
// F                  Specialized Magnetization Vectors
// ____________________________________________________________________________


MagVec MagVec::Mx() const
  {
  MagVec M(*this);
  for(int i=0; i<NComps(); i++)
    {
    M.put(1.0, 3*i);
    M.put(0.0, 3*i+2);
    M.put(0.0, 3*i+3);
    }
  return M;
  }

MagVec MagVec::My() const
  {
  MagVec M(*this);
  for(int i=0; i<NComps(); i++)
    {
    M.put(0.0, 3*i);
    M.put(1.0, 3*i+2);
    M.put(0.0, 3*i+3);
    }
  return M;
  }

MagVec MagVec::Mz() const
  {
  MagVec M(*this);
  for(int i=0; i<NComps(); i++)
    {
    M.put(0.0, 3*i);
    M.put(0.0, 3*i+2);
    M.put(1.0, 3*i+3);
    }
  return M;
  }

MagVec MagVec::MxVec(int NC)
  {
  MagVec M(NC);
  for(int i=0; i<M.NComps(); i++)
    { M.put(1.0, 3*i); M.put(0.0, 3*i+2); M.put(0.0, 3*i+3); }
  return M;
  }

MagVec MagVec::MyVec(int NC)
  {
  MagVec M(NC);
  for(int i=0; i<M.NComps(); i++)
    { M.put(0.0, 3*i); M.put(1.0, 3*i+2); M.put(0.0, 3*i+3); }
  return M;
  }

MagVec MagVec::MzVec(int NC)
  {
  MagVec M(NC);
  for(int i=0; i<M.NComps(); i++)
    { M.put(0.0, 3*i); M.put(0.0, 3*i+2); M.put(1.0, 3*i+3); }
  return M;
  }

MagVec MagVec::MxVec() const { return Mx(); }
MagVec MagVec::MyVec() const { return My(); }
MagVec MagVec::MzVec() const { return Mz(); }


// ____________________________________________________________________________
// F                Magnetization Vector Evolution Functions
// ____________________________________________________________________________
    
/* These functions are decidedly Bloch based. That is, they are set to evolve 
   the magnetization vector according to the phenomenological Bloch equations.
   The user must supply a 3
*/

#endif								// MagVec.cc
