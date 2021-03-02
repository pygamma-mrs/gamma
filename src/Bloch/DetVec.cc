/* DetVec.cc ****************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Bloch Detection Vector 				Implementation 	**
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
** The class DetVec defines a Bloch detection vector. Such vectors      **
** are used when determining detectable magnetization values from       **
** magnetization vectors that evolve according to the phenomenological  **
** Bloch equations. A magnetization vector as specified herein involves **
** N individual sub-vectors, each of these sub-vectors has three        **
** magnetization components. The vector with N components takes the     **
** blocked form                                                         **
**                   |M> = ||M > |M > |M > ..... |M   >>                **
**                            0    1    2          N-1                  **
**                                                                      **
** & each sub-vector is                                                 **
**                        |M > = |M  , M   , M  >                       **
**                          i      ix   iy    iz                        **
**                                                                      **
** A Bloch detection vector is a vector whcih selects of various        **
** components from the magnetization vector. It will (should) have a    **
** corresponding blocked structure.                                     **
**                                                                      **
**                   <D| = <<D | <D | <D | ..... <D   ||                **
**                            0    1    2          N-1                  **
**                                                                      **
** As an example one might consider detection of x magnetization. The   **
** compoennts of |M> are generally given by                             **
**                                                                      **
**           <3*i|M> = M     <3*i+1|M> = M     <3*i+2|M> = M            **
**                      ix                iy                iz          **
**                                                                      **
** so the Bloch detection vector will have components                   **
**                                                                      **
**        <D |3*i> = 1     <D |3*i+1> = 0     <D |3*i+2> = 0            **
**          x                y                  z                       **
**                                                                      **
** Bloch detection vectors are simply GAMMA row vectors set to interact **
** with Bloch magnetization vectors (MagVec). They add additional       **
** functionality beyond row vectors (row_vector) that make simulations  **
** involving the phenomenological Bloch equations integrate smoothly    **
** into the rest of the platform.                                       **
**                                                                      **
*************************************************************************/

#ifndef   DetVec_cc_			// Is file already included?
#  define DetVec_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Bloch/DetVec.h>		// Includes the interface 
#include <Basics/Gutils.h>              // Include GAMMA errors/queries
#include <Basics/StringCut.h>		// Include Gdec function
#include <HSLib/HSLibIF.h>		// Include Hilbert space module
#include <Level2/acquire1D.h>		// Inlcude 1D rapid acquisitions
#include <Level1/coord.h>		// Include GAMMA coordinates

#ifdef _MSC_VER                         // If we are using MSVC++
 #pragma warning (disable : 4786)       //   Kill STL namelength warnings
#endif

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                   Bloch Detection Vector Error Handling
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


void DetVec::BDVerror(int eidx, int noret) const
  {
  std::string hdr("Bloch Detection Vector");
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
    case 60: msg = std::string("Cant Determine X Detection");	// (60)
             GAMMAerror(hdr,msg,noret);  break;
    case 61: msg = std::string("Cant Determine Y Detection");	// (61)
             GAMMAerror(hdr,msg,noret);  break;
    case 62: msg = std::string("Cant Determine Z Detection");	// (62)
             GAMMAerror(hdr,msg,noret);  break;
    case 63: msg = std::string("Cant Determine All Sub-Vectors");	// (63)
             GAMMAerror(hdr,msg,noret);  break;
    default: GAMMAerror(hdr, eidx, noret); break;		// (-1)
    }
  }

void DetVec::BDVerror(int eidx, const std::string& pname, int noret) const
  {
  std::string hdr("Bloch Detection Vector");
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

volatile void DetVec::BDVfatal(int eidx) const
  {  
  BDVerror(eidx, 1);				// Normal non-fatal error
  if(eidx) BDVerror(0);				// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

volatile void DetVec::BDVfatal(int eidx, const std::string& pname) const
  {  
  BDVerror(eidx, pname, 1);			// Normal non-fatal error
  if(eidx) BDVerror(0);				// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii                   Bloch Detection Vector Checking Functions
// ____________________________________________________________________________


bool DetVec::CheckNorms(const std::vector<double>& norms, bool warn) const
  {
  int nn = norms.size();			// Number of norms
  int nc = NComps();				// Number of sub-vectors
  if(nn != nc)					// Insure number of norms
    {						// is proper for this vector
    if(warn)
      {
      BDVerror(10, 1);				// Dimensioning mismatch
      BDVerror(15, 1);				// Improper # of norms
      }
    return false;
    }
  for(int i=0; i<nn; i++)			// Insure no norms are
    if(norms[i] < 0)				// negative
      {
      if(warn) BDVerror(27, 1);			// Negative norm set
      return false;
      }
  return true;
  }

bool DetVec::CheckRange(int cmp, bool warn) const
  {
  int nc = NComps();				// Number of sub-vectors
  if(cmp >= nc)
    {
    if(warn)
      {
      BDVerror(10, 1);				// Dimensioning problem
      BDVerror(16, 1);				// Subvector out of range
      }
    return false;
    }
  return true;
  }

// ____________________________________________________________________________
// iii                Bloch Detection Vector Setup Functions
// ____________________________________________________________________________

/* These are protected functions because they allow specific aspects of the
   vector to be set up without worrying about consistency! Basically what we
   must determine are

                          1.) Number of sub-vectors
                          2.) Which sub-vectors are detected
                          3.) Which components are detected

   Parameter NMagVecs is used indicate how many sub-vectors are defined.
   If that is NOT found, it is assumed that there is only 1 sub-vector.     */

bool DetVec::SetVector(const ParameterSet& pset, int pfx, bool warn)
  {
  ParameterSet  subpset;                        // Working parameter set
  if(pfx != -1) subpset = pset.strip(pfx);      // Get only params with [#]
  else          subpset = pset;                 // Or use full pset

  int nsv;					// For # of sub-vectors
  if(!GetNVects(subpset, nsv, false))		// If # sub-vectors not found
    {						// assume that theres only 1
    double Dx, Dy, Dz;				//   Detection components
    if(GetMxMyMz(subpset,Dx,Dy,Dz,-1,warn))	//   Look for Dx, Dy, Dz
      {						//     If found, we know D
      *this = DetVec(1);			//     Set for 1 sub-vector
      put(Dx, 0);				//     Set Dx at <0|M>
      put(Dy, 1);				//     Set Dy at <1|M>
      put(Dz, 2);				//     Set Dz at <2|M>
      return true;				//     All finished
      }
    else					//   No MagVects, No Dx,Dy,Dz
      {						//   so don't know how we work
      if(warn) BDVerror(56, 1);			//     Warn no # mag vects
      if(warn) BDVerror(40, 1);			//     Warn not enough params
      return false;				//     Return we failed
      }
    }
  if(SetSubVects(subpset,nsv,warn)) return true;// Set up nsv sub-vectors

  if(warn)					// If here, we cannot find
    {						// all of the sub-vectors
    BDVerror(63, 1);				//   Cannot get all sub-vectors
    BDVerror(40, 1);				//   Warn insufficient params
    }
  return false;					//   Return we failed
  }

bool DetVec::GetNVects(const ParameterSet& pset, int& nm, bool warn) const
  {
  std::string pstate, pname("N");		// Parameter name
  ParameterSet::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);                      // Pix in param. list for NSpins
  if(item != pset.end())                        // If it was found in the list
    {
    (*item).parse(pname,nm,pstate);		//   Retrieve # of spins
    return true;				//   Return we found it
    }
  if(warn) BDVerror(102, pname, 1);
  return false;
  }

bool DetVec::SetSubVects(const ParameterSet& pset, int N, bool warn)
  {
  DetVec BDV(N);					// New vector
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
        BDVerror(103, Gdec(i), 1);		//     No components found 
        }
      return false;
      }
    BDV.put(Mx, 3*i);				// Set Mx ith subvector
    BDV.put(My, 3*i+1);				// Set My ith subvector
    BDV.put(Mz, 3*i+2);				// Set Mz ith subvector
    }
  *this = BDV;
  return true;					// Set all N sub-vecotrs
  }

bool DetVec::GetCoord(const ParameterSet& pset, coord& Pt, int i, bool warn) const
  {
  std::string sufx;					// Name suffix
  if(i != -1) sufx = "(" + Gdec(i) + ")";	// Use suffix if idx != -1
  std::string pname = "MCoord" + sufx; 		// Parameter name
  ParameterSet::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);                      // Pix in list for MCoord(i)
  std::string pstate, pname1;                        // Temp std::string
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
    BDVerror(54, 1);
    BDVerror(102, pname, 1);
    }
  Pt = UnitZ;
  return false;
  }

bool DetVec::GetMxMyMz(const ParameterSet& pset,
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
      { BDVerror(60,1); BDVerror(102,pname,1); }
    return false;				//   Return our failure
    }
  std::string pstate; 				// Temp string
  (*item).parse(pname,Mx,pstate);		// Retrieve Mx value
  pname = "My" + sufx;				// Parameter name
  item = pset.seek(pname);                      // Pix in list for My(i)
  if(item == pset.end())                        // If My(i) not is pset, we
    {						// have failed miserably
    if(warn) 					//   Issue warnings if desired
      { BDVerror(61,1); BDVerror(102,pname,1); }
    return false;				//   Return our failure
    }
  (*item).parse(pname,My,pstate);		// Retrieve My value
  pname = "Mz" + sufx;				// Parameter name
  item = pset.seek(pname);                      // Pix in list for Mz(i)
  if(item == pset.end())                        // If Mz(i) not is pset, we
    {						// have failed miserably
    if(warn) 					//   Issue warnings if desired
      { BDVerror(62,1); BDVerror(102,pname,1); }
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

DetVec::DetVec(int N) : row_vector(3*N)
  {							// For N sub-vectors
  for(int i=0; i<N; i++)				// Fill each sub-vector
    {							// with {1,0,0} by
    put(1,3*i);						// Default (detect Mx)
    put(0,3*i+1);
    put(0,3*i+2);
    }
  }

DetVec::DetVec(const DetVec&     BD) : row_vector(BD) {}
DetVec::DetVec(const row_vector& RV) : row_vector(RV) {}

// ----------------------------------------------------------------------------
//              Constructors Using Individual Bloch Detection Vectors
// ----------------------------------------------------------------------------

DetVec::DetVec(double Dx, double Dy, double Dz) : row_vector(3)
  { put(Dx, 0); put(Dy, 1); put(Dz, 2); }

DetVec::DetVec(const coord& D) : row_vector(3)
  { put(D.x(), 0); put(D.y(), 1); put(D.z(), 2); }

DetVec::DetVec(double Dx1,double Dy1,double Dz1,double Dx2,double Dy2,double Dz2)
       : row_vector(6)
  { put(Dx1,0); put(Dy1,1); put(Dz1,2); put(Dx2,3); put(Dy2,4); put(Dz2,5); }

DetVec::DetVec(const coord& D1, const coord& D2) : row_vector(6)
  { 
  put(D1.x(),0); put(D1.y(),1); put(D1.z(),2);
  put(D2.x(),3); put(D2.y(),4); put(D2.z(),5); 
  }

DetVec::DetVec(const std::vector<coord>& Ds) : row_vector(3*Ds.size())
  { 
  int N = Ds.size();
  for(int i=0; i<N; i++)
    {	
    put(Ds[i].x(), 3*i);
    put(Ds[i].y(), 3*i+1);
    put(Ds[i].z(), 3*i+2);
    }
  }

// ----------------------------------------------------------------------------
//              Constructors Of Common Detection Vectors
// ----------------------------------------------------------------------------


DetVec::DetVec(int N, char x) : row_vector(3*N)
  {
  complex dx=0, dy=0, dz=0;				// 
  switch(x)
    {
    default:
    case 'x': dx=complex1;               break;
    case 'y': dy=complex1;               break;
    case 'z': dz=complex1;               break;
    case '+': dx=complex1; dy= complexi; break;
    case '-': dx=complex1; dy=-complexi; break;
    }
  for(int i=0; i<N; i++)				// Fill each sub-vector
    {							// with {1,0,0} by
    put(dx,3*i);					// Default (detect Mx)
    put(dy,3*i+1);
    put(dz,3*i+2);
    }							// For N sub-vectors
  }							// For N sub-vectors


// ----------------------------------------------------------------------------
//                        Assignment and Destruction
// ----------------------------------------------------------------------------

// DetVec& DetVec::operator= (const DetVec& BDV)		// Inherited
// DetVec::~DetVec () 						// Inherited

// ____________________________________________________________________________
// B                     Bloch Detection Vector Access
// ____________________________________________________________________________

//int DetVec::size()					// Inherited
int DetVec::NComps() const { return size()/3; }

// ----------------------------------------------------------------------------
//                    Detecton (Sub)Vector Components
// ----------------------------------------------------------------------------

double DetVec::Dx(int cmp) const { return getRe(3*cmp);   }
double DetVec::Dy(int cmp) const { return getRe(3*cmp+1); }
double DetVec::Dz(int cmp) const { return getRe(3*cmp+2); }

void DetVec::Dx(int cmp, double dx) { put(dx, 3*cmp);   }
void DetVec::Dy(int cmp, double dy) { put(dy, 3*cmp+1); }
void DetVec::Dz(int cmp, double dz) { put(dz, 3*cmp+2); }

double DetVec::x(int    cmp) const { return getRe(3*cmp);   }
double DetVec::y(int    cmp) const { return getRe(3*cmp+1); }
double DetVec::z(int    cmp) const { return getRe(3*cmp+2); }
double DetVec::norm(int cmp) const
  {
  double DX = x(cmp);
  double DY = y(cmp);
  double DZ = z(cmp);
  return sqrt(DX*DX + DY*DY + DZ*DZ);
  }

double DetVec::theta(int cmp) const
  { coord pt(x(cmp), y(cmp), z(cmp)); return pt.theta(); }

double DetVec::phi(int cmp) const
  { coord pt(x(cmp), y(cmp), z(cmp)); return pt.phi(); }


// ----------------------------------------------------------------------------
//                    Detection (Sub)Vector Norms
// ----------------------------------------------------------------------------

std::vector<double> DetVec::Norms() const
  {
  std::vector<double> Nrms; 
  int nm = NComps();
  for(int i=0; i<nm; i++)
    Nrms.push_back(Norm(i));
  return Nrms;
  }

void DetVec::Norms(const std::vector<double>& Ns)
  {
  if(!CheckNorms(Ns))  BDVfatal(35);		// Check norm count matches
  int nm = NComps();
  for(int i=0; i<nm; i++)
    Norm(Ns[i], i);
  }

double DetVec::Norm(int cmp) const
  {
  if(!CheckRange(cmp)) BDVfatal(104, Gdec(cmp));
  int     I = 3*cmp;				// Index of component
  double Dx = getRe(I);				// Get Dx of component cmp
  double Dy = getRe(I+1);			// Get Dy of component cmp
  double Dz = getRe(I+2);			// Get Dz of component cmp
  return sqrt(Dx*Dx + Dy*Dy + Dz*Dz);		// Return norm of component
  }

void DetVec::Norm(double nv, int cmp) 
  {
  if(!CheckRange(cmp)) BDVfatal(105, Gdec(cmp));
  double sf =  nv/Norm(cmp);			// Scaling to adjust norm
  int I = 3*cmp;				// Index of component
  put(get(I)*sf, I);				// Set Dx of component cmp
  put(get(I+1)*sf, I+1);			// Set Dy of component cmp
  put(get(I+2)*sf, I+2);			// Set Dz of component cmp
  }

// ____________________________________________________________________________
// C              Bloch Detection Vector Parameter Set Functions
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//       Functions To Make A Bloch Detection Vector From A Parameter Set
// ----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                   Parameter Set From Bloch Detection Vector
//-----------------------------------------------------------------------------

DetVec::operator ParameterSet( ) const
  { ParameterSet pset; pset += *this; return pset; }

void operator+= (ParameterSet& pset, const DetVec& BDV)
  { BDV.PSetAdd(pset); }

bool DetVec::PSetAdd(ParameterSet& pset, int pfx) const
  {
  std::string prefx;						// Parameter prefix
  if(pfx != -1)						// Use prefix if pfx
    prefx = std::string("[")+Gdec(pfx)+std::string("]");		// is NOT -1
  int NC = NComps();					// Get # subvectors
  std::string pname = prefx+std::string("NDetVecs");		// # of subvectors
  std::string pstate("No. of Detection SubVectors");		// Parameter statement
  SinglePar par(pname, NC, pstate);			// # subvects parameter
  pset.push_back(par);					// Add param to pset
  coord MSV;						// Working coordinate
  for(int i=0; i<NC; i++)				// Loop sub-vectors
    {
    pname = prefx + "MCoord(" + Gdec(i) + ")";		//   Coordinate name
    pstate = std::string("Detection SubVector Comps");	//   Coordinate state
    MSV = coord(getRe(3*i), getRe(3*i+1), getRe(3*i+2));//   Coordiante
    par = MSV.param(pname, pstate);			//   Coordinate param
    pset.push_back(par);				//   Add parameter
    }
  return true;
  }
 
//-----------------------------------------------------------------------------
//            Parameter Set File From Bloch Detection Vector
//-----------------------------------------------------------------------------

        // Input                DV       : Detection vector (this)
        //                      filename: Output file name
        //                      ofstr   : Output file stream
        //                      pfx     : Vector prefix (default -1)
        //                      warn    : Warning level
        // Output               none    : Vector is written as a parameter
        //                                to file or output file stream

bool DetVec::write(const std::string& filename, int pfx, int warn) const
  {
  std::ofstream ofstr(filename.c_str());     // Open filename for input
  if(!write(ofstr, pfx, warn?1:0))      // If file bad then exit
    {
    BDVerror(40, filename, 1);           // Filename problems
    if(warn>1) BDVfatal(20);              // Fatal error
    return false;
    }
  ofstr.close();                        // Close it now
  return true;
  }

bool DetVec::write(std::ofstream& ofstr, int pfx, int warn) const
  {
  ParameterSet pset;                    // Declare a parameter set
  PSetAdd(pset, pfx);                   // Add in interaction parameters
  if(!pset.write(ofstr, warn?1:0))      // Use parameter set to write
    {
    if(warn)
      {
      BDVerror(22, 1);                    // Problems writing to filestream
      if (warn>1) BDVfatal(23);           // Fatal error
      }
    return false;
    }
  return true;
  }

// ____________________________________________________________________________
// D                  Bloch Detection Vector Input Functions
// ____________________________________________________________________________

bool DetVec::read(const std::string& filein, int idx, int warn)
  {
  ParameterSet pset;                    // Declare a parameter set
  if(!pset.read(filein, warn?1:0))      // Read in pset from file
    {
    if(warn)
      {
      BDVerror(40, filein, 1);		// Problems with file filein
      if(warn>1) BDVfatal(101,filein);	// Can't read from file filein
      else       BDVerror(101,1);
      }
    return false;
    }
  return read(pset, idx, warn);
  }

bool DetVec::read(const ParameterSet& pset,int idx, int warn)
  {
  if(SetVector(pset, idx, warn?true:false))	// Try and set the vector
    return true;
  if(warn)                                      // Looks like we can't read
    {                                           // the system at all
    std::string sl("");				//   String for index
    if(idx != -1) sl = std::string(Gdec(idx));	//   Set it if one exists
                  BDVerror(77, sl, 1);
    if(warn > 1)  BDVfatal(5);
    else          BDVerror(78, sl);
    }
  return false;
  }

std::string DetVec::ask_read(int argc, char* argv[], int argn)
  {
  std::string filename;				// Name of input file
  std::string M("\n\tDetection vector filename? ");	// Query we will ask if
  query_parameter(argc, argv, argn,M,filename);	// Get/ask for filename
  read(filename);                               // Read system from filename
  return filename;				// Return filename used
  }

std::string DetVec::ask_read(int argc, char* argv[], int argn, const std::string& def)
  {
  std::string M = "\n\tDetection vector filename"	// Query we will ask if
           + std::string(" [") + def + std::string("]? ");// it is needed
  std::string filename = def;                        // Name of spin system file
  ask_set(argc,argv,argn,M,filename);		// or ask for it
  read(filename);                               // Read system from filename
  return filename;                              // Return filename used
  }
 
// ____________________________________________________________________________
// E                      SYSTEM STANDARD I/O FUNCTIONS
// ____________________________________________________________________________

        // Input                BDV      : A magnetization vector (this)
        //                      ostr    : Output stream
        // Output               none    : Vector info is sent
        //                                into the output stream

std::vector<std::string> DetVec::printStrings() const
  {
  std::vector<std::string> PStrings;		// Strings we will return
  int NS = NComps();				// Get number components
  std::string hdr("Bloch Detection Vector");	// Output header
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

                                 Bloch Detection Vector

                     <Di|     Dx       Dy       Dz    ||Di|| 
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

  std::string MXC("    Dx    ");			// 2nd column header
  std::string MYC("    Dy    ");			// 3rd column header
  std::string MZC("    Dz    ");			// 4th column header
  std::string MNC("  ||Di||  ");			// 6th column header

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

std::ostream& DetVec::print(std::ostream& ostr, int np) const
  {
  int NS = NComps();				// Get number components
  std::string hdr("Bloch Detection Vector");	// Output header
  int hl;					// Header length
  if(!NS)					// If there are no vectors
    {						// just write system is empty
    hdr = "Empty " + hdr;
    hl = hdr.length();
    ostr << "\n" << std::string(40-hl/2, ' ') << hdr;
    return ostr;
    }

/*          Output Some Printed Lines Which Will Look Like The Following

                                 Bloch Detection Vector

                     <Di|     Dx       Dy       Dz    ||Di|| 
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

  std::string MXC("    Dx    ");			// 2nd column header
  std::string MYC("    Dy    ");			// 3rd column header
  std::string MZC("    Dz    ");			// 4th column header
  std::string MNC("  ||Di||  ");			// 6th column header

  int    cgap = 1;                              // Gap between columns
  std::string cspc = std::string(cgap, ' ');              // String for gap

  std::vector<std::string> chdrs;                         // Column headers
  chdrs.push_back(std::string("<Di|"));		// 1st column header
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

std::ostream& operator<< (std::ostream& ostr, const DetVec& BDV)
  { return BDV.print(ostr); }

// ____________________________________________________________________________
// F                  Specialized Bloch Detection Vectors
// ____________________________________________________________________________


DetVec DetVec::Dx() const
  {
  DetVec D(*this);
  for(int i=0; i<NComps(); i++)
    {
    D.put(1.0, 3*i);
    D.put(0.0, 3*i+2);
    D.put(0.0, 3*i+3);
    }
  return D;
  }

DetVec DetVec::Dy() const
  {
  DetVec D(*this);
  for(int i=0; i<NComps(); i++)
    {
    D.put(0.0, 3*i);
    D.put(1.0, 3*i+2);
    D.put(0.0, 3*i+3);
    }
  return D;
  }

DetVec DetVec::Dz() const
  {
  DetVec D(*this);
  for(int i=0; i<NComps(); i++)
    {
    D.put(0.0, 3*i);
    D.put(0.0, 3*i+2);
    D.put(1.0, 3*i+3);
    }
  return D;
  }

DetVec DetVec::DxVec(int NC)
  {
  DetVec D(NC);
  for(int i=0; i<D.NComps(); i++)
    { D.put(1.0, 3*i); D.put(0.0, 3*i+2); D.put(0.0, 3*i+3); }
  return D;
  }

DetVec DetVec::DyVec(int NC)
  {
  DetVec D(NC);
  for(int i=0; i<D.NComps(); i++)
    { D.put(0.0, 3*i); D.put(1.0, 3*i+2); D.put(0.0, 3*i+3); }
  return D;
  }

DetVec DetVec::DzVec(int NC)
  {
  DetVec D(NC);
  for(int i=0; i<D.NComps(); i++)
    { D.put(0.0, 3*i); D.put(0.0, 3*i+2); D.put(1.0, 3*i+3); }
  return D;
  }

DetVec DetVec::DxVec() const { return Dx(); }
DetVec DetVec::DyVec() const { return Dy(); }
DetVec DetVec::DzVec() const { return Dz(); }


// ____________________________________________________________________________
// F                Bloch Detection Vector Evolution Functions
// ____________________________________________________________________________
    
/* These functions are decidedly Bloch based. That is, they are set to evolve 
   the magnetization vector according to the phenomenological Bloch equations.
   The user must supply a 3
*/

#endif						// DetVec.cc
