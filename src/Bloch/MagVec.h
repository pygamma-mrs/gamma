/* MagVec.h *****************************************************-*-c++-*-
**									**
**	                        G A M M A				**
**								 	**
**	Bloch Magnetization Vector 		Interface		**
**						 			**
**      Copyright (c) 2002                                              **
**      Dr. Scott A. Smith						**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**						 			**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**								 	**
** Description							 	**
**								 	**
** The class MagVec defines a Bloch magnetization vector. Such vectors	**
** are meant to evolve in time as described by the phenomenological	**
** Bloch equations. A magnetization vector as specified herein involves **
** N individual sub-vectors, each of these sub-vectors has three	**
** magnetization components. The vector with N components takes the 	**
** blocked form								**
**		     |M> = ||M > |M > |M > ..... |M   >			**
**                            0    1    2          N-1			**
**								 	**
** & each sub-vector is							**
**                        |M > = |M  , M   , M  >		 	**
**                          i      ix   iy    iz			**
**								 	**
** The components of |M> are then given by				**
**								 	**
**           <3*i|M> = M     <3*i+1|M> = M     <3*i+2|M> = M		**
**                      ix                iy                iz		**
**								 	**
** a static Bo field, an applied rf-field (B1), simplistic relaxation	**
** (T1 & T2), and possibly exchange between magnetizaiton vectors.	**
**									**
** The magnetization vector has ONLY 1 type of evolution implicitly	**
** defined, that dictated by the phenomenological Bloch equations.	**
**									**
**			     -Gt					**
**                 |M(t)> = e	|M(0)-M	  > + |M   >			**
**                                     inf      inf			**	
**									**
** where the evolution matrix G and the initinite time vector M		**
** must be specified appropriately.                            inf	**
**									**
*************************************************************************/

#ifndef   MagVec_h_			// Is file already included?
#  define MagVec_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <string>			// Include libstdc++ strings
#include <vector>			// Include libstdc++ STL vectors
#include <Matrix/col_vector.h>		// Include GAMMA column vectors
#include <Basics/ParamSet.h>		// Include GAMMA parameter sets
#include <Level1/coord.h>		// Include GAMMA coordinates

class MagVec : public col_vector
{

friend class BlochSys;
  
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                    Magnetization Vector Error Handling
// ____________________________________________________________________________


        // Input                sys     : Spin system (this)
        //                      ei	: Error index
        //                      nr	: Flag for linefeed (0=linefeed)
        //                      pn	: string in message
     
         void MVerror(int ei,                        int nr=0) const;
         void MVerror(int ei, const std::string& pn, int nr=0) const;
volatile void MVfatal(int ei) const;
volatile void MVfatal(int ei, const std::string& pn) const;

// ____________________________________________________________________________
// ii                   Magnetization Vector Checking Functions
// ____________________________________________________________________________

bool CheckNorms(const  std::vector<double>&  Ns, bool warn=true) const;
bool CheckRange(int cmp,                    bool warn=true) const;

// ____________________________________________________________________________
// iii                     Magnetization Vector Setup Functions
// ____________________________________________________________________________

/* These are protected functions because they allow specific aspects of the
   vector to be set up without worrying about consistency! Basically what we
   must determine are

                          1.) Number of sub-vectors
                          2.) Mx, My, and Mz for each sub-vector

   Parameter NMagVecs is used indicate how many sub-vectors are defined.
   If that is NOT found, it is assumed that there is only 1 sub-vector.     */

bool SetVector(const   ParameterSet& P, int pfx=-1, bool W=true);
bool GetNVects(const   ParameterSet& P, int&     N, bool W=true) const;
bool SetSubVects(const ParameterSet& P, int      N, bool W=true);
bool GetCoord(const  ParameterSet& P, 
                                     coord& pt, int idx=-1, bool W=true) const;
bool GetMxMyMz(const   ParameterSet& P, 
            double& Mx, double& My, double& Mz, int idx=-1, bool W=true) const;

// ----------------------------------------------------------------------------
//         These Read All Single Magnetization Vector Parameters
// ----------------------------------------------------------------------------


//bool GetVect(const ParameterSet& pset, int i, double& v, Isotope& I,
//             double& R1, double& R2, coord& Pt, int& Sp, bool warn=true) const;

// ----------------------------------------------------------------------------
//         These Read All Magnetization Vector Pair Parameters
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

// ____________________________________________________________________________
// A                   SPIN SYSTEM CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

///Center Spin System Algebraic
///F_list MagVec		- Constructor
///F_list ~			- Destructor
///F_list =			- Assignment	

// ----------------------------------------------------------------------------
//                          Simple Constructors
// ----------------------------------------------------------------------------

MSVCDLC MagVec(int   N=0);
MSVCDLC MagVec(const MagVec& MV);
MSVCDLC MagVec(const col_vector& CV);

// ----------------------------------------------------------------------------
//              Constructors Using Individual Magnetization Vectors
// ----------------------------------------------------------------------------

MSVCDLC MagVec(double Mx, double My, double Mz);
MSVCDLC MagVec(const coord& M);
MSVCDLC MagVec(double Mx1,double My1,double Mz1,double Mx2,double My2,double Mz2);
MSVCDLC MagVec(const coord& M1, const coord& M2);
MSVCDLC MagVec(const std::vector<coord>& Ms);

// ----------------------------------------------------------------------------
//                       Assignment and Destruction
// ----------------------------------------------------------------------------

//        ~MagVec();
//MagVec& operator= (const MagVec& MV);

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

//MagVec operator-  ()                 const;
MSVCDLL MagVec  operator+  (const MagVec& M1) const;
MSVCDLL MagVec& operator+= (const MagVec& M1);
MSVCDLL MagVec  operator-  (const MagVec& M1) const;
MSVCDLL MagVec& operator-= (const MagVec& M1);

// ____________________________________________________________________________
// C                    Magnetization Vector Access
// ____________________________________________________________________________

//int size()   const;			// Inherited
MSVCDLL int NComps() const;

// ----------------------------------------------------------------------------
//                  Magnetization (Sub)Vector Components
// ----------------------------------------------------------------------------

MSVCDLL double Mx(int cmp) const;
MSVCDLL double My(int cmp) const;
MSVCDLL double Mz(int cmp) const;

MSVCDLL void   Mx(int cmp, double mx);
MSVCDLL void   My(int cmp, double my);
MSVCDLL void   Mz(int cmp, double mz);

MSVCDLL double x(int     cmp=0) const;
MSVCDLL double y(int     cmp=0) const;
MSVCDLL double z(int     cmp=0) const;
MSVCDLL double norm(int  cmp=0) const;
MSVCDLL double theta(int cmp=0) const;
MSVCDLL double phi(int   cmp=0) const;

// ----------------------------------------------------------------------------
//                     Magnetization (Sub)Vector Norms
// ----------------------------------------------------------------------------

MSVCDLL std::vector<double> Norms() const;
MSVCDLL void                Norms(const std::vector<double>& Ns);

MSVCDLL double Norm(int i) const;
MSVCDLL void   Norm(double nv, int i);

// ____________________________________________________________________________
// C              Magnetization Vector Parameter Set Functions
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//       Functions To Make A Magnetization Vector From A Parameter Set
// ----------------------------------------------------------------------------

/*
virtual int write(const std::string &filename, int idx=-1, int warn=2) const;
virtual int write(std::ofstream& ofstr, int idx=-1, int warn=2) const; 
*/

//-----------------------------------------------------------------------------
//                   Parameter Set From Magnetization Vector
//-----------------------------------------------------------------------------

MSVCDLL             operator ParameterSet( ) const;
MSVCDLL friend void operator+= (ParameterSet& pset, const MagVec& MV);
MSVCDLL bool PSetAdd(ParameterSet& pset, int pfx=-1)   const;

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

MSVCDLL bool write(const std::string& filename, int pfx=-1, int warn=2) const;
MSVCDLL bool write(std::ofstream& ofstr,        int pfx=-1, int warn=2) const;

// ____________________________________________________________________________
// D                  Magnetization Vector Input Functions
// ____________________________________________________________________________

///F_list read		- Read magnetization vector from disk file
///F_list ask_read	- Ask for file, read magnetization vector from file

	// Input		M	 : Magnetization vector (this)
	// 			filename : Input filename
	//			idx	 : Parameter index value used for
	//				   prefix [#] in input names
	//			argc	 : Number of arguments
	//			argv     : Vector of argc arguments
	//			argn     : Argument index
	//			warn	 : Warning output level
	//					0 = no warnings
	//					1 = warnings
	//				       >1 = fatal warnings
	// Output		none	 : Vector is filled with values in
	//				   parameters read from file
	// Output		string   : The parameter argn of array argc
	//				   is used to supply a filename
	//				   from which the spin system is read
	//				   If the argument argn is not in argv,
	//				   the user is asked to supply a filename
        //                                 The set filename is returned 
	// Note			 	 : The file should be an ASCII file
	//				   containing recognized M parameters
	// Note			 	 : The vector is modifed (filled)


MSVCDLL bool read(const std::string& fn,    int idx=-1, int warn=2);
MSVCDLL bool read(const ParameterSet& pset, int idx=-1, int warn=2);

MSVCDLL std::string ask_read(int argc, char* argv[], int argn);
MSVCDLL std::string ask_read(int argc, char* argv[], int argn,
                                                       const std::string& def); 

// ____________________________________________________________________________
// E                Magnetization Vector Standard I/O Functions
// ____________________________________________________________________________

MSVCDLL std::vector<std::string> printStrings() const;
MSVCDLL std::ostream&            print(std::ostream& out, int np=20) const;
MSVCDLL friend std::ostream& operator<<(std::ostream& out, const MagVec& M);

// ____________________________________________________________________________
// F                  Specialized Magnetization Vectors
// ____________________________________________________________________________

MSVCDLL        MagVec Mx()    const;
MSVCDLL        MagVec My()    const;
MSVCDLL        MagVec Mz()    const;
MSVCDLL static MagVec MxVec(int NC);
MSVCDLL static MagVec MyVec(int NC);
MSVCDLL static MagVec MzVec(int NC);
MSVCDLL        MagVec MxVec() const;
MSVCDLL        MagVec MyVec() const;
MSVCDLL        MagVec MzVec() const;

// ____________________________________________________________________________
// G                Magnetization Vector Evolution Functions
// ____________________________________________________________________________


};

#endif						// MagVec.h
