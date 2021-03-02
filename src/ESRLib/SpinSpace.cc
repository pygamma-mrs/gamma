/* SpinSpace.cc *************************************************-*-c++-*-
**									**
** 	                       G A M M A				**
**									**
**	Spin Hilbert Space                           Implementation	**
**					                                **
**	Copyright (c) 2000				                **
**	Scott A. Smith                                                  **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**							 		**
**      $Header: 
**								 	**
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** Class SpinSpace defines a quantity which contains a count of the	    **
** number of spins and their isotope types.  As it is assumed that	    **
** complex spin system classes will be derived from SpinSpace, the	    **
** following functions have been declared virtual.                      **
**                                                                      **
**	~	The SpinSpace destructor                                         **
**	print	The SpinSpace printing function                              **
**                                                                      **
** This allows any any derived classes to overwrite the	defaults of 	**
** these functions.	                                                    **
**                                                                      **
*************************************************************************/

#ifndef   SpinSpace_cc_				// Is the file already included?
#  define SpinSpace_cc_				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation			// Then this is the implementation
#  endif

#include <HSLib/SpinSpace.h>			// Include the interface
#include <Basics/Gutils.h>			// Include GAMMA errors
#include <Basics/StringCut.h>			// Include GAMMA Gdec function

int    SpinSpace::_warn  = 2;			// Set warnings active by default
string SpinSpace::DefIso = string("1H");	// Set default isotope type to 1H


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                     CLASS SPIN SPACE ERROR HANDLING
// ____________________________________________________________________________

/*      Input               HSSp    : Hilbert spin space (this)
                            eidx    : Error index
                            noret   : Flag for linefeed (0=linefeed)
                            pname   : string in message
        Output              none    : Error message output
                                      Program execution stopped if fatal     */
                                      
void SpinSpace::error(int eidx, int noret) const
  {
  string hdr("Spin Hilbert Space");
  switch (eidx)
    {
    case 0: GAMMAerror(hdr, 0, noret); break;	// Program Aborting        (0)
    case 3: GAMMAerror(hdr, 3, noret); break;	// !Construct From Pset    (3)
    case 4: GAMMAerror(hdr, 4, noret); break;	// !Construct From File    (4)
    case 5: GAMMAerror(hdr, 5, noret); break;	// !Write To Param File    (5)
    case 6: GAMMAerror(hdr, 6, noret); break;	// !Write To Filestream    (6)
    case 7: GAMMAerror(hdr, 7, noret); break;	// !Output Anything	   (7)
    case 10: GAMMAerror(hdr,string("Accessed Spin Index Out of Range"),noret);
             break;                             //                         (10)
    case 11: GAMMAerror(hdr,string("Accessed Spin Pair Out of Range"),noret);
             break;                             //                         (11)
    case 12: GAMMAerror(hdr,string("Requested Isotope NULL System"), noret);
             break;                             //                         (12)
    case 13: GAMMAerror(hdr,string("Cannot Access Spin Isotope Type"), noret);
             break;                             //                         (13)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  if(!noret) cout << ".\n";
  }

void SpinSpace::error(int eidx, const string& pname, int noret) const
  {
  string hdr("Spin Hilbert Space");
  string msg;
  switch(eidx)
    {
    case 1:  GAMMAerror(hdr, 1, pname, noret); break;	// File Problems   (1)
    case 2:  GAMMAerror(hdr, 2, pname, noret); break;	// !Read SingPar   (2)
    case 3:  msg = string("Warning - Unspecified Isotope Type(s) Set To ");
             GAMMAerror(hdr,msg+pname,noret);  break;	//		   (3)
    case 4:  msg = string("Warning - No Isotope Type Information: ");
             GAMMAerror(hdr,msg+pname,noret);  break;	//		   (4)
    default: GAMMAerror(hdr, -11, pname, noret); break;	// Unknown Error   (-1)
    }
  }
 
volatile void SpinSpace::fatality(int eidx) const
  {  
  error(eidx, 1);				// Normal non-fatal error
  if(eidx) error(0);				// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                   SPIN SPACE CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________
 

SpinSpace::SpinSpace(int ns)
  {
  nspins   = ns;				// Set the number of spins
  sysname  = string("");			// There is no system name
  if(ns>0)					// If there are spins
    {
    Isotope IsoTmp(DefIso);			//   The default isotope type
    for(int i=0; i<ns; i++)			//   Now we initialize all of
      {						//   the isotopes & spin flags
      spinflags.push_back(1);			//	Set spin flag to TRUE
      Isotopes.push_back(IsoTmp);		//      Set spin type to DefIso
      }
    int hs = HS();				// Get system Hilbert space
    bsmx = matrix(hs, hs, i_matrix_type);	// Set a default basis (matrix)
    }
  }


SpinSpace::SpinSpace(const SpinSpace &sys)
  {
  nspins    = sys.nspins;			// Copy number of spins
  sysname   = sys.sysname;			// Copy the system name
  spinflags = sys.spinflags;			// Copy the spin flags
  Isotopes  = sys.Isotopes;			// Copy the spin isotopes
  bsmx      = sys.bsmx;				// Copy the basis matrix
  }


SpinSpace::~SpinSpace() { }


SpinSpace& SpinSpace::operator=(const SpinSpace& sys)
  {
  if(this == &sys) return *this;		// Avoid self copying
  nspins    = sys.nspins;			// Copy the number of spins
  spinflags = sys.spinflags;			// Copy the spin flags
  Isotopes  = sys.Isotopes;			// Copy the spin isotopes
  sysname   = sys.sysname;			// Copy the system name
  bsmx      = sys.bsmx;				// Copy the basis matrix
  return *this;
  }

// ____________________________________________________________________________
// B                              COMPARISONS
// ____________________________________________________________________________
 
/* Systems are considered equal if # spins and spin isotopes types match.
   This is irrespective of how system spin flags are set.                    */

int SpinSpace::operator==(const SpinSpace& sys) const
  {
  if(this == &sys) return 1;			// True if sys is this!
  if(nspins != sys.spins()) return 0;		// Fail if # spins mismatch
  for(int i=0; i<nspins; i++)			// Insure all spin types equal
    if(Isotopes[i] != sys.Isotopes[i]) return 0;
  return 1;
  }
 
int SpinSpace::operator!=(const SpinSpace &sys) const
  {
  if(nspins != sys.spins()) return 1;		// True if # spins mismatch
  for(int i=0; i<nspins; i++)			// Insure not all spin types
    if (Isotopes[i] != sys.Isotopes[i]) 	// types are equal
      return 1;
  return 0;
  }
 
// ____________________________________________________________________________
// C                   BASIC SPIN SPACE MANIPULATIONS
// ____________________________________________________________________________
 
int SpinSpace::spins() const { return HSs.size(); }
 
int SpinSpace::HS() const
  {
  int hs=1;
  for(int i=0; i<nspins; i++)
    hs *= (Isotopes[i]).HS();
  return hs;
  }

int SpinSpace::HS(int spin) const
  { check_spin(spin); return Isotopes[spin].HS(); }

// ____________________________________________________________________________
// D             SPIN ANGULAR MOMENTUM AND ISOTOPE MANIUPLATIONS
// ____________________________________________________________________________

const Isotope& SpinSpace::isotope(int spin) const
  { check_spin(spin); return Isotopes[spin]; }

double SpinSpace::weight(int spin) const
  { check_spin(spin); return Isotopes[spin].weight(); }

string SpinSpace::symbol(int spin) const
  { check_spin(spin); return Isotopes[spin].symbol(); }

double SpinSpace::qn(int spin) const
  { check_spin(spin); return Isotopes[spin].qn(); }

double SpinSpace::qn() const
  { 
  double q = 0;
  for(int spin=0; spin<nspins; spin++)
    q += Isotopes[spin].qn();
   return q;
  }

string SpinSpace::element(int spin) const
  { check_spin(spin); return Isotopes[spin].element(); }

double SpinSpace::gamma(int spin) const
  { check_spin(spin); return Isotopes[spin].gamma(); }

double SpinSpace::gamma(const string& iso) const
  { return (Isotope(iso)).gamma(); }
 
void SpinSpace::isotope(int spin, const string& symbol)

        // Input                int     : Spin label [0, nspins-1]
        // 			symbol	: Isotope label (e.g. 1H, 13C)
	// Output		void    : System spin is changed to symbol

  { isotope(spin,Isotope(symbol)); }

 
void SpinSpace::isotope(int spin, const Isotope& Iso)

        // Input                int     : Spin label [0, nspins-1]
        // 			Iso     : Isotope type (e.g. 1H, 13C)
	// Output		void    : System spin is changed to Iso

  {
  if(!check_spin(spin,1)) fatality(13); 	// Insure spin exists
  double Ival = Isotopes[spin].qn();		// I value of current type
  Isotopes[spin] = Iso;				// Set new isotope type
  if(Ival != Iso.qn())				// If Hilbert space change
    bsmx = matrix(HS(), HS(), i_matrix_type);	// construct new basis (matrix)
  }


string SpinSpace::momentum (int spin) const

        // Input                int     : spin label [0, nspins-1]
        // Output               char *  : momentum (e.g. 1/2, 3/2)
        // Note 			  default is 1/2 (proton)

  { check_spin(spin); return Isotopes[spin].momentum(); }


string SpinSpace::momentum() const
  {
  double d=qn();
  if (int(d)==d) return Gdec(int(d));
  else           return string(Gdec(int(2*d)))+string("/2");
  }

 
row_vector SpinSpace::qState(int state) const
 
        // Input                sys	: Spin system (this)
        // 			int	: Number of state (0<=state<HS())
        // Output               vx	: Row vector whose elements contain
	//				  the quantum number of each spin

  {
  row_vector states(nspins);			// Vector, element per spin
  double mz;					// Spin mz value
  int hs; 					// For spin Hilbert space
  for(int spin=nspins-1; spin>=0; spin--)	// Loop over spins
    {
    mz = Isotopes[spin].qn();			//   Spin mz value
    hs = Isotopes[spin].HS();			//   Spin hilbert space
    states.put(complex(mz- state%hs), spin);	//   Set element for spin
    state /= hs;
    }
  return states;
  }


matrix SpinSpace::qStates() const
 
        // Input                sys	: Spin system (this)
        // Output               matrix  : Matrix whose rows contain the
	//				  information of base function qState

  {
  matrix states(HS(), nspins);		// States, rows=bfs, cols=spins
  int hs, tHS = HS();			// Hilbert space dimension
  int bs = tHS;				// Blocks with qn have size bs
  int i, state;
  for(int spin=0; spin<nspins; spin++)	// Loop over the spins
    {
    hs = HS(spin);			//   This spins Hilbert space
    bs /= hs;				//   Blocks with same quantum #
    double qns = qn(spin);		//   HS = 2*qn+1
    int pos=0; 
    while(pos<tHS)
      for(state=0; state<hs; state++)
        for(i=0; i<bs; i++, pos++)
          states.put(complex(qns-state),pos,spin);
    }
  return states;
  }


double SpinSpace::qnState(int state) const

        // Input                sys	: Spin system (this)
        // 			int	: Number of state (0<=state<HS())
        // Output               double  : Quantum number of a state

  {
  double qn = 0;
  for(int spin=nspins-1; spin>=0; spin--)
    {
    qn += Isotopes[spin].qn()- (state%Isotopes[spin].HS());
    state /= Isotopes[spin].HS();
    }
  return qn;
  }


col_vector SpinSpace::qnStates() const
 
        // Input                sys	: Spin system (this)
        // 			int	: Number of state (0<=state<HS())
        // Output               vx	: Column vector whose HS() elements
	//				  contain quantum #s of all states

  {
  col_vector qnv(HS());
  qnv.put(complex(qn()),0);
  int bs=1;
  for(int spin=nspins-1; spin>=0; spin--)
    {
    int hs=HS(spin);
    int pos=bs;
    for(int q=1; q<hs; q++)
      for(int p=0; p<bs; p++, pos++)
	qnv.put(qnv.get(pos-bs)-1, pos);
    bs *= hs;
    }
  return qnv;
  }


row_vector SpinSpace::qnDist() const
  {
  row_vector dist(int(2*qn()+1), complex0);	// Vector for output
  int size=1;					// Current # used slots in dist 
  int length=1;
  double tmp;
  dist.put(complex1,0);				// Start @ 1 for this coherence
  for(int spin=0; spin<nspins; spin++)
    {
    length = int(2*qn(spin))+1;
    size += length-1;
    for (int i=size-1; i>=0; i--)
      {
      tmp = 0;
      for(int j=0; j<length; j++)
      if((i-j) >= 0)
        tmp += Re(dist.get(i-j));
      dist.put(complex(tmp),i);
      }
    }
  return dist;
  }


row_vector SpinSpace::CoherDist() const
  {
  row_vector dist(int(4*qn()+1), complex0);	// Vector for output
  row_vector qndist=qnDist();			// Get distribution of levels
  int qnt=int(2*qn());
  for (int i=0; i<=qnt; i++)
    for (int j=0; j<=qnt; j++)
      dist.put( dist.get(qnt+j-i) + qndist.get(i)*qndist.get(j), (qnt+j-i) );
  return dist;
  }

 
/*    Function   Arguments Output                     Comments
   ------------- --------- ------  --------------------------------------------
    homonuclear    ---      bool   True if all spins in system homonuclear
   heteronuclear   ---      bool   True if any spin types s in system differ
     spinhalf      ---      bool   True if all spin in system have I=1/2
     electrons     ---      int    Returns the number of electron in system
      nepair       i,j      bool   True if spin i & j are a nucleus/e- pair
      enpair       i,j      bool   True if spin i & j are a nucleus/e- pair
     pairidx       i,j      int    Index (dipolar) for spin pair i&j
     isotopes      ---      int    Number of unique isotopes in system
     isotopes       i      string  Isotope type of (isotope) index i 
     isotopes     string    bool   True if isotope of type input in system   */

 
bool SpinSpace::homonuclear() const
  {
  Isotope iso = isotope(0);
  for(int i=1; i<nspins; i++)
    if(isotope(i) != iso) return false;
  return true;
  }

bool SpinSpace::heteronuclear() const
  {
  Isotope iso = isotope(0);
  for(int i=1; i<nspins; i++)
    if (isotope(i) != iso) return true;
  return false;
  }

bool SpinSpace::spinhalf() const
  {
  for(int i=0; i<nspins; i++)
    if(qn(i) != 0.5) return false;
  return false;
  }

int SpinSpace::electrons() const
  {
  int ne = 0;
  for(int i=0; i<nspins; i++)
    if(name(i) == string("Electron")) ne++;
  return ne;
  }

bool SpinSpace::nepair(int i, int j) const { return enpair(i,j); }
bool SpinSpace::enpair(int i, int j) const
  {
  string el("e-");
  if(symbol(i) != el)
     { if(symbol(j) == el) return true; }
  else if(symbol(j) != el) return true;
  return false;
  }

int SpinSpace::pairidx(int i, int j) const
  {
  if(i==j) return -1;			// No i,i pair index
  if(i>nspins-1 || i<0) return -1;	// No index, outside range
  if(j>nspins-1 || j<0) return -1;	// No index, outside range
  int k=0;
  if(i>j)				// Always use j<i indices
    {
    k = j;
    j = i;
    i = k;
    k = 0;
    }
  for(int r=0; r<i; r++)
    k += nspins-1-r;
  k += j - i - 1; 
  return k;
  }
 

int SpinSpace::isotopes() const
  {
  if(!nspins) return 0;			// No isotopes if no spins
  Isotope *Is = NULL;			// Array of isotopes
  Is = new Isotope[nspins];		// Allocate isotopes array
  int niso=0; 				// Counter for isotopes
  Is[niso++] = isotope(0);		// Set 1st isotope & count
  bool fnd;				// Flag if isotope found
  int i=0, j=0;				// Spin/isotope indices
  for(i=1; i<nspins; i++)		// Loop spins in system
    {
    fnd = false;			//   Assume spin is unique type
    for(j=0; j<niso && !fnd; j++)	//   Now search previous types
      if(isotope(i) == Is[j]) fnd=true; //   and see if it is unique
    if(!fnd) Is[niso++] = isotope(i);   //   If it is unique, update
    }					//   the total isotope count
  delete [] Is;				// Don't need this anymore
  return niso;				// Return number of isotopes
  }


string SpinSpace::isotopes(int idx) const
  {
  if(!nspins)				// Can't do anything if no spins
    { error(10,1); fatality(12); }	// Out of range, no spins present
  if(!idx) return symbol(0);		// Done if 1st type requested
  int niso = 0; 			// Number of found isotopes
  int i=0, j=0;				// Working spin/isotope index
  Isotope* Isos = NULL;			// Array of found isotopes
  Isos = new Isotope[nspins];		// Array of isotopes
  Isos[niso++] = isotope(i);		// Store 1st spins type
  bool old;				// Flag if isotope known
  for(i=1; i<nspins; i++)		// Loop the remaining spins
    {
    old = false;			// Assume spin i type new
    for(j=0; j<niso && !old; j++)	// Compare with known types
      {					// and see if it matches.
      if(isotope(i) == Isos[j])		// If so, its just of a type
      old = true; 			// we already know about
      }
    if(!old)				// If it is a new type, see
      {					// if its the one requested
      if(niso == idx)			// then we'll return it after
         {				// we deallocate the array
         delete [] Isos;		// of isotopes we've used
         return symbol(i);
         }
      Isos[niso++] = isotope(i);        // If not one requested, store
      }					// and increment the count
    }
  error(10);				// Spin index is out of range
  fatality(12);				// No spins present
  return symbol(0);			// This shouldn't be reached
  }


bool SpinSpace::isotopes(const string& I) const
  {
  for(int i=0; i<nspins; i++)           // Loop all spins
    if(symbol(i) == I) return true;     // Return true if any of type I
  return false;
  }

// ____________________________________________________________________________
// E                          SPIN FLAG FUNCTIONS
// ____________________________________________________________________________

/* For each spin in the system there exists an On-Off "flag" which can be
   used in functions that act on a set of chosen spins.  They do not have
   any other quality that affects the internal workings of the class. 

   Function     Argument(s)   Output                 Comments
  ---------  ---------------  ------  -----------------------------------------
  flags             TF         void   Sets all spin flags to TF
  flag             i,TF        void   Sets spin i's flag to TF
  flag            iso,TF       void   Sets flags for spins of type iso to TF
  flag              i           int   Gets flag for spin i
  get_flags      int* FLGS     void   FLGS[i] set to Spin i flag
  set_flags    int* FLGS,TF    void   FLGS[i] set to TF for i=[0,nspins)
  set_flags   int* FLGS,i,TF   void   FLGS[i] set to TF if spin i exists
  set_flag   int* FLGS,iso,TF  void   FLGS[i] set to TF if spin i type iso
  set_flag   int* FLGS,iso,TF  void   FLGS[i] set to TF if spin i type iso   */

void SpinSpace::SetFlags(bool TF)
  { for (int i=0; i<nspins; i++) spinflags[i]=TF; }

void SpinSpace::SetFlag(int spin, bool TF)
  { check_spin(spin); spinflags[spin] = TF; }

void SpinSpace::SetFlags(const string& isoin, bool TF)
  { Isotope Iso(isoin); SetFlags(Iso,TF); }

void SpinSpace::SetFlags(const Isotope& Iso, bool TF)
  { for(int i=0; i<nspins; i++) if(Isotopes[i]==Iso) spinflags[i]=TF; }

flagvec SpinSpace::GetFlags() const { return spinflags; }

bool SpinSpace::GetFlag(int i) const
{ if(check_spin(i,2)) return spinflags[i]?true:false; return false; }

flagvec SpinSpace::GetFlags(bool TF) const
  { 
  flagvec FlagsCopy = spinflags;
  for(int i=0; i<nspins; i++) FlagsCopy[i] = TF;
  return FlagsCopy;
  }

flagvec SpinSpace::GetFlags(int spin, bool TF, bool DefTF) const
  { 
  flagvec FlagsCopy=spinflags;
  for(int i=0; i<nspins; i++) FlagsCopy[i] = DefTF;
  if(check_spin(spin,2)) FlagsCopy[spin]=TF;
  return FlagsCopy;
  }

flagvec SpinSpace::GetFlags(const string& isoin, bool TF, bool DefTF) const
  {
  flagvec FlagsCopy=spinflags;
  for(int i=0; i<nspins; i++)
    {
    FlagsCopy[i] = DefTF;
    if(symbol(i)==isoin) FlagsCopy[i]=TF;
    }
  return FlagsCopy;
  }

flagvec SpinSpace::GetFlags(const Isotope& isoin, bool TF, bool DefTF) const
  {
  flagvec FlagsCopy=spinflags;
  for(int i=0; i<nspins; i++)
    {
    FlagsCopy[i] = DefTF;
    if(Isotopes[i]==isoin) FlagsCopy[i]=TF;
    }
  return FlagsCopy;
  }


// ______________________________________________________________________
// F                        SPIN SPACE NAME
// ______________________________________________________________________
 
inline void SpinSpace::name(const string& name) { sysname = name; }

inline const string& SpinSpace::name(int i) const
  { if(i<0) return sysname; check_spin(i); return Isotopes[i].name(); }

void SpinSpace::warnings(int warnf) { _warn = warnf; }

	// Input		sys	: A spin system (this)
	//			warnf	: Warning flag
	// Output		none	: The funciton sets the intermal
	//				  value of _warn to warnf
	// Note				: _warnf determines wheteher some 
	//				  non-fatal problems issue a warning
	//				  message or not


int SpinSpace::warnings() const { return _warn; }

	// Input		sys	: A spin system (this)
	// Output 		warnf	: Value of warning flag


string SpinSpace::IsoDefault() { return DefIso; }
 
        // Input                sys     : A spin system (this)
        // Output               none    : The function returns a string
        //                                of the default isotope type
 

void SpinSpace::IsoDefault(const string& DI)

        // Input                sys     : A spin system (this)
        //                      DI      : string for an isotope type
        // Output               none    : The function sets the intermal
        //                                default isotope type to DI

  { DefIso = DI; }
// sosi - must check that DI is a valid isotope here


// ____________________________________________________________________________
// G                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//          Functions To Make A Parameter Set From A Spin System
// ----------------------------------------------------------------------------

SpinSpace::operator ParameterSet( ) const

	// Input		ss    : A spin system (this)
	// Output		pset  : Parameter set with
	//			        only spin system parameters

  { 
  ParameterSet pset;			// A Null parameter set
  pset += *this; 			// Add in spin system parameters
  return pset;
  }


void operator+= (ParameterSet& pset, const SpinSpace& ss) {ss.PSetAdd(pset);}

	// Input		ss    : A spin system
	//  			pset  : A parameter set
	// Output		pset  : Parameter set with spin system
	//			        parameters added to it


void SpinSpace::PSetAdd(ParameterSet& pset, int idx) const

	// Input		ss	: A spin system
        //                      pset	: Parameter set
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        // Output               void    : Spin system parameters are
        //                                are added to the parameter set
        //                                with interaction prefix idx
        // Note                         : The parameter names & types
        //                                output here MUST match those used
        //                                in setting spin systems
        //                                from parameters sets

  {
  string pname;
  string pdata;
  string pstate;
  SinglePar par;
  string prefx;                                 // Parameter prefix
  if(idx != -1)                                 // Only use suffix if idx
    prefx = string("[")+Gdec(idx)+string("]");   // is NOT -1

   pname = prefx + string("SysName");		// Add Spin System Name
   pstate = string("Name of the Spin System");
   pdata = name();
   par = SinglePar(pname, 2, pdata, pstate);
   pset.push_back(par);

   pname = prefx + string("NSpins");		// Add Number of Spins
   pstate = string("Number of Spins in the System");
   pdata = Gdec(spins());
   par = SinglePar(pname, 1, pdata, pstate);
   pset.push_back(par);

   int ns = spins();
   pstate = string("Spin Isotope Type");	// Add Isotope Labels
   for(int i=0; i < ns; i++)
     {
     pname = prefx + string("Iso(");
     pname += Gdec(i);
     pname += string(")");
     pdata = symbol(i);
     par = SinglePar(pname, 2, pdata, pstate);
     pset.push_back(par);
     }
   return;
   } 


// ----------------------------------------------------------------------------
//           Functions To Make A Spin System From A Parameter Set
// ----------------------------------------------------------------------------

 
int SpinSpace::getSpins(const ParameterSet& pset, int warn) const
 
        // Input                sys	: A spin system (this)
        //                      pset    : A parameter set
	//			warn	: Warning level
	//					0 = no warnings
	//					1 = warnings
	//				       >1 = fatal warnings
        // Output               ns      : Number of spins specified by 
        //                                parameter NSpins in pset
        // Note                         : This does NOT set the number of spins
	// Note				: Return -1 if # of spins not found
 
  {
  int ns = 0;					// Assume no spins
  string pstate, pname("NSpins");		// Parameter name
  list<SinglePar>::const_iterator item;		// A pix into parameter list
  item = pset.seek(pname);		 	// Pix in param. list for NSpins
  if(item != pset.end()) 			// Retrieve the number of spins
    (*item).parse(pname,ns,pstate);		// if it was found in the list
  else
    {
    if(warn)
      {
      error(2, pname, 1); 			// Can't read in NSpins
      if(warn > 1)  fatality(3);		// This is a fatal error
      }
    ns = -1;
    }
  return ns;
  }  

 
void SpinSpace::setName(const ParameterSet& pset)
 
        // Input                sys	: A spin system (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system name is
        //                                set from parameter in pset
        // Note                         : This is not required
 
  {
  string pname("SysName");                     // System parameter name
  string sval, pstate;                          // System parameter value
  list<SinglePar>::const_iterator item;		// A pix into parameter list
  item = pset.seek(pname);		 	// Pix in pset for SysName
  if(item != pset.end()) 			// Retrieve the system name
    {
    (*item).parse(pname,sval,pstate);
    name(sval);
    }
  }


void SpinSpace::setIs(const ParameterSet& pset)

        // Input                sys 	: A spin system (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system isotope types
        //                                are set from parameters in pset
	// Note				: Any isotopes not specified 
	//				  will be set to 1H

  {
  string sval, pstate;   	                // System parameter value
  string* NoISet;	   			// string of I's not found
  NoISet = new string[nspins];
  string Ipname, Ist("Iso(");			// Start of I parameter names
  string Ifi(")");				// Finish of I parameter names
  list<SinglePar>::const_iterator item;		// A pix into parameter list
  int i, count=0;				// Counter of isotopes
  for(i=0; i<nspins; i++)			// Retrieve the isotope types
    {                                           // Isotopes not specified are
    Ipname = Ist + Gdec(i) + Ifi;                // Set parameter name
    item = pset.seek(Ipname);		 	// Pix in pset for Ipname
    if(item != pset.end()) 			// Retrieve the parameter
      {
      (*item).parse(Ipname,sval,pstate);
      isotope(i, sval);
      } 
    else
      {
      isotope(i, DefIso);			// Set to default isotope type
      NoISet[count] = Ipname;			// Store par name not found
      count++;					// Count isotopes not specified
      } 
    }   
  int j=0;					// Counter for linefeed
  if(count && _warn)				// Output warnings if desired
    {						// that isotope type(s) haven't 
    for(i=0; i<count; i++)			// been specified in pset
      {
      j++;
      if(j==1) error(4, NoISet[i], 1);	         // Warning that no I set
      else
        {
        cout << "; " << NoISet[i];              // Nor has this I been set
        if(j>=3) j=0;
        }
      }  
    error(3, DefIso);				// Say we've set type to DefIso
    }
  delete [] NoISet;
  }
     

void SpinSpace::operator= (const ParameterSet& pset) { setSsys(pset); }

	// Input		sys      : Spin system (this)
	// 			pset     : A parameter set
	// Output		none	 : Spin system filled with
	//				   parameters n pset
	// Note				 : Three things are gleaned from
	//				   the parameter set for a SpinSpace
	//				   1.) The number of spins
	//				   2.) An isotope type for each spin
	//				   3.) An optional spin system name
	// Note				 : Functions which place a SpinSpace
	//				   into a parameter set must contain
	//				   add the information read here


// ----------------------------------------------------------------------------
//      Functions To Output Spin System To ASCII From A Parameter Set
// ----------------------------------------------------------------------------


int SpinSpace::write(const string &filename, int idx, int warn) const
 
	// Input		ss	: Spin system (base)
        //                      filename: Output file name
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning level
	// Output		none	: Spin system is written as a 
	//				  parameter set to file filename

  {
  if(!spins()) return 1;		// Nothing if no spins
  ofstream ofstr(filename.c_str());	// Open filename for input
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!write(ofstr, idx, w2))            // If file bad then exit
    {
    if(warn)
      {  
      error(1, filename);		// Problems with file
      if(warn>1) fatality(5);		// Fatal error
      }  
    return 0;
    }
  ofstr.close();                        // Close it now
  return 1;
  }


int SpinSpace::write(ofstream& ofstr, int idx, int warn) const

	// Input		ss	: Spin system (base)
        //                      ofstr   : Output file stream
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning level
        // Output               none    : Spin system is written as a
        //                                parameter set to output filestream
        // Note                         : This depends on function PSetAdd!
 
  {
  if(!spins()) return 1;                // Nothing if no spins
  if(!ofstr.good())                     // If file bad then exit
    {
    if(warn) error(6);			//      Problems with filestream
    if(warn > 1) fatality(7);		//      It's a fatal error
    return 0;
    }
  ParameterSet pset;			// Declare a parameter set
  PSetAdd(pset, idx);                   // Add system to parameter set
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!pset.write(ofstr, w2))            // Use parameter set to write
    {
    if(warn)
      {
      error(6, 1);			// Problems writing to filestream
      if(warn>1) fatality(7);		// Fatal error
      }  
    return 0;
    }
  return 1;
  }


// ____________________________________________________________________________
//                          SPIN SPACE INPUT FUNCTIONS
// ____________________________________________________________________________


int SpinSpace::read(const string& filename, int idx, int warn)

	// Input		sys      : Spin system (this)
	// 			filename : Input filename
	//			idx	 : Parameter index value used for
	//				   prefix [#] in input names
	//			warn	 : Warning output level
	//					0 = no warnings
	//					1 = warnings
	//				       >1 = fatal warnings
	// Output		TF	 : Spin system is filled with
	//				   parameters read from file
	//				   TRUE if read is successful
	// Note			 	 : The file should be an ASCII file
	//				   containing recognized sys parameters

  {
  ParameterSet pset;			// Declare a parameter set
  if(!pset.read(filename, warn?1:0))	// Try and read in pset
    {					// If we cannot read the file then
    if(warn)				// we'll issue warnings as desired
      {
      error(1, filename);		//	Problems with file
      if(warn>1) fatality(4);		// 	This is a fatal problem
      if(warn>1) error(4);		// 	Or maybe it isn't so bad...
      }
    return 0;				// Well flag we didn't read & exit
    }
  return read(pset, idx, warn);		// Fill up SpinSpace with parameters
  }


int SpinSpace::read(const ParameterSet& pset, int idx, int warn)

	// Input		sys      : Spin system (this)
	// 			pset	 : A parameter set
	//			idx	 : Parameter index value used for
	//				   prefix [#] in input names
	//			warn	 : Warning output level
	//					0 = no warnings
	//					1 = warnings
	//				       >1 = fatal warnings
	// Output		TF	 : Spin system is filled with
	//				   parameters in pset
	//				   TRUE if system filled properly

  {
  int TF = setSsys(pset, idx, warn?1:0);	// Use overload to read
  if(!TF)					// If setSsys didn't handle
    {						// the system read from pset
    if(warn)					// then we'll issue some
      {
      error(8, 1);				//    Problems with pset
      if(warn>1) fatality(4);			//    This is a fatal problem
      if(warn>1) error(4);			//    Or maybe it isn't so bad..
      }
    return 0;
    }
  return TF;
  }


string SpinSpace::ask_read(int argc, char* argv[], int argn)

	// Input		sys     : A basic spin system (this)
	//			argc	: Number of arguments
	//			argv    : Vector of argc arguments
	//			argn    : Argument index
	// Output		filename: The parameter argn of array argc
	//				  is used to supply a filename
	//				  from which the spin system is read
	//				  If the argument argn is not in argv,
	//				  the user is asked to supply a filename
	//				  The filename is returned
	// Note			 	: The file should be an ASCII file
	//				  containing recognized sys parameters
	// Note			 	: The spin system is modifed (filled)

  {
  string filename;				// Name of spin system file  
  query_parameter(argc, argv, argn,		// Get filename from command
       "\n\tSpin system filename? ", filename);	// Or ask for it
  read(filename);		           	// Read system from filename
  return filename;				// Give back the filename
  }

 
// ____________________________________________________________________________
//                             DEFAULT BASIS FUNCTIONS
// ____________________________________________________________________________


matrix SpinSpace::get_basis() const { return bsmx; }

	// Input		sys	: Spin system
	// Output		bs	: The default basis defined by
	//				  the spin system Hilbert space
	// Note				: The return is matrix rather than
	//				  basis as the latter class is not
	//				  known to GAMMA at this level

// ____________________________________________________________________________
//                           SPIN SPACE MAPPING FUNCTIONS
// ____________________________________________________________________________

/* These functions map the system's single spin & spin pair basis functions as
   well as their coherences to those of the full system.  This is quite useful
   in switching between reduced and full Hilbert spaces.	   	     */

//-----------------------------------------------------------------------------
//                            Single Spin Mappings
//-----------------------------------------------------------------------------


matrix SpinSpace::BasisMap1() const

	// Input		sys      : Spin system (this)
	// Output		bmap	 : A ns x hs array of basis mappings.
	//				   For each spin is a row the length
	//				   of the full Hilbert space.  Each
	//				   column (system basis function) will
	//				   contain the spin's corresponding
	//				   basis function # in the REAL part.
	//				   to the systems basis function #
	//				   The imaginary part will contain the
	//				   total basis Fz - spin's mz 

  {
  int NS = spins();                     // Number of spins
  int FHS = HS();			// Full spin Hilbert space
  matrix QS = qStates();                // Get quantum values
  col_vector QNS = qnStates();          // Get Fz values
  int I;

//           This is For The Spin Pair Sub-Space

  SpinSpace sys1(1);                     // Default for spin pair
  int hs;                               // Hilbert space
  matrix qs;                            // Spin flips (mz) all transitions
  col_vector qns;                       // Vector of Fz totals

  matrix bmap(NS,FHS);			// Basis mapping array
  int BF, bf, match;                    // Basis indices, match flag
  double delFz;                         // Delta Fz for basis functions
  for(I=0; I<NS; I++)                   // Loop over the spins
    {
    sys1.isotope(0, symbol(I));         // Set I isotope type
    qs = sys1.qStates();                // Sub-space states
    qns = sys1.qnStates();              // Sub-space Fz totals
    hs = sys1.HS();
    for(BF=0; BF<FHS; BF++)              // Loop basis functions
      {  
      match = 0;                        // No basis mapping
      for(bf=0; bf<hs && !match; bf++)
        {
        if(QS.getRe(BF,I) != qs.getRe(bf,0))
          bmap.put(0,I,BF);
        else
          { 
          delFz = QNS.getRe(BF) - qns.getRe(bf);
          bmap.put(complex(bf, delFz),I,BF);
          match++;
          }
        }
      }
    }    
  return bmap;
  }


matrix SpinSpace::TransitionMap1() const

        // Input                sys      : Spin system (this)
        // Output               tmap     : An ns x hs array of coherence maps.
        //                                 For each spin is a row the length
        //                                 of the full Liouville space.  Each
        //                                 column (system coherences) will
        //                                 contain the spin's corresponding
        //                                 coherence # in the REAL part.

  {
  matrix QS = qStates();			// Get quantum values
  matrix bmap = BasisMap1();    	        // Single spin basis map
 
  int NS = spins();				// Number of spins       
  int FHS = HS();				// Full spin Hilbert space
  int LS = FHS*FHS;				// Full spin Liouville space
  matrix tmap(NS, LS);				// Transition map
  int H,I;					// Spin indices, full space
  int L,M,N;					// Basis indices, full space
  int m,n,hs;					// Basis indices, sub space
  double delMz;
  SpinSpace sys1(1);
  for(I=0; I<NS; I++)				// Loop over the spins, I
    {
    sys1.isotope(0, symbol(I));			// Set I isotope type
    hs = sys1.HS();				// Get sub-Hilbert space
    for(L=0, M=0; M<FHS; M++)			// Loop all transitions
      {
      m = int(bmap.getRe(I,M));			// Get map of m -> M
      for(N=0; N<FHS; N++, L++)
        {
        n = int(bmap.getRe(I,N));		// Get map of n -> N
        delMz = bmap.getIm(I,M) - bmap.getIm(I,N);
        if(delMz) tmap.put(-1, I, L);		// No map if delMz != 0
        else					// To insure mapping we must
          {					// check all spin flips now
          int map=1;				//   Assume map is O.K.
          for(H=0; H<NS && map; H++)
            {
            if(H!=I && QS.getRe(N,H) != QS.getRe(M,H))
              {
              map = 0;
              tmap.put(-1, I, L);
              }
            }
         if(map)
            {
            tmap.put(complex(m,n),I,L);
            }
          }
        }    
      }
    }    
  return tmap;
  }  


//-----------------------------------------------------------------------------
//                              Spin Pair Mapping
//-----------------------------------------------------------------------------


matrix SpinSpace::BasisMap2() const

	// Input		sys      : Spin system (this)
	// Output		bmap	 : A nd x hs array of basis mappings.
	//				   Each spin pair is a row the length
	//				   of the full Hilbert space.  Each
	//				   column (system basis function) will
	//				   contain the spin pair's corresponding
	//				   basis function # in the REAL part.
	//				   The imaginary part will contain the
	//				   total basis Fz - spin pair's mz 

  {
  int NS = spins();                     // Number of spins
  int FHS = HS();			// Full spin Hilbert space
  int ND = spinpairs();                 // Number of dipoles
  matrix QS = qStates();                // Get quantum values
  col_vector QNS = qnStates();          // Get Fz values
  int I,J,K;
 
//                   This is For The Spin Pair Sub-Space
 
  SpinSpace sys2(2);                     // Default for spin pair
  int hs;                               // Hilbert space
  matrix qs;                            // Spin flips (mz) all transitions
  col_vector qns;                       // Vector of Fz totals
 
  matrix bmap(ND,FHS);			// Basis mapping array
  int BF, bf, match;                    // Basis indices, match flag
  double delFz;                         // Delta Fz for basis functions
  for(I=0, K=0; I<NS-1; I++)            // Loop over the dipoles
    {
    sys2.isotope(0, symbol(I));         // Set I isotope type
    for(J=I+1; J<NS; J++, K++)          // (and their spin pairs)
      {
      sys2.isotope(1, symbol(J));       // Set S isotope type
      qs = sys2.qStates();              // Sub-space states
      qns = sys2.qnStates();            // Sub-space Fz totals
      hs = sys2.HS();
      for(BF=0; BF<FHS; BF++)		// Loop basis functions
        {
        match = 0;                      // No basis mapping
        for(bf=0; bf<hs && !match; bf++)
          {
          if(QS.getRe(BF,I) != qs.getRe(bf,0))
            bmap.put(0,K,BF);
          else if(QS.getRe(BF,J) != qs.getRe(bf,1))
            bmap.put(0,K,BF);
          else
            {
            delFz = QNS.getRe(BF) - qns.getRe(bf);
            bmap.put(complex(bf, delFz),K,BF);
            match++;
            }
          }
        }
      }
    }
  return bmap;
  }


matrix SpinSpace::TransitionMap2() const

        // Input                sys      : Spin system (this)
        // Output               tmap     : An nd x hs array of coherence maps 
        //                                 Each spin pair has a row the length
        //                                 of the full Liouville space.  Each
        //                                 column (system coherences) will
        //                                 contain the spin pair's corresponding
        //                                 coherence # in the REAL part.
	// Note				 : If there is no mapping, the indices
	//				   will be set to -1,0

  {
  matrix QS = qStates();                // Get quantum values
  matrix bmap = BasisMap2();            // Dipolar basis map
 
  int NS = spins();                     // Number of spins
  int ND = spinpairs();                 // Number of dipoles
  int FHS = HS();			// Full spin Hilbert space
  int LS = FHS*FHS;			// Full spin Liouville space
  matrix tmap(ND, LS, complex(-1));	// Transition map
  int H,I,J,K;                          // Spin indices, full space
  int L,M,N;                            // Basis indices, full space
  int m,n,hs;				// Basis indices, sub space
  double delMz;
  SpinSpace sys2(2);
  for(I=0, K=0; I<NS-1; I++)            // Loop over the dipoles, K
    {
    sys2.isotope(0, symbol(I));         // Set I isotope type
    for(J=I+1; J<NS; J++, K++)          // (and their spin pairs)
      {
      sys2.isotope(1, symbol(J));       // Set S isotope type
      hs = sys2.HS();                   // Get sub-Hilbert space
      for(L=0, M=0; M<FHS; M++)		// Loop all transitions
        {
        m = int(bmap.getRe(K,M));       // Get map of m -> M
        for(N=0; N<FHS; N++, L++)
          {
          n = int(bmap.getRe(K,N));     // Get map of n -> N
          delMz = bmap.getIm(K,M) - bmap.getIm(K,N);
          if(delMz) tmap.put(-1, K, L);	// No map if delMz != 0
          else                          // To insure mapping we must
            {                           // check all spin flips now
            int map=1;                  //   Assume map is O.K.
            for(H=0; H<NS && map; H++)
              {
              if(H!=I && H!=J && QS.getRe(N,H) != QS.getRe(M,H))
                {
                map = 0;
                tmap.put(-1, K, L);
                }
              }
            if(map)
              {
              tmap.put(complex(m,n),K,L);
              }
            }
          }    
        }    
      }
    }    
  return tmap;
  }  


// ____________________________________________________________________________
//                            STANDARD I/O FUNCTIONS
// ____________________________________________________________________________


ostream& SpinSpace::print(ostream& ostr) const

	// Input		sys      : Spin system (this)
	//			ostr	 : An output stream
	// Output		ostr	 : The output stream modified by
	//				   the spin system parameters

  {
  string s;
  if(sysname.length()) 				// Output the system name
    ostr << "System   : " << sysname << "\n";	// if there is one
  ostr << "Spin     :";				// Output the spin labels
  int i;
  for(i=0; i<nspins; i++) ostr << Gdec("%10d",i);
  ostr << "\nIsotope  :";			// Output the Isotope types
  for(i=0; i<nspins; i++)
    {
    s = symbol(i);
    ostr << string(10-s.length(), ' ') << s;
    }
  ostr << "\nMomentum :";			// Out the spin angular momentum
  for (i=0; i<nspins; i++)
    {
    s = momentum(i);
    ostr << string(10-s.length(), ' ') << s;
    }
  ostr << "\n";
  return ostr;
  }


ostream& operator<< (ostream& ostr, const SpinSpace& sys)

	// Input		ostr	 : An output stream
	// Input		sys      : Spin system (this)
	// Output		ostr	 : The output stream modified by
	//				   the spin system parameters

  { sys.print(ostr); return ostr; }


vector<string> SpinSpace::printstrings() const

	// Input		sys      : Spin system (this)
        // Output               SV       : string vector containing strings
	//				   which can be used for printing
	//				   the spin system

  {
  vector<string> StrArray;			// Make a string vector
  int i;					// Array index
  string Stemp, s;				// Temporary strings
  if(sysname.length()) 				// Store the system name
    {
    Stemp = string("System   : ") + sysname;
    StrArray.push_back(Stemp);
    }
  Stemp = string("Spin     :");			// String is for spin labels
  for(i=0; i<nspins; i++) Stemp += Gdec("%10d",10);
  StrArray.push_back(Stemp);			// Store spin labels
  Stemp = string("Isotope  :");			// Next string is for the
  for(i=0; i<nspins; i++) 			// the spin isotope types
    {
    s = symbol(i);
    Stemp += string(10-s.length(), ' ') + s;
    }
  StrArray.push_back(Stemp); 			// Store the isotope types
  Stemp = string("Momentum :");			// Next string is for the
  for(i=0; i<nspins; i++)			// the spin angular momentum
    {
    s = momentum(i);
    Stemp += string(10-s.length(), ' ') + s;
    }
  StrArray.push_back(Stemp);			// Store angular momenta
  return StrArray;				// Return the string array
  }

#endif						// SpinSpace.cc
