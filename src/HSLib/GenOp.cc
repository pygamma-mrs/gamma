/* GenOp.cc *****************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      General Operator                            Implementation      **
**                                                                      **
**      Copyright (c) 1990, 1991, 1992                                  **
**      Scott Smith                                                     **
**      Eidgenoessische Technische Hochschule                           **
**      Labor fuer physikalische Chemie                                 **
**      8092 Zuerich / Switzerland                                      **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
**  Class gen_op defines general operators for GAMMA in C++.  It        **
**  contains all of all the algebraic operations (+, -, *, /), most     **
**  useful complex functions (exp, log) and input/output routines.      **
**                                                                      **
*************************************************************************/

#ifndef   GenOp_cc_			// Is file already included?
#  define GenOp_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// this is the implementation
#  endif

#include <GamGen.h>			// Include OS specific stuff
#include <HSLib/GenOp.h>		// Include the interface
#include <Basics/Gutils.h>		// Include GAMMA std errors
#include <Basics/StringCut.h>	 	// Include GAMMA Gdec function

//# define OpDebug 1			// Define this for debugging

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

int gen_op::DBPr     = 499;			// Default Basis Priority 
int gen_op::EBPr     = 500;			// EigenBasis Priority
int gen_op::MaxReps = 3;			// Maximum Operator Reps

// ____________________________________________________________________________
// i                CLASS GENERAL OPERATOR ERROR HANDLING
// ____________________________________________________________________________

/*      Input               Op      : General operator (this)
                            eidx    : Error index
                            noret   : Flag for linefeed (0=linefeed)
                            pname   : string in message
        Output              none    : Error message output
                                      Program execution stopped if fatal

The following error messages use the defaults set in the Gutils package
 
               Case                          Error Message
 
 NO PNAME       (0)                     Program Aborting.....
    		(9)                     Problems During Construction
                default                 Unknown Error 
 WITH PNAME     (1)                     Problems With File PNAME
                (2)                     Cannot Read Parameter PNAME
    		(3)  			Invalid Use Of Function PNAME
                default                 Unknown Error - PNAME                */  
 
void gen_op::GenOperror(int eidx, int noret) const
  {
  std::string hdr("General Operator");
  std::string msg;
  switch (eidx)
    {
    case 6:  GAMMAerror(hdr,"Bad Internal Component Access",noret);break;//(6)
    case 7:  GAMMAerror(hdr,"Accessing Empty Operator",noret);    break;// (7)
    case 20: GAMMAerror(hdr,"Unable To Perform Addition",noret);  break;// (20)
    case 21: GAMMAerror(hdr,"Cannot Perform Subtraction",noret);  break;// (21)
    case 22: GAMMAerror(hdr,"Unable To Do Multiplication",noret); break;// (22)
    case 39: GAMMAerror(hdr,"Bad Matrix/Basis Mixing",noret);     break;// (39)
    case 40: GAMMAerror(hdr,"Trouble Mixing Two Operators",noret);break;// (40)
    case 41: GAMMAerror(hdr,"Trouble Mixing With Matrix",noret);  break;// (41)
    case 42: GAMMAerror(hdr,"Problems During Assignment",noret);  break;// (42)
    case 50: GAMMAerror(hdr,"Rectangular Array Use",noret);       break;// (50)
    case 51: GAMMAerror(hdr,"Mixing Mismatched Dimensions",noret);break;// (51)
    case 56: GAMMAerror(hdr,"Element Access Out Of Range",noret); break;// (56)
    case 57: GAMMAerror(hdr,"Min. Representations Limit 3",noret);break;// (57)
    case 67: GAMMAerror(hdr,"Cannot Access Element",noret);       break;// (67)
    case 77: GAMMAerror(hdr,"Problems Setting Basis",noret);      break;// (77)
    case 92: GAMMAerror(hdr,"Representation Sizes Differ",noret); break;// (92)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }  

void gen_op::GenOperror(int eidx, const std::string& pname, int noret) const
  {
  std::string hdr("General Operator");
  std::string msg;
  switch(eidx)
    {
    case 5: msg = std::string("Bad Use Of ") + pname + std::string(" Function ");
             GAMMAerror(hdr,msg+pname,noret);  break;		       // (5)
    case 9 :msg = std::string("Construction With ") + pname + std::string("OpReps");
             GAMMAerror(hdr,msg+pname,noret);  break;		       // (9)
    default: GAMMAerror(hdr, eidx, pname, noret); break;// Unknown Error  (-1)
    }
  }
     
volatile void gen_op::GenOpfatality(int eidx) const
  {
  GenOperror(eidx, 1);				// Normal non-fatal error
  if(eidx) GenOperror(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

 
volatile void gen_op::GenOpfatality(int eidx, const std::string& pname) const
  {
  GenOperror(eidx, pname, 1);			// Normal non-fatal error
  if(eidx) GenOperror(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii               GENERAL OPERATOR COMMON INTERNAL FUNCTIONS
// ____________________________________________________________________________

/* The next two functions provide access to the first operator representation
   in the operator list. These are necessary because fo the way the STL handles
   vectors and the fact that we have mutable functions.

   In this class

   1. WBR, EBR, and DBR are mutable pointers to genoprep
   2. STL vector functions end() & begin() are (const) iterators into the
      vector of operator representations vector<genoprep> for the operator
   3. *(end()) is a (const) genoprep
   4. &(end()) is the (const) address of genoprep
   5. *(&end()) is a (const) pointer to genoprep
   6. Obegin()  will return a pointer to the 1st Op representation
   7. ObeginC() will return a constant pointer to the 1st Op representation  */
 
      genoprep* gen_op::Obegin()        { return &(*(begin())); }
const genoprep* gen_op::ObeginC() const { return &(*(begin())); }


/* These functions help in manipulating the operator representations.  They
   individually don't worry about Op integrity so they are PRIVATE and must
   be used with caution. Op integrity demands that if any representations
   exist then WBR should point to one of the, whereas if no representations
   exist WBR,EBR,&DBR should all be NULL.

 Function              Action                          Precautions
 --------  --------------------------------- -----------------------------
  ZeroOp   Zeroes OpReps, Nulls All Pointers This should always be fine
  setNULL  Nulls All Pointers                OpReps Can Exist, But NO WBR!
  AddRep   Adds New OpRep, Sets WBR EBR etc. Pre-Set WBR,DBR,EBR as needed
  AddRepM  Constant AddRep! Mutable!         As AddRep, but very tricky.
  GetIndex Returns Index of Specified OpRep  Returns -1 if no OpRep found    */

void gen_op::ZeroOp()        { if(size()) clear(); setNULL(); }
void gen_op::setNULL()const  { WBR=NULL; DBR=NULL; EBR=NULL; }

void gen_op::AddRep(const gen_op& Op1)
  {
  if(!(Op1.WBR)) return;			// Do nothing if Op1 NULL
  AddRep(*(Op1.WBR));				// Use function overload
  }


void gen_op::AddRep(const genoprep& OpRep)
  {
# ifdef OpDebug
  std::cout << "\n\tAdding New Rep To ";
  if(OpName.length()) std::cout << OpName << " ";
  std::cout << "Operator";
# endif
  int DBI = GetIndex(DBR); 			// Index of DBR
  int EBI = GetIndex(EBR); 			// Index of EBR
  push_back(OpRep);				// Add OpRep to rep
  std::vector<genoprep>::iterator item;		// An iterator into the list
  item = end()-1;				// The last in the list
  WBR = &(*item);				// This will be the WBR now
  if(DBR && DBI>=0) { DBR=Obegin(); DBR+=DBI; }	// Reset DBR or see if there
  else                check_DBR();		// is a DBR if none yet exists
  if(EBR && EBI>=0) { EBR=Obegin(); EBR+=EBI; }	// Reset EBR or see if there
  else                check_EBR();		// is an EBR if none yet exists
  }						// (No Check of Rep Limits)					

void gen_op::AddRepM(const genoprep& OpRep) const
  {
# ifdef OpDebug
  std::cout << "\n\tAdding New Rep To Mutable ";
  if(OpName.length()) std::cout << OpName << " ";
  std::cout << "Operator";
# endif
  int DBI = GetIndex(DBR); 				// Initial index of DBR
  int EBI = GetIndex(EBR); 				// Initial index of EBR
  gen_op* OpTmp = const_cast<gen_op*>(this);		// Cast away const
  (*OpTmp).push_back(OpRep);				// Add OpRep to rep
  std::vector<genoprep>::iterator item;			// The last in the list
  item = (*OpTmp).end()-1;				// Set WBR to this rep
  WBR = &(*item);					// This is the WBR now
  if(DBR && DBI>=0) { DBR=(*OpTmp).Obegin(); DBR+=DBI; }// Reset DBR or check 
  else                check_DBR();			// for DBR if isnt one
  if(EBR && EBI>=0) { EBR=(*OpTmp).Obegin(); EBR+=EBI; }// Reset EBR or check
  else                check_EBR();			// for EBR if isnt one
  }							// (Dont Chk Rep Limits)					

int gen_op::GetIndex(const genoprep* OpRep) const
  {
  if(!size()) return -1;				// -1 if no reps
  std::vector<genoprep>::const_iterator item;		// Iterator into vector
  item = begin();					// First in the vector
  int i=0, idx=-1;					// Temp counters
  while(item!=end() && idx==-1)
    {
    if(&(*item) == OpRep) idx = i;			//   Save indx if found
    item++;						//   Point to next
    i++;						//   Update counter
    }
  return idx;						// This is rep index
  }

// ____________________________________________________________________________
// iii           INTERNAL OPERATOR REPRESENTATION MANIPULATIONS
// ____________________________________________________________________________
 
// *******************  Special Basis Representation Checks *******************
 
/* These functions look for a particular OpRep and will set the WBR to that
   representation if found.  They do not alter the operator in any other way.
   Note that these are constant functions because only the OpRep that WBR
   points to is changed, i.e. WBR may be set to point to a new OpRep but the
   OpRep isn't altered (except to maybe change a priority value).  This is
   possible because WBR was declared mutable in the class definition.

   Functions  Argument              Result                       Caution
   ---------  --------  --------------------------------  ---------------------
   check_DBR    Op1     If Op1 DBR=WBR, sets Op DBR=WBR   NULL WBR's,DBR's Bad?
   check_EBR    Op1     If Op1 WBR=EBR, sets Op EBR=WBR   NULL WBR's,EBR's Bad?
   check_DBR     -      If Op WBR tests DBR, Set DBR=WBR  WBR MUST NOT BE NULL
   check_EBR   cutoff   If Op WBR tests EBR, Set EBR=WBR  WBR MUST NOT BE NULL

   Why are these private? For one thing, it isn't mandatory that a DBR or
   EBR exists for any particular operator (same for WBR in a NULL operator)
   so accessing OpReps through these pointers could lead to big trouble. But
   also, these will over write any existing OpRep pointer without checking.
   Use caution! Mostly these are called after construction, assignment, or
   an operation that leaves the operator with only 1 OpRep, WBR.             */

void gen_op::check_DBR(const gen_op &Op1) const {if(Op1.WBR==Op1.DBR) DBR=WBR;}
void gen_op::check_EBR(const gen_op &Op1) const {if(Op1.WBR==Op1.EBR) EBR=WBR;}

void gen_op::check_DBR() const
  {						// See if WBR basis is default
  if(WBR->RepBs.isDefaultBasis())		// (provided in class basis)
    {
    DBR = WBR;					// Set DBR equal to WBR
    if(DBPr>WBR->RepPty) WBR->RepPty=DBPr;	// Set DBR priority >= WBR
    }
  }

void gen_op::check_EBR(double cutoff) const
  {				
  if(WBR->RepMx.test_type(d_matrix_type,cutoff) // Test if WBR mx is diagonal
                               == d_matrix_type)
    {
    WBR->RepMx.set_type(d_matrix_type);		// Set the matrix diagonal
    EBR = WBR;					// Set EBR equal to WBR
    if(EBPr>WBR->RepPty) WBR->RepPty=EBPr; 	// Set WBR priortty >= EBR
    }
  }
 
// ********************* Finding Specific Representations *********************

/* Functions  Argument              Result                       Caution
   ---------  --------  --------------------------------  ---------------------
    FindRep      bs     Pointer to OpRep with basis bs     Return NULL !found
    SetRep       bs     Set WBR to OpRep with basis bs     Return 0 if unable
    LP_rep       -      Pointer to OpRep with low prior.         -----
 
 Note that LPR, the OpRep of lowest priority, is deleted in delete_rep. We hope 
 that higher priority reps are stored at list top from IBR and lower priority
 reps at the bottom of the list. WBR is preferentially NOT chosen as LPR.    */
 
const genoprep* gen_op::FindRep(const basis &bs) const
  {
  const genoprep* REP = NULL;		// Don't know proper OpRep
  std::vector<genoprep>::const_iterator item;// Iterator into list
  item = begin();			// Begin with 1st OpRep
  const genoprep* PTR = ObeginC();	// Begin with 1st OpRep
  while(item!=end() && !REP)		// Loop through OpReps
    {
    if(PTR->RepBs == bs) REP=PTR;	//   Set pointer if bases match
    PTR++;
    item++;
    }
  return REP;
  }

int gen_op::SetRep(const basis &bs) const
  {
  if(!size()) return 0;				// Do nothing if no reps
  int TF = 0;					// Haven't found OpRep
  gen_op* OpTmp = const_cast<gen_op*>(this);	// Cast away const
  std::vector<genoprep>::const_iterator item;	// Iterator into list
  item = (*OpTmp).begin();
  genoprep* PTR = (*OpTmp).Obegin();		// Begin with 1st OpRep
  while((item != ((*OpTmp).end())) && !TF)	// Loop through OpReps
    {
    if(PTR->RepBs == bs) 			//   See if RepBs is bs
      {						//   If so, PTR is the OpRep
      WBR=PTR;					//   we want as our WBR, so
      TF++;					//   set it & flag we found it
      }
    item++;
    PTR++;					// Move to next OpRep
    }
  return TF;
  }

const genoprep* gen_op::LP_rep() const
  {
  const genoprep* LPR = NULL;		// Don't know proper OpRep
  int Pty = 500000;			// Start with large priority
  std::vector<genoprep>::const_iterator item;// Iterator into list
  item = begin();
  const genoprep* PTR = ObeginC();	// Begin with 1st OpRep
  while(item != end())			// Loop through OpReps
    {
    if((PTR->RepPty)<Pty && PTR!=WBR)	// Lower priority rep found
      {					// That is NOT the WBR
      Pty = PTR->RepPty;		//   Store lowest priority
      LPR = PTR;			//   Pointer to low priority rep
      }
    PTR++;
    item++;
    }
  return LPR;
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A              CLASS GENERAL OPERATOR CONSTRUCTORS/DETRUCTOR
// ____________________________________________________________________________

 
/* These constructors set up a new operator.  There are several ways to make
   an operator as listed below:
 
             Input Arguments                Resulting Operator
          ---------------------    -----------------------------------------
                   -               Empty Operator, No Representations
                  mx               Op With 1 Representation, mx in DBR
                  SOp              Op With 1 Representation, SOp's mx in DBR
                mx1,mx2            Op With 1 Representation, {mx1,mx2}
                 mx,bs             Op With 1 Representation, {mx,bs} 
          vector<mx>,vector<bs>    Op With 1 Representation, {MXCmp,BSCmp} 
               *mx,n,*bs           Op With 1 Representation, {MXCmp,BSCmp} 
               N,*mx,*bs           Op With N Representations {mx[i],bs[i]} 
 
   Above, MXCmp & BSCmp are larger composite space arrays formed by taking 
   smaller space arrays/bases and placing them on the diagonal of the    
   larger array/basis to result in a block-diagonal form.  This supports
   use of multi-system spin systems in GAMMA.                                */

gen_op::gen_op():std::vector<genoprep>() {setNULL();}// Absolutely NO representions

gen_op::gen_op(const matrix &mx) : std::vector<genoprep>()
  {
  if(!OpCheck(mx, 1)) GenOpfatality(9);		// Insure array is square
  setNULL(); 					// Insure no rep pointers
  if(!mx.cols()) return; 			// Empty Op if no array
  AddRep(genoprep(mx, basis(mx.rows()), DBPr));	// else set mx to 1st OpRep
  }

gen_op::gen_op(const spin_op& SOp) : std::vector<genoprep>()
  { 
  matrix mx = SOp.get_mx();
  *this = gen_op(mx);
  }
//  { *this = gen_op(SOp.get_mx()); }

gen_op::gen_op(const matrix &mx1, const matrix &mx2) : std::vector<genoprep>()
  {
  if(!OpCheck(mx1,mx2,1)) GenOpfatality(9);	// Insure arrays match
  setNULL(); 					// Insure no rep pointers
  if(!mx1.cols()) return; 			// Empty Op if no array
  AddRep(genoprep(mx1, basis(mx2), DBPr));	// Set matrix to 1st OpRep
  }

gen_op::gen_op(const matrix &mx, const basis &bs) : std::vector<genoprep>()
  {
  if(!OpCheck(mx,bs,1)) GenOpfatality(9);	// Insure mx & bs match
  setNULL(); 					// Insure no rep pointers
  if(!mx.cols()) return; 			// Empty Op if no array
  AddRep(genoprep(mx, bs, DBPr));		// Set matrix to 1st OpRep
  }

gen_op::gen_op(const gen_op& Op1) : std::vector<genoprep>()
  { setNULL(); OpName = Op1.OpName; AddRep(Op1); }
 
gen_op::gen_op(const std::vector<matrix>& mxc, const std::vector<matrix>& bsc)
       :std::vector<genoprep>()
  {
  int dimM = 0;					// Total Op dimension
  int nc = mxc.size();				// Number of components
  int* ncd;
  ncd = new int[nc];				// Array for dimensions
  int i=0;
  for(i=0; i<nc; i++)				// Loop input components
    {
    dimM += mxc[i].rows();			//   Adjust Op dimension
    if(!(OpCheck(mxc[i], bsc[i]),1))		//   Make sure bsc & mxc
      GenOpfatality(9);				//   properly dimensioned
    }
  matrix mx(dimM, dimM, i_matrix_type);		// For our Op rep matrix
  matrix bs(dimM, dimM, i_matrix_type);		// For our Op rep basis
  int pinblk = 0;				// Start at <0|mx|0>,<0|bs|0>
  for(i=0; i<nc; i++)				// Loop over our components
    {
    ncd[i] = mxc[i].rows();			// Set component dimension
    mx.put_block(pinblk, pinblk, mxc[i]);	// Put mxc component into mx 
    bs.put_block(pinblk,pinblk,bsc[i]);		// Put bsc component into bs
    pinblk += mxc[i].rows();			// Where next component goes
    }						// i.e. <pb|mx|pb> & <pb|bs|pb>
  AddRep(genoprep(mx,basis(bs,nc,ncd),DBPr));	// Set first Op representation
						// in the composite space
  delete [] ncd;				// Remove array of dimensions
  }


 
        // Input                nc  : Number of components in multi-sys
        //                      mxc : Operator matrices associated with comps
        //                      bsc : Basis matrices associated with components
        // Output               Op  : General operator (this) containing
        //                            mxc as diagonal blocks of mx, and bsc
        //                            as diagonal blocks of bs
	// Note			    : bsc array can be empty to trigger use
	//			      of a default basis (I matrix)
 
gen_op::gen_op(matrix* mxc, int nc, matrix* bsc) : std::vector<genoprep>()
  {
  int dimM = 0;					// Total Op dimension
  int* ncd;					// Array for dimensions
  ncd = new int[nc];
  int i=0;
  for(i=0; i<nc; i++)				// Loop input components
    {
    dimM += mxc[i].rows();			//   Adjust Op dimension
    if(bsc && !(OpCheck(mxc[i], bsc[i]),1))	//   Make sure bsc & mxc
      GenOpfatality(9);				//   properly dimensioned
    }
  matrix mx(dimM, dimM, i_matrix_type);		// For our Op rep matrix
  matrix bs(dimM, dimM, i_matrix_type);		// For our Op rep basis
  int pinblk = 0;				// Start at <0|mx|0>,<0|bs|0>
  for(i=0; i<nc; i++)				// Loop over our components
    {
    ncd[i] = mxc[i].rows();			// Set component dimension
    mx.put_block(pinblk, pinblk, mxc[i]);	// Put mxc component into mx 
    if(bsc) bs.put_block(pinblk, pinblk, bsc[i]);// Put bsc component into bs
    pinblk += mxc[i].rows();			// Where next component goes
    }						// i.e. <pb|mx|pb> & <pb|bs|pb>
  AddRep(genoprep(mx,basis(bs,nc,ncd),DBPr));	// Set first Op representation
  delete [] ncd;
  }


gen_op::gen_op(int N, matrix* mxs, matrix* bss) : std::vector<genoprep>()

        // Input                N    : Number of representations
        //                      mxs  : Vector of operator rep. matrices
        //                      bss  : Vector of operator rep. bases
        // Output               Op   : General operator (this) constructed
        //                             with N representations from {mxs,bss}
        // Note                      : First will be set to WBR.
        // Note                      : If Op1 is NULL Op will be NULL.
	// Note			     : This makes an operator with  multiple 
	//			       representations.  It does NOT make
	//			       a single operator in a composite space

  {
  if(N<0) GenOpfatality(9, Gdec(N));		// Negative # of reps, die;
  if(!mxs[0].cols()) return; 			// Empty Op if no array
  setNULL(); 					// Insure no rep pointers
  AddRep(genoprep(mxs[0],basis(bss[0]),DBPr));	// Start first representation
  for(int i=1; i<N; i++)			// Now add in other reps
    {
    if(!OpCheck(mxs[i], bss[i], 1))		// Insure acceptable pair
      GenOpfatality(9);				// else unable to construct
    if(dim() != mxs[i].rows())			// Insure space is same as other
      { 					// representations
      GenOperror(92);				//   OpRep sizes differ!
      GenOpfatality(9); 			//   Can't construct
      }
    AddRep(genoprep(mxs[i],basis(bss[i]),DBPr));	// Add this representation
    }
  }

gen_op::~gen_op () { }			// Everyone deletes themselves


// ____________________________________________________________________________
// B               OPERATOR FUNCTIONS, OPERATOR WITH OPERATOR
// ____________________________________________________________________________

/*     Operator    Arguments                           Result 
   --------------- ---------   ------------------------------------------------ 
          +         Op1,Op2    Op1+Op2 in WBR of Op1
          +=        Op1,Op2    Op1 = Op1+Op2 in WBR of Op1
          -            Op1     -Op1 in WBR of Op1
          -         Op1,Op2    Op1-Op2 in WBR of Op1
          -=        Op1,Op2    Op1 = Op1-Op2 in WBR of Op1
          *         Op1,Op2    Op1*Op2 in WBR of Op1
          *=        Op1,Op2    Op1 = Op1*Op2 in WBR of Op1
          &=        Op1,Op2    Op1 = Op2*Op1 in WBR of Op1
*/

gen_op gen_op::operator+ (const gen_op &Op2) const
  {
  if(!WBR) return Op2; 				// Op1 NULL, result is Op2
  if(!Op2.WBR) return *this;			// Op1 NULL, result is Op1
  if(!OpCheck(Op2,1))				// Check dimensions
    {
    GenOperror(40, 1); 				// Error during Op-Op operation 
    GenOpfatality(20); 				// Error during addition
    }
  Op2.Op_base(*this);				// Put Op2 in WBR of Op1
  return gen_op(WBR->RepMx + Op2.WBR->RepMx, WBR->RepBs);// Return Op1 + Op2
  }

gen_op & gen_op::operator+= (const gen_op &Op1)
  {
  if(!Op1.WBR) return (*this);			// Don't add if Op1 is NULL
  if(!WBR) { *this=Op1; return (*this); }	// If Op is NULL, result is Op1
  if(!OpCheck(Op1, 1))			// Insure Op & Op1 compatible
    {
    GenOperror(40, 1);		 	// Error during Op-Op operation 
    GenOpfatality(20);			// Error during addition
    }
  Op_base(Op1);				// Insure Op in basis of Op1
  setOnlyWBR();				// Insure only one representation
  WBR->RepMx += Op1.WBR->RepMx;		// Add arrays in mutual WBR basis
  OpName = std::string("");		// No longer has a name
  if(EBR)				// Check if eigenbasis destroyed
  { if(Op1.EBR==NULL)
    { EBR = NULL;	          	//   Can't be EBR if Op1 had no EBR 
    }
    else
    { if(Op1.EBR->RepBs != EBR->RepBs)//  NM
        EBR = NULL;
    }
  }
  // sosiz - Must Adjust The Representation Priority!
  return (*this);
  }

gen_op gen_op::operator- (const gen_op &Op2) const
  {
  if(!Op2.WBR) return *this;			// If NULL Op2, return Op1
  if(!WBR)					// If NULL Op1, return -Op2
    return gen_op((Op2.WBR->RepMx).operator-(), Op2.WBR->RepBs);
  if(!OpCheck(Op2,1))				// Check Op1-Op2 dimensions
    {
    GenOperror(40, 1);   			// Error in Op-Op workings
    GenOperror(21);				// Can't do the subtraction
    }
  Op2.Op_base(*this);				// Put Op2 in WBR of Op1
  return gen_op(WBR->RepMx - Op2.WBR->RepMx, WBR->RepBs);// Return Op1 - Op2
  }

gen_op gen_op::operator-() const
  {
  if(!WBR) return *this;				// If no WBR return ourself
  gen_op NegOp((WBR->RepMx).operator-(), WBR->RepBs);	// Set up -Op1
  if(OpName.length()) NegOp.name("Negated " + OpName);	// Set new name
  return NegOp;
  }

/*
gen_op operator - (const gen_op &Op1)
  {
  if(!Op1.WBR) return Op1;			// Just return if NULL Op1
  return gen_op((Op1.WBR->RepMx).operator-(), Op1.WBR->RepBs);	// Else return -Op1
  }
*/

gen_op & gen_op::operator -= (const gen_op &Op1)
  {
  if(!Op1.WBR) return (*this);			// Do nothing if Op1 NULL
  if(!WBR) { *this = -Op1; return (*this); }	// If No Op, result is -Op1 
  if(!OpCheck(Op1,1))			// Check dimensions
    {
    GenOperror(40, 1);		 	// Error during Op-Op operation 
    GenOperror(21); 			// Can't do the subtraction
    }
  Op_base(Op1);				// Put Op into WBR of Op1
  setOnlyWBR();				// Put Op exclusively in WBR
  WBR->RepMx -= Op1.WBR->RepMx;	
  if(EBR)				// Check if eigenbasis destroyed
  { if(Op1.EBR==NULL)
      {
      EBR = NULL;
      WBR->RepPty = Op1.WBR->RepPty;
      }
    else if(Op1.EBR->RepBs != EBR->RepBs)
      {
      EBR = NULL;
      WBR->RepPty = Op1.WBR->RepPty;
      }
  }
  OpName = std::string("");
  return (*this);
  }


// sosi				Need EBR adjust like *= ?
gen_op gen_op::operator* (const gen_op& Op) const
  {
  if(!Op.WBR || !WBR) 				// If Op1 or Op NULL 
    return gen_op();				// Return a NULL result
  if(!OpCheck(Op,1))				// Check dimensions
    {
    GenOperror(40, 1); 				// Error in Op-Op mix
    GenOpfatality(22);				// Can't do multiplication
    }
  Op.Op_base(*this);				// Put Op in WBR of Op1
  gen_op PdtOp(WBR->RepMx * Op.WBR->RepMx, WBR->RepBs);
  PdtOp.OpName = std::string("");
  return PdtOp;
  }

gen_op & gen_op::operator *= (const gen_op &Op1)
  {
  if(!Op1.WBR || !WBR) 			// If Op or Op1 NULL, result NUll
    { *this=gen_op(); return (*this); }
  if(!OpCheck(Op1,1))			// Check dimensions
    {
    GenOperror(40, 1);			// Error Op-Op use
    GenOpfatality(22);			// Can't do multiplication
    }
  Op1.Op_base(*this);			// Put Op1 into WBR of Op
  setOnlyWBR();				// Delete all reps by WBR of Op
  WBR->RepMx *= Op1.WBR->RepMx;		// Multiply operator matrices
  if(EBR)				// Check if eigenbasis destroyed
    {
    if(Op1.EBR==NULL)
      {
      EBR = NULL;
      WBR->RepPty = Op1.WBR->RepPty;
      }
    else if(Op1.EBR->RepBs != EBR->RepBs)
      {
      EBR = NULL;
      WBR->RepPty = Op1.WBR->RepPty;
      }
    }
  OpName = std::string("");

  return (*this);
  }

gen_op & gen_op::operator &= (const gen_op &Op1)
  {
  if(!Op1.WBR || !WBR) 			// If Op or Op1 NULL, result NUll
    { *this=gen_op(); return (*this); }
  if(!OpCheck(Op1,1))			// Check dimensions
    {
    GenOperror(40, 1); 			// Error during Op-Op operation 
    GenOpfatality(22);			// Can't do multiplication
    }
  Op_base(Op1);				// Put Op1 into WBR of Op
  setOnlyWBR();				// Delete all reps by WBR of Op
  WBR->RepMx = Op1.WBR->RepMx*WBR->RepMx;	// Apply "reverse" mutiplication
  OpName = std::string("");

  return (*this);
  }



	// Input		Op1  : General operator.
	// 			Op   : General operator (this).
	// Return		Op   : Operator which is a copy of
	//			       the input operator, Op = Op1.
	// Note		             : Result EXCLUSIVELY in WBR of Op1

void gen_op::operator= (const gen_op& Op1)
  {
  if(this == &Op1) return; 			// Avoid self copy
  OpName = Op1.OpName;				// Copy any oprator name
/*
  std::vector<genoprep>::operator= (Op1);	// Copy all representations
  DBR = Op1.DBR;				// Copy pointer to DBR
  EBR = Op1.EBR;				// Copy pointer to EBR
  WBR = Op1.WBR;				// Copy pointer to WBR
*/
// sosi replaced these two with last 4 lines
  ZeroOp();					// Zero current Op
  AddRep(Op1);
  }


/* These functions blend operators and spin operators.  Their mixing demands
   it to be done in the default basis, so the results will be exclusively in
   DBR. These functions just use the Op-mx functions.                        */

//void gen_op::operator= (const spin_op &SOp) { (*this) = SOp.get_mx(); }
void gen_op::operator= (const spin_op &SOp)
 {
 matrix mx = SOp.get_mx();
 (*this) = mx; 
 }

gen_op operator + (const gen_op& Op1, const spin_op& SOp)
  { gen_op Op(Op1); Op += SOp; return Op; }

gen_op & gen_op::operator += (const spin_op &SOp)
//  { matrix mx = SOp.get_mx(); (*this) += mx; }
  {
  matrix mx = SOp.get_mx();
  (*this) += mx;

  return (*this);
  }

gen_op operator - (const gen_op& Op1, const spin_op& SOp)
  { gen_op Op(Op1); Op -= SOp; return Op; }

gen_op & gen_op::operator -= (const spin_op &SOp) 
{
  (*this) -= SOp.get_mx();

  return (*this);
}


// ____________________________________________________________________________
// C                OPERATOR FUNCTIONS, OPERATOR WITH MATRIX
// ____________________________________________________________________________

/* These functions blend operators with matrices. The arrays must be square
   and match the dimension of the operartor with which it is interacting. All
   arrays are assumed equal to an operator in the default basis. As such, the
   result of matrix-operator interactions is always in the default basis. 

       Operator    Arguments                           Result 
   --------------- ---------   ------------------------------------------------ 
          +          Op,mx     Return Op1 = Op+mx            in DBR of Op
          +          mx,Op     Return Op1 = mx+Op            in DBR of Op
          +=         Op,mx     Set Op = Op+mx                in DBR of Op
          -          Op,mx     Return Op1 = Op-mx            in DBR of Op
          -          mx,Op     Return Op1 = mx-Op            in DBR of Op
          -=         Op,mx     Set Op = Op-mx                in DBR of Op
          *          Op,mx     Return Op1 = Op*mx            in DBR of Op
          *          mx,Op     Return Op1 = mx*Op            in DBR of Op
          *=         Op,mx     Set Op = Op*mx                in DBR of Op
          &=         Op,mx     Set Op = mx*Op                in DBR of Op
          =          Op,mx     Set Op = mx                   in DBR          */

/*
gen_op gen_op::operator+ (const matrix& mx) const
  {
  if(!WBR) return gen_op(mx);			// Return mx if NULL Op1
  if(!OpCheck(mx,1))				// Insure mx, Op compatible
    {
    GenOperror(41, 1); 				// Error during Op-mx operation 
    GenOpfatality(20);				// Can't do the addition
    }
  set_DBR();					// Put Op1 in default basis
  gen_op AddOp(WBR->RepMx+mx, WBR->RepBs);	// Add the arrays
  AddOp.OpName = std::string("");		// No name now
  return AddOp;					// Return sum
  }

gen_op operator + (const matrix& mx,  const gen_op& Op1) { return Op1 + mx; }

*/

gen_op operator + (const gen_op &Op1, const matrix &mx)  { return mx + Op1; }
gen_op operator + (const matrix &mx,  const gen_op &Op1)
  {
  if(!Op1.WBR) return gen_op(mx);		// Return mx if NULL Op1
  if(!Op1.OpCheck(mx,1))			// Insure mx, Op compatible
    {
    Op1.GenOperror(41, 1); 			// Error during Op-mx operation 
    Op1.GenOpfatality(20);			// Can't do the addition
    }
  Op1.set_DBR();				// Put Op1 in default basis
  return gen_op(Op1.WBR->RepMx+mx, Op1.WBR->RepBs);	// Form return operator
  }

void gen_op::operator += (const matrix& mx)
  {
  if(!WBR) { *this = mx; OpName=""; return; }	// If NULL, result is mx
  if(!OpCheck(mx, 1))				// Check Op, mx dimensions
    {
    GenOperror(41, 1); 				// Error during Op-mx operation 
    GenOpfatality(20);				// Can't do the addition
    }
  set_DBR();					// Set Op to DBR
  setOnlyWBR();					// Set Op EXCLUSIVELY to DBR
  WBR->RepMx += mx;				// Matrix addition
  EBR = NULL;					// Remove prior reference to EBR
  check_EBR();					// Perhaps WBR is now EBR
  if(!EBR) WBR->RepPty = DBPr;			// If not set, set DBPr priority
  OpName = std::string("");
  }

/*
gen_op gen_op::operator- (const matrix& mx) const
  {
  if(!WBR) return gen_op(mx.operator-()); 	// If Op1 NULL, ret. negated mx
  if(!OpCheck(mx, 1))				// Check Op1-mx compatibility
    {
    GenOperror(41, 1); 				// Error during Op-mx operation 
    GenOpfatality(21);				// Can't do the subtraction
    }
  set_DBR();					// Put Op1 in default basis
  gen_op DiffOp(WBR->RepMx-mx, WBR->RepBs);	// Add the arrays
  DiffOp.OpName = std::string("");		// No name now
  return DiffOp;				// Return difference
  }
*/

gen_op operator - (const gen_op& Op1, const matrix &mx)
  {
  if(!Op1.WBR) return gen_op(mx.operator-()); 	// If Op1 NULL, ret. negated mx
  if(!Op1.OpCheck(mx, 1))			// Check Op1-mx compatibility
    {
    Op1.GenOperror(41, 1); 			// Error during Op-mx operation 
    Op1.GenOpfatality(21);			// Can't do the subtraction
    }
  Op1.set_DBR();				// Put Op1 in default basis
  return gen_op(Op1.WBR->RepMx-mx, Op1.WBR->RepBs);
  }



	// Input		mx   : A matrix
	// 			Op1  : General operator.
	// Return		Op   : Operator difference of input
	//			       matrix and operator
	//			       Op = mx - Op1
	// Note		             : Result EXCLUSIVELY in DBR
	
gen_op operator - (const matrix &mx, const gen_op &Op1)
  {
  if(!Op1.WBR) return gen_op(mx);	// Op1 NULL, Op is matrix
  if(!Op1.OpCheck(mx,1))		// Check Op1-mx compatibility
    {
    Op1.GenOperror(41, 1);		// Error during Op-mx operation 
    Op1.GenOpfatality(21);		// Can't do the subtraction
    }
  Op1.set_DBR();			// Put Op1 in default basis
  return gen_op(mx-Op1.WBR->RepMx, Op1.WBR->RepBs);
  }
	


	// Input		mx   : A matrix
	// 			Op   : General operator (this).
	// Return		Op   : Operator input minus input matrix
	//			       Op = Op - mx
	// Note		             : Result EXCLUSIVELY in DBR

void gen_op::operator -= (const matrix& mx)
  {
  if(!WBR) { *this=mx.operator-(); return; } 	// Op NULL, set to -mx
  if(!OpCheck(mx,1))			// Check dimensions
    {
    GenOperror(41, 1); 			// Error during Op-mx operation 
    GenOpfatality(21);			// Can't do the subtraction
    }
  set_DBR();				// Set Op to DBR
  setOnlyWBR();				// Set Op exclusively to DBR
  WBR->RepMx -= mx;			// Matrix subtraction
  EBR = NULL;				// Remove prior reference to EBR
  check_EBR();				// Check possibility WBR is now EBR
  if(!EBR) WBR->RepPty = DBPr;		// If not set, set DBPr priority
  }



	// Input		Op1  : General operator.
	//			mx   : A matrix
	// Return		Op   : Operator product of the input
	//			       operator and matrix
	//			       Op =  Op1 * mx.
	// Note			     : Result EXCLUSIVELY in DBR
	// Note		             : Order matters - Op1*mx != mx*Op1

gen_op operator * (const gen_op& Op1, const matrix& mx)
  {
  if(!Op1.WBR) return Op1;		// If Op1 NULL, return NULL Op
  if(!mx.cols()) return gen_op();	// If mx NULL, return NULL Op
  if(!Op1.OpCheck(mx,1))		// Check dimensions
    {
    Op1.GenOperror(41, 1); 		// Error during Op-mx operation 
    Op1.GenOpfatality(22);		// Can't do multiplication
    }
  Op1.set_DBR();			// Put Op1 in DBR
  return gen_op(Op1.WBR->RepMx*mx, Op1.WBR->RepBs);
  }
  


	// Input		mx   : A matrix
	//			Op1  : General operator.
	// Return		Op   : Operator product of the input
	//			       matrix and operator
	//			       Op =  mx * Op1
	// Note			     : Result EXCLUSIVELY in DBR
	// Note		             : Order matters - Op1*mx != mx*Op1

gen_op operator * (const matrix& mx, const gen_op& Op1)
  {
  if(!Op1.WBR) return Op1;		// If Op1 NULL, return NULL Op
  if(!mx.cols()) return gen_op();	// If mx NULL, return NULL Op
  if(!Op1.OpCheck(mx,1))		// Check dimensions
    {
    Op1.GenOperror(41, 1); 		// Error during Op-mx operation 
    Op1.GenOpfatality(22);		// Can't do multiplication
    }
  Op1.set_DBR();			// Put Op1 in DBR
  return gen_op(mx*Op1.WBR->RepMx, Op1.WBR->RepBs);
  }
	


	// Input		mx   : A matrix
	// 			Op   : General operator (this).
	// Return		Op   : Operator input times input matrix
	//			       Op = Op * mx
	// Note		             : Result EXCLUSIVELY in DBR

void gen_op::operator *= (const matrix &mx)
  {
  if(!WBR) return;			// If NULL Op, nothing happens
  if(!mx.cols()) { ZeroOp(); return; }	// If NULL mx, zero current Op
  if(!OpCheck(mx,1))			// Check dimensions
    {
    GenOperror(41, 1);			// Error during Op-mx operation 
    GenOpfatality(22);			// Can't do multiplication
    }
  set_DBR();				// Set Op to DBR
  setOnlyWBR ();			// Set Op exclusively to DBR
  WBR->RepMx *= mx;			// Matrix multiplication
  EBR = NULL;				// Remove prior reference to EBR
  check_EBR();				// Check possibility WBR is now EBR
  if(!EBR) WBR->RepPty = DBPr;		// If not set, set DBPr priority
  }
	

	// Input		mx   : A matrix
	// 			Op   : General operator (this).
	// Return		Op   : Operator multiplied by input matrix
	//			       Op = mx * Op
	// Note		             : Result EXCLUSIVELY in DBR

void gen_op::operator &= (const matrix& mx)
  {
  if(!WBR) return;			// If NULL Op, nothing happens
  if(!mx.cols()) { ZeroOp(); return; }	// If NULL mx, zero current Op
  if(!OpCheck(mx,1))			// Check dimensions
    {
    GenOperror(41, 1); 			// Error during Op-mx operation 
    GenOpfatality(22);			// Can't do multiplication
    }
  set_DBR();				// Set Op to DBR
  setOnlyWBR ();			// Set Op exclusively to DBR
  WBR->RepMx *= mx*WBR->RepMx;		// Matrix multiplication
  EBR = NULL;				// Remove prior reference to EBR
  check_EBR();				// Check possibility WBR is now EBR
  if(!EBR) WBR->RepPty = DBPr;		// If not set, set DBPr priority
  }



	// Input		mx   : A matrix
	// 			Op   : General operator (this).
	// Return		Op   : Operator which is a copy of
	//			       the input matrix, Op = mx.
	// Note			     : Op is EXCLUSIVELY in DBR

void gen_op::operator= (const matrix& mx)
  {
  ZeroOp();					// Zero current Op
  if(!mx.rows()) return;			// Done if no matrix
  if(!OpCheck(mx,1))				// Insure mx is square
    {
    GenOperror(41, 1); 				// Error in Op-mx operation 
    GenOpfatality(42);				// Assignment problems
    }
  genoprep OpRep(mx,basis(mx.rows()),DBPr);	// New OpRep we will add
  AddRep(OpRep);				// Add new OpRep as WBR/DBR
  }

// ____________________________________________________________________________
// D                  OPERATOR FUNCTIONS, OPERATOR WITH SCALAR
// ____________________________________________________________________________

/* These function alter an operator by an input scalar value, typically by
   multiplication or division.  The result will be a new operator that is in
   the working basis of the input operator ONLY. 

       Operator    Arguments                           Result 
   --------------- ---------   ------------------------------------------------ 
          *          Op,z      Return Op1 = z*Op             in WBR of Op
          *          z,Op      Return Op1 = z*Op             in WBR of Op
          *          Op,d      Return Op1 = d*Op             in WBR of Op
          *          d,Op      Return Op1 = d*Op             in WBR of Op
          *=         Op,z      Set Op = z*Op                 in WBR of Op
          *=         Op,d      Set Op = d*Op                 in WBR of Op
          /          Op,z      Return Op1 = Op/z             in WBR of Op
          /          Op,d      Return Op1 = Op/d             in WBR of Op
          /=         Op,z      Set Op = Op/z                 in WBR of Op
          /=         Op,d      Set Op = Op/d                 in WBR of Op    */

gen_op operator * (const gen_op&  Op1, const complex& z)   { return z*Op1; }
gen_op operator * (const complex& z,   const gen_op&  Op1)
  {
  if(!Op1.WBR) return Op1;			// Return Null if Op1 is NULL
  return gen_op(Op1.WBR->RepMx*z, Op1.WBR->RepBs);// Multiply WBR matrix only
  }

gen_op operator * (const gen_op& Op1, double d) { return  (complex)d * Op1; }
gen_op operator * (double d, const gen_op& Op1) { return  (complex)d * Op1; }

gen_op & gen_op::operator *= (const complex& z)
{ 
	if(WBR) 
	{ 
		setOnlyWBR(); 
		WBR->RepMx *= z; 
	} 
	return *this;
}

gen_op & gen_op::operator *= (double r)
{ 
	if(WBR) 
	{ 
		setOnlyWBR(); 
		WBR->RepMx *= (complex)r; 
	} 
	return *this;
}

gen_op operator / (const gen_op& Op1, const complex& z)
  {
  if(!Op1.WBR) return Op1;			// Return NULL if Op1 NULL
  return gen_op(Op1.WBR->RepMx/z, Op1.WBR->RepBs);
  }

gen_op operator / (const gen_op &Op1, double r) { return Op1/(complex)r; }

gen_op & gen_op::operator /= (const complex& z)
{ 
  if(WBR) 
  { 
    setOnlyWBR(); 
    WBR->RepMx /= z; 
  } 
  return (*this);
}

gen_op & gen_op::operator /= (double r)
{ 
  if(WBR) 
  { 
    setOnlyWBR(); 
    WBR->RepMx /= (complex)r; 
  } 
  return (*this);
}

                                                                      
// ____________________________________________________________________________
// E                       COMPLEX OPERATOR FUNCTIONS
// ____________________________________________________________________________
 

complex gen_op::det() const { return (WBR)?complex0:(WBR->RepMx).det(); }

	// Input		Op   : General operator (this)
        // Return		z    : Complex number, determinant of Op
	// Note			     : Performed in WBR of Op


complex gen_op::trace() const { return (!WBR)?complex0:(WBR->RepMx).trace(); }

	// Input		Op   : General operator (this)
        // Return		z    : Complex number, trace of Op
	// Note			     : Performed in WBR of Op


complex gen_op::trace(const gen_op &Op1) const

	// Input		Op   : General operator (this)
	// 			Op1  : General operator
        // Return		z    : Complex, trace of (Op*Op1)
	// Note			     : Performed in WBR of Op

  {
  if(!WBR || !Op1.WBR) return complex0;		// Zero if Op or Op1 NULL
  if(!OpCheck(Op1,1)) GenOpfatality(3,"trace"); // Check Op-Op1 compatibility
  const genoprep* REP;				// An operator representation
  REP = Op1.FindRep(WBR->RepBs);		// Seek Op1 rep. in Op's WBR
  if(REP) return (WBR->RepMx).trace(REP->RepMx);// If exists, output Tr{Op*Op1}
  REP = FindRep(Op1.WBR->RepBs);		// Seep Op rep. in Op1's WBR
  if(REP) return (Op1.WBR->RepMx).trace(REP->RepMx);// If exists, output Tr{Op*Op1}
/*
if(REP)
{
std::cout << "\n\tFound Matching Op-Op for Trace....";
status(1);
Op1.status(1);
std::cout << "\n\t\tTr(Op1*Op2) = " << (WBR->RepMx).trace(REP->RepMx);
std::cout << "\n\t\tTr(Op2*Op1) = " << (REP->RepMx).trace(WBR->RepMx);
// sosi - we should now look for any matching representation between Op & Op1
//        before doing any calculations......
}
std::cout << "\n\tMust Generate New Base For Trace....";
std::cout.flush();
status(1);
Op1.status(1);
*/
  set_DBR();					// Put Op into its DBR
  Op1.set_DBR();				// Put Op1 into its DBR
  return (DBR->RepMx).trace(Op1.DBR->RepMx);	// Now we take the trace...
  }


complex gen_op::proj(const gen_op &Op2, int normf) const

	// Input		Op1  : General operator (this)
	// 			Op2  : General operator
	// 			normf: Value of normalization
	//			       if zero (default) it is computed
        // Return		z    : Complex, projection of Op1 on Op2
	//				      t		    t
	//				tr(Op2 * Op1)/tr(Op2  * Op2)

  {
  if(!WBR || !Op2.WBR) return complex0;	// Zero if Op1 or Op2 NULL
  if(!OpCheck(Op2,1))			// Check Op1 - Op2 compatibility
    GenOpfatality(3, "proj"); 		// Op projection error
  complex znum,zden=1;
  Op_base(Op2);				// Set Op1 to basis of Op2
  gen_op Op2adj = adjoint(Op2);		// Store Op2 adjoint
  znum=(Op2adj.WBR->RepMx).trace(WBR->RepMx); // Get <Op2|Op1>
  if(normf == 0)			// Compute denominator if not given
    zden = (Op2adj.WBR->RepMx).trace(Op2.WBR->RepMx);
  if(norm(zden) < 1.e-15)		// Check for denominator zero
    GenOpfatality(3, "proj"); 		// Op projection error
  return (znum/zden);			// All O.K., return projection
  }


	// Input		Op   : General operator (this)
        // Return		int  : Op Liouville space dimension

int gen_op::dim()    const { return (!WBR)?0:WBR->RepMx.cols(); }
int gen_op::HS()     const { return (!WBR)?0:WBR->RepMx.cols(); }
int gen_op::LS()     const { return (!WBR)?0:WBR->RepBs.dim_LS(); }
int gen_op::dim_LS() const { return (!WBR)?0:WBR->RepBs.dim_LS(); }



	// Input		Op1   : General operator (this)
        // Return		Op    : exponential of Op1
	//				Op = exp(Op1)
        // Note			      : Computed in EBR of Op1
        // Note			      : Returns I matrix for Null Op

gen_op gen_op::exp() const
  {
  int hs = dim();				// Operator Hilbert space
  if(!WBR)					// If we have no reps then
    {						// return an identity matrix
    if(!hs) GenOpfatality(3, "exp"); 		// or Op exponential error
    return gen_op(matrix(hs,hs,i_matrix_type)); // if there is no dimension
    }
  set_EBR();					// Put Op in EBR
  gen_op ExpOp(*this); 				// Copy Op in diagonal form
  complex z;
  for(int i=0; i<hs; i++) 			// Exponentiate all elements
    {  						// (of which are only diags)
    z = (WBR->RepMx).get(i,i);
    (ExpOp.WBR->RepMx).put(z.Zexp(),i,i);  
    }  
  return ExpOp;
  }     



        // Input                Op    : Operator (this)
        //                      t     : Exponential factor
        //                      cutoff: Exponential factor roundoff
        // Return               ExpOp : Exponential of Op
        //                              ExpOp = exp(t*Op)
        // Note                       : Exponential output in EBR of Op
        // Note                       : Op's EBR is generated herein
        // Note                       : Value of t is considered 0 if
        //                              it's magnituded is less than cutoff

gen_op gen_op::exp(const complex& t, double cutoff) const
  {
  int hs = dim();				// Operator Hilbert space
  if(!WBR && !hs) GenOpfatality(3,"exp"); 	// Op exponential error
  if(!WBR || norm(t) < fabs(cutoff))
    return gen_op(matrix(hs, hs, i_matrix_type));
  set_EBR();					// Put Op in EBR
  gen_op ExpOp(*this); 
  complex z;
  for(int i=0; i<hs; i++)
    {  
    z = t * ((WBR->RepMx).get(i,i));
    (ExpOp.WBR->RepMx).put(z.Zexp(),i,i);  
    }  
  return ExpOp;
  }
 
	// Input		Op   : General operator (this)
        // Return		int  : Op Liouville space dimension


	// Input		Op1   : General operator (this)
        // Return		Op    : exponential of Op1
	//				Op = exp(Op1)
        // Note			      : Computed in same base as Op1
        // Note			      : Returns I matrix for Null Op

gen_op gen_op::expm() const
  {
  int hs = dim();				// Operator Hilbert space
  if(!WBR)					// If we have no reps then
    {						// return an identity matrix
    if(!hs) GenOpfatality(3, "exp"); 		// or Op exponential error
    return gen_op(matrix(hs,hs,i_matrix_type)); // if there is no dimension
    }
  gen_op ExpOp(*this); 				// Copy Op in diagonal form
  matrix mx1 = this->get_mx();
  ExpOp.put_mx(mx1.expm());
  return ExpOp;
  }     



        // Input                Op    : Operator (this)
        //                      t     : Exponential factor
        //                      cutoff: Exponential factor roundoff
        // Return               ExpOp : Exponential of Op
        //                              ExpOp = exp(t*Op)
        // Note                       : Exponential output in same base as Op
        // Note                       : Value of t is considered 0 if
        //                              it's magnituded is less than cutoff

gen_op gen_op::expm(const complex& t, double cutoff) const
  {
  int hs = dim();				// Operator Hilbert space
  if(!WBR && !hs) GenOpfatality(3,"exp"); 	// Op exponential error
  if(!WBR || norm(t) < fabs(cutoff))
    return gen_op(matrix(hs, hs, i_matrix_type));
  gen_op ExpOp(*this); 
  matrix mx1 = this->get_mx()*t;
  ExpOp.put_mx(mx1.expm());
  return ExpOp;
  }
 

        // Input                Op      : Operator (this)
        //                      power   : Exponential power
        // Return               Op^n    : Op taken to the nth power
        // Note                         : Output in EBR of Op
        // Note                         : Op's EBR is generated herein
 
gen_op gen_op::Pow(int power) const
  {
  int hs = dim();				// Operator Hilbert space
  if(!hs && !WBR) GenOpfatality(3,"Pow"); 	// Op power error
  complex z, zpow(power);  			// Working complex numbers
  set_EBR();					// Put Op in EBR
  gen_op PowOp(*this);				// Copy the (diagonal) operator
  for(int i=0; i<hs; i++)			// Loop over the elements
    {						// and take the power
    z = PowOp.WBR->RepMx(i,i);
    PowOp.WBR->RepMx.put(pow(z,zpow), i, i);
    }
  return PowOp;
  }  

gen_op tensor_product(const gen_op &Op1, const gen_op &Op2)

        // Input                Op1  : General operator.
        //                      Op2  : General operator.
        // Return               Op   : Operator tensor product of the two input
        //                             operators, Op =  Op1 X Op2.
        // Note                      : Order matters - Op1XOp2 != Op2XOp1
	// Note			     : No basis checking is done here
	//			       nor any dimension checking
	// sosi			     - What is this used for? Floquet?

  {
  if(!Op2.WBR || !Op1.WBR) return gen_op();	// Op1 or Op2 NULL, return NULL
    return gen_op (tensor_product(Op1.WBR->RepMx,Op2.WBR->RepMx),
                   tensor_product(Op1.WBR->RepBs,Op2.WBR->RepBs) );
  }  



        // Input                Op1  : General operator
        // Return               Op   : logarithms of Op1
        //                             Op = log(Op1)
        // Note                      : Computed in EBR of Op1

gen_op log(const gen_op &Op1)
  {
  if(!Op1.WBR) Op1.GenOpfatality(3,"log"); 	// Op log error if NULL
  Op1.set_EBR();				// Put Op in EBR
  gen_op LogOp(Op1);
  for(int i=0; i < Op1.WBR->RepMx.rows(); i++)
    LogOp.WBR->RepMx.put(log(Op1.WBR->RepMx(i,i)), i, i);
  return LogOp;
  }
  

	// Input		Op1  : General operator (this)
	// 			Op2  : General operator
        // Return		Op   : General operator from similarity
	//			       transform of Op2 by Op1 
	// 			       Op = Op1 * Op2 * [Op1]
        // Note			     : Op EXCLUSIVELY in WBR of Op1
        // Note			     : Op2 switched to WBR of Op1

gen_op gen_op::sim_trans(const gen_op& Op2) const
  {
  if(!WBR || !Op2.WBR) return gen_op();	// Result NULL if Op1 or Op2 NULL
  if(!OpCheck(Op2,1))			// Check Op1-Op2 dimensions
    GenOpfatality(3,"sim_trans");	// Op sim_trans error
  Op2.Op_base((*this));			// Put Op2 in WBR of Op1
  return gen_op(WBR->RepMx			// Op = Op1*Op2*(Op1)t
          *times_adjoint(Op2.WBR->RepMx,WBR->RepMx),WBR->RepBs);
  }


void gen_op::sim_trans_ip(const gen_op &Op1)

	// Input		Op   : General operator (this)
	// 			Op1  : General operator
        // Return		None : Similarity transform of Op by Op1 
	//					            t
	// 			       Op = Op1 * Op * [Op1]
        // Note			     : Op EXCLUSIVELY in WBR of Op1

  {
  if(!WBR) return;			// Result NULL if Op NULL
  if(!Op1.WBR) *this = gen_op();	// Result NULL if Op1 NULL
  if(!OpCheck(Op1,1))			// Check dimensions
    GenOpfatality(3,"sim_trans_ip");	// Op sim_trans error
  Op1.Op_base(*this);			// Put Op1 into WBR of Op
  setOnlyWBR ();			// Delete all reps by WBR of Op
  WBR->RepMx = Op1.WBR->RepMx			// Op = Op1*Op*(Op1)t
	    * times_adjoint(WBR->RepMx,Op1.WBR->RepMx);
  }


//gen_op gen_op::adjoint()

	// Input		Op1   : General operator (this)
        // Return		Op    : Hermitian adjoint of Op1
	//                          	          *t
	// 			        Op = [Op1]  (transpose & conjugate)
        // Note                       : EXCLUSIVELY in the WBR of Op

//{
//  if (WBR)			// Check for NULL Op1
				// Construct during return
//    return gen_op(adjoint(WBR->RepMx), WBR->RepBs);
//  else
//    return gen_op(*this);		// Op1 NULL, return NULL
//}


row_vector gen_op::eigvals() const

	// Input		Op   : General operator (this)
        // Return		vx   : Vector of eigenvalues

  {
  int hs=dim();				// Get Op dimension
  row_vector evals(hs);			// Array for return
  set_EBR();				// Set Op into eigenbasis
  for(int i=0; i<hs; i++)		// Fill vector with Op diag
    evals.put(get(i,i),i);
  return evals;
  }
  

void gen_op::eigvals(double* vx, int sort) const

	// Input		Op   : General operator (this)
	//			vx   : Double vector
	//			sort : Flag for eigenvalue sorting 
	//			        0 - default
	//				1 - sort on reals
	//			       -1 - sort on imaginaries
	//				# - sort on norms
        // Return		void : vx filled with eigenvalues

  {
sort=0;
//  double max, maxe, vxj=0;
//  int maxi;
  set_EBR();
  for(int i=0; i<dim(); i++)
    vx[i] = Re(get(i,i));
//  if(sort)			// Sort the eigenvalues if needed
//    {
//    for(i=0; i<dim()-1; i++)
//      {
//      if(sort == 1)
//        max = Re(vx[i]);
//      else if(sort == -1)
//        max = Im(vx[i]);
//      else
//        max = norm(vx[i]);
//      maxe = vx[i];
//      maxi = i;
//      for(int j=i+1; j<dim(); j++)
//        {
//        if(sort == 1)
//          vxj = Re(vx[j]);
//        else if(sort == -1)
//          max = Im(vx[j]);
//        else
//          max = norm(vx[j]);
//        if(vxj > max)
//          {
//          maxe = vx[j];
//          maxi = j;
//          if(sort == 1)
//            max = Re(vx[i]);
//          else if(sort == -1)
//            max = Im(vx[i]);
//          else
//            max = norm(vx[i]);
//          }
//        }
//      if(maxi != i)
//        {
//        vx[maxi] = vx[i];
//        vx[i] = maxe;
//        }
//      }
//    }
  return;
  }
  

// These functions are now marked for deletion since they have been
// replaced by their member function counterparts!
// Sosi 2/21/92

complex det(const gen_op& Op)                        { return Op.det(); }
complex trace(const gen_op& Op)                      { return Op.trace(); }
complex trace(const gen_op& Op, const gen_op &Op1)   { return Op.trace(Op1); }
complex proj(const gen_op &Op,  const gen_op &Op1, int N)
                                                     { return Op.proj(Op1,N); }
int dim(const gen_op &Op)                            { return Op.dim(); }
gen_op exp(const gen_op &Op)                         { return Op.exp(); }
gen_op pow(const gen_op &Op, int power)              { return Op.Pow(power); }
gen_op sim_trans(const gen_op& Op,const gen_op& X)   { return Op.sim_trans(X); }

gen_op adjoint(const gen_op &Op1)

	// Input		Op1   : General operator
        // Return		Op    : Hermitian adjoint of Op1
	//                          	          *t
	// 			        Op = [Op1]  (transpose & conjugate)
        // Note                       : EXCLUSIVELY in the WBR of Op

  {
  if(!Op1.WBR) return Op1;		// Just return for NULL Op1
  return gen_op(adjoint(Op1.WBR->RepMx), Op1.WBR->RepBs);
  }

// ____________________________________________________________________________
// F                      OPERATOR COMPONENT MANIPULATIONS
// ____________________________________________________________________________

/*     Function    Arguments                           Result 
   --------------- ---------   ------------------------------------------------ 
       get_mx         ---      Returns the matrix of current OpRep (WBR)
     get_matrix       ---      Returns the matrix of current OpRep (WBR)
       put_mx          mx      Sets WBR matrix to be mx, other OpReps deleted
     put_matrix        mx      Sets WBR matrix to be mx, other OpReps deleted
       get_bs         ---      Returns the basis of current OpRep (WBR)
     get_basis        ---      Returns the basis of current OpRep (WBR)
       put_bs          bs      Sets WBR basis to be bs, other OpReps deleted
     put_basis         bs      Sets WBR basis to be bs, other OpReps deleted
       (i,j)         int,int   Returns <i|Op|j> in current OpRep (WBR) 
        get          int, int  Returns <i|Op|j> in current OpRep (WBR) 
        put          int, int  Sets <i|Op|j> in WBR, others deleted

  Note that Op(i,j) returns a copy of the operator element, NOT the element.
  Thus, code such as "Op(i,j) = z;" will not work!                           */

matrix gen_op::get_mx()     const { return get_matrix(); }
matrix gen_op::get_matrix() const { return (!WBR)?matrix(0,0):WBR->RepMx; }

void gen_op::put_mx(const matrix &mx) { put_matrix(mx); }
void gen_op::put_matrix(const matrix &mx)
  {
  if(!WBR) *this = gen_op(mx);		// If NULL, Op is mx (= checks mx)
  setOnlyWBR();				// Delete all reps but WBR
  if(!OpCheck(mx, WBR->RepBs, 1))	// Insure acceptable mx,bs pair
    {
    GenOperror(6,1);			// Trouble accessing Op internals
    GenOpfatality(39);			// Bad mx,bs mixture
    }
  WBR->RepMx = mx;			// Replace Op WBR matrix
  EBR = NULL;				// Set EBR to be NULL
  check_EBR();				// Try and find an EBR too
  }

// -------------------------- Basis Manipulations -----------------------------

basis gen_op::get_bs()    const { return get_basis(); }
basis gen_op::get_basis() const { return (!WBR)?basis(1):WBR->RepBs; }

void gen_op::put_bs(const basis &bs) { put_basis(bs); }
void gen_op::put_basis(const basis &bs)
  {
  if(!WBR)				// Check for NULL Op
    {
    *this = gen_op();			// Return NULL Op
    GenOperror(6,1);			// Trouble accessing Op internals
    GenOperror(7,1);			// Accessing empty Op
    GenOperror(77);			// Error setting Op basis
    }
  setOnlyWBR();				// Delete all reps but WBR
  WBR->RepBs = bs;			// Replace Op WBR basis
  check_DBR();				// Perhaps this is DBR?
  }

// -------------------- Individual Element Manipulations ----------------------

void gen_op::put(const complex &z, int row, int col)
  {
  if(!WBR)				// Check for NULL Op
    {
    GenOperror(7,1);			// Dealing with NULL Op
    GenOpfatality(67);			// Cannot access element
    }
  if(!OpCheck(row,col))			// Check element indices
    {
    GenOperror(6,1);			// Trouble accessing Op internals
    GenOpfatality(67);			// Cannot access element
    }
  setOnlyWBR();				// Delete any other representations
  if(row!=col) EBR = NULL;		// Cannot be EBR if off-diagonal
  WBR->RepMx.put(z,row,col);		// Set the element value
  }

complex gen_op::operator() (int i, int j) const { return get(i,j); }
complex gen_op::get(int row, int col)     const
  {
  if(!WBR) return complex0;		// Nothing if NULL Op
  if(!OpCheck(row,col))			// Check element indices
    {
    GenOperror(6,1);			// Trouble accessing Op internals
    GenOpfatality(67);			// Cannot access element
    }
  return (WBR->RepMx.get(row,col));	// Output the element value
  }

// ____________________________________________________________________________
// G                     OPERATOR AUXILIARY FUNCTIONS
// ____________________________________________________________________________

std::string gen_op::name() const               { return OpName; }
void        gen_op::name(const std::string& n) { OpName = n; }

void gen_op::bsname(const std::string& bn)
  { if(WBR) (WBR->RepBs).name(bn); }

int gen_op::exists() const { return WBR?1:0; }

	// Input		Op    : General operator (this)
	// Output		int   : TRUE if operator not NULL



        // Input                Op     : General operator (this)
        // Output               vect   : Vector of aligned Op elements
	// Note 		       : Added in support of multi-sys

col_vector gen_op::superket() const
  {
  basis bs = (*this).get_basis();
  matrix mx = (*this).get_mx();
  int dimS = bs.dim_LS();
  col_vector mS(dimS);
  int mSptr=0, pinblock, j, k;
  complex elt;
  for(int i=0; i<bs.sub_N(); i++)
    {
    pinblock = bs.sub_anchor(i);
    for(j=pinblock; j<pinblock+bs.sub_dim(i); j++)
      for(k=pinblock; k<pinblock+bs.sub_dim(i); k++)
        {
        elt = mx.get(j,k);
        mS.put(elt, mSptr);
	mSptr++;
	}	
     }
  return mS;
  }



	// Input		: Operator (this) with preserved basis
	//			  structure
	//			: Superket that had been formed from Op
	// Output		: General Operator (this) with updated mx
	// Note 		       : Added in support of multi-sys

void gen_op::desuperket(const col_vector& mS)
  {
  basis bs=(*this).get_basis();
  matrix mX=(*this).get_mx();
  int mSptr=0, pinblock, i, j, k;
  complex elt;
  for (i=0; i<bs.sub_N(); i++)
    {
      pinblock = bs.sub_anchor(i);
      for (j=pinblock; j<pinblock+bs.sub_dim(i); j++)
        for (k=pinblock; k<pinblock+bs.sub_dim(i); k++)
          {
            elt = mS.get(mSptr);
            mX.put(elt,j,k);
            mSptr++;
          }
     }
  setOnlyWBR ();
  (*this).put_mx(mX);
  }

gen_op gen_op::project_sub(int ic) const

	// Input		Op  : General operator (this)
	//			int : Multi_sys component ic
	// Output		    : Operator's ic component in
	//			      the subspace of ic
	// Note 		    : Added in support of multi-sys

  {
  basis bS = get_basis();			// Get basis of Op
  int nc = bS.sub_N();				// Number of components
  if((ic<0) || (ic>=nc)) return gen_op();	// Null Op if bad component
  if(nc == 1) return gen_op(*this);		// If 1 component, return Op
  int ssi = bS.sub_anchor(ic);			// Index where subspace ic starts
  int n = bS.sub_dim(ic);			// Dimension of subspace ic
  matrix mX = get_mx();				// Get matrix of Op
  matrix bSx = bS.U();
//  matrix bSx = (const matrix&)bS;		// Get basis as a matrix
  matrix mxic = mX.get_block(ssi,ssi,n,n);	// Clip out Op mx this component
  matrix bsx = bSx.get_block(ssi,ssi,n,n);	// Clip out Op basis this coponent
  basis bsic(bsx);				// Form basis this component
  return gen_op(mxic, bsic);			// Return operator for this component
  }  

// ____________________________________________________________________________
// H                  OPERATOR REPRESENTATION MANIPULATIONS
// ____________________________________________________________________________

/* These functions are quick checks to see of a particular reperesentation
   exists and/or whether the current represenation is the working OpRep.  The
   functions check for EBR & DBR, an eigenbasis and the default basis OpReps.*/

int gen_op::test_EBR() const { return (EBR != NULL); }
int gen_op::test_DBR() const { return (DBR != NULL); }
int gen_op::in_EBR()   const { return (EBR == WBR && WBR); }
int gen_op::in_DBR()   const { return (DBR == WBR && WBR); }


	// Input		Op    : General operator (this)
	// Output		void  : Outputs operator status

void gen_op::status(int pf) const
  {
  int nreps = size();
  std::cout << "\n\tCurrent Number of Representations: " << nreps;
  int spaces = 0;
  std::vector<genoprep>::const_iterator item;// Iterator into list
  item = begin(); 			// The 1st OpRep in list
  const genoprep* REP = ObeginC();
  for(int i=0; i<nreps; i++)
    {
    std::cout << "\n\tRep " << i+1 << "- priority = " << REP->RepPty;
    spaces++;
    if(REP->RepPty > 9)  spaces++;
    if(REP->RepPty > 99) spaces++;
    for(int j=spaces; j<3; j++) std::cout << " "; 
    if(REP == WBR) std::cout << " (WBR)";
    if(REP == DBR) std::cout << " (DBR)";
    if(REP == EBR) std::cout << " (EBR)";
    std::cout << "\n\t     - matrix = " << REP->RepMx.refs() << " refs";
    switch((REP->RepMx).stored_type())
      {
      case i_matrix_type: std::cout << ", Identity Matrix";  break;
      case d_matrix_type: std::cout << ", Diagonal Matrix";  break;
      case h_matrix_type: std::cout << ", Hermitian Matrix"; break;
      case n_matrix_type: std::cout << ", Complex Matrix";   break;
      default:            std::cout << ", Generic Matrix";   break;
      }
    switch((REP->RepMx).stored_type())
      {
      case h_matrix_type: 				             break;
      case i_matrix_type:
      case d_matrix_type:
      case n_matrix_type:
        if((REP->RepMx).test_hermitian())
	  std::cout << ", Hermitian"; 
        else                              
	  std::cout << ", Non-Hermitian";
	break;
      default:                            
	std::cout << ", Non-Hermitian";
	break;
      }

    std::cout << "\n\t     - basis  = " << (REP->RepBs).refs() << " refs";
//    switch((((const matrix&)REP->RepBs)).stored_type())
    switch(((REP->RepBs).U()).stored_type())      {
      case i_matrix_type: std::cout << ", Identity Matrix";  break;
      case d_matrix_type: std::cout << ", Diagonal Matrix";  break;
      case h_matrix_type: std::cout << ", Hermitian Matrix"; break;
      case n_matrix_type: std::cout << ", Complex Matrix";   break;
      default:            std::cout << ", Generic Matrix";   break;
      }
//    switch(((const matrix&)(REP->RepBs)).stored_type())
    switch(((REP->RepBs).U()).stored_type())      {
      case h_matrix_type: 				             break;
      case i_matrix_type:
      case d_matrix_type:
      case n_matrix_type:
        if((REP->RepMx).test_hermitian())
	  std::cout << ", Hermitian";
        else                             
	  std::cout << ", Non-Hermitian";
	break;
      default:                           
	std::cout << ", Non-Hermitian";
	break;
      }
    REP++;
    item++;
    }
  if(pf)
    {
    REP = ObeginC();
    std::string rn;
    for(int k=0; k<nreps; k++)
      {
      rn = "Representation " + Gdec(k);
      std::cout << "\n\n" << CenterString(rn) << "\n";
      REP->print(std::cout);
      REP++;
      }
    }
  std::cout << "\n";
  }

// *************  Set Operator to a Specific Representation *************

/* These functions are will put an operator into a particular basis.  A new
   OpRep will be formed if the representation doesn't already exist.  The
   Op_base functions set "this" into the basis/Op.WBR input as an argument.
   A cutoff can be specified that will "check" the validity of an EBR.       */
	

	// Input		Op    : General operator (this)
	// Output		none  : Op WBR set to DBR 
	// Note			      : If no DBR exists for Op then
	//				one is produced with NEW
	//				default basis produced.

void gen_op::set_DBR() const
  {
  if(!WBR) return;				// Nothing if NULL Op
  if(WBR == DBR) return;			// Nothing if WBR already DBR
  if(DBR) { WBR=DBR; return; }			// If DBR exists, set WBR to DBR
  matrix mx=WBR->RepBs.convert_back(WBR->RepMx);// Convert WBR mx to DBR mx
  genoprep OpRep(mx,basis(mx.rows()),DBPr);	// New Op representation (DBR)
  AddRepM(OpRep);				// Add new representation
  SetLimits(MaxReps);				// Impose rep limit on Op
  }

	// Input		Op   : General operator (this)
	// Output		none : Op WBR set to EBR 
	// Note			     : Will generate an EBR if none exists
	// Note			     : Assumes the MaxReps limit always
	//			       always allows for EBR

void gen_op::set_EBR() const
  {
  if(!WBR) return;				// Do nothing if NULL Op
  if(WBR == EBR) return;			// Do nothing if Op WBR is EBR
  if(EBR != NULL) { WBR=EBR; return; }		// EBR exists, set to WBR & Done
  set_DBR();					// First place Op in DBR
  if(EBR != NULL) { WBR=EBR; return; }		// Remote Chance DBR made EBR?
  basis bs0 = DBR->RepBs;			// This is the default basis
  matrix mxd, mxev;				// Arrays for diagonalization
  diag(DBR->RepMx, mxd, mxev);			// Diagonalize Op in DBR
  basis bs(bs0, mxev);				// Construct eigenbasis
  genoprep OpRep(mxd,bs,EBPr);			// EBR Op representaton
  AddRepM(OpRep);				// Add EBR to Op reps
  SetLimits(MaxReps);				// Impose rep limit on Op
  }

	// Input		Op     : General operator(this)
	// 			Op1    : General operator
	//			cutoff : Cutoff for checking EBR
	//				 Default values is 1.e-12
        // Return		none   : Op put into basis of Op1
	// 			         Op = Op1.bs * Op * [Op1.bs]
	// Note			       : It is assumed that Op & Op1
	//			         share the same default basis	
	//				 & indeed this is enforced to conserve
	//				 potential future basis checking
	// Note 		       : This is one of the few 
	//                               places where we must insure there
	//                               are not too many representations
	//                               (call to function SetLimits)

void gen_op::Op_base(const gen_op &Op1, double cutoff) const
  {
  if(WBR == NULL) return;	 		// A NULL Op is in any base
  if(Op1.WBR == NULL)				// Check for NULL Op1
    {
    GenOperror(7,1);				// Dealing with NULL Op
    GenOpfatality(77);				// Problems setting basis
    }
  if(!OpCheck(Op1,1)) GenOpfatality(77);	// Problems setting basis
  if(WBR->RepBs == Op1.WBR->RepBs) return;	// Exit if Op in Op1 basis
  if(SetRep(Op1.WBR->RepBs)) return;		// Have OpRep in Op1 WBR?
# ifdef OpDebug
  std::cout << "\n\tSetting Basis Of ";
  if(OpName.length()) std::cout << OpName << " ";
  std::cout << "Operator";
  if(Op1.name().length()) std::cout << " To That Of "
                      << Op1.name() << " Operator";
# endif
  set_DBR();					// Put Op into its DBR
  if(Op1.DBR) WBR->RepBs = Op1.DBR->RepBs; 	// Insure Op-Op1 DBRs bases same
  if(Op1.WBR != Op1.DBR)			// If DBR isn't OpRep wanted
    {						// (most likely it isn't) 
    matrix mx=Op1.WBR->RepBs.convert(WBR->RepMx);// Op array into basis of Op1
    basis bs = Op1.WBR->RepBs;			// Basis will be Op1 basis
    AddRepM(genoprep(mx,bs,1));			// Add OpRep, default priority
    check_EBR(cutoff);				// See if rep is EBR
    SetLimits(MaxReps);				// Impose rep limit on Op
// sosi
// xxxxx this is new, check it works & look at EBR priority?
    }
  }



	// Input		Op   : General operator(this)
	// 			bs   : Basis
        // Return		none : Op put into basis of Op1
	// 			       Op = bs * Op * [bs]t

void gen_op::Op_base(const basis& bs) const
  {
  if(WBR == NULL) return; 			// Exit if NULL Op
  if(!OpCheck(bs,1)) GenOpfatality(77);		// Problems setting basis
  if (WBR->RepBs == bs) return; 		// If Op basis is bs were done
  if(SetRep(bs)) return;			// Have OpRep in basis bs?
  set_DBR();					// Set Op into its DBR
  if (DBR->RepBs == bs) return; 		// Maybe DBR is basis we want?
  matrix mx = bs.convert(DBR->RepMx);		// Set new rep to basis
  AddRepM(genoprep(mx,bs,1));			// Add OpRep, default priority
  SetLimits(MaxReps);				// Impose rep limit on Op2
  }



	// Input		Op   : General operator (this)
	// Output		none : Deletes all reps but WBR of Op
	// Note			     : Does nothing if Op is NULL
	// Note			     : This routine seems to be too complex.
	//			       Since use of "erase" seems to move
	//			       around OpRep pointers in the list,
	//			       we have to track WBR explicitly.....!?
	//			       Perhaps this is a compiler bug in STL?

void gen_op::setOnlyWBR ()
  {
  if(!WBR) return;			// If NULL Op do nothing
  if(size() == 1) return;		// Exit if only WBR
  int EKeep=0, DKeep=0;			// Flags to maintain EBR,DBR
  if(EBR == WBR) EKeep=1; 		// We must save EBR=WBR
  if(DBR == WBR) DKeep=1; 		// We must save DBR=WBR
  int WBRpos=0, nOpReps=size();		// WBR position, # OpReps
  genoprep* REP = Obegin(); 		// The 1st OpRep in list
  for(int i=0; i<nOpReps; i++, REP++)	// Loop reps and find WBR
    if(REP == WBR) WBRpos = i;		// position OpRep list
  std::vector<genoprep>::iterator item;	// Iterator into list
  item = begin(); 			// The 1st OpRep in list
  erase(item, item+WBRpos);		// Remove all before WBR
  if(size() > 1) 			// Remove any OpReps after WBR
    {
    item = begin()+1; 			// The 2nd OpRep in list
    erase(item, item+size()-1);		// Remove all but 1st one
    }
  WBR = Obegin();			// Reset WBR to 1st OpRep
  (EKeep)?EBR=WBR:EBR=NULL; 		// Set EBR to WBR as needed
  (DKeep)?DBR=WBR:DBR=NULL;		// Set DBR to WBR as needed
  }


// ******************* Representation Limits and Priorities *******************

void gen_op::Op_priority(int pty) { WBR->RepPty = pty; }

	// Input		Op   : General operator (this).
	// 			pty  : Priority value
	// Output		none : Assigns Op WBR priority
	// Note			     : Higher pty value implies rep
	//			       preferentially maintained



	// Input		Op   : General operator (this).
	// 			limit: Representation limit number 
	// Output		none : Insures number of Op reps
	//			       does not exceed limit
	// Note                      : This very simple function is quite
	//                             critical. It keeps the number of
	//                             operator representations from
	//			       becoming infinitely large!

void gen_op::SetLimits(int limit) const
  {
  LimCheck(limit);				// Insure value reasonable
  if(!WBR) return; 				// Don't do anything if NULL Op
  if(int(size()) <= limit) return;		// Do nothing limit not exceeded
# ifdef OpDebug
  std::cout << "\n\tLimiting No. Of ";
  if(OpName.length()) std::cout << OpName << " ";
  std::cout << "Operator Representations";
  std::cout << "\n\tCurrent Representations: " << size();
  std::cout << ", Maximum Allowed: " << limit; 
# endif
  gen_op* OpTmp = const_cast<gen_op*>(this);	// Cast away const (mutable!)
  genoprep* REP = NULL;				// Pointer to OpRep
  genoprep* LP  = NULL;				// Another pointer to OpRep
  int PtyMin;					// For OpRep priority checking
  std::vector<genoprep>::iterator item;		// Iterator into list
  std::vector<genoprep>::iterator lpi;
  while(int(size()) > limit)			// Since # OpReps>limit we have
    {						// to delete some of them.
    REP  = (*OpTmp).Obegin();			//   Get at 1st OpRep (genoprep*) 
    item = (*OpTmp).begin();			//   Same but vector iterator
    if(REP==WBR) { REP++; item++; } 		//   Must insure that NOT WBR
    PtyMin = REP->RepPty;			//   Set initial priority value
    LP = REP; 					//   Set intiial OpRep pointer 
    lpi = item;
    while(item != (*OpTmp).end())		//   Loop rest of OpRep list
      {						//   & find low priority OpRep
      if(REP->RepPty<PtyMin && REP!=WBR)	//     Reset low priority rep
        { 					//     if necessary
        PtyMin = REP->RepPty;			//	 Lowest Oprep priority
        LP = REP;				//	 Lowest priority OpRep
        lpi = item;				//       Pointer to low one
        }
      REP++;					//   Go to next OpRep in list
      item++;
      }
    (*OpTmp).erase(lpi);			// Delete lowest priority OpRep
    }						// See if we must delete another
// sosi - the use of erase has now been found to change
//        the OpRep pointers. Does this mess up because of it?  Might
//        have to do some bookkeeping explicitly for this to work then....
  }

// ____________________________________________________________________________
// I                  CLASS OPERATOR EQUALITY & INEQUALITY
// ____________________________________________________________________________

/* These functions check whether two operators are equivalent.  They are taken
   to be equal if their default basis matrix representations are the same. If
   class gen_op is protecting its data properly, there is no way to generate
   two operators having the same DBR and not being equal.  Note that this then
   assumes that the number of representations may be unequal even though the
   operators are taken as equal.

           Input                Op	: First operator (this)
                                Op1	: Second operator
           Output               TF	: True if Op = Op1 (==)
				          True if Op != Op2 (!=)            */
 
// thse are now in section M
/*
int gen_op::operator!= (const gen_op& Op1) const { return !(*this==Op1); }
int gen_op::operator== (const gen_op& Op1) const
  {
//std::cout << "\n\tOperator used to call == " << *this;
  if(!WBR) return (!Op1.WBR);		// Return TRUE if both NULL
//std::cout << "\n\tOperator compared with in == " << Op1;
  if(!Op1.WBR) return 0;		// Return FALSE if Op1 NULL & not Op
//std::cout << "\nThis is the gen_op::== operator, do I reach here?";
  const genoprep* REP = NULL;		// Don't know if any OpRep's match
  const genoprep* PTR = Op1.begin();	// Begin with 1st OpRep of Op1
  while(PTR!=end())			// Loop through OpReps of Op1
    {
    REP = FindRep(PTR->RepBs);		// Look for OpRep in Op with same basis
    if(REP)				// If we have found match, then the Ops
//{
//if(PTR->RepMx == REP->RepMx)
//std::cout << "\n\tTWO Ops FOUND TO BE EQUAL";
      return(PTR->RepMx == REP->RepMx); // are s equal if their matrices match
//}
    PTR++;
    }
  set_DBR();				// No matching OpReps, so as a last 
  Op1.set_DBR();			// resort we'll see if DBR's match
  return (DBR->RepMx == Op1.DBR->RepMx);// Check DBR representations
// sosiw need to insure Op & Op1 return in original representation....
  }
*/

// ____________________________________________________________________________
// J                        CLASS OPERATOR CHECKS
// ____________________________________________________________________________

/* This set of functions insures that the Op/mx/bs is properly square and/or
   of compatible dimension to mix with another Op/mx/bs.  All return TRUE if
   dimensions look OK.  The functions take a warning level flag which will
   be used to decide what happens under a FAIL.  If the flag is 0 then nothing
   will result.  If the flags is 1 a non-fatal error is output.  If the flag
   is >1 a fatal error is output and program execution stopped.              */
 
int gen_op::OpCheck(const gen_op &Op1, int warn) const
  {
  if(dim() == Op1.dim()) return 1;      // Check dimension equivalence
  if(warn)
    {
    if(warn > 1) GenOpfatality(51);     // Dimensioning mismatch!
    else         GenOperror(51,1);
    }
  return 0;
  }
 
bool gen_op::OpCheck(const matrix& mx, int warn) const
  {
  if(mx.cols()!=mx.rows())		// If array is not square it cannot
    {					// interact with with GenOp
    if(warn)				//   Output warnings, maybe quit
      {
      if(warn > 1) GenOpfatality(50);
      else         GenOperror(50,1);
      }
    return false;
    }
  return true;
  }  


int gen_op::OpCheck(const basis &bs, int warn) const
  {
  if(WBR->RepMx.cols() == bs.dim()) return 1; // OK if dimensions match
  if(warn)
    {
    if(warn > 1) GenOpfatality(50);
    else         GenOperror(50,1);
    }
  return 0;
  }  


int gen_op::OpCheck(const matrix &mx1, const matrix &mx2, int warn) const
  {
  int TF  = OpCheck(mx1,1); 		// Insure square matrix mx1
      TF *= OpCheck(mx2,1); 		// Insure square matrix mx2
      TF *= (mx1.cols() == mx2.cols());	// Inxure mx1 & mx2 same size
  if(TF) return 1; 
  if(warn)
    {
    if(warn > 1) GenOpfatality(51);
    else         GenOperror(51,1);
    }
  return 0;
  }  

int gen_op::OpCheck(const matrix &mx, const basis &bs, int warn) const
  {
  int TF  = OpCheck(mx,1); 		// Insure square matrix mx1
      TF *= (mx.cols() == bs.dim());	// Insure mx & bs same size
  if(TF) return 1; 
  if(warn)
    {
    if(warn > 1) GenOpfatality(51);
    else         GenOperror(51,1);
    }
  return 0;
  }  


/* This set of functions insures that the Op limits are met, Op element access
   is within matrix representation bounds, etc.  These will issue non-fatal
   warnings and continue on.                                                 */

int gen_op::OpCheck(int row, int col, int warn) const
  {
  int TF=1;
  if(row<0 || row>=dim()) TF*=0;	// Check row index
  if(col<0 || col>=dim()) TF*=0;	// Check col index
  if(TF) return 1;
  if(warn)
    {
    if(warn > 1) GenOpfatality(56);
    else         GenOperror(56,1);
    }
  return 0;
  }


int gen_op::LimCheck(int limit, int warn) const
  {
  if(limit >= 3) return 1;			// Insure minimum # reps is 3
  if(warn)
    {
    if(warn > 1) GenOpfatality(57);
    else         GenOperror(57,1);
    }
  return 0;
  }

// ____________________________________________________________________________
// K                    CLASS OPERATOR I/O FUNCTIONS
// ____________________________________________________________________________

// ------------------------ ASCII Output Functions ----------------------------

/*              Input           Op   : Operator (this)
                                ostr : Output ASCII file stream
                                full : Flag for amount of output
                Return          void : Op is sent to the output stream
                Note                 : Op WBR representation only            */  

std::ostream& gen_op::print(std::ostream& ostr, int full) const
  {
  if(!WBR)
    {
    std::string emp("Empty");
    if(OpName.length()) emp += std::string(" "); 
    ostr << CenterString(emp + OpName + " Operator") << "\n";
    }
  else
    {
    ostr << CenterString(OpName + " Operator") << "\n";
    ostr << *WBR;
    }
  if(full) ostr << "\n";
  return   ostr;
  }
 
std::ostream& operator<< (std::ostream& ostr, const gen_op &Op) {return Op.print(ostr);}

std::ostream& gen_op::eigenvalues(std::ostream& ostr, int rc, int ncol) const
  {
  if(!dim()) return ostr;			// Bail if no eigenvalues
  row_vector ev = eigvals();			// Get vector of eigenvalues
  return ev.printcols(ostr,ncol,rc);		// Print vector elements
  }
  

// ------------------------ Binary Output Functions ---------------------------

/*              Input           Op   : Operator (this)
                                fn   : Input binary filename
                                fp   : File stream (pointing at Op)
                Return          void : Op is written to either the
                                       specified file or filestream.
                Note                 : Output format is partially set by
				       class matrix (matrix typing)  
		Note		     : NULL operators will NOT be output 
                Note                 : Only WBR representation is output     */  

void gen_op::write(const std::string& fn) const
  {
  std::ofstream fp;					// Construct a file stream
  fp.open(fn.c_str(),std::ios::out|std::ios::binary);	// Open the file for output
  write(fp);					// Write Op to file stream
  fp.close();					// Close file stream
  }

std::ofstream& gen_op::write(std::ofstream& fp) const
  { if(WBR) (*WBR).write(fp); return fp; }

// ------------------------ Binary Input Functions ----------------------------

/*              Input           Op   : Operator (this)
                                fn   : Input binary filename
                                fp   : File stream (pointing at Op)
                                bs   : An optional basis in which to
                                       force the operator into so that many
                                       read operators may share one basis!
                Return          void : Op is read in from either the
                                       specified file or filestream.
                Note                 : Read Op will reside in only 1 basis   */

void gen_op::read(const std::string& fn, const basis& bs)
  {
  read(fn);				// Use function overload
  if(get_basis() == bs) put_basis(bs);	// Set basis to bs if match up
  }

void gen_op::read(const std::string& fn)
  {
  std::ifstream fp;					// Construct a file stream
  fp.open(fn.c_str(),std::ios::in|std::ios::binary);	// Open file for reading
  read(fp);					// Read Op, use funct. overload
  fp.close();					// Close the file stream
  }

std::ifstream& gen_op::read(std::ifstream &fp)
  {
  ZeroOp();				// Remove any existing OpReps
  genoprep OpRep;			// An operator reperesentation
  OpRep.read(fp);			// Read the operator matrix
  AddRep(OpRep);			// Construct the operator
  return fp;
  }

// ____________________________________________________________________________
// L                     CLASS OPERATOR TESTING FUNCTIONS
// ____________________________________________________________________________

/* Since GAMMA operators take the brunt of the work in most MR simulations, it
   is essential that they are in proper working order. As such, I have included
   some means to test various essential aspects of them. In particular, two
   basic aspects are the ability to handle their diagonalization and multiple
   representation (OpRep) tracking. Some of these tests simply parallel the
   functionality provided in class matrix.
 
     Function    Arguments                         Result 
 --------------- --------- ------------------------------------------------ 
 TestEigenSystem    int    Diagonalize the operator & check its eigensystem
 TestReps           int    Check operator equivalance over representations    */

 
double gen_op::TestEigenSystem(int pf)
  {
  set_DBR();					// Put Op in default basis
  matrix mxDBR = get_matrix();			// Get DBR OpRep matrix
  set_EBR();					// Put Op in eigenbasis
  matrix mxEBR = get_matrix();			// Get EBR OpRep matrix
  basis  bsEBR = get_basis();			// Get EBR OpRep basis
  double bsdev = norm(bsEBR.TestBasis(pf));	// Do basis tests U*inv(U)=1
  matrix mx = bsEBR.convert_back(mxEBR);	// Put EBR OpRep back into DBR
  matrix delmx = mx - mxDBR;			// Difference between them
  double mxdev = norm(delmx.maxZ());		// Max difference between them
  if(pf)
    {
    std::cout << "\n\t\tLargest Deviation From Op - S*D*Sm1: " << mxdev;
    }
  return gmax(bsdev, mxdev);
  }

/*
double gen_op::TestReps(int pf)
  {
// sosi
  }
*/

// ----------------------------------------------------------------------------
//                     Matrix Functionality For GenOp
// ----------------------------------------------------------------------------

/* This allows GenOp to act like matrix when using these testing functions.
   The difference is that these act on the current basis representation.

     Function    Output                       Description
   ------------  ------  ------------------------------------------------------
   is_symmetric   bool    TF if <i|mx|j>-<j|mx|i>  < d       (def. d GMxCut)
   is_hermitian   bool    TF if <i|mx|j>-<j|mx|i>* < d       (def. d GMxCut)
   is_unitary     bool    TF if inv(mx) == adjoint mx, CPU intensive
   is_real        bool    TF if Im(<i|mx|j>) < d for all i,j (def. d GMxCut)
   is_imaginary   bool    TF if Re(<i|mx|j>) < d for all i,j (def. d GMxCut)
   is_complex     bool    TF if is_real && is_imaginary      (def. d GMxCut)
   is_zero        bool    TF if ||<i|mx|j>|| < d for all i,j (def. d GMxCut)
   is_diagonal    bool    TF if ||<i|mx|j>|| < d for all i!=j(TRUE,d GMxCut)
   is_square      bool    TF if rows_ == cols_                               */

bool gen_op::is_symmetric(const double d) const { return WBR->RepMx.is_symmetric(d); }
bool gen_op::is_hermitian(const double d) const { return WBR->RepMx.is_hermitian(d); }
bool gen_op::is_unitary(const   double d) const { return WBR->RepMx.is_unitary(d);   }
bool gen_op::is_real(const      double d) const { return WBR->RepMx.is_real(d);      }
bool gen_op::is_imaginary(const double d) const { return WBR->RepMx.is_imaginary(d); }
bool gen_op::is_complex(const   double d) const { return WBR->RepMx.is_complex(d);   }
bool gen_op::is_zero(const      double d) const { return WBR->RepMx.is_zero(d);      }
bool gen_op::is_diagonal(const  double d) const { return WBR->RepMx.is_diagonal(d);  }
bool gen_op::is_square()                  const { return true;                }

 
// ____________________________________________________________________________
// M                 CLASS OPERATOR CONTAINER SUPPORT FUNCTIONS
// ____________________________________________________________________________

/* Aside for providing basic tests as to whether to operators are equvalent or
   not, these operators are necessary if any STL container classes are to
   be used based on operators (e.g. list<gen_op> or vector<gen_op>)          */

bool gen_op::operator== (const gen_op& Op) const
  {
  if(!WBR) return (!Op.WBR);		// Return TRUE if both NULL
  if(!Op.WBR) return false;		// Return FALSE if Op NULL & not Op
  const genoprep* REP = NULL;		// Don't know if any OpRep's match
  const genoprep* PTR = Op.ObeginC();	// Begin with 1st OpRep of Op
  std::vector<genoprep>::const_iterator item;// Iterator into list
  item = Op.begin();
  while(item != Op.end())		// Loop through OpReps of Op
    {
    PTR = &(*item);
    REP = FindRep(PTR->RepBs);		// Look for OpRep in Op with same basis
    if(REP)				// If we have found match, then the Ops
      return(PTR->RepMx == REP->RepMx); // are equal if their matrices match
    item++;
    }
  basis thisbs = WBR->RepBs; 		// No matching OpReps, so as a last 
  basis Opbs   = Op.WBR->RepBs; 	// resort we must see if DBR's match
  set_DBR();				// So, store initial bases, set both
  Op.set_DBR();				// into teir DBR, compare, then reset
  bool TF=(DBR->RepMx == Op.DBR->RepMx);// them back into their original bases
  return TF;				// before returning the result
  }

bool gen_op::operator!=(const gen_op& Op) const { return !(*this==Op); }
bool gen_op::operator<(const  gen_op& Op) const { return (size()<Op.size()); }
bool gen_op::operator>(const  gen_op& Op) const { return (size()>Op.size()); }


#endif						// GenOp.cc
