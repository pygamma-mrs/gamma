/* MultiExch.cc *************************************************-*-c++-*-
**									**
**                                G A M M A 				**
**									**
**      Multiple System Exchange 			Implementation	**
**									**
**      Copyright (c) 1995						**
**      Nikolai Skrynnikov						**
**      Dr. Scott A. Smith						**
**      National High Magnetic Field Laboratory				**
**      1800 E. Paul Dirac Drive					**
**      Tallahassee Florida, 32306-4005					**
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
** This module of function supports the multi_sys, the GAMMA class 	**
** handling mulitple spin systems.  The routines herein	generally	**
** involve such a spin system and builds common non-mutual exchange 	**
** superoperators, in this case in a direct product Liouville space of	**
** the systems involved.						**
**									**
*************************************************************************/

#ifndef   MultiExch_cc			// Is the file already included?
#  define MultiExch_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <MultiSys/MultiExch.h>		// Include out header file
#include <MultiSys/MultiSys.h>		// Include multi_sys spin systems
#include <MultiSys/MultiLib.h>		// Include the D_basis function
#include <HSLib/SpinOpCmp.h>		// Include the Fx
#include <Level2/MutExch.h>		// Include the Kex function

using std::vector;			// Using libstdc++ STL vectors
using std::ostream;			// Using libstdc++ output streams
using std::string;			// Using libstdc++ strings

/*************************************************************************
**									**
**			      Function Xnm				**
**									**
** The purpose of this function is to calcualte a superoperator that	**
** represents non-mutual exchange between the spin system components 	**
** in a given multi_sys spin system.  The input spin system, msys, 	**
** contains a set of defined non-mutual exchange processes, each of	**
** which contains the spin systems exchanging, their exchange rates,	**
** and their spin<-->spin mappings.					**
**									**
** The general principle is to compute the element of exchange matrix	**
** Xex between								**
**									**
**	        LHS				    RHS			**
** |m1, m2,...mN><m'1, m'2,...m'N|  &  |n1, n2,...nL><n'1, n'2,...n'L|	**
**									**
** Some of the spins from the left hand side (LHS) are exchanging with	**
** their counnterparts from the right. For example 2nd spin from the	**
** left can be mapped by exchange onto 5th spin from the right. The	**
** exchange must preserve the spin isotope type and its quantum state: 	**
**									**
**                        m2 = n5, m'2 = n'5				**
**									**
** Some other spins are not mapped between left and right. For example,	**
** 1st spin from the left and 7th spin from the right (they are left 	**
** for some third molecule). These spins must be in |m><m| state both 	**
** on the right and on the left, i.e.					**
**									**
**                       m1 = m'1 and n7 = n'7				**
**									**
** If these conditions are fullfilled, the element of Xex exits.	**
**									**
** In the basis of conveniently normalized spin operators (multiplied	**
** by population of respective species) the exchange matrix assumes 	**
** the following form:							**
**									**
**                  	[   -k1      k2     ]				**
**		        |   	            |				**
**			[    k1     -k2     ]				**
**									**
**									**
** where k factors multiply the whole blocks between (or within) 	**
** exchanging components and are related to each other as		**
**									**
**                     p1 k1 = p2 k2 = p3 k3 = ...			**
**									**
** Herein we will normalize the populations, that is to say that the	**
** exchange matrix for p1=10 and p2=5 will be the same as that for	**
** p1=2 and p2=1. This is done by finding out the component with the	**
** smallest population (p   ) then using				**
**                       min                                            **
**				[        ]				**
**	                   k  = | p /p	 | K				**
**                          i   [  j  min]                              **
**									**
** to determine k1 and k2 for an exchange between components i and j	**
** having the exchange rate K.						**
**									**
** For more information on the theory, see 				**
**									**
**   R.O.Kuhne, T.Schaffhauser, A.Wokaun, R.R.Ernst JMR 35, 39 (1979)	**
**									**
** In addition, blocks have normalization associated with HS dimension	**
** of the component spin systems.					**
**									**
** N. Skrynnikov, Oct. 1995						**
**									**
*************************************************************************/


	// Input		msys    : A multi_spin spin system
	// Output		XL	: Non-mutual exchange superoperator
	// Note				: This superoperator is returned
	//				  in the default (product) basis

// ____________________________________________________________________________
//                   Full Non-Mutual Exchange Superoperator
// ____________________________________________________________________________

/* The non-mutual exchange superoperator is determined from the sum over all
   non-mutual exchange processes defined in the system. Thus we have
                          
                                 processes
                                    ---
                                    \
                               X =  /    X
                                    ---   p
                                     p                                        */

super_op Xnm(const multi_sys& sys)
  {
  int   dim = sys.LS();				// Determine overall LS dim
  matrix Xmx(dim,dim,complex0);			// Superop matrix to host Xex
  for(int p=0; p<sys.NExProcs(); p++)		// Loop over exchange processes
    Xmx += Xnmp(sys, p);			// Add this exchange process
  basis Dbs = D_basis(sys);	        	// Default basis (composite LS)
  return super_op(Xmx,Dbs);			// Return as superoperator
  }

// ____________________________________________________________________________
//                   Single Process Non-Mutual Exchange Superoperator
// ____________________________________________________________________________

/* This function generates the non-mutual exchange superoperator matrix for a
   specified exchange process p. The calculation is performed in blocks where
   each block corresponds to two components in exchange. Square blocks about
   the diagonal is a component with itself and will be zero unless it is
   involved in the exchange.  If the latter is the case then only diagonal
   elements in the block will may be non-zero. Off-diagonal (maybe rectangular)
   blocks are for two different components in exchange. If the two components
   do not exchange in process p then the block is left empty. If they do 
   exchange the elements in the block are calculated.

   Blocks for component A & B are of dimension LSA x LSB whtere LSA is the
   Liouville space of the component spin system A. Given the base global 
   indices of each component in Liouville space, IAo and JBo respectively, the
   block of elements for A exchanging with B is show below

                     <IAo|Xp|JBo>    ------    <IAo|Xp|JBo+LSB-1>
                          |                           |
                          |                           |
                          |                           |
                 <IAo+LSA-1|Xp|JBo>  ------ <IAo+LSA-1|Xp|JBo+LSB-1>         */

matrix Xnmp(const multi_sys& msys, int p)
  {
  ExchProc Pro = msys.ExProc(p);		// This exchange process
  double   K   = Pro.Kex();			// Process exchange rate (1/sec)
  if(K == 0) 					// We are done if there
    return matrix(msys.LS(),msys.LS(),0);	// is no exchange!

  int ncmps    = msys.NComps();			// # comps. in the system
  //int ncl      = msys.NCompsLHS(p);		// # comps. on left,  ex. proc p
  //int ncr      = msys.NCompsLHS(p);		// # comps. on right, ex. proc p
  //int comp0    = Pro.LHSComp(0);		// 0th component on lhs
  double pmin  = msys.popmin();                 // Smallest component population
  int    dim   = msys.LS();			// Determine overall LS dim
  basis Dbs    = D_basis(msys);	        	// Default basis
  matrix Xp(dim, dim, 0.0);			// Start with zero superop

  double popI, popJ;
  vector<int> lsds  = msys.LSs();		// Array of component LS dims.
  int cmpI, cmpJ;				// Liouville basis indices			
  int Jo=0, Jend=0;				// Global Xnmp column indices
  int Io=0, Iend=0, I, J;			// Global Xnmp row indices
  bool cmpIex=false, cmpJex=false;		// Flag if comps. I/J exchange
  double KI, KJ;
  for(cmpI=0,Io=0,Iend=0; cmpI<ncmps; cmpI++)	// Loop Xnmp row components
    {
    Io      = Iend;                             //  Start of left comp. block
    Iend   += lsds[cmpI];                       //  End   of left comp. block
    cmpIex  = Pro.CompInLHS(cmpI);              //  Flag if comp. exchanging
    if(cmpIex)
      {
      popI = msys.pop(cmpI);                    //  Population of left comp.
      Jend = 0;
      for(cmpJ=0; cmpJ<ncmps; cmpJ++)
        {
        Jo      = Jend;                         //  Start of rght comp. block
        Jend   += lsds[cmpJ];                   //  End   of rght comp. block
        cmpJex = Pro.CompInRHS(cmpJ);           //  Flag if comp. exchanging
        if(cmpJex)
          {
          popJ = msys.pop(cmpJ);		//  Population of rt comp.
          KI = popJ*K/pmin;
          KJ = popI*K/pmin;
          for(I=Io; I<Iend; I++)		//  Fill the <I|Xp|I> elements
             Xp.put(Xp.get(I,I) + KI, I, I);	//  with scaled K values
          Xnmpblk(msys, Pro, Xp, KJ,		//    Generate cross-terms in 
                cmpI, cmpJ, Io, Iend, Jo, Jend);//    block cmpI w/ block cmpJ
          Xnmpblk(msys, Pro, Xp, KI,		//    Generate cross-terms in 
                cmpJ, cmpI, Jo, Jend, Io, Iend);//    block cmpI w/ block cmpJ
          for(J=Jo; J<Jend; J++)		//  Fill the <J|Xp|J> elements
             Xp.put(Xp.get(J,J) + KJ, J, J);	//  with scaled K values
          }					//     Next column component
        Jo = Jend;				//     Start of next col comp.
        }					//     Move to next col comp.
      }						//     Move to next col comp.
    Io = Iend;					//   Start of next row comp.
    }						//   Move  to next row comp.
  return Xp;					// Return process exchange mx
  }

// ____________________________________________________________________________
//                        Exchange Superoperator Block
// ____________________________________________________________________________

/* This functions fills all the elements of a single block in a non-mutual
   exchange superoperator. The superoperator will be associated with multiple
   spin system msys and for a single exchange process Pro. Elements in the
   superoperator Xp will be filled for the block between cmpI and cmpJ.
   The exchange rate K is modified from the process Pro exchange rate to 
   account for any differences in the component populations.

   Each two component spin systems in msys, cmpI and cmpJ, define a block of
   Xp (in spin Liouville space). The row length (height) is the spin Liouville
   space of the 1st component, cmpI.  The column length (width) is the spin
   Liouville space of the 2nd component, cmpJ. The start of this block will
   be <Io|Xp|Jo> and the end of the block is <Io+Iend-1|Xp|Jo+Jend-1>.       */

void Xnmpblk(const multi_sys& msys, const ExchProc& Pro, matrix& Xp,
              double K, int cmpI, int cmpJ, int Io, int Iend, int Jo, int Jend)
  {
  int I, J;					// Global Xp indices
  row_vector braI, ketI;			// for |I><I| basis function
  row_vector braJ, ketJ;			// for |J><J| basis function
  int hsnorm;					// Normalization by dim of HS
  bool match;					// Flag if non-zero exch. elem.
  for(I=Io; I<Iend; I++)			//   Loop <I|Xp|J> rows
    { 						//   for process p 
    braI = LS_qState_bra(msys, I);		//   Essence of |ketI><braI|
    ketI = LS_qState_ket(msys, I);		//   form superbra & superket
    for(J=Jo; J<Jend; J++) 			//     Loop <I|Xp|J> columns
      { 					//     for process p 
      braJ  = LS_qState_bra(msys, J);		//     Essence of |ketJ><braJ|
      ketJ  = LS_qState_ket(msys, J);		//     form superbra & superket
      match = Xnmpelem(msys,Pro,braI,ketI,	//     See if element has exch.
                    braJ,ketJ,cmpI,cmpJ,hsnorm);//     (set hsnorm if it is)
      if(match) 				// Exchange will swap I & J
        {					// so we must set <I|Xp|J>
        Xp.put(-K/double(hsnorm),I,J); 		// to normalized rate good
        }					//   for any non-mutual exch. 
      }						// Next column of Xp in comp.
    }						// Next column component
  }

// ____________________________________________________________________________
//               Exchange Superoperator Off-Diagonal Element
// ____________________________________________________________________________

/* This function checkes whether a particular element in a non-mutual exchange
   superoperator may be non-zero. This exchange supoeroperator is associated
   with the multiple spin system msys and the partiuclar exchange process Pro.
   The element <I|Xp|J> is specified by two Hilbert space bras and kets, i.e.
   <I|Xp|J> = <|i><i||Xp||j><j||. These are input as braI, ketI, braJ, ketJ
   respectively. Each off-diagonal element must be associated with exchange
   between two component spin systems cmpI and cmpJ (where cmpI != cmpJ).
   The function returns false if <I|Xp|J> is zero for exchange process Pro. If
   the element is non-zero, the function returns true and sets the value of
   the element normalization, hsnorm, based on spin Hilbert spaces.         

   To decide if an element is non-zero a comparison is made between the spin
   quantum states of two components spins. The Hilbert space bra and ket of
   component A is given by

         |A> = |Am  Am  ... Am  >   and   <A| = <Am  Am  ... Am  |
                  z1  z2      zn                   z1  z2      zn
  
   where n is the number of spins in A. An element of Xp can be considered to
   be
                    '   '   '          ''   ''   ''      '''  '''  '''
     <|A>|B>|C>...<A |<B |<C |...|Xp||A  >|B  >|C  >...<A   |B   |C   |...|

   where the primes indicate that the associated mz values are different. Now
   consider an element <I|Xp|J> which happens to line in the area where say
   component B exchanges with component D. We know that, within the block 
   associated with B & D, the spin quantum numbers will match for any other 
   components (since they are unchanging within the block). Thus we might 
   write
                                           '        '
                         <I|Xp|J> == <|B><B||Xp|D><D||

   That is, only the spin quantum numbers for the exchange components need to
   be considered. Expanding over the individual spin quantum numbers we have

                                 '   '       '                      '      '
 <I|Xp|J> = <|Bm  Bm  .. Bm  ><Bm  Bm  ... Bm  ||Xp||Dm  .. Dm  ><Dm  .. Dm  |>
                z1  z2     zn    z1  z2      zn        z1     zl    z1     zl

   An exchange process will potentially swap spins on the left with spins on
   the right. 
*/

bool Xnmpelem(const multi_sys& msys, const ExchProc& Pro,
       const row_vector& braI, const row_vector& ketI,
       const row_vector& braJ, const row_vector& ketJ,
                                        int cmpI, int cmpJ, int& hsnorm)
  {
  hsnorm = 1;					// Starting normalization
  bool passive;					// Flag if spin in exchange
  double qIb = 0, qIk = 0;			// Mz quantum numbers for sI
  double qJb = 0, qJk = 0;			// Mz quantum numbers for sJ
  int sI, sJ;					// Spin indices
  int ncI = braI.cols();			// Number spins in braI
  int ncJ = braJ.cols();			// Number spins in braJ

/* ----------------------------------------------------------------------------
   Before Calling This Function We Insured That Some Spins In Component cmpI
   Were In Exchange With Spins In Componet cmpJ. Any Spin sI In Component cmpI
   May Have At Most 1 Exchange Spin Partner, sJ, In Component cmpJ. Here We
   Loop Over All Spins sI In cmpI And Find Their Corresponding Spin In cmpJ,
   If Any. If Spin sI Has No Exchange Partner, It Cannot Change Spin States In
   The Basis Function Coherence, (<cmpI, mzsI| == |cmpI, mzsI|), Else The 
   Element Must Be Zero. If Spin sI Is Passive In The Exchange Then It Will
   Contribute To The Norm Of Any Exchange Matrix Element. Alternatively, If We
   Find An Exchanging Spin Pair We Insure That The Assosicated Mz States Are
   Swapped Or Else The Element Must Be Zero.
   --------------------------------------------------------------------------*/

  for(sI=0; sI<ncI; sI++)			// Loop spins in braI & ket I
    {						// (i.e. spins in component I)
    passive = true;				//   Assume spin sI passive
    qIb = braI.getRe(sI);			//   Mz of spin sI in braI <i|
    qIk = ketI.getRe(sI);			//   Mz of spin sI in ketI |i>
    for(sJ=0; sJ<ncJ && passive; sJ++)		//   Loop spins in braJ & ket J
      { 					//   (i.e. spins in compon. J)
      if(Pro.mapped(cmpI,sI,cmpJ,sJ))		//     If partners in exchange
        {					//     cmpI sI <--> cmpJ sJ
        passive = false;
        qJb = braJ.getRe(sJ);			//        Mz spin sJ, braJ <j|
        qJk = ketJ.getRe(sJ);			//        Mz spin sJ, ketJ |j>
        if((qIb!=qJb) || (qIk!=qJk))		//        Unequal quantum state
          return false;				//        element must be zero
        }
      }
    if(passive)					// If cmpI spin sI not exchange
      {						// with any spin in comp. cmpJ
      if(qIb!=qIk) return false;		//   No element if not |m><m|
      hsnorm *= (msys.Comp(cmpI)).HS(sI); 	//   Else a passive spin will
      } 					//   contribute to norm factor 
    }

// ----------------------------------------------------------------------------
// At This Point We Know That For All Spins In Components cmpI And cmpJ This
// Exchange Superoperator Element Has Properly Swapped Spin Quantum Numbers In
// Its Two Basis Functions. We Have Also Set A Proper Normalization Factor.
// ----------------------------------------------------------------------------

  for(sJ=0; sJ<ncJ; sJ++) 			// Loop spins in braJ & ketJ
    {						// (i.e. spins in component J)
    qJb = braJ.getRe(sJ);			//   Mz of spin sJ in braJ
    qJk = ketJ.getRe(sJ);			//   Mz of spin sJ in ketJ
    passive = true;				//   Assume spin sJ is passive
    for(sI=0; sI<ncI && passive; sI++)		//   Loop spins in braI & ketI
      {						//   (i.e. spins in comp I)
      if(Pro.mapped(cmpJ,sJ,cmpI,sI))		//     See if sJ & sI exchange
        passive = false;			//     & flag non-passive if so
      }
    if(passive)					//   If sJ is not exchanging,
      if(qJb != qJk) return false; 		//   No element if not |m><m|
    }
  return true;
  }

// ----------------------------------------------------------------------------
//  This Function Mimics Xnm Above But Just Outputs Details Of The Calculation
// ----------------------------------------------------------------------------

void Xnm(ostream& ostr, const multi_sys& sys)
  {
  string hdr("Non-Mutual Exchange Superoperator Generation");
  int len = hdr.length();
  string spcr(40-len/2, ' ');
  ostr << "\n\n" << spcr << hdr;
  ostr << "\n"   << spcr << string(len, '=');
  ostr << "\n";
  ostr << "\n\t* Liouville Space:    " << sys.LS(); 
  ostr << "\n\t* Exchange Processes: " << sys.NExProcs();
  ostr << "\n\t* System Components:  " << sys.NComps();
  for(int p=0; p<sys.NExProcs(); p++)
    {
    ostr << "\n\t* Calculating Xnmp For Process: " << p;
    Xnmp(ostr, sys, p);
    }
  }

void Xnmp(ostream& ostr, const multi_sys& msys, int p)
  {
  ExchProc Pro = msys.ExProc(p);		// This exchange process
  double   K   = Pro.Kex();			// Process exchange rate (1/sec)
  int      ls  = msys.LS();			// Liouville space dimension
  if(K == 0) 					// We are done if there
    {
    ostr << "\n\t - Zero Exchange Rate Specified!";
    ostr << "\n\t - Process Array Is A Zero " 
         << ls << "x" << ls << " Matrix";
    return;
    }
  int ncmps = msys.NComps();			// # comps. in the system
  int ncl   = msys.NCompsLHS(p);		// # comps. on left,  ex. proc p
  int ncr   = msys.NCompsLHS(p);		// # comps. on right, ex. proc p
  //int comp0 = Pro.LHSComp(0);			// 0th component on lhs
  double pmin  = msys.popmin();			// Smallest component population
  ostr << "\n\t - Exchange Rate Is "                       << K << "/sec";
  ostr << "\n\t - Components In Exchange On The Left Is "  << ncl;
  ostr << "\n\t - Components In Exchange On The Right Is " << ncr;
  ostr << "\n\t - Population Of Smallest Component Is "    << pmin;

  double popI, popJ;
  vector<int> lsds  = msys.LSs();		// Array of component LS dims.
  int cmpI, cmpJ;				// Liouville basis indices			
  int Jo=0, Jend=0;				// Global Xnmp column indices
  int Io=0, Iend=0;				// Global Xnmp row indices
  bool cmpIex=false, cmpJex=false;		// Flag if comps. I/J exchange
  for(cmpI=0; cmpI<ncmps; cmpI++)		// Loop left side components
    {
    Io      = Iend;				//  Start of left comp. block
    Iend   += lsds[cmpI];			//  End   of left comp. block
    cmpIex  = Pro.CompInLHS(cmpI);		//  Flag if comp. exchanging
    if(cmpIex)
      {
      popI = msys.pop(cmpI);			//  Population of left comp.
      Jend = 0;					
      for(cmpJ=0; cmpJ<ncmps; cmpJ++)
        {
        Jo      = Jend;				//  Start of rght comp. block
        Jend   += lsds[cmpJ];			//  End   of rght comp. block
        cmpJex = Pro.CompInRHS(cmpJ);		//  Flag if comp. exchanging
        if(cmpJex)
          {
          popJ = msys.pop(cmpJ);
          ostr << "\n\t - Treating Exchange Between Component "
               << cmpI << ", Population " << msys.pop(cmpI);
          ostr << "\n\t   And Component "
               << cmpJ << ", Population " << msys.pop(cmpJ);
	  ostr << "\n\t - Four Blocks Of Exchange Mx To Be Filled";
	  ostr << "\n\t   1. Comp " << cmpI << " Disspation";
	  ostr << "\n\t      Rows " << Io << " To " << Iend;
	  ostr << "\n\t      Cols " << Io << " To " << Iend;
	  ostr << "\n\t      Effective K (k1) Value " << (popJ/pmin)*K;
          Xnmpdblk(ostr, msys, (popJ/pmin)*K, Io, Iend);
	  ostr << "\n\t   2. Comp " << cmpI << " Exchange With Comp " << cmpJ;
	  ostr << "\n\t      Rows " << Io << " To " << Iend;
	  ostr << "\n\t      Cols " << Jo << " To " << Jend;
	  ostr << "\n\t      Effective K (-k2) Value " << -(popI/pmin)*K;
          Xnmpblk(ostr, msys, Pro, (popI/pmin)*K, cmpI, cmpJ, Io, Iend, Jo, Jend);
	  ostr << "\n\t   3. Comp " << cmpJ << " Exchange With Comp " << cmpI;
	  ostr << "\n\t      Rows " << Jo << " To " << Jend;
	  ostr << "\n\t      Cols " << Io << " To " << Iend;
	  ostr << "\n\t      Effective K (-k1) Value " << -(popJ/pmin)*K;
          Xnmpblk(ostr, msys, Pro, (popJ/pmin)*K, cmpJ, cmpI, Jo, Jend, Io, Iend);
	  ostr << "\n\t   4. Comp " << cmpJ << " Disspation";
	  ostr << "\n\t      Rows " << Jo << " To " << Jend;
	  ostr << "\n\t      Cols " << Jo << " To " << Jend;
	  ostr << "\n\t      Effective K (k2) Value " << (popI/pmin)*K;
          Xnmpdblk(ostr, msys, (popI/pmin)*K, Jo, Jend);
          }
        Jo = Jend;				//   Start next rgt. comp. blk
        }
      }
    Io = Iend;					//   Start next lft. comp. blk
    }

/*
  double pop0  = msys.pop(comp0);	        // population of 0th component
  ostr << "\n\t - Base Population Is " << pop0;
  Jo=0, Jend=0;					// Global Xnmp column indices
  Io=0, Iend=0, I;				// Global Xnmp row indices
  cmpIex=false, cmpJex=false;			// Flag if comps. I/J exchange
  for(cmpI=0,Io=0,Iend=0; cmpI<ncmps; cmpI++)	// Loop Xnmp row components
    {
    Iend += lsds[cmpI];				// Index at end of row block
    cmpIex = Pro.involves(cmpI);		// Flag if comp. exchanging
    ostr << "\n\t - Generating Blocks Involving Component "
         << cmpI << ", Population " << msys.pop(cmpI);
    ostr << "\n\t - Block Spans Rows " << Io << " To " << Iend-1;
    if(!cmpIex)
      ostr << "\n\t - Component Not In Exchange, Skipping.....";
    for(cmpJ=0,Jo=0,Jend=0; 			// Loop Xnmp colunm components
                  cmpJ<ncmps && cmpIex; cmpJ++)	// if cmpI is exchanging
      {
      Jend   += lsds[cmpJ];			// Index at column end of block
      cmpJex  = Pro.involves(cmpJ);		// Flag if comp. exchanging
      popJ    = msys.pop(cmpJ);			// Pop. of component exchanging
//if(cmpI == cmpJ) kexJo = popJ*K/pmin;
//else             kexJo = pop0*K/pmin;
      kexJo   = (pop0*K)/popJ;			// Normalized exchange rate
      ostr << "\n\t   . Calculating Exchange"
           << " With Component " << cmpJ << ", Population " << popJ;
      ostr << "\n\t   . Block Spans Rows " << Io << " To "  << Iend-1
           << " And Columns "              << Jo << " To "  << Jend-1;
      ostr << "\n\t   . Scaled Exchange Rate Is "     << kexJo << "/sec";
      if(cmpI == cmpJ)				// 1) For "mutual exchange"
        {					//    only diagonal decay terms,
        ostr << "\n\t   . Block Is Diagonal, Only Diagonal Elements NonZero";
        bool nr   = true;
        int  nc   = 0;
        for(I=Io; I<Iend; I++)
          {
          if(nr) { ostr << "\n\t     "; nr=false; }
          ostr << "  <" << I << "|Xnmp|" << I << "> = " << kexJo;
          nc++;
          if(nc >= 4) { nc=0; nr=true; }
          }
        Jo = Jend;				//    Set J start to next col comp.
        }
      else if(!cmpJex)				// 2) If component J does not
        { 					//    exchange block is zero
        ostr << "\n\t   . Block " << cmpJ << " Is Not Exchanging, Left Zero";
        Jo = Jend;
        }
      else if(!Pro.mixes(cmpI, cmpJ)) 		// 3) If components I & J don't
        { 					//    exchange block is zero
        ostr << "\n\t   . Blocks " << cmpI << " & " << cmpJ
             << " Do Not Exchanging, Left Zero";
        Jo = Jend;
        }
      else					// 4) Components in I & J are
        {					//    exchanging
        Xnmpblk(ostr, msys, Pro, kexJo, 	//    Generate cross-terms in 
                cmpI, cmpJ, Io, Iend, Jo, Jend);//    block cmpI w/ block cmpJ
        Jo = Jend;				//     Start of next col comp.
        }					//     Next column component
      }						//     Move to next col comp.
    Io = Iend;					//   Start of next row comp.
    }						//   Move  to next row comp.
*/
  }

void Xnmpdblk(ostream& ostr, const multi_sys& msys, double K, int Io, int Iend)
  {
  ostr << "\n\t      Block Is Diagonal, Only Diagonal Elements NonZero";
  bool nr   = true;
  int  nc   = 0;
  int I;
  for(I=Io; I<Iend; I++)
    {
    if(nr) { ostr << "\n\t     "; nr=false; }
    ostr << "  <" << I << "|Xnmp|" << I << "> = " << K;
    nc++;
    if(nc >= 4) { nc=0; nr=true; }
    }
  }

void Xnmpblk(ostream& ostr, const multi_sys& msys, const ExchProc& Pro,
              double K, int cmpI, int cmpJ, int Io, int Iend, int Jo, int Jend)
  {
  int I, J;					// Global Xp indices
  row_vector braI, ketI;			// for |I><I| basis function
  row_vector braJ, ketJ;			// for |J><J| basis function
  int hsnorm;					// Normalization by dim of HS
  bool match;					// Flag if non-zero exch. elem.
  bool nr   = true;
  int  nc   = 0;
  ostr << "\n\t      Block Is Off-Diagonal, Displaying NonZero Elements";
  for(I=Io; I<Iend; I++)			//   Loop <I|Xp|J> rows
    { 						//   for process p 
    braI = LS_qState_bra(msys, I);		//   Essence of |ketI><braI|
    ketI = LS_qState_ket(msys, I);		//   form superbra & superket
    for(J=Jo; J<Jend; J++) 			//     Loop <I|Xp|J> columns
      { 					//     for process p 
      braJ  = LS_qState_bra(msys, J);		//     Essence of |ketJ><braJ|
      ketJ  = LS_qState_ket(msys, J);		//     form superbra & superket
      match = Xnmpelem(msys,Pro,braI,ketI,	//     See if element has exch.
                    braJ,ketJ,cmpI,cmpJ,hsnorm);//     (set hsnorm if it is)
      if(match) 				// Exchange will swap I & J
        {					// so we must set <I|Xp|J>
        if(nr) { ostr << "\n\t     "; nr=false; }
        ostr << "  <" << I << "|Xnmp|" << J << "> = " << -K/double(hsnorm);
        nc++;
        if(nc >= 4) { nc=0; nr=true; }
        }					//   for any non-mutual exch. 
      }						// Next column of Xp in comp.
    }						// Next column component
  }

// sosik

// ----------------------------------------------------------------------------
//     This is an older function for generating Xnm, left for comparison.
// ----------------------------------------------------------------------------

super_op XXnm(const multi_sys& mixture)
  {
  int   dim = mixture.LS();			// Determine overall LS dim
  basis Dbs = D_basis(mixture);	        	// Default basis
  matrix Xmx(dim, dim, 0.0);			// Superop matrices to host Xex
  int i, j, k;			
  int ic1=0, ic2=0;				// Multi_sys components
  int s1, s2;					// Spin numbers within comp's
  int match, mover;				// Flags
  double q1b = 0, q1k = 0;			// Mz quantum numbers of s1, s2
  double q2b = 0, q2k = 0;			// Mz quantum numbers of s1, s2
  double pop2, kex2;
  int spinHS, norm;				// Normalization by dim of HS
  ExchProc Pro;					// A single exchange process
  sys_dynamic dsys;
  row_vector bra1, ket1, bra2, ket2;
  for(i=0; i<mixture.NExProcs(); i++)		// Loop over exchange processes
    {
    Pro = mixture.ExProc(i);			// Get exchange process i
    matrix Xel(dim, dim, 0.0);			// Start with zero superop
    int ic0 = Pro.LHSComp(0);			// 0th component on lhs
    double pop0 = mixture.pop(ic0);	        // population of 0th component
    double kex0 = Pro.Kex();		        // k refers to 0th comp on lhs
    if(kex0 == 0) continue;			// Don't contribute if K=0
    for(j=0; j<dim; j++)			// Loop rows and columns of
      { 					// Xel for process i and set
      for(k=0; k<dim; k++)  			// <j|Xel|k> elements
	{
        ic1 = Dbs.which_sub_LS(j);		// Index of component 1
        ic2 = Dbs.which_sub_LS(k);	        // Index of component 2 
        pop2 = mixture.pop(ic2);
        kex2 = (pop0*kex0)/pop2;

        if(ic1 == ic2)				// Component with itself 
	  if(Pro.involves(ic1) )		// It's exchanging component
            {
            if(j == k) Xel.put(kex2, j, j); 	// Only diagonal decay terms
            continue;
            }
         if(!Pro.mixes(ic1, ic2)) continue;	// Not in exchange
						// Generate cross-terms
         bra1 = LS_qState_bra(mixture, j);	// |ket1><bra1| form superbra
         ket1 = LS_qState_ket(mixture, j);	//
         bra2 = LS_qState_bra(mixture, k);	//
         ket2 = LS_qState_ket(mixture, k);	// |ket2><bra2| form superket
         norm = 1;
         match = 1;
         for(s1=0; s1<bra1.cols(); s1++)	// parse bra1;
	   {
           mover = 0;				// flag for exchanging spin s1
           for(s2=0; s2<bra2.cols(); s2++)	// look for s2 <-> s1 
	     {
             q1b = bra1.getRe(s1);		// Mz of spin s1 in bra1
             q1k = ket1.getRe(s1);		// Mz of spin s1 in ket1
             q2b = bra2.getRe(s2);		// Mz of spin s2 in bra2
             q2k = ket2.getRe(s2);		// Mz of spin s2 in ket2
	     if(Pro.mapped(ic1,s1,ic2,s2))	// if partners in exchange
               {
               mover = 1;			// s1 in exchange with s2
               if((q1b!=q2b) || (q1k!=q2k))	// different quantum state
                 {
		 match = 0;			// do not click
                 break;
                 }
               }
             }
           if(!match) break;
           if(!mover)			    // s1 is not exchanging
             {
             if(q1b!=q1k)                  // it's not |m><m|
               {
               match = 0;
               break;
               }
             dsys = mixture.Comp(ic1);	// spin that is in superbra
             spinHS = dsys.HS(s1);		// but not in superket
             norm *= spinHS;			// contributes into norm 
             } 
           }
	 if(!match) continue;
         for(s2=0; s2<bra2.cols(); s2++)	// now parse bra2 to find
           {					// local non-movers
           q2b = bra2.getRe(s2);		// Mz of spin s2 in bra2
           q2k = ket2.getRe(s2);		// Mz of spin s2 in ket2

           mover = 0;				// flag for exchanging s2
           for(s1=0; s1<bra1.cols(); s1++)     
	     if(Pro.mapped(ic2, s2, ic1, s1))	  
	       mover = 1;			// this is mover
             if(!mover)				// s2 is not exchanging
               if(q2b!=q2k)			// it is not |m><m|
		 {
                 match = 0;
                 break;
                 }
           }			
 
         if(!match) continue;
         else
           {
	   kex2 /= double(norm);
           Xel.put(-kex2,j,k);			// Normalization good for any 
           }					// non-mutual exchange. 
         }
       }
     Xmx += Xel;				// Add this exchange process'
     }						// contribution to Xnm
   return super_op(Xmx,Dbs);
   }









/*************************************************************************
**									**
**			      Function Xm				**
**									**
** The purpose of this function is to calculate a superoperator that	**
** represents all mutual exchanges within all components   		**
** in a given multi_sys spin system.  The input spin system, msys, 	**
** contains a number of components, each (possibly) with a set of 	**
** defined mutual exchange processes, each of which contains the spins 	**
** exchanging and their exchange rates					**
**									**
** Jacco van Beek, 13-03-2009						**
**									**
*************************************************************************/


super_op Xm(const multi_sys& msys)
  {
  super_op LOp;

  int nc = msys.NComps();			// Number of sub-spaces
  int cmp, ls=0;				// Indicies for sub-spaces
  matrix I;
  matrix *mxc, *bsc;					// Array Liouv. sub-matrices
  int *ncd;					// Array Liouv. sub-space dims
  mxc = new matrix[nc];
  bsc = new matrix[nc];
  ncd = new int[nc];
  gen_op H;
  sys_dynamic sysd;
  super_op Lex;

// 		This Is Composite Liouville Space Superoperator
//	      (Now We Must Deal With All Sub-Spaces Individually)
  for(cmp=0; cmp<nc; cmp++)			// Loop subspaces (components)
    {
    sysd = msys.Comp(cmp); 			//extract the current component
    H = Fx(sysd); 				//generate an arbitrary gen_op for its basis
    Lex = Kex(sysd, H.get_basis()); 		//mutual exchange operator in default Liouville basis
    
    mxc[cmp] = Lex.get_mx();			//store mutual exchange matrix, this compon.
    bsc[cmp] = H.get_basis().get_mx();		//store Hilbert space basis matrix, this compon.
    ncd[cmp] = mxc[cmp].rows();			//Store subspace dimension
    ls += ncd[cmp];				//Track total size of Liouville space
    }

  LOp = super_op(mxc, nc, bsc);
  delete [] mxc;
  delete [] bsc;
  delete [] ncd;
  return LOp;
  }

#endif							// MultiExch.cc
