/* SpinOpCmp.cc *************************************************-*-c++-*-
**									**
**                           G A M M A					**
**									**
**      Composite Spin Operators                   Implementation	**
**									**
**      Copyright (c) 1994						**
**      S.A. Smith							**
**      National High Magnetic Field Laboratory				**
**      1800 E. Paul Dirac Drive					**
**      Tallahassee Florida, 32306-4005					**
**									**
**      Z.L. Madi, T. Levante						**
**      Eidgenoessische Technische Hochschule				**
**      Labor fur physikalische Chemie					**
**      8092 Zurich / Switzerland					**
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
**  Description								**
**									**
** Module comp_spin_op contains functions which produce spin operators	**
** for a composite spin system.  These will span a spin Hilbert space,	**
** defined defined by a list of spins specified by some spin system 	**
** (class spin_sys).							**
**									**
** The composite space matrix representation is in the natural	basis	**
** constructed from the direct product of single spin bases.		**
**								 	**
**  See class spin_sys for spin system function and structure.		**
**  See module single_spin_op for single spin operators.		**
**  See class spin_op (near death) for spin operator workings.		**
**  See class gen_op (future) if spin_op has been removed!		**
**								 	**
*************************************************************************/

#ifndef _SpinOpCmp_cc_			// Is file already included?
#define _SpinOpCmp_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

#include <HSLib/SpinOpCmp.h>		// Include the header file
#include <Basics/Gutils.h>		// Include GAMMA error handling
#include <HSLib/SpinOp.h>		// Knowledge of spin operators
#include <HSLib/SpinOpSng.h>		// Knowledge of single spin operators 

// ____________________________________________________________________________
// i              COMPOSITE SPIN OPERATORS ERROR FUNCTIONS
// ____________________________________________________________________________

        // Input 		eidx    : Error index
        //                      noret   : Flag for return (0=linefeed)
        // Output               none    : Error message
        //				  Stops Execution (fatal)
 
void SOpCmperror(int eidx, int noret)
  {
  std::string hdr("Composite Spin Operator");
  std::string msg;
  switch(eidx)
    {
    case 0: GAMMAerror(hdr, 0, noret); break;   // Program Aborting        (0)
    case 20:msg = std::string("Polarization Operator Ia Undefined for Spins I>1/2");
            GAMMAerror(hdr,msg,noret);	break;	//			   (20)				 
    case 21:msg = std::string("Polarization Operator Ib Undefined for Spins I>1/2");
            GAMMAerror(hdr,msg,noret);	break;	//			   (21)				 
    case 22: msg = std::string("Use Polarization Op. Function Ipol for Spins I>1/2");
            GAMMAerror(hdr,msg,noret);	break;	//			   (22)				 
    case 23: msg = std::string("Bad Quantum Value for Polarization Operator");
            GAMMAerror(hdr,msg,noret);	break;	//			   (23)				 
    case 24: msg = std::string("Too Many Characters in Fpdt Operator Name");
            GAMMAerror(hdr,msg,noret);	break;	//			   (24)				 
    case 25: msg = std::string("Truncating Fpdt Name To Proper Length");
            GAMMAerror(hdr,msg,noret);	break;	//			   (25)				 
    case 26: msg = std::string("Unknown Character in Fpdt Operator Name");
            GAMMAerror(hdr,msg,noret);	break;	//			   (26)				 
    case 27: msg = std::string("Setting Unknown Character in Fpdt Name To 0");
            GAMMAerror(hdr,msg,noret);	break;	//			   (27)				 
    default:GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error   (-1)
    }
  }

void volatile SOpCmpfatal(int eidx)
  {
  SOpCmperror(eidx, 1);
  if(eidx) SOpCmperror(0);
  GAMMAfatal();					// Clean exit from program
  }
     

// ____________________________________________________________________________
// A            SINGLE SPIN SPIN OPERATORS, COMPOSITE SYSTEMS
// ____________________________________________________________________________

/* The spin operators returned here reside in the composite Hilbert space
   of the spin system given in the function call even though they operate
   only on a single spin!  For example, Ix for a single proton would be
   represented by a 2x2 array whereas that same protons Ix when part of
   a three proton system would by 8x8 where 8 is the system Hilbert space.

		Input		sys	: A spin system
				spin	: A spin index [0, ns)
		Output		SOp	: A spin operator in the full
					  Hilbert space of sys               */

spin_op Iu(const spin_sys& sys, int spin, int type)
  {
  int hs, ns=sys.spins();		// Number of system spins
  matrix* spr = new matrix[ns];		// Array of sub-space Ie's
  for(int i=0; i<ns; i++)		// Loop over all spins
    {
    hs = sys.HS(i);			//	Spin i Hilbert space
    if(spin == i)			//	If this is the active spin
      {					//	we must make proper sub-matirx
      switch(type)
        {
        default:
        case 0: spr[i]=Ie(hs); break;	//	For Ie operator
        case 1: spr[i]=Ix(hs); break;	//	For Ix operator
        case 2: spr[i]=Iy(hs); break;	//	For Iy operator
        case 3: spr[i]=Iz(hs); break;	//	For Iz operator
        case 4: spr[i]=Ip(hs); break;	//	For I+ operator
        case 5: spr[i]=Im(hs); break;	//	For I- operator
        case 6:				//	For Ialpha operator
          if(hs != 2) 			// 	Error: Ia on I>1/2 spin
            {				//	Should use Ipol instead
            SOpCmperror(20, 1);		//      Cannot have Ia on I>1/2
            SOpCmpfatal(22);		//	Died on polarization Op
            }
          spr[i]=matrix(hs,hs,0,d_matrix_type);
          (spr[i]).put(1.0, 0,0);
          break;
        case 7:				//	For Ibeta operator
          if(hs != 2) 			// 	Error: Ib on I>1/2 spin
            {				//	Should use Ipol instead
            SOpCmperror(21, 1);		//      Cannot have Ib on I>1/2
            SOpCmpfatal(22);		//	Died on polarization Op
            }
          spr[i]=matrix(hs,hs,0,d_matrix_type);
          (spr[i]).put(1.0, 1,1);
          break;
        }
      }
    else spr[i] = Ie(hs);		// If not active spin, just Ie
    }
  spin_op SOp(ns, spr);
  if(spr) delete [] spr;
  return SOp;
  }

spin_op Ie(const spin_sys& sys, int spin) { return Iu(sys, spin, 0); }
spin_op Ix(const spin_sys& sys, int spin) { return Iu(sys, spin, 1); }
spin_op Iy(const spin_sys& sys, int spin) { return Iu(sys, spin, 2); }
spin_op Iz(const spin_sys& sys, int spin) { return Iu(sys, spin, 3); }
spin_op Ip(const spin_sys& sys, int spin) { return Iu(sys, spin, 4); }
spin_op Im(const spin_sys& sys, int spin) { return Iu(sys, spin, 5); }
spin_op Ia(const spin_sys& sys, int spin) { return Iu(sys, spin, 6); }
spin_op Ib(const spin_sys& sys, int spin) { return Iu(sys, spin, 7); }

spin_op Ipol(const spin_sys& sys, double m, int spin)
  {
  int hs;
  double I;
  int ns = sys.spins();			// Number of system spins
  matrix* spr = new matrix[ns];		// Array of sub-space Ie's
  for (int j=0; j<ns; j++)
    {
    hs = sys.HS(j);			// Get spin Hilbert space
    if(spin == j)			// If this is the active spin
      {					// must get Ipol for spin 
      spr[j]= matrix(hs, hs,0,d_matrix_type);
      I = sys.qn(j); 			//   Spin's Iz value
      if(fabs(m)>I) SOpCmpfatal(23);	//   Bad spin quantum number
      int i = int(I-m);			//   Index for polarization
      (spr[j]).put(1, i,i);		//   Set I[m] polarization
      }
    else spr[j] = Ie(hs);		// For non-active spin, use Ie
    }     
  spin_op SOp(ns,spr);			// Consruct spin operator
  if(spr) delete [] spr;		// Delete local subspace mxs
  return SOp;				// Return the spin op
  }    

/* The Fu Functions Below Just Return The Iu Resuts But Are Added Into The
   Module For Convenience (so either I or F names will work for 1 spin)      */

spin_op Fe(const spin_sys& sys, int spin) { return spin_op(Ie(sys, spin)); }
spin_op Fx(const spin_sys& sys, int spin) { return spin_op(Ix(sys, spin)); }
spin_op Fy(const spin_sys& sys, int spin) { return spin_op(Iy(sys, spin)); }
spin_op Fz(const spin_sys& sys, int spin) { return spin_op(Iz(sys, spin)); }
spin_op Fp(const spin_sys& sys, int spin) { return spin_op(Ip(sys, spin)); }
spin_op Fm(const spin_sys& sys, int spin) { return spin_op(Im(sys, spin)); }
spin_op Fa(const spin_sys& sys, int spin) { return spin_op(Ia(sys, spin)); }
spin_op Fb(const spin_sys& sys, int spin) { return spin_op(Ib(sys, spin)); }
spin_op Fpol(const spin_sys& sys, double m, int spin)
                                          { return spin_op(Ipol(sys,m,spin));}

// ____________________________________________________________________________
// B                  SPIN OPERATORS OVER THE SPIN SYSTEM
// ____________________________________________________________________________

/* The spin operators returned by the functions below reside in the composite
   Hilbert space of the spin system given in the function call and they operate
   on every spin in that system.  That is to say, that
                        
                  spins                          type
                   ---
             F  =  \     F          e ---> Identity    p ---> F+
              u    /      ui        x ---> Fx          m ---> F-
                   ---              y ---> Fy          a ---> Falpha
                    i               z ---> Fz          b ---> Fbeta
                                    pol ---> Fpol

		Input		sys	: A spin system
				type	: An operator type 
		Output		SOp	: A spin operator in the full
					  Hilbert space of sys               */

spin_op Fe(const spin_sys & sys) { return Faxis(sys, 'e'); }
spin_op Fx(const spin_sys & sys) { return Faxis(sys, 'x'); }
spin_op Fy(const spin_sys & sys) { return Faxis(sys, 'y'); }
spin_op Fz(const spin_sys & sys) { return Faxis(sys, 'z'); }
spin_op Fp(const spin_sys & sys) { return Faxis(sys, 'p'); }
spin_op Fm(const spin_sys & sys) { return Faxis(sys, 'm'); }
spin_op Fa(const spin_sys & sys) { return Faxis(sys, 'a'); }
spin_op Fb(const spin_sys & sys) { return Faxis(sys, 'b'); }
spin_op Fpol(const spin_sys & sys, double m) { return Fpol_gen_new(sys, m); }

// ____________________________________________________________________________
// C               SPIN OPERATORS OVER A SPIN ISOTOPE TYPE
// ____________________________________________________________________________

/* The spin operators returned by the functions below reside in the composite
   Hilbert space of the spin system given in the function call. They operate
   exclusively on spins of a particular isotope type. That is to say, that
                        
           spins                                        type
            ---         [ 1 if iso == iso
      F  =  \     F   * |         i         e ---> Identity    p ---> F+
       u    /      ui   [ 0 if iso != iso   x ---> Fx          m ---> F-
            ---                   i         y ---> Fy          a ---> Falpha
             i                              z ---> Fz          b ---> Fbeta
                                                   pol ---> Fpol

		Input		sys	: A spin system
				type	: An operator type                   */

spin_op Fe(const spin_sys& sys, const std::string& iso)
  { return Faxis(sys, iso, 'e'); }

spin_op Fx(const spin_sys& sys, const std::string& iso)
  { return Faxis(sys, iso, 'x'); }

spin_op Fy(const spin_sys& sys, const std::string& iso)
  { return Faxis(sys, iso, 'y'); }

spin_op Fz(const spin_sys& sys, const std::string& iso)
  { return Faxis(sys, iso, 'z'); }

spin_op Fp(const spin_sys& sys, const std::string& iso)
  { return Faxis(sys, iso, 'p'); }

spin_op Fm(const spin_sys& sys, const std::string& iso)
  { return Faxis(sys, iso, 'm'); }

spin_op Fa(const spin_sys& sys, const std::string& iso)
  { return Faxis(sys, iso, 'a'); }

spin_op Fb(const spin_sys& sys, const std::string& iso)
  { return Faxis(sys, iso, 'b'); }

spin_op Fpol(const spin_sys& sys, double m, const std::string& iso)
  { return Fpol_gen_new(sys, iso, m); }


// ____________________________________________________________________________
// D                 SPIN OPERATORS OVER FLAGGED SPINS 
// ____________________________________________________________________________

/* The spin operators returned here reside in the composite Hilbert space of
   the spin system given in the function and they operate only on spins in
   that system whose spin flags have been set to TRUE (see class spin_sys)
 
           spins                                        type
            ---         [ 1 if flag(i)!=0
      F  =  \     I   * |                   e ---> Identity    p ---> F+
       u    /      ui   [ 0 if flag(i)=0    x ---> Fx          m ---> F-
            ---                             y ---> Fy          a ---> Falpha
             i                              z ---> Fz          b ---> Fbeta
                                                   pol ---> Fpol
   
    The last two functions, Fb_sp & Fpol_sp, take the dummy argument iso for 
    arument input compatibility with other F*_sp functions.  Some compilers
    insist that all arguments are used, these functions do use the string
    iso AFTER the function return, but only to avoid compiler warnings.     */

spin_op Fe(const spin_sys& sys, const flagvec& sflags)
                                           { return Faxis(sys, sflags, 'e'); }
spin_op Fx(const spin_sys& sys, const flagvec& sflags)
                                           { return Faxis(sys, sflags, 'x'); }
spin_op Fy(const spin_sys& sys, const flagvec& sflags)
                                           { return Faxis(sys, sflags, 'y'); }
spin_op Fz(const spin_sys& sys, const flagvec& sflags)
                                           { return Faxis(sys, sflags, 'z'); }
spin_op Fp(const spin_sys& sys, const flagvec& sflags)
                                           { return Faxis(sys, sflags, 'p'); }
spin_op Fm(const spin_sys& sys, const flagvec& sflags)
                                           { return Faxis(sys, sflags, 'm'); }
spin_op Fa(const spin_sys& sys, const flagvec& sflags)
                                           { return Faxis(sys, sflags, 'a'); }
spin_op Fb(const spin_sys& sys, const flagvec& sflags)
                                           { return Faxis(sys, sflags, 'b'); }
spin_op Fpol(const spin_sys& sys, double m, const flagvec& sflags)
                                      { return Fpol_gen_new(sys, sflags, m); }

spin_op Fe_sp(const spin_sys& sys) { return Faxis_sp(sys, 'e'); }
spin_op Fx_sp(const spin_sys& sys) { return Faxis_sp(sys, 'x'); }
spin_op Fy_sp(const spin_sys& sys) { return Faxis_sp(sys, 'y'); }
spin_op Fz_sp(const spin_sys& sys) { return Faxis_sp(sys, 'z'); }
spin_op Fp_sp(const spin_sys& sys) { return Faxis_sp(sys, 'p'); }
spin_op Fm_sp(const spin_sys& sys) { return Faxis_sp(sys, 'm'); }
spin_op Fa_sp(const spin_sys& sys) { return Faxis_sp(sys, 'a'); }
spin_op Fb_sp(const spin_sys& sys) { return Faxis_sp(sys, 'b'); }
spin_op Fpol_sp(const spin_sys& sys, double m) { return Fpol_gen(sys, m); }

// ____________________________________________________________________________
// E                     GENERIC SPIN OPERATOR FUNCTIONS
// ____________________________________________________________________________

/* The spin operators returned here reside in the composite Hilbert space of
   the spin system given in the function.  These are "generic" functions in
   that other functions use them to construct the desired operator and this
   conserves redundant code.                                                 */

// ------------- Associated With A Cartesian Axis Or Specific Type ------------

/* These functions generate spin operators that are a specific type (or a
   specific axis).  The operator type is specified by a single input "type".
   Valid characters and their meaning are as follows:

      character: x      y      z      p      m      e      a       b
      operator:  Ix     Iy     Iz     I+     I-     Ie   Ialpha  Ibeta
 
         OverLoads of Function where Type of Second Argument Changes
 
           OL Function    Type    Label  Result
           -- --------    ----    ----   ------------------------------
           A. Faxis       int     spin   Spin index, only spin affected
           B. Faxis       string  iso    Isotope label, spins of type iso hit
           C. Faxis       ---            All spins in sys affected
           D. Faxis_sp    ---            Spins with flags TRUE affected
           E. Faxis_sp    vector  flags  Spins with flags TRUE affected
 
    Note that some of these spin operators will have a specialized matrix
    structure. Currently this is NOT handled well by class spin_op so the
    structures must be enforced here.  That can be removed when spin_op
    works better with GAMMA matrix classes - sosi.                           */

spin_op Faxis(const spin_sys& sys, int spin, char axis)
  {
  spin_op SOp;
  switch(axis)
    {
    case 'x': SOp = Ix(sys,spin); break;				// Fx
    case 'y': SOp = Iy(sys,spin); break;				// Fy
    case 'z': SOp = Iz(sys,spin); break;			 	// Fz
    case 'e': SOp = Ie(sys,spin); break;				// Fe
    case 'p': SOp = Ip(sys,spin); break;				// F+
    case 'm': SOp = Im(sys,spin); break;				// F-
    case 'a': SOp = Ia(sys,spin); break;				// Falp
    case 'b': SOp = Ib(sys,spin); break;				// Fbet
    }
  return SOp;
  }

spin_op Faxis(const spin_sys & sys, const std::string& isoin, char axis)
  { return Faxis(sys, sys.GetFlags(isoin,1,0), axis); }

spin_op Faxis(const spin_sys & sys, char axis)
  {
  spin_op SOp;
  int ns = sys.spins();
  int spin = 0;
  switch(axis)
    {
    case 'x': for(spin=0; spin<ns; spin++) SOp += Ix(sys,spin); break;	// Fx
    case 'y': for(spin=0; spin<ns; spin++) SOp += Iy(sys,spin); break;	// Fy
    case 'z': for(spin=0; spin<ns; spin++) SOp += Iz(sys,spin); break; 	// Fz
    case 'e': for(spin=0; spin<ns; spin++) SOp += Ie(sys,spin);	break;	// Fe
    case 'p': for(spin=0; spin<ns; spin++) SOp += Ip(sys,spin); break;	// F+
    case 'm': for(spin=0; spin<ns; spin++) SOp += Im(sys,spin); break;	// F-
    case 'a': for(spin=0; spin<ns; spin++) SOp += Ia(sys,spin); break;	// Falp
    case 'b': for(spin=0; spin<ns; spin++) SOp += Ib(sys,spin); break; 	// Fbet
    }
  SOp.FaxisStruct(axis);		// Force Fz,Fe,Fz,&Fb diagonal
  return SOp;
  }

spin_op Faxis_sp(const spin_sys & sys, char axis)
  {
  spin_op SOp;
  int ns = sys.spins();
  int i = 0;
  switch(axis)
    {
    case 'x': for(i=0;i<ns;i++) if(sys.GetFlag(i)) SOp+=Ix(sys,i); break; // Fx
    case 'y': for(i=0;i<ns;i++) if(sys.GetFlag(i)) SOp+=Iy(sys,i); break; // Fy
    case 'z': for(i=0;i<ns;i++) if(sys.GetFlag(i)) SOp+=Iz(sys,i); break; // Fz
    case 'e': for(i=0;i<ns;i++) if(sys.GetFlag(i)) SOp+=Ie(sys,i); break; // Fe
    case 'p': for(i=0;i<ns;i++) if(sys.GetFlag(i)) SOp+=Ip(sys,i); break; // F+
    case 'm': for(i=0;i<ns;i++) if(sys.GetFlag(i)) SOp+=Im(sys,i); break; // F-
    case 'a': for(i=0;i<ns;i++) if(sys.GetFlag(i)) SOp+=Ia(sys,i); break; // Fa
    case 'b': for(i=0;i<ns;i++) if(sys.GetFlag(i)) SOp+=Ib(sys,i); break; // Fb
    }
  SOp.FaxisStruct(axis);		// Force Fz,Fe,Fz,&Fb diagonal
  return SOp;
  }

spin_op Faxis(const spin_sys& sys, const flagvec& flags, char axis)
  {
  spin_op SOp;
  int ns = sys.spins();
  int i=0;
  switch (axis)
    {
    case 'x': for(i=0; i<ns; i++) if(flags[i]) SOp+=Ix(sys,i); break; // Fx
    case 'y': for(i=0; i<ns; i++) if(flags[i]) SOp+=Iy(sys,i); break; // Fy
    case 'z': for(i=0; i<ns; i++) if(flags[i]) SOp+=Iz(sys,i); break; // Fz
    case 'e': for(i=0; i<ns; i++) if(flags[i]) SOp+=Ie(sys,i); break; // Fe
    case 'p': for(i=0; i<ns; i++) if(flags[i]) SOp+=Ip(sys,i); break; // F+
    case 'm': for(i=0; i<ns; i++) if(flags[i]) SOp+=Im(sys,i); break; // F-
    case 'a': for(i=0; i<ns; i++) if(flags[i]) SOp+=Ia(sys,i); break; // Fa
    case 'b': for(i=0; i<ns; i++) if(flags[i]) SOp+=Ib(sys,i); break; // Fb
    }
  SOp.FaxisStruct(axis);		// Force Fz,Fe,Fz,&Fb diagonal
  return SOp;
  }

// --------------------------- Polarization Operators -------------------------

spin_op Fpol_gen_new(const spin_sys & sys, double m)

	// Input		sys  : Spin system (base)
	// 			m    : Spin quantum value
	// Return		SOp  : Spin polarization operator which
	//			       is Fpol for spins flagged
	// Note		             : SOp initialized to zero at creation.
	// Note		             : Includes all spins in sys

  {
  int ns = sys.spins();
  spin_op SOp;
  for(int spin=0; spin<ns; spin++)
    SOp += Ipol(sys, m, spin);
  SOp.FaxisStruct('z');			// Force Fpol diagonal
  return SOp;
  }

spin_op Fpol_gen_new(const spin_sys & sys, const std::string& isoin, double m)

	// Input		sys  : Spin system (base)
	//			isoin: Label for an isotope type
	// 			m    : Spin quantum value
	// Return		SOp  : Spin polarization operator which
	//			       is Fpol over spins of type isoin
	// Note		             : SOp initialized to zero at creation.
	// Note		             : Spins included if of type isoin

  { return Fpol_gen_new(sys, sys.GetFlags(isoin,1,0), m); }


spin_op Fpol_gen_new(const spin_sys& sys, const flagvec& flags, double m)

	// Input		sys  : Spin system (base)
	//			flags: Flags for which spins to include
	// 			m    : Spin quantum value
	// Return		SOp  : Spin polarization operator which
	//			       is Fpol for spins flagged
	// Note		             : SOp initialized to zero at creation.
	// Note		             : Spins included if flags = TRUE

  {
  int ns = sys.spins();
  spin_op SOp;
  for(int spin=0; spin<ns; spin++)
    if(flags[spin]) SOp += Ipol(sys,m,spin);
  SOp.FaxisStruct('z');				// Force Fpol diagonal
  return SOp;
  }


spin_op Fpol_gen(const spin_sys & sys, double m)

	// Input		sys  : Spin system (base)
	// 			m    : Spin quantum value
	// Return		SOp  : Spin polarization operator which
	//			       is Fpol for spins flagged
	// Note		             : SOp initialized to zero at creation.
	// Note		             : Spins included if flags = TRUE

  {
  int ns = sys.spins();
  spin_op SOp;
  for (int spin=0; spin<ns; spin++)
    if (sys.GetFlag(spin)) SOp += Ipol(sys,m,spin);
  SOp.FaxisStruct('z');			// Force Fpol diagonal
  return SOp;
  }

// ------------------------------ Product Operators ---------------------------


spin_op Ipdt(const spin_sys &sys, std::string name)

	// Input		sys  : Spin system (base)
	// 			name : Spin operator name, contains a
	//			       single character per spin in
	//			       the system.
	// Return		SOp  : Spin operator product requested
	// Note		             : Characters allowed in name

	// 				x - Ix		e - identity
	// 				y - Iy		p - I+
	// 				z - Iz		m - I-
	// 				0 - identity	a - Ialpha
	// 				1 - Ix		b - Ibeta
	// 				2 - Iy		+ - I+
	// 				3 - Iz		- - I-

	// Note		             : Name examples are -
	//			       "+-x" = I+(0)I(1)Ix(2)
	//			       "y0z" = Iy(0)Iz(2)

  {
  spin_op Iproduct = Ie(sys,0);		// Initialize base operator to unity
  int ns = sys.spins();			// Number of spins in sys
  int nOps = name.length();		// Number of characters in nOps
  if(nOps > ns)
    {
    SOpCmperror(24, 1);			//   Too many characters in name
    SOpCmperror(25);			//   Truncating Fpdt name
    nOps = ns;				//   Only use ns characters of nOps
    }
  char letter;
  for(int spin=0; spin<nOps; spin++)	// Loop through all operators and
    {					// form the appropriate product
    letter = name[spin];
    if(letter=='x' || letter== '1')     Iproduct *= Ix(sys,spin);
    else if(letter=='y' || letter=='2') Iproduct *= Iy(sys,spin);
    else if(letter=='z' || letter=='3') Iproduct *= Iz(sys,spin);
    else if(letter=='p' || letter=='+') Iproduct *= Ip(sys,spin);
    else if(letter=='m' || letter=='-') Iproduct *= Im(sys,spin);
    else if(letter=='a')                Iproduct *= Ia(sys,spin);
    else if(letter=='b')                Iproduct *= Ib(sys,spin);
    else if(letter!='e' && letter!='0')
      {
      SOpCmperror(26, 1);		//   Unknown character in name
      SOpCmperror(27);			//   Setting character to 0
      }
    }
  return Iproduct;
  }


spin_op Fpdt(const spin_sys &sys, std::string name) { return Ipdt(sys, name); }


#endif 							// SpinOpCmp.cc
