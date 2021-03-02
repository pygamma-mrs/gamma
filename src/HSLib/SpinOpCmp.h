/* SpinOpCmp.h **************************************************-*-c++-*-
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      Composite Spin Operators                   Interface		**
**                                                      `               **
**      Copyright (c) 1999                                              **
**      S.A. Smith							**
**      National High Magnetic Field Laboratory				**
**      1800 E. Paul Dirac Drive					**
**      Tallahassee Florida, 32306-4005					**
**									**
**      Z.L. Madi, T. Levante						**
**      Eidgenoessische Technische Hochschule				**
**      Labor fur physikalische Chemie					**
**      8092 Zurich / Switzerland					**
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** Module SpinOpCmp contains functions which produce spin operators 	**
** for a composite spin system.  These will span a spin Hilbert space,	**
** defined by a list of spins specified by some spin system (see class	**
** spin_sys).								**
**									**
**  The composite space matrix representation is in the natural	basis	**
**  constructed from the direct product of single spin bases.		**
**								 	**
**  See class spin_sys for spin system function and structure.		**
**  See module SpinOpSng for single spin operators.			**
**  See class SpinOp (near death) for spin operator workings.		**
**                                                                      **
*************************************************************************/

#ifndef SpinOpCmp_h_			// Is file already included?
#  define SpinOpCmp_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <HSLib/SpinOp.h>		// Include spin operators
#include <HSLib/SpinSys.h>		// Include base spin systems
#include <string>			// Include stdlibc++ strings

// ____________________________________________________________________________
// i              COMPOSITE SPIN OPERATORS ERROR FUNCTIONS
// ____________________________________________________________________________

        // Input                eidx    : Error index
        //                      noret   : Flag for return (0=linefeed)
        // Output               none    : Error message
        //				  Stops Execution (fatal)

void          SOpCmperror(int eidx, int noret=0);
void volatile SOpCmpfatal(int eidx);

// ____________________________________________________________________________
// A              SINGLE SPIN SPIN OPERATORS, COMPOSITE SYSTEMS
// ____________________________________________________________________________
 
/* The spin operators returned here reside in the composite Hilbert space
   of the spin system given in the function call even though they operate
   only on a single spin!  For example, Ix for a single proton would be
   represented by a 2x2 array whereas that same protons Ix when part of
   a three proton system would by 8x8 where 8 is the system Hilbert space.
 
                Input           sys     : A spin system
                                spin    : A spin index [0, ns)
                Output          SOp     : A spin operator in the full
                                          Hilbert space of sys               */  
 
MSVCDLL spin_op Iu(const   spin_sys &sys, int spin, int type);
MSVCDLL spin_op Ie(const   spin_sys &sys, int spin);
MSVCDLL spin_op Iz(const   spin_sys &sys, int spin);
MSVCDLL spin_op Ix(const   spin_sys &sys, int spin);
MSVCDLL spin_op Iy(const   spin_sys &sys, int spin);
MSVCDLL spin_op Ip(const   spin_sys &sys, int spin);
MSVCDLL spin_op Im(const   spin_sys &sys, int spin);
MSVCDLL spin_op Ia(const   spin_sys &sys, int spin);
MSVCDLL spin_op Ib(const   spin_sys &sys, int spin);
MSVCDLL spin_op Ipol(const spin_sys &sys, double m, int spin);

/* The Fu Functions Below Just Return The Iu Resuts But Are Added Into The
   Module For Convenience (so either I or F names will work for 1 spin)      */  
 
MSVCDLL spin_op Fe(const spin_sys &sys, int spin);
MSVCDLL spin_op Fx(const spin_sys &sys, int spin);
MSVCDLL spin_op Fy(const spin_sys &sys, int spin);
MSVCDLL spin_op Fz(const spin_sys &sys, int spin);
MSVCDLL spin_op Fp(const spin_sys &sys, int spin);
MSVCDLL spin_op Fm(const spin_sys &sys, int spin);
MSVCDLL spin_op Fa(const spin_sys &sys, int spin);
MSVCDLL spin_op Fb(const spin_sys &sys, int spin);
MSVCDLL spin_op Fpol(const spin_sys &sys, double m, int spin);

// ____________________________________________________________________________
// B                SPIN OPERATORS OVER THE SPIN SYSTEM
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
 
                Input           sys     : A spin system
                                type    : An operator type
                Output          SOp     : A spin operator in the full
                                          Hilbert space of sys               */

MSVCDLL spin_op Fe(const spin_sys &sys);
MSVCDLL spin_op Fx(const spin_sys &sys);
MSVCDLL spin_op Fy(const spin_sys &sys);
MSVCDLL spin_op Fz(const spin_sys &sys);
MSVCDLL spin_op Fp(const spin_sys &sys);
MSVCDLL spin_op Fm(const spin_sys &sys);
MSVCDLL spin_op Fa(const spin_sys &sys);
MSVCDLL spin_op Fb(const spin_sys &sys);
MSVCDLL spin_op Fpol(const spin_sys &sys, double m);

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
 
                Input           sys     : A spin system
                                type    : An operator type                   */

MSVCDLL spin_op Fe(const spin_sys &sys, const std::string& iso);
MSVCDLL spin_op Fx(const spin_sys &sys, const std::string& iso);
MSVCDLL spin_op Fy(const spin_sys &sys, const std::string& iso);
MSVCDLL spin_op Fz(const spin_sys &sys, const std::string& iso);
MSVCDLL spin_op Fp(const spin_sys &sys, const std::string& iso);
MSVCDLL spin_op Fm(const spin_sys &sys, const std::string& iso);
MSVCDLL spin_op Fa(const spin_sys &sys, const std::string& iso);
MSVCDLL spin_op Fb(const spin_sys &sys, const std::string& iso);
MSVCDLL spin_op Fpol(const spin_sys &sys, double m, const std::string& iso);

// ____________________________________________________________________________
// D                      SPIN OPERATORS OVER FLAGGED SPINS
// ____________________________________________________________________________

/* The spin operators returned here reside in the composite Hilbert space of
   the spin system given in the function and they operate only on spins in
   that system whose spin flags have been set to TRUE (see class spin_sys)

           spins                                        type
            ---         [ 1 if flag(i)!=0
      F  =  \     I   * |         i         e ---> Identity    p ---> F+
       u    /      ui   [ 0 if flag(i)=0    x ---> Fx          m ---> F-
            ---                   i         y ---> Fy          a ---> Falpha
             i                              z ---> Fz          b ---> Fbeta
                                                   pol ---> Fpol
 
    The last two functions, Fb_sp & Fpol_sp, take the dummy argument iso for
    arument input compatibility with other F*_sp functions.  Some compilers
    insist that all arguments are used, these functions do use the string
    iso AFTER the function return, but only to avoid compiler warnings.     */
 
MSVCDLL spin_op Fe(const spin_sys& sys, const flagvec& sflags);
MSVCDLL spin_op Fx(const spin_sys& sys, const flagvec& sflags);
MSVCDLL spin_op Fy(const spin_sys& sys, const flagvec& sflags);
MSVCDLL spin_op Fz(const spin_sys& sys, const flagvec& sflags);
MSVCDLL spin_op Fp(const spin_sys& sys, const flagvec& sflags);
MSVCDLL spin_op Fm(const spin_sys& sys, const flagvec& sflags);
MSVCDLL spin_op Fa(const spin_sys& sys, const flagvec& sflags);
MSVCDLL spin_op Fb(const spin_sys& sys, const flagvec& sflags);
MSVCDLL spin_op Fpol(const spin_sys& sys, double m, const flagvec& sflags);
 
MSVCDLL spin_op Fe_sp(const spin_sys& sys);
MSVCDLL spin_op Fx_sp(const spin_sys& sys);
MSVCDLL spin_op Fy_sp(const spin_sys& sys);
MSVCDLL spin_op Fz_sp(const spin_sys& sys);
MSVCDLL spin_op Fp_sp(const spin_sys& sys);
MSVCDLL spin_op Fm_sp(const spin_sys& sys);
MSVCDLL spin_op Fa_sp(const spin_sys& sys);
MSVCDLL spin_op Fb_sp(const spin_sys& sys);
MSVCDLL spin_op Fpol_sp(const spin_sys& sys, double m);

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
           B. Faxis       string  I      Isotope label, spins of type iso hit
           C. Faxis       ---            All spins in sys affected
           D. Faxis_sp    ---            Spins with flags TRUE affected
           E. Faxis_sp    int*    flags  Spins with flags TRUE affected      */

MSVCDLL spin_op Faxis(const    spin_sys& sys, int spin,             char axis);
MSVCDLL spin_op Faxis(const    spin_sys& sys, const std::string& I, char axis);
MSVCDLL spin_op Faxis_sp(const spin_sys &sys,                       char axis);
MSVCDLL spin_op Faxis(const    spin_sys &sys,                       char axis);
MSVCDLL spin_op Faxis(const    spin_sys &sys, const flagvec& flags, char axis);

// --------------------------- Polarization Operators -------------------------

MSVCDLL spin_op Fpol_gen_new(const spin_sys& sys,                      double m);
MSVCDLL spin_op Fpol_gen_new(const spin_sys& sys,const std::string& I, double m);
MSVCDLL spin_op Fpol_gen_new(const spin_sys& sys,const flagvec& flags, double m);
MSVCDLL spin_op Fpol_gen(const     spin_sys& sys,                      double m);

// ------------------------------ Product Operators ---------------------------

MSVCDLL spin_op Ipdt(const spin_sys &sys, std::string name);
MSVCDLL spin_op Fpdt(const spin_sys &sys, std::string name);

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

#endif							// SpinOpCmp.h
