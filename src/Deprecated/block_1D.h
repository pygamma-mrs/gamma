/*************************************************************************
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      Block_1D                                      Interface		**
**                                                                      **
**      Copyright (c) 1990, 1991, 1992                                  **
**      Scott Smith                                                     **
**      Eidgenoessische Technische Hochschule                           **
**      Labor fuer physikalische Chemie                                 **
**      8092 Zurich / Switzerland                                       **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/


/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** This class specifies the structure and properties of a 1-dimensional **
** block of simulated (or any other) data.  The data block is a row     **
** vector in RAM memory and has an associated parameter set.  Functions **
** are provided for basic algebraic manipulations, specification of     **
** parameters, and general I/O.  Many external programs deal with data  **
** sets coupled to a block of parameters and this class allows users    **
** to more readily support I/O to these programs.  One example would be **
** the Felix module which provides functions for input/output of GAMMA  **
** data to/from the program Felix from Molecular Simulations (MSI),     **
** nee Biosym, nee Hare Research, Inc.                                  **
**                                                                      **
*************************************************************************/

#ifndef   Gblock_1D_h_			// Is this file already included?
#  define Gblock_1D_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/row_vector.h>		// Include base class row_vector 
#include <Basics/ParamSet.h>		// Include base class ParameterSet 

class block_1D : public row_vector, public ParameterSet
  {

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                    CLASS BLOCK_1D ERROR HANDLING
// ____________________________________________________________________________


void Blk1DError(int eidx, int noret=0) const;

	// Input                BLK	: A 1D data block
        //                      eidx    : Error index
        //                      noret   : Flag for linefeed (0=linefeed)
	// Output               none    : Error message output

 

void Blk1DError(int eidx, int i, int j, int noret=0) const;

        // Input                BLK     : A 1D data block
        //                      eidx    : Error index   
	//			i,j	: Block indices
        //                      noret   : Flag for linefeed (0=linefeed)
	// Output               none    : Error message output


volatile void Blk1DFatal(int eidx);

        // Input                BLK     : A 1D data block
        //                      eidx    : Error index   
	// Output		none : Stops Execution & Error Message Output


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A                       BLOCK_1D CONSTRUCTORS / DESTRUCTOR
// ____________________________________________________________________________
	

MSVCDLC block_1D ( ) {};

	// Input			none : 
	// Output			none : creates an empty block_1D
        // Note                              : only pointers created,
        //                                     no array space allocated


MSVCDLC block_1D (int dim);

	// Input			dim  : integer
	// Output			none : creates an empty block_1D
	//				       of dimension specified by dim


MSVCDLC block_1D (const block_1D& BLK):row_vector(BLK), ParameterSet(BLK) {};

	// Input			BLK  : block_1D
	// Output			     : returns a new block_1D equivalent
	//				       to input block_1D BLK


MSVCDLC block_1D (const row_vector& vec):row_vector(vec) {};

	// Input			vec  : A row vector
	// Output			BLK  : Returns a new block_1D with
	//				       data equal to input block_1D BLK

MSVCDLC ~block_1D ( ) {};

	// Input			none : 
	// Output			none : destroys block_1D



// ______________________________________________________________________
//                CLASS BLOCK_1D FUNCTIONS, BLOCK WITH BLOCK
// ______________________________________________________________________


MSVCDLL friend block_1D operator + (const block_1D& BLK1, const block_1D& BLK2);

	// Input		BLK1, BLK2 : two block_1Ds
	// Return		BLK  : block_1D which is the addition
        //	                       of the two input block_1Ds,
	//			       this = BLK =  BLK1 + BLK2


MSVCDLL void operator += (const block_1D& BLK1);

	// Input		BLK1 : a block_1D
	// Return		BLK  : block_1D which is the input
        //	                       input block_1D BLK1 added to itself
	//			       BLK = BLK + BLK1


MSVCDLL friend block_1D operator - (const block_1D& BLK1, const block_1D& BLK2);

	// Input		BLK1 : a data block_1D
	//			BLK2 : a data block_1D
	// Return		BLK  : block_1D which is the subtraction
        //	                       of the two input block_1Ds,
	//			       BLK = BLK1 - BLK2


MSVCDLL friend  block_1D operator - (const block_1D& BLK1);

	// Input		BLK1  : a block_1D
	// Return		BLK: block_1D which is the negative
        //	                     of the input block_1D, BLK = - BLK1


MSVCDLL void operator -= (const block_1D& BLK1 );

	// Input		BLK1    : a block_1D
	// Return		BLK  : block_1D which is the input
        //	                       input block_1D BLK1 added to itself
	//			       BLK = BLK - BLK1


MSVCDLL void operator = (const block_1D& BLK);

	// Input		this,BLK  : two block_1D
	// Return		copies BLK to this, and destroys the old
        //                      contest of this block_1D


MSVCDLL void operator = (const matrix& mx);

	// Input		this,BLK  : two block_1D
	// Return		copies BLK to this, and destroys the old
        //                      contest of this block_1D



// ********** Basic Function Definitions: Block_1D with Scalar **********


MSVCDLL friend block_1D operator * (const block_1D& BLK1, const complex& z);
MSVCDLL friend block_1D operator * (const block_1D& BLK1, double d);

	// Input		BLK1  : A block_1D
        //                      z/d   : A complex number/double
	// Return		BLK: block_1D which is the multlipication
        //	                     of the block_1D by the complex number
        //                           BLK = z * BLK1 or BLK = d * BLK


MSVCDLL friend block_1D operator * (const complex& z, const block_1D& BLK1);
MSVCDLL friend block_1D operator * (double d, const block_1D& BLK1);

        // Input                z  : a complex number
	//      		BLK1  : a block_1D
	// Return		BLK: block_1D which is the multlipication
        //	                     of the block_1D by the complex number
        //                           BLK = z * BLK1


MSVCDLL friend block_1D operator / (const block_1D& BLK1, const complex& z);

	// Input		BLK1  : a block_1D
        //                      z  : a complex number
	// Return		BLK: block_1D which is the multlipication
        //	                     of the block_1D by the inverted complex
	//			     number. BLK = (1/z) * BLK1


// ******* Basic Function Definitions: Block_1D with Parameter Set ******


//    friend void operator += (block_1D& BLK, const ParameterSet& pset);

	// Input		BLK  : a block_1D
        //                      pset : a parameter set
	// Return		BLK  : block_1D has the paramter set pset added


// ______________________________________________________________________
//                           DATA SET FUNCTIONS
// ______________________________________________________________________


MSVCDLL friend void product(block_1D &BLK1, const block_1D &BLK2);

	// Input		BLK1 : data block
	// 			BLK2 : data block
	// Return		none : BLK1 multiplied point by point by BLK2


// ______________________________________________________________________
//                         PARAMETER SET FUNCTIONS
// ______________________________________________________________________


MSVCDLL void set_defaults(int level);

	// Input		BLK  : data block (this)
	//			level: level to which defaults are set
	// Return		none : BLK parameter set defaults set


// ______________________________________________________________________
//                 CLASS BLOCK_1D CHECKS AND ERROR MESSAGES
// ______________________________________________________________________


MSVCDLL int Blk1DCheck(int pt) const;

	// Input		BLK  : block_1D
	// 			pt   : point (integer)
	// Output		int  : insures that point pt is a valid
        //			       point for a block_1D. 0 if bad


MSVCDLL int Blk1DCheck(const block_1D &BLK2) const;

	// Input		BLK1    : block_1D
	// 	 		BLK2    : block_1D
	// Output		int     : insures that the two block_1Ds
	//			          can computationally intermix
	//			          0 if incompatible pair


// ____________________________________________________________________________
//                       CLASS BLOCK_1D I/O FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                         I/O On Block Parameter Set
// ----------------------------------------------------------------------------

MSVCDLL void write_pset(const std::string& filename) const;
 
        // Input                BLK  : A 1D data block (this)   
        //                      file : Name of disk file to produce
        // Return               none : A file is produced on disk containing
        //                             the parameter set in BLK
 
 
MSVCDLL friend std::ostream& operator<< (std::ostream& ostr, const block_1D& BLK);

	// Input		ostr : string, error message
	// 			BLK  : block_1D
	// Return		     : stream, prints block_1D


MSVCDLL friend void print_pset(const block_1D& BLK);

	// Input		BLK  : a block_1D
	// Return		none : BLK parameters printed to std I/O


MSVCDLL friend void print_dset(block_1D& BLK);
// sosi - make constant when vector classes and matrices are constant

	// Input		BLK  : a block_1D
	// Return		none : BLK data vector printed to std I/O

  };

extern void product(block_1D &, const block_1D &);

#endif 						/* __CLASS_BLOCK_1D_H__ */

