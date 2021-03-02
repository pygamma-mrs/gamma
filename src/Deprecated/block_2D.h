/* block_2D.h ***********************************-*-c++-*-
**							**
**                     G A M M A 			**
**							**
**	Block_2D 	      Interface definition	**
**							**
**	Copyright (c) 1990, 1991, 1992			**
**	Scott Smith					**
**	Eidgenoessische Technische Hochschule		**
**	Labor fuer physikalische Chemie			**
**	8092 Zuerich / Switzerland			**
**							**
**      $Header: $
**							**
**      Modifications:					**
**							**
*********************************************************/

/*********************************************************
**							**
** 	Description					**
**							**
**	The class block_2D defines the attributes	**
**	of a 2-Dimensional data block quantity,		**
**	assoicated functions, & I/O routines		**
**							**
*********************************************************/

#ifndef   Gblock_2D_h_		// Is file already included?
#  define Gblock_2D_h_ 1	// If no, then remember it
#  if defined(GAMPRAGMA)	// Using the GNU compiler?
#    pragma interface		// this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/matrix.h>
#include <Basics/ParamSet.h>

class block_2D : public matrix, public ParameterSet
  {

// ----------------------------------------------------------------------
// ------------------------ PRIVATE FUNCTIONS ---------------------------
// ----------------------------------------------------------------------

// ______________________________________________________________________
//                   CLASS BLOCK_2D ERROR HANDLING
// ______________________________________________________________________


  friend void block_2D_error(int error);

	// Input		error: Error Flag
	// Output		none : Error Message Output
 

  friend void block_2D_error(int error, int i, int j);

	// Input		error: Error Flag
	// Output		none : Error Message Output


 friend volatile void block_2D_fatality(int error);

	// Input		error: Error Flag
	//			i,j  : Labels
	// Output		none : Stops Execution & Error Message Output


// ----------------------------------------------------------------------
// ------------------------- PUBLIC FUNCTIONS ---------------------------
// ----------------------------------------------------------------------

  public:

// ______________________________________________________________________
//                      BLOCK_2D CONSTRUCTORS / DESTRUCTOR
// ______________________________________________________________________
	

MSVCDLC block_2D ( ) {};

	// Input			none : 
	// Output			none : creates an empty block_2D
        // Note                              : only pointers created,
        //                                     no array space allocated


MSVCDLC block_2D (int dim1, int dim2);

	// Input			dim  : integer
	// Output			none : creates an empty block_2D
	//				       of dimension specified by dim


MSVCDLC block_2D (const block_2D& BLK):matrix(BLK), ParameterSet(BLK) {};

	// Input			BLK  : block_2D
	// Output			     : returns a new block_2D equivalent
	//				       to input block_2D BLK


MSVCDLC block_2D (const matrix& mx):matrix(mx), ParameterSet() { set_defaults(1); };

	// Input			mx   : matrix
	// Output			     : returns a new block_2D equivalent
	//				       to input matrix mx


MSVCDLC ~block_2D ( ) {};

	// Input			none : 
	// Output			none : destroys block_2D



// ______________________________________________________________________
//                CLASS BLOCK_2D FUNCTIONS, BLOCK WITH BLOCK
// ______________________________________________________________________

MSVCDLL friend block_2D operator + (block_2D& BLK1 , block_2D& BLK2);

	// Input		BLK1, BLK2 : two block_2Ds
	// Return		BLK  : block_2D which is the addition
        //	                       of the two input block_2Ds,
	//			       this = BLK =  BLK1 + BLK2



MSVCDLL void operator += (block_2D& BLK1);

	// Input		BLK1 : a block_2D
	// Return		BLK  : block_2D which is the input
        //	                       input block_2D BLK1 added to itself
	//			       BLK = BLK + BLK1


MSVCDLL friend block_2D operator - (block_2D& BLK1, block_2D& BLK2);

	// Input		BLK1 : a data block_2D
	//			BLK2 : a data block_2D
	// Return		BLK  : block_2D which is the subtraction
        //	                       of the two input block_2Ds,
	//			       BLK = BLK1 - BLK2


MSVCDLL friend  block_2D operator - (block_2D& BLK1);

	// Input		BLK1  : a block_2D
	// Return		BLK: block_2D which is the negative
        //	                     of the input block_2D, BLK = - BLK1


MSVCDLL void operator -= (block_2D& BLK1);

	// Input		BLK1    : a block_2D
	// Return		BLK  : block_2D which is the input
        //	                       input block_2D BLK1 added to itself
	//			       BLK = BLK - BLK1


MSVCDLL void operator = (const block_2D& BLK);

	// Input		this,BLK  : two block_2D
	// Return		copies BLK to this, and destroys the old
        //                      contest of this block_2D


MSVCDLL void operator = (const matrix& mx);

	// Input		this,BLK  : two block_2D
	// Return		copies BLK to this, and destroys the old
        //                      contest of this block_2D



// ********** Basic Function Definitions: D_Block with Scalar **********


MSVCDLL friend block_2D operator * (block_2D& BLK1, complex& z);

	// Input		BLK1  : a block_2D
        //                      z  : a complex number
	// Return		BLK: block_2D which is the multlipication
        //	                     of the block_2D by the complex number
        //                           BLK = z * BLK1


MSVCDLL friend block_2D operator * (complex& z, block_2D& BLK1);

        // Input                z  : a complex number
	//      		BLK1  : a block_2D
	// Return		BLK: block_2D which is the multlipication
        //	                     of the block_2D by the complex number
        //                           BLK = z * BLK1


MSVCDLL friend block_2D operator / (block_2D& BLK1, complex& z);

	// Input		BLK1  : a block_2D
        //                      z  : a complex number
	// Return		BLK: block_2D which is the multlipication
        //	                     of the block_2D by the inverted complex
	//			     number. BLK = (1/z) * BLK1


// ******* Basic Function Definitions: D_Block with Parameter Set ******


MSVCDLL void operator += (const ParameterSet& pset);

	// Input		BLK  : a block_2D
        //                      pset : a parameter set
	// Return		BLK  : block_2D has the paramter set pset added


// ______________________________________________________________________
//                           DATA SET FUNCTIONS
// ______________________________________________________________________


// ______________________________________________________________________
//                         PARAMETER SET FUNCTIONS
// ______________________________________________________________________


MSVCDLL void set_defaults(int level);

	// Input		BLK  : data block (this)
	//			level: level to which defaults are set
	// Return		none : BLK parameter set defaults set


// ______________________________________________________________________
//                 CLASS BLOCK_2D CHECKS AND ERROR MESSAGES
// ______________________________________________________________________


MSVCDLL friend int block_2D_check (block_2D& BLK, int pt);

	// Input		BLK  : block_2D
	// 			pt   : point (integer)
	// Output		int  : insures that point pt is a valid
        //			       point for a block_2D. 0 if bad


MSVCDLL friend int block_2D_check (block_2D &BLK1, block_2D &BLK2);

	// Input		BLK1    : block_2D
	// 	 		BLK2    : block_2D
	// Output		int     : insures that the two block_2Ds
	//			          can computationally intermix
	//			          0 if incompatible pair


// ______________________________________________________________________
//                   CLASS BLOCK_2D I/O FUNCTIONS
// ______________________________________________________________________


MSVCDLL friend void write_pset ( block_2D& BLK, char* filename );

	// Input		BLK      : a block_2D
	//			filename :
	// Return		none     : BLK parameter set written to
	//				   disk file "filename"

 
MSVCDLL friend std::ostream& operator<< (std::ostream& ostr, block_2D& BLK);

	// Input		ostr : string, error message
	// 			BLK  : block_2D
	// Return		     : stream, prints block_2D


MSVCDLL friend void print_pset (block_2D& BLK);

	// Input		BLK  : a block_2D
	// Return		none : BLK parameters printed to std I/O


MSVCDLL friend void print_dset ( block_2D& BLK );

	// Input		BLK  : a block_2D
	// Return		none : BLK data vector printed to std I/O

};

#endif /* __CLASS_BLOCK_2D_H__ */

