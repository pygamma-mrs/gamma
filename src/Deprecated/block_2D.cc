/* block_2D.cc **********************************-*-c++-*-
**							**
** 	            G A M M A				**
**							**
**	Block_2D                   Implementation 	**
**						 	**
**	Copyright (c) 1990, 1991, 1992		 	**
**	Scott Smith				 	**
**	Eidgenoessische Technische Hochschule	 	**
**	Labor fuer physikalische Chemie		 	**
**	8092 Zurich / Switzerland		 	**
**						 	**
**      $Header: $
**						 	**
** Modifications:					**
**						 	**
*********************************************************/

/*********************************************************
**						 	**
** 	Description				 	**
**						 	**
**	The Class Block_2D specifies the structure 	**
**	and properties of a 2-dimensional block of	**
**	simulated (or any other) data.  The data	**
**	block is a matrix in RAM memory and has		**
**	an associated parameter set.  Functions are	**
**	provided for basic algebraic manipulations, 	**
**	parameter specifications, and I/O. 	 	**
**						 	**
*********************************************************/

#ifndef   block_2D_cc_			// Is file already included?
#  define block_2D_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// this is the implementation
#  endif

#include <Deprecated/block_2D.h>
#include <assert.h>
#include <iostream>
#include <stdlib.h>

using namespace std;


// ----------------------------------------------------------------------
// ------------------------ PRIVATE FUNCTIONS ---------------------------
// ----------------------------------------------------------------------

// ______________________________________________________________________
//                   CLASS BLOCK_2D ERROR HANDLING
// ______________________________________________________________________


void block_2D_error(int error)
  {
  std::cout << "\nClass Block_2D: ";
  switch (error)
    {
    case 0:							// (0)
      std::cout << "Program Aborting.";
      break;
    case 1:							// (1)
      std::cout << "Mixing of Data Blocks of Differing Dimensions";
      break;
    case 2:							// (2)
      std::cout << "Problems with Block_2D - Block_2D Operation";
      break;
    case 5:							// (5)
      std::cout << "Block_2D Data Point NOT Accessible";
      break;
    case 6:							// (6)
      std::cout << "Accessed Data Point EXCEEDS Size of Block_2D";
      break;
    case 50:							// (50)
      std::cout << "Error in Block_2D - Block_2D Addition";
      break;
    case 51:							// (51)
      std::cout << "Error in Block_2D - Block_2D Subtraction";
      break;
    default:
      std::cout << "Unkown Error, number " << error;
    }
  std::cout << ".\n";
  }


void block_2D_error(int error, int i, int j)
     
	// Input		error : Error Flag
	//			i,j   : Labels
	// Output		none  : Error Message Output

  {
  std::cout << "\nClass Block_2D: ";
  switch (error)
    {
    case 100:							// (100)
      std::cout << "Size = " << i << " and Current Point Accessed = " << j;
      break;
    case 101:							// (101)
      std::cout << "Block Data Size = " << i << " mixed with Block Data Size " << j;
      break;
    default:
      std::cout << "Unkown Error, number " << error;
    }
  std::cout << ".\n";
  }


volatile void block_2D_fatality(int i)

	// Input		none :
	// Output		none : Stops Execution & Error Message Output

  {
  block_2D_error (i);
  if(i)
    std::cout << "\nClass Block_2D: Program Aborting.\n";
  exit(-1);
  }

// ----------------------------------------------------------------------
// ------------------------- PUBLIC FUNCTIONS ---------------------------
// ----------------------------------------------------------------------

// ______________________________________________________________________
//                 CLASS BLOCK_2D CONSTRUCTION, DESTRUCTION
// ______________________________________________________________________

	
block_2D::block_2D(int dim1, int dim2):matrix(dim1,dim2), ParameterSet()

	// Input			dim1 : row dimension
	// 				dim2 : column dimension
	// Output			none : creates an empty block_2D
	//				       of dimension specified

  { set_defaults(1); }			// Set Block Parameter Defaults


// ______________________________________________________________________
//                CLASS BLOCK_2D FUNCTIONS, BLOCK WITH BLOCK
// ______________________________________________________________________

	
block_2D operator + (block_2D &BLK1 , block_2D &BLK2)

	// Input		BLK1, BLK2 : two block_2Ds
	// Return		BLK        : block_2D which is the addition
        //	                             of the two input block_2Ds,
	//			             this = BLK =  BLK1 + BLK2

  {
  block_2D BLK3(BLK1);
  BLK3 += BLK2;
  return BLK3;
  }


void block_2D::operator += (block_2D &BLK2)

	// Input		BLK1 : a block_2D
	// Return		BLK  : block_2D which is the input
        //	                       input block_2D BLK1 added to itself
	//			       BLK = BLK + BLK1

  {
  if(BLK2.rows() && BLK2.cols())		// if BLK1 has data
    {
    if(!rows() || !cols())			// if BLK1 has no data
      matrix::operator=(BLK2);
    else
      {
      if(block_2D_check(*this, BLK2))		// insure dimensions match
        matrix::operator+= (BLK2);		// add data matrices
      else
        {
        block_2D_error(2);			// error in BLK with BLK operation
        block_2D_fatality(50);			// fatality during BLK + BLK
        }
      }
    }
// sosi -kludge
//  ParameterSet::operator+= (BLK2);			// add BLK2 parameters
  return;
  }


block_2D operator - (block_2D &BLK1, block_2D &BLK2)

	// Input		BLK1 : a data block_2D
	//			BLK2 : a data block_2D
	// Return		BLK  : block_2D which is the subtraction
        //	                       of the two input block_2Ds,
	//			       BLK = BLK1 - BLK2

  {
  block_2D BLK3(BLK1);
  BLK3 -= BLK2;
  return BLK3;
  }

/* *** This function is currently endlessly recursive. */
block_2D operator - (block_2D &BLK1)

	// Input		BLK1  : a block_2D
	// Return		BLK   : block_2D which is the negative
        //	                        of the input block_2D, BLK = - BLK1

  {
  block_2D BLK(BLK1);

#ifdef _MSC_VER
	// *DCT* This call leads to infinite recursion and will not compile on windows.
	cout << "\nFATAL: Call to this deprecated funtion leads to infinite recursion.\n";
	assert(0);
#else
  BLK.matrix::operator= (-BLK1);	// matrix negation
  BLK.ParameterSet::operator= (BLK1);		// copy parameter set BLK1 into BLK
#endif

  return BLK;
  }


void block_2D::operator -= (block_2D &BLK2)

	// Input		BLK1 : a block_2D
	// 			BLK2 : a block_2D
	// Return		none : block_2D BLK1 is modified to have
        //	                       input block_2D BLK2 added to itself
	//			       BLK1 = BLK1 - BLK2

  {
  if(BLK2.rows() && BLK2.cols())		// if BLK1 has data
    {
    if(!rows() || !cols())			// if BLK1 has no data
      matrix::operator= (-BLK2);		// matrix negation
    else
      {
      if(block_2D_check((*this), BLK2))		// insure dimensions match
        matrix::operator-= (BLK2);		// subtract matrices
      else
        {
        block_2D_error(2);			// error in BLK with BLK operation
        block_2D_fatality(51);			// fatality during BLK - BLK
        }
      }
    }
  }


void block_2D::operator = (const block_2D &BLK)

	// Input		this,BLK  : two block_2D
	// Return		copies BLK to this, and destroys the old
        //                      contents of this block_2D

  {
  matrix::operator= (BLK);	// Copy BLK matrix
  ParameterSet::operator= (BLK);	// Copy BLK parameter set
  }


void block_2D::operator = (const matrix& mx)

	// Input		BLK	: A 2D data block (this)
	//			mx	: A  matrix
	// Return		copies BLK to this, and destroys the old
        //                      contents of this block_2D

  {
  matrix::operator= (mx);	// Copy BLK matrix
  ParameterSet::operator= (ParameterSet());	// Copy BLK parameter set
  }


// ______________________________________________________________________
//                CLASS BLOCK_2D FUNCTIONS, BLOCK WITH SCALAR
// ______________________________________________________________________

	
block_2D operator * (block_2D &BLK, complex &z)

	// Input		BLK  : a data block
	//			z    : a complex number
	// Return		BLK1 : a new data block which is the input
	//			       block multiplied by scalar z

  {
  block_2D BLK1(BLK);
  BLK1 *= z;			// scalar-matrix multiplication
  return BLK1;
  }


block_2D operator * (complex &z, block_2D &BLK1)

	// Input		BLK  : a data block
	//			z    : a complex number
	// Return		BLK1 : a new data block which is the input
	//			       block multiplied by scalar z

  {
  block_2D BLK(BLK1);
  BLK *= z;			// scalar-matrix multiplication
  return BLK;
  }


block_2D operator / (block_2D &BLK1, complex &z)

	// Input		BLK  : a data block
	//			z    : a complex number
	// Return		BLK1 : a new data block which is the input
	//			       block divided by scalar z

  {
  block_2D BLK(BLK1);
  BLK /= z;			// scalar-matrix divsion
  return BLK;
  }


// ______________________________________________________________________
//             CLASS BLOCK_2D FUNCTIONS, BLOCK WITH PARAMETER SET
// ______________________________________________________________________


void block_2D::operator += (const ParameterSet &pset)

	// Input		BLK  : a data block
	//			pset : a parameter set
	// Return		none : the input data block is returned
	//			       with parameters from pset added
     
  { 
// sosi - kludge
//ParameterSet::operator+= (pset); }		// add parameters in pset
}


// ______________________________________________________________________
//                  CLASS BLOCK_2D DATA SET FUNCTIONS
// ______________________________________________________________________


// ______________________________________________________________________
//                  CLASS BLOCK_2D PARAMETER SET FUNCTIONS
// ______________________________________________________________________

// *****************  Entire Parameter Set Manipulations ****************

void block_2D::set_defaults(int level)

	// Input		BLK  : data block(this)
	//			level: level to which defaults are set
	// Return		none : BLK parameter set defaults set

  {
  if(level == 0) clear(); 		// Zero any current parameters
  std::string pname("Dim1");
  std::string pstate("Data Set Row Dimension of Block");
  SinglePar par(pname, rows(), pstate);
  push_back(par); 
  pname = std::string("Dim2");
  pstate = std::string("Data Set Column Dimension of Block");
  par = SinglePar(pname, rows(), pstate);
  push_back(par); 
  }


// ______________________________________________________________________
//                        CLASS BLOCK_2D CHECKS
// ______________________________________________________________________


int block_2D_check(block_2D &BLK, int pt)

  {
  int i = 1;
  if(pt >= BLK.rows())			// Insure point not beyond matrix
    {
    block_2D_error(6);			// Point exceeds block data size
    block_2D_error(100, BLK.rows(), pt);
    i = 0;
    }
  return (i);
  }


int block_2D_check(block_2D &BLK1, block_2D &BLK2)

  {
  int i = 1;
  if (BLK1.rows() != BLK2.rows())	// Dimension equivalence check
    {
    block_2D_error(1);			// Dimension mismatch
    block_2D_error(101, (BLK1.rows()), (BLK2.rows()));
    i = 0;
    }
  return (i);
  }


// ______________________________________________________________________
//                   CLASS BLOCK_2D I/O FUNCTIONS
// ______________________________________________________________________

 
void write_pset(block_2D &BLK, char* filename)

	// Input		BLK  : data block
	//			char*: name for disk files
	// Return		none : a file is written to disk containing
	//			       the parameter set in BLK

  {
//  write_pset((ParameterSet)BLK, filename);
  ((ParameterSet)BLK).write(filename);
  }


std::ostream& operator << (std::ostream& ostr, block_2D &BLK)

	// Input		BLK  : data block
	// Return		none : parameter set contained in BLK sent
	//			       to standard I/O	

  { return ostr << (matrix)BLK; }


void print_pset(block_2D &BLK)

	// Input		BLK  : data block
	// Return		none : parameter set contained in BLK sent
	//			       to standard I/O	

  {
  std::cout << "\n\n\t\t\t\tCurrent Block_2D Parameter Set\n\n";
  std::cout << ParameterSet(BLK);
  }


void print_dset(block_2D &BLK)

	// Input		BLK  : data block
	// Return		none : BLK data block is set to standard output

  {
  std::cout << "\n\n\t\t\t\tCurrent Block_2D Data Vector\n\n";
  int nrows = BLK.rows();
  int ncols = BLK.cols();
  int flag = 0;
  for(int i=0; i<nrows; i++)
    for(int j=0; j<ncols; j++)
    {
      std::cout << i+1 << ", " << j+1 << ".\t" << BLK(i,j);
      flag = !flag;
      if (flag)
	std::cout << "\t";
      else
	std::cout << "\n";
    }
  }


#endif /* __CLASS_BLOCK_2D_CC__ */
