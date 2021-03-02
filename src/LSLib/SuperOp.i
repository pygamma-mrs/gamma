/* SuperOp.i */
// Swig interface file

%{
#include "LSLib/SuperOp.h"
%}

%include "std_vector.i"
%include "std_string.i"

%feature("autodoc", "1" );

%rename(__add__)  super_op::operator+ const;
%rename(__iadd__) super_op::operator+=;

%rename(__sub__)  super_op::operator- const;
%rename(__isub__) super_op::operator-=;
%rename(__neg__)  super_op::operator- ();

%rename(__mul__)  super_op::operator* const;
%rename(__imul__) super_op::operator*=;
%rename(__iand__) super_op::operator &=;

%rename(__idiv__) super_op::operator/=;

%rename(__eq__)   super_op::operator== const;


class super_op;  // Forward declaration so following functions will compile

// These declarations were added here,
// as the "friend" versions were not considered declarations
// but rather as defining the relationship to the super_op class.

super_op left(const gen_op& Op); 	// LOp*Op1 = Op*Op1

//super_op left(const matrix& mx, const basis& bs);
//super_op left(const matrix& mx); 

super_op right(const gen_op& Op); 	// LOp*Op1 = Op1*Op

//super_op right(const matrix& mx); 	// LOp*Op1 = Op1*mx
//super_op right(const matrix& mx, const basis& bs);

super_op Hsuper(const gen_op& Heff);
 
super_op U_transform(const gen_op& Op);

//super_op U_transform(const matrix& mx);

super_op commutator(const gen_op& Op); 	// LOp*Op1 = [Op, Op1]

//super_op commutator(const matrix& mx);	// LOp*mx = [mx, mx1] 

super_op d_commutator(const gen_op& Op, const complex& z = complex1);

///super_op d_commutator(const matrix& mx);

super_op d_commutator(const gen_op& Op1, const gen_op& Op2);
super_op d_commutator(const gen_op& Op1, const gen_op& Op2,
                                                             const complex& z);
//super_op d_commutator (const matrix& mx1, const matrix& mx2);



class super_op
{

public:


super_op();

//super_op(const matrix& mx);
//super_op(const matrix& mx,  const matrix& bs);
//super_op(const matrix& mx,  const basis& bs);
//super_op(const std::vector<matrix>& mxc, const std::vector<matrix>& bsc);
//super_op(matrix* mxc, int nc, matrix* bsc=NULL);

super_op(const super_op& LOp1);
super_op(const gen_op& Op1, const gen_op& Op2);

~super_op( );

super_op  operator +  (const super_op& LOp1) const;
super_op& operator += (const super_op& LOp1);
super_op  operator -  (const super_op& LOp1) const;
super_op  operator -  () const;
super_op& operator -= (const super_op& LOp1);

super_op  operator *  (const super_op& LOp1) const;
super_op& operator *= (const super_op& LOp1);
super_op& operator &= (const super_op& LOp1);


//friend gen_op operator * (const super_op& LOp, const gen_op& Op1);
//friend gen_op operator * (const gen_op& Op1, const super_op& LOp);

//friend super_op  operator *  (const super_op& LOp1, const complex& z);
//friend super_op  operator *  (const complex& z,     const super_op& LOp1);
//friend super_op  operator *  (const super_op& LOp1, double d);
//friend super_op  operator *  (double d,             const super_op& LOp1);

super_op& operator *= (const complex& z);
super_op& operator *= (double d);

//friend super_op  operator /  (const super_op& LOp1, const complex& z);
//friend super_op  operator /  (const super_op& LOp1, double d);

super_op& operator /= (const complex& z);
super_op& operator /= (double d);
/*
friend super_op left(const gen_op& Op); 	// LOp*Op1 = Op*Op1
friend super_op left(const matrix& mx, const basis& bs);
friend super_op left(const matrix& mx); 
friend super_op right(const gen_op& Op); 	// LOp*Op1 = Op1*Op
friend super_op right(const matrix& mx); 	// LOp*Op1 = Op1*mx
friend super_op right(const matrix& mx, const basis& bs);

friend super_op commutator(const gen_op& Op); 	// LOp*Op1 = [Op, Op1]
friend super_op commutator(const matrix& mx);	// LOp*mx = [mx, mx1] 

friend super_op d_commutator(const gen_op& Op, const complex& z);
friend super_op d_commutator(const matrix& mx);
friend super_op d_commutator(const gen_op& Op1, const gen_op& Op2);
friend super_op d_commutator(const gen_op& Op1, const gen_op& Op2,
                                                             const complex& z);
friend super_op d_commutator (const matrix& mx1, const matrix& mx2);
 
friend super_op U_transform(const gen_op& Op);
friend super_op U_transform(const matrix& mx);

friend super_op project(const gen_op& Op);

friend super_op project(const matrix& mx);
*/

super_op exp() const;

super_op exp(const complex& t, double cutoff=1.e-12) const;


/*
friend super_op exp(const super_op& LOp1);

friend super_op exp(const super_op& LOp1, const complex& t);

friend super_op pow(const super_op& LOp, int power);
*/

void set_EBR() const;

void set_HBR() const;

void set_DBR() const;

void LOp_base(const super_op& LOp1) const;

void LOp_Hbase(const super_op& LOp1, int warn=0) const;

void LOp_base(const    gen_op& Op) const;

void SetHSBaseOf(const gen_op& Op) const;


int HS()   const;			// Operator Hilbert space
int size() const;			// Operator Liouville space
int dim()  const;			// Operator Liouville space
int LS()   const;			// Operator Liouville space
  
void eigenvalues(int nc = 4, int ri=0) const;

//matrix Mx()     const;
//matrix get_mx() const;

//void put_mx(const matrix& mx1);

basis Bs()       const;
basis get_basis() const;

void put_basis(const basis& Hbs);

basis LBs()        const;
basis get_Lbasis() const;

void put_Lbasis(const basis& Lbs);

complex operator() (int row, int col) const;

void put(int row, int col, const complex& z);

complex get(int row, int col) const;

bool checkLOp(const super_op& LOp1,                 int warn=2) const;
bool checkLOp(const gen_op& Op,                     int warn=2) const;

//bool checkLOp(const matrix& mx,                     int warn=2) const;
//bool checkLOp(const matrix& mx1, const matrix& mx2, int warn=2) const;
//bool checkLOp(const matrix& mx,  const basis& bs,   int warn=2) const;

bool checkLOp(int row,           int col,           int warn=2) const;

void status() const;

int operator == (const super_op& LOp1);

int below(double d) const;

 
//std::ostream& print(std::ostream& ostr, int flag=0) const;
//friend std::ostream& operator << (std::ostream& ostr, const super_op& LOp);

void           write(const std::string& fn)  const;

//std::ofstream& write(std::ofstream& fp) const;


void           read(const std::string& fn);
void           read(const std::string& fn, const gen_op&   Op);
void           read(const std::string& fn, const super_op& LOp1);

//std::ifstream& read(std::ifstream& fp);

//friend super_op Hsuper(const gen_op& Heff);

//friend super_op HsuperX(const gen_op& Heff);

}; 

