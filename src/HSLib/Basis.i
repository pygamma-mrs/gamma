/* basis.i */
// Swig interface file.

%{
#include "HSLib/Basis.h"
%}

%include "std_vector.i"
%include "std_string.i"

%feature("autodoc", "1" );

%rename(__eq__)  basis::operator== const;
%rename(__ne__)  basis::operator!= const;

%rename(__assign__) basis::operator=;


class basis : private matrix
{ 

public:
  
basis();
basis(int dim);
basis(const std::vector<int> dims);

//basis(const matrix& mx, int nc=1, int* ncd=NULL);
//basis(const basis&  bs, const matrix& mx);
basis(const basis&  bs);

       ~basis();

basis& operator=(const basis& bs);
//basis& operator=(const matrix& mx);

int         size() const;
int         dim()  const;

std::string name() const;
void        name(const std::string& nm);

int dim_LS() const;

int sub_N() const;

int sub_dim(int ic) const;

int sub_anchor(int ic) const;

int sub_anchor_LS (int ic) const;

int which_sub_LS(int i) const;

//friend basis defLSbasis(const basis& bs);
//friend basis defbasis(const basis& bs);
 
//matrix U()                            const;
//matrix get_matrix()                   const;
//matrix get_mx()                       const;
//matrix convert(const      matrix& mx) const;
//matrix convert_back(const matrix& mx) const; 


//friend basis tensor_product(const basis& bs1, const basis& bs2);

bool operator==(const basis& bs2) const;
bool operator!=(const basis& bs2) const;

bool isDefaultBasis() const;


int refs() const;

bool check(const basis& bs1) const;

//       std::ostream& print(std::ostream& ostr, int full=1) const;
//friend std::ostream& operator<<  (std::ostream& ostr, const basis& bs);

void           write(const std::string& fn) const;

//std::ofstream& write(std::ofstream& fp)     const;

void           read(const std::string& fn);

//std::ifstream& read(std::ifstream&     fp);

double TestBasis(int pf=0) const;

};


