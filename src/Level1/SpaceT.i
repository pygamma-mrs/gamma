/* SpaceT.i */

%{
#include "Level1/SpaceT.h"
%}

%feature("autodoc", "1" );

%rename(__assign__) space_T::operator=;

class space_T;

space_T A1(double x, double y, double z);

space_T A1(coord &pt);

space_T A1(row_vector &vx);

complex A1(double x, double y, double z, int m, int l=1);

complex A1(coord &pt, int m, int l=1);

complex A11(double x, double y, double z, int m);

space_T SphA1(complex plus, complex zero, complex minus);

space_T SphA1(coord &pt);

space_T A2(double Aiso, double delzz, double eta,
			double alpha, double beta, double gamma);
space_T A2(coord &Tcomps);

space_T A2(coord &Tcomps, coord &Tangles);

space_T A2(const matrix &mx, double prec);

complex A2(int l, int m, double Aiso, double delzz, double eta);

complex A2(int l, int m, const matrix& mx);

complex A20(int m, double Aiso, double delzz, double eta);

complex A20(int m, const matrix &mx);

complex A21(int m, double Aiso, double delzz, double eta);

complex A21(int m, const matrix &mx);

complex A22(int m, double Aiso, double delzz, double eta);

complex A22(int m, const matrix &mx);

complex T_comp(const space_T &SphT, int L, int M);

space_T T_mult(const space_T &SphT1, const space_T &SphT2);

space_T T_mult(const space_T &SphT1, int l1, const space_T &SphT2, int l2);

space_T T_mult(const space_T &SphT1, int l1,
					 const space_T &SphT2, int l2, int L);
complex T_mult(const space_T &SphT1, int l1,
				const space_T &SphT2, int l2, int L, int M);

space_T T_rot(space_T &SphT1, double alpha, double beta, double gamma);

void T_rot(int num, space_T* SphT, space_T* SphTrot, double alpha, double beta, double gamma);

 void T_rot(int num, space_T* SphT, space_T* SphTrot, matrix* D);

complex T_rot(space_T &SphT1, int l, int m,
				 double alpha, double beta, double gamma);

class space_T
  {

  int rank;		
  row_vector **vx;	
  coord EA;		
  coord PAS_EA;		
  coord PAS;		

void SphTerror(int eidx, int noret=0) const;

void SphTerror(int eidx, const std::string& pname, int noret=0) const;

volatile void SphTfatality(int eidx) const;

volatile void SphTfatality(int eidx, const std::string& pname) const;
         
void updatePAS();

public:

space_T();

space_T(const space_T& SphT);

space_T(const space_T& SphT, int l);

space_T(const SinglePar& par);

virtual ~space_T();

virtual space_T& operator = (const space_T &SphT);

coord PASys() const;

coord PASys_EA() const;

double iso() const;

void iso(double Aiso);

double delz() const;

void delz(double delzz);

double eta() const;

void eta(double ETA);

double alpha() const;

double beta() const;

double gamma() const;

int exists() const;

int exists(int l) const;

int Rank() const;

complex component(int L, int M) const;

double Ccomponent(int r, int c=0) const;

space_T rotate(double alpha, double beta, double gamma) const;

space_T rotate(const coord &EA) const;

complex rotate(int l, int m, double alpha, double beta, double gamma) const;

complex rotate(int l, int m, const coord &EA) const;

SinglePar param(const std::string& pname) const;

SinglePar param(const std::string& pname, const std::string& pstate) const;

//operator ParameterSet();

virtual void write(const std::string &filename);

virtual void read(const std::string &filename);

};
