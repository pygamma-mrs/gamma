/* coord.h */
// Swig interface file.

%{
#include "Level1/coord.h"
%}

%include "std_string.i"

%rename(__add__)  coord::operator+ const;
%rename(__iadd__) coord::operator+=;
%rename(__sub__)  coord::operator- const;
%rename(__isub__) coord::operator-=;

%rename(__mul__)  coord::operator* const;
%rename(__imul__) coord::operator*=;
%rename(__div__)  coord::operator/ const;
%rename(__idiv__) coord::operator/=;

%rename(__eq__)  coord::operator== const;
%rename(__ne__)  coord::operator!= const;
%rename(__lt__)  coord::operator<  const;
%rename(__gt__)  coord::operator>  const;

%rename(__assign__) coord::operator=;



matrix Rmx(double alpha, double beta, double gamma);

class coord
{
    

double cx;				// First ordinate
double cy;				// Second ordinate
double cz;				// Third ordinate
				// Output format flags
static int         PTotype;  		//   Output type
static int         PTscience; 	//   Scientific notation
static int         PTdigits;		//   Digits per point
static int         PTprecise;		//   Digits after decimal
static std::string PTform;		//   Output format
static coord       DefCoord;		// Default coordinate
static double      OrdCutoff;		// Zero cutoff for ordinate

  
private:

void PTerror(int eidx,                        int nr=0) const;
void PTerror(int eidx, const std::string& pn, int nr=0) const;

volatile void PTfatal(int eidx)                                  const;
volatile void PTfatal(int eidx, const std::string& pn)           const;

bool SetPtCartesian(const   ParameterSet& pset, int indx, int warn=1);
bool SetPtSpherical(const   ParameterSet& pset, int indx, int warn=1);
bool SetPtCylindrical(const ParameterSet& pset, int indx, int warn=1);


public:

coord( );
coord(double xx, double yy=0, double zz=0);
coord(const coord& pt1);
coord(const ParameterSet& pset, int idx=-1, int warn=2);
coord(const SinglePar& par);
~coord( );

coord& operator= (const coord& pt);

double get(int i)    const;
double x()           const;
void   x(double xx);
double y()           const;
void   y(double yy);
double z()           const;
void   z(double zz);
void   xyz(double xx, double yy, double zz);
void   xyz(const coord& pt1);
double norm() const;
 
double Rad() const;
double Rad(const coord& pt2) const;
//friend double Rad(double x, double y, double z);

double theta() const;
double theta(const coord& pt2) const;

//friend double theta(double x, double y, double z);

double phi()                 const;
double phi(const coord& pt2) const;

//friend double phi(double x, double y, double z);

void invert();

static matrix Rz(double phi, int rad=0);
static matrix Rx(double theta, int rad=0);
static matrix Ry(double theta, int rad=0);
 
coord xrotate(double theta, int rad=0) const;
coord yrotate(double theta, int rad=0) const;
coord zrotate(double phi, int rad=0) const;

matrix Ralpha(double alpha, int rad) const; 
matrix Rbeta(double  beta,  int rad) const;
matrix Rgamma(double gamma, int rad) const;

matrix REuler(double alpha, double beta, double gamma, int rad) const;

/*
friend matrix Rmx1(double alpha);
friend matrix Rmx2(double beta);
friend matrix Rmx3(double gamma);
friend matrix Rmx(double alpha, double beta, double gamma);
friend matrix Rmx(coord EA);
*/

coord rotate(double alpha, double beta, double gamma);

coord rotate(coord& EA);

coord trans_x(double    delx);
void  trans_x_ip(double delx);
coord trans_y(double    dely);
void  trans_y_ip(double dely);
coord trans_z(double    delz);
void  trans_z_ip(double delz);

coord translate(double    delx, double dely=0, double delz=0);
void  translate_ip(double delx, double dely=0, double delz=0);

coord translate(const    coord& del) const;
void  translate_ip(const coord& del);


coord operator +  (const coord& del) const;
coord operator -  (const coord& del) const;
coord&  operator += (const coord& del);
coord&  operator -= (const coord& del);

coord  operator *  (double r) const;
coord&   operator *= (double r);
coord  operator /  (double r) const;
coord&   operator /= (double r);

bool operator==(const coord& pt) const;
bool operator!=(const coord& pt) const;
bool operator>(const  coord& pt) const;
bool operator<(const  coord& pt) const;


/*
friend double Rad(const coord& pt1, const coord& pt2);
friend double theta(const coord& pt1, const coord& pt2);
friend double phi(const coord& pt1, const coord& pt2);
friend coord cdvect(const coord& pt1, const coord& pt2);
friend coord  operator *  (double r, const coord& pt1);
*/
 
//friend coord operator * (const matrix& mx, const coord& pt);
 
SinglePar param(const std::string& pname) const;
SinglePar param(const std::string& pname, const std::string& ps) const;

int read(const std::string& filein, int indx, int warn=1);
int read(const ParameterSet&  pset, int indx, int warn=1);

static int length();

//friend void coord_setf (int otype, int science, int digits, int precise);
//friend void coord_getf (int& otype, int& science, int& digits, int& precise);
//std::ostream& print(std::ostream& ostr) const;
//friend std::ostream& operator << (std::ostream& ostr, const coord& pt);
//friend std::istream& operator >> (std::istream& istr, coord& pt);

coord Cart2Sph(int rad=1) const;
coord Sph2Cart(int rad=1) const;
coord Cart2Cyl(int rad=1) const;
coord Cyl2Cart(int rad=1) const;
coord Sph2Cyl(int  rad=1) const;
coord Cyl2Sph(int  rad=1) const;
 
static coord getDefCoord();
static void  setDefCoord(const coord& dpt);
 
static void SetCutoff(double co=-1);

};

extern const coord UnitX;		// Unitx  = 1i + 0j + 0k
extern const coord UnitY;		// Unity  = 0i + 1j + 0k
extern const coord UnitZ;		// Unitz  = 0i + 0j + 1k
extern const coord coord0;		// coord0 = 0i + 0j + 0k

