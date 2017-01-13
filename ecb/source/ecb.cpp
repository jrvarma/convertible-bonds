/*
    ecb.cpp for Valuation of Convertible Bonds
    Copyright (C) 2001, 2004  Prof. Jayanth R. Varma, jrvarma@iimahd.ernet.in,
    Indian Institute of Management, Ahmedabad 380 015, INDIA

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program (see file COPYING); if not, write to the 
    Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
    Boston, MA  02111-1307  USA
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "stddef.h"
using namespace std;

#include "newton.h"
#include "lattice.h"
#include "OpenXMLParser.h"

//MyLattice 

class MyNode
{
public:
  double s;  // stock price
  double vB; // coupon and redemption part of bond value 
  double vO; // option value or option part of convertible
  double v;  // = vB + vO = total value of convertible
  char  ch;  // flag to denote action (exercise, call or put)
};

typedef lattice<MyNode, timemap> MyLattice;
// the nodes of the lattice are of the class MyNode defined above
// the conversion from lattice stages to real time is the standard timemap defined in
// lattice.h where the nodes are equally spaced at distances of dt

//forward declarations
class lattice_params;
class bond; 
class Out;
class Input;
class Results;
class Params;
class option;
class stock;
class MyError{public: int n;};

typedef enum {continuous, annual, semi_annual, price} Compounding;

/* function declarations */

double f_of_kd(const double x, MyLattice &L);
// compute bond value as function of kd (i.e. set kd=x, compute the bond value using 
// the lattice and return this bond value)
double f_of_sigma(const double x, MyLattice &L);
// compute option/bond value as function of sigma (i.e. set sigma=x, compute the value using 
// the lattice and return this value)
FILE *getfile(int argc, char **argv);
// open and return the file given as command line argument. 
// If there is no command line argument returns stdin (standard input)
// This means that the input can be piped into the program if necessary
void Vega(MyLattice &L, Results& res);
// compute vega
void binomial(MyLattice &L, Results& res);
// main computation routine which initialises the binomial lattice and rolls it back
// to compute the option values
void init_lattice(MyLattice& L);
// initialise binomial lattice
void ytm(MyLattice &L, Results& res);
// compute ytm
void print_params(void);
// print out lattice parameters (u, d, p etc.)
void node_value(struct nodevalue& v, const int n);
// compute the option/bond value at a node
void rollback(MyLattice& L, Results& res);
// roll back lattice computing values at time t-1 from values at time t
void compute(MyLattice& L, Results& res);
// compute option/bond values and greeks from the fully folded lattice
void print_results(const Results& res);
// print option/bond values and other results
MyLattice& data(FILE *infile);
// read input file, set up lattice, bond, stock and other classes and return the lattice
void print_data(const Input& input);
// print input data
void print_lattice(const MyLattice& L);
// print the full lattice
char *rate_type(const Compounding type);
// return a text representation of the given Compounding type
double effv_rate(const double rate, const Compounding type);
// return equivalent continuously compounded rate
double convert_rate(const double rate, const Compounding type);
// convert continuously compounded rate to equivalent rate at required Compounding type
int enum_value(const AdvXMLParser::Attribute &elem, const char* str_list[]);
  //convert a text string into an enumerated value using the given list of strings
char *dbl2str(const double x);
// convert double to text representation
char *rate2str(const double x);
// convert continuously compounded rate to text representation of equivalent rate
// at required Compounding type using global variable pOut->ratetype


/*global variables */
// these are read only variables
const timemap *TMptr;              // time map associated with the lattice L
const Out *pOut;                   // desired outputs
const Params *pParams;             // computational parameters
const stock *pS;                   // stock details
const bond *pB;                    // bond details (f_of_kd uses const cast to modify this)
const lattice_params *pLP;         // lattice parameters (f_of_sigma uses const cast to modify)
AdvXMLParser::Document xmldoc("Output"); // output xml document
char dblformatstd[] = "%.6g";
char dblformatcompact[] = "%.3g";
char dblformatnarrow[] = "%.2g";
char *dblformat = dblformatstd;
ostringstream MyErrMsg;
MyError MyErr;


/* classes*/

class lattice_params  // defines the lattice parameters
{
public:
  double u, d, p, sigma, RminusQ;
  void set_params(const double sigma0, const double RminusQ0);
};

void lattice_params::set_params(const double sigma0, const double RminusQ0)
  // we use a default constructor and then use set_params for actual construction
  // set_params allows us to reset the object with destroying and recreating it
{
  sigma = sigma0;                     // volatility
  RminusQ = RminusQ0;                 // risk free rate minus continuous dividend yield
  double dt = TMptr->dt();            // time interval between successive stages
  double a = exp(RminusQ*dt);         // risk neutral growth in stock price during interval dt
  u = exp(sigma*sqrt(dt));            // upmove 
  d = 1.0/u;                          // downmove d=1/u defines recombining lattice
  p = (a - d) / (u - d);              // probability of upmove
}


class stock  // defines the stock details
{
public:
  double su0, // unadjusted stock price
    su,       // adjusted for country discount
    s,        // further adjusted for PV of known dividends
    q,        // continuous dividend yield
    *Div,     // known dividend amounts
    *xD,      // known dividend ex-div dates
    discount; // country discount on stock price
  int ndiv;   // number of known dividends
  void set_params(const double sdot_div, const double r);
  double pvdiv(const double t, const double r) const;
    
private:
  void pv0(const double sdot_div, const double r);
};

void stock::set_params(const double sdot_div, const double r)
{
  su = su0*(1-discount);        // adjust country discount
  pv0(sdot_div, r);             // compute PV at time zero of each known dividends
  s = su - pvdiv(0, r);         // adjust for PV of all known dividends
}
void stock::pv0(const double sdot_div, const double r)
{
  for (int i = 0; i < ndiv; i++)
    Div[i] *= exp((sdot_div - r)*xD[i]);   // convert to base currency and discount
  // Div is now pv at time 0
}

double stock::pvdiv(const double t, const double r) const
  // apart from being called with t=0 during set_params, this method is called 
  // repeatedly during lattice computations for different values of t
{
  double sum = 0;
  for (int i = 0; i < ndiv; i++)
    if (xD[i] >= t)
      sum += Div[i];      // sum is PV at 0 of all dividends receivable at or after time t
  return (sum*exp(r*t));  // convert to PV at time t
}

struct pv{
  double price;
  double duration;
}; // used to return PV and duration of the straight bond

struct nodevalue
{
  double st, // stock price at this node
    liveval, // value if option is kept alive
    deadval, // value if option is exercised now
    callval, // intermediate variable representing value if bond is called
    refval,  // value considering possible exercise and possible calls
    putval,  // intermediate variable representing value if bond is put
    bestval, // value considering possible exercise, possible calls & possible puts
    bondval, // coupon and redemption part of value at this node
    optval, // option part of value at this node
    totval,   // total value at this node including coupon if any
    coup;    // coupon payment at this node
  bool exercise, // is the option exercised?
    called,      // is the bond called?
    converted,   // is it converted (exercised or converted on call)
    put_it;      // is the bond put?
  char ch;       // characted representing action taken at this node (exercise, call, put)
};  // used in computing value of option/convertible at a node

class soft_call
{
  // details of soft call option to issuer
  // one instance for each soft call option
public:
  bool callable(const double st, const int n) const
  {
    bool time_ok =  (n  <= TMptr->n(time_max) && n  >= TMptr->n(time_min));
    // are we within the time window of the soft call (time_min to time_max)?
    bool price_ok = (st >= price);
    // is the stock price above the floor price for the soft call ?
    return(time_ok && price_ok); // bond is called if both conditions are met
  }
  void set (const double tmin, const double tmax, const double st_price, 
	    const double redeem_value)
    // we use the default constructor
    // and then use this method to set private variables 
  {
    time_min = tmin; time_max = tmax; price = st_price; redeem =redeem_value;
  }
  // methods to read private variables
  double Time_min() const {return time_min;} // earliest time of call
  double Time_max() const {return time_max;} // latest time of call
  double Price() const {return price;}       // floor stock price for call
  double Redeem() const {return redeem;}     // redemption value if called
private:    
  double time_min, time_max, price, redeem; 
};


class put_option
{
  // details of put option 
  // one instance for each put option
public:
  void set(const double tmin, const double tmax, const double yld_or_price, 
	   const Compounding yld_type)  
    // we use the default constructor
    // and then use this method to set private variables 
  {
    time_min = tmin; time_max = tmax; yield = yld_or_price;
    type = yld_type;
  }
  bool puttable(const int n) const
  {
    return n  <= TMptr->n(time_max) && n  >= TMptr->n(time_min);
    // bond is puttable if we are within the time window (time_min to time_max)
  }
  double value(const bond &b, const int n) const; // return put value
  // return private variables
  double Time_min() const {return time_min;} // earliest time of put
  double Time_max() const {return time_max;} // latest time of put
  double Yield() const {return yield;}       // yield to put or put value
  Compounding Type() const {return type;}    // compounding type of the put yield
private:
  double time_min, time_max, yield; 
  Compounding type;
};

class option
{
  // this class is used through the derived class bond
public:
  double R,         // risk free rate
    x,              // exercise price
    te,             // earliest exercise (0 for American, tl for European)
    tl;             // option maturity
  bool convertible, // is this a convertible or a stock option
    call,           // is it a call option?
    put;            // is it a put option?
  int ne,           // lattice time corresponding to te
    nl;             // lattice time corresponding to tl
};

class bond : public option
{
public:
  //static constants
  static const int disc_at_rf = 0;
  static const int disc_at_kd = 1;
  static const int disc_prob_wt = 2;
  static const int disc_val_wt = 3;
  //public variables
  //this class allows user to directly modify these variables
  //but user must call set_params after such modifications
  double coupon, redeem, bond_price, maturity, cpn_prd;
  int cpn_freq;
  double sdot_coupon, sdot_redeem, sdot_o_price, sdot_p_price, sdot_o_redeem, s_redeem;
  double kd, wrf, wkd;
  int disc_rate;
  //constructors and set up 
  bond(void) : _no(0), _np(0) {}
  bond(int no, int np);
  void set_no(int no);
  void set_np(int np);
  void set_params();
  //return private variables
  int no(void) const {return _no;}
  int np(void) const {return _np;}
  class soft_call& co(const int i) const {return _co[i];} 
  // this is const but user may modify _co[i]
  class put_option& po(const int i) const {return _po[i];}
  // this is const but user may modify _po[i]
  double straight() const {return _straight;};
  //other methods
  double cpn(const int i) const;
  struct pv fstraight(const double kd) const;
  void put_call_value(struct nodevalue &v, const int n) const;
private:
  class soft_call *_co;
  class put_option *_po;
  int _no, _np;
  double _straight;
};
    
bond::bond(int no, int np) : _no(no), _np(np)
{
  _co = new soft_call[no];
  _po = new put_option[np];
}    

void bond::set_no(int no)
{
  if(_no > 0)
    delete [] _co;
  _no = no;
  if (_no > 0)
    _co = new soft_call[no];
}    
    
void bond::set_np(int np)
{
  if(_np > 0)
    delete [] _po;
  _np = np;
  if(_np > 0)
    _po = new put_option[np];
}    

struct pv bond::fstraight(const double kd) const
{
  double wkd = exp(kd*cpn_prd);
  struct pv f;
  f.price = s_redeem*wkd;
  f.duration = maturity*s_redeem*wkd;
  for (double t = maturity; t >= 0; t-=cpn_prd){
    double coup = (t == 0) ? 0.0 : coupon*cpn_prd*exp(sdot_coupon*t);
    f.price = coup + f.price/wkd;
    f.duration = t*coup + f.duration/wkd;
  }
  f.duration = f.duration/f.price;
  return f;
}


double bond::cpn(const int i) const
  //calculate coupons at lattice time i
{
  double dt = TMptr->dt();
  double _cpn = 0;
  for (double t = maturity; t > 0; t -= cpn_prd){
    if (i == TMptr->n(t))
      _cpn += coupon*cpn_prd*exp(sdot_coupon*t-kd*(t/dt - i));
  }
  return _cpn;
}

void bond::put_call_value(struct nodevalue &v, const int n) const
  //calculate option value considering calls and puts
  //on entry v.deadval, v.liveval, v.exercise and v.st are set
{
  v.called = false;
  double t = n*TMptr->dt();
  v.refval = (v.exercise) ? v.deadval : v.liveval;
  for (int i = 0; i < _no; i++){
    v. callval = x*co(i).Redeem()*exp(sdot_o_price*t);
    // value to bondholders if they redeem instead of converting
    v.callval = max(v.deadval, v.callval);
    // bondholders choose their best option - redeem or convert
    if (co(i).callable(v.st, n) && v.callval < v.liveval){
      v.called = true; // this call is worth exercising
      v.refval = min(v.refval, v.callval); // is this call better than other calls seen so far?
    }
  }
  v.converted = v.called && (v.refval == v.deadval); 
  // investor converted instead of redeeming
  v.put_it = false;
  v.bestval = v.refval;
  for (int i = 0; i < _np; i++){
    bool time_ok = po(i).puttable(n);
    if (time_ok){
      v.putval = po(i).value(*this, n);
      v.put_it = v.put_it || (time_ok && v.putval > v.refval);
    }
    if(v.put_it) v.bestval = max(v.bestval, v.putval);
  }
  if (v.put_it){
    v.bondval =  v.bestval + v.coup;
    v.optval = 0;
    }else if (v.converted || v.exercise){
      v.optval = v.refval;
      v.bondval = v.coup;
    }else{
      v.bondval += v.coup; // leave v.optval unchanged and add coupon to v.bondval
  }
  v.ch = (v.put_it) ? 'P' : (v.called) ? 'C' : (v.exercise) ? 'E' : ' ';
}

void bond::set_params(void)
{
  const double dt = TMptr->dt();
  //discount conversion flows at risk free rate
  wrf = exp(R*dt);
  //discount bond flows at bond yield
  wkd = exp(kd*dt);
  ne = TMptr->n(te); // convert te to lattice time
  nl = TMptr->n(tl); // convert tl to lattice time
  if (convertible){
    s_redeem = exp(sdot_redeem*maturity)*redeem; // convert to base currency
    cpn_prd = 1.0/cpn_freq;
    _straight = fstraight(kd).price; // compoute straight bond value
  }
}

inline double put_option::value(const bond &b, const int n) const
{
  double t = TMptr->t(n); // convert lattice time to real time
  if(type == price){
    //convert to base currency
    return b.x*(yield*exp(b.sdot_p_price*t)/100);
  }else{
    double tv_c = 0;
    for (double s = b.maturity; s > 0; s -= b.cpn_prd){
      if (TMptr->n(s) <= n)
	tv_c += b.coupon*b.cpn_prd*exp(b.sdot_coupon*s+yield*(t-s));
      //value today of all past coupons
    }
    return b.x*(exp(yield*t) - tv_c); // amount needed now to achieve desired yield to put
  }
}

class Results
{
public:
  double delta, Gamma, theta1, theta2, theta3, vega; // option greeks
  double f, f42, fl1, fl2; // intermediate variables to calculate greeks and option value
  double optval, straight, bondval; // value of options, of staright bond and of convertible
  double ytm; 
  double x, rf, price; // exercise/conversion price, risk free rate and bond price
  // above variables copied here for convenience
  const char *YTMerrmsg; // error if any in computing ytm
  int YTMiter; // number of iterations used in computing ytm
};

class Params
{
  //computational parameters
public:
  bool find_ytm;         // is ytm to be computed?
  bool ignore_put;       // are put options to be ignored
  bool ignore_call;      // are call options to be ignored
  double toler;          // tolerance for computing ytm
  int max_iter;          // max number of iterations for computing ytm
  int find_vega;         // is vega to be computed?
  static const int default_toler_reciprocal = 100; // default value 
  // this is actually default_toler = 0.01 but ISO C++ forbids non int
  // member constants
  static const int default_max_iter = 20;     // default value 
};

class Out
{
public:
  char opt_name[20];
  bool data, // print input data
    params,  // print lattice parameters
    lattice, // print full lattice
    option,  // print option results
    bond;    // print bond valuation results
  int silent;// internal variable to supress printing during ytm/vega iterations
  Compounding ratetype; // compounding type to be used while printing interest rates
};

int main(int argc, char *argv[])
{
  static Results res;
  MyErrMsg << "<?xml version = \"1.0\" encoding = \"UTF-8\"?>" << endl << "<Error>" << endl;
  int hdr_len = MyErrMsg.str().size();
  try{
    MyLattice& L = data(getfile(argc, argv));
    binomial(L, res);
    Out& out = *(const_cast<Out *>(pOut));
    if (pParams->find_vega){
      out.silent++;
      Vega(L, res);
      out.silent--;
    }
    if(pParams->find_ytm){
      out.silent++;
      ytm(L, res);
      out.silent--;
    }
    if (!pOut->silent && (pOut->option || (pB->convertible && pOut->bond)))
      print_results(res);
    string strXML = xmldoc.GenerateXML();
    cout << strXML << endl;
  }
  catch(std::bad_alloc){
    MyErrMsg << "Insufficient memory" << endl;
    throw;
  }
  catch(...){
    //Write XML file for error message and exit
    if(MyErrMsg.str().size() == hdr_len)
      MyErrMsg << "Unknown (Internal) Error" << endl;
    cout << MyErrMsg.str() << endl << "Exiting\n</Error>";
  }
}


void ytm(MyLattice& L, Results& res)
{

  if (pB->disc_rate == bond::disc_at_rf){
    res.YTMerrmsg = "YTM cannot be calculated when discounting at risk free rate";
    res.YTMiter = 0;
    return;
  }
  newton_solver<MyLattice> ns(f_of_kd, &L); // initialise the newton solver
  ns.xmin = -1;                             // interest rate >= -1
  ns.target = pB->bond_price;               // set target for solver
  ns.max_iter = pParams->max_iter;          // set max number of iterations
  ns.tolerance = pParams->toler;            // set tolerance
  res.ytm = ns.solve(pB->kd, res.bondval);  // solve
  res.YTMerrmsg = ns.err_msg();             // store error message if any
  res.YTMiter = ns.iter_no();               // store number of iterations
}

void Vega(MyLattice& L, Results& res)
{
  numdiff<MyLattice> nd(f_of_sigma, &L);
  //initialise the numerical differentiator
  nd.type = (pParams->find_vega == 1) ? 
    numdiff_forward : numdiff_central;
  //set numerical differentiation method
  //use defaults for other parameters of numerical differentiator
  res.vega = nd.deriv(res.f, pLP->sigma);
  // differentiate
}

double f_of_kd(const double x, MyLattice& L)
  //function used by newton solver to compute ytm
{
  static Results res; // use this instead of global res
  bond * _pB = const_cast<bond *>(pB); // allow us to modify it
  double old_kd = _pB->kd;
  _pB->kd = x;                         // change kd
  _pB->set_params();                   // set_params changes other paramaters accordingly
  binomial(L, res);                    // process lattice
  _pB->kd = old_kd;                    // restore  kd
  _pB->set_params();                   // set_params changes other paramaters accordingly
  return res.bondval;                  // return the bond value for kd=x
}

double f_of_sigma(const double x, MyLattice& L)
{
  //function used by numerical differentiator to compute vega
  static Results res; // use this instead of global res
  double old_sigma = pLP->sigma;
  (const_cast<lattice_params *>(pLP))->set_params(x, pLP->RminusQ);
  //modify sigma after casting away const-ness
  binomial(L, res); // process the lattice
  (const_cast<lattice_params *>(pLP))->set_params(old_sigma, pLP->RminusQ);
  //restore old sigma after casting away const-ness
  return res.f;     // return the lattice option value for sigma=x
}

void binomial(MyLattice& L, Results& res)
{
  
  if (!pOut->silent && pOut->params) print_params(); // print lattice parameters if needed
  init_lattice(L); // initialise the lattice
  while(L.n() > 0){
    if (!pOut->silent && pOut->lattice) print_lattice(L);
    rollback(L, res); // roll back the lattice
  }
  if (!pOut->silent && pOut->lattice) print_lattice(L);
  compute(L, res); // compute the results at the end
}

void init_lattice(MyLattice& L)
{
  L.rewind(); // put us in the last stage of the lattice
  int n = L.n(); 
  const double d = pLP->d;
  const double u = pLP->u;
  const double s = pS->s;
  //fill in the stock prices at the last stage of the lattice
  L[n/2].s = (2*(n/2) == n) ? s : s*u;
  for (int i = n/2 + 1; i <= n; i ++)
    L[i].s = L[i-1].s*d*d;
  for (int i = n/2 - 1; i >= 0; i --)
    L[i].s = L[i+1].s*u*u;
  struct nodevalue v;
  v.coup = (!pB->convertible) ? 0: pB->x*pB->cpn(n);// coupon at this node
  for (int i = 0; i <= n; i++){
    v.st = L[i].s+pS->pvdiv(TMptr->t(n), pB->R);
    //stock price in lattice excludes the PV of known future dividends
    //we add this back to get the true market price of the stock
    v.liveval = v.bondval = (!pB->convertible) ? 0 : pB->x*pB->s_redeem;
    v.optval = 0;
    //if option is not exercised, value is 
    //0 for stock option and bond redemption value for convertible
    node_value(v, n);  //compute the value at this node
    L[i].vB = v.bondval;
    L[i].vO = v.optval;
    L[i].v = v.optval + v.bondval;
    L[i].ch = v.ch;             //store the action flag
  }
}

void rollback(MyLattice& L, Results& res)
{
  if (L.n() == 4) res.f42 = L[2].v; /* This value is used for theta estimate */
  L.swap(); //prepare to rollback by marking current nodes as "old" and new nodes as current
  int n = L.n(); // current lattice time
  struct nodevalue v;
  const double p = pLP->p;
  v.coup = (!pB->convertible) ? 0: pB->x*pB->cpn(n); // coupon at this node
  for (int i = 0; i <= n; i++){
    L[i].s = L.upmove(i).s*pLP->d; // compute stock price at this node
    v.bondval = (p*L.upmove(i).vB + (1-p)*L.dnmove(i).vB)/pB->wkd;
    v.optval =  (p*L.upmove(i).vO + (1-p)*L.dnmove(i).vO)/pB->wrf;
    v.liveval = v.bondval + v.optval;
    //live value is discounted expected value of value at lattice time n+1
    v.st = L[i].s+pS->pvdiv(TMptr->t(n), pB->R);
    //stock price in lattice excludes the PV of known future dividends
    //we add this back to get the true market price of the stock
    node_value(v, n);  //compute the value at this node
    L[i].vB = v.bondval;
    L[i].vO = v.optval;
    L[i].v = v.optval + v.bondval;
    L[i].ch = v.ch;            //store the action flag
  }
}

double deadvalue(const double st, const bool last)
{
  double value;

  if (pB->put)
    value = pB->x - st;
  else if (pB->call)
    value = st - pB->x;
  if(pB->convertible)
    value = st;
  if (last) value = max(0.0, value);
  return value;
}

void node_value(struct nodevalue& v, const int n)
{
  v.deadval = deadvalue(v.st, TMptr->N() == n);
  v.exercise = (n >= pB->ne && n <= pB->nl && v.deadval > v.liveval);
  if(v.exercise){
    v.optval = v.deadval;
    v.bondval = 0;
  }// else values set in rollback are left unchanged
  if(!pB->convertible){
    v.ch = (v.exercise) ? 'E' : ' ';
  }else{
    pB->put_call_value(v, n); // compute value taking call and put into account
  }
}

void compute(MyLattice& L, Results& res)
{
  const double dt = TMptr->dt();
  const double s = pS->s;
  const double sigma = pLP->sigma;
  const double RminusQ = pLP->RminusQ;
  res.f = L[0].v;
  res.delta = (L.nPlus2(0).v - L.nPlus2(2).v)/(L.nPlus2(0).s - L.nPlus2(2).s);
  res.Gamma = ( (L.nPlus2(0).v - L.nPlus2(1).v)/(L.nPlus2(0).s - L.nPlus2(1).s) -
		(L.nPlus2(1).v - L.nPlus2(2).v)/(L.nPlus2(1).s - L.nPlus2(2).s) ) /
    (0.5*(L.nPlus2(0).s - L.nPlus2(2).s));
  res.theta1 = (L.nPlus2(1).v - L[0].v)/(2*dt);
  if (TMptr->N() >= 4) res.theta2 = (res.f42 - L[0].v)/(4*dt);
  res.theta3 = RminusQ*res.f - RminusQ*s*res.delta - 0.5*sigma*sigma*s*s*res.Gamma;
  if(pB->convertible){
    res.x = pB->x;
    res.rf = pB->R;
    res.price = pB->bond_price;
    res.bondval = res.f*100/res.x; 
    res.straight= 100*pB->straight();
    res.optval = res.bondval - res.straight;
  }
}

char *rate_type(const Compounding type)
  //return a text representation of the compounding type
{
  static char s[15];
  switch(type){
  case continuous:  strcpy(s, "continuous)"); break;
  case annual:      strcpy(s, "annual"); break;
  case semi_annual: strcpy(s, "semi-annual"); break;
  case price:       strcpy(s, "price"); break;
  }
  return(s);
}

double effv_rate(const double rate, const Compounding type)
  //convert from given compounding type to continuously compounded equivalent rate
{
  switch(type){
  case continuous:	return(rate);
  case annual:		return(log(1.0+rate));
  case semi_annual:	return(2.0*log(1.0+rate/2.0));
  default:       	return(rate);
  }
}

double convert_rate(const double rate, const Compounding type)
  //convert from continuously compounded rate to equivalent rate of given compounding type
{
  switch(type){
  case continuous:	return(rate);
  case annual:		return(exp(rate)-1.0);
  case semi_annual:	return(2.0*exp(rate/2.0)-2.0);
  default:       	return(rate);
  }
}

struct Input
//define a set of variables used while reading input
//they are not directly used in the actual computation
//but are transformed into something else
{
  Compounding R_type, q_type, sdot_type, kd_type;
  double R0, q0, kd_0;
  int N0, BondMatYears, BondMatDays;
  double sdot_div;
};

using namespace AdvXMLParser;

int enum_value(const Attribute &elem, const char* str_list[])
  //convert a text string into an enumerated value using the given list of strings
{
  int j = -1;
  const char *str = elem.GetValue().c_str();
  for(int i = 0; *str_list[i]; i++){ //the empty string is a sentinel value
    if (strcmp(str,str_list[i]) == 0){
      j = i;
      break;
    }
  }
  return j;
}

//these variables are used by the function enum_value above
//each array must have "" as its last element 
//this is a sentinel to prevent an infinite loop in the function
const char* CompoundingList[] = {"continuous", "annual", "semiAnnual", ""};
const char* DiscountTypeList[] = {"riskFree", "bondYield", "probWeighted", 
				  "probAndValueWeighted", ""};
const char* VegaTypeList[] = {"none", "forwardDifferences", "centralDifferences", ""};
const char* PutValueTypeList[] = {"continuousYieldPercent", "annualYieldPercent", 
				  "semiAnnualYieldPercent", "amount", ""};
const char* BoolList[] = {"false", "true", ""};
const char* optionTypeList[] = {"call", "put", ""};

MyLattice& data(FILE *infile)
{
  //these classes are created here and global read only pointers are set up to point to these
  static Input in;
  static Out out;
  static Params params;
  static stock S;
  static bond B;
  static lattice_params LP;
  //this is returned
  static MyLattice L;
  
  //set up global read only pointers 
  pOut = &out;
  pParams = &params;
  pS = &S;
  pB = &B;
  pLP = &LP;
  //read input file into buffer
  long nSize = 20000;
  char* pBuffer = OpenXmlFile(infile, nSize);
  try{
    Parser parser;
    //parse the input XML file
    auto_ptr<Document> pDocument(parser.Parse(pBuffer, nSize));
    {
      //is there a convertible root element in the XML file?
      const Element& test = pDocument->GetRoot()("convertible");
      B.convertible = ! test.IsEmpty();
    }
    //for convertibles the root but one element is convertible, else it is American
    const char* rootstr = (!B.convertible) ? "American" : "convertible";
    const Element& root = pDocument->GetRoot()(rootstr);
    const char *outputs, *stock, *computation; // these are assigned below
    if (!B.convertible){
      B.call = (enum_value(root("Option")["type"], optionTypeList) == 0);
      B.put = ! B.call;
      stock = "Stock"; outputs = "Outputs"; computation = "Computation";
    }else{
      B.put = false; B.call = true; 
      stock = "stock"; outputs = "outputs"; computation = "computation";
    }
    in.N0 << root(computation)["numberOfLatticePoints"];
    int N = (in.N0 < 2) ? 2 : in.N0;
    //enum_value converts false/true to 0/1 and the bool cast does the rest of the work
    out.data = (bool) enum_value(root(outputs)["inputData"], BoolList);
    out.params = (bool) enum_value(root(outputs)["latticeParameters"], BoolList);
    out.lattice = (bool) enum_value(root(outputs)["fullLattice"], BoolList);
    out.option = (bool) enum_value(root(outputs)["optionResults"], BoolList);
    if(B.convertible)
      out.bond = (bool) enum_value(root(outputs)["bondResults"], BoolList);
    out.ratetype = (Compounding)enum_value(root(outputs)["yieldCompounding"], 
					CompoundingList);
    S.su0 << root(stock)["currentPrice"];
    if(B.convertible)
      S.discount << root(stock)["countryDiscountPercent"];
    else 
      S.discount = 0;
    S.discount /= 100;
    double sigma;
    sigma << root(stock)["volatilityPercent"];
    sigma /= 100;
    if(!B.convertible){
      B.x << root("Option")["exercisePrice"];
      B.te << root("Option")["earliestExercise"];
      B.tl << root("Option")["maturity"];
      in.R0 << root("riskFreeRatePercent")["rate"];
      in.R_type = (Compounding) enum_value(root("riskFreeRatePercent")["compounding"],
					CompoundingList);
    }else{
      B.x << root("conversion")["conversionPrice"];
      B.te << root("conversion")["earliestConversion"];
      if (!root("conversion")["maturity"].IsEmpty())
	B.tl << root("conversion")["maturity"];
      in.R0 << root("bond")("riskFreeRatePercent")["rate"];
      in.R_type = (Compounding) enum_value(root("bond")("riskFreeRatePercent")["compounding"],
					CompoundingList);
    }
    in.R0 /= 100;
    B.R = effv_rate(in.R0, in.R_type);
    {
      const Element& elem = root(stock)("dividends");
      S.ndiv = count_children(elem);
      S.Div = new double[S.ndiv];
      S.xD = new double[S.ndiv];
      ConstIterator<Element> it = elem.Begin();
      ConstIterator<Element> itEnd = elem.End();
      for(int i = 0; it < itEnd; ++it, i++){
	S.Div[i] << (*it)["amount"];
	S.xD[i] << (*it)["year"];
      }
      in.q0 << elem["yieldPercent"];
      in.q0 /= 100;
      in.q_type = (Compounding)enum_value(elem["yieldCompounding"], CompoundingList);
    }
    S.q = effv_rate(in.q0/(1-S.discount), in.q_type);
    // q0 is the yield on the undiscounted price
    params.find_vega = enum_value(root(outputs)["vega"], VegaTypeList);
    if(B.convertible){
      in.BondMatYears << root("bond")("maturity")["years"];
      in.BondMatDays << root("bond")("maturity")["days"];
      B.maturity = in.BondMatYears + in.BondMatDays/365.0;
      if (root("conversion")["maturity"].IsEmpty())
	B.tl = B.maturity;
    }else
      B.maturity = B.tl;
    L = *(new MyLattice(N, B.maturity));
    TMptr = L.timemap_ptr();
    LP.set_params(sigma, B.R - S.q);
    if(B.convertible){
      B.coupon << root("bond")("coupon")["rate"];
      B.coupon /= 100;
      B.cpn_freq << root("bond")("coupon")["frequency"];
      in.kd_0 << root("bond")("requiredYieldPercent")["rate"];
      in.kd_0 /= 100;
      in.kd_type = (Compounding)enum_value(root("bond")("requiredYieldPercent")
					["compounding"], CompoundingList);
      B.kd = effv_rate(in.kd_0, in.kd_type);
      B.redeem << root("bond")["redemptionValue"];
      B.redeem /= 100;
      B.bond_price << root("bond")["currentPrice"];
      params.find_ytm = (bool) enum_value(root(outputs)["ytm"], BoolList);
      if(root(outputs)["ytmToler"].IsEmpty())
          params.toler = 1.0/params.default_toler_reciprocal;
      else
          params.toler << root(outputs)["ytmToler"];
      if(root(outputs)["maxIter"].IsEmpty())
          params.max_iter = params.default_max_iter;
      else
          params.max_iter << root(outputs)["maxIter"];
      {
	const Element& elem = root("callProvisions");
	B.set_no(count_children(elem));
	ConstIterator<Element> it = elem.Begin();
	ConstIterator<Element> itEnd = elem.End();
	for(int i = 0; it < itEnd; ++it, i++){
	  double tmin, tmax, price, redeem;
	  tmin << (*it)["minTime"];
	  tmax << (*it)["maxTime"];
	  price << (*it)["minPrice"];
	  redeem << (*it)["redemptionValue"];
	  redeem /= 100;
	  B.co(i).set(tmin, tmax, price, redeem);
	}
	params.ignore_call = (bool) enum_value(elem["ignore"], BoolList);
      }
      {
	const Element& elem = root("putProvisions");
	B.set_np(count_children(elem));
	ConstIterator<Element> it = elem.Begin();
	ConstIterator<Element> itEnd = elem.End();
	for(int i = 0; it < itEnd; ++it, i++){
	  double tmin, tmax, yield;
	  Compounding type;
	  tmin << (*it)["minTime"];
	  tmax << (*it)["maxTime"];
	  yield << (*it)["putValue"];
	  type = (Compounding) enum_value((*it)["valueType"], PutValueTypeList);
	  if(type != price){
	    yield /= 100;
	    yield = effv_rate(yield, type);
	  }
	  B.po(i).set(tmin, tmax, yield, type);
	}
	params.ignore_put = (bool) enum_value(elem["ignore"], BoolList);
      }
      in.sdot_div << root("currencyAppreciation")["dividends"];
      B.sdot_coupon << root("currencyAppreciation")["coupons"];
      B.sdot_redeem << root("currencyAppreciation")["redeem"];
      B.sdot_o_price << root("currencyAppreciation")["callPrice"];
      B.sdot_o_redeem << root("currencyAppreciation")["callRedeem"];
      B.sdot_p_price << root("currencyAppreciation")["putValue"];
      in.sdot_div /= 100;
      B.sdot_coupon /= 100;
      B.sdot_redeem /= 100;
      B.sdot_o_price /= 100;
      B.sdot_o_redeem /= 100;
      B.sdot_p_price /= 100;
      in.sdot_type = (Compounding)enum_value(root("currencyAppreciation")["compounding"], 
					  CompoundingList);
      if(params.ignore_put) B.set_np(0);
      if(params.ignore_call) B.set_no(0);
      B.disc_rate = enum_value(root(computation)["latticeDiscounting"], 
			       DiscountTypeList);
    }else{
      in.sdot_div = B.sdot_coupon = B.sdot_redeem = B.sdot_o_price = B.sdot_p_price
	= B.sdot_o_redeem = 0;
    } // if B.convertible
    in.sdot_div = effv_rate(in.sdot_div, in.sdot_type);
    B.sdot_coupon = effv_rate(B.sdot_coupon, in.sdot_type);
    B.sdot_redeem = effv_rate(B.sdot_redeem, in.sdot_type);
    B.sdot_o_price = effv_rate(B.sdot_o_price, in.sdot_type);
    B.sdot_p_price = effv_rate(B.sdot_p_price, in.sdot_type);
    B.sdot_o_redeem = effv_rate(B.sdot_o_redeem, in.sdot_type);
    S.set_params(in.sdot_div, B.R);
    B.set_params();
    out.ratetype = (Compounding)
      enum_value(root(outputs)["outputInterestRateCompounding"], CompoundingList);
    Element& xmlroot = xmldoc.GetRoot();
    xmlroot["Output_Interest_Rate_Compounding"] << rate_type(pOut->ratetype);
  }
  catch(ParsingException e){
    // Parsing error
    //*(MyErr.msg) << "XML Parsing error at line " << e.GetLine() << " in input file" 
    //		 << endl;
    throw MyErr;
  }
  delete[] pBuffer;
  if(!pOut->silent && pOut->data) 
    print_data(in);
  {
    // MyLattice M(6,6);
    //cerr << "M.n is " << M.n() << endl;
  }
  return L;
}

char *dbl2str(const double x)
  //convert a double into a text string using the global format variable dblformat
{
  static char s[100];
  sprintf(s, dblformat, x);
  return s;
}

char *rate2str(const double x)
  //convert an interest rate into the desired output compounding type pOut->ratetype
  //and then convert it into a text string using dbl2str
{
  double y = convert_rate(x, pOut->ratetype);
  return dbl2str(100.0*y);
}

void print_data(const Input &in)
{
  Element& xmlroot = xmldoc.GetRoot();
  Element& xmldata = xmlroot("Input_Data", 0);
  xmldata["type"] << ((pB->convertible) ? "Call (Convertible)" : ((pB->put) ? "Put" : "Call"));
  {
    Element& xmlouts = xmldata("Outputs", 0);
    xmlouts["Input_Data"] << pOut->data;
    xmlouts["Lattice_Parameters"] << pOut->params;
    xmlouts["Full_Lattice"] << pOut->lattice;
    xmlouts["Option_Results"] << pOut->option;
    if(pB->convertible)xmlouts["Bond_Results"] << pOut->bond;
  }{
    Element& xmlcompute = xmldata("Computational_Parameters", 0);
    xmlcompute["No_of_Lattice_Points"] << in.N0;
    xmlcompute["Vega_Computation"] << VegaTypeList[pParams->find_vega];
    if(pB->convertible){
      xmlcompute["YTM_Computation"] << pParams->find_ytm;
      xmlcompute["YTM_Tolerance"] << dbl2str(pParams->toler);
      xmlcompute["YTM_Iteration_Limit"] << dbl2str(pParams->max_iter);
      xmlcompute["Lattice_Discounting"] << DiscountTypeList[pB->disc_rate];
    }
  }{
    Element& xmlstock = xmldata("Stock", 0);
    xmlstock["Unadjusted_Stock_Price"] << dbl2str(pS->su0);
    if(pS->discount != 0){
      if(pB->convertible)
	xmlstock["Country-Discount_on_Stock_Price_Percent"] << dbl2str(pS->discount*100);
      xmlstock["Discounted_Stock_Price"] << dbl2str(pS->su);
    }
    xmlstock["Dividend_Yield_Percent"] << rate2str(pS->q);
    Element& Kdiv = xmlstock("Known_Dividends", 0);
    for (int i = 0; i < pS->ndiv; i++){
      Kdiv("Dividend", i)["Amount"] << pS->Div[i];
      Kdiv("Dividend", i)["Date"] << pS->xD[i];
    }
  }{
    Element& xmloption = xmldata("Option", 0);
    xmloption["Exercise_Price"] << dbl2str(pB->x);
    xmloption["Earliest_Exercise"] << dbl2str(pB->te);
    xmloption["Latest_Exercise"] << dbl2str(pB->tl);
    xmloption["Volatility_Percent"] << dbl2str(pLP->sigma*100);
    xmloption["Riskfree_Rate_Percent"] << rate2str(pB->R);
  }
  if(pB->convertible){
    {
      Element& xmlbond = xmldata("Bond");
      xmlbond["Residual_Bond_Maturity_in_Years"] << dbl2str(pB->maturity);
      xmlbond["Coupon_Rate_Percent"] << dbl2str(pB->coupon*100);
      xmlbond["Coupon_Frequency_per_Year"] << pB->cpn_freq;
      xmlbond["Required_Yield_Percent"] << rate2str(pB->kd);
      xmlbond["Redemption_Value_Percent_of_Face_Value"] << dbl2str(pB->redeem*100);
      if(pParams->find_ytm)
	xmlbond["Bond_Price"] << dbl2str(pB->bond_price);
    }{
      Element& xmlcall = xmldata("Call_Provisions");
      if(pParams->ignore_call)
	xmlcall["Options_to_issuer_ignored"] << "true";
      else
	for (int i = 0; i < pB->no(); i++){
	  Element& temp = xmlcall("Call", i);
	  temp["Min_Time"] << dbl2str(pB->co(i).Time_min());
	  temp["Max_Time"] << dbl2str(pB->co(i).Time_max());
	  temp["Price"] << dbl2str(pB->co(i).Price());
	  temp["Redeem"] << dbl2str(pB->co(i).Redeem()*100);
	}
    }{
      Element& xmlput = xmldata("Put_Provisions");
      if(pParams->ignore_put)
	xmlput["Put_options_ignored"] << "true";
      else
	for (int i = 0; i < pB->np(); i++){
	  Element& temp = xmlput("Put", i);
	  temp["Min_Time"] << dbl2str(pB->po(i).Time_min());
	  temp["Max_Time"] << dbl2str(pB->po(i).Time_max());
	  if(pB->po(i).Type() == price)
	    temp["Price"] << dbl2str(pB->po(i).Yield());
	  else
	    temp["Yield_to_Put_Percent"] << rate2str(pB->po(i).Yield());
	}
    }{
      Element& xmlcurrency = xmldata("Expected_Currency_Appreciation", 0);
      xmlcurrency["Dividends"] << rate2str(in.sdot_div);
      xmlcurrency["Coupons"] << rate2str(pB->sdot_coupon);
      xmlcurrency["Redemption"] << rate2str(pB->sdot_redeem);
      xmlcurrency["Call_Price_Limits"] << rate2str(pB->sdot_o_price);
      xmlcurrency["Call_Redemption_Values"] << rate2str(pB->sdot_o_redeem);
      xmlcurrency["Put_Prices"] << rate2str(pB->sdot_p_price);
    }
  } // if(pB->convertible)
}

void print_params()
{
  Element& xmlroot = xmldoc.GetRoot();
  Element& xmlparams = xmlroot("Lattice_Paramaters", 0);
  xmlparams["No_of_Lattice_Points"] << TMptr->N();
  xmlparams["Present_Value_of_Discrete_Dividends"] << dbl2str(pS->pvdiv(0, pB->R));
  xmlparams["Stock_Price_Adjusted_for_Dividends"] << dbl2str(pS->s);
  xmlparams["Continously_Compounded_Riskfree_Rate"] << dbl2str(pB->R);
  xmlparams["Continously_Compounded_Dividend_Rate"] << dbl2str(pS->q);
  xmlparams["Continously_Compounded_Net_Growth_Rate"] << dbl2str(pLP->RminusQ);
  xmlparams["u"] << dbl2str(pLP->u);
  xmlparams["d"] << dbl2str(pLP->d); 
  xmlparams["p"] << dbl2str(pLP->p);
}


void print_lattice(const MyLattice& L)
{
  static Element *pTree;
  static char s[2] = "";

  int n = L.n();
  int N = TMptr->N();
  if(n == N){
    Element& xmlroot = xmldoc.GetRoot();
    pTree = &(xmlroot("Full_Lattice", 0));
    if(pB->convertible)
      (*pTree)["Bond_Redemption_Value"] << dbl2str(pB->x*pB->s_redeem);
  }
  double dt = TMptr->dt();
  double t = TMptr->t(n);
  Element& stage = (*pTree)("Stage", N-n);
  stage["Step"] << n;
  stage["Time"] << dbl2str(t);
  if (pS->ndiv)
    stage["Present_Value_of_Dividends"] << dbl2str(pS->pvdiv(t, pB->R));
  if(pB->convertible)
    stage["Coupons"] << dbl2str(pB->x*pB->cpn(n));
  for (int i = 0; i <= n; i++){
    Element& node = stage("Node", i);
    node["S"] <<  dbl2str(L[i].s);
    if (pB->convertible){
      node["VB"] << L[i].vB;
      node["VO"] << L[i].vO;
    }
    node["V"] << L[i].v;
    switch(L[i].ch){
    case 'P' : node["Action"] << "Put"; break;
    case 'E' : node["Action"] << "Exercised"; break;
    case 'C' : node["Action"] << "Called"; break;
    }
  }
}

void print_results(const Results& res)
{
  double t = 0;
  double dt = TMptr->dt();

  Element& xmlroot = xmldoc.GetRoot();
  if (!pB->convertible){
    Element& xmlres = xmlroot("Option_Results", 0);
    xmlres["Option_Value"] << dbl2str(res.f);
    xmlres["Delta"] << dbl2str(res.delta);
    xmlres["Gamma"] << dbl2str(res.Gamma);
    xmlres["Theta_Indirect"] << dbl2str(res.theta3);
    xmlres["Theta_Direct_2dt"] << dbl2str(res.theta1);
    xmlres["Theta_Direct_4dt"] << dbl2str(res.theta2);
    if (pParams->find_vega) 
      xmlres["Vega"] << dbl2str(res.vega);
    return;
  }
  if(pOut->option){
    Element& xmlres = xmlroot("Option_Results", 0);
    xmlres["Bond_Convertible_into_One_Share_Face_Value"] << dbl2str(res.x);
    xmlres["Bond_Convertible_into_One_Share_Straight_Value"] << dbl2str(res.straight*res.x/100);
    xmlres["Bond_Convertible_into_One_Share_Option_Value"] << dbl2str(res.optval*res.x/100);
    xmlres["Delta"] << dbl2str(res.delta);
    xmlres["Gamma"] << dbl2str(res.Gamma);
    if (pParams->find_vega) 
      xmlres["Vega"] << dbl2str(res.vega);
  }
  Element& xmlres = xmlroot("Bond_Value_Results", 0);
  xmlres["Bond_Value"] << dbl2str(res.bondval);
  xmlres["Straight_Bond_Value"] << dbl2str(res.straight);
  xmlres["Option_Value"] << dbl2str(res.optval);
  if (pParams->find_ytm){
    Element& xmlytm = xmlroot("YTM", 0);
    xmlytm["Bond_Price"] << dbl2str(res.price);
    if(strlen(res.YTMerrmsg) != 0) {
      xmlytm["Error_Message"] << res.YTMerrmsg;
      xmlytm["Iteration_No."] << res.YTMiter;
      xmlytm["Bond_Value"] << dbl2str(res.bondval);
      xmlytm["Error"] << dbl2str(res.bondval - res.price);
    }
    if(res.ytm < res.rf)
      xmlytm["YTM_Meaningless"] << "true";
    xmlytm["YTM_Percent"] << rate2str(res.ytm);
    int oas = (int) (10000*(convert_rate(res.ytm, pOut->ratetype) - 
			    convert_rate(res.rf, pOut->ratetype)) + 0.5);
    xmlytm["OAS_Basis_Points"] << oas;
  }
}

FILE *getfile(int argc, char **argv)
{

  FILE *infile;

  if (argc < 2){
    // MyErrMsg << "No Input File Specified"  << endl;
    return stdin;
  }
  if (NULL == (infile = fopen(argv[1],"rb"))){
    MyErrMsg << "Unable to open " << argv[1] << endl;
    throw MyErr;
  }
  return infile;
}
