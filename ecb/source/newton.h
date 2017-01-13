/*
    newton.h, part of ecb, a program for Valuation of Convertible Bonds
    Copyright (C) 2001  Prof. Jayanth R. Varma, jrvarma@iimahd.ernet.in,
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

  static const int numdiff_central = 0;
  static const int numdiff_forward = 1;
  static const int numdiff_backward = -1;

template <class Z> class numdiff 
// This class is used to compute the numerical derivative of a function f with respect to x
// i.e., to compute d/dx of f(x,Z).
// The function f could depend on several other variables (Z) and may involve a very complex
// computational process. Z may be a simple structure containing Z1, Z2, etc, but, in general, 
// Z can be an arbitrary class.
// The easiest way to use numdiff for a complicated function is to create a wrapper function 
// f(x, Z) in which all other variables are encapsulated in the class Z. 
// The wrapper function could set global variables using the values of members of
// Z. It may call other subroutines using members of Z as arguments. It may carry out complex 
// computations. But at the end, the wrapper function must return f(x, Z)
{
 public:
  static const int central = 0;
  static const int forward = 1;
  static const int backward = -1;
  //parameters that can be set by user
  double c1;  // parameter to set delta x as max(c1, c2*x);
  double c2;  // parameter to set delta x as max(c1, c2*x);
  int type;   // whether central, forward or backward differences
  //pointer to function f(x,Z)
  double (*fp)(const double, Z&);
  //pointer to Z
  Z *pZ;
  //constructor requires pointer to f and Z
  //other parameters are set to default values
  numdiff (double (*fp0)(const double, Z&), Z *p) : 
    fp(fp0), pZ(p) {type = central; c1 = c2 = 0.001;};
  //method deriv computes the numerical derivative
  double deriv(const double f, const double x);
};


template <class Z> double numdiff<Z>::deriv(const double f, const double x)
{
  double dx = max(c1, c2*x);
  switch (type){
  case central : {
    double f1 = (*fp)(x-0.5*dx, *pZ);
    double f2 = (*fp)(x+0.5*dx, *pZ);
    return (f2 - f1) / dx;
  }
  case forward : {
    double f2 = (*fp)(x+dx, *pZ);
    return (f2 - f) / dx;
  }
  case backward : {
    double f1 = (*fp)(x-dx, *pZ);
    return (f - f1) / dx;
  }
  default : {  // this is only to prevent a compiler warning
    double f2 = (*fp)(x+dx, *pZ);
    return (f2 - f) / dx;
  }
  }
}

template <class Z> class newton_solver
// This class is used to solve for x* such that f(x,Z) equals a target value.
// i.e., to find x* such that f(x*,Z)=target.
// x may be subject to lower and upper bounds
// The function f could depend on several other variables (Z) and may involve a very complex
// computational process. Z may be a simple structure containing Z1, Z2, etc, but, in general, 
// Z can be an arbitrary class.
// The easiest way to use newton_solver for a complicated function is to create a wrapper 
// function f(x, Z) in which all other variables are encapsulated in the class Z. 
// The wrapper function could set global variables using the values of members of
// Z. It may call other subroutines using members of Z as arguments. It may carry out complex 
// computations. But at the end, the wrapper function must return f(x, Z)
// If the user can provide an analytic derivative of f, newton_solver can make use of it
// Otherwise, it uses numdiff to compute a numerical derivative.
{
 public:
  //the following are part of the problem definition
  double target;          // the target value f(x,Z) = target (default is 0)
  double xmax;            // upper bound on x (default is 1e30 or approx infinity)
  double xmin;            // lower bound on x (default is 1e30 or approx -infinity)
  // the following parameters can be set by user to finetune the algorithm
  int max_iter;           // maximum number of iterations
  double tolerance;       // tolerance i.e. |f(x,Z) - target| < tolerance
  double step_limit;      // used to keep x away from its lower and upper bound
  double step_reduction;  // factor by which newton step is reduced to prevent uphill step
  double min_improve;     // used in definition of uphill direction
  int max_step_reductions;// number of times step is reduced to avoid uphill direction
  int numdiff_type;       // central, forward or backward differences in numdiff
  double c1;              // parameter to set delta x as max(c1, c2*x) in numdiff;
  double c2;              // parameter to set delta x as max(c1, c2*x) in numdiff;
  int debug;              // controls printing out intermediate results for debugging
  //these are outputs
  int iter_no(void)       // the number of iterations
    {return iter;};
  const char *err_msg()         // returns error message if any
    {return _err_msg;}
  //constructor for analytic differentiation
  //requires pointer to f(x,Z), g=d/dx of f(x,Z) and Z
  newton_solver(double (*fp0)(const double, Z&), double (*gp0)(const double, Z&), Z *p)
    : fp(fp0), gp(gp0), pZ(p) {numerical = false; init();};
  //constructor for numerical differentiation
  //requires pointer to f(x,Z) and Z
  newton_solver(double (*fp0)(const double, Z&), Z *p)
    : fp(fp0), pZ(p) {numerical = true; init();};
  //this is the solver. Requires starting point x0 and f0=f(x0,Z)
  double solve(const double x0, const double f0);
	
 private:
  double (*fp)(const double, Z&);
  double (*gp)(const double, Z&);
  double eval_error(const double x){return (*fp)(x, *pZ) - target;}  
  double eval_g(const double x){return (*gp)(x, *pZ);}
  Z *pZ;
  bool numerical;
  double error;
  double x;
  double g;
  int iter;
  void init();
  const char *_err_msg;
};

template <class Z> void newton_solver<Z>::init()
     //set defaults parameter values
{
  max_iter = 100;
  tolerance = 1e-6;
  target = 0;
  step_limit = 0.9;
  step_reduction = 0.5;
  min_improve = 0.5;
  xmax = 1e30;
  xmin = -1e30;
  max_step_reductions = 10;
  debug = 0;
  numdiff_type = numdiff<Z>::central;
  c1 = 0.001;
  c2 = 0.001;
}

template <class Z> double newton_solver<Z>::solve (const double x0, const double f0){
  numdiff<Z> nd(fp, pZ);
  nd.c1 = c1;
  nd.c2 = c2;
  nd.type = numdiff_type;
  _err_msg="";
  x = x0;
  double error_new = error = f0 - target; 
  iter = 0;
  while (iter <= max_iter && fabs(error) > tolerance){ // main loop
    if(debug > 0)
      cout << "Iteration " << iter << ": error =" << error << ".  x = " << x << endl;
    //compute the gradient g
    if (numerical)
      g = nd.deriv(target+error, x);
    else
      g = eval_g(x);
    double step = -error/g; // newton step
    //reduce the step to keep x away from its lower and upper bounds
    if(step < step_limit*(xmin-x)) step = step_limit*(xmin-x);
    if(step > step_limit*(xmax-x)) step = step_limit*(xmax-x);
    if(debug > 1)
      cout << "Error= " << error << ". g=" << g << ". Try step=" << step << endl;
    int count = 0;
    do{
      if(debug > 1 && count > 0)
	cout << "Failed. Reduction in error = " << fabs(error) - fabs(error_new) 
	     << ". Reqd = " << min_improve*step*g << endl;
      if (count > 0)
        step *= step_reduction;
      error_new = eval_error(x+step); 
      if(debug > 1)
	cout << "x = " << x+step << ". f = " << error_new+target << endl;
      count++;
    } while (fabs(error) - fabs(error_new) < min_improve*step*g && count <= max_step_reductions);
    //predicted improvement is step times gradient
    //if actual improvement is much below this, the step is too large
    if(count > max_step_reductions){
      //the gradient appears to point uphill
      _err_msg = "Uphill gradient";
      return (x);
    }
    iter++;
    x += step;
    error = error_new;
  }
  if(fabs(error) > tolerance){
    _err_msg = "Iteration count exceeded";
    return(x);
  }
  return x;
}

