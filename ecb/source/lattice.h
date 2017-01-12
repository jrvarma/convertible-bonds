/*
    lattice.h, part of ecb, a program for Valuation of Convertible Bonds
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
class timemap
// the purpose of this class is to convert from lattice time to real time and back
// this implements a simple conversion which simply rescales by the factor dt
// real_time t = (lattice_time n) * dt
{
 private:
  int _N;
  double _T, _dt;
    
 public:
  virtual double t(const int n) const {return n*_dt;}
  //convert lattice time to real time
  virtual int n(const double t) const 
    {return (int) ( (floor(t/_dt+0.5) > _N) ?  _N : floor(t/_dt+0.5));} 
  //{return (int) ( (rint(t/_dt) > _N) ?  _N : rint(t  /_dt));} 
  //convert real time to lattice time
  int N(void) const {return _N;}        // maximum lattice time
  double T(void) const {return _T;}     // maximum real time
  double dt(void) const {return _dt;}   // scaling factor
  timemap(void) : _T(0), _N(0), _dt(1) {}
  timemap(const timemap &x) : _T(x._T), _N(x._N), _dt(x._dt){}
  timemap(const int N, const double T) {_T = T; _N = N; _dt = _T/_N;}
};

class bad_subscript
{
 public:
  int low, high, actual;
  bad_subscript(int l, int h, int a): low(l), high(h), actual(a) {}
};

template <class X, class M> class lattice
// this template depends on two classes
// the node class X
// and the timemap class M
{
 private:
  X *v1, *v2, *v3; // three arrays of nodes
  int _n;          // current lattice time
  M *pTM;          // timemap
  void init();
    
 public:
  lattice(void) : _n(0) {pTM = new timemap();} // empty lattice
  lattice(const int N, const double T) {pTM = new timemap(N, T); init();}
  //lattice time runs from 0-N and real time from 0-T
  lattice(const lattice& x);
  int n(void) const {return _n;}
  M *timemap_ptr(void) const {return pTM;}
  int rewind() {_n = pTM->N();}
  void swap() {X *temp = v1; v1 = v3; v3 = v2; v2 = temp; if (_n > 0) _n--;}
  void redim(const int N, const double T);
  X& operator [] (int i) const 
    {if (0 <= i && i <= pTM->N()) return v1[i]; else throw(bad_subscript(0,pTM->N(),i));}
  X& nPlus1(int i) const 
    {if (0 <= i && i <= pTM->N()) return v2[i]; else throw(bad_subscript(0,pTM->N(),i));}
  X& nPlus2(int i) const 
    {if (0 <= i && i <= pTM->N()) return v3[i]; else throw(bad_subscript(0,pTM->N(),i));}
  X& upmove(int i) const 
    {if (0 <= i && i <= pTM->N()) return v2[i]; else throw(bad_subscript(0,pTM->N(),i));}
  X& dnmove(int i) const 
    {if (0 <= i && i <  pTM->N()) return v2[i+1]; else throw(bad_subscript(0,pTM->N(),i));}
};

template <class X, class M> void lattice<X, M>::init()
{
  _n = pTM->N();
  if (_n > 0){
    v1 = new X[_n+1];
    v2 = new X[_n+1];
    v3 = new X[_n+1];
  }
}   

template <class X, class M> lattice<X, M>::lattice(const lattice& x) 
{
  pTM = new timemap(*(x.pTM));
  init();
  for (int i = 0; i <= pTM->N(); i++) {
    v1[i] = x.v1[i];
    v2[i] = x.v2[i];
    v3[i] = x.v3[i];
  }
}

template <class X, class M> void lattice<X, M>::redim(const int N, const double T)
{
  if (pTM->N() > 0){
    delete [] v1;
    delete [] v2;
    delete [] v3;
  }
  delete pTM;
  pTM = new timemap(N, T);
  init(n);
}    
