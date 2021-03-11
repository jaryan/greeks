#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

SEXP BSput (SEXP S, SEXP X, SEXP b, SEXP r, SEXP t, SEXP v)
{

  /* possible add vectorization to this to make it faster, though
     the R code is plenty fast enough */

  int i = 0;
  /*for(i=0; i < length(S); i++) {*/

  double _S = REAL(S)[i];
  double _X = REAL(X)[i];
  double _b = REAL(b)[i];
  double _r = REAL(r)[i];
  double _t = REAL(t)[i];
  double _v = REAL(v)[i];

  double d1, d2;
  double put;
  double delta, gamma, vega, theta, rho;
  double vanna, charm, speed, zomma;
  double colour, DvegaDtime, vomma, dualdelta, dualgamma;

  SEXP ans;

  d1 = (log(_S/_X) + (_r - _b + (_v*_v)/2.0) * _t)/(_v * sqrt(_t));
  d2 = d1 - _v * (sqrt(_t));

  put  = exp(-_r * _t) * _X * pnorm(-d2,0,1,1,0) -
         exp(-_b * _t) * _S * pnorm(-d1,0,1,1,0);

  delta = -exp(-_b * _t) * pnorm(-d1,0,1,1,0);
  vega  = _S * exp(-_b * _t) * dnorm(d1,0,1,0) * sqrt(_t);
  theta = -exp(-_b * _t) * (_S * dnorm(d1,0,1,0) * _v)/
          (2.0 * sqrt(_t)) + (_r * _X * exp(-_r * _t) * pnorm(-d2,0,1,1,0));
  rho   = -_X * _t * exp(-_r * _t) * pnorm(-d2,0,1,1,0);
  gamma = exp(-_b * _t) * ( dnorm(d1,0,1,0)/( _S * _v * sqrt(_t)) );
  vanna = -exp(-_b * _t) * dnorm(d1,0,1,0) * (d2/_v);
  charm = _b * exp(-_b * _t) * pnorm(-d1,0,1,1,0) +
          exp(-_b * _t) * dnorm(d1,0,1,0) *
          ( (2.0 * (_r - _b) * _t - d2 * _v * sqrt(_t))/(2.0 * _t * _v * sqrt(_t)) );
  speed = -(gamma/_S) * (d1/(_v * sqrt(_t)) + 1.0);
  zomma = gamma * ( (d1*d2 - 1.0)/_v );
  colour = -exp(-_b * _t) * (dnorm(d1,0,1,0)/(2.0*_S*_t*_v*sqrt(_t))) *
           (2.0 * _b * _t + 1.0 +
             ((2.0 * (_r-_b) * _t - d2 * _v * sqrt(_t))/
               (_v * sqrt(_t))) * d1);
  DvegaDtime = _S * exp(-_b * _t) * dnorm(d1,0,1,0) * sqrt(_t) *
               (_b + ((_r-_b)*d1)/(_v*sqrt(_t))-((1.0+d1*d2)/(2.0*_t)));
  vomma = vega * ((d1*d2) / _v);
  dualdelta = exp(-_r * _t) * pnorm(-d2,0,1,1,0);
  dualgamma = exp(-_r * _t) * ( dnorm(d2,0,1,0)/(_X*_v*sqrt(_t)) );

  PROTECT(ans = allocVector(REALSXP, 15));
  double *_ans = REAL(ans);

  _ans[0] = put;
  _ans[1] = delta;
  _ans[2] = gamma;
  _ans[3] = vega;
  _ans[4] = theta;
  _ans[5] = rho;
  _ans[6] = vanna;
  _ans[7] = charm;
  _ans[8] = zomma;
  _ans[9] = speed;
  _ans[10] = colour;
  _ans[11] = DvegaDtime;
  _ans[12] = vomma;
  _ans[13] = dualdelta;
  _ans[14] = dualgamma;
  /*
  ans = coerceVector(ans, VECSXP);
  */


  SEXP GreeksNames = allocVector(STRSXP,15);
  PROTECT(GreeksNames);
  SET_STRING_ELT(GreeksNames,0,mkChar("value"));
  SET_STRING_ELT(GreeksNames,1,mkChar("delta"));
  SET_STRING_ELT(GreeksNames,2,mkChar("gamma"));
  SET_STRING_ELT(GreeksNames,3,mkChar("vega"));
  SET_STRING_ELT(GreeksNames,4,mkChar("theta"));
  SET_STRING_ELT(GreeksNames,5,mkChar("rho"));
  SET_STRING_ELT(GreeksNames,6,mkChar("vanna"));
  SET_STRING_ELT(GreeksNames,7,mkChar("charm"));
  SET_STRING_ELT(GreeksNames,8,mkChar("zomma"));
  SET_STRING_ELT(GreeksNames,9,mkChar("speed"));
  SET_STRING_ELT(GreeksNames,10,mkChar("colour"));
  SET_STRING_ELT(GreeksNames,11,mkChar("DvegaDtime"));
  SET_STRING_ELT(GreeksNames,12,mkChar("vomma"));
  SET_STRING_ELT(GreeksNames,13,mkChar("dualdelta"));
  SET_STRING_ELT(GreeksNames,14,mkChar("dualgamma"));
  setAttrib(ans, R_NamesSymbol, GreeksNames);
  setAttrib(ans, R_ClassSymbol, mkString("greeks"));

  UNPROTECT(2);
  return(ans);
}

SEXP BScall (SEXP S, SEXP X, SEXP b, SEXP r, SEXP t, SEXP v)
{
  double _S = REAL(S)[0];
  double _X = REAL(X)[0];
  double _b = REAL(b)[0];
  double _r = REAL(r)[0];
  double _t = REAL(t)[0];
  double _v = REAL(v)[0];

  double d1, d2, call;
  SEXP ans;

  d1 = (log(_S/_X) + (_r - _b + (_v*_v)/2.0) * _t)/(_v * sqrt(_t));
  d2 = d1 - _v * (sqrt(_t));

  call = (_S * exp(-_b * _t) * pnorm(d1,0,1,1,0)) -
         (_X * exp(-_r * _t) * pnorm(d2,0,1,1,0));

  PROTECT(ans = allocVector(REALSXP, 1));
  REAL(ans)[0] = call;
  UNPROTECT(1);
  return(ans);
}
