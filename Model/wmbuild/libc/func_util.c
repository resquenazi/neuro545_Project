/*****************************************************************************/
/*                                                                           */
/*  func_util.c                                                              */
/*  wyeth bair                                                               */
/*  caltech                                                                  */
/*  02/03/94                                                                 */
/*                                                                           */
/*  For miscellaneous functions, and 1d sampled arrays.                      */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <math.h>

#include "my_util.h"

/**************************************-**************************************/
/*                                                                           */
/*                                 FUNC_ALPHA                                */
/*                                                                           */
/*  Alpha function.                                                          */
/*                                                                           */
/*  - Peak value is 1.0.                                                     */
/*  - Time to peak is t0 + 1/alpha.                                          */
/*  - Used to model the time course of post-synaptic currents.               */
/*                                                                           */
/*****************************************************************************/
double func_alpha(x,t0,alpha)
     double x,t0,alpha;
{
  if (x < t0)
    return 0.0;
  else
    return (alpha*exp(1.0) * (x-t0) * exp(-alpha*(x-t0)));
}
/**************************************-**************************************/
/*                                                                           */
/*                                 FUNC_BOXCAR                               */
/*                                                                           */
/*  Boxcar function.                                                         */
/*  - Value 1.0 from [t0,t0+width).                                          */
/*                                                                           */
/*****************************************************************************/
double func_boxcar(x,t0,width)
     double x,t0,width;
{
  if ((x >= t0)&&(x < (t0+width)))
    return 1.0;
  else
    return 0.0;
}
/**************************************-**************************************/
/*                                                                           */
/*                                FUNC_DIFF_EXP                              */
/*                                                                           */
/*  Difference of exponentials (DOE).                                        */
/*  tau1 - slow, falling phase                                               */
/*  tau2 - fast, rising phase                                                */
/*                                                                           */
/*****************************************************************************/
double func_diff_exp(x,t0,tau1,tau2)
     double x,t0,tau1,tau2;
{
  if (x >= t0)
    return exp((t0-x)/tau1) - exp((t0-x)/tau2);
  else
    return 0.0;
}
/**************************************-**************************************/
/*                                                                           */
/*                                FUNC_GAMMA_N3                              */
/*                                                                           */
/*  Gamma function, n=3. (p78 Papoulis, Prob, RV & Stoch. Proc. 3rd Ed.)     */
/*                                                                           */
/*****************************************************************************/
double func_gamma_n3(x,s)
     double x,s;
{
  if (x < 0.0)
    return 0.0;
  else
    return (x*x*x/(s*s*s) * exp(-x/s));
}
/**************************************-**************************************/
/*                                                                           */
/*                                FUNC_INV_GAMMA                             */
/*                                                                           */
/*                               by Blake Richards                           */
/*                                                                           */
/*  Inverse gamma distribution.                                              */
/*  Important note: due to some weird historical quirk, the function gamma() */
/*  in glibc was actually the log of the gamma function. In c99, they put in */
/*  a function tgamma() (true gamma). Originally I tried to use this, but it */
/*  doesn't work without '-std=c99' in the compile command. I put this in    */
/*  makelib.ss, but that threw up a bunch of warnings and errors. Thus, I    */
/*  decided to simply take the exponent of the "gamma" function to get the   */
/*  function I wanted. - Blake, 08/02/07                                     */
/*                                                                           */
/*****************************************************************************/
double func_inv_gamma(x,alpha,beta)
     double x,alpha,beta;
{
  exit_error("FUNC_INV_GAMMA","gamma() is deprecated?");
  /*
  if (x <= 0.0)
    return 0.0;
  else
    return ((pow(beta,alpha)/exp(gamma(alpha)))*pow(x,-alpha-1)*exp(-beta/x));
  */
  return 0.0;
}

/**************************************-**************************************/
/*                                                                           */
/*                                FUNC_GAUSSIAN                              */
/*                                                                           */
/*  Gaussian distribution.                                                   */
/*                                                                           */
/*  - Area is 1.                                                             */
/*  - Mean and peak at mu.                                                   */
/*  - Variance is sigma^2.                                                   */
/*                                                                           */
/*****************************************************************************/
double func_gaussian(x,mu,sigma)
     double x,mu,sigma;
{
  return (1.0/(sqrt(2.0*M_PI)*sigma)*exp(-0.5*((x-mu)/sigma)*((x-mu)/sigma)));
}
/**************************************-**************************************/
/*                                                                           */
/*                              FUNC_GAUSSIAN_ONE                            */
/*                                                                           */
/*  Gaussian shaped function.                                                */
/*                                                                           */
/*  - Max is 1.                                                              */
/*  - Mean and peak at mu.                                                   */
/*  - Variance is sigma^2.                                                   */
/*                                                                           */
/*****************************************************************************/
double func_gaussian_one(x,mu,sigma)
     double x,mu,sigma;
{
  return exp(-0.5*((x-mu)/sigma)*((x-mu)/sigma));
}
/**************************************-**************************************/
/*                                                                           */
/*                               FUNC_2D_GAUSSIAN                            */
/*                                                                           */
/*  2D Gaussian distribution.                                                */
/*                                                                           */
/*  - Area is 1 if 'normflag' = 1  (Wyeth - check this)                      */
/*  - Mean and peak at m1,m2.                                                */
/*  - Variance is s1^2 and s2^2.                                             */
/*                                                                           */
/*****************************************************************************/
double func_2d_gaussian(x,y,m1,m2,s1,s2,normflag)
     double x,y,m1,m2,s1,s2;
     int normflag;
{
  if (normflag)
    return (1.0/(2.0*M_PI*s1*s2) *
	    exp(-0.5*(((x-m1)/s1)*((x-m1)/s1) + ((y-m2)/s2)*((y-m2)/s2))));
  else
    return exp(-0.5*(((x-m1)/s1)*((x-m1)/s1) + ((y-m2)/s2)*((y-m2)/s2)));
}
/**************************************-**************************************/
/*                                                                           */
/*                                FUNC_LOGISTIC                              */
/*                                                                           */
/*  Logistic function.                                                       */
/*                                                                           */
/*  - 0 at -inf, 1 at +inf, 0.5 at 0.                                        */
/*                                                                           */
/*****************************************************************************/
double func_logistic(x,beta)
     double x,beta;
{
  return (1.0/(1.0 + exp(-2.0*beta*x)));
}
/**************************************-**************************************/
/*                                                                           */
/*                                FUNC_MAXWELL                               */
/*                                                                           */
/*  Maxwell distribution (p78 Papoulis, Prob, RV & Stoch. Proc. 3rd Ed.)     */
/*                                                                           */
/*  Peaks at sqrt(2)*s.                                                      */
/*                                                                           */
/*****************************************************************************/
double func_maxwell(x,s)
     double x,s;
{
  if (x < 0.0)
    return 0.0;
  else
    return ((sqrt(2.0/M_PI)/s) * x*x/(s*s) * exp(-0.5*x*x/(s*s)));
}
/**************************************-**************************************/
/*                                                                           */
/*                             FUNC_PHOTORECEPTOR                            */
/*                                                                           */
/*****************************************************************************/
double func_photoreceptor(x,s)
     double x,s;
{
  if (x < 0.0)
    return 0.0;
  else
    return ((sqrt(2.0/M_PI)/s) * x*x*x*x/(s*s*s*s) * exp(-0.5*x/(s)));
}
/**************************************-**************************************/
/*                                                                           */
/*                          FUNC_PRIMATE_PHOTOCURRENT                        */
/*                                                                           */
/*  From Schnapf et al. 1990, Eqn. 7.                                        */
/*                                                                           */
/*****************************************************************************/
double func_primate_photocurrent(x,tau_r,tau_d,tau_p,phi)
     double x,tau_r,tau_d,tau_p,phi;
{
  double tr,tr3,td,tp,phirad;

  tr = x/tau_r;
  tr3 = tr*tr*tr;
  td = x/tau_d;
  tp = x/tau_p;
  phirad = phi * M_PI/180.0;

  if (x < 0.0)
    return 0.0;
  else
    return tr3/(1.0+tr3) * exp(-td*td) * cos(2.0*M_PI*tp + phirad);
}
/**************************************-**************************************/
/*                                                                           */
/*                                 FUNC_QUICK                                */
/*                                                                           */
/*  sigmoidal function.                                                      */
/*                                                                           */
/*  The "Quick" function is the integral of a Weibull distribution.  We      */
/*  often use this to fit psychometric functions.                            */
/*                                                                           */
/*    p(x) = d - (d-c) exp (-(x/a)^b)                                        */
/*                                                                           */
/*  where                                                                    */
/*    p - performance                                                        */
/*    a - threshold                                                          */
/*    b - slope                                                              */
/*    c - guessing probability                                               */
/*    d - asymptotic probability                                             */
/*                                                                           */
/*****************************************************************************/
double func_quick(a,b,c,d,x)
     double a,b,c,d,x;
{
  return (d-(d-c)*exp(-pow((x/a),b)));
}
/**************************************-**************************************/
/*                                                                           */
/*                                FUNC_SIGMOID                               */
/*                                                                           */
/*   A rather slow sigmoid in the tails.                                     */
/*                                                                           */
/*   x/(1+|x|)                                                               */
/*                                                                           */
/*****************************************************************************/
double func_sigmoid(x)
     double x;
{
  if (x >= 0.0)
    return x/(1.0 + x);
  else
    return x/(1.0 - x);
}
/**************************************-**************************************/
/*                                                                           */
/*                              FUNC_SIGMOID_LOG                             */
/*                                                                           */
/*                              by Blake Richards                            */
/*                                                                           */
/*   A logistic sigmoid function.                                            */
/*                                                                           */
/*   s*(1/(1 + exp(-(x + h)/T)) - 0.5) + m                                   */
/*                                                                           */
/*****************************************************************************/
double func_sigmoid_log(x,s,m,T,h)
     double x;
     double s;
     double m;
     double T;
     double h;
{
  float result;

  result = s*(1.0/(1.0+exp(-(x+h)/T)) - 0.5) + m;
  
  // to prevent nan being returned
  if (isnan(result) && x > 0)
    return s-m;
  else if (isnan(result) && x < 0)
    return -m; 
  else
    return result;
}
/**************************************-**************************************/
/*                                                                           */
/*                              FUNC_SIGMOID_TANH                            */
/*                                                                           */
/*  A rather straight sigmoid - linear in a broad range in the center,       */
/*  then bends quickly at the asymptotes.                                    */
/*                                                                           */
/*****************************************************************************/
double func_sigmoid_tanh(x,c,mu)
     double x,c,mu;
{
  float result,ep,en;

  ep = exp(c*(x-mu)/2.0);
  en = exp(-c*(x-mu)/2.0);
  result = (ep - en)/(ep + en);
  
  // Blake edited this, to prevent nan being returned
  if (isnan(result) && x > 0)
    return 1.0;
  else if (isnan(result) && x < 0)
    return -1.0; 
  else
    return result;
}
/**************************************-**************************************/
/*                                                                           */
/*                              FUNC_ONE_SIDED_EXP                           */
/*                                                                           */
/*  - Peak is 1 at t0.                                                       */
/*                                                                           */
/*****************************************************************************/
double func_one_sided_exp(x,t0,tau)
     double x,t0,tau;
{
  if (x < t0)
    return 0.0;
  else
    return (exp(-(x-t0)/tau));
}
/**************************************-**************************************/
/*                                                                           */
/*                               FUNC_LOGNORMAL                              */
/*                                                                           */
/*  Lognormal function. (p78 Papoulis, Prob, RV & Stoch. Proc. 3rd Ed.)      */
/*                                                                           */
/*****************************************************************************/
double func_lognormal(x,s)
     double x,s;
{
  double xos;

  xos = x/s;
  if (xos <= 0.0)
    return 0.0;
  else
    return 1.0/xos * exp(-0.5 * log(xos)*log(xos));
}
/**************************************-**************************************/
/*                                                                           */
/*                                ALPHA_FARRAY                               */
/*                                                                           */
/*****************************************************************************/
float *alpha_farray(t0,alpha,ampl,n)
     float t0,alpha,ampl;
     int n;
{
  int i;
  float *data;

  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = ampl*(float)func_alpha((double)i,(double)t0,(double)alpha);
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                                BOXCAR_FARRAY                              */
/*                                                                           */
/*****************************************************************************/
float *boxcar_farray(t0,width,ampl,n)
     float t0,width,ampl;
     int n;
{
  int i;
  float *data;

  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = ampl*(float)func_boxcar((double)i,(double)t0,(double)width);
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                                DIFF_EXP_FARRAY                            */
/*                                                                           */
/*****************************************************************************/
float *diff_exp_farray(t0,a,tau1,tau2,n)
     float t0,a,tau1,tau2;
     int n;
{
  int i;
  float *data;

  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = a * (float)func_diff_exp((double)i,(double)t0,(double)tau1,
				       (double)tau2);
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                               GAMMA_N3_FARRAY                             */
/*                                                                           */
/*****************************************************************************/
float *gamma_n3_farray(m,s,ampl,n)
     float m,s,ampl;
     int n;
{
  int i;
  float *data;

  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = ampl*(float)func_gamma_n3((double)i-(double)m,(double)s);
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                            DIFF_INV_GAMMA_FARRAY                          */
/*                                                                           */
/*                              by Blake Richards                            */
/*                                                                           */
/*  A difference of inverse gamma distribution functions.                    */
/*                                                                           */
/*****************************************************************************/
float *diff_inv_gamma_farray(scale,alpha,beta,posamp,negamp,n)
     float scale,alpha,beta,posamp,negamp;
     int n;
{
  int i, peak;
  float x;
  float *data;
  
  data = (float *)myalloc(n*sizeof(float));
  peak = -1;
  
  // create the first gamma function, then subtract the second
  for(i=0; i<n; i++){
    x = scale*i;
    data[i] = (float)posamp*func_inv_gamma((double)x,(double)alpha,
					   (double)beta);
    
    // check for the peak
    if (i > 0 && peak == -1){
      if (data[i] < data[i-1]){
	peak = i-1;
      }
    }
    
    // if passed the peak, subtract the second function
    if(peak != -1 && i > peak){
      x = scale*(i - peak);
      data[i] -= (float)negamp*func_inv_gamma((double)x,(double)alpha,
					      (double)beta); 
    }
  }
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                               COSINE_FARRAY                               */
/*                                                                           */
/*****************************************************************************/
float *cosine_farray(t0,freq,phase,ampl,n)
     float t0,freq,phase,ampl;
     int n;
{
  int i;
  float *data;

  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = ampl*cos(freq*((double)i-t0) + phase);
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                               GAUSSIAN_FARRAY                             */
/*                                                                           */
/*****************************************************************************/
float *gaussian_farray(mu,sigma,ampl,n)
     float mu,sigma,ampl;
     int n;
{
  int i;
  float *data;

  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = ampl*(float)func_gaussian((double)i,(double)mu,(double)sigma);
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GAUSSIAN_ONE_FARRAY                            */
/*                                                                           */
/*****************************************************************************/
float *gaussian_one_farray(mu,sigma,ampl,n)
     float mu,sigma,ampl;
     int n;
{
  int i;
  float *data;

  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = ampl*(float)func_gaussian_one((double)i,(double)mu,
					    (double)sigma);
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                             FUNC_LOGISTIC_FARRAY                          */
/*                                                                           */
/*****************************************************************************/
float *func_logistic_farray(m,beta,ampl,n)
     float m,beta,ampl;
     int n;
{
  int i;
  float *data;

  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = ampl*(float)func_logistic((double)i-(double)m,(double)beta);
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                                MAXWELL_FARRAY                             */
/*                                                                           */
/*****************************************************************************/
float *maxwell_farray(m,s,ampl,n)
     float m;         // Time origin
     float s;         // scale parameter
     float ampl;      // amplitude
     int n;           // length of array to return
{
  int i;
  float *data;

  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = ampl*(float)func_maxwell((double)i-(double)m,(double)s);
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                            FUNC_PHOTORECPTOR_FARRAY                       */
/*                                                                           */
/*****************************************************************************/
float *func_photoreceptor_farray(m,s,ampl,n)
     float m,s,ampl;
     int n;
{
  int i;
  float *data;

  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = ampl*(float)func_photoreceptor((double)i-(double)m,(double)s);
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                        FUNC_PRIMATE_PHOTOCURRENT_FARRAY                   */
/*                                                                           */
/*****************************************************************************/
float *func_primate_photocurrent_farray(m,ampl,tr,td,tp,phi,n)
     float m,ampl,tr,td,tp,phi;
     int n;
{
  int i;
  float *data;

  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = ampl*(float)func_primate_photocurrent((double)i-(double)m,
						    (double)tr,(double)td,
						    (double)tp,(double)phi);
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                               FUNC_QUICK_FARRAY                           */
/*                                                                           */
/*  Set s=0.01 and n=100 to get points from 0.00 to 0.99.                    */
/*                                                                           */
/*****************************************************************************/
float *func_quick_farray(s,a,b,c,d,n)
     float s,a,b,c,d;
     int n;
{
  int i;
  float *data;

  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = (float)func_quick(a,b,c,d,(double)((float)i*s));
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                             FUNC_LOGNORMAL_FARRAY                         */
/*                                                                           */
/*****************************************************************************/
float *func_lognormal_farray(m,s,ampl,n)
     float m,s,ampl;
     int n;
{
  int i;
  float *data;

  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = ampl*(float)func_lognormal((double)i-(double)m,(double)s);
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                            ONE_SIDED_EXP_FARRAY                           */
/*                                                                           */
/*****************************************************************************/
float *one_sided_exp_farray(t0,tau,ampl,n)
     float t0,tau,ampl;
     int n;
{
  int i;
  float *data;

  //printf("tau = %f\n",tau);
  //printf("ampl = %f\n",ampl);

  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = ampl*(float)func_one_sided_exp((double)i,(double)t0,(double)tau);
  return data;
}
