/*****************************************************************************/
/*                                                                           */
/*  min_util.c                                                               */
/*                                                                           */
/*  Minimization                                                             */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "my_util.h"
#include "iarray_util.h"
#include "farray_util.h"

#define GLIMIT 100.0                 // Max mag. for one parabolic fit step
#define TINY 1.0e-20                 // Small number

#define GOLDRAT 1.61803399           // Golden ratio, 1.6180339887...
#define GRATIO 0.381966              // Golden ratio
#define MAX_ITER  99                 // Maximum iterations
#define MAX_ITER_PR  222             // Maximum iterations
#define EPS_0     1.0e-10            // Stay away from zero minimum
#define SHFT(a,b,c,d) (a)=(b); (b)=(c); (c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))


#define TOLERANCE  2.0e-4

// GLOBALS
int     min_util_glob_n;             // Number of dimensions
float  *min_util_glob_x;             // Current coordinates
float  *min_util_glob_v;             // Current vector
float (*min_util_glob_f)(float []);  // Funciton being minimized

/**************************************-**************************************/
/*                                                                           */
/*                               MINU_TEST_FUNC                              */
/*                                                                           */
/*  Simple functions to test the min computations.                           */
/*                                                                           */
/*****************************************************************************/
float minu_test_func(x)
  float *x;
{
  float r;
  float x0,x1,c;

  if (1){
    //
    //  f(x1,x2) = 1.0 + (x1-c)^2 * x2^2;
    //
    c = 3.0;
    x0 = x[0];
    x1 = x[1];
    r = 1.0 + (x0-c)*(x0-c) + x1*x1;
  }else{
    ;
  }

  //printf("  f(%f,%f) = %f\n",x[0],x[1],r);

  return r;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MINU_TEST_FUNCD                              */
/*                                                                           */
/*  Simple functions to test the min computations.                           */
/*                                                                           */
/*****************************************************************************/
int minu_test_funcd(x,d)
  float *x;
  float *d;  // Returned derivative values
{
  float c;

  if (1){
    //  f(x1,x2) = (x1-c)^2 * x2^2;
    c = 3.0;
    d[0] = 2.0*(x[0]-c);  // d/dx1
    d[1] = 2.0*x[1];  // d/dx2
  }else{
    d[0] = 0.0;
    d[1] = 0.0;
  }
  //printf("d[0] = %f\n",d[0]);

  return 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                                  MINU_FA                                  */
/*                                                                           */
/*  Return the value of the global function evaluated at 'a' times the       */
/*  current vector added to the current point.                               */
/*                                                                           */
/*****************************************************************************/
float minu_fa(float a)  // This format needed to make argument work
{
  int i;
  float f,*xt;

  //printf("______a = %f\n",a);
  //exit(0);

  //
  //  WYETH - THIS PROC SHOULD BE PUT IN-LINE ?
  //

  xt = get_farray(min_util_glob_n);

  for(i=0;i<min_util_glob_n;i++)
    xt[i] = min_util_glob_x[i] + a * min_util_glob_v[i];

  //printf("(minu_fa) a = %f ==>  x0,x1 = %f %f\n",a,xt[0],xt[1]);

  f = min_util_glob_f(xt);

  myfree(xt);

  return f;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MINU_INIT_RANGE                              */
/*                                                                           */
/*  Find an initial range that contains a minimum.                           */
/*                                                                           */
/*****************************************************************************/
void minu_init_range(ax,bx,cx,fa,fb,fc,f)
     float *ax;  // Initial point;  Return an endpoint of the range
     float *bx;  // Initial point;  Return a lower point within the range
     float *cx;  //                 Return an endpoint of the range
     float *fa;
     float *fb;
     float *fc;
     float (*f)(float);
{
  float ulim,u,r,q,fu,dum;

  //printf("ax = %f  bx=%f \n",*ax,*bx);

  *fa = (*f)(*ax);
  *fb = (*f)(*bx);
  if (*fb > *fa){
    //temp = fa;    fa = fb;    fb = temp;  // Swap the y-values
    //temp = ax;    ax = bx;    bx = temp;  // Swap the x-values
    SHFT(dum,*ax,*bx,dum);
    SHFT(dum,*fb,*fa,dum);
  }

  *cx = (*bx) + GOLDRAT*(*bx - *ax);       // First guess for c
  *fc = (*f)(*cx);
  while(*fb > *fc){
    r = (*bx - *ax) * (*fb - *fc);
    q = (*bx - *cx) * (*fb - *fa);
    // WYETH - Replace FMAX with something simpler
    u = (*bx) - ((*bx - *cx)*q - (*bx - *ax)*r) /
      (2.0 * SIGN(FMAX(fabs(q-r),TINY),q-r));
    ulim = (*bx) + GLIMIT * (*cx - *bx);

    if ((*bx - u) * (u - *cx) > 0.0){  // Go between 'b' and 'c'
      fu = (*f)(u);
      if (fu < *fc){        // Found min between 'b' and 'c'
	*ax = (*bx);
	*bx = u;
	*fa = (*fb);
	*fb = fu;
	return;
      }else if (fu > *fb){  // Found min between 'a' and 'u'
	*cx = u;
	*fc = fu;
	return;
      }
      u = (*cx) + GOLDRAT * (*cx - *bx);  // Parabolic fit no good
      fu = (*f)(u);
    }else if ((*cx-u) * (u - ulim) > 0.0){
      fu = (*f)(u);
      if (fu < *fc){
	SHFT(*bx,*cx,u,*cx + GOLDRAT*(*cx - *bx));
	SHFT(*fb,*fc,fu,(*f)(u));
      }
    }else if (((u - ulim) * (ulim - *cx)) >= 0.0){
      u = ulim;
      fu = (*f)(u);
    }else{
      u = (*cx) + GOLDRAT * (*cx - *bx);
      fu = (*f)(u);
    }
    SHFT(*ax,*bx,*cx,u);
    SHFT(*fa,*fb,*fc,fu);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               MINU_BRENT_1D                               */
/*                                                                           */
/*  Return the value of 'f' at its minimum, using Richard Brent's method.    */
/*                                                                           */
/*****************************************************************************/
float minu_brent_1d(ax,bx,cx,f,tlrnc,xmin)
     float ax;          // The triplet (ax,bx,cx) contains a minimum, because
     float bx;          //   ax < bx < cx, and f(bx) < f(ax) and f(cx)
     float cx;          //
     float (*f)(float); // Minimize this function
     float tlrnc;       // To within this tolerance
     float *xmin;       // x-coord at found minimum
{
  float a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;

  // WYETH
  int i;

  if (ax < cx){
    a = ax;
    b = cx;
  }else{
    a = cx;
    b = ax;
  }

  x = w = v = bx;
  fw = fv = fx = (*f)(x);

  e = 0.0;
  for(i=0;i<MAX_ITER;i++){
    xm = 0.5 * (a+b);
    tol2 = 2.0 * (tol1 = tlrnc*fabs(x)+EPS_0);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))){      // Are we done?
      *xmin = x;
      return fx;
    }
    if (fabs(e) > tol1){
      r = (x-w) * (fx-fv);
      q = (x-v) * (fx-fw);
      p = (x-v)*q - (x-w)*r;
      q = 2.0*(q-r);
      if (q > 0.0)
	p = -p;
      q = fabs(p);
      etemp = e;
      e = d;
      if (fabs(p) >= fabs(0.5*q*etemp) || (p <= q*(a-x)) || (p >= q*b-x)){
	d = GRATIO * (e = (x >= xm ? a-x : b-x));
      }else{
	d = p/q;
	u = x+d;
	if ((u-a < tol2) || ((b-u) < tol2))
	  d = SIGN(tol1,xm-x);
      }
    }else{
      d = GRATIO * (e = (x >= xm ? a-x : b-x));
    }
    u = (fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));

    fu = (*f)(u);  // Evaluate the function

    if (fu <= fx){
      if (u >= x)
	a = x;
      else
	b = x;
      SHFT(v,w,x,u);
      SHFT(fv,fw,fx,fu);
    }else{
      if (u < x)
	a = u;
      else
	b = u;
      if ((fu <= fw) || (w == x)){
	v = w;
	w = u;
	fv = fw;
	fw = fu;
      }else if ((fu <= fv) || (v == x) || (v == w)){
	v = u;
	fv = fu;
      }
    }
  }

  printf("  *** (MINU_BRENT_1D)  Maximum iterations exceeded\n");

  *xmin = x;
  return fx;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 MINU_LINE                                 */
/*                                                                           */
/*  Return the value of the function 'func' at the minimum that is found     */
/*  by searching along the vector direction 'vd' from the point 'x'.         */
/*  'x' will be replaced by the coordinates of the minimum, and 'vd' will    */
/*  be replaced by the displacement from the initial to final 'x'.           */
/*                                                                           */
/*****************************************************************************/
float minu_line(x,v,n,func)
     float *x;                // [n] starting point
     float *v;                // [n] vector direction
     int n;                   // Number of dimensions
     float (*func)(float *);  // Function to be minimized
{
  float xx,xmin,fx,fb,fa,bx,ax;
  float minu_fa(float x);

  // WYETH 
  int i;
  float min_val;

  //
  //  Store values into globals
  //
  min_util_glob_n = n;
  min_util_glob_x = (float *)myalloc(n*sizeof(float));
  min_util_glob_v = (float *)myalloc(n*sizeof(float)); 
  min_util_glob_f = func;
  for(i=0;i<n;i++){
    min_util_glob_x[i] = x[i];
    min_util_glob_v[i] = v[i];
  }

  ax = 0.0;
  xx = 1.0;
  minu_init_range(&ax,&xx,&bx,&fa,&fx,&fb,minu_fa);
  //printf("INIT_RANGE: (%f, %f, %f)  f(x):  %f, %f, %f\n",ax,xx,bx,fa,fx,fb);
  min_val = minu_brent_1d(ax,xx,bx,minu_fa,TOLERANCE,&xmin);

  //
  //  Set final values
  //
  for(i=0;i<n;i++){
    v[i] *= xmin;
    x[i] += v[i];
  }

  myfree(min_util_glob_x);
  myfree(min_util_glob_v);

  return min_val;
}
/**************************************-**************************************/
/*                                                                           */
/*                                  MINU_PR                                  */
/*                                                                           */
/*  An implementation of the Polak-Ribiere method.                           */
/*                                                                           */
/*****************************************************************************/
void minu_pr(x,n,tlrnc,func,dfunc,riter,rval)
     float *x;                // [n] starting point
     int n;                   // Number of dimensions
     float tlrnc;             // Tolerance
     float (*func)(float *);  // Function to be minimized
     int (*dfunc)(float *, float *);  // Deriv of func
     //void (*dfunc)(float *, float *);  // Deriv of func
     int   *riter;            // Return number of iterations
                              //   -1 indicates "Too many iterations"
                              //   -2 indicates NaN was encountered
     float *rval;             // Return value at min
{
  int i,j;
  int dflag;
  float gg,gam,fp,dgg;
  float *g,*h,*v;

  // WYETH - Why not store global min_util_glob_n  here?

  g = (float *)myalloc(n*sizeof(float));
  h = (float *)myalloc(n*sizeof(float));
  v = (float *)myalloc(n*sizeof(float));

  //for(i=0;i<n;i++)
  //printf("x[%d] = %f\n",i,x[i]);

  fp = (*func)(x);
  dflag = (*dfunc)(x,v);    // Get the vector 'v' of derivatives
  for(j=0;j<n;j++){
    g[j] = -v[j];
    v[j] = h[j] = g[j];
    //v[j] *= -1.0;  // WYETH - USE THIS
    //h[j] = g[j] = v[j];
  }
  for(i=0;i<MAX_ITER_PR;i++){
    *riter = i;
    *rval = minu_line(x,v,n,func);
    if (isnan(x[0])){
      //printf("_______________IS NAN - EXITING\n");
      *riter = -2;  // NaN was encountered
      myfree(g);
      myfree(h);
      myfree(v);
      return;
    }
    if (2.0*fabs(*rval-fp) <= tlrnc*(fabs(*rval)+fabs(fp)+EPS_0)){
      // Most typical return
      myfree(g);
      myfree(h);
      myfree(v);
      //printf("   fp   = %f\n",fp);
      //printf("   rval = %f\n",*rval);
      //printf("(minu_pr) done 01\n");
      return;
    }
    fp = (*func)(x);
    dflag = (*dfunc)(x,v);
    if (dflag < 0){  // The derivative routine reports an error
      // Note, dflag value should be 1 (success) or >= -10 (error flag)
      printf("(min_util.c)  dflag = %d\n",dflag);
      *riter = dflag;
      myfree(g);
      myfree(h);
      myfree(v);
      return;
    }
    dgg = gg = 0.0;
    for(j=0;j<n;j++){
      gg += g[j]*g[j];
      // dgg += v[j]*v[j];   // Use for Fletcher Reeves method
      dgg += (v[j] + g[j]) * v[j];
    }
    if (gg == 0.0){
      // If the gradient is zero, stop
      myfree(g);
      myfree(h);
      myfree(v);
      printf("(minu_pr) gradient is zero\n");
      return;
    }
    gam = dgg/gg;
    for(j=0;j<n;j++){
      g[j] = -v[j];
      v[j] = h[j] = g[j] = gam*h[j];
    }
  }
  *riter = -1;
  printf("  *** MINU_PR  Too many iterations\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                               MINU_TEST_FIT                               */
/*                                                                           */
/*****************************************************************************/
void minu_test_fit()
{
  int n,niter;
  float minval,tlrnc,*x;
  float (*func)(float *);
  int (*dfunc)(float *, float *);

  tlrnc = 0.00001;

  n = 2;  // Number of dimensions
  x = (float *)myalloc(n*sizeof(float));
  x[0] = 10.0;
  x[1] =  7.0;

  func  = minu_test_func;
  dfunc = minu_test_funcd;
  minu_pr(x,n,tlrnc,func,dfunc,&niter,&minval);

  printf("    Number of iterations:  %d\n",niter);
  printf("    Minimum value:  %f\n",minval);
  printf("    x0,x1 = %f %f\n",x[0],x[1]);

  printf("\n");
  exit(0);
}
