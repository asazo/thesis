#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#include<stdio.h>
#include<math.h>
#include<mex.h>

double H(double x);
double IH(double x);
double I2H(double x);
double I3H(double x);
int wrap(int j, int n);

/*************************************************************************/
void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[])
{

  /* **********************************************************************
   * MATLAB Interface:
   * **********************************************************************
   * vf = ls2vf3D(x,y,z,ls,Z,n1,n2,n3)
   * **********************************************************************
   * This routine returns a volume-of-fluid style representation for an
   * interface from an input of a level-set style representation.
   * **********************************************************************
   * INPUT:
   *  x = x-coords of pixels at which level-set representation is given.
   *      One dimensional array. SHOULD BE int32!!!
   *  y = y-coords of pixels at which level-set representation is given.
   *      One dimensional array. SHOULD BE int32!!!
   *  z = z-coords of pixels at which level-set representation is given.
   *      One dimensional array. SHOULD BE int32!!!
   *  ls = Level set values at the x, y & z coordinates. 1D array.
   *  Z = Workspace variable; should be all -1's. Size should be the same
   *      as the 3D computational grid.
   *  n1 = First dimension of computational grid.
   *  n2 = Second dimension of computational grid.
   *  n3 = Third dimension of computational grid.
   * **********************************************************************
   */

  int *x, *y, *z;
  int i, j, k, p, ind, xk, yk, zk, E, W, N, S, U, D, Npix, m, n1, n2, n3, test;
  double *ls, *Z, *vf;
  double mx, my, mz, val;
  
  x = (int *) mxGetData(prhs[0]); /* Array of x coords. */
  y = (int *) mxGetData(prhs[1]); /* Array of y coords. */
  z = (int *) mxGetData(prhs[2]); /* Array of z coords. */
  ls = (double *) mxGetData(prhs[3]); /* Level set values. */
  Z = (double *) mxGetData(prhs[4]); /* Workspace. */

  Npix = mxGetM(prhs[0]);  /* Number of pixels currently in the subset. */
  n1 = (int) mxGetScalar(prhs[5]); /* Dimension of grid. */
  n2 = (int) mxGetScalar(prhs[6]); /* Dimension of grid. */
  n3 = (int) mxGetScalar(prhs[7]); /* Dimension of grid. */

  //printf("n1 n2 n3 = %i %i %i\n",n1,n2,n3);
  //printf("Npix = %i\n",Npix);

  /* Allocate memory and get pointer for output variables. */
  plhs[0] = mxCreateDoubleMatrix(Npix,1,mxREAL);	
  vf = mxGetPr(plhs[0]);
  
  /* Populate Z with level set values: */
  for (k=0;k<Npix;k++) { /* Loop over pixels. */
    xk = x[k] - 1;
    yk = y[k] - 1;
    zk = z[k] - 1;
    ind = xk + n1*yk + n1*n2*zk;
    Z[ind] = ls[k];
    /* printf("%f\n",ls[k]); */
  }
  
  /* Real action: */
  for (k=0;k<Npix;k++){ /* Loop over pixels. */
    xk = x[k] - 1;
    yk = y[k] - 1;
    zk = z[k] - 1;
    val = ls[k];
    ind = xk + n1*yk + n1*n2*zk;

    E = wrap(xk+1,n1) + n1*yk + n1*n2*zk;
    W = wrap(xk-1,n1) + n1*yk + n1*n2*zk;
    N = xk + n1*wrap(yk+1,n2) + n1*n2*zk;
    S = xk + n1*wrap(yk-1,n2) + n1*n2*zk;
    U = xk + n1*yk + n1*n2*wrap(zk+1,n3);
    D = xk + n1*yk + n1*n2*wrap(zk-1,n3);
        
    mx = 0.5*(Z[E]-Z[W]);
    my = 0.5*(Z[N]-Z[S]);
    mz = 0.5*(Z[U]-Z[D]);

    test = 0;
    if (fabs(mx)<1e-6) test = test + 1;
    if (fabs(my)<1e-6) test = test + 2;
    if (fabs(mz)<1e-6) test = test + 4;    

    switch (test) {

      case 0 :
        //printf("case 0\n");
        vf[k] = I3H(mx+my+mz+val) - I3H(mx+my-mz+val) -
                I3H(mx-my+mz+val) + I3H(mx-my-mz+val) -
                I3H(-mx+my+mz+val) + I3H(-mx+my-mz+val) +
                I3H(-mx-my+mz+val) - I3H(-mx-my-mz+val);
        vf[k] = vf[k] / (mx*my*mz);
        break;

      case 1 :
        //printf("case 1\n");
        vf[k] = I2H(my+mz+val) - I2H(my-mz+val) +
                I2H(-my-mz+val) - I2H(-my+mz+val);
        vf[k] = 2*vf[k]/(my*mz);
        break;

      case 2 :
        //printf("case 2\n"); 
        vf[k] = I2H(mx+mz+val) - I2H(mx-mz+val) +
                I2H(-mx-mz+val) - I2H(-mx+mz+val);
        vf[k] = 2*vf[k]/(mx*mz);
        break;

      case 4 :
        //printf("case 4\n");
        vf[k] = I2H(mx+my+val) - I2H(mx-my+val) -
                I2H(-mx+my+val) + I2H(-mx-my+val);
        vf[k] = 2*vf[k]/(mx*my);
        break;

      case 3 :
        //printf("case 3\n");
        vf[k] = IH(mz+val) - IH(-mz+val);
        vf[k] = 4*vf[k]/mz;
        break;

      case 5 :
        //printf("case 5\n");
        vf[k] = IH(my+val) - IH(-my+val);
        vf[k] = 4*vf[k]/my;
        break;

      case 6 :
        //printf("case 6\n");
        vf[k] = IH(mx+val) - IH(-mx+val);
        vf[k] = 4*vf[k]/mx;
        break;

      case 7 :
        //printf("case 7\n");
        vf[k] = 8*H(val);
        break;

      } /* end switch */
    
    vf[k] = vf[k] / 8.0;
    
    if ( max(max(max(max(max(max(Z[E],Z[W]),Z[N]),Z[S]),val),Z[U]),Z[D]) < 0 ) { 
      vf[k] = 0;
    }
    if ( min(min(min(min(min(min(Z[E],Z[W]),Z[N]),Z[S]),val),Z[U]),Z[D]) > 0 ) {
      vf[k] = 1.0;
    }
    
  } /* for k */

  /* Clean up Z: */
  for (k=0;k<Npix;k++){
    xk = x[k]-1;
    yk = y[k]-1;
    zk = z[k]-1;
    ind = xk + n1*yk + n1*n2*zk;
    Z[ind] = -1;
  }

}

double H(double x) 
{
  if (x <= 0) return 0.0;
  else return 1.0;
}
  
double IH(double x) 
{
  if (x<=0) return 0.0;
  else return x;
}
  
double I2H(double x) 
{
  if (x<=0) return 0.0;
  else return 0.5*x*x;
}

double I3H(double x)
{
  if (x<=0) return 0.0;
  else return x*x*x/6.0;
}

int wrap(int j, int n)
{
  if (j<0) return (n-1);
  else {
    if (j>n-1) return 0;
  }
  return j;
}
