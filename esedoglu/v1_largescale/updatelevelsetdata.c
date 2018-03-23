#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#include<stdio.h>
#include<math.h>
#include<mex.h>

double surfacetension(double ori1, double ori2);

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[])
{
 /* 
  *  updatelevelsetdata(presence,grains,ID,ORI);
  *
  *  CAUTION: MODIFIES MATLAB INPUT *IN PLACE*.
  */

  mxArray *indices, *grainlevsetvals, *grainconvvals, *locs;
  double *pindices, *plocs, *pgrainlevsetvals, *pgrainconvvals, *id, *ori;
  double sum, mink, st;
  int N,i,j,k,ell,dims,nograins,gind,gind2,idk,idell;
  double temp[100], phi[100], minphi[100];
  
  dims = mxGetM(prhs[0]); /* Number of pixels. */
  N = mxGetM(prhs[1]);    /* Number of grains. */  

  id = (double *) mxGetData(prhs[2]);  /* List of grain IDs. */
  ori = (double *) mxGetData(prhs[3]); /* Grain orientations. */
  
  for (j=0;j<dims;j++){ /* Loop over pixels. */
    indices = mxGetCell(prhs[0],j); /* Grains near this pixel. */
    pindices = (double *) mxGetData(indices);
    locs = mxGetCell(prhs[0],2*dims+j); /* Location of pixel in grain's data. */
    plocs = (double *) mxGetData(locs);
    nograins = mxGetN(indices); /* Number of grains near this pixel. */
    
    for (k=0;k<nograins;k++){ /* Loop over grains. */
      gind = (int) pindices[k]; /* Index of grain in list of all grains. */
      i = (int) plocs[k]-1; /* Location of pixel within grain's data. */
      grainconvvals = mxGetCell(prhs[1],2*N+gind-1);
      pgrainconvvals = (double *) mxGetData(grainconvvals);
      temp[k] = pgrainconvvals[i];
    }

    /* These lines implement the redistribution step in Esedoglu-Otto algorithm: */
    /* Form the "phi" functions: */
    for (k=0;k<nograins;k++){
      gind = (int) pindices[k]; /* Index of grain in list of all grains. */
      sum = 0.0;
      idk = (int) id[gind-1]; /* id of the k-th grain in the local list. */
      for (ell=0;ell<nograins;ell++){
        if (ell != k) {
          gind = (int) pindices[ell]; /* Index of grain in list of all grains. */
          idell = (int) id[gind-1]; /* id of the ell-th grain in the local list. */
          st = surfacetension(ori[idk-1],ori[idell-1]);
          sum = sum + st*temp[ell];
        }
      }
      phi[k] = sum;
    }
    
    /* Minimization over the "phi" functions involved in forming level set functions: */
    for (k=0;k<nograins;k++){
      mink = 1e100;
      for (ell=0;ell<nograins;ell++){
        if (ell != k) {
          mink = min( mink , phi[ell] );
        }
      }
      minphi[k] = mink;
    }    

    /* Form the level set functions: */
    for (k=0;k<nograins;k++){
      gind = (int) pindices[k]; /* Index of grain in list of all grains. */
      i = (int) plocs[k]-1; /* Location of pixel within grain's data. */
      grainlevsetvals = mxGetCell(prhs[1],N+gind-1);
      pgrainlevsetvals = (double *) mxGetData(grainlevsetvals);

      pgrainlevsetvals[i] = (minphi[k] - phi[k]);
      if (nograins==1) pgrainlevsetvals[i] = temp[k];
    }
    
  } /* (for j). Loop over pixels ends. */
    
}

double surfacetension(double ori1, double ori2)
{
  double ang1, minang, st;
  if (ori2>ori1) ang1 = 6.28 - ori2 + ori1;
  else ang1 = 6.28 - ori1 + ori2;
  minang = min( ang1, fabs(ori1-ori2) );
  
  /* Read-Shockley with Brandon angle */
  /*
  if (minang>1.5) st = 1;
  else st = minang / 1.5 * ( 1-log(minang/1.5) );
  */

  /* Equal surface tensions: */  
  st = 1.0;
  
  return st;
}
