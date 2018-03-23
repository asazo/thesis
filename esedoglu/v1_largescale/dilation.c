#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#include<stdio.h>
#include<math.h>
#include<mex.h>

struct grid3dsubset {
  int *x; /* List of element indices. */
  int *y; /* Ditto. */
  int *z; /* Ditto. */
  int *bx;/* List of bdry element indices. */
  int *by;/* Ditto. */
  int *bz;/* Ditto. */
  int *p; /* Stores position of pixels in element list. */
  int n1; /* Innermost dimension of grid. */
  int n2; /* Middle dimension of grid. */
  int n3; /* Outer dimension of grid. */
  int N; /* Number of elements. */
  int bN; /* Number of boundary elements. */
  int maxN; /* Max number of elements allowed. */
};

struct intarray1d {
  int *D; /* Points to beginning of data. */
  int N;  /* Number of elements. */
};

void create_grid3dsubset(struct grid3dsubset *in, int n1, int n2, int n3, int maxN, int *W);
void destroy_grid3dsubset(struct grid3dsubset *in);
void add_element_grid3dsubset(struct grid3dsubset *in, int x, int y, int z);
int isbdry(struct grid3dsubset *in, int x, int y, int z);
void updatebdry_grid3dsubset(struct grid3dsubset *in);
void grow_grid3dsubset(struct grid3dsubset *in, int w);
void create_intarray1d(struct intarray1d *in, int N);
void destroy_intarray1d(struct intarray1d *u);
void create_grid3dsubset_nocalloc(struct grid3dsubset *in, int n1, int n2, int n3, 
				  int maxN, int *W, int *longarray);
void destroy_grid3dsubset_nocalloc(struct grid3dsubset *in);

/*************************************************************************/
void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[])
{

  /* **********************************************************************
   * MATLAB Interface:
   * **********************************************************************
   * [xo,yo,zo] = dilation(x,y,z,w,W);
   * **********************************************************************
   * Grows (dilates) the subset of grid by w pixels outwards.
   * Periodic boundary conditions used. 4 point neighborhood (E, W, N, S).
   * The vectors x, y, and z should be COLUMN VECTORS of type int32!
   * W is a workspace array of (-1)'s the same size as computational grid 
   * and of type int32.
   * **********************************************************************
   */
  
  int *x, *y, *z, w, *W, *longarray;
  const int *dims;
  double *xout, *yout, *zout;
  int i, j, k, ind, numpix, n1, n2, n3;
  struct grid3dsubset S;

  /* Input variables: */
  x = (int *) mxGetData(prhs[0]); /* Array of x coords. */
  y = (int *) mxGetData(prhs[1]); /* Array of y coords. */
  z = (int *) mxGetData(prhs[2]); /* Array of z coords. */
  w = (int) mxGetScalar(prhs[3]); /* Width of growth. */
  W = (int *) mxGetData(prhs[4]); /* Workspace. */
  longarray = (int *) mxGetData(prhs[5]); /* Workspace. */

  numpix = mxGetM(prhs[0]);  /* Number of pixels currently in the subset. */
  dims = mxGetDimensions(prhs[4]);
  n1 = (int) dims[0]; /* Dimension of grid. */
  n2 = (int) dims[1]; /* Dimension of grid. */
  n3 = (int) dims[2]; /* Dimension of grid. */

  /* Maximum anticipated number of points in the subset is set here: */
  create_grid3dsubset_nocalloc(&S, n1, n2, n3, 20000000, W, longarray);
  
  for (j=0;j<numpix;j++){ /* Fill the subset with points from list of coords. */
    add_element_grid3dsubset(&S,x[j]-1,y[j]-1,z[j]-1);
  }

  updatebdry_grid3dsubset(&S); /* Boundary data structure made consistent. */
  grow_grid3dsubset(&S,w); /* Real work happens here. */

  /* Allocate memory and get pointer for output variables. */
  plhs[0] = mxCreateDoubleMatrix(S.N,1,mxREAL);	
  xout = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(S.N,1,mxREAL);
  yout = mxGetPr(plhs[1]);
  plhs[2] = mxCreateDoubleMatrix(S.N,1,mxREAL);
  zout = mxGetPr(plhs[2]);
  
  for (i=0;i<S.N;i++){ /* Fill output variables with coordinates. */
    *(xout+i) = S.x[i]+1;
    *(yout+i) = S.y[i]+1;
    *(zout+i) = S.z[i]+1;
  }

  destroy_grid3dsubset_nocalloc(&S);
  
}
/*************************************************************************/
void create_grid3dsubset(struct grid3dsubset *in, int n1, int n2, int n3, int maxN, int *W)
{
  
  /* Allocate memory. Set to empty subset. */
  (*in).x = (int *) calloc(maxN,sizeof(int));
  (*in).y = (int *) calloc(maxN,sizeof(int));
  (*in).z = (int *) calloc(maxN,sizeof(int));
  (*in).bx = (int *) calloc(maxN,sizeof(int));
  (*in).by = (int *) calloc(maxN,sizeof(int));
  (*in).bz = (int *) calloc(maxN,sizeof(int));
  (*in).n1 = n1;
  (*in).n2 = n2;
  (*in).n3 = n3;
  (*in).p = W; /* Workspace array that is provided by caller is used. */
  (*in).N = 0;
  (*in).bN = 0;
  (*in).maxN = maxN;
      
}

/*************************************************************************/
void create_grid3dsubset_nocalloc(struct grid3dsubset *in, int n1, int n2, int n3, 
				  int maxN, int *W, int *longarray)
{

  /* Allocate memory. Set to empty subset. */
  (*in).x = longarray;
  (*in).y = longarray + maxN;
  (*in).z = longarray + 2*maxN;
  (*in).bx = longarray + 3*maxN;
  (*in).by = longarray + 4*maxN;
  (*in).bz = longarray + 5*maxN;
  (*in).n1 = n1;
  (*in).n2 = n2;
  (*in).n3 = n3;
  (*in).p = W; /* Workspace array that is provided by caller is used. */
  (*in).N = 0;
  (*in).bN = 0;
  (*in).maxN = maxN;

}
/*************************************************************************/
void destroy_grid3dsubset(struct grid3dsubset *in)
{
  int i, n2, n3, ind, x, y, z;
      
  /* Clear the workspace array: */
  n2 = (*in).n2;
  n3 = (*in).n3;
  for (i=0;i<(*in).N;i++) {
    x = (*in).x[i];
    y = (*in).y[i];
    z = (*in).z[i];
    ind = x*n2*n3 + y*n3 + z;
    (*in).p[ind] = -1;
  }

  free( (*in).x );
  free( (*in).y );
  free( (*in).z );
  free( (*in).bx );
  free( (*in).by );
  free( (*in).bz );
  
}
/*************************************************************************/
void destroy_grid3dsubset_nocalloc(struct grid3dsubset *in)
{
  int i, n2, n3, ind, x, y, z;
      
  /* Clear the workspace array: */
  n2 = (*in).n2;
  n3 = (*in).n3;
  for (i=0;i<(*in).N;i++) {
    x = (*in).x[i];
    y = (*in).y[i];
    z = (*in).z[i];
    ind = x*n2*n3 + y*n3 + z;
    (*in).p[ind] = -1;
  }
  
}
/*************************************************************************/
void add_element_grid3dsubset(struct grid3dsubset *in, int x, int y, int z)
{
  /* Adds element to grid subset. Note that boundary is NOT updated. */
  int n1,n2,n3,N,ind;  
  N = (*in).N;
  
  n1 = (*in).n1;
  n2 = (*in).n2;
  n3 = (*in).n3;

  if (N == 0) {/* If subset currently empty... */
    (*in).x[0] = x;
    (*in).y[0] = y;
    (*in).z[0] = z;
    ind = x*n2*n3 + y*n3 + z;
    (*in).p[ind] = 0;
    ((*in).N)++;
  }
  else { /* If subset already has elements... */
    ind = x*n2*n3 + y*n3 + z;
    if ( (*in).p[ind] == -1 ) { /* If new loc. not already in subset... */
      (*in).p[ind] = N; /* Stores position of new element in the list. */
      (*in).x[N] = x; /* Add new x coord to list of x coords. */
      (*in).y[N] = y; /* Add new y coord to list of y coords. */
      (*in).z[N] = z; /* Add new y coord to list of z coords. */    
      ((*in).N)++;
    } /* end if */
  } /* end else. */
}
/*************************************************************************/
int isbdry(struct grid3dsubset *in, int x, int y, int z)
{
  int neighmin, flag, N, S, E, W, U, D, n1, n2, n3, ind;

  n1 = (*in).n1;
  n2 = (*in).n2;
  n3 = (*in).n3;

  E = x+1;
  if (E == n1) E = 0;
  W = x-1;
  if (W == -1) W = n1-1;
  N = y+1;
  if (N == n2) N = 0;
  S = y-1;
  if (S == -1) S = n2-1;
  U = z+1;
  if (U == n3) U = 0;
  D = z-1;
  if (D == -1) D = n3-1;

  neighmin = 0;
  ind = E*n2*n3 + y*n3 + z;
  neighmin = min( neighmin, (*in).p[ind] );
  ind = W*n2*n3 + y*n3 + z;
  neighmin = min( neighmin, (*in).p[ind] );
  ind = x*n2*n3 + N*n3 + z;
  neighmin = min( neighmin, (*in).p[ind] );
  ind = x*n2*n3 + S*n3 + z;
  neighmin = min( neighmin, (*in).p[ind] );
  ind = x*n2*n3 + y*n3 + U;
  neighmin = min( neighmin, (*in).p[ind] );
  ind = x*n2*n3 + y*n3 + D;  
  neighmin = min( neighmin, (*in).p[ind] );
    
  if (neighmin < 0) flag = 1;
  else flag = 0;

  return flag;
}
/*************************************************************************/
void updatebdry_grid3dsubset(struct grid3dsubset *in)
{
  int i,x,y,z,N,bN;
  
  N = (*in).N; /* Number of elements in grid subset. */
  bN = 0;  

  for (i=0;i<N;i++){
    x = (*in).x[i];
    y = (*in).y[i];
    z = (*in).z[i];
    if ( isbdry(in,x,y,z) == 1 ) {
      (*in).bx[bN] = x;
      (*in).by[bN] = y;
      (*in).bz[bN] = z;
      bN++;
    }
  }
  
  (*in).bN = bN;
}
/*************************************************************************/
void grow_grid3dsubset(struct grid3dsubset *in, int w)
{
  int i,j,k,n1,n2,n3,t,x,y,z,bx,by,bz,pos,count,bN,ind,E,W,N,S,U,D,maxN;
  struct intarray1d candidate_x, candidate_y, candidate_z;

  n1 = (*in).n1;
  n2 = (*in).n2;
  n3 = (*in).n3;
  maxN = (*in).maxN;

  create_intarray1d(&candidate_x,maxN);
  create_intarray1d(&candidate_y,maxN);
  create_intarray1d(&candidate_z,maxN);

  for (t=0;t<w;t++){

    count = 0; /* Counts how many new points added to the subset
		              on current pass. */
   
    for (i=0;i<(*in).bN;i++){ /* Go through current boundary points. */
      bx = (*in).bx[i];
      by = (*in).by[i];
      bz = (*in).bz[i];

      /* Add neighbors of point (bx,by) if they aren't already in subset.
       * Four-neighborhood is used.
       */

      /* East neighbor: */
      E = bx+1;
      if (E==n1) E=0;
      ind = E*n2*n3 + by*n3 + bz;
      if ( (*in).p[ind] < 0 ) {
        pos = (*in).N;
        (*in).x[pos] = E;
        (*in).y[pos] = by;
	(*in).z[pos] = bz;
        (*in).N++;
        (*in).p[ind] = pos;
        candidate_x.D[count] = E;
        candidate_y.D[count] = by;
	candidate_z.D[count] = bz;
        count++;
      }
      
      /* West neighbor: */
      W = bx-1;
      if (W==-1) W=n1-1;
      ind = W*n2*n3 + by*n3 + bz;
      if ( (*in).p[ind] < 0 ) {
        pos = (*in).N;
        (*in).x[pos] = W;
        (*in).y[pos] = by;
	(*in).z[pos] = bz;
        (*in).N++;
        (*in).p[ind] = pos;
        candidate_x.D[count] = W;
        candidate_y.D[count] = by;
	candidate_z.D[count] = bz;
        count++;
      }
      
      /* North neighbor: */
      N = by+1;
      if (N==n2) N=0;
      ind = bx*n2*n3 + N*n3 + bz;
      if ( (*in).p[ind] < 0 ) {
        pos = (*in).N;
        (*in).x[pos] = bx;
        (*in).y[pos] = N;
	(*in).z[pos] = bz;
        (*in).N++;
        (*in).p[ind] = pos;
        candidate_x.D[count] = bx;
        candidate_y.D[count] = N;
	candidate_z.D[count] = bz;
        count++;
      }
      
      /* South neighbor: */
      S = by-1;
      if (S==-1) S=n2-1;
      ind = bx*n2*n3 + S*n3 + bz;
      if ( (*in).p[ind] < 0 ) {
        pos = (*in).N;
        (*in).x[pos] = bx;
        (*in).y[pos] = S;
	(*in).z[pos] = bz;
        (*in).N++;
        (*in).p[ind] = pos;
        candidate_x.D[count] = bx;
        candidate_y.D[count] = S;
	candidate_z.D[count] = bz;
        count++;
      }

      /* Up neighbor: */
      U = bz+1;
      if (U==n3) U=0;
      ind = bx*n2*n3 + by*n3 + U;
      if ( (*in).p[ind] < 0 ) {
        pos = (*in).N;
        (*in).x[pos] = bx;
        (*in).y[pos] = by;
	(*in).z[pos] = U;
        (*in).N++;
        (*in).p[ind] = pos;
        candidate_x.D[count] = bx;
        candidate_y.D[count] = by;
	candidate_z.D[count] = U;
        count++;
      }
      
      /* Down neighbor: */
      D = bz-1;
      if (D==-1) D=n3-1;
      ind = bx*n2*n3 + by*n3 + D;
      if ( (*in).p[ind] < 0 ) {
        pos = (*in).N;
        (*in).x[pos] = bx;
        (*in).y[pos] = by;
	(*in).z[pos] = D;
        (*in).N++;
        (*in).p[ind] = pos;
        candidate_x.D[count] = bx;
        candidate_y.D[count] = by;
	candidate_z.D[count] = D;
        count++;
      }
      
    } /* i loop (through current bdry points) ends. */
    
    /* Check the newly added points; new bdry is a subset of those. */
    bN = 0;
    for (i=0;i<count;i++){
      x = candidate_x.D[i];
      y = candidate_y.D[i];
      z = candidate_z.D[i];
      if ( isbdry(in,x,y,z) == 1 ) {
        (*in).bx[bN] = x;
        (*in).by[bN] = y;
	(*in).bz[bN] = z;
        bN++;
      }
    }
    (*in).bN = bN;
        
    
  } /* t loop ends. */
  destroy_intarray1d(&candidate_x);
  destroy_intarray1d(&candidate_y);
  destroy_intarray1d(&candidate_z);
}
/*************************************************************************/
void create_intarray1d(struct intarray1d *in, int N)
{
  (*in).D = (int *) calloc(N,sizeof(int));
  (*in).N = N;
}
/*************************************************************************/
void destroy_intarray1d(struct intarray1d *u)
{
  free( (*u).D );
}
/*************************************************************************/
