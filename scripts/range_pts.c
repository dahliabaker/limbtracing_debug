///////////////////////////////////////////////////////////////////////////
// Flash Lidar Measurement
// Author: Jay and Steve
// Modified by: Ann 
// 03/27/2016
//          EDIT: now calculated in sensor frame.
///////////////////////////////////////////////////////////////////////////

#include    "mex.h"
#include	<math.h>
#include	<stdarg.h>
#include	<stdio.h>
#include	<stdlib.h>
#include    <string.h>

///////////////////
// -- Outputs -- //
///////////////////

double *range_Ptr;

//////////////////
// -- Inputs -- //
//////////////////

double *pts, *grid_uvec;
double *verts, *nhats, *facets;
int    numPts, n_facets, n_verts;

/////////////////
// -- Locals -- //
/////////////////

double r, p[3], grid_loc[3];
int    hit, iv, ii, jj, kk, mm, Aind, Bind, Cind;
double A_loc[3], B_loc[3], C_loc[3], nhat_loc[3], lpt[3], lvec[3];

/////////////////////
// -- Functions -- //
/////////////////////

//

void inverse ( double answer[], double x[], double y[], double n[] );
int in_triangle_3D( double *A, double *B, double *C, double *P, double *nhat);
double dot ( double a[3], double b[3]);
double line_intersect_plane ( double *l, double *nhat, double *l0, double *p0, double *p);

/***************************************************************
 * main program
 ***************************************************************/
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
        
{
    
    ///////////////////////////////////////////////////////////////////////
    
    /////////////////////////
    // -- Assign Inputs -- //
    /////////////////////////
    
    // - Scalar
    
    numPts               = (int)mxGetScalar(prhs[0]);   // 
    n_facets             = (int)mxGetScalar(prhs[1]);   //
    n_verts              = (int)mxGetScalar(prhs[2]);   //
    
    // - Vector and Matrix
    
    grid_uvec            = mxGetPr(prhs[3]);
    pts                  = mxGetPr(prhs[4]);       // these are sensor grid
    facets               = mxGetPr(prhs[5]);        // 
    verts                = mxGetPr(prhs[6]);       // (change to verts..?)
    nhats                = mxGetPr(prhs[7]);       // 
    
    
    
    ///////////////////////////////////////////////////////////////////////
    
    //////////////////////////
    // -- Assign Outputs -- //
    //////////////////////////
    
    plhs[0]  = mxCreateDoubleMatrix(numPts, 1, mxREAL);
    range_Ptr = mxGetPr(plhs[0]);
    
    
    for (ii=0; ii<numPts; ii++) {
            grid_loc[0] = pts[3*ii];
            grid_loc[1] = pts[3*ii + 1];
            grid_loc[2] = pts[3*ii + 2];
            
            // Add pre-processing here that checks to see if it is within a bubble of the center to stop from
            // having to search all the facets if we know it won't intersect

            // See if this point/ray sees the shape model
            for ( kk=0; kk<n_facets; kk++) {
                
                Aind = (int)facets[3*kk];
                Bind = (int)facets[3*kk+1];
                Cind = (int)facets[3*kk+2];                
                
                for (mm=0; mm<3; mm++) {
                    A_loc[mm] = verts[3* Aind + mm]; // this to verts?
                    B_loc[mm] = verts[3* Bind + mm];
                    C_loc[mm] = verts[3* Cind + mm];
                    nhat_loc[mm] = nhats[3*kk + mm];
                }
    
                // find intersection with plane of triangle
                r = line_intersect_plane ( grid_uvec, nhat_loc, grid_loc, A_loc, p);
            
                // check if intersection is inside triangle
                hit =  in_triangle_3D( A_loc, B_loc, C_loc, p, nhat_loc);
                
                if (hit && ( (range_Ptr[ii] == 0.0) || ( fabs(r) < range_Ptr[ii] ) )) {                    
                    
                    //Save range 
                    range_Ptr[ii] = fabs(r);     
               
                }
            }

            if (range_Ptr[ii] == 0.0) {
                range_Ptr[ii] = mxGetNaN();
            }
    }
    

    
    
    
    return;
    
} // For main

///////////////////////////////////////////////////////////////////////////

////// SUBFUNCTIONS //////

void inverse ( double answer[6], double x[3], double y[3], double n[3] )
{
    double denom;
    
    denom = x[0]*(y[1]*n[2]-n[1]*y[2]) - y[0]*(x[1]*n[2]-n[1]*x[2]) + n[0]*(x[1]*y[2]-y[1]*x[2]);
    
    answer[0] = (y[1]*n[2]-n[1]*y[2])/denom;
    answer[1] = (y[2]*n[0]-n[2]*y[0])/denom;
    answer[2] = (y[0]*n[1]-n[0]*y[1])/denom;
    answer[3] = (x[2]*n[1]-n[2]*x[1])/denom;
    answer[4] = (x[0]*n[2]-n[0]*x[2])/denom;
    answer[5] = (x[1]*n[0]-n[1]*x[0])/denom;
    //    answer[6] = (y[2]*x[1]-x[2]*y[1])/denom;
    //    answer[7] = (y[0]*x[2]-x[0]*y[2])/denom;
    //    answer[8] = (y[1]*x[0]-x[1]*y[0])/denom;
    
    return;
}

double dot ( double a[3], double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    
}

int in_triangle_3D( double *A, double *B, double *C, double *P, double *nhat)
{
    // In triangle stuff
    double v0[3], v1[3], v2[3], answer[6];
    double uv[2];
    int i;
    
    for (i=0;  i<3; i++){
        v0[i] = C[i] - A[i];
        v1[i] = B[i] - A[i];
        v2[i] = P[i] - A[i];
    }
    
    inverse( answer, v1, v0, nhat );
    
    uv[0] = answer[0]*v2[0] + answer[1]*v2[1] + answer[2]*v2[2];
    uv[1] = answer[3]*v2[0] + answer[4]*v2[1] + answer[5]*v2[2];
    
    return ((uv[0] >= -1e-10) && (uv[1] >= -1e-10) && (uv[0] + uv[1] <= (1+1e-10))); 
}

double line_intersect_plane ( double *l, double *nhat, double *l0, double *p0, double *p)
{
    double denom, d, numer, diff[3];
    int ii;
    
    denom = dot(l, nhat);
    
    if ((denom < 1e-10) && (denom > -1e-10)) {
        return 0;
    } else {
        for (ii=0; ii<3; ii++) {
            diff[ii] = p0[ii] - l0[ii];
        }
        numer = dot(nhat, diff);
        d = numer/denom;
        for (ii=0; ii<3; ii++) {
            p[ii] = d*l[ii] + l0[ii];
        }
        return d;
    }
}


///////////////////////////////////////////////////////////////////////////