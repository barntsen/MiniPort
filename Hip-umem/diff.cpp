#include "hip/hip_runtime.h"
// Differentiator object.
#include "diff.h"
#include "util.h"
#include "gpu.h"

extern "C" 
{
#include <stdlib.h>
#include <stdio.h>
}

// Internal functions
 void __global__ dxminus(float *w,int l,float *A, float *dA, float dx, int nx, int ny);
 void __global__ dxplus(float *w,int l,float *A, float *dA, float dx, int nx, int ny);
 void __global__ dyminus(float *w,int l,float *A, float *dA, float dx, int nx, int ny);
 void __global__ dyplus(float *w,int l,float *A, float *dA, float dx, int nx, int ny);

//DiffNew creates a new differentiator object.
struct diff * DiffNew( int l)
{
  struct diff * Diff; 
  int i,j,k;

  Diff = (struct diff *)GpuNew(sizeof(struct diff));
  Diff->lmax = 8;

  if(l < 1){l=1;}
  if(l > Diff->lmax){l=Diff->lmax;}

  Diff->l = l;
  Diff->coeffs = (float *) GpuNew(sizeof(float)*Diff->lmax*Diff->lmax);
  Diff->w      = (float*) GpuNew(sizeof(float)*l);

  // Load coefficients
  for (j=0; j<Diff->lmax; j=j+1){
    for (i=0; i<Diff->lmax; i=i+1){
      Diff->coeffs[idx2(Diff->lmax,i,j)] = 0.0;
    } 
  }
  // l=1
  Diff->coeffs[idx2(Diff->lmax,0,0)] = 1.0021;

  // l=2
  Diff->coeffs[idx2(Diff->lmax,1,0)] = 1.1452;
  Diff->coeffs[idx2(Diff->lmax,1,1)] = -0.0492;
  
  // l=3
  Diff->coeffs[idx2(Diff->lmax,2,0)] = 1.2036;
  Diff->coeffs[idx2(Diff->lmax,2,1)] = -0.0833;
  Diff->coeffs[idx2(Diff->lmax,2,2)] = 0.0097;

  // l=4
  Diff->coeffs[idx2(Diff->lmax,3,0)] = 1.2316;
  Diff->coeffs[idx2(Diff->lmax,3,1)] = -0.1041;
  Diff->coeffs[idx2(Diff->lmax,3,2)] = 0.0206;
  Diff->coeffs[idx2(Diff->lmax,3,3)] = -0.0035;

  // l=5
  Diff->coeffs[idx2(Diff->lmax,4,0)] = 1.2463;
  Diff->coeffs[idx2(Diff->lmax,4,1)] = -0.1163;
  Diff->coeffs[idx2(Diff->lmax,4,2)] = 0.0290;
  Diff->coeffs[idx2(Diff->lmax,4,3)] = -0.0080;
  Diff->coeffs[idx2(Diff->lmax,4,4)] = 0.0018;

  // l=6
  Diff->coeffs[idx2(Diff->lmax,5,0)] = 1.2542;
  Diff->coeffs[idx2(Diff->lmax,5,1)] = -0.1213;
  Diff->coeffs[idx2(Diff->lmax,5,2)] = 0.0344;
  Diff->coeffs[idx2(Diff->lmax,5,3)] = -0.017;
  Diff->coeffs[idx2(Diff->lmax,5,4)] = 0.0038;
  Diff->coeffs[idx2(Diff->lmax,5,5)] = -0.0011;

  // l=7
  Diff->coeffs[idx2(Diff->lmax,6,0)] = 1.2593;
  Diff->coeffs[idx2(Diff->lmax,6,1)] = -0.1280;
  Diff->coeffs[idx2(Diff->lmax,6,2)] = 0.0384;
  Diff->coeffs[idx2(Diff->lmax,6,3)] = -0.0147;
  Diff->coeffs[idx2(Diff->lmax,6,4)] = 0.0059;
  Diff->coeffs[idx2(Diff->lmax,6,5)] = -0.0022;
  Diff->coeffs[idx2(Diff->lmax,6,6)] = 0.0007;

  // l=8
  Diff->coeffs[idx2(Diff->lmax,7,0)] = 1.2626;
  Diff->coeffs[idx2(Diff->lmax,7,1)] = -0.1312;
  Diff->coeffs[idx2(Diff->lmax,7,2)] = 0.0412;
  Diff->coeffs[idx2(Diff->lmax,7,3)] = -0.0170;
  Diff->coeffs[idx2(Diff->lmax,7,4)] = 0.0076;
  Diff->coeffs[idx2(Diff->lmax,7,5)] = -0.0034;
  Diff->coeffs[idx2(Diff->lmax,7,6)] = 0.0014;
  Diff->coeffs[idx2(Diff->lmax,7,7)] = -0.0005;


  for(k=0;k<l;k=k+1)
  {
    Diff->w[k] = Diff->coeffs[idx2(Diff->lmax,l-1,k)];
  }

  return(Diff);
}
// Dxminus computes the backward derivative in the x-direction.
//
//  Parameters:
//  Diff           : Diff object 
//  float  A       : Input 2D array
//  float dx       : Sampling interval
//  float dA       : Output array 
//
//  The output array, dA, contains the derivative for each point computed
//  as:
//  dA[i,j] = (1/dx) sum_{k=1}^l w[k](A[i+(k-1)dx,j]-A[(i-kdx,j]
//
//  w[k] are weights and l is the length of the differentiator.
//  (see DiffNew for the definitions of these)
void DiffDxminus(struct diff *Diff , float *A, float *dA, float dx, int nx, int ny)
{
  int l;
  float *w;

  l= Diff->l;
  w = Diff->w;
  // Kernel call 
  hipLaunchKernelGGL(dxminus, NBLOCKS, NTHREADS, 0, 0, w,l,A,dA,dx,nx,ny);
  GpuError();
}

// dxminus computes backward derivative in x-direction
 __global__ void dxminus(float *w,int l,float *A, float *dA, float dx, int nx, int ny)
{
  int i,j,k;
  float sum;
  int p;
  int N;

  //
  // Left border (1 <i < l+1)
  //

//  for(j=0;j<ny;j=j+1)
//  {
//    for(i=0;i<l;i=i+1)
//    {

  N=l*ny; // No of processors
  for(p=blockIdx.x*blockDim.x + threadIdx.x; p<N; p+=blockDim.x*gridDim.x)
  {
    i = p%(l);
    j = p/l;

    sum=0.0;
    for(k=1; k<i+1; k=k+1)
    {
      sum = -w[k-1]*A[idx2(nx,i-k,j)] + sum; 
    }
    for(k=1; k<l+1; k=k+1)
    {
      sum = w[k-1]*A[idx2(nx,i+(k-1),j)] +sum; 
    }
    dA[idx2(nx,i,j)] = sum/dx;
  }
//    } 
//  }

  //
  // Outside border area 
  //
  //for(j=0; j<ny; j=j+1)
  //{
  //  for(i=l; i<nx-l; i=i+1)
  //  {

  N=(nx-2*l)*ny; // No of processors
  for(p=blockIdx.x*blockDim.x + threadIdx.x; p<N; p+=blockDim.x*gridDim.x)
  {
    i = l+p%(nx-2*l);
    j = p/(nx-2*l);

    sum=0.0;
    for(k=1; k<l+1; k=k+1)
    {
      sum = w[k-1]*(-A[idx2(nx,i-k,j)]+A[idx2(nx,i+(k-1),j)]) +sum; 
    }
    dA[idx2(nx,i,j)] = sum/dx;
  }
   // }
  //} 

  //
  // Right border 
  //
//  for(j=0; j<ny; j=j+1)
//  {
//    for(i=nx-l; i<nx; i=i+1)
//    {

  N=l*ny; // No of processors
  for(p=blockIdx.x*blockDim.x + threadIdx.x; p<N; p+=blockDim.x*gridDim.x)
  {
    i = (nx-l) + p%(l);
    j = p/l;

    sum = 0.0;
    for(k=1; k<l+1; k=k+1)
    {
      sum = -w[k-1]*A[idx2(nx,i-k,j)] + sum;
    }

    for(k=1; k<(nx-i+1); k=k+1)
    {
      sum = w[k-1]*A[idx2(nx,i+(k-1),j)] + sum;
    }
    dA[idx2(nx,i,j)] = sum/dx;
  }
  //  }
 // }

}

// Dxplus computes the forward derivative in the x-direction.
//
// Arguments:
//  Diff: Diff object 
//  float  A       : Input 2D array
//  float dx       : Sampling interval
//  float dA       : Output array 
//
//  The output array, dA, contains the derivative for each point computed
//  as:
//  dA[i,j] = (1/dx) sum_{k=1}^l w[k](A[i+kdx,j]-A[(i-(k-1)dx,j]
//
//  w[k] are weights and l is the length of the differentiator.
//  (see DiffNew for the definitions of these)
//------------------------------------------------------------------------------
void DiffDxplus(struct diff *Diff , float *A, float *dA, float dx, int nx, int ny)
{
  int l;
  float *w;

  l= Diff->l;
  w = Diff->w;
  // Kernel call 
  hipLaunchKernelGGL(dxplus, NBLOCKS, NTHREADS, 0, 0, w,l,A,dA,dx,nx,ny);
  GpuError();
}

__global__  void dxplus(float *w, int l, float *A, float *dA, float dx, int nx, int ny)
{
  int i,j,k;
  float sum;
  int N,p;

  // Left border

  //for(j=0; j<ny; j=j+1)
  //{
  //  for(i=0; i<l; i=i+1)
  //  {

  N=l*ny; // No of processors
  for(p=blockIdx.x*blockDim.x + threadIdx.x; p<N; p+=blockDim.x*gridDim.x)
  {
    i = p%(l);
    j = p/l;

    sum=0.0;
    for(k=1; k<i+2; k=k+1)
    {
      sum = -w[k-1]*A[idx2(nx,i-(k-1),j)] + sum; 
    }
    for(k=1; k<l+1; k=k+1)
    {
      sum = w[k-1]*A[idx2(nx,i+k,j)] +sum; 
    }
    dA[idx2(nx,i,j)] = sum/dx;
  }
   // }
//  } 
  //
  // Between left and right border
  //
  //for(j=0; j<ny; j=j+1)
  //{
  //  for(i=l; i<nx-l; i=i+1)
   // {
  N=(nx-2*l)*ny; // No of processors
  for(p=blockIdx.x*blockDim.x + threadIdx.x; p<N; p+=blockDim.x*gridDim.x)
  {
    i = l+p%(nx-2*l);
    j = p/(nx-2*l);
      sum=0.0;
      for(k=1; k<l+1; k=k+1)
      {
        sum = w[k-1]*(-A[idx2(nx,i-(k-1),j)]+A[idx2(nx,i+k,j)]) +sum; 
      }
      dA[idx2(nx,i,j)] = sum/dx;
   }
   // }
  //} 

  //
  // Right border 
  //
  //for(j=0; j<ny; j=j+1)
  //{
    //for(i=nx-l; i<nx; i=i+1)
    //{
    //N=l*ny-1; // No of processors
    N=l*ny-1; // No of processors
    for(p=blockIdx.x*blockDim.x + threadIdx.x; p<N; p+=blockDim.x*gridDim.x)
    {
      i = (nx-l) + p%(l);
      j = p/l;

      sum = 0.0;
      for(k=1; k<l+1; k=k+1)
      {
        sum = -w[k-1]*A[idx2(nx,i-(k-1),j)] + sum;
      }

      for(k=1; k<nx-i; k=k+1)
      {
        sum = w[k-1]*A[idx2(nx,i+k,j)] + sum;
      }
      dA[idx2(nx,i,j)] = sum/dx;
    }
    //}
  //}
}
//------------------------------------------------------------------------------
// Dyminus computes the backward derivative in the y-direction
//
// Arguments:
//  Diff: Diff object 
//  float  A       : Input 2D array
//  float dx       : Sampling interval
//  float dA       : Output array 
//
//  The output array, dA, contains the derivative for each point computed
//  as:
//  dA[i,j] = (1/dx) sum_{k=1}^l w[k](A[i,j+(k-1)dx]-A[i,j-kdx,j]
//
//  w[k] are weights and l is the length of the differentiator.
//  (see DiffNew for the definitions of these)
//------------------------------------------------------------------------------
void DiffDyminus(struct diff *Diff , float *A, float *dA, float dx, int nx, int ny)
{
  int l;
  float *w;

  l= Diff->l;
  w = Diff->w;
  // Kernel call 
  hipLaunchKernelGGL(dyminus, NBLOCKS, NTHREADS, 0, 0, w,l,A,dA,dx,nx,ny);
  GpuError();
}

void __global__ dyminus(float *w, int l,float *A, float *dA, float dx, int nx, int ny){
  int i,j,k;
  float sum;
  int N,p;

  // Left border 

//  for(j=0; j<l; j=j+1)
//  {
//    for(i=0; i<nx; i=i+1)
//    {

  N=l*nx; // No of processors
  for(p=blockIdx.x*blockDim.x + threadIdx.x; p<N; p+=blockDim.x*gridDim.x)
  {
    i = p%(nx);
    j = p/nx;

    sum=0.0;
    for(k=1; k<j+1; k=k+1)
    {
      sum = -w[k-1]*A[idx2(nx,i,j-k)] + sum; 
    }
    for(k=1; k<l+1; k=k+1)
    {
      sum = w[k-1]*A[idx2(nx,i,j+(k-1))] +sum; 
    }
    dA[idx2(nx,i,j)] = sum/dx;
  }
    //}
  //} 
  //
  // Outside border area 
  //
  //for(j=l; j<ny-l; j=j+1)
  //{
  //  for(i=0; i<nx; i=i+1)
  //  {
  N=(ny-2*l)*nx; // No of processors
  for(p=blockIdx.x*blockDim.x + threadIdx.x; p<N; p+=blockDim.x*gridDim.x)
  {
    i = p%nx;
    j = l+p/nx;
    sum=0.0;
    for(k=1; k<l+1; k=k+1)
    {
      sum = w[k-1]*(-A[idx2(nx,i,j-k)]+A[idx2(nx,i,j+(k-1))]) +sum; 
    }
    dA[idx2(nx,i,j)] = sum/dx;
  }
    //}
  //} 

  //
  // Right border 
  //
  //for(j=ny-l; j<ny; j=j+1)
  //{
  //  for(i=0; i<nx; i=i+1)
  //  {
    N=(l)*nx; // No of processors
    for(p=blockIdx.x*blockDim.x + threadIdx.x; p<N; p+=blockDim.x*gridDim.x)
    {
      i = p%nx;
      j = ny-l +p/nx;

      sum = 0.0;
      for(k=1; k<l+1; k=k+1)
      {
        sum = -w[k-1]*A[idx2(nx,i,j-k)] + sum;
      }

      for(k=1; k<(ny-j+1); k=k+1)
      {
        sum = w[k-1]*A[idx2(nx,i,j+(k-1))] + sum;
      }
      dA[idx2(nx,i,j)] = sum/dx;
    }
   // }
  //}
}
// Dyplus computes the forward derivative in the x-direction
//
// Arguments:
//  Diff: Diff object 
//  float  A       : Input 2D array
//  float dx       : Sampling interval
//  float dA       : Output array 
//
//  The output array, dA, contains the derivative for each point computed
//  as:
//  dA[i,j] = (1/dx) sum_{k=1}^l w[k](A[i+kdx,j]-A[(i-(k-1)dx,j]
//
//  w[k] are weights and l is the length of the differentiator.
//  (see DiffNew for the definitions of these)
//------------------------------------------------------------------------------
void DiffDyplus(struct diff *Diff , float *A, float *dA, float dx, int nx, int ny)
{
  int l;
  float *w;

  l= Diff->l;
  w = Diff->w;
  // Kernel call 
  hipLaunchKernelGGL(dyplus, NBLOCKS, NTHREADS, 0, 0, w,l,A,dA,dx,nx,ny);
  GpuError();
}
__global__ void dyplus(float *w, int l, float *A, float *dA, float dx, int nx, int ny)
{
  int i,j,k;
  float sum;
  int N,p;

  //
  // Left border (1 <i < l+1)
  //

  // Left border
  //for(j=0; j<l; j=j+1)
  //{
    //for(i=0; i<nx; i=i+1)
    //{
  N=l*nx; // No of processors
  for(p=blockIdx.x*blockDim.x + threadIdx.x; p<N; p+=blockDim.x*gridDim.x)
  {
    i = p%(nx);
    j = p/nx;
      sum=0.0;
      for(k=1; k<j+2; k=k+1)
      {
        sum = -w[k-1]*A[idx2(nx,i,j-(k-1))] + sum; 
      }
      for(k=1; k<l+1; k=k+1)
      {
        sum = w[k-1]*A[idx2(nx,i,j+k)] +sum; 
      }
      dA[idx2(nx,i,j)] = sum/dx;
   }
 //   }
 // } 

  //
  // Between left and right border
  //
  //for(j=l; j<ny-l; j=j+1)
  //{
  //  for(i=0; i<nx; i=i+1)
  //  {
  N=(ny-2*l)*nx; // No of processors
  for(p=blockIdx.x*blockDim.x + threadIdx.x; p<N; p+=blockDim.x*gridDim.x)
  {
    i = p%nx;
    j = l+p/nx;

      sum=0.0;
      for(k=1; k<l+1; k=k+1)
      {
        sum = w[k-1]*(-A[idx2(nx,i,j-(k-1))]+A[idx2(nx,i,j+k)]) +sum; 
      }
      dA[idx2(nx,i,j)] = sum/dx;
  } 
  //  } 
  //}
  //
  // Right border 
  //
  //for(j=ny-l; j<ny; j=j+1)
  //{
  //  for(i=0; i<nx; i=i+1)
  //  {
  N=l*nx-1; // No of processors
  for(p=blockIdx.x*blockDim.x + threadIdx.x; p<N; p+=blockDim.x*gridDim.x)
  {
      i = p%nx;
      j = ny-l +p/nx;

      sum = 0.0;
      for(k=1; k<l+1; k=k+1)
      {
        sum = -w[k-1]*A[idx2(nx,i,j-(k-1))] + sum;
      }

      for(k=1; k<ny-j; k=k+1){
        sum = w[k-1]*A[idx2(nx,i,j+k)] + sum;
      }
      dA[idx2(nx,i,j)] = sum/dx;
   }
   // }
  //}
}
