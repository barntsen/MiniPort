// Ac2d methods

// Imports
extern "C"
{
#include <stdlib.h>
#include <stdio.h>
}
#include "diff.h"
#include "rec.h"
#include "src.h"
#include "model.h"
#include "ac2d.h"
#include "util.h"
#include "gpu.h"
#include "hip/hip_runtime.h"

// Internal functions

  void Ac2dvx(struct ac2d *Ac2d, struct model *Model);
  void Ac2dvy(struct ac2d *Ac2d, struct model *Model);
  void Ac2dstress(struct ac2d *Ac2d, struct model *Model, int sx, int sy, float amp);
  __global__ void ac2dvx(float *Rho, float *exx, float *vx, float *thetax, float *Drhox,
                         float *Eta1x, float *Eta2x, float Dt, int nx, int ny);
  __global__ void ac2dvy(float *Rho, float *eyy, float *vy, float *thetay, float *Drhoy,
                         float *Eta1y, float *Eta2y, float Dt, int nx, int ny);
  __global__ void ac2dstress(float Dt, float *Kappa, float *p, float *gammax,
                            float *gammay, float *Dkappax, float * Dkappay,
                            float *exx, float *eyy,
                            float *Alpha1x, float *Alpha2x, float *Alpha1y, float *Alpha2y, 
                            int nx, int ny, int sx, int sy, float amp);
// Public functions

// Ac2dNew creates a new Ac2d object
//
// Parameters:
//   - Model : Model object
//
// Return    :Ac2d object  
  struct ac2d *Ac2dNew(struct model *Model)
{
  struct ac2d *Ac2d;
  int i,j;
  int Nx,Ny;

  Nx = Model->Nx;
  Ny = Model->Ny;
  
  Ac2d = (struct ac2d*)GpuNew(sizeof(struct ac2d)); 
  Ac2d->p=(float*)GpuNew(sizeof(float)*Nx*Ny);
  Ac2d->vx=(float*)GpuNew(sizeof(float)*Nx*Ny);
  Ac2d->vy=(float*)GpuNew(sizeof(float)*Nx*Ny);
  Ac2d->exx=(float*)GpuNew(sizeof(float)*Nx*Ny);
  Ac2d->eyy=(float*)GpuNew(sizeof(float)*Nx*Ny);
  Ac2d->gammax=(float*)GpuNew(sizeof(float)*Nx*Ny);
  Ac2d->gammay=(float*)GpuNew(sizeof(float)*Nx*Ny);
  Ac2d->thetax=(float*)GpuNew(sizeof(float)*Nx*Ny);
  Ac2d->thetay=(float*)GpuNew(sizeof(float)*Nx*Ny);

  hipMalloc((void**)&(Ac2d->dp),sizeof(float)*Nx*Ny);
  hipMalloc((void**)&(Ac2d->dvx),sizeof(float)*Nx*Ny);
  hipMalloc((void**)&(Ac2d->dvy),sizeof(float)*Nx*Ny);
  hipMalloc((void**)&(Ac2d->dexx),sizeof(float)*Nx*Ny);
  hipMalloc((void**)&(Ac2d->deyy),sizeof(float)*Nx*Ny);
  hipMalloc((void**)&(Ac2d->dgammax),sizeof(float)*Nx*Ny);
  hipMalloc((void**)&(Ac2d->dgammay),sizeof(float)*Nx*Ny);
  hipMalloc((void**)&(Ac2d->dthetax),sizeof(float)*Nx*Ny);
  hipMalloc((void**)&(Ac2d->dthetay),sizeof(float)*Nx*Ny);
  
  for (i=0; i<Nx; i=i+1){ 
    for (j=0; j<Ny; j=j+1){ 
      Ac2d->p[idx2(Nx,i,j)]       = 0.0;
      Ac2d->vx[idx2(Nx,i,j)]      = 0.0;
      Ac2d->vy[idx2(Nx,i,j)]      = 0.0;
      Ac2d->exx[idx2(Nx,i,j)]     = 0.0;
      Ac2d->eyy[idx2(Nx,i,j)]     = 0.0;
      Ac2d->gammax[idx2(Nx,i,j)]  = 0.0;
      Ac2d->gammay[idx2(Nx,i,j)]  = 0.0;
      Ac2d->thetax[idx2(Nx,i,j)]  = 0.0;
      Ac2d->thetay[idx2(Nx,i,j)]  = 0.0;
      Ac2d->ts = 0;
    }
  }
  hipMemcpy(Ac2d->dp, Ac2d->p,sizeof(float)*Nx*Ny,hipMemcpyHostToDevice);
  hipMemcpy(Ac2d->dvx, Ac2d->vx,sizeof(float)*Nx*Ny,hipMemcpyHostToDevice);
  hipMemcpy(Ac2d->dvy, Ac2d->vy,sizeof(float)*Nx*Ny,hipMemcpyHostToDevice);
  hipMemcpy(Ac2d->dexx, Ac2d->exx,sizeof(float)*Nx*Ny,hipMemcpyHostToDevice);
  hipMemcpy(Ac2d->deyy, Ac2d->eyy,sizeof(float)*Nx*Ny,hipMemcpyHostToDevice);
  hipMemcpy(Ac2d->dgammax, Ac2d->gammax,sizeof(float)*Nx*Ny,hipMemcpyHostToDevice);
  hipMemcpy(Ac2d->dgammay, Ac2d->gammay,sizeof(float)*Nx*Ny,hipMemcpyHostToDevice);
  hipMemcpy(Ac2d->dthetax, Ac2d->thetax,sizeof(float)*Nx*Ny,hipMemcpyHostToDevice);
  hipMemcpy(Ac2d->dthetay, Ac2d->thetay,sizeof(float)*Nx*Ny,hipMemcpyHostToDevice);
  return(Ac2d);
}
// Ac2dSolve computes the solution of the acoustic wave equation.
// The acoustic equation of motion are integrated using Virieux's (1986) stress-velocity scheme.
// (See the notes.tex file in the Doc directory).
// 
//     vx(t+dt)   = dt/rhox d^+x[ sigma(t)] + dt fx + vx(t)
//                + thetax D[1/rhox]
//     vy(t+dt)   = dt/rhoy d^+y sigma(t) + dt fy(t) + vy(t)
//                + thetay D[1/rhoy]
//
//     dp/dt(t+dt) = dt Kappa[d^-x dexx/dt + d-y deyy/dt + dt dq/dt(t) 
//                 + dt [gammax Dkappa + gammay Dkappa]
//                 + p(t)
//     dexx/dt     =  d^-_x v_x 
//     deyy/dt     =  d^-_z v_y 
//
//     gammax(t+dt) = alpha1x gammax(t) + alpha2x dexx/dt 
//     gammay(t+dt) = alpha1y gammay(t) + alpha2y deyy/dt 
//
//     thetax(t+dt) = eta1x thetax(t) + eta2x d^+x p
//     thetay(t+dt) = eta1y thetay(t) + eta2y d^+y p
//  
//  Parameters:  
//    Ac2d : Solver object
//    Model: Model object
//    Src  : Source object
//    Rec  : Receiver object
//    nt   : Number of timesteps to do starting with current step  
//    l    : The differentiator operator length
int Ac2dSolve(struct ac2d *Ac2d, struct model *Model, struct src *Src, 
              struct rec *Rec,int nt,int l)
{
  int sx,sy;         // Source x,y-coordinates 
  struct diff *Diff; // Differentiator object
  int ns,ne;         // Start stop timesteps
  int Nx,Ny;         // Model sizes
  int i,k;           // Loop indices

  float perc,oldperc; // Percentage finished current and old
  int iperc;          // Percentage finished

  Nx = Model->Nx;
  Ny = Model->Ny;

  Diff = DiffNew(l);  // Create differentiator object

  oldperc=0.0;
  ns=Ac2d->ts;         //Get current timestep 
  ne = ns+nt;         
  for(i=ns; i<ne; i=i+1){

    // Compute spatial derivative of stress
    // Use exx and eyy as temp storage
    DiffDxplus(Diff,Ac2d->dp,Ac2d->dexx,Model->Dx,Nx,Ny); // Forward differentiation x-axis
    Ac2dvx(Ac2d,Model);                        // Compute vx
    DiffDyplus(Diff,Ac2d->dp,Ac2d->deyy,Model->Dx,Nx,Ny); // Forward differentiation y-axis
    Ac2dvy(Ac2d,Model);                        // Compute vy

    DiffDxminus(Diff,Ac2d->dvx,Ac2d->dexx,Model->Dx,Nx,Ny); //Compute exx     
    DiffDyminus(Diff,Ac2d->dvy,Ac2d->deyy,Model->Dx,Nx,Ny); //Compute eyy   


    // Add source
    float amp=0.0;
    for (k=0; k<Src->Ns;k=k+1)
    {
      sx=Src->Sx[k];
      sy=Src->Sy[k];
      //Ac2d->p[idx2(Nx,sx,sy)] = Ac2d->p[idx2(Nx,sx,sy)]
      amp = Model->Dt*(Src->Src[i]/(Model->Dx*Model->Dx)) ; 
    }

    // Update stress
      Ac2dstress(Ac2d,Model,sx,sy,amp);  

    // Print progress
    perc=1000.0*(i)/(ne-ns-1);
    if(perc-oldperc >= 10.0){
      iperc=perc/10;
      if(iperc%10==0){
        fprintf(stderr,"%d\n",iperc);
        fflush(stderr);
      }
      oldperc=perc;
   }

    // Record data
    if(Rec->recon == OK){
      hipMemcpy(Ac2d->p, Ac2d->dp, sizeof(float)*Nx*Ny, hipMemcpyDeviceToHost);
      RecReceiver(Rec,i,Ac2d->p,Nx,Ny); 
    }

    // Record snapshot
    if(Rec->snpon == OK){
      RecSnap(Rec,i,Ac2d->p,Nx,Ny);
    }
  }
  return(OK);
}
// Ac2vx computes the x-component of particle velocity
//
// Parameters:
//   Ac2d : Solver object 
//   Model: Model object
void Ac2dvx(struct ac2d *Ac2d, struct model *Model)
{
  int nx,ny;

  nx = Model->Nx;
  ny = Model->Ny;
  
  ac2dvx<<< NTHREADS, NBLOCKS>>>(Model->dRho, Ac2d->dexx, Ac2d->dvx, 
                                   Ac2d->dthetax, Model->dDrhox,Model->dEta1x, 
                                   Model->dEta2x, Model->Dt, nx, ny);
  GpuError();
                       
}
__global__ void ac2dvx(float *Rho, float *exx, float *vx, float *thetax, float *Drhox,
                       float *Eta1x, float *Eta2x, float Dt, int nx, int ny)
{
  int N,p;
  int i,j;

  // The derivative of stress in x-direction is stored in exx
  // Scale with inverse density and advance one time step

  N=ny*nx; // No of processors
  for(p=blockIdx.x*blockDim.x + threadIdx.x; p<N; p+=blockDim.x*gridDim.x)
  {
    i = p%(nx);
    j = floorf(p/nx);

        vx[idx2(nx,i,j)] = Dt*(Rho[idx2(nx,i,j)])*exx[idx2(nx,i,j)] 
                               + vx[idx2(nx,i,j)]
                               + Dt*thetax[idx2(nx,i,j)]*Drhox[idx2(nx,i,j)];

        thetax[idx2(nx,i,j)]  = Eta1x[idx2(nx,i,j)]*thetax[idx2(nx,i,j)] 
                                    + Eta2x[idx2(nx,i,j)]*exx[idx2(nx,i,j)];
  }
}
// Ac2vy computes the y-component of particle velocity
//
// Parameters:
//   Ac2d : Solver object 
//   Model: Model object
void Ac2dvy(struct ac2d *Ac2d, struct model *Model)
{
  int nx,ny;

  nx = Model->Nx;
  ny = Model->Ny;
  
  ac2dvy<<< NTHREADS, NBLOCKS>>>(Model->dRho, Ac2d->deyy, Ac2d->dvy, 
                                   Ac2d->dthetay, Model->dDrhoy,Model->dEta1y, 
                                   Model->dEta2y, Model->Dt, nx, ny);
  GpuError();
                       
}
__global__ void ac2dvy(float *Rho, float *eyy, float *vy, float *thetay, float *Drhoy,
                       float *Eta1y, float *Eta2y, float Dt, int nx, int ny)
{
  int N,p;
  int i,j;

  // The derivative of stress in x-direction is stored in exx
  // Scale with inverse density and advance one time step

  N=ny*nx; // No of processors
  for(p=blockIdx.x*blockDim.x + threadIdx.x; p<N; p+=blockDim.x*gridDim.x)
  {
    i = p%(nx);
    j = floorf(p/nx);

        vy[idx2(nx,i,j)] = Dt*(Rho[idx2(nx,i,j)])*eyy[idx2(nx,i,j)] 
                               + vy[idx2(nx,i,j)]
                               + Dt*thetay[idx2(nx,i,j)]*Drhoy[idx2(nx,i,j)];

        thetay[idx2(nx,i,j)]  = Eta1y[idx2(nx,i,j)]*thetay[idx2(nx,i,j)] 
                                    + Eta2y[idx2(nx,i,j)]*eyy[idx2(nx,i,j)];
  }
}

// Ac2dstress computes acoustic stress 
//
// Parameters:
//   Ac2d : Solver object 
//   Model: Model object
void Ac2dstress(struct ac2d *Ac2d, struct model *Model, int sx, int sy, float amp)
{
  int nx, ny;

  nx = Model->Nx;
  ny = Model->Ny;

  ac2dstress<<<NTHREADS,NBLOCKS>>>(Model->Dt, Model->dKappa,Ac2d->dp, Ac2d->dgammax, 
                                  Ac2d->dgammay, Model->dDkappax, Model->dDkappay,
                                  Ac2d->dexx, Ac2d->deyy,
                                  Model->dAlpha1x,Model->dAlpha2x, Model->dAlpha1y, 
                                  Model->dAlpha2y, nx,ny,sx,sy,amp);
  GpuError();
}
 __global__ void ac2dstress(float Dt, float *Kappa, float *p, float *gammax,
                            float *gammay, float *Dkappax, float * Dkappay,
                            float *exx, float *eyy,
                            float *Alpha1x, float *Alpha2x, float *Alpha1y, float *Alpha2y, 
                            int nx, int ny, int sx, int sy,float amp)
 { 
  int N,pno;
  int i,j;

  N=ny*nx; // No of processors
  for(pno=blockIdx.x*blockDim.x + threadIdx.x; pno<N; pno+=blockDim.x*gridDim.x)
  {
    i = pno%(nx);
    j = floorf(pno/nx);

      p[idx2(nx,i,j)] = Dt*Kappa[idx2(nx,i,j)]
                             *(exx[idx2(nx,i,j)] + eyy[idx2(nx,i,j)]) 
                             + p[idx2(nx,i,j)]
                             + Dt*(gammax[idx2(nx,i,j)]*Dkappax[idx2(nx,i,j)]
                             +gammay[idx2(nx,i,j)]*Dkappay[idx2(nx,i,j)]);

      gammax[idx2(nx,i,j)] = Alpha1x[idx2(nx,i,j)]*gammax[idx2(nx,i,j)] 
                           + Alpha2x[idx2(nx,i,j)]*exx[idx2(nx,i,j)];
      gammay[idx2(nx,i,j)] = Alpha1y[idx2(nx,i,j)]*gammay[idx2(nx,i,j)] 
                                 + Alpha2y[idx2(nx,i,j)]*eyy[idx2(nx,i,j)];
      if((i==sx) && (j==sy))
      {
        p[idx2(nx,i,j)] = p[idx2(nx,i,j)]+amp;
      }
  }
}
