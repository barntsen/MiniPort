// Methods for the model struct
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "model.h" // Model struct definition
#include "util.h"
#include "gpu.h"

struct model *Modelmaxwell(float *vp, float *rho, float *Q, 
                      float Dx, float Dt, float W0, int Nb, int Nx, int Ny);

struct model *Modelsls(float *vp, float *rho, float *Q, 
                      float Dx, float Dt, float W0, int Nb,int Nx, int Ny);

int Modeld(float *d, float dx, int nb, int n);


// ModelNew creates a new model.
//
// Parameters: 
//
//  - vp :  P-wave velocity model
//  - rho:  Density 
//  - Q  :  Q-values
//  - Dx :  Grid interval in x- and y-directions
//  - Dt :  Modeling time sampling interval
//  - W0 :  Q-model peak angular frequency
//  - Nb :  Width of border attenuation zone (in grid points)
//  - Rheol : Type of Q-model. Rheol=MAXWELL (Maxwell solid)
//                             Rheol=SLS     (Standard linear solid)
//
// Return:  Model structure
//
// ModelNew creates the parameters needed by the Ac2d object
// to perform 2D acoustic modeling.
// For the details of the MAXWELL or SLS type models
// see the comments in Modelmaxwell and Modelsls.


struct model *ModelNew(float *vp, float *rho, float *Q, 
                      float Dx, float Dt, float W0, int Nb, int Rheol, int Nx, int Ny)
{
  struct model *m;

  if(Rheol == MAXWELL)
  {
    m = Modelmaxwell(vp, rho, Q, Dx, Dt, W0, Nb, Nx, Ny); 
  } else if(Rheol == SLS)
  {
    m= Modelsls(vp, rho, Q, Dx, Dt, W0, Nb, Nx, Ny);
  }
  else
  {
    fprintf(stderr,"Uknown Q-model\n"); 
    // Bailing out
    exit(ERR);
  } 
  
  return(m);
}
// Modelmaxwell creates a new model.
//
// Parameters: 
//
//  - vp :  P-wave velocity model
//  - rho:  Density 
//  - Q  :  Q-values
//  - Dx :  Grid interval in x- and y-directions
//  - Dt :  Modeling time sampling interval
//  - W0 :  Q-model peak angular frequency
//  - Nb :  Width of border attenuation zone (in grid points)
//
// Return:  Model structure
//
// ModelNew creates the parameters needed by the Ac2d object
// to perform 2D acoustic modeling.
// The main parameters are density $\rho$ and bulk modulus $\kappa$ which are
// calculated from the wave velocity and density.
// In addition are the visco-elastic coefficients $\alpha_1$, $\alpha_2$ ,
// $\eta_1$  and $\eta_2$ computed.
//
// The model is defined by several 2D arrays, with the x-coordinate
// as the first index, and the y-coordinate as the second index.
// A position in the model (x,y) maps to the arrays as [i,j]
// where x=Dx*i, y=Dx*j
// The absorbing boundaries is comparable to the CPML method
// but constructed using a visco-elastic medium with
// relaxation specified by a standard-linear solid, while 
// a time dependent density which uses a standard-linear solid
// relaxation mechanism.
//
//                     Nx                Outer border        
//    |----------------------------------------------|
//    |           Qmin=1.1                           |
//    |                                              |
//    |           Qmax=Q(x,y=Dx*Nb)     Inner border |
//    |      ----------------------------------      |
//    |      |                                |      |
//    |      |                                |      | Ny
//    |      |      Q(x,y)                    |      |
//    |      |                                |<-Nb->|
//    |      |                                |      |
//    |      |                                |      |
//    |      ----------------------------------      |
//    |                                              |
//    |                                              |
//    |                                              |
//    |-----------------------------------------------
//
//    Fig 1: Organisation of the Q-model.
//           The other arrays are organised in the same way.
//
// The Boundary condition is implemented by using a strongly
// absorbing medium in a border zone with width Nb.
// The border zone has the same width both in the horizontal
// and vertical directions.
// The medium in the border zone has a Q-value of Qmax
// at the inner bondary (taken from the Q-model) and
// the Q-value is gradualy reduced to Qmin at the outer boundary.
//
//  In the finit-edifference solver we use the Maxwell 
//  solid to implement time dependent 
//  bulk modulus and density.
//  The Maxwell solid model uses
//  one parameter, tau0.
//  tau0 is related to the Q-value by
// (See the notes.tex in the Doc directory for the equations.)
//  
//    taue(Q0) = Q(W0)/W0
//
//  Q0 is here the value for Q at the frequency W0.
//
//  The coeffcients needed by the solver methods in the Ac2d object are
//    alpha1x =  exp(d_x/Dt)exp(tau0x),                                  \\
//    alpha2x =  - dx Dt/tau0x
//    alpha1y =  exp(d_x/Dt)exp(tau0y),                                  \\
//    alpha2y =  -dx Dt/tau0y
//    eta1x   =  exp(d_x/Dt)exp(tau0x),                                  \\
//    eta2x   =  -dx Dt/tau0x
//    eta1y   =  exp(d_x/Dt)exp(tau0y),                                  \\
//    eta2y   =  -dx Dt/tau0y
//
// tau0 is interpolated between the values given by the Q-value 
// Qmax at the inner border of the model and the Qmin at the outer border. 
// For the interpolation we just assume that the relaxation times
// varies proportionaly with the square of the distance from
// the inner border, according to
//
//   tau0x(x) = tau0xmin + (tau0xmax-tau0xmin)*d(x)
//   tau0y(x) = tau0ymin + (tau0xmax-tau0ymin)*d(y)
//                       
// where 
//
//   d(x) = (x/L)^2
//
// x is the distance from the outer border, while
// L is the length of the border.
// We also have
//
//   tau0xmax = tau0(Qmax)
//   tau0xmin = tau0(Qmin)
//   tau0ymax = tau0(Qmax)
//   tau0ymin = tau0(Qmin)
//
// Here Qmin= 1.1, while Qmax is equal to the value 
// of Q at the inner border.
// 
struct model *Modelmaxwell(float *vp, float *rho, float *Q, 
                      float Dx, float Dt, float W0, int Nb, int Nx, int Ny)
{
  struct model *Model; // Object to instantiate

  // Smoothing parameters
  float Qmin, Qmax;       // Minimum and Maximum Q-values in boundary zone
  float tau0min,tau0max;  // Taue values corresponding to Qmin and Qmax

  // Relaxation times
  float tau0x, tau0y;     

  float argx;            // Temp variabels
  float argy;            // Temp variables
  int i,j;               // Loop indices

  Model= (struct model*)GpuNew(sizeof(struct model));
  Model->Dx = Dx;
  Model->Dt = Dt;
  Model->Nx = Nx;
  Model->Ny = Ny;
  Model->Nb = Nb;
  Model->W0 = W0;
  Nx = Model->Nx;
  Ny = Model->Ny;
  Model->Rho   = (float*)GpuNew(sizeof(float)*Nx*Ny);
  Model->Q     = (float*)GpuNew(sizeof(float)*Nx*Ny);
  Model->Kappa = (float*)GpuNew(sizeof(float)*Nx*Ny);

  // The following parameters are the change in the 
  // bulk modulus caused by visco-elasticity
  // A separate factor is used for the x- and y-directions
  // due to tapering
  Model->Dkappax = (float*)GpuNew(sizeof(float)*Nx*Ny);
  Model->Dkappay = (float*)GpuNew(sizeof(float)*Nx*Ny);
  Model->Drhox   = (float*)GpuNew(sizeof(float)*Nx*Ny);
  Model->Drhoy   = (float*)GpuNew(sizeof(float)*Nx*Ny);

  // Coeffcients used for updating memory functions
  Model->Alpha1x   =  (float*)GpuNew(sizeof(float)*Nx*Ny);
  Model->Alpha1y   =  (float*)GpuNew(sizeof(float)*Nx*Ny);
  Model->Alpha2x   =  (float*)GpuNew(sizeof(float)*Nx*Ny);
  Model->Alpha2y   =  (float*)GpuNew(sizeof(float)*Nx*Ny);
  Model->Eta1x   =  (float*)GpuNew(sizeof(float)*Nx*Ny);
  Model->Eta1y   =  (float*)GpuNew(sizeof(float)*Nx*Ny);
  Model->Eta2x   =  (float*)GpuNew(sizeof(float)*Nx*Ny);
  Model->Eta2y   =  (float*)GpuNew(sizeof(float)*Nx*Ny);

  // Tapering (profile) functions for
  // the x- and y-directions.
  Model->dx   =  (float*)GpuNew(sizeof(float)*Nx);
  Model->dy   =  (float*)GpuNew(sizeof(float)*Ny);

  // Store the model
  for(j=0; j<Ny;j=j+1){
    for(i=0; i<Nx;i=i+1){
      Model->Kappa[idx2(Nx,i,j)] = rho[idx2(Nx,i,j)]*vp[idx2(Nx,i,j)]*vp[idx2(Nx,i,j)];
      Model->Rho[idx2(Nx,i,j)]   = 1.0/rho[idx2(Nx,i,j)];
      Model->Q[idx2(Nx,i,j)]       = Q[idx2(Nx,i,j)];
    }
  }

  //Compute 1D profile functions
    Modeld(Model->dx, Model->Dx, Model->Nb,Nx);
    Modeld(Model->dy, Model->Dx, Model->Nb,Ny);
 
  // Compute relaxation times
  for(j=0; j<Ny;j=j+1){
    for(i=0; i<Nx;i=i+1){

      // Compute relaxation times corresponding to Qmax and Qmin
      // Note that we compute the inverse
      // of tau0, and use the same
      // name for the inverse, tau0=1/tau0.
      
      Qmin = 1.1;  // MinimumQ-value at the outer boundaries:       
      tau0min = Qmin/Model->W0;
      tau0min = 1.0/tau0min;
      Qmax  = Model->Q[idx2(Nx,Nb,j)];
      tau0max = Qmax/Model->W0;
      tau0max = 1.0/tau0max;

      // Interpolate tau0 in x-direxction
      tau0x = tau0min + (tau0max-tau0min)*Model->dx[i];

      Qmax  = Model->Q[idx2(Nx,i,Nb)];
      tau0max = Qmax/Model->W0;
      tau0max = 1.0/tau0max;

      // Interpolate tau0 in y-direxction
      tau0y = tau0min + (tau0max-tau0min)*Model->dy[j];

      // In the equations below the relaxation time tau0 
      // is inverse (1/tau0)
      // Compute alpha and eta coefficients
      argx = Model->dx[i];
      argy = Model->dy[j];
      // An extra tapering factor of exp(-(x/L)**2)
      // is used to taper some coefficeints 
      Model->Alpha1x[idx2(Nx,i,j)]   = expf(-argx)*expf(-Model->Dt*tau0x);
      Model->Alpha1y[idx2(Nx,i,j)]   = expf(-argy)*expf(-Model->Dt*tau0y);
      Model->Alpha2x[idx2(Nx,i,j)]   = -Model->Dt*tau0x;
      Model->Alpha2y[idx2(Nx,i,j)]   = -Model->Dt*tau0y;
      Model->Eta1x[idx2(Nx,i,j)]     = expf(-argx)*expf(-Model->Dt*tau0x);
      Model->Eta1y[idx2(Nx,i,j)]     = expf(-argy)*expf(-Model->Dt*tau0y);
      Model->Eta2x[idx2(Nx,i,j)]     = -Model->Dt*tau0x;
      Model->Eta2y[idx2(Nx,i,j)]     = -Model->Dt*tau0y;
 
      // For the Maxwell solid Dkappa = kappa and Drho = 1/rho
      // to comply with the solver algorithm in ac2d.e
      Model->Dkappax[idx2(Nx,i,j)]   = Model->Kappa[idx2(Nx,i,j)];
      Model->Dkappay[idx2(Nx,i,j)]   = Model->Kappa[idx2(Nx,i,j)];
      Model->Drhox[idx2(Nx,i,j)]     = Model->Rho[idx2(Nx,i,j)];
      Model->Drhoy[idx2(Nx,i,j)]     = Model->Rho[idx2(Nx,i,j)];
    }
  }

  return(Model);
}
// Modelsls creates a new model with Standard Linear Solid Q
//
// Parameters: 
//
//  - vp :  P-wave velocity model
//  - rho:  Density 
//  - Q  :  Q-values
//  - Dx :  Grid interval in x- and y-directions
//  - Dt :  Modeling time sampling interval
//  - W0 :  Q-model peak angular frequency
//  - Nb :  Width of border attenuation zone (in grid points)
//
// Return:  Model structure
//
// ModelNew creates the parameters needed by the Ac2d object
// to perform 2D acoustic modeling.
// The main parameters are density $\rho$ and bulk modulus $\kappa$ which are
// calculated from the wave velocity and density.
// In addition are the visco-elastic coefficients $\alpha_1$, $\alpha_2$ ,
// $\eta_1$  and $\eta_2$ computed.
//
// The model is defined by several 2D arrays, with the x-coordinate
// as the first index, and the y-coordinate as the second index.
// A position in the model (x,y) maps to the arrays as [i,j]
// where x=Dx*i, y=Dx*j
// The absorbing boundaries is comparable to the CPML method
// but constructed using a visco-elastic medium with
// relaxation specified by a standard-linear solid, while 
// a time dependent density which uses a standard-linear solid
// relaxation mechanism.
//
//                     Nx                Outer border        
//    |----------------------------------------------|
//    |           Qmin=1.1                           |
//    |                                              |
//    |           Qmax=Q(x,y=Dx*Nb)     Inner border |
//    |      ----------------------------------      |
//    |      |                                |      |
//    |      |                                |      | Ny
//    |      |      Q(x,y)                    |      |
//    |      |                                |<-Nb->|
//    |      |                                |      |
//    |      |                                |      |
//    |      ----------------------------------      |
//    |                                              |
//    |                                              |
//    |                                              |
//    |-----------------------------------------------
//
//    Fig 1: Organisation of the Q-model.
//           The other arrays are organised in the same way.
//
// The Boundary condition is implemented by using a strongly
// absorbing medium in a border zone with width Nb.
// The border zone has the same width both in the horizontal
// and vertical directions.
// The medium in the border zone has a Q-value of Qmax
// at the inner bondary (taken from the Q-model) and
// the Q-value is gradualy reduced to Qmin at the outer boundary.
//
//  In the finit-edifference solver we use the standard
//  linear solid to implement time dependent 
//  bulk modulus and density.
//  The standard linear solid model uses
//  two parameters, $\tau_{sigma}$ and $\tau_{\epsilon}$.
//  These are related to the Q-value by
// (See the notes.tex in the Doc directory for the equations.)
//  
//    taue(Q0) = tau0/Q0(\sqrt{Q^2_0+1} +1\right)
//    taus(Q0) = tau0/Q0(\sqrt{Q^2_0+1} +1\right)
//
//  Q0 is here the value for Q at the frequency W0.
//
//  The coeffcients needed by the solver methods in the Ac2d object are
//    alpha1x =  exp(d_x/Dt)exp(tausx),                                  \\
//    alpha2x =  dx Dt/tauex
//    alpha1y =  exp(d_x/Dt)exp(tausy),                                  \\
//    alpha2y =  dx Dt/tauey
//    eta1x   =  exp(d_x/Dt)exp(tausx),                                  \\
//    eta2x   =  dx Dt/tauex
//    eta1y   =  exp(d_x/Dt)exp(tausy),                                  \\
//    eta2y   =  dx Dt/tauey
//
// Relaxation times are interpolated between the values given by the Q-value 
// Qmax at the inner border of the model and the Qimin at the outer border. 
// For the interpolation we just assume that the relaxation times
// varies proportionaly with the square of the distance from
// the inner border, according to
//
//   tausx(x) = tausxmin + (tausxmax-tausxmin)*d(x)
//   tausy(x) = tausymin + (tausxmax-tausymin)*d(y)
//   tauex(x) = tauexmin + (tauexmax-tauexmin)*d(x)
//   tauey(x) = taueymin + (tausymax-tausymin)*d(y)
//                       
// where 
//
//   d(x) = (x/L)^2
//
// x is the distance from the outer border, while
// L is the length of the border.
// We also have
//
//   tausxmax = taus(Qmax)
//   tausxmin = taus(Qmin)
//   tausymax = taus(Qmax)
//   tausymin = taus(Qmin)
//   tauexmax = taue(Qmax)
//   tauexmin = taue(Qmin)
//   taueymax = taue(Qmax)
//   taueymin = taue(Qmin)
//
// Here Qmin= 1.1, while Qmax is equal to the value 
// of Q at the inner border.
struct model *Modelsls(float *vp, float *rho, float *Q, 
                      float Dx, float Dt, float W0, int Nb,int Nx, int Ny)
{
  struct model *Model; // Object to instantiate

  float tau0;         // Relaxation time at Peak 1/Q-value
  
  // Smoothing parameters
  float Qmin, Qmax;       // Minimum and Maximum Q-values in boundary zone
  float tauemin,tauemax;  // Taue values corresponding to Qmin and Qmax
  float tausmin,tausmax;  // Taus values corresponding to Qmin and Qmax

  // Relaxation times
  float tausx, tausy;     
  float tauex, tauey;

  float argx;            // Temp variabels
  float argy;            // Temp variables
  int i,j;               // Loop indices

  Model= (struct model*)GpuNew(sizeof(struct model));
  Model->Dx = Dx;
  Model->Dt = Dt;
  Model->Nx = Nx;
  Model->Ny = Ny;
  Model->Nb = Nb;
  Model->W0 = W0;
  Nx = Model->Nx;
  Ny = Model->Ny;
  Model->Rho     =  (float*)GpuNew(sizeof(float)*Nx*Ny); // Density
  Model->Q       =  (float*)GpuNew(sizeof(float)*Nx*Ny); // Density
  Model->Kappa   = (float*)GpuNew(sizeof(float)*Nx*Ny); // Density

  // The following parameters are the change in the 
  // bulk modulus caused by visco-elasticity
  // A separate factor is used for the x- and y-directions
  // due to tapering
  Model->Dkappax = (float*)GpuNew(sizeof(float)*Nx*Ny);  
  Model->Dkappay = (float*)GpuNew(sizeof(float)*Nx*Ny);  
  Model->Drhox     = (float*)GpuNew(sizeof(float)*Nx*Ny);  
  Model->Drhoy     = (float*)GpuNew(sizeof(float)*Nx*Ny);  

  // Coeffcients used for updating memory functions
  Model->Alpha1x   =  (float*)GpuNew(sizeof(float)*Nx*Ny);
  Model->Alpha1y   =  (float*)GpuNew(sizeof(float)*Nx*Ny);
  Model->Alpha2x   =  (float*)GpuNew(sizeof(float)*Nx*Ny);
  Model->Alpha2y   =  (float*)GpuNew(sizeof(float)*Nx*Ny);
  Model->Eta1x   =  (float*)GpuNew(sizeof(float)*Nx*Ny);
  Model->Eta1y   =  (float*)GpuNew(sizeof(float)*Nx*Ny);
  Model->Eta2x   =  (float*)GpuNew(sizeof(float)*Nx*Ny);
  Model->Eta2y   =  (float*)GpuNew(sizeof(float)*Nx*Ny);

  // Tapering (profile) functions for
  // the x- and y-directions.
  Model->dx   =  (float*)GpuNew(sizeof(float)*Nx);
  Model->dy   =  (float*)GpuNew(sizeof(float)*Ny);

  // Store the model
  for(j=0; j<Ny;j=j+1){
    for(i=0; i<Nx;i=i+1){
      Model->Kappa[idx2(Nx,i,j)] = rho[idx2(Nx,i,j)]*vp[idx2(Nx,i,j)]*vp[idx2(Nx,i,j)];
      Model->Rho[idx2(Nx,i,j)]   = 1.0/rho[idx2(Nx,i,j)];
      Model->Q[idx2(Nx,i,j)]       = Q[idx2(Nx,i,j)];
    }
  }

  //Compute 1D profile functions
    Modeld(Model->dx, Model->Dx, Model->Nb,Nx);
    Modeld(Model->dy, Model->Dx, Model->Nb,Ny);
 
  // Compute relaxation times
  for(j=0; j<Ny;j=j+1){
    for(i=0; i<Nx;i=i+1){
      tau0 = 1.0/Model->W0;   // Relaxation time corresponding to absorption top
      Qmin = 1.1;            // MinimumQ-value at the outer boundaries

      // Compute relaxation times corresponding to Qmax and Qmin
      tauemin = (tau0/Qmin)*(sqrtf(Qmin*Qmin+1.0)+1.0);
      tauemin = 1.0/tauemin;
      tausmin = (tau0/Qmin)*(sqrtf(Qmin*Qmin+1.0)-1.0);
      tausmin = 1.0/tausmin;

      Qmax  = Model->Q[idx2(Nx,Nb,j)];
      // Note that we compute the inverse
      // of relaxation times, and use the same
      // name for the inverses, taus=1/taus.
      // In all formulas below this section we
      // work with the inverse of the relaxation times.
      tauemax = (tau0/Qmin)*(sqrtf(Qmax*Qmax+1.0)+1.0);
      tauemax = 1.0/tauemax;
      tausmax = (tau0/Qmin)*(sqrtf(Qmax*Qmax+1.0)-1.0);
      tausmax = 1.0/tausmax;
      tauex = tauemin + (tauemax-tauemin)*Model->dx[i];
      tausx = tausmin + (tausmax-tausmin)*Model->dx[i];
      Qmax  = Model->Q[idx2(Nx,i,Nb)];
      tauemax = (tau0/Qmin)*(sqrtf(Qmax*Qmax+1.0)+1.0);
      tauemax = 1.0/tauemax;
      tausmax = (tau0/Qmin)*(sqrtf(Qmax*Qmax+1.0)-1.0);
      tausmax = 1.0/tausmax;

      // Interpolate relaxation times 
      tauey = tauemin + (tauemax-tauemin)*Model->dy[j];
      tausy = tausmin + (tausmax-tausmin)*Model->dy[j];

      // In the equations below the relaxation times taue and taus
      // are inverses (1/taue, 1/taus)
      // Compute alpha and eta coefficients
      argx = Model->dx[i];
      argy = Model->dy[j];
      // An extra tapering factor of exp(-(x/L)**2)
      // is used to taper some coefficeints 
      Model->Alpha1x[idx2(Nx,i,j)]   = expf(-argx)*expf(-Model->Dt*tausx);
      Model->Alpha1y[idx2(Nx,i,j)]   = expf(-argy)*expf(-Model->Dt*tausy);
      Model->Alpha2x[idx2(Nx,i,j)]   = Model->Dt*tauex;
      Model->Alpha2y[idx2(Nx,i,j)]   = Model->Dt*tauey;
      Model->Eta1x[idx2(Nx,i,j)]     = expf(-argx)*expf(-Model->Dt*tausx);
      Model->Eta1y[idx2(Nx,i,j)]     = expf(-argy)*expf(-Model->Dt*tausy);
      Model->Eta2x[idx2(Nx,i,j)]     = Model->Dt*tauex;
      Model->Eta2y[idx2(Nx,i,j)]     = Model->Dt*tauey;
 
      // Compute the change in moduli due to
      // visco-ealsticity (is equal to zero for the elastic case)
      Model->Dkappax[idx2(Nx,i,j)]   = Model->Kappa[idx2(Nx,i,j)]
                             *(1.0-tausx/tauex);
      Model->Dkappay[idx2(Nx,i,j)]   = Model->Kappa[idx2(Nx,i,j)]
                             *(1.0-tausy/tauey);
      Model->Drhox[idx2(Nx,i,j)]     = (Model->Rho[idx2(Nx,i,j)])
                             *(1.0-tausx/tauex);
      Model->Drhoy[idx2(Nx,i,j)]     = (Model->Rho[idx2(Nx,i,j)])
                             *(1.0-tausy/tauey);
    }
  }

  return(Model);
}
//
// Modelstability checks velocity model for stability.
// 
// Parameters:
//       
//     - Model : Model object
//
// Return      : Stability index
float ModelStability(struct model *Model)
{
  int nx,ny;
  int i,j;
  float vp,stab;

  nx = Model->Nx;
  ny = Model->Ny;
  for(j=0; j<ny; j=j+1){
    for(i=0; i<nx; i=i+1){
      vp = sqrtf(Model->Kappa[idx2(nx,i,j)]*Model->Rho[idx2(nx,i,j)]);
      stab = (vp*Model->Dt)/Model->Dx;
      if(stab > 1.0/sqrtf(2.0)){
        fprintf(stderr,"Stability index too large! %f \n", stab);
        fprintf(stderr,"vp: %f \n", vp);
        fprintf(stderr,"dt: %f \n", Model->Dt);
        fprintf(stderr,"dx: %f \n", Model->Dx);
      }
    }
  }

  return(stab);
}
// Modeld creates a 1D profile function tapering the left
// and right borders. 
// 
// Parameters:
//
//   d  : Input 1D float array
//   dx : Grid spacing
//   nb : Width of border zone   
//   
//   Return: OK if no error, ERR in all other cases.
int Modeld(float *d, float dx, int nb, int n){
  int i;

  for(i=0; i<n; i=i+1){
    d[i]=1.0;
  }

  // Taper left border
  for(i=0; i<nb;i=i+1){
      d[i] = d[i]*((i*dx)/(nb*dx))
                 *((i*dx)/(nb*dx));
  }

  // taper right border
  for(i=n-1-nb; i<n;i=i+1){
      d[i] = d[i]*((n-1-i)*dx)/(nb*dx)
                 *((n-1-i)*dx)/(nb*dx);
  }

  return(OK);
}
