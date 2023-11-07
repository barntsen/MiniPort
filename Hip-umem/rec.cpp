// Rec object

// Imports
#include <stdlib.h>
#include <stdio.h> 
#include "gpu.h"
#include "rec.h"
#include "util.h"
// RecNew is the constructor for receiver objects.
//
// Arguments:
//   Model:  Model object
//   rx:     Integer array with position of receivers in the 
//           x-direction (gridpoints)
//   ry:     Integer array with position of receivers in the 
//           y-direction (gridpoints)
//   nt:     No of time samples in the receiver data
//   resamp: Resample factor relative to the modelling time sample interval
//   file:   File name for snap shots
//
//  Returns: Receiver object  
//----------------------------------------------------------------------------
struct rec *RecNew(int *rx, int *ry, int nt, 
                  int resamp, int sresamp, char *file, int nr)
{
  struct rec *Rec;
  int i,j;

  Rec = (struct rec *)GpuNew(sizeof(struct rec));
  Rec->nr = nr;
  Rec->rx = rx;
  Rec->ry = ry;
  Rec->nt = nt;
  Rec->p = (float*) GpuNew((Rec->nr*Rec->nt)*sizeof(float));
  for(j=0; j<Rec->nr; j=j+1)
  {
    for(i=0; i<Rec->nt; i=i+1){
      Rec->p[idx2(Rec->nt,i,j)]=0.0;
    }
  }
  Rec->resamp = resamp;
  Rec->sresamp = sresamp;
  Rec->pit = 0;
  if(Rec->sresamp > 0){
    Rec->fd = fopen(file,"w");
  }
  
  return(Rec);
}  
// RecReciver records data at the receiver
//
// Arguments: 
//  Rec:    : Receiver object
//  it      : Current time step
//  p:      : Pressure data at time step no it
//
// Returns  : Integer (OK or ERR)
//-----------------------------------------------------------------------------
int RecReceiver(struct rec *Rec,int it, float *p, int Nx, int Ny)
{
  int pos;
  int ixr,iyr;

  if(Rec->pit > Rec->nt-1){return(ERR);}

  if(it%Rec->resamp == 0){
    for (pos=0;pos<Rec->nr; pos=pos+1){  
      ixr=Rec->rx[pos];
      iyr=Rec->ry[pos];
      Rec->p[idx2(Rec->nt,Rec->pit,pos)] = p[idx2(Nx,ixr,iyr)];       
    } 
    Rec->pit = Rec->pit+1;
  }
  return(OK);
}
// Recsave stores receiver recording on file
//
// Arguments: 
//  Rec:    : Receiver object
//  file    : Output file name
//
// Returns  : Integer (OK or ERR)
//-----------------------------------------------------------------------------
int RecSave(struct rec *Rec, char *file)
{
  FILE *fd;
  int n;

  fd = fopen(file,"w");
  n = Rec->nr*Rec->nt;
  fwrite((void*)Rec->p,sizeof(float),n,fd);
  fclose(fd);

  return(OK);
}
// RecSnap records snapshots
//
// Arguments: 
//  Rec:    : Receiver object
//  it      : Current time step       
//  snp     : Pressure data
// Returns  : Integer (OK or ERR)
int RecSnap(struct rec *Rec,int it, float *snp, int Nx, int Ny)
{
  int n;
  char *tmp;
  int nb;
  
  if (Rec->sresamp <= 0){
    return(OK);
  }
  n = Nx*Ny;
  if(it%Rec->sresamp == 0){
    nb=fwrite((void*)snp,sizeof(float),n,Rec->fd);
  }
  return(OK);
}
