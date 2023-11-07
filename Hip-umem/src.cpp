// Src object

// Imports
#include<stdlib.h>
#include<stdio.h>
#include "src.h"
#include "util.h"
#include "math.h"
#include "gpu.h"

// Sricker creates a Ricker wavelet
int Srcricker(float *src,    float t0, float f0, int nt, float dt);

// SrcNew creates a new source object
struct src *SrcNew(float *source, int nsrc, int *sx, int *sy, int ns)
{
  int i;
  struct src *Src;
  Src = (struct src*)GpuNew(sizeof(struct src));
  Src->Src = (float*)GpuNew(sizeof(int)*nsrc);
  for (i=0; i< nsrc; i=i+1){
    Src->Src[i] = source[i];
  }

  Src->Sx = sx;
  Src->Sy = sy;
  Src->Ns = ns;
  
  
  return(Src);
}

// SrcDel deletes a source object
int SrcDel(struct src *Src)
{
  free(Src); 
  return(OK);
}
// Ricker pulse
int Srcricker(float *source, float t0, float f0, int nt, float dt)
{
  float t;
  float w0;
  float arg;
  int i;

  for(i=0; i<nt; i=i+1){
    t = i*dt-t0;
    w0 = 2.0*3.14159*f0;
    arg = w0*t; 
    source[i] = (1.0-0.5*arg*arg)*exp(-0.25*arg*arg);
  }
  return(OK);
}
