#include "gpu.h"
extern "C" {
#include <stdio.h>
}
#include "hip/hip_runtime.h"

// GpuNew -- Allocate unified memory on host and gpu
void *GpuNew(int n);
void GpuDelete(void *f);
void GpuErrorCheck(char *s);

void * GpuNew(int n){
  void * mem;
  
  mem = malloc(n*sizeof(float));
  
  return(mem);
}

//GpuDelete -- Delete unified memory on host and gpu
void GpuDelete(void *f){
  free(f);
}

//GpuError -- Check for gpu errors
void GpuError(){
  hipDeviceSynchronize();
  hipError_t cerr;
  cerr = hipGetLastError();
  if(cerr != hipSuccess){
    fprintf(stderr,"%s\n",hipGetErrorString(cerr));
    exit(1);
  }
}
