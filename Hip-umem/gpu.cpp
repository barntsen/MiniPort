#include "hip/hip_runtime.h"

#include "gpu.h"
extern "C" {
#include <stdio.h>
}

// GpuNew -- Allocate unified memory on host and gpu
void *GpuNew(int n);
void GpuDelete(void *f);
void GpuErrorCheck(char *s);
void * GpuNew(int n){
  void *f;
  hipError_t cerr;
  cerr = hipMallocManaged(&f,n);
  if(cerr != hipSuccess){
    fprintf(stderr,"GpuNew:%s\n ", hipGetErrorString(cerr)) ;
    exit(1);
  }
  return(f);
}

//GpuDelete -- Delete unified memory on host and gpu
void GpuDelete(void *f){
  hipError_t cerr;

  cerr=hipFree(f);
  if(cerr != hipSuccess){
    fprintf(stderr,"GpuDelete:%s\n ", hipGetErrorString(cerr)) ;
    exit(1);
  }
  cerr=hipDeviceSynchronize();
  if(cerr != hipSuccess){
    fprintf(stderr,"GpuDelete:%s\n ", hipGetErrorString(cerr)) ;
    exit(1);
  }
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
