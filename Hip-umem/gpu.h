//#define NTHREADS 1024
//#define NBLOCKS  1024

extern int NBLOCKS;
extern int NTHREADS;

void *GpuNew(int n);
void GpuDelete(void *f);
void GpuError();




