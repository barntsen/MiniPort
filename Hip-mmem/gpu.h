//#define NTHREADS 1024
//#define NBLOCKS  256000
extern int NTHREADS;
extern int NBLOCKS;

void *GpuNew(int n);
void GpuDelete(void *f);
void GpuError();




