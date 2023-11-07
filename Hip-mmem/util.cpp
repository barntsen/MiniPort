
#include<time.h>
// Clock returns the time in nano seconds
float Clock()
{
  struct timespec tp;
  clock_gettime(CLOCK_MONOTONIC, &tp);
  return (float)((double)tp.tv_sec + (double)tp.tv_nsec*1.0e-9) ;
}
