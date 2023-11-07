// Differentiator exports

//diff object
struct diff {
  int  l;        // Differentiator length
  int lmax;      // Max differentiator length
  float *coeffs; 
  float *w;     // Differentiator weights 
};
//Methods:
struct diff *DiffNew(int l);
void DiffDxminus(struct diff *Diff, float *A, float *dA, float dx, int nx, int ny);
void DiffDyminus(struct diff *Diff, float *A, float *dA, float dx, int nx, int ny);
void DiffDxplus(struct diff *Diff, float *A, float *dA, float dx,  int nx, int ny);
void DiffDyplus(struct diff *Diff, float *A, float *dA, float dx,  int nx, int ny);
