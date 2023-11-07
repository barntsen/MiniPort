// rec is the receiver object

struct rec 
{
    int nr; // No of receivers
    int *rx;       // Receiver x-postions
    int *ry;       // Receiver y-postions 
    FILE *fd;      // Snapshot output file descriptor
    int nt;        // No of time samples
    float *p;      // Pressure p[i,j] time sample no j at
                   // position no i 
    int Nx;        // Size of p 1st dim
    int Ny;        // Size of p 2nd dim
    float *wrk;    // Work array
    int   resamp;  // Resample factor for receivers
    int   sresamp; // Resample factor for snapshots
    int recon;     // =OK Record data =OFF No data recorded
    int snopn;     // =OK Record snapshots =OFF Do not record snapshots
    int pit;       // Next time sample to be recorded
};

// RecNew is the constructor for receivers.
// The return value is a Rec object
struct rec *RecNew(int * rx, int * ry, int nt, int resamp, int sresamp,
                  char* file, int nr);

// RecReciver records data at the receivers
int RecReceiver(struct rec *Rec,int it, float *p, int Nx, int Ny);

// Recsave stores receiver recording on file
int RecSave(struct rec *Rec , char *file);

// RecSnap records snapshots
int RecSnap(struct rec *Rec,int it, float * snp, int Nx, int Ny);

