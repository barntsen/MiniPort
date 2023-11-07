// Ac2D is the acoustic modeling object

struct ac2d {
  float * p;     // Stress 
  float * vx;     // x-component of particle velocity
  float * vy;     // y-component of particle velocity
  float * exx;    // time derivative of strain x-component
  float * eyy;    // time derivative of strain y-component
  float * gammax; // Gamma functions in x-direction
  float * gammay; // Gamma function in y-direction
  float * thetax; // Theta function in x-direction
  float * thetay; // Theta function in y-direction
  float * dp;     // Stress 
  float * dvx;    // Device x-component of particle velocity
  float * dvy;    // Device y-component of particle velocity
  float * dexx;   // Device time derivative of strain x-component
  float * deyy;   // Device time derivative of strain y-component
  float * dgammax;// Device gamma in x-direction
  float * dgammay;// Device gamma in y-direction
  float * dthetax;// Device Theta in x-direction
  float * dthetay;// Device Theta in y-direction
  int ts;         // Timestep no
};

// Ac2dNew creates a new Ac2d object
  struct ac2d *Ac2dNew(struct model *Model);

// Ac2dSolve computes the pressure at the next timestep
   int Ac2dSolve(struct ac2d *Ac2d, struct model *Model,struct src *Src, 
                 struct rec *Rec,int nt, int l);
