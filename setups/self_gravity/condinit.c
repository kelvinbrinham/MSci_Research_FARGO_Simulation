#include "fargo3d.h"

double Sigma(double r) {
  return SIGMA0*pow(r/R0,-SIGMASLOPE)*exp(-r/SIGMACUTOFF) ; 
}

double cs(double r) {
  real cs0 = INITIALTOOMREQ * 
    M_PI * G * SIGMA0 * pow(SIGMACUTOFF, 1.5-SIGMASLOPE) * exp(-1) ;

  return cs0 * pow(r/SIGMACUTOFF, -0.25);
}

void InitDensPlanet() {

  int i,j,k;
  real *field;

  field = Density->field_cpu;
    
  boolean GhostInclude = TRUE;
  
  int begin_k =(GhostInclude ? 0 : NGHZ);
  int end_k = Nz+2*NGHZ-begin_k;
  int begin_j =(GhostInclude ? 0 : NGHY);
  int end_j = Ny+2*NGHY-begin_j;
  int begin_i = (GhostInclude ? 0 : NGHX);
  int end_i = Nx+2*NGHX-begin_i;

  for (k = begin_k; k<end_k; k++) {
    for (j = begin_j; j<end_j; j++) {
      for (i = begin_i; i<end_i; i++) {
	field[l] = 
	  Sigma(Ymed(j)) * (1 + NOISE*(drand48()-.5)) * (1 + 0.01*sin(i*6*M_PI/NX)) ;
      }
    }
  }
}

void InitSoundSpeedPlanet() {

  int i,j,k;
  real *field;
  real dr, dz;
  real r, z, H, r0, rho_o, t, omega, vk;
  real rho;
  FILE *fo;
  real *d;
  real *e;

  field = Energy->field_cpu;
  d = Density->field_cpu;

  boolean GhostInclude = TRUE;
  
  int begin_k =(GhostInclude ? 0 : NGHZ);
  int end_k = Nz+2*NGHZ-begin_k;
  int begin_j =(GhostInclude ? 0 : NGHY);
  int end_j = Ny+2*NGHY-begin_j;
  int begin_i = (GhostInclude ? 0 : NGHX);
  int end_i = Nx+2*NGHX-begin_i;

  
  for (k = begin_k; k<end_k; k++) {
    for (j = begin_j; j<end_j; j++) {
      for (i = begin_i; i<end_i; i++) {	
	r = Ymed(j);
	vk = sqrt(G*MSTAR/r);
	field[l] = cs(r) ;
#ifdef ADIABATIC
	field[l] = field[l]*field[l]*d[l]/(GAMMA*(GAMMA-1.0));
#endif
      }
    }
  }    
}

real ConstructSequence (real *u, real *v, int n)
{
  int i;
  real lapl=0.0;
  for (i = 1; i < n; i++)
    u[i] = 2.0*v[i]-u[i-1];
  for (i = 1; i < n-1; i++) {
    lapl += fabs(u[i+1]+u[i-1]-2.0*u[i]);
  }
  return lapl;
}



void InitVazimPlanet() {

  INPUT(Pot) ;

  int i,j,k;
  real *field;
  real dr, dz;
  real r, z, H, r0, rho_o, t;
  real rho;
  FILE *fo;
  real vt, vk;
  real *d;
  real *vr;
  real *e;
  real *phi;


  field = Vx->field_cpu;
  vr = Vy->field_cpu;
  e = Energy->field_cpu;
  phi = Pot->field_cpu;
  d   = Density->field_cpu;

  int pitch = Pitch_cpu;

  boolean GhostInclude = FALSE;
  
  int begin_k =(GhostInclude ? 0 : NGHZ);
  int end_k = Nz+2*NGHZ-begin_k;
  int begin_j =(GhostInclude ? 0 : NGHY);
  int end_j = Ny+2*NGHY-begin_j;
  int begin_i = (GhostInclude ? 0 : NGHX);
  int end_i = Nx+2*NGHX-begin_i;


  // Construct vphi to balance forces
  real* GlobalRmed = malloc( (NY)   * sizeof(real));
  real* Radii      = malloc( (NY+1) * sizeof(real));
  real* vt_int     = malloc( (NY)   * sizeof(real));
  real* GlobalP    = malloc( (NY)   * sizeof(real));
  real* GlobalS    = malloc( (NY)   * sizeof(real));
  real* GlobalPot  = malloc( (NY)   * sizeof(real));
  real* vt_cent    = malloc( (NY+1) * sizeof(real));
  
  
  // Read in the global grid
  FILE* radii_file;
  char radii_filename[200];
  sprintf(radii_filename,"%sdomain_y.dat",OUTPUTDIR);
  radii_file = fopen(radii_filename,"r");
  for (i=0; i<NGHY; i++) {
    real temp;
    fscanf(radii_file,"%lg",&temp);
  }
  for (i=0; i< NY+1; i++) {
    fscanf(radii_file,"%lg",Radii+i);
  }
  fclose(radii_file);

  real cs0 = INITIALTOOMREQ * 
    M_PI * G * SIGMA0 * pow(SIGMACUTOFF, 1.5-SIGMASLOPE) * exp(-1) ;

  // Compute global arrays 
  for (i=0; i<NY; i++) {
    GlobalRmed[i] = 0.5*(Radii[i+1] + Radii[i]);
    if (i > 0)
      GlobalS[i] = 0.5*(Sigma(GlobalRmed[i]) + Sigma(GlobalRmed[i-1])) ;
    GlobalP[i] = Sigma(GlobalRmed[i])* cs(GlobalRmed[i])*cs(GlobalRmed[i]) ;
#ifdef ADIABATIC
    GlobalP[i] /= GAMMA ;
#endif
#ifdef SELF_GRAVITY
    GlobalPot[i] = -G*MSTAR/GlobalRmed[i] ;
#endif
  }

#ifdef SELF_GRAVITY
  // Gather the total potential
  real *localpot = malloc( NY * sizeof(real)) ;
  int  *offsets = malloc( CPU_Number * sizeof(int)) ;
  int  *sizes   = malloc( CPU_Number * sizeof(int)) ;

  // Get the size of the data on each processor:
  MPI_Allgather(&Ny, 1, MPI_INT, sizes, 1, MPI_INT, MPI_COMM_WORLD);
  // Compute the offsets
  offsets[0] = 0 ;
  //sizes[0] *= sizeof(real) ;
  for (i=1; i < CPU_Number; i++) {
    offsets[i] = offsets[i-1] + sizes[i-1] ;
    //sizes[i] *= sizeof(real) ;
  }

  // Fill the local potential array (use azimuthal average)
  i = k = 0;
  for (j = 0; j < Ny; j++) {
    localpot[j] = 0;
    for (i = 0;  i < Nx; i++) 
      localpot[j] += phi[l]/Nx ;
  }
  // Get the global potential array
  MPI_Allgatherv(localpot, sizes[CPU_Rank], MPI_DOUBLE,
		 GlobalPot, sizes, offsets, MPI_DOUBLE, MPI_COMM_WORLD);

  free(sizes) ;
  free(offsets) ;
  free(localpot);
#endif


  vt_int[0]=0;
  for (i = 1; i < NY; i++) {
    vt_int[i]= 
      ((GlobalP[i] - GlobalP[i-1])/GlobalS[i] + GlobalPot[i] - GlobalPot[i-1]) /
      (GlobalRmed[i]-GlobalRmed[i-1]);
	  
    vt_int[i] = sqrt(vt_int[i]*Radii[i])-Radii[i]*OMEGAFRAME;
  }

  real t1 = vt_cent[0] = vt_int[1]+.75*(vt_int[1]-vt_int[2]);
  real r1 = ConstructSequence (vt_cent, vt_int, NY);
  vt_cent[0] += .25*(vt_int[1]-vt_int[2]);
  real t2 = vt_cent[0];
  real r2 = ConstructSequence (vt_cent, vt_int, NY);
  
  vt_cent[0] = t1-r1/(r2-r1)*(t2-t1);
  ConstructSequence (vt_cent, vt_int, NY);

  vt_cent[NY]=vt_cent[NY-1];

  free(GlobalPot);
  free(GlobalS);
  free(GlobalP);
  free(vt_int);
  free(Radii);
  free(GlobalRmed);


  for (k = begin_k; k<end_k; k++) {
    for (j = begin_j; j<end_j; j++) {
      for (i = begin_i; i<end_i; i++) {
	real r = Ymin(j);
	real omega = sqrt(G*MSTAR/r/r/r);

	real vt = vt_cent[j+Y0-NGHY];

	field[l] = vt*(1.+ASPECTRATIO*NOISE*(drand48()-.5));
	vr[l] = vt*ASPECTRATIO*NOISE*(drand48()-.5);
      }
    } 
  }
  free(vt_cent);
  
}

void CondInit() {
  
  INPUT(Pot);
  OUTPUT(Pot) ;
  OUTPUT(Density);
  OUTPUT(Energy);
  OUTPUT(Vx);
  OUTPUT(Vy);

  int i,j,k;
  int index;
  real vt;
  real *field;
  real *rho;
  real *v1;
  real *v2;
  real *e;

#ifdef PLANETS
  Sys = InitPlanetarySystem(PLANETCONFIG);
  ListPlanets();
  if(COROTATING)
    OMEGAFRAME = GetPsysInfo(FREQUENCY);
  else
#endif
    OMEGAFRAME = OMEGAFRAME;

  
  InitDensPlanet ();
  InitSoundSpeedPlanet ();

#ifdef SELF_GRAVITY
  // Compute the potential so that we can use it to compute vphi
  FARGO_SAFE(Potential()); 
  FARGO_SAFE(SelfGravityPotential());
#endif

  InitVazimPlanet ();
}
