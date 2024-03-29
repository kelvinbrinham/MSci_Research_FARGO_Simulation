#include "fargo3d.h"

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
	field[l] = SIGMA0*pow((Ymed(j)/R0),-SIGMASLOPE)*(1.0+NOISE*(drand48()-.5));
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
	field[l] = ASPECTRATIO * pow(Ymed(j)/R0, FLARINGINDEX) * vk; //sqrt(G*MSTAR/Ymed(j))
#ifdef ADIABATIC
	field[l] = field[l]*field[l]*d[l]/(GAMMA-1.0);
#endif
      }
    }
  }    
}

void InitVazimPlanet() {

  int i,j,k;
  real *field;
  real dr, dz;
  real r, z, H, r0, rho_o, t;
  real rho;
  FILE *fo;
  real vt, omega;
  real *vr;
  real *cs;

  field = Vx->field_cpu;
  vr = Vy->field_cpu;
  cs = Energy->field_cpu;
    
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
	omega = sqrt(G*MSTAR/r/r/r);
	vt = omega*r*sqrt(1.0+pow(ASPECTRATIO,2)*pow(r/R0,2*FLARINGINDEX)*
			  (2.0*FLARINGINDEX - 1.0 - SIGMASLOPE));

	vt -= OMEGAFRAME*r;

	field[l] = vt*(1.+ASPECTRATIO*NOISE*(drand48()-.5));
	vr[l] = r*omega*ASPECTRATIO*NOISE*(drand48()-.5);
      }
    }
  }    
}

void CondInit() {
  
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
  InitVazimPlanet ();

#ifdef DUST
  real K0 = 1.0 / STOKES  ;    
  real rho0 = SIGMA0 ;

  int iDust ;
  for (iDust = 0; iDust < NDUST; iDust++) {
    for (k=0; k<Nz+2*NGHZ; k++) {
      for (j=0; j<Ny+2*NGHY; j++) {
	for (i=NGHX; i<Nx+NGHX; i++) {
	  Dust_Density[iDust]->field_cpu[l] = DUST_TO_GAS*Density->field_cpu[l];
	  Dust_Vx[iDust]->field_cpu[l] = Vx->field_cpu[l] ;
	  Dust_Vy[iDust]->field_cpu[l] = 0 ;
	}
      }
    }
#ifdef EPSTEIN_2D
    Drag_Coeff[iDust] = K0 * pow(10, iDust) / rho0 ;
#elif (defined STOKES_NUMBER)
    Drag_Coeff[iDust] = K0 * pow(10, iDust) ;
#else
#error Drag law not supported
#endif
  }
#endif 
}
