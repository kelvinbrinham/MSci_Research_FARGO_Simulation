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
  int begin_i = 0;
  int end_i = Nx;

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
  int begin_i = 0;
  int end_i = Nx;
  
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

real ConstructSequence (u, v, n)
     real *u, *v;
     int n;
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

  int i,j,k;
  real *field;
  real dr, dz;
  real r, z, H, r0, rho_o, t;
  real rho;
  FILE *fo;
  real vt, omega;
  real *vr;
  real *cs;
  real nu;

  field = Vx->field_cpu;
  vr = Vy->field_cpu;
  cs = Energy->field_cpu;
    
  boolean GhostInclude = FALSE;
  
  int begin_k =(GhostInclude ? 0 : NGHZ);
  int end_k = Nz+2*NGHZ-begin_k;
  int begin_j =(GhostInclude ? 0 : NGHY);
  int end_j = Ny+2*NGHY-begin_j;
  int begin_i = 0;
  int end_i = Nx;

  real* GlobalRmed = malloc ( (NY) * sizeof(real));
  real* Radii = malloc( (NY+1) * sizeof(real));
  real* vt_int = malloc(NY*sizeof(real));
  real* GLOBAL_SOUNDSPEED = malloc( (NY)*sizeof(real));
  real* Sigma = malloc( (NY) *sizeof(real));
  real* vt_cent = malloc( (NY+1) * sizeof(real));

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

  for (i=0; i<NY; i++) {
    GlobalRmed[i] = 0.5*(Radii[i+1] + Radii[i]);
    GLOBAL_SOUNDSPEED[i] = ASPECTRATIO * pow(GlobalRmed[i]/R0, FLARINGINDEX) * sqrt(G*MSTAR/GlobalRmed[i]);
    Sigma[i] = SIGMA0*pow(( GlobalRmed[i] /R0),-SIGMASLOPE);
  }

  vt_int[0]=0;
  for (i = 1; i < NY; i++) {
    vt_int[i]=(GLOBAL_SOUNDSPEED[i]*GLOBAL_SOUNDSPEED[i]*Sigma[i]- \
	       GLOBAL_SOUNDSPEED[i-1]*GLOBAL_SOUNDSPEED[i-1]*Sigma[i-1])/	\
      (.5*(Sigma[i]+Sigma[i-1]))/(GlobalRmed[i]-GlobalRmed[i-1])+	\
      G*(1.0/GlobalRmed[i-1]-1.0/GlobalRmed[i])/(GlobalRmed[i]-GlobalRmed[i-1]);
	  
    vt_int[i] = sqrt(vt_int[i]*Radii[i])-Radii[i]*OMEGAFRAME;
  }
  real t1 = vt_cent[0] = vt_int[1]+.75*(vt_int[1]-vt_int[2]);
  real r1 = ConstructSequence (vt_cent, vt_int, NY);
  vt_cent[0] += .25*(vt_int[1]-vt_int[2]);
  real t2 = vt_cent[0];
  real r2 = ConstructSequence (vt_cent, vt_int, NY);
  t1 = t1-r1/(r2-r1)*(t2-t1);
  vt_cent[0] = t1;
  ConstructSequence (vt_cent, vt_int, NY);
  vt_cent[NY]=vt_cent[NY-1];

  free(vt_int);
  free(GLOBAL_SOUNDSPEED);
  free(Sigma);
  free(GlobalRmed);
  free(Radii);



  for (k = begin_k; k<end_k; k++) {
    for (j = begin_j; j<end_j; j++) {
      for (i = begin_i; i<end_i; i++) {
	r = Ymin(j);
	omega = sqrt(G*MSTAR/r/r/r);

	vt = vt_cent[j+Y0-NGHY];

	field[l] = vt*(1.+ASPECTRATIO*NOISE*(drand48()-.5));
	if (ALPHA != 0) {
	  const real aspectratiosqd = pow(ASPECTRATIO,2.)*pow(r/R0,2*FLARINGINDEX);
	  nu = ALPHA * aspectratiosqd * r*r * omega;
	}
	else
	  nu = NU;
	vr[l] = -3. * (1 - SIGMASLOPE + 2 * FLARINGINDEX) * nu / r;
      }
    }
  } 

  free(vt_cent);
   
}

inline boolean IsPlanetCloserToStar (int iplanet, real radius) {
  char error_message[200];
  if (iplanet >= Sys->nb) {
    sprintf(error_message,"Attempting to access planet %i, but there are only %i in this simulation!",iplanet,Sys->nb);
    prs_error(error_message);
  }
  real x = Sys->x[iplanet];
  real y = Sys->y[iplanet];
  real radius_planet = x*x + y*y;
  if (radius_planet < radius)
    return TRUE;
  else
    return FALSE;
}

void SetZeroInsidePlanet(Field* field, int iplanet) {
  int i,j,k;
  real *array = field->field_cpu;
  real radius;  
  boolean GhostInclude = TRUE;
  
  int begin_k =(GhostInclude ? 0 : NGHZ);
  int end_k = Nz+2*NGHZ-begin_k;
  int begin_j =(GhostInclude ? 0 : NGHY);
  int end_j = Ny+2*NGHY-begin_j;
  int begin_i = 0;
  int end_i = Nx;

  for (k = begin_k; k<end_k; k++) {
    for (j = begin_j; j<end_j; j++) {
      for (i = begin_i; i<end_i; i++) {
        radius = Ymed(j);
	if (!IsPlanetCloserToStar(iplanet, radius))
	  array[l] = 0.;
      }
    }
  }

}

void SetZeroOutsidePlanet(Field* field, int iplanet) {
  int i,j,k;
  real *array = field->field_cpu;
  real radius;  
  boolean GhostInclude = TRUE;
  
  int begin_k =(GhostInclude ? 0 : NGHZ);
  int end_k = Nz+2*NGHZ-begin_k;
  int begin_j =(GhostInclude ? 0 : NGHY);
  int end_j = Ny+2*NGHY-begin_j;
  int begin_i = 0;
  int end_i = Nx;

  for (k = begin_k; k<end_k; k++) {
    for (j = begin_j; j<end_j; j++) {
      for (i = begin_i; i<end_i; i++) {
        radius = Ymed(j);
	if (IsPlanetCloserToStar(iplanet, radius))
	  array[l] = 0.;
      }
    }
  }

}

void InitDust()
{

  int i,j,k, iDust;
    
  boolean GhostInclude = TRUE;
  
  int begin_k =(GhostInclude ? 0 : NGHZ);
  int end_k = Nz+2*NGHZ-begin_k;
  int begin_j =(GhostInclude ? 0 : NGHY);
  int end_j = Ny+2*NGHY-begin_j;
  int begin_i = 0;
  int end_i = Nx;

  real STOKES0 = 1.0 ;

  for(i = 0; i < NDUST; ++i)
    Drag_Coeff[i] = (1. / STOKES0) * pow(10, i) ;

  double kDrag = (1 - 2*FLARINGINDEX + SIGMASLOPE) ;
  for (k = begin_k; k<end_k; k++) {
    for (j = begin_j; j<end_j; j++) {
      for (i = begin_i; i<end_i; i++) {

	for (iDust = 0; iDust < NDUST; iDust++)
	  {
	    Dust_Density[iDust]->field_cpu[l] =  0.1*Density->field_cpu[l] ;
	    Dust_Vx[iDust]->field_cpu[l] = Vx->field_cpu[l];
	    Dust_Vy[iDust]->field_cpu[l] = Vy->field_cpu[l] ;
	  }
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
  InitDust() ;

}
