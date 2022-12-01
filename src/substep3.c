//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
#include "cool.h"
#include <assert.h>
//<\INCLUDES>

void SubStep3_cpu (real dt) {

//<USER_DEFINED>
  INPUT(Energy);
#ifdef COOLING
  INPUT(Density) ;
#endif
#ifdef X
  INPUT(Vx_temp);
#endif
#ifdef Y
  INPUT(Vy_temp);
#endif
#ifdef Z
  INPUT(Vz_temp);
#endif
  OUTPUT(Energy);
//<\USER_DEFINED>

//<EXTERNAL>
  real* e   = Energy->field_cpu;
#ifdef COOLING
  real* d   = Density->field_cpu;
#endif
#ifdef X
  real* vx  = Vx_temp->field_cpu;
#endif
#ifdef Y
  real* vy  = Vy_temp->field_cpu;
#endif
#ifdef Z
  real* vz  = Vz_temp->field_cpu;
#endif
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = XIP; 
  int size_y = Ny+2*NGHY-1;
  int size_z = Nz+2*NGHZ-1;
//<\EXTERNAL>

//<INTERNAL>
  int i; //Variables reserved
  int j; //for the topology
  int k; //of the kernels
  int ll;
#ifdef X
  int llxp;
#endif
#ifdef Y
  int llyp;
#endif
#ifdef Z
  int llzp;
#endif
  real term;
  real div_v;
#ifdef COOLING
  real r;
  real t_cool ;
#ifdef RADIATIVE_2D_COOLING
  real dens ;
  real T ;
  real kappa ;
  real tau ;
#endif
#endif
//<\INTERNAL>
  
//<CONSTANT>
// real xmin(Nx+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real Sxj(Ny+2*NGHY);
// real Syj(Ny+2*NGHY);
// real Szj(Ny+2*NGHY);
// real Sxk(Nz+2*NGHZ);
// real Syk(Nz+2*NGHZ);
// real Szk(Nz+2*NGHZ);
// real InvVj(Ny+2*NGHY);
// real BETA_COOLING_TIME(1);
// real BETA_COOLING_SLOPE(1);
// real SIGMACUTOFF(1);
// real GAMMA(1);
//<\CONSTANT>

//<MAIN_LOOP>
  
  i = j = k = 0;
  
#ifdef Z
  for(k=0; k<size_z; k++) {
#endif
#ifdef Y
    for(j=0; j<size_y; j++) {
#endif
#ifdef X
      for(i=0; i<size_x; i++) {
#endif
//<#>

	ll = l;
#ifdef X
	llxp = lxp;
#endif
#ifdef Y
	llyp = lyp;
#endif
#ifdef Z
	llzp = lzp;
#endif
	div_v = 0.0;
#ifdef X
	div_v += (vx[llxp]-vx[ll])*SurfX(j,k);
#endif
#ifdef Y
	div_v += (vy[llyp]*SurfY(j+1,k)-vy[ll]*SurfY(j,k));
#endif
#ifdef Z
	div_v += (vz[llzp]*SurfZ(j,k+1)-vz[ll]*SurfZ(j,k));
#endif
	term = 0.5 * dt * (GAMMA - 1.) * div_v * InvVol(j,k);
	e[ll] *= (1.0-term)/(1.0+term);

#if (defined ADIABATIC && defined COOLING)

#if   defined CARTESIAN
	r = sqrt(XC*XC + YC*YC + ZC*ZC) ;
#elif defined CYLINDRICAL
	r = sqrt(ymed(j)*ymed(j) + zmed(k)*zmed(k)) ;
#elif defined SPHERICAL
	r = Ymed(j) ;
#endif

#ifdef BETA_COOLING
	// Cooling time is a constant multiple of the orbital time
	t_cool = BETA_COOLING_TIME * sqrt(r*r*r/(G*MSTAR)) ;
	
	if (BETA_COOLING_SLOPE != 0)
	  t_cool *= pow(r/SIGMACUTOFF, BETA_COOLING_SLOPE) ;
#elif RADIATIVE_2D_COOLING
	T = (GAMMA-1) * (e[ll]/d[ll]) ;
	dens = d[ll] / sqrt(2*M_PI * GAMMA*T * r*r*r/(G*MSTAR)) ;	
	T /= R_MU 

	kappa = OPACITY(dens, T) ;
	tau = kappa * d[ll] ;
	t_cool = 
	  3*((GAMMA-1)/R_MU) * (sig*tau + 1 / kappa) / (16 * STEFANK * T*T*T) ;
#endif
	e[ll] /= 1 + dt / t_cool ;

#endif

//<\#>
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
//<\MAIN_LOOP>
}
