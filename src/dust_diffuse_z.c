// DATE: 5 - Feb - 2015
// AUTHOR: richard_booth

//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void Diffuse_dust_z_cpu (real dt, Field *DustDensity, Field *DustVz) {
  
//<USER_DEFINED>
  INPUT(Density);
  INPUT(Energy);
  INPUT(DustDensity) ;
  OUTPUT(DustDensity) ;
#ifdef Z
  INPUT(DustVz) ;
  OUTPUT(DustVz) ;
#endif
//<\USER_DEFINED>

//<EXTERNAL>
  real* rho_g = Density->field_cpu;
  real* rho_d = DustDensity->field_cpu;
  real* e     = Energy->field_cpu;
#ifdef Z
  real* vz      = DustVz->field_cpu;
#endif
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY-1;
  int size_z = Nz+2*NGHZ-1;
  real dx = Dx;
//<\EXTERNAL>

//<INTERNAL>
  int i; //Variables reserved
  int j; //for the topology
  int k; //of the kernels
#ifdef Z
  int ll;
  int llzm;
  real D;
  real vdiff;
  real cs2;
  real C_l;
  real C_r;
  real r;
#endif
//<\INTERNAL>


//<CONSTANT>
// real xmin(Nx+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real ALPHA(1);
// real NU(1);
// real SMALLDUSTDENSITY(1);
// real GAMMA(1);
// real SCHMIDT(1);
//<\CONSTANT>

//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for(k=1; k<size_z; k++) {
#endif
#ifdef Y
    for(j=1; j<size_y; j++) {
#endif
#ifdef X
      for(i=0; i<size_x; i++) {
#endif
//<#>
#ifdef Z
	ll = l;
	llzm = lzm;
#if (defined(ISOTHERMAL) || defined(POLYTROPIC))
	cs2 = 0.5 * (e[ll] + e[llzm]) ;
	cs2 = cs2*cs2;
#endif
#ifdef ADIABATIC
	cs2 = 0.5*GAMMA*(GAMMA-1)*(e[ll]/rho_g[ll] + e[llzm]/rho_g[llzm]);
#endif		   
#ifdef ALPHAVISCOSITY
#ifdef CARTESIAN
	r = sqrt(XC*XC + YC*YC + ZC*ZC) ;
#else
	r = ymed(j) ;
#endif	
	D = ALPHA*cs2*sqrt(r*r*r/(G*MSTAR)) / SCHMIDT;
#else
	D = NU / SCHMIDT ;
#endif
	
	C_l = rho_d[llzm] / rho_g[llzm] ;
	C_r = rho_d[ll] / rho_g[ll] ;

#ifndef SPHERICAL
	vdiff =  - D * (log(C_r) - log(C_l)) / (zmed(k) - zmed(k-1));
#else
	vdiff =  - D * (log(C_r) - log(C_l)) / (ymed(j)*(zmed(k) - zmed(k-1)));
#endif
	// Flux Limit
	vdiff /= (1 +vdiff*vdiff/cs2) ;

	vdiff += vz[ll] ;

	if (C_l <= SMALLDUSTDENSITY && !(vdiff < 0))
	  vdiff = 0 ;
	else if (C_r <= SMALLDUSTDENSITY && !(vdiff > 0))
	  vdiff = 0 ;

	vz[ll] = vdiff ;
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
