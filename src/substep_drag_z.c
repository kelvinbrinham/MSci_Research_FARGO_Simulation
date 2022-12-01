// DATE: 21 - Nov - 2014
// AUTHOR: richard_booth

//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
#include "drag_law.h"
#include "drag_update.h"
//<\INCLUDES>

void SubStep_drag_z_cpu (real dt, int iDustSpecies, real DragCoeff) {

//<USER_DEFINED>
  INPUT(Pot);
  
  INPUT(Dust_Density[iDustSpecies]) ;
#ifdef X
  INPUT(Dust_Vx[iDustSpecies]) ;
#endif
#ifdef Z
  INPUT(Dust_Vz[iDustSpecies]) ;
  OUTPUT(Dust_Vz[iDustSpecies]) ;
#endif
  
  INPUT(Density);
  INPUT(Energy);
#ifdef Z
  INPUT(Vz);
  INPUT(Vz_temp);
  OUTPUT(Vz_temp);
#endif
//<\USER_DEFINED>

//<EXTERNAL>
  real* pot = Pot->field_cpu;
  real* rho_g = Density->field_cpu;
  real* e     = Energy->field_cpu;
  real* rho_d = Dust_Density[iDustSpecies]->field_cpu;
#ifdef X
  real* vx    = Dust_Vx[iDustSpecies]->field_cpu;
#endif
#ifdef Z
  real* vz      = Dust_Vz[iDustSpecies]->field_cpu;
  real* vz_gas_old = Vz->field_cpu ;
  real* vz_gas_new = Vz_temp->field_cpu ;
#endif
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
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
#ifdef Z
  int llzm;
#endif
  real vphi;
  real r ;
  real cs ;
  DragData gas  ;
  DragData dust ;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real OMEGAFRAME(1);
// real VERTICALDAMPING(1);
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
#endif
#ifdef X
	llxp = lxp;
#endif
#ifdef Z
	llzm = lzm;

#ifdef CARTESIAN
	r = sqrt(XC*XC + YC*YC + Zmin(k)*Zmin(k));
#else
	r = ymed(j) ;
#endif
#if (defined(ISOTHERMAL) || defined(POLYTROPIC))
	cs = 0.5 * (e[ll] + e[llzm]) ;
#endif
#ifdef ADIABATIC
	cs = sqrt(0.5*GAMMA*(GAMMA-1)*(e[ll]/rho_g[ll] + e[llzm]/rho_g[llzm]));
#endif
	
	gas.density = 0.5*(rho_g[ll] + rho_g[llzm]) ;
	gas.vel = vz_gas_new[ll] ;
	gas.accel = (vz_gas_new[ll] - vz_gas_old[ll]) / dt ;

	dust.Kdrag = drag_constant(DragCoeff, r, gas.density, cs);
	dust.density = 0.5*(rho_d[ll] + rho_d[llzm]) ;
	dust.vel = vz[ll] ;

	dust.accel = 0;

#ifdef CARTESIAN
#ifdef POTENTIAL
	dust.accel -= (pot[ll]-pot[llzm]) / (zmed(k)-zmed(k-1));
#endif
#endif

#ifdef CYLINDRICAL
#ifdef POTENTIAL
	dust.accel -= (pot[ll]-pot[llzm]) / (zmed(k)-zmed(k-1));
#endif
#endif

#ifdef SPHERICAL
	vphi = .25*(vx[ll] + vx[llxp] + vx[llzm] + vx[llxp-stride]);
	vphi += ymed(j)*sin(zmin(k))*OMEGAFRAME;

	dust.accel += vphi*vphi*cos(zmin(k))/(sin(zmin(k))*ymed(j));

#ifdef POTENTIAL
	dust.accel -= (pot[ll]-pot[llzm]) / (ymed(j)*(zmed(k)-zmed(k-1)));
#endif
#endif	
	dust.vel += dust.accel * dt;

	drag_update(dt, &gas, &dust, 1) ;

	vz_gas_new[ll] = gas.vel ;
	vz[ll] = dust.vel ;
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
