//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
#include "drag_law.h"
#include "drag_update.h"
//<\INCLUDES>

void SubStep_drag_x_cpu (real dt, int iDustSpecies, real DragCoeff) {

//<USER_DEFINED>
  INPUT(Pot);
  
  INPUT(Dust_Density[iDustSpecies]);
#ifdef X
  INPUT(Dust_Vx[iDustSpecies]);
  OUTPUT(Dust_Vx[iDustSpecies]);
#endif

  INPUT(Density);
  INPUT(Energy);
#ifdef X
  INPUT(Vx);
  INPUT(Vx_temp);
  OUTPUT(Vx_temp);
#endif
//<\USER_DEFINED>

//<EXTERNAL>
  real* pot = Pot->field_cpu;
  real* rho_g = Density->field_cpu;
  real* e     = Energy->field_cpu;
  real* rho_d = Dust_Density[iDustSpecies]->field_cpu;
#ifdef X
  real* vx_gas_old  = Vx->field_cpu;
  real* vx_gas_new  = Vx_temp->field_cpu;
  real* vx      = Dust_Vx[iDustSpecies]->field_cpu ;
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
  int ll;
#ifdef X
  int llxm;
#endif
#ifdef Y
  int llyp;
#endif
#ifdef Z
  int llzp;
#endif 
  real r;
  real cs;
  DragData gas  ;
  DragData dust ;
//<\INTERNAL>

//<CONSTANT>
// real GAMMA(1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
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
#ifdef X
	ll = l;
	llxm = lxm;

#ifdef CARTESIAN
	r = sqrt(Xmin(i)*Xmin(i) + YC*YC + ZC*ZC);
#else
	r = ymed(j) ;
#endif
#if (defined(ISOTHERMAL) || defined(POLYTROPIC))
	cs = 0.5 * (e[ll] + e[llxm]) ;
#endif
#ifdef ADIABATIC
	cs = sqrt(0.5*GAMMA*(GAMMA-1)*(e[ll]/rho_g[ll] + e[llxm]/rho_g[llxm]));
#endif

	gas.density  = 0.5*(rho_g[ll] + rho_g[llxm]) ;
	gas.vel = vx_gas_new[ll] ;
	gas.accel = (vx_gas_new[ll] - vx_gas_old[ll]) / dt ;

	dust.Kdrag = drag_constant(DragCoeff, r, gas.density, cs);
	dust.density = 0.5*(rho_d[ll] + rho_d[llxm]) ;
	dust.vel = vx[ll] ;

	dust.accel = 0 ;	
#ifdef POTENTIAL
	dust.accel -= (pot[ll]-pot[llxm])/zone_size_x(j,k);
#endif
	dust.vel += dust.accel * dt ;

	drag_update(dt, &gas, &dust, 1) ;

	vx_gas_new[ll] = gas.vel ;
	vx[ll] = dust.vel ;
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
