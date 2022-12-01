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

void SubStep_drag_y_cpu (real dt, int iDustSpecies, real DragCoeff) {
  
//<USER_DEFINED>
  INPUT(Pot);  

  INPUT(Dust_Density[iDustSpecies]) ;
#ifdef X
  INPUT(Dust_Vx[iDustSpecies]) ;
#endif
#ifdef Y
  INPUT(Dust_Vy[iDustSpecies]);
  OUTPUT(Dust_Vy[iDustSpecies]);
#endif
#ifdef Z
  INPUT(Dust_Vz[iDustSpecies]) ;
#endif

  INPUT(Density);
  INPUT(Energy) ;
#ifdef Y
  INPUT(Vy);
  INPUT(Vy_temp);
  OUTPUT(Vy_temp);
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
#ifdef Y
  real* vy_gas_old  = Vy->field_cpu ;
  real* vy_gas_new  = Vy_temp->field_cpu ;
  real* vy      = Dust_Vy[iDustSpecies]->field_cpu;
#endif
#ifdef Z
  real* vz      = Dust_Vz[iDustSpecies]->field_cpu;
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
  int llxp;
#endif
#ifdef Y
  int llym;
#endif
#ifdef Z
  int llzp;
#endif
#ifndef CARTESIAN
  real vphi;
#endif
#ifdef SHEARINGBOX
  real vm1, vm2;
#endif
#ifdef SPHERICAL
  real vzz;
#endif 
  real r;
  real cs;
  DragData gas;
  DragData dust;
//<\INTERNAL>


//<CONSTANT>
// real xmin(Nx+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real OMEGAFRAME(1);
// real OORTA(1);
// real VERTICALDAMPING(1);
// real GAMMA(1);
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
#ifdef Y
	ll = l;
#ifdef X
	llxp = lxp;
#endif
#ifdef Y
	llym = lym;
#endif
#ifdef Z
	llzp = lzp;
#endif
    
#ifdef CARTESIAN
	r = sqrt(XC*XC + ymin(j)*ymin(j) + ZC*ZC) ;
#else
	r = ymin(j) ;
#endif
#if (defined(ISOTHERMAL) || defined(POLYTROPIC))
	cs = 0.5 * (e[ll] + e[llym]) ;
#endif
#ifdef ADIABATIC
	cs = sqrt(0.5*GAMMA*(GAMMA-1)*(e[ll]/rho_g[ll] + e[llym]/rho_g[llym]));
#endif
	
	gas.density = 0.5*(rho_g[ll] + rho_g[llym]) ;
	gas.vel = vy_gas_new[ll] ;
	gas.accel = (vy_gas_new[ll] - vy_gas_old[ll]) / dt ;

	dust.Kdrag = drag_constant(DragCoeff, r, gas.density, cs);
	dust.density = 0.5*(rho_d[ll] + rho_d[llym]) ;
	dust.vel = vy[ll] ;

	dust.accel = 0;
	
#ifdef CARTESIAN
#ifdef SHEARINGBOX 
	vm1 = vx[ll]+vx[llxp]; 
	vm2 = vx[llym]+vx[llxp-pitch]; 
	dust.accel += (.5*OMEGAFRAME*(vm1+vm2) -
		       4.0*OORTA*OMEGAFRAME*ymin(j)); 
#endif 
#endif

#ifdef CYLINDRICAL
	vphi = .25*(vx[ll] + vx[llxp] + vx[llym] + vx[llxp-pitch]);
	vphi += ymin(j)*OMEGAFRAME;

	dust.accel += vphi*vphi/ymin(j) ;
#endif

#ifdef SPHERICAL
	vphi =  .25*(vx[ll] + vx[llxp] + vx[llym] + vx[llxp-pitch]);
	vphi += ymin(j)*sin(zmed(k))*OMEGAFRAME;

	vzz = .25*(vz[ll] + vz[llzp]  + vz[llym] + vz[llzp-pitch]);
	dust.accel += (vphi*vphi + vzz*vzz)/ymin(j) ;
#endif

#ifdef POTENTIAL
	dust.accel -= (pot[ll]-pot[llym]) /(ymed(j)-ymed(j-1));
#endif

	dust.vel += dust.accel * dt;

	drag_update(dt, &gas, &dust, 1) ;

	vy_gas_new[ll] = gas.vel ;
	vy[ll] = dust.vel ;
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
