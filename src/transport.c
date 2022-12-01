#include "fargo3d.h"

// CHANGES:
// AUTHOR: richard_booth
// DATE: 22 - Nov - 2014
//  - Modified routines to allow density / velocity to be specified
//  - Added a parameter specifying whether to the fluid is a gas, which controls
//    whether the internal energy is transported

static void (*__VanLeerX) (Field *, Field *, Field *, real);

void VanLeerX(Field *Dens, Field *DensStar, Field *Vx_t, real dt) {
  FARGO_SAFE(VanLeerX_a(Dens));
  FARGO_SAFE(VanLeerX_b(dt, Dens, DensStar, Vx_t));
}

void TransportX(Field *Q, Field *Qs, Field* Dens, Field *Vx_t, real dt) { 
  FARGO_SAFE(DivideByRho(Q,Dens));
  __VanLeerX(DivRho, Qs, Vx_t, dt);
  FARGO_SAFE(UpdateX (dt, Q, Qs, Vx_t));
}

void TransportY(Field *Q, Field *Qs, Field* Dens, Field *Vy_t, real dt) {
  FARGO_SAFE(DivideByRho(Q, Dens));
  FARGO_SAFE(VanLeerY_a(DivRho));
  FARGO_SAFE(VanLeerY_b(dt, DivRho, Qs, Vy_t));

  FARGO_SAFE(UpdateY (dt, Q, Qs, Vy_t));
}

void TransportZ(Field *Q, Field *Qs, Field* Dens, Field *Vz_t, real dt) {
  FARGO_SAFE(DivideByRho(Q, Dens));
  FARGO_SAFE(VanLeerZ_a(DivRho));
  FARGO_SAFE(VanLeerZ_b(dt, DivRho, Qs, Vz_t));
  FARGO_SAFE(UpdateZ (dt, Q, Qs, Vz_t));
}

void X_advection (Field *Dens, Field *Vx_t, real dt, boolean isGas) {
  int i;
#ifdef X
  __VanLeerX(Dens, DensStar, Vx_t, dt);
  TransportX(Mpx, Qs, Dens, Vx_t, dt);
  TransportX(Mmx, Qs, Dens, Vx_t, dt);
#endif
#ifdef Y
  TransportX(Mpy, Qs, Dens, Vx_t, dt);
  TransportX(Mmy, Qs, Dens, Vx_t, dt);
#endif
#ifdef Z
  TransportX(Mpz, Qs, Dens, Vx_t, dt);
  TransportX(Mmz, Qs, Dens, Vx_t, dt);
#endif
  if (isGas)
    {
#ifdef ADIABATIC
      TransportX(Energy, Qs, Dens, Vx_t, dt);
#endif
    }
  TransportX(Dens, Qs, Dens, Vx_t, dt);
}

void transport(real dt, Field* Dens, 
	       Field *Vx_t, Field *Vy_t, Field *Vz_t,
	       Field *Vx_new, Field *Vy_new, Field *Vz_new,
	       boolean isGas){
  int i;
#ifdef X
  FARGO_SAFE(momenta_x(Dens, Vx_t));
#endif
#ifdef Y
  FARGO_SAFE(momenta_y(Dens, Vy_t));
#endif
#ifdef Z
  FARGO_SAFE(momenta_z(Dens, Vz_t));
#endif

#ifdef DUST_DIFFUSE
  if (!isGas)
    {
#ifdef X
      FARGO_SAFE(Diffuse_dust_x(dt, Dens, Vx_t)) ;
#endif
#ifdef Y
      FARGO_SAFE(Diffuse_dust_y(dt, Dens, Vy_t)) ;
#endif
#ifdef Z
      FARGO_SAFE(Diffuse_dust_z(dt, Dens, Vz_t)) ;
#endif
    }
#endif


#ifdef X
#ifndef STANDARD
    FARGO_SAFE(ComputeVmed(Vx_t)); 
#endif
#endif


#ifdef Z
  FARGO_SAFE(VanLeerZ_a(Dens));
  FARGO_SAFE(VanLeerZ_b(dt, Dens, DensStar, Vz_t));
#ifdef X
  TransportZ(Mpx, Qs, Dens, Vz_t, dt);
  TransportZ(Mmx, Qs, Dens, Vz_t, dt);
#endif
#ifdef Y
  TransportZ(Mpy, Qs, Dens, Vz_t, dt);
  TransportZ(Mmy, Qs, Dens, Vz_t, dt);
#endif
#ifdef Z
  TransportZ(Mpz, Qs, Dens, Vz_t, dt);
  TransportZ(Mmz, Qs, Dens, Vz_t, dt);
#endif
  if (isGas)
    {
#ifdef ADIABATIC
      TransportZ(Energy, Qs, Dens, Vz_t, dt);
#endif
    }
  TransportZ(Dens, Qs, Dens, Vz_t, dt);
#endif

#ifdef Y
  FARGO_SAFE(VanLeerY_a(Dens));
  FARGO_SAFE(VanLeerY_b(dt, Dens, DensStar, Vy_t));

#ifdef X  
  TransportY(Mpx, Qs, Dens, Vy_t, dt);
  TransportY(Mmx, Qs, Dens, Vy_t, dt);
#endif
#ifdef Y
  TransportY(Mpy, Qs, Dens, Vy_t, dt);
  TransportY(Mmy, Qs, Dens, Vy_t, dt);
#endif
#ifdef Z
  TransportY(Mpz, Qs, Dens, Vy_t, dt);
  TransportY(Mmz, Qs, Dens, Vy_t, dt);
#endif

  if (isGas)
    {
#ifdef ADIABATIC
      TransportY(Energy, Qs, Dens, Vy_t, dt);
#endif
    }
  TransportY(Dens, Qs, Dens, Vy_t, dt);
#endif

#ifdef X
#ifdef STANDARD
  __VanLeerX = VanLeerX;
  X_advection (Dens, Vx_t, dt, isGas);
#else // FARGO algorithm below
  FARGO_SAFE(ComputeResidual(dt, Vx_t, Vx_new));
  __VanLeerX = VanLeerX;
  X_advection (Dens, Vx_new, dt, isGas); // Vx => variable residual
  //__VanLeerX= VanLeerX;
  __VanLeerX= VanLeerX_PPA;
  X_advection (Dens, Vx_t, dt, isGas); // Vx_temp => fixed residual @ given r. This one only is done with PPA
  __VanLeerX = VanLeerX;
  AdvectSHIFT(Mpx, Nshift);
  AdvectSHIFT(Mmx, Nshift);
#ifdef Y
  AdvectSHIFT(Mpy, Nshift);
  AdvectSHIFT(Mmy, Nshift);
#endif
#ifdef Z
  AdvectSHIFT(Mpz, Nshift);
  AdvectSHIFT(Mmz, Nshift);
#endif
  if (isGas)
    {
#ifdef ADIABATIC
      AdvectSHIFT(Energy, Nshift);
#endif
    }
  AdvectSHIFT(Dens, Nshift);
#endif
#endif

  FARGO_SAFE(NewVelocity_x(Dens, Vx_new));
  FARGO_SAFE(NewVelocity_y(Dens, Vy_new));
  FARGO_SAFE(NewVelocity_z(Dens, Vz_new));
}
