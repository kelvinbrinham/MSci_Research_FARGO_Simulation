#include "fargo3d.h"

void CondInit() {
  int i,j,k;
  real* rho = Density->field_cpu;
  real* e = Energy->field_cpu;
  int index;
#ifdef Z
  real* v1 = Vz->field_cpu;
  int dim1 = Nz + 2*NGHZ;
#define Q1 (Zmed(i) - ZMIN)/(ZMAX - ZMIN)
#endif
#ifdef Y
  real* v1 = Vy->field_cpu;
  int dim1 = Ny + 2*NGHY;
#define Q1 (Ymed(i) - YMIN)/(YMAX - YMIN)
#endif
#ifdef X
  real* v1 = Vx->field_cpu;
  int dim1 = Nx;
#define Q1 (Xmed(i) - XMIN)/(XMAX - XMIN)
#endif
  
  for (i = 0; i<dim1; i++) {
    e[i]   = 1.0/(GAMMA-1.0);
    rho[i] = 1.0;
    v1[i]  = 0.0;
    if (Q1 > 0.5) {
      rho[i] = 0.125;
#ifdef ADIABATIC
      e[i] = 0.1/(GAMMA-1.0);
#endif
    }
  }

#ifdef DUST
  real K0 = 0.1;
  int iDust ;
  for (iDust = 0; iDust < NDUST; iDust++) {
    for (i = 0; i<dim1; i++) {
      Dust_Vz[iDust]->field_cpu[i] = 0 ;
      Dust_Density[iDust]->field_cpu[i] = rho[i] * DUST_TO_GAS;
    }
    Drag_Coeff[iDust] = K0 * pow(10, iDust) ;
  }
  
#endif
  
}
