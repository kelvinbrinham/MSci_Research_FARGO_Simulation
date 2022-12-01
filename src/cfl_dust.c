//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void cfl_dust_cpu(int iDustSpecies) { 

//<USER_DEFINED>
  OUTPUT(DensStar);
#ifdef X
  INPUT(Dust_Vx[iDustSpecies]);
  INPUT2D(VxMed);
#endif
#ifdef Y
  INPUT(Dust_Vy[iDustSpecies]);
#endif
#ifdef Z
  INPUT(Dust_Vz[iDustSpecies]);
#endif
  INPUT2DINT(Dust_Zmin[iDustSpecies]);
  INPUT2DINT(Dust_Zmax[iDustSpecies]);
//<\USER_DEFINED>

//<EXTERNAL>
  real* dtime = DensStar->field_cpu;
#ifdef X
  real* vx = Dust_Vx[iDustSpecies]->field_cpu;
  real* vxmed = VxMed->field_cpu;
#endif
#ifdef Y
  real* vy = Dust_Vy[iDustSpecies]->field_cpu;
#endif
#ifdef Z
  real* vz = Dust_Vz[iDustSpecies]->field_cpu;
#endif
  int* kmin = Dust_Zmin[iDustSpecies]->field_cpu;
  int* kmax = Dust_Zmax[iDustSpecies]->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+NGHY;
  int size_z = Nz+NGHZ;
  real dx = Dx;
  int pitch2d = Pitch2D;
  real nu = NU;
//<\EXTERNAL>

//<INTERNAL>
  real dtmin = 1e30;
  int i;
  int j;
  int k;
  int ll;
  int llxp;
  int llyp;
  int llzp;
  real cfl2=0.0;
  real cfl3=0.0;
  real cfl4=0.0;
  real vxx, vxxp;
  real soundspeed;
  real soundspeed2;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real GAMMA(1);
// real CFL(1);
//<\CONSTANT>

//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for (k=NGHZ; k<size_z; k++) {
#endif
#ifdef Y
    for (j=NGHY; j<size_y; j++) {
#endif
#ifdef X
      for (i=0; i<size_x; i++) {
#endif
//<#>
	// Skip cells that we have identified as boundaries
	if (k >= kmin[l2D_XY] && k < kmax[l2D_XY]) {
	  ll = l;
	  llxp = lxp;
	  llyp = lyp;
	  llzp = lzp;
#ifdef X
#ifdef STANDARD
	  vxx = vx[ll];
	  vxxp= vx[llxp];
#else
	  vxx = vx[ll] - vxmed[l2D];
	  vxxp= vx[llxp] - vxmed[l2D];
#endif
#endif

#ifdef X
	  cfl2 = (max2(fabs(vxx),fabs(vxxp)))/zone_size_x(j,k);
#endif
#ifdef Y
	  cfl3 = (max2(fabs(vy[ll]),fabs(vy[llyp])))/zone_size_y(j,k);
#endif
#ifdef Z
	  cfl4 = (max2(fabs(vz[ll]),fabs(vz[llzp])))/zone_size_z(j,k);
#endif


	  dtime[ll] = CFL/sqrt(cfl2*cfl2 + cfl3*cfl3 + cfl4*cfl4);

	  //	  if (dtime[ll] <= dtmin) {
	  //	    dtmin = dtime[ll];
	  //	    INSPECT_REAL (cfl5_a);
	  //	    INSPECT_REAL (cfl5_b);
	  //	    INSPECT_REAL (cfl5_c);
	  //	    INSPECT_INT (i);
	  //	    INSPECT_INT (j);
	  //	    INSPECT_INT (k);
	  //	  }
	}
//<\#>
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
//<\MAIN_LOOP>
#endif
//<LAST_BLOCK>
  cfl_b();
//<\LAST_BLOCK>

}
