#include "fargo3d.h"

#define MIN_DUST_DENSITY    (1e-10)
#define DUST_DENSITY_FLOOR  (1e-20)

void boundary_dust_zmin_internal(int iDustSpecies, real DragCoeff) ;
void set_zmin_boundary_condition(int iDustSpecies, real DragCoeff) ;

void set_zmin_boundary_location(int iDustSpecies) ;
void set_zmax_boundary_location(int iDustSpecies) ;
void smooth_zmin_boundary_location(int iDustSpecies) ;
void smooth_zmax_boundary_location(int iDustSpecies) ;


void CondInit() {
  int i,j,k;
  real *v1;
  real *v2;
  real *v3;
  real *e;
  real *rho;
  real h;
  
  real omega;
  real r, r3;

  rho = Density->field_cpu;
  e   = Energy->field_cpu;
  v1  = Vx->field_cpu;
  v2  = Vy->field_cpu;
  v3  = Vz->field_cpu;

  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      h = ASPECTRATIO*Ymed(j);
      r = Ymed(j);
      r3 = r*r*r;
      omega = sqrt(G*MSTAR/(r3));
      for (i=NGHX; i<Nx+NGHX; i++) {
	v2[l] = v3[l] = 0.0;
	v1[l] = omega*r;
#ifdef CYLINDRICAL
	rho[l] = SIGMA0*pow(r/R0,-SIGMASLOPE)*exp(-pow(Zmed(k)/h,2.0)/2.0)/(ZMAX-ZMIN);
#else
	real xi = SIGMASLOPE+1.+FLARINGINDEX;
	real beta = 1.-2*FLARINGINDEX;
	real h = ASPECTRATIO*pow(r/R0,FLARINGINDEX);
	if (FLARINGINDEX == 0.0) {
	  rho[l] = SIGMA0/sqrt(2.0*M_PI)/(R0*ASPECTRATIO)*pow(r/R0,-xi)* \
	    pow(sin(Zmed(k)),-beta-xi+1./(h*h));
	} else {
	  rho[l] = SIGMA0/sqrt(2.0*M_PI)/(R0*ASPECTRATIO)*pow(r/R0,-xi)* \
	    pow(sin(Zmed(k)),-xi-beta)*					\
	    exp((1.-pow(sin(Zmed(k)),-2.*FLARINGINDEX))/2./FLARINGINDEX/(h*h));
	}
	v1[l] *= sqrt(pow(sin(Zmed(k)),-2.*FLARINGINDEX)-(beta+xi)*h*h);
	v1[l] -= OMEGAFRAME*r*sin(Zmed(k));
#endif
#ifdef ISOTHERMAL
	e[l] = h*sqrt(G*MSTAR/r);
#else
	e[l] = rho[l]*h*h*G*MSTAR/r/(GAMMA-1.0);
#endif
      }
    }
  }

#ifdef DUST
  real K0 = 100.;
  real rho0 = SIGMA0/(sqrt(2.0*M_PI)*(R0*ASPECTRATIO)) ;
  real cs0 = ASPECTRATIO*R0 / sqrt(G*MSTAR/(R0*R0*R0)) ;

  int iDust ;
  for (iDust = 0; iDust < NDUST; iDust++) {
    for (k=0; k<Nz+2*NGHZ; k++) {
      for (j=0; j<Ny+2*NGHY; j++) {
	for (i=NGHX; i<Nx+NGHX; i++) {
	  Dust_Density[iDust]->field_cpu[l] = DUST_TO_GAS * rho[l] ;
	  Dust_Vx[iDust]->field_cpu[l] = v1[l] ;
	  Dust_Vy[iDust]->field_cpu[l] = v2[l] ;
	  Dust_Vz[iDust]->field_cpu[l] = v3[l] ;
	}
      }
    }
    Drag_Coeff[iDust] = K0 * pow(10, iDust) / (rho0 * cs0);   
  }

  boundary_zmin_dust = boundary_dust_zmin_internal ;
#endif  
}

void boundary_dust_zmin_internal(int iDustSpecies, real DragCoeff) {
  // Find the new Z-location of the boundary and smooth twice
  set_zmin_boundary_location(iDustSpecies) ; 
  smooth_zmin_boundary_location(iDustSpecies);
  smooth_zmin_boundary_location(iDustSpecies);
  // Zero-out the region outside the new domain
  set_zmin_boundary_condition(iDustSpecies,DragCoeff) ;
}


void set_zmin_boundary_condition(int iDustSpecies, real DragCoeff) {

//<USER_DEFINED>
  INPUT(Dust_Density[iDustSpecies]);
  INPUT(Dust_Vx[iDustSpecies]);
  INPUT(Dust_Vy[iDustSpecies]);
  INPUT(Dust_Vz[iDustSpecies]);
  OUTPUT(Dust_Density[iDustSpecies]);
  OUTPUT(Dust_Vx[iDustSpecies]);
  OUTPUT(Dust_Vy[iDustSpecies]);
  OUTPUT(Dust_Vz[iDustSpecies]);

  INPUT2DINT(Dust_Zmin[iDustSpecies]);

  INPUT(Vx) ;
//<\USER_DEFINED>

//<INTERNAL>
  int i;
  int j;
  int k;
  int jact;
  int jgh;
  int kact;
  int kgh;
  int lgh;
  int lghs;
  int lact;
  int lacts;
  int lacts_null;
  real rho_mid;
  int k0;
  int kb ;
  int ll;
//<\INTERNAL>

//<EXTERNAL>
  real* gas_vx  = Vx->field_cpu;
  real* dust_density = Dust_Density[iDustSpecies]->field_cpu;
  real* dust_vx = Dust_Vx[iDustSpecies]->field_cpu;
  real* dust_vy = Dust_Vy[iDustSpecies]->field_cpu;
  real* dust_vz = Dust_Vz[iDustSpecies]->field_cpu;
  int* kmin = Dust_Zmin[iDustSpecies]->field_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = NGHZ;
  int nx = Nx;
  int ny = Ny;
  int nz = Nz;
  int nghy = NGHY;
  int nghz = NGHZ;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  real dx = Dx;
  real r0 = R0;
  real aspectratio = ASPECTRATIO;
  real flaringindex = FLARINGINDEX;
  real sigmaslope = SIGMASLOPE;
  real omegaframe = OMEGAFRAME;
//<\EXTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
//<\CONSTANT>

  //<MAIN_LOOP>

  i = j = k = 0;
#ifdef Y
    for(j=0; j<size_y; j++) {
#endif
#ifdef X
      for(i=0; i<size_x; i++) {
#endif
	// Find the vertical boundary:	
	k = nghz ;
	k0 = kmin[l2D_XY];

	// Reset the values outside the boundary to make life easy:
	for (k = 0 ; k < k0; k++) {
	  ll = l ;
	  dust_density[ll] = DUST_DENSITY_FLOOR;
	  dust_vz[ll] = 0;
	  dust_vy[ll] = 0;
	  dust_vx[ll] = gas_vx[ll] ;
	}
	// Extra edge for the z-velocity
	k = k0 ;
	dust_vz[l] = 0;

#ifdef X
      }
#endif
#ifdef Y
    }
#endif
//<\MAIN_LOOP>

}


void set_zmin_boundary_location(int iDustSpecies) {


//<USER_DEFINED>
  INPUT(Dust_Density[iDustSpecies]);
  OUTPUT2DINT(Dust_Zmin[iDustSpecies]);

  INPUT(Density);
//<\USER_DEFINED>

//<INTERNAL>
  int i;
  int j;
  int k;
  real rho_mid;
//<\INTERNAL>

//<EXTERNAL>
  real* density = Density->field_cpu;
  real* dust_density = Dust_Density[iDustSpecies]->field_cpu;
  int* kmin = Dust_Zmin[iDustSpecies]->field_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = NGHZ;
  int nx = Nx;
  int ny = Ny;
  int nz = Nz;
  int nghy = NGHY;
  int nghz = NGHZ;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
//<\EXTERNAL>

//<CONSTANT>
//<\CONSTANT>

//<MAIN_LOOP>

  i = j = k = 0;
#ifdef Y
    for(j=0; j<size_y; j++) {
#endif
#ifdef X
      for(i=0; i<size_x; i++) {
#endif
	// Find the vertical boundary:	
	k = nz + nghz;
	rho_mid = density[l] ;
	for (; k >  nghz; k--) 
	  if (dust_density[l] < MIN_DUST_DENSITY*rho_mid*DUST_TO_GAS)
	    break ;
	// Add a buffer a buffer zone of one cell.
	//   This should help the boundary move outwards when the density
	//   is increasing by collecting mass in the cell. Consider increasing
	//   the number of buffer cells if mass diffusing across the boundary
	//   is preventing expansion.
	kmin[l2D_XY] = MAX(k-1, nghz) ;
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
//<\MAIN_LOOP>
}




void smooth_zmin_boundary_location(int iDustSpecies) {


//<USER_DEFINED>
  INPUT2DINT(Dust_Zmin[iDustSpecies]);
  OUTPUT2DINT(Dust_Zmin[iDustSpecies]);
//<\USER_DEFINED>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  int llxm;
  int llxp;
//<\INTERNAL>

//<EXTERNAL>
  int* kmin = Dust_Zmin[iDustSpecies]->field_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = NGHZ;
  int nx = Nx;
  int ny = Ny;
  int nz = Nz;
  int nghy = NGHY;
  int nghz = NGHZ;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
//<\EXTERNAL>

//<CONSTANT>
//<\CONSTANT>

//<MAIN_LOOP>

  i = j = k = 0 ;
#ifdef Y
  for(j=0; j<size_y; j++) {
#endif
#ifdef X
    for(i=0; i<size_x; i++) {
#endif
      ll = l ;  
      llxm = lxm;
      llxp = lxp;
	  
      if (j > 0) {
	kmin[ll] = MIN(kmin[ll], kmin[llxm-pitch]+1) ;
	kmin[ll] = MIN(kmin[ll], kmin[ll  -pitch]+1) ;
	kmin[ll] = MIN(kmin[ll], kmin[llxp-pitch]+1) ;
      }

      kmin[ll] = MIN(kmin[ll], kmin[llxm]+1) ;
      kmin[ll] = MIN(kmin[ll], kmin[llxp]+1) ;

      if (j+1<size_y) {
	kmin[ll] = MIN(kmin[ll], kmin[llxm+pitch]+1) ;
	kmin[ll] = MIN(kmin[ll], kmin[ll  +pitch]+1) ;
	kmin[ll] = MIN(kmin[ll], kmin[llxp+pitch]+1) ;
      }
#ifdef X
    }
#endif
#ifdef Y
  }   
#endif
  //<\MAIN_LOOP>
}


void set_zmax_boundary_location(int iDustSpecies) {


//<USER_DEFINED>
  INPUT(Dust_Density[iDustSpecies]);
  OUTPUT2DINT(Dust_Zmax[iDustSpecies]);

  INPUT(Density);
//<\USER_DEFINED>

//<INTERNAL>
  int i;
  int j;
  int k;
  real rho_mid;
//<\INTERNAL>

//<EXTERNAL>
  real* density = Density->field_cpu;
  real* dust_density = Dust_Density[iDustSpecies]->field_cpu;
  int* kmax = Dust_Zmax[iDustSpecies]->field_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = NGHZ;
  int nx = Nx;
  int ny = Ny;
  int nz = Nz;
  int nghy = NGHY;
  int nghz = NGHZ;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
//<\EXTERNAL>

//<CONSTANT>
//<\CONSTANT>

//<MAIN_LOOP>

  i = j = k = 0;
#ifdef Y
    for(j=0; j<size_y; j++) {
#endif
#ifdef X
      for(i=0; i<size_x; i++) {
#endif
	// Find the vertical boundary:	
	k = nghz;
	rho_mid = density[l] ;
	for (; k < nz+nghz; k++) 
	  if (dust_density[l] < MIN_DUST_DENSITY*rho_mid*DUST_TO_GAS)
	    break ;
	// Add a buffer a buffer zone of one cell.
	//   This should help the boundary move outwards when the density
	//   is increasing by collecting mass in the cell. Consider increasing
	//   the number of buffer cells if mass diffusing across the boundary
	//   is preventing expansion.
	kmax[l2D_XY] = MIN(k+1, nz+nghz);
	k = kmax[l2D_XY];
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
//<\MAIN_LOOP>
}




void smooth_zmax_boundary_location(int iDustSpecies) {


//<USER_DEFINED>
  INPUT2DINT(Dust_Zmax[iDustSpecies]);
  OUTPUT2DINT(Dust_Zmax[iDustSpecies]);
//<\USER_DEFINED>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  int llxm;
  int llxp;
//<\INTERNAL>

//<EXTERNAL>
  int* kmax = Dust_Zmax[iDustSpecies]->field_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = NGHZ;
  int nx = Nx;
  int ny = Ny;
  int nz = Nz;
  int nghy = NGHY;
  int nghz = NGHZ;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
//<\EXTERNAL>

//<CONSTANT>
//<\CONSTANT>

//<MAIN_LOOP>

  i = j = k = 0 ;
#ifdef Y
  for(j=0; j<size_y; j++) {
#endif
#ifdef X
    for(i=0; i<size_x; i++) {
#endif
      ll = l ;  
      llxm = lxm;
      llxp = lxp;
	  
      if (j > 0) {
	kmax[ll] = MAX(kmax[ll], kmax[llxm-pitch]-1) ;
	kmax[ll] = MAX(kmax[ll], kmax[ll  -pitch]-1) ;
	kmax[ll] = MAX(kmax[ll], kmax[llxp-pitch]-1) ;
      }

      kmax[ll] = MAX(kmax[ll], kmax[llxm]-1) ;
      kmax[ll] = MAX(kmax[ll], kmax[llxp]-1) ;

      if (j+1<size_y) {
	kmax[ll] = MAX(kmax[ll], kmax[llxm+pitch]-1) ;
	kmax[ll] = MAX(kmax[ll], kmax[ll  +pitch]-1) ;
	kmax[ll] = MAX(kmax[ll], kmax[llxp+pitch]-1) ;
      }
#ifdef X
    }
#endif
#ifdef Y
  }   
#endif
  //<\MAIN_LOOP>
}

