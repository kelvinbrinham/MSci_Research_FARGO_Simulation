#include "fargo3d.h"

void boundary_dust_zmax_internal(int iDustSpecies, real DragCoeff) ;
void set_boundary_location(int iDustSpecies) ;
void smooth_boundary_location(int iDustSpecies) ;
void set_boundary_condition(int iDustSpecies, real DragCoeff) ;


void CondInit() {
  int i,j,k;
  real* rho = Density->field_cpu;
  real* e = Energy->field_cpu;
  int index;

  real* vz = Vz->field_cpu;

  boolean GhostInclude = TRUE ;
  
  int begin_k =(GhostInclude ? 0 : NGHZ);
  int end_k = Nz+2*NGHZ-begin_k;
  int begin_j =(GhostInclude ? 0 : NGHY);
  int end_j = Ny+2*NGHY-begin_j;
  int begin_i = (GhostInclude ? 0 : NGHX);
  int end_i = Nx+2*NGHX-begin_i;

  for (j = begin_j; j<end_j; j++) {
    for (i = begin_i; i<end_i; i++) {	
      real R = Ymed(j) ;
      real h = ASPECTRATIO * pow(R/R0, FLARINGINDEX) ;
      real cs = h * sqrt(G*MSTAR/R) ;

      k = begin_k ;
      e[l] = cs ;
      rho[l] = 1 ;
      vz[l] = 0 ;
      
      real Sigma = 0 ;
      if (!GhostInclude) 
	Sigma = rho[l] * (Zmin(k+1) - Zmin(k)) ;
      
      for (k = begin_k + 1; k<end_k; k++) {
	real r_p = sqrt(R*R + Zmed(k)*Zmed(k)) ;
	real r_m = sqrt(R*R + Zmed(k-1)*Zmed(k-1)) ;
	
	real d_phi = (R/r_m - R/r_p) / (2*h*h) ;

	rho[l] = rho[lzm] * (1 - d_phi) / (1 + d_phi) ;       

	e[l]  = cs;
	vz[l] = 0 ;

	if (k >= NGHZ)
	  Sigma += rho[l] * (Zmin(k+1) - Zmin(k)) ;
      }
      for (k = begin_k; k<end_k; k++) 
	rho[l] *= SIGMA0 * pow(R/R0, SIGMASLOPE) / Sigma ;
    }
  }

#ifdef DUST
  real K0 = 100.;
  int iDust ;
  for (iDust = 0; iDust < NDUST; iDust++) {
    for (k = begin_k; k<end_k; k++) {
      for (j = begin_j; j<end_j; j++) {
	for (i = begin_i; i<end_i; i++) {	
	  real Sigma = SIGMA0 * pow(Ymed(j)/R0, SIGMASLOPE);
	  Dust_Vz[iDust]->field_cpu[l] = 0 ;
	  if (k == NGHZ)
	    Dust_Density[iDust]->field_cpu[l] = 
	      Sigma * DUST_TO_GAS / (Zmin(k+1) - Zmin(k));
	  else
	    Dust_Density[iDust]->field_cpu[l] = Sigma * DUST_TO_GAS *1e-100 ;
	}
      }
    }

    i=begin_i; j=begin_j; k=begin_k;
    Drag_Coeff[iDust] = K0 * pow(10, iDust) / (rho[l]*e[l]) ;
  }

  boundary_zmax_dust = boundary_dust_zmax_internal ;
#endif
  
}


void boundary_dust_zmax_internal(int iDustSpecies, real DragCoeff) {
  // Find the new Z-location of the boundary and smooth twice
  set_boundary_location(iDustSpecies) ; 
  smooth_boundary_location(iDustSpecies);
  smooth_boundary_location(iDustSpecies);
  // Zero-out the region outside the new domain
  set_boundary_condition(iDustSpecies,DragCoeff) ;
}


void set_boundary_condition(int iDustSpecies, real DragCoeff) {

//<USER_DEFINED>
  INPUT(Dust_Density[iDustSpecies]);
#ifdef X
  INPUT(Dust_Vx[iDustSpecies]);
  OUTPUT(Dust_Vx[iDustSpecies]);
  INPUT(Vx) ;
#endif
#ifdef  Y
  INPUT(Dust_Vy[iDustSpecies]);
  OUTPUT(Dust_Vy[iDustSpecies]);
#endif

  INPUT(Dust_Vz[iDustSpecies]);
  OUTPUT(Dust_Vz[iDustSpecies]);

  OUTPUT(Dust_Density[iDustSpecies]);
  INPUT2DINT(Dust_Zmax[iDustSpecies]);

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
  real* dust_density = Dust_Density[iDustSpecies]->field_cpu;
#ifdef X
  real* gas_vx  = Vx->field_cpu;
  real* dust_vx = Dust_Vx[iDustSpecies]->field_cpu;
#endif
#ifdef Y
  real* dust_vy = Dust_Vy[iDustSpecies]->field_cpu;
#endif
  real* dust_vz = Dust_Vz[iDustSpecies]->field_cpu;
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
	k0 = kmax[l2D_XY];

	// Reset the values outside the boundary to make life easy:
	for (k = k0 ; k < nz+2*nghz; k++) {
	  ll = l ;
	  dust_density[ll] = 1e-20;
	  dust_vz[ll] = 0;
#ifdef Y
	  dust_vy[ll] = 0;
#endif
#ifdef X
	  dust_vx[ll] = gas_vx[ll] ;
#endif
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


void set_boundary_location(int iDustSpecies) {


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
	  if (dust_density[l] < 1e-10*rho_mid*DUST_TO_GAS)
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




void smooth_boundary_location(int iDustSpecies) {


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

