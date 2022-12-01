#include "fargo3d.h"
void boundaries() {

#ifdef DUST
  int iDust ;
#endif

  if (!PERIODICZ) {
#ifdef Z
    if(Gridd.bc_down) {
      boundary_zmin();
#ifdef DUST
      for (iDust=0; iDust < NDUST; iDust++)
	boundary_zmin_dust(iDust, Drag_Coeff[iDust]) ;   
#endif
    }
    if(Gridd.bc_up) {
      boundary_zmax();
#ifdef DUST
      for (iDust=0; iDust < NDUST; iDust++)
	boundary_zmax_dust(iDust, Drag_Coeff[iDust]) ;   
#endif
    }
#endif
  }

  if (!PERIODICY) {
#ifdef Y
    if(Gridd.bc_left) {
      boundary_ymin();
#ifdef DUST
      for (iDust=0; iDust < NDUST; iDust++)
	boundary_ymin_dust(iDust, Drag_Coeff[iDust]) ;   
#endif
    }
    if(Gridd.bc_right) {
      boundary_ymax();
#ifdef DUST
      for (iDust=0; iDust < NDUST; iDust++)
	boundary_ymax_dust(iDust, Drag_Coeff[iDust]) ;   
#endif
    }
#endif
  }

#ifdef GHOSTSX 
  Fill_GhostsX();
#endif
}

