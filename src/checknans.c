#include "fargo3d.h"

int CheckNansField (Field *f) {
  int i,j,k;
  for(k=0; k<Nz+2*NGHZ; k++) {
    for(j=0; j<Ny+2*NGHY; j++) {
      for(i=0; i<Nx; i++) {
	if (isnan(f->field_cpu[l])) {
	  return l;
	}
      }
    }
  }
  return -1;
}

void CheckNans (char *string){
  static int count=0;
  int i;
  Field *g;
  g = ListOfGrids;
  while (g != NULL) {
    if (*(g->owner) == g)
      INPUT(g) ;
    if ((i = CheckNansField (g)) >= 0) {
      printf ("Found NaNs in grid %s at location (i=%d, j=%d, k=%d)\n",	\
	      g->name, i%Nx, (i/Nx)%(Ny+2*NGHY), i/(Nx*(Ny+2*NGHY)));
      printf ("position: after call to %s\n", string);
      prs_exit(1) ;
    }
    g = g->next;
  }
}
