
#ifndef _FARGO_DRAG_LAW
#define _FARGO_DRAG_LAW

#ifdef GPU
#define HOST   __host__
#define DEVICE __device__
#else
#define HOST  
#define DEVICE
#endif

#if defined(STOKES_NUMBER)
HOST DEVICE 
static inline real drag_constant(real coeff, real r, real rho_g, real c_s) {
  return coeff * sqrt((G * MSTAR) / (r*r*r)) ;
}
#elif defined(EPSTEIN)
HOST DEVICE 
static inline real drag_constant(real coeff, real r, real rho_g, real c_s) {
  return coeff * rho_g * c_s ;
}
#elif defined(EPSTEIN_2D)
HOST DEVICE 
static inline real drag_constant(real coeff, real r, real rho_g, real c_s) {
  return coeff * rho_g * sqrt((G * MSTAR) / (r*r*r)) ;
}
#else
HOST DEVICE 
static inline real drag_constant(real coeff, real r, real rho_g, real c_s) {
  return coeff ;
}
#endif

#endif// _FARGO_DRAG_LAW
