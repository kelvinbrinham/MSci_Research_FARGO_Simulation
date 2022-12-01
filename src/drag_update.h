
#ifndef _FARGO_DRAG_UPDATE
#define _FARGO_DRAG_UPDATE

#ifdef DUST

#ifdef GPU
#define HOST   __host__
#define DEVICE __device__
#else
#define HOST  
#define DEVICE
#endif

// Semi-implicit update including time-dependence when back-reaction can be
// neglected
HOST DEVICE 
static inline void __test_particle_drag(real dt, const DragData* in_gas,
					DragData* dust, int nDust) {
  DragData gas = *in_gas ; // Local copy to avoid pointer aliasing.
  
  int iDust ;
  for (iDust = 0; iDust < nDust; iDust++) {
    real a_rel = dust[iDust].accel - gas.accel;
    real v_rel = dust[iDust].vel   - gas.vel    - a_rel * dt;

    real f_Drag_v = + exp  (- dt * dust[iDust].Kdrag) ;
    real f_Drag_a = - expm1(- dt * dust[iDust].Kdrag) / dust[iDust].Kdrag ;
    
    dust[iDust].vel = gas.vel + v_rel*f_Drag_v + a_rel*f_Drag_a ;
  }
}

// Short-friction time solution when back-reaction can be neglected. I.e. dust
// velocity is the terminal velocity.
HOST DEVICE 
static inline void __test_particle_sf(real dt, const DragData* in_gas,
				      DragData* dust, int nDust) {
  DragData gas = *in_gas ; // Local copy to avoid pointer aliasing.
  
  int iDust ;
  for (iDust = 0; iDust < nDust; iDust++) {
    real a_rel = dust[iDust].accel - gas.accel;

    dust[iDust].vel = gas.vel + a_rel / dust[iDust].Kdrag ;
  }
}

// Semi-implicit update including time-dependence and back-reaction
// Only works for one dust species
HOST DEVICE 
static inline void __feedback_drag(real dt, DragData* gas,
				   DragData* dust, int nDust) {

  int iDust = 0 ;
  
  // Compute centre of mass quantities:
  real v_com = gas->density*gas->vel   + dust[iDust].density*dust[iDust].vel;
  real a_com = gas->density*gas->accel + dust[iDust].density*dust[iDust].accel;
  real d_tot = gas->density            + dust[iDust].density ;

  v_com /= d_tot ;
  a_com /= d_tot ;

  // Compute delta v
  real a_rel = dust[iDust].accel - gas->accel;
  real v_rel = dust[iDust].vel   - gas->vel   - a_rel * dt;
  real eps = dust[iDust].density / d_tot ;

  real K = dust[iDust].Kdrag / (1 - eps) ;
  
  real f_Drag_v = + exp  (- dt * K) ;
  real f_Drag_a = - expm1(- dt * K) / K;

  dust[iDust].vel = v_rel*f_Drag_v + a_rel * f_Drag_a ;

  // Compute the gas velocity
  real gas_vel = v_com - eps * dust[iDust].vel ;

  // Convert delta_v to v
  dust[iDust].vel += gas_vel ;

  // Return the gas velocity
  gas->vel = gas_vel ;
}

// Short-friction time solution including back-reaction. I.e. dust velocity is
// the terminal velocity.
HOST DEVICE 
static inline void __feedback_sf(real dt, DragData* gas,
				 DragData* dust, int nDust) {
  int iDust ;  
  
  // Compute centre of mass quantities:
  real v_com = gas->density * gas->vel ;
  real a_com = gas->density * gas->accel ;
  real d_tot = gas->density ;
  
  for (iDust = 0; iDust < nDust; iDust++) {
    v_com += dust[iDust].density * dust[iDust].vel ;
    a_com += dust[iDust].density * dust[iDust].accel ;
    d_tot += dust[iDust].density ;
  }
  v_com /= d_tot ;
  a_com /= d_tot ;
  
  // Compute the delta Vs, and gas vel
  real gas_vel = v_com ;
  for (iDust = 0; iDust < nDust; iDust++) {
    dust[iDust].vel = (dust[iDust].accel - a_com) / dust[iDust].Kdrag ;
    gas_vel -= dust[iDust].vel * (dust[iDust].density / d_tot) ;
  }
  
  // Convert delta_v to v
  for (iDust = 0; iDust < nDust; iDust++)
    dust[iDust].vel += gas_vel ;

  // Return the new gas velocity
  gas->vel = gas_vel ;  
}

// Select desired drag implementation
#if (defined TEST_PARTICLE)

#if (defined FEEDBACK)
#error "Only one of TEST_PARTICLE and FEEDBACK may be defined"
#endif

#define MAX_DUST_SPECIES 0

#ifdef TERMINAL_VELOCITY
#define drag_update __test_particle_sf
#else
#define drag_update __test_particle_drag
#endif

#elif (defined FEEDBACK)

#define MAX_DUST_SPECIES 1

#ifdef TERMINAL_VELOCITY
#define drag_update __feedback_sf
#else
#define drag_update __feedback_drag
#endif

#else // NOT FEEBACK

#error "One of FEEDBACK or TEST_PARTICLE must be specified"

#endif // CLOSE NOT FEEDBACK

#else // IFDEF DUST
#define drag_update 
#endif // IFDEF DUST


#endif// _FARGO_DRAG_UPDATE
