
#include "fargo3d.h"


typedef struct SG_data_cpu {  
#ifdef SELF_GRAVITY
  // These will exist on only one the CPU or GPU
  real* source;
  real* kernel;
  real* pot;

#ifdef FFTW
  fftw_plan fwd_plan ;
  fftw_plan bwd_plan ;
#endif
  struct sg_comm src_comm ;
  struct sg_comm pot_comm ;  

  
  real soft ; // Softening parameter (softening length / R)
  real norm ; // Normalization for fft work
  real delta_u; // u-space separation between radial grid cells
  real delta_t; // phi-space separation between radial grid cells
  real Rmid0;   // Middle of first radial grid cell
#endif
} SGData ;
static SGData sg_data_cpu ;
static FFTGrid fft_grid ;


// Prototypes for internal functions
void sgkernel();
void sgsource() ;
void sgpot() ;


void SelfGravityInit_cpu() {
#ifndef GPU
#ifndef FFTW
  fprintf(stderr, "Error: self-gravity requires compilation with FFTW\n") ;
  prs_exit(1);
  
#elif defined(SELF_GRAVITY)
  /* Self-gravity only works in 2D cylindrical co-ordinates for now */

#ifndef CYLINDRICAL
#error Self-gravity requires cylindrical co-ordinates
#endif

#if !(defined (X) && defined(Y))
#error X and Y must be active for self-gravity
#endif

#ifdef Z
#error Self-gravity does not work with Z dimensions
#endif
  
  if (!((toupper(*SPACING) == 'L') && (toupper(*(SPACING+1)) == 'O'))) {
    fprintf(stderr, "Error: self-gravity requires a logarithmic grid\n");
    prs_exit(1) ;
  }

  // Set up the forward and backward plan
  if (CPU_Rank == 0) 
    printf("Computing FFT plan...\n") ;
  fftw_mpi_init() ;

  ptrdiff_t local_n0, local_0_start ;
  ptrdiff_t local_size = fftw_mpi_local_size_2d(2*NY, NX/2+1, MPI_COMM_WORLD,
						&local_n0, &local_0_start) ;

  // Make sure that we have enough space for boundaries
  local_size = MAX(local_size, (NY+2*NGHY)*(NX/2+1)) ;
  sg_data_cpu.source = fftw_alloc_real(2*local_size);
  sg_data_cpu.kernel = fftw_alloc_real(2*local_size);
  sg_data_cpu.pot    = fftw_alloc_real(2*local_size);

  if (sg_data_cpu.source == NULL || sg_data_cpu.kernel == NULL ||
      sg_data_cpu.pot == NULL) {
    fprintf(stderr, "Error allocating fftw arrays for self-gravity\n");
    prs_exit(1);
  }

  sg_data_cpu.fwd_plan =
    fftw_mpi_plan_dft_r2c_2d(2*NY, NX,
			     sg_data_cpu.source,
			     (fftw_complex*) sg_data_cpu.source,
			     MPI_COMM_WORLD,
			     FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT) ;

  sg_data_cpu.bwd_plan =
    fftw_mpi_plan_dft_c2r_2d(2*NY, NX,
			     (fftw_complex*) sg_data_cpu.pot, sg_data_cpu.pot,
			     MPI_COMM_WORLD,
			     FFTW_MEASURE | FFTW_MPI_TRANSPOSED_IN) ;
			   
  // Store grid sizes
  fft_grid.Ny = Nx ;
  fft_grid.stride_y = 2*(NX/2 + 1);
  fft_grid.local_Nx = local_n0 ;
  fft_grid.local_i_start = local_0_start ;
  fft_grid.local_size = 2*local_size ;

  // Save the additional constants
  sg_data_cpu.soft = SELFGRAVITYSOFTENING;
  sg_data_cpu.delta_u = log(YMAX/YMIN) / NY ;
  sg_data_cpu.delta_t = 2*M_PI/NX ;
  sg_data_cpu.Rmid0 = (YMIN + YMIN*exp(sg_data_cpu.delta_u))/2 ;
  sg_data_cpu.norm =
    sg_data_cpu.delta_u * sg_data_cpu.delta_t * sg_data_cpu.Rmid0 / (2*NX*NY);
  
  // TODO:
  //   Work out who we need to grab data from
  
  // Find out the FFTW domains of each processor
  int *fft_ends  = malloc(sizeof(int) * CPU_Number + 1) ;

  fft_ends[0] = 0 ;
  int end = fft_grid.local_i_start + fft_grid.local_Nx ;
  MPI_Allgather(&end, 1, MPI_INT, fft_ends+1, 1, MPI_INT, MPI_COMM_WORLD) ;

  // Now get the Mesh ends of each process
  int *mesh_ends = malloc(sizeof(int) * CPU_Number + 1) ;
  MPI_Allgather(&Ny, 1, MPI_INT, mesh_ends+1,  1, MPI_INT, MPI_COMM_WORLD) ;

  // From the mesh-ends construct the starts / ends
  int i;
  int tot = 0 ;
  for (i = 0; i < CPU_Number; ++i) {
    mesh_ends[i] = tot ;
    tot += mesh_ends[i+1];
  }
  mesh_ends[CPU_Number] = tot ;

  // Fill the communicators for the potential source terms (density)
  struct sg_comm comm ;
  comm.sends = malloc(sizeof(struct sg_domain) * CPU_Number) ;
  comm.recvs = malloc(sizeof(struct sg_domain) * CPU_Number) ;

  // Work out where we send our data:
  comm.Nsend = 0;
  for (i = 0; i < CPU_Number; i++) {
    // Overlap range:
    int start = MAX(fft_ends[i],   mesh_ends[CPU_Rank]) ;
    int end   = MIN(fft_ends[i+1], mesh_ends[CPU_Rank+1]) ;

    if (end > start) {
      struct sg_domain *send = &comm.sends[comm.Nsend++] ;
      send->rank   = i ;
      send->start  = start - mesh_ends[CPU_Rank];
      send->end    = end   - mesh_ends[CPU_Rank];
    }
  }

  // Work out where our data is recieved from:
  comm.Nrecv = 0 ;
  for (i = 0; i < CPU_Number; i++) {
    // Overlap range:
    int start = MAX(fft_ends[CPU_Rank],   mesh_ends[i]) ;
    int end   = MIN(fft_ends[CPU_Rank+1], mesh_ends[i+1]) ;

    if (end > start) {
      struct sg_domain *recv = &comm.recvs[comm.Nrecv++] ;
      recv->rank   = i ;
      recv->start  = start - fft_ends[CPU_Rank];
      recv->end    = end   - fft_ends[CPU_Rank];
    }
  }
  sg_data_cpu.src_comm = comm ;

  // Now compute the transfers for the return, including ghost cells in radial
  // direction
  comm.sends = malloc(sizeof(struct sg_domain) * (CPU_Number + 1)) ;
  comm.recvs = malloc(sizeof(struct sg_domain) * (CPU_Number + 1)) ;

  // Work out where we send our data:
  comm.Nsend = 0 ;
  for (i = 0; i < CPU_Number; i++) {
    // Overlap range:
    int start = MAX(fft_ends[CPU_Rank],   mesh_ends[i]   - NGHY) ;
    int end   = MIN(fft_ends[CPU_Rank+1], mesh_ends[i+1] + NGHY) ;

    if (end > start) {
      struct sg_domain *send = &comm.sends[comm.Nsend++] ;
      send->rank   = i ;
      send->start  = start - fft_ends[CPU_Rank];
      send->end    = end   - fft_ends[CPU_Rank];
    }
  }

  // Work out where our data is recieved from:
  comm.Nrecv = 0;
  for (i = 0; i < CPU_Number; i++) {
    // Overlap range:
    int start = MAX(fft_ends[i],   mesh_ends[CPU_Rank]   - NGHY) ;
    int end   = MIN(fft_ends[i+1], mesh_ends[CPU_Rank+1] + NGHY) ;

    if (end > start) {
      struct sg_domain *recv = &comm.recvs[comm.Nrecv++] ;
      recv->rank   = i ;
      recv->start  = start - (mesh_ends[CPU_Rank] - NGHY) ;
      recv->end    = end   - (mesh_ends[CPU_Rank] - NGHY);
    }
  }

  // Finally handle the edge case of periodic boundary conditions:
  if (CPU_Rank == CPU_Number-1) {   
    struct sg_domain *send = &comm.sends[comm.Nsend++] ;
    send->rank  = 0 ;
    send->start = fft_ends[CPU_Number] - NGHY - fft_ends[CPU_Rank] ;
    send->end   = fft_ends[CPU_Number]        - fft_ends[CPU_Rank] ;
  }
  if (CPU_Rank == 0) {
    struct sg_domain *recv = &comm.recvs[comm.Nrecv++] ;
    recv->rank   = CPU_Number - 1 ;
    recv->start  = 0 ;
    recv->end    = NGHY ;
  }  
  sg_data_cpu.pot_comm = comm ;
  
  free(fft_ends) ;
  free(mesh_ends) ;
  
  
  // Initialize the kernel
  sgkernel() ;
#endif
#else
  fprintf(stderr, "Error: self-gravity is running on the GPU\n");
  prs_exit(1) ;
#endif
}

// Compute the potential due to self-gravity
void SelfGravityPotential_cpu() {
#ifndef GPU
  sgsource() ;
  sgpot() ;
#else
  fprintf(stderr, "Error: self-gravity is running on the GPU\n");
  prs_exit(1) ;
#endif 
}


///============================================================================
// Internal functions

#ifndef GPU

// Setup the gravity kernel
void sgkernel() {
#ifdef SELF_GRAVITY

  fftw_plan kern_plan =
    fftw_mpi_plan_dft_r2c_2d(2*NY, NX,
			     sg_data_cpu.kernel,
			     (fftw_complex*)sg_data_cpu.kernel,
			     MPI_COMM_WORLD,
			     FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT) ;

  int ii, jj ;
  
  int size_x = fft_grid.local_Nx;
  int size_y = fft_grid.Ny ;

  int i_start = fft_grid.local_i_start ;
  
  real du   = sg_data_cpu.delta_u ;
  real dphi = sg_data_cpu.delta_t ;

  real eps = sg_data_cpu.soft*sg_data_cpu.soft ;

  int stride_fft = fft_grid.stride_y ;
  
  for (ii=0; ii < size_x; ii++) 
    for (jj=0; jj < size_y; jj++) {
      // u runs over range [-umax, umax), where umax = log(YMAX/YMIN).
      real u = du * (ii + i_start) ;
      if (ii + i_start >= NY) u -= 2*NY * du ;
      
      real phi = dphi * jj ;
      
      sg_data_cpu.kernel[ii*stride_fft + jj] = 
	- G * sg_data_cpu.norm / sqrt(2*(cosh(u) - cos(phi)) + eps*exp(u)) ;
    }

  // Do the FFT of the kernel
  fftw_execute(kern_plan) ;

  fftw_destroy_plan(kern_plan) ;  
#endif
}


// Compute the fft of the source term (density)
void sgsource() {
#ifdef SELF_GRAVITY
  INPUT(Density);

  real *rho = Density->field_cpu ;

  int size_x = fft_grid.local_Nx;
  int size_y = fft_grid.Ny ;
  int local_i_start = fft_grid.local_i_start ;
  int stride_fft = fft_grid.stride_y ;

  double du = sg_data_cpu.delta_u ;


  // Gather the density from the parent processor
  struct sg_comm comm = sg_data_cpu.src_comm ;
  
  MPI_Request* reqsts = malloc((comm.Nsend+comm.Nrecv)*sizeof(MPI_Request));
  MPI_Status* stats = malloc((comm.Nsend+comm.Nrecv)*sizeof(MPI_Status));
  
  // Post the recieves first:
  int iproc;
  for (iproc=0; iproc < comm.Nrecv; iproc++) {
    struct sg_domain recv = comm.recvs[iproc] ;

    MPI_Irecv(sg_data_cpu.source + stride_fft*recv.start,
	      stride_fft*(recv.end-recv.start)*sizeof(real), MPI_BYTE,
	      recv.rank, 42, MPI_COMM_WORLD, reqsts+iproc) ;
  }

  // Setup the data and send it:
  for (iproc=0; iproc < comm.Nsend; iproc++) {
    struct sg_domain send = comm.sends[iproc] ;

    int i, j, k = 0 ;
    int ii, jj ;
    for (ii = 0; ii < NX; ii++) 
      for (jj = send.start; jj < send.end; jj++) {
	i = ii + NGHX ;
	j = jj + NGHY ;	
	sg_data_cpu.pot[jj*stride_fft + ii] = rho[l] ;
      }

    MPI_Isend(sg_data_cpu.pot + send.start*stride_fft,
	      stride_fft*(send.end-send.start)*sizeof(real), MPI_BYTE,
	      send.rank, 42, MPI_COMM_WORLD, reqsts+iproc+comm.Nrecv) ;
  }

  MPI_Waitall(comm.Nsend + comm.Nrecv, reqsts, stats) ;

  free(stats) ;
  free(reqsts) ;

  int i, j ;  
  for (i = 0; i < size_x; i++) {
    if (i + local_i_start < NY) {     
      double r3 = exp(1.5 * du * (i + local_i_start)) ;
      for (j = 0; j < size_y; j++)
	sg_data_cpu.source[i*stride_fft + j] *= r3 ;
    }
    else {
      for (j = 0; j < size_y; j++)
	sg_data_cpu.source[i*stride_fft + j] = 0 ;
    }
  }

  // Do the FFT of the source term
  fftw_execute(sg_data_cpu.fwd_plan) ;
#endif
}


// Compute the potential via the convolution of source and kernel
void sgpot() {
#ifdef SELF_GRAVITY
  INPUT(Pot);
  OUTPUT(Pot);

  int i, j, k ;
  int size = fft_grid.local_size ;
  int stride_fft = fft_grid.stride_y ;

  // The convolution is just a multiplication. Here we handle the fact that the
  // the fourier transforms of the source / kernel are complex.
  //
  // Note:
  //   Here the FFT data is actually transposed. However, we don't care since we
  //   operate element-wise anyway
  real *src    = sg_data_cpu.source ;
  real *kern   = sg_data_cpu.kernel ;
  real *pot   = sg_data_cpu.pot ;
  for (i = 0; i < size; i +=2) {
    int i_real = i ;
    int i_imag = i+1;

    pot[i_real] = (src[i_real]*kern[i_real] - src[i_imag]*kern[i_imag]);
    pot[i_imag] = (src[i_real]*kern[i_imag] + src[i_imag]*kern[i_real]);
  }

  fftw_execute(sg_data_cpu.bwd_plan);


  // Now send the data back to the host processor
  struct sg_comm comm = sg_data_cpu.pot_comm ;
  
  MPI_Request* reqsts = malloc((comm.Nsend+comm.Nrecv)*sizeof(MPI_Request));
  MPI_Status* stats = malloc((comm.Nsend+comm.Nrecv)*sizeof(MPI_Status));
  
  // Post the recieves first:
  int iproc;
  for (iproc=0; iproc < comm.Nrecv; iproc++) {
    struct sg_domain recv = comm.recvs[iproc] ;

    MPI_Irecv(sg_data_cpu.source + stride_fft*recv.start,
	      stride_fft*(recv.end-recv.start)*sizeof(real), MPI_BYTE,
	      recv.rank, 42, MPI_COMM_WORLD, reqsts+iproc);
  }

  // Setup the data and send it:
  for (iproc=0; iproc < comm.Nsend; iproc++) {
    struct sg_domain send = comm.sends[iproc] ;

    MPI_Isend(sg_data_cpu.pot + send.start*stride_fft,
	      stride_fft*(send.end-send.start)*sizeof(real), MPI_BYTE,
	      send.rank, 42, MPI_COMM_WORLD, reqsts+iproc+comm.Nrecv);
  }

  MPI_Waitall(comm.Nsend + comm.Nrecv, reqsts, stats) ;

  free(stats) ;
  free(reqsts) ;


  // Make a pointer to source as it now contains the potential:
  real *sg_pot = sg_data_cpu.source ;
  pot = Pot->field_cpu ;  

  // Now renormalize the potential and add back to the potential arrays
  i = j = k = 0 ;
  real Rmid0 = sg_data_cpu.Rmid0 ;
  for (j = 0; j < Ny + 2*NGHY; j++) {
    double sqrt_r = sqrt(Rmid0/ymed(j)) ;

    for (i = NGHX; i < Nx + NGHX; i++)  {
      int ii = i - NGHX ;
      pot[l] += sg_pot[j*stride_fft + ii] * sqrt_r ;
    }
  }

#ifdef GHOSTSX
  for (k = 0; j < Nz + 2*NGHZ; k++) {
    for (j = 0; j < Ny + 2*NGHY; j++) {
      for (i = 0; i < NGHX; i++) {
	int lghost1 = l;
	int lcopy1 = l + NX; 
	int lcopy2 = lghost1 + NGHX;
	int lghost2 = lcopy2 + NX; 

	pot[lghost1] = pot[lcopy1] ;
	pot[lghost2] = pot[lcopy2] ;
      }
    }
  }
#endif
#endif
}

#endif
