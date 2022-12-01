#define __GPU
#define __NOPROTO

#include <cufft.h>
#include <cuComplex.h>

#include "fargo3d.h"
 
typedef struct SG_data_gpu {  
#ifdef SELF_GRAVITY
  // These will exist on only one the CPU or GPU
  real* source;
  real* kernel;
  real* pot;

#ifdef GPU
  cufftHandle fwd_cuplan ;
  cufftHandle bwd_cuplan ;
#endif
  
  real soft ; // Softening parameter (softening length / R)
  real norm ; // Normalization for fft work
  real delta_u; // u-space separation between radial grid cells
  real delta_t; // phi-space separation between radial grid cells
  real Rmid0;   // Middle of first radial grid cell
#endif
} SGData ;
static SGData sg_data_gpu ;
static FFTGrid fft_grid ;

#ifdef GPU

#define ymin(i) ymin_s[(i)]
CONSTANT(real, ymin_s, 7686);

#ifdef FLOAT
typedef cufftComplex complex ;
#else
typedef cufftDoubleComplex complex ;
#endif

#endif

extern "C" void prs_exit(int) ;
extern "C" void CheckNans (char *string) ;

// Prototypes for internal functions
void sgkernel_gpu();
void sgsource_gpu() ;
void sgpot_gpu() ;

// Compute the potential due to self-gravity
extern "C" void SelfGravityPotential_gpu() {
#ifdef GPU
  sgsource_gpu() ;
  sgpot_gpu() ; 
#endif
}


extern "C" void SelfGravityInit_gpu() {
#ifdef GPU
#if defined(SELF_GRAVITY)
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

#ifdef FLOAT
  if (cufftPlan2d(&(sg_data_gpu.fwd_cuplan), 2*NY, NX, CUFFT_R2C) != CUFFT_SUCCESS) 
#else
  if (cufftPlan2d(&(sg_data_gpu.fwd_cuplan), 2*NY, NX, CUFFT_D2Z) != CUFFT_SUCCESS) 
#endif
    {

      fprintf(stderr, "Error: failed to create fwd cufft_cuplan for SG") ;
      prs_exit(1) ;
    }

#ifdef FLOAT
  if (cufftPlan2d(&(sg_data_gpu.bwd_cuplan), 2*NY, NX, CUFFT_C2R) != CUFFT_SUCCESS) 
#else
  if (cufftPlan2d(&(sg_data_gpu.bwd_cuplan), 2*NY, NX, CUFFT_Z2D) != CUFFT_SUCCESS) 
#endif
    {
      fprintf(stderr, "Error: failed to create bwd cufft_cuplan for SG") ;
      prs_exit(1) ;
    }


  // Make sure that we have enough space for boundaries  
  size_t local_size ;
  cufftGetSize(sg_data_gpu.fwd_cuplan, &local_size) ;

  cudaMalloc((void**)&(sg_data_gpu.source), local_size) ;
  check_errors("SG_alloc");
  cudaMalloc((void**)&(sg_data_gpu.kernel), local_size) ;
  check_errors("SG_alloc");
  cudaMalloc((void**)&(sg_data_gpu.pot), local_size) ;
  check_errors("SG_alloc");

			   
  // Store grid sizes
  fft_grid.Ny = NX ;
  fft_grid.stride_y = (NX/2 + 1);
  fft_grid.local_Nx = 2*NY ;

  // Save the additional constants
  sg_data_gpu.soft = SELFGRAVITYSOFTENING;
  sg_data_gpu.delta_u = log(YMAX/YMIN) / NY ;
  sg_data_gpu.delta_t = 2*M_PI/NX ;
  sg_data_gpu.Rmid0 = (YMIN + YMIN*exp(sg_data_gpu.delta_u))/2 ;
  sg_data_gpu.norm =
    sg_data_gpu.delta_u*sg_data_gpu.delta_t * sg_data_gpu.Rmid0/(2*NX*NY);
  
  // Initialize the kernel
  sgkernel_gpu() ;
  printf("Initialised self-gravity\n") ;
#endif
#endif
}

///============================================================================
// Internal functions

#ifdef GPU

__global__  void __sgkernel(SGData sg_data, FFTGrid fft_grid)
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  int j = threadIdx.y + blockIdx.y * blockDim.y;

  if (i < fft_grid.local_Nx && j < fft_grid.Ny) {
    real du   = sg_data.delta_u ;
    real dphi = sg_data.delta_t ;

    real eps = sg_data.soft*sg_data.soft ;

    int stride_fft = fft_grid.Ny ;
  
    // u runs over range [-umax, umax), where umax = log(YMAX/YMIN).
    real u = du * i ;
    if (i >= fft_grid.local_Nx/2) u -= fft_grid.local_Nx * du ;
      
    real phi = dphi * j ;
      
     sg_data.source[i*stride_fft + j] = 
       - G * sg_data.norm / sqrt(2*(cosh(u) - cos(phi)) + eps*exp(u)) ;
  }
}


// Setup the gravity kernel
void sgkernel_gpu() {
#ifdef SELF_GRAVITY
  dim3 block (BLOCK_X, BLOCK_Y);
  dim3 grid ((fft_grid.local_Nx+block.x-1)/block.x,
	     (fft_grid.Ny+block.y-1)/block.y) ;

  __sgkernel<<<grid,block>>>(sg_data_gpu, fft_grid) ;
  check_errors("__sgkernel") ;
      
  // Do the FFT of the kernel
#ifdef FLOAT
  cufftExecR2C(sg_data_gpu.fwd_cuplan, sg_data_gpu.source, (complex*) sg_data_gpu.kernel);
#else
  cufftExecD2Z(sg_data_gpu.fwd_cuplan, sg_data_gpu.source, (complex*) sg_data_gpu.kernel);
#endif
  check_errors("sgkernel_fft") ;
#endif
}

__global__ void __sgsource(SGData sg_data, FFTGrid fft_grid, 
			   real* rho, int pitch, int stride) {

  int ii = threadIdx.x + blockIdx.x * blockDim.x;
  int jj = threadIdx.y + blockIdx.y * blockDim.y;
  
  int stride_fft = fft_grid.Ny ;
  int Nx = fft_grid.local_Nx / 2 ; // Must be NY
  double du = sg_data.delta_u ;

  if (jj < fft_grid.Ny) {
    if (ii < Nx) {
      int i = 0, j=0, k = 0 ;
      i = jj + NGHX ; // Transposed
      j = ii + NGHY ;
    
      sg_data.pot[ii*stride_fft + jj] = rho[l] * exp(1.5 * du * ii) ;	
    }
    else {
      sg_data.pot[ii*stride_fft + jj] = 0 ;
    }
  }
}

// Compute the fft of the source term (density)
void sgsource_gpu() {
#ifdef SELF_GRAVITY
  INPUT(Density);

  dim3 block (BLOCK_X, BLOCK_Y);
  dim3 grid ((fft_grid.local_Nx+block.x-1)/block.x,
	     (fft_grid.Ny+block.y-1)/block.y) ;

  __sgsource<<<grid,block>>>(sg_data_gpu, fft_grid, Density->field_gpu,	
			     Pitch_gpu, Stride_gpu) ;
  check_errors("__sgsource") ;
      
  // Do the FFT of the source term
#ifdef FLOAT
  cufftExecR2C(sg_data_gpu.fwd_cuplan, sg_data_gpu.pot, (complex*) sg_data_gpu.source) ;
#else
  cufftExecD2Z(sg_data_gpu.fwd_cuplan, sg_data_gpu.pot, (complex*) sg_data_gpu.source) ;
#endif
  check_errors("sgsource_fft") ;
#endif
}

__global__ void __sgpot1(SGData sg_data, FFTGrid fft_grid) {

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  int j = threadIdx.y + blockIdx.y * blockDim.y;
  
  complex* src  = (complex*) sg_data.source ;
  complex* kern = (complex*) sg_data.kernel ;

  if (i < fft_grid.local_Nx && j < fft_grid.stride_y) {
    int idx = i*fft_grid.stride_y + j ;
    complex result ;
    result.x = src[idx].x*kern[idx].x - src[idx].y*kern[idx].y ;
    result.y = src[idx].x*kern[idx].y + src[idx].y*kern[idx].x ;
    //src[idx] = kern[idx] ;
    src[idx] =  result ;
  }

}

__global__ void __sgpot2(SGData sg_data, FFTGrid fft_grid,
			 real* pot, int size_x, int size_y,
			 int pitch, int stride) {
  
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  int j = threadIdx.y + blockIdx.y * blockDim.y;
  int k = 0 ;

  if (i < size_x && j < size_y) {

    int ii = i - NGHX ;
    int jj = j - NGHY ;
    if (jj < 0) jj += fft_grid.local_Nx ;

#ifdef GHOSTSX
    if (ii <  0) ii += fft_grid.Ny ;
    if (ii > fft_grid.Ny) ii -= fft_grid.Ny ;
#endif

    int stride_fft = fft_grid.Ny ;
    real Rmid0 = sg_data.Rmid0 ;
    real sqrt_r = sqrt(sg_data.Rmid0/ymed(j)) ;

    pot[l] += sg_data.pot[jj*stride_fft + ii] * sqrt_r ;
  }
}

// Compute the potential via the convolution of source and kernel
void sgpot_gpu() {
#ifdef SELF_GRAVITY
  INPUT(Pot);
  OUTPUT(Pot);

#ifdef BIGMEM
#define ymin_d &Ymin_d
#endif
  CUDAMEMCPY(ymin_s, ymin_d, sizeof(real)*(Ny+2*NGHY+1), 0, \
	     cudaMemcpyDeviceToDevice);

  // Multiply source by kernel
  dim3 block (BLOCK_X, BLOCK_Y);
  dim3 grid ((fft_grid.local_Nx+block.x-1)/block.x,
	     (fft_grid.stride_y+block.y-1)/block.y) ;

  __sgpot1<<<grid,block>>>(sg_data_gpu, fft_grid) ;
  check_errors("__sgpot1") ;

  // Do the FFT of the source term
#ifdef FLOAT
  cufftExecC2R(sg_data_gpu.bwd_cuplan, (complex*) sg_data_gpu.source, sg_data_gpu.pot) ;
#else
  cufftExecZ2D(sg_data_gpu.bwd_cuplan, (complex*) sg_data_gpu.source, sg_data_gpu.pot) ;
#endif
  check_errors("sgpot_fft") ;

  dim3 grid2 ((Nx+2*NGHX+block.x-1)/block.x,
	      ((Ny+2*NGHY)+block.y-1)/block.y) ;

  __sgpot2<<<grid2,block>>>(sg_data_gpu, fft_grid, Pot->field_gpu, 
			    Nx + 2*NGHX, Ny+2*NGHY,
			    Pitch_gpu, Stride_gpu) ;
  check_errors("__sgpot2") ;
#endif
}


#endif // GPU
