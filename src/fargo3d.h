//System includes
#ifdef GPU
#include <cuda.h>
#include <driver_functions.h> //for some CUDA structs
#include <cuda_runtime_api.h>
#endif

#ifndef __APPLE__
#include <malloc.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>
#include <time.h>

#include "define.h"
#include "types_def.h"
#include "fondam.h"

#define NVCC_PROBLEMS //Nvcc is problematic with double references.

#if  !(defined(__NOPROTO) && defined(NVCC_PROBLEMS))
#include "prototypes.h"
#else
extern "C" void Input_GPU(Field *, int, const char *);
extern "C" void Output_GPU(Field *, int, const char *);
extern "C" void Input2D_GPU(Field2D *, int, const char *);
extern "C" void Output2D_GPU(Field2D *, int, const char *);
extern "C" void Input2DInt_GPU(FieldInt2D *, int, const char *);

extern "C" void Output2DInt_GPU(FieldInt2D *, int, const char *);

extern "C" void check_errors(const char*);

extern "C" real reduction_full_SUM (Field *, int, int, int, int);
extern "C" void reduction_SUM_cpu (Field *, int, int, int, int);
extern "C" real reduction_full_MIN (Field *, int, int, int, int);
extern "C" void reduction_MIN_cpu (Field *, int, int, int, int);

extern "C" void cfl_b(void);
#endif

// MPI includes
#ifdef PARALLEL
#include <mpi.h>
#else
#include "mpi_dummy.h"
#endif 

// FFT includes
#ifdef FFTW

#ifdef GPU
#include <cufft.h>
#ifdef PARALLEL
#error "FFT only works on CPU or 1 GPU"
#endif

 
#else // not GPU

#include <fftw3-mpi.h>
#ifndef PARALLEL
#error Self-gravity with fft requires MPI
#endif

#endif
#endif 

#include "param.h"

#ifndef __LOCAL // ONLY VAR.C HAS __LOCAL
#include "structs.h"
#include "global_ex.h"
#else
#include "structs.h"
#include "global.h"
#endif
