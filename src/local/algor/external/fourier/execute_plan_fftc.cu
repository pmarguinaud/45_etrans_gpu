#define cufftSafeCall(err) __cufftSafeCall(err, __FILE__, __LINE__)
#include <cufft.h>
#include <stdio.h>

static const char *_cudaGetErrorEnum (cufftResult error)
{
  switch (error)
    {
#define cr(x) case CUFFT_##x: return #x
      cr (SUCCESS);
      cr (INVALID_PLAN);
      cr (ALLOC_FAILED);
      cr (INVALID_TYPE);
      cr (INVALID_VALUE);
      cr (INTERNAL_ERROR);
      cr (EXEC_FAILED);
      cr (SETUP_FAILED);
      cr (INVALID_SIZE);
      cr (UNALIGNED_DATA);
#undef cr
    }

  return "UNKNOWN";
}

static inline void __cufftSafeCall (cufftResult err, const char * file, const int line)
{
  if( CUFFT_SUCCESS != err) 
    {
      fprintf (stderr, "CUFFT error at %s:%d\n", file, line);
      fprintf (stderr, "CUFFT error %d %s\n", err, _cudaGetErrorEnum (err)); 
      cudaDeviceReset (); 
    }
 }

extern "C"
void
#ifdef TRANS_SINGLE
execute_plan_fftc_ (cufftHandle * PLANp, int * ISIGNp, cufftComplex * data)
#else
execute_plan_fftc_ (cufftHandle * PLANp, int * ISIGNp, cufftDoubleComplex * data)
#endif
{
  cufftHandle plan = *PLANp;
  int ISIGN = *ISIGNp;
  
/*if (cudaDeviceSynchronize() != cudaSuccess){
  	fprintf(stderr, "Cuda error: Failed to synchronize\n");
  	return;	
}*/
  
  if (ISIGN== -1)
    {
#ifdef TRANS_SINGLE
    cufftSafeCall(cufftExecR2C(plan, (cufftReal*)data, data));
#else
    cufftSafeCall(cufftExecD2Z(plan, (cufftDoubleReal*)data, data));
#endif
    }
  else if (ISIGN== 1)
    {
#ifdef TRANS_SINGLE
    cufftSafeCall(cufftExecC2R(plan, data, (cufftReal*)data));
#else
    cufftSafeCall(cufftExecZ2D(plan, data, (cufftDoubleReal*)data));
#endif
    }
  else 
    {
      abort();
    }

// cudaDeviceSynchronize();

//if (cudaDeviceSynchronize() != cudaSuccess){
//	fprintf(stderr, "Cuda error: Failed to synchronize\n");
//	return;	
//}


}

