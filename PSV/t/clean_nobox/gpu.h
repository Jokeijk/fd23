/* Author Dunzhu Li  dli@caltech.edu 
*/
#include<cuda.h>
#include<cstdio>
#define safecall(call) do{\
  cudaError_t err = call ;\
  if (cudaSuccess != err){\
	fprintf(stderr, "cuda error at %s:%d, %s\n",\
		__FILE__, __LINE__, cudaGetErrorString(err));\
  }\
}while(0)

#define CUT_CHECK_ERROR(errorMessage) do {                           \
  cudaThreadSynchronize();                                           \
  cudaError_t err = cudaGetLastError();                              \
  if( cudaSuccess != err) {                                          \
	fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",\
		errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
	exit(EXIT_FAILURE);                                              \
  }} while (0)
#ifdef USE_RESTRICT
#define userestrict __restrict__
#else
#define userestrict
#endif
