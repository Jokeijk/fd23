#include <stdio.h>



__device__ __constant__ char d_coef[2];
 
__global__ 
void add_arrays(char *a, int *b) 
{
	a[threadIdx.x] = d_coef[0];
}


/*
char g_coef[2]={'a','b'};
void pre()
{
  cudaMemcpyToSymbol(d_coef,g_coef,sizeof(char)*2);
}
*/

void pre();
#define N 7
int main()
{
    pre();
	char a[N] = "Hello ";
	int b[N] = {15, 10, 6, 0, -11, 1,0};
  
	char *ad;
	int *bd;
	const int csize = N*sizeof(char);
	const int isize = N*sizeof(int);
 
	// print the contents of a[]
	//printf("%s", a);
 
	// Allocate and Transfer memory to the device
	cudaMalloc( (void**)&ad, csize ); 
	cudaMalloc( (void**)&bd, isize ); 
	
	cudaMemcpy( ad, a, csize, cudaMemcpyHostToDevice ); 
	cudaMemcpy( bd, b, isize, cudaMemcpyHostToDevice ); 
	
	// Perform the array addition
	dim3 dimBlock( N  );  
	dim3 dimGrid ( 1  );
	add_arrays<<<dimGrid, dimBlock>>>(ad, bd);
	
	// Copy the Contents from the GPU
	cudaMemcpy( a, ad, csize, cudaMemcpyDeviceToHost ); 
	cudaFree( ad );
	
	// Display the results
	printf("%s\n", a);
	return EXIT_SUCCESS;
}
