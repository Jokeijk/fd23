#include <stdio.h>
const int N = 16; 
const int blocksize = 16; 

void test2();

__constant__ char d_coef[5][4];
char g_coef[5][4];

void set(){ 
  g_coef[0][1]='k';
  cudaMemcpyToSymbol(d_coef,g_coef,sizeof(char)*20);
}

__global__ 
void hello(char *a, int *b) 
{
    a[threadIdx.x] = d_coef[0][1];
}
int test1()
{
    char a[N] = "Hello \0\0\0\0\0\0";
    int b[N] = {15, 10, 6, 0, -11, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    char *ad;
    int *bd;
    const int csize = N*sizeof(char);
    const int isize = N*sizeof(int);
    printf("%s", a);
    cudaMalloc( (void**)&ad, csize ); 
    cudaMalloc( (void**)&bd, isize ); 
    cudaMemcpy( ad, a, csize, cudaMemcpyHostToDevice ); 
    cudaMemcpy( bd, b, isize, cudaMemcpyHostToDevice ); 
    
    dim3 dimBlock( blocksize, 1 );
    dim3 dimGrid( 1, 1 );
    hello<<<dimGrid, dimBlock>>>(ad, bd);
    cudaMemcpy( a, ad, csize, cudaMemcpyDeviceToHost ); 
    cudaFree( ad );
    cudaFree( bd );
    
    printf("%s\n", a);
    return EXIT_SUCCESS;
}

int main()
{
  set();
  test1();
  test2();
}
