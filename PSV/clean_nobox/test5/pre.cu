extern __device__ __constant__ char d_coef[2];

char g_coef[2]={'a','b'};

void pre()
{
  cudaMemcpyToSymbol(d_coef,g_coef,sizeof(char)*2);
}
