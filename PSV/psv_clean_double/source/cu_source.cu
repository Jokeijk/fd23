#include"../gpu.h"
__global__ void convolve(double *f1,double *f2,double* out)
{
  //n1 is short, n2 is long
  //output size is n2
  extern  __shared__ double s[];  //size is 2*n1
  int k=blockDim.x*(blockIdx.x)+threadIdx.x;
  s[threadIdx.x]=f2[(int)(k-blockDim.x)];
  s[threadIdx.x+blockDim.x]=f2[k];
  __syncthreads();

  double sum=0;
  for(int i=0;i<blockDim.x;i++){
	sum += f1[i]*s[threadIdx.x+i+1];
  }
  out[k]=sum;
}
void cu_con(int N1,double *x,int N2,double *f2)
{
  static bool firstrun=true;
  static double *d_f1,*d_f2,*d_out;

  int n1=((N1-1)/32+1)*32;
  int n2=((N2-1)/n1+1)*n1;

  if(firstrun){
	firstrun=false;
	double * f1=new double[n1];
	for(int i=0;i<n1;i++) f1[i]=0.0;
	for(int i=0;i<N1;i++) f1[n1-1-i]=x[i];

	safecall(cudaMalloc((void**)&(d_f1),sizeof(double)*(n1)));
	safecall(cudaMemcpy( d_f1 , f1 , sizeof(double)*n1,cudaMemcpyHostToDevice));

	safecall(cudaMalloc((void**)&(d_f2),sizeof(double)*(n2+n1)));
	safecall(cudaMemset(d_f2, 0,sizeof(double)*(n2+n1)));
	d_f2+=n1;

	safecall(cudaMalloc((void**)&(d_out),sizeof(double)*n2));
  }
  dim3 nblocks((N2-1)/n1+1,1);
  dim3 blocksize(n1,1);
  safecall(cudaMemcpy( d_f2 , f2 , sizeof(double)*N2,cudaMemcpyHostToDevice));
  convolve<<<nblocks,blocksize,2*n1*sizeof(double)>>>(d_f1,d_f2,d_out);
  safecall(cudaMemcpy( f2 , d_out ,sizeof(double)*N2,cudaMemcpyDeviceToHost));
}
void cu_con(double *f2,double *x,int N2,int N1,double *y)
{
  static bool firstrun=true;
  static double *d_f1,*d_f2,*d_out;

  int n1=((N1-1)/32+1)*32;
  int n2=((N2-1)/n1+1)*n1;

  if(firstrun){
	firstrun=false;
	double * f1=new double[n1];
	for(int i=0;i<n1;i++) f1[i]=0.0;
	for(int i=0;i<N1;i++) f1[n1-1-i]=x[i];

	safecall(cudaMalloc((void**)&(d_f1),sizeof(double)*(n1)));
	safecall(cudaMemcpy( d_f1 , f1 , sizeof(double)*n1,cudaMemcpyHostToDevice));

	safecall(cudaMalloc((void**)&(d_f2),sizeof(double)*(n2+n1)));
	safecall(cudaMemset(d_f2, 0,sizeof(double)*(n2+n1)));
	d_f2+=n1;

	safecall(cudaMalloc((void**)&(d_out),sizeof(double)*n2));
  }
  dim3 nblocks((N2-1)/n1+1,1);
  dim3 blocksize(n1,1);
  safecall(cudaMemcpy( d_f2 , f2 , sizeof(double)*N2,cudaMemcpyHostToDevice));
  convolve<<<nblocks,blocksize,2*n1*sizeof(double)>>>(d_f1,d_f2,d_out);
  safecall(cudaMemcpy( f2 , d_out ,sizeof(double)*N2,cudaMemcpyDeviceToHost));
}

/*
int main()
{
  int N1=10;
  int N2=200;
  int n1=((N1-1)/32+1)*32;
  int n2=((N2-1)/n1+1)*n1;
  double *f1, *f2,*out;

  f1=new double[n1];
  for(int i=0;i<n1;i++) f1[i]=0.0;
  for(int i=0;i<N1;i++) f1[i]=(double)(i+1);

  f2=new double[n2];
  for(int i=0;i<n2;i++) f2[i]=0.0;
  for(int i=0;i<N2;i++) f2[i]=(double)(i+1);

  cu_con(N1,f1,N2,f2);
  for(int i=0;i<10;i++){
	printf("%f\n",f2[i]);
  }
  return;
  
  double *d_f1, *d_f2,*d_out;

  safecall(cudaMalloc((void**)&(d_f1),sizeof(double)*(n1)));
  safecall(cudaMemcpy( d_f1 , f1 , sizeof(double)*n1,cudaMemcpyHostToDevice));

  safecall(cudaMalloc((void**)&(d_f2),sizeof(double)*(n2+n1)));
  safecall(cudaMemset(d_f2, 0,sizeof(double)*(n2+n1)));
  d_f2+=n1;
  safecall(cudaMemcpy( d_f2 , f2 , sizeof(double)*n2,cudaMemcpyHostToDevice));

  safecall(cudaMalloc((void**)&(d_out),sizeof(double)*n2));

  dim3 nblocks((N2-1)/n1+1,1);
  dim3 blocksize(n1,1);

  convolve<<<nblocks,blocksize,2*n1*sizeof(double)>>>(d_f1,d_f2,d_out);

  out=new double[n2];
  safecall(cudaMemcpy( out , d_out , sizeof(double)*n2,cudaMemcpyDeviceToHost));

  for(int i=0;i<10;i++){
	printf("%f\n",out[i]);
  }
  
  for(int i=0;i<n1;i++) f1[i]=0.0;
  for(int i=0;i<N1;i++) f1[i]=(double)(i+1);
  cu_con(f2,f1,N2,N1,out);
  for(int i=0;i<10;i++){
	printf("%f\n",f2[i]);
  }

  CUT_CHECK_ERROR("error");
}
*/
