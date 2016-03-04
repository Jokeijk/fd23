#include"record.h"
#include"../gpu.h"
#include"../field/field.h"
#include"../material/material.h"
#include"../param/param.h"
#include<fstream>
#include<iostream>
#include<iterator>
#include<cstdio>
#include<cstdlib>

void RECORD::init_line(int total_rec, int ix0, int iz0, int idx, int idz, const char * output)
{
	nrec=total_rec*2;
	rx=new int[nrec];
	rz=new int[nrec];
	g_rec=new float[nrec];
	for(int i=0;i<nrec/2;i++){
	  rx[2*i]=ix0+idx*i;
	  rz[2*i]=iz0+idz*i;

	  rx[2*i+1]=rx[2*i]-1;
	  rz[2*i+1]=rz[2*i];
	}
	sprintf(recordfile,"%s",output);
	fid.open(recordfile,std::ios::binary);

}

void RECORD::init_list(const char *sta_file,const char * output){
	int Tx[100000];
	int Tz[100000];
	nrec=0;
	std::ifstream file(sta_file);
	if(!file.is_open())
	{
	  fprintf(stderr,"%s Not exist!",sta_file);
	  exit(-1);
	}

	std::istream_iterator<int> eos;
	std::istream_iterator<int> iit(file);

	while (true) { 
	  Tz[nrec]=*iit;
	  iit++;
	  Tx[nrec]=*iit;
	  iit++;

	  printf("nrec is %s %d %d %d \n",sta_file,nrec,Tz[nrec], Tx[nrec]);
	  nrec += 1;

	  if(iit==eos) break;

	  if(nrec>=100000){
		fprintf(stderr,"Too much record.\n");
		exit(-1);
	  }
	}
	file.close();

	nrec=nrec*2;

	g_rec=new float[nrec];
	rx=new int[nrec];
	rz=new int[nrec];

	for(int i=0;i<nrec/2;i++){
	  rx[2*i]=Tx[i];
	  rz[2*i]=Tz[i];

	  rx[2*i+1]=rx[2*i]-1;
	  rz[2*i+1]=rz[2*i];
	}
	sprintf(recordfile,"%s",output);
	fid.open(recordfile,std::ios::binary);
}

void RECORD::cu_init_record(int k,PARAM &param)
{
  int* trx=new int[nrec];
  int* trz=new int[nrec];
  int* tind_rec=new int[nrec];

  int ind=0;
  for (int i=0;i<nrec;i++){
	if(rz[i] >= param.a_nz1[k] && rz[i] <= param.a_nz2[k] )
	{
	  trx[ind]=rx[i];
	  trz[ind]=rz[i];
	  tind_rec[ind]=i;
	  ind=ind+1;
	}
  }
  lrec[k]=ind;
  safecall(cudaMalloc((void**)&d_rx[k],sizeof(int)*lrec[k]));
  safecall(cudaMalloc((void**)&d_rz[k],sizeof(int)*lrec[k]));
  safecall(cudaMemcpy(d_rx[k],trx,sizeof(int)*lrec[k],cudaMemcpyHostToDevice));
  safecall(cudaMemcpy(d_rz[k],trz,sizeof(int)*lrec[k],cudaMemcpyHostToDevice));
  ind_rec[k]=new int[lrec[k]];
  for(int i1=0;i1<lrec[k];i1++){
	ind_rec[k][i1]=tind_rec[i1];
  }
  safecall(cudaMalloc((void**)&d_recbuf[k],sizeof(float)*lrec[k]));
  g_d_recbuf[k]=new float[lrec[k]];

  delete[] trx;
  delete[] trz;
  delete[] tind_rec;
}

void RECORD::close()
{
  fid.close();
}



int cpy_backhost(int deviceid,PARAM & param,FIELD & fld,FIELD & d_fld)
{
  int nz1,nz2,tnz;
  int nx=fld.nx;
  nz1=param.a_nz1[deviceid];
  nz2=param.a_nz2[deviceid];
  tnz=nz2-nz1+1;

  safecall(cudaMemcpy(fld.Txx+nz1*nx,d_fld.Txx+nz1*nx,sizeof(float)*nx*tnz,cudaMemcpyDeviceToHost));
  safecall(cudaMemcpy(fld.Txz+nz1*nx,d_fld.Txz+nz1*nx,sizeof(float)*nx*tnz,cudaMemcpyDeviceToHost));
  safecall(cudaMemcpy(fld.Tzz+nz1*nx,d_fld.Tzz+nz1*nx,sizeof(float)*nx*tnz,cudaMemcpyDeviceToHost));
  safecall(cudaMemcpy(fld.U  +nz1*nx,d_fld.U  +nz1*nx,sizeof(float)*nx*tnz,cudaMemcpyDeviceToHost));
  safecall(cudaMemcpy(fld.W  +nz1*nx,d_fld.W  +nz1*nx,sizeof(float)*nx*tnz,cudaMemcpyDeviceToHost));

  CUT_CHECK_ERROR("cpy_backhost");
  return 0;
}

/* record */
__global__ void cudarecord_list(int nx,float* d_U,
	int nrec,int *rx, int*rz,float *rec)
{
  //record according to index
  int ir=blockIdx.x*blockDim.x+threadIdx.x;
  if (ir <nrec){
	rec[ir]=d_U[rx[ir]+rz[ir]*nx];
  }
}

void RECORD::cu_flush(int deviceid){
  if(deviceid==0){
	fid.write((char*)g_rec,sizeof(float)*nrec);
  }
}

void RECORD::cu_record(int deviceid,int nx,float* U)
{
  int myrec=lrec[deviceid];

  dim3 nblocks(myrec/256+1,1);
  dim3 blocksize(256,1);

  cudarecord_list<<<nblocks,blocksize>>>(nx,U,
	  myrec,d_rx[deviceid],d_rz[deviceid],d_recbuf[deviceid]);

  cudaMemcpy(g_d_recbuf[deviceid],d_recbuf[deviceid],sizeof(float)*myrec,cudaMemcpyDeviceToHost);
  for(int i=0;i<myrec;i++){
	int g_ind=ind_rec[deviceid][i];
	g_rec[g_ind] = g_d_recbuf[deviceid][i];
  }
  CUT_CHECK_ERROR("record");
}
