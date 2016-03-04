/* Author Dunzhu Li  dli@caltech.edu 
*/
#include"slide.h"
#include"../gpu.h"
#include<fstream>
#include"../param/param.h"

SLIDE::SLIDE(PARAM & param)
{
  int nx=param.nx;
  int nz=param.nz;
  float before=param.before;
  float after =param.after ;
  float dt    =param.dt;
  char*timefile=param.timefile;

  usepunch=param.usepunch;
  current_step=-1;
  nxblock=nx/BDIMX;
  nzblock=nz/BDIMY;
  first_arrive_step=NULL;
  last_step=NULL;

  if(usepunch!=1)   return;

  std::ifstream fid;
  fid.open(timefile,std::ios::binary);

  float *buf=new float[BDIMY*nx];
  for(int i=0;i<BDIMY*nx;i++) buf[i]=0.0;



  int *firstbuf=new int[nxblock*nzblock];
  int *lastbuf=new int[nxblock*nzblock];

  for(int k=0;k<nzblock;k+=1){
	//read BDIMY ROW one time

	fid.read((char*)buf,sizeof(float)*BDIMY*nx);
	for(int i=0;i<nxblock;i+=1){

	  float minvalue=1E28;
	  float maxvalue=-1E28;
	  for(int kk=0;kk<BDIMY;kk++){
		for(int ii=0;ii<BDIMX;ii++){
		  int ind=kk*nx+i*BDIMX+ii;
		  minvalue=min(minvalue,buf[ind]);
		  maxvalue=max(maxvalue,buf[ind]);
		}
	  }

	  firstbuf[k*nxblock+i]=int( (minvalue-before)/dt );
	  lastbuf[k*nxblock+i]=int( (maxvalue+after)/dt );

	}
  }

  fid.close();
  safecall(cudaMalloc((void**)&(first_arrive_step),sizeof(int)*nxblock*nzblock));
  safecall(cudaMemcpy(first_arrive_step,firstbuf, sizeof(int)*nxblock*nzblock,cudaMemcpyHostToDevice));

  safecall(cudaMalloc((void**)&(last_step),sizeof(int)*nxblock*nzblock));
  safecall(cudaMemcpy(last_step,lastbuf, sizeof(int)*nxblock*nzblock,cudaMemcpyHostToDevice));

  delete[]buf;
  delete[]firstbuf;
  delete[]lastbuf;
}
