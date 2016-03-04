#include"point_source.h"
#include"../param/param.h"
#include"../field/field.h"
#include"../gpu.h"
#include"../param/const.h"
void cu_point_source(int deviceid,int it,int lsrc,float *src, PARAM &param, FIELD & fld)
{
  float temp;
  int nx=param.nx;
  int nz1=param.a_nz1[deviceid];
  int nz2=param.a_nz2[deviceid];
  int zs=param.zs;
  int xs=param.xs;


  double delta=param.dip*PI/180.0;
  double lambda=param.rake*PI/180.0;
  double phi=(param.strike-param.azimuth)*PI/180.0;

  float Mxy = sin(delta)*cos(lambda)*cos(2.0*phi)
   	         +0.5*sin(2.0*delta)*sin(lambda)*sin(2.0*phi);

  float Mzy = -(  cos(delta)*cos(lambda)*sin(phi) 
	            - cos(2.0*delta)*sin(lambda)*cos(phi) );


  int ixs,izs;

  float factor=1.0/PI/1.4142135623731/(param.h*param.h);
  Mxy *= factor;
  Mzy *= factor;

  /* Note stress has been normalized by dt,h and so on */


  /* Mxy */
  ixs=xs;  izs=zs;
  if(izs>=nz1 && izs<=nz2 && it<lsrc ){ 
	safecall(cudaMemcpy(&temp,fld.Tx+izs*nx+ixs,sizeof(float),cudaMemcpyDeviceToHost));
	temp -= 0.5*src[it]*Mxy;
	safecall(cudaMemcpy(fld.Tx+izs*nx+ixs,&temp,sizeof(float),cudaMemcpyHostToDevice));
  }

  ixs=xs-1;  izs=zs;
  if(izs>=nz1 && izs<=nz2 && it<lsrc ){ 
	safecall(cudaMemcpy(&temp,fld.Tx+izs*nx+ixs,sizeof(float),cudaMemcpyDeviceToHost));
	temp -= 0.5*src[it]*Mxy;
	safecall(cudaMemcpy(fld.Tx+izs*nx+ixs,&temp,sizeof(float),cudaMemcpyHostToDevice));
  }

  /* Mzy */
  ixs=xs;  izs=zs;
  if(izs>=nz1 && izs<=nz2 && it<lsrc ){ 
	safecall(cudaMemcpy(&temp,fld.Tz+izs*nx+ixs,sizeof(float),cudaMemcpyDeviceToHost));
	temp -= 0.5*src[it]*Mzy;
	safecall(cudaMemcpy(fld.Tz+izs*nx+ixs,&temp,sizeof(float),cudaMemcpyHostToDevice));
  }

  ixs=xs;  izs=zs+1;
  if(izs>=nz1 && izs<=nz2 && it<lsrc ){ 
	safecall(cudaMemcpy(&temp,fld.Tz+izs*nx+ixs,sizeof(float),cudaMemcpyDeviceToHost));
	temp -= 0.5*src[it]*Mzy;
	safecall(cudaMemcpy(fld.Tz+izs*nx+ixs,&temp,sizeof(float),cudaMemcpyHostToDevice));
  }

}

void cu_point_source_p(int deviceid,int it,int lsrc,float *src, PARAM &param, FIELD & fld)
{
  float temp;
  int nx=param.nx;
  int nz1=param.a_nz1[deviceid];
  int nz2=param.a_nz2[deviceid];
  int zs=param.zs;
  int xs=param.xs;


  float delta=param.dip*PI/180.0;
  float lambda=param.rake*PI/180.0;
  float phi=(param.strike-param.azimuth)*PI/180.0;

  float Mxy = sin(delta)*cos(lambda)*cos(2.0*phi)
   	         +0.5*sin(2.0*delta)*sin(lambda)*sin(2.0*phi);

  float Mzy = -(  cos(delta)*cos(lambda)*sin(phi) 
	            - cos(2.0*delta)*sin(lambda)*cos(phi) );


  int ixs,izs;

  float factor=1.0/PI/1.4142135623731/(param.h*param.h)*param.dt/param.h;
  Mxy *= factor;
  Mzy *= factor;

  if(it>=lsrc){
	it=lsrc-1;
  }


  /* Mxy */
  ixs=xs;  izs=zs;
  if(izs>=nz1 && izs<=nz2 && it<lsrc ){ 
	safecall(cudaMemcpy(&temp,fld.Tx+izs*nx+ixs,sizeof(float),cudaMemcpyDeviceToHost));
	temp -= 0.5*src[it]*Mxy;
	safecall(cudaMemcpy(fld.Tx+izs*nx+ixs,&temp,sizeof(float),cudaMemcpyHostToDevice));
  }

  ixs=xs-1;  izs=zs;
  if(izs>=nz1 && izs<=nz2 && it<lsrc ){ 
	safecall(cudaMemcpy(&temp,fld.Tx+izs*nx+ixs,sizeof(float),cudaMemcpyDeviceToHost));
	temp -= 0.5*src[it]*Mxy;
	safecall(cudaMemcpy(fld.Tx+izs*nx+ixs,&temp,sizeof(float),cudaMemcpyHostToDevice));
  }

  /* Mzy */
  ixs=xs;  izs=zs;
  if(izs>=nz1 && izs<=nz2 && it<lsrc ){ 
	safecall(cudaMemcpy(&temp,fld.Tz+izs*nx+ixs,sizeof(float),cudaMemcpyDeviceToHost));
	temp -= 0.5*src[it]*Mzy;
	safecall(cudaMemcpy(fld.Tz+izs*nx+ixs,&temp,sizeof(float),cudaMemcpyHostToDevice));
  }

  ixs=xs;  izs=zs+1;
  if(izs>=nz1 && izs<=nz2 && it<lsrc ){ 
	safecall(cudaMemcpy(&temp,fld.Tz+izs*nx+ixs,sizeof(float),cudaMemcpyDeviceToHost));
	temp -= 0.5*src[it]*Mzy;
	safecall(cudaMemcpy(fld.Tz+izs*nx+ixs,&temp,sizeof(float),cudaMemcpyHostToDevice));
  }

  /* -Mxy */
  ixs=xs-1;  izs=zs;
  if(izs>=nz1 && izs<=nz2 && it<lsrc ){ 
	safecall(cudaMemcpy(&temp,fld.Tx+izs*nx+ixs,sizeof(float),cudaMemcpyDeviceToHost));
	temp += 0.5*src[it]*Mxy;
	safecall(cudaMemcpy(fld.Tx+izs*nx+ixs,&temp,sizeof(float),cudaMemcpyHostToDevice));
  }

  ixs=xs-1-1;  izs=zs;
  if(izs>=nz1 && izs<=nz2 && it<lsrc ){ 
	safecall(cudaMemcpy(&temp,fld.Tx+izs*nx+ixs,sizeof(float),cudaMemcpyDeviceToHost));
	temp += 0.5*src[it]*Mxy;
	safecall(cudaMemcpy(fld.Tx+izs*nx+ixs,&temp,sizeof(float),cudaMemcpyHostToDevice));
  }

  /* -Mzy */
  ixs=xs-1;  izs=zs;
  if(izs>=nz1 && izs<=nz2 && it<lsrc ){ 
	safecall(cudaMemcpy(&temp,fld.Tz+izs*nx+ixs,sizeof(float),cudaMemcpyDeviceToHost));
	temp += 0.5*src[it]*Mzy;
	safecall(cudaMemcpy(fld.Tz+izs*nx+ixs,&temp,sizeof(float),cudaMemcpyHostToDevice));
  }

  ixs=xs-1;  izs=zs+1;
  if(izs>=nz1 && izs<=nz2 && it<lsrc ){ 
	safecall(cudaMemcpy(&temp,fld.Tz+izs*nx+ixs,sizeof(float),cudaMemcpyDeviceToHost));
	temp += 0.5*src[it]*Mzy;
	safecall(cudaMemcpy(fld.Tz+izs*nx+ixs,&temp,sizeof(float),cudaMemcpyHostToDevice));
  }

}
