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


#ifdef EXPLOSION_SOURCE
  float Mxx=1;
  float Mzz=1;
  float Mxz=0;
#else

  float delta=param.dip*PI/180.0;
  float lambda=param.rake*PI/180.0;
  float phi=(param.strike-param.azimuth)*PI/180.0;

  float Mxx = -( sin(delta)*cos(lambda)*sin(2.0*phi) + 
	  sin(2.0*delta)*sin(lambda)*sin(phi)*sin(phi) );

  float Mxz = -( cos(delta)*cos(lambda)*cos(phi) + 
	  cos(2.0*delta)*sin(lambda)*sin(phi) );

  float Mzz = sin(2.0*delta)*sin(lambda);
#endif

  int ixs,izs;

  float factor=1.0/PI/1.4142135623731/(param.h*param.h);
  Mxx *= factor;
  Mxz *= factor;
  Mzz *= factor;

  /* Note stress has been normalized by dt,h and so on */


  /* Mxx */
  ixs=xs;  izs=zs;
  if(izs>=nz1 && izs<=nz2 && it<lsrc ){ 
	safecall(cudaMemcpy(&temp,fld.Txx+izs*nx+ixs,sizeof(float),cudaMemcpyDeviceToHost));
	temp -= src[it]*Mxx;
	safecall(cudaMemcpy(fld.Txx+izs*nx+ixs,&temp,sizeof(float),cudaMemcpyHostToDevice));
  }


  /* Mzz */
  ixs=xs;   izs=zs;
  if(zs>=nz1 && zs<=nz2 && it<lsrc ){ 
	safecall(cudaMemcpy(&temp,fld.Tzz+izs*nx+ixs,sizeof(float),cudaMemcpyDeviceToHost));
	temp -= src[it]*Mzz;
	safecall(cudaMemcpy(fld.Tzz+izs*nx+ixs,&temp,sizeof(float),cudaMemcpyHostToDevice));
  }
  /* Mxz */
  for(ixs=xs-1;ixs<=xs;ixs++){
	for(izs=zs;izs<=zs+1;izs++){
	  if(izs>=nz1 && izs<=nz2 && it<lsrc ){ 
		safecall(cudaMemcpy(&temp,fld.Txz+izs*nx+ixs,sizeof(float),cudaMemcpyDeviceToHost));
		temp -= src[it]*Mxz*0.25;
		safecall(cudaMemcpy(fld.Txz+izs*nx+ixs,&temp,sizeof(float),cudaMemcpyHostToDevice));
	  }
	}
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

#ifdef EXPLOSION_SOURCE
  float Mxx=1;
  float Mzz=1;
  float Mxz=0;
#else

  float delta=param.dip*PI/180.0;
  float lambda=param.rake*PI/180.0;
  float phi=(param.strike-param.azimuth)*PI/180.0;

  float Mxx = -( sin(delta)*cos(lambda)*sin(2.0*phi) + 
	  sin(2.0*delta)*sin(lambda)*sin(phi)*sin(phi) );

  float Mxz = -( cos(delta)*cos(lambda)*cos(phi) + 
	  cos(2.0*delta)*sin(lambda)*sin(phi) );

  float Mzz = sin(2.0*delta)*sin(lambda);
#endif

  int ixs,izs;

  float factor=1.0/PI/1.4142135623731/(param.h*param.h)*param.dt/param.h;
  Mxx *= factor;
  Mxz *= factor;
  Mzz *= factor;

  /* Note stress has been normalized by dt,h and so on */

  if(it>=lsrc){
	it=lsrc-1;
  }


  /* Mxx */
  ixs=xs;  izs=zs;
  if(izs>=nz1 && izs<=nz2 && it<lsrc ){ 
	safecall(cudaMemcpy(&temp,fld.Txx+izs*nx+ixs,sizeof(float),cudaMemcpyDeviceToHost));
	temp -= src[it]*Mxx;
	safecall(cudaMemcpy(fld.Txx+izs*nx+ixs,&temp,sizeof(float),cudaMemcpyHostToDevice));
  }

  ixs=xs-1;  izs=zs;
  if(izs>=nz1 && izs<=nz2 && it<lsrc ){ 
	safecall(cudaMemcpy(&temp,fld.Txx+izs*nx+ixs,sizeof(float),cudaMemcpyDeviceToHost));
	temp += src[it]*Mxx;
	safecall(cudaMemcpy(fld.Txx+izs*nx+ixs,&temp,sizeof(float),cudaMemcpyHostToDevice));
  }


  /* Mzz */
  ixs=xs;   izs=zs;
  if(zs>=nz1 && zs<=nz2 && it<lsrc ){ 
	safecall(cudaMemcpy(&temp,fld.Tzz+izs*nx+ixs,sizeof(float),cudaMemcpyDeviceToHost));
	temp -= src[it]*Mzz;
	safecall(cudaMemcpy(fld.Tzz+izs*nx+ixs,&temp,sizeof(float),cudaMemcpyHostToDevice));
  }

  ixs=xs-1;   izs=zs;
  if(zs>=nz1 && zs<=nz2 && it<lsrc ){ 
	safecall(cudaMemcpy(&temp,fld.Tzz+izs*nx+ixs,sizeof(float),cudaMemcpyDeviceToHost));
	temp += src[it]*Mzz;
	safecall(cudaMemcpy(fld.Tzz+izs*nx+ixs,&temp,sizeof(float),cudaMemcpyHostToDevice));
  }
  /* Mxz */
  for(ixs=xs-1;ixs<=xs;ixs++){
	for(izs=zs;izs<=zs+1;izs++){
	  if(izs>=nz1 && izs<=nz2 && it<lsrc ){ 
		safecall(cudaMemcpy(&temp,fld.Txz+izs*nx+ixs,sizeof(float),cudaMemcpyDeviceToHost));
		temp -= src[it]*Mxz*0.25;
		safecall(cudaMemcpy(fld.Txz+izs*nx+ixs,&temp,sizeof(float),cudaMemcpyHostToDevice));

		safecall(cudaMemcpy(&temp,fld.Txz+izs*nx+ixs-1,sizeof(float),cudaMemcpyDeviceToHost));
		temp += src[it]*Mxz*0.25;
		safecall(cudaMemcpy(fld.Txz+izs*nx+ixs-1,&temp,sizeof(float),cudaMemcpyHostToDevice));
	  }
	}
  }
}
