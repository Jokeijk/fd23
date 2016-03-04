#include"material.h"
#include<cmath>
#include<cstdio>
#include<cstdlib>
#include"../model/util.h"
#include"../model/model.h"
#include<algorithm>
#include<fstream>
#include"../param/const.h"
#include"../param/param.h"
using namespace std;
#define ALLOC_MAT \
  BV=new float[nx*nz];\
 MUX=new float[nx*nz];\
 MUZ=new float[nx*nz];\
 std::fill_n(BV,nx*nz,0);\
 std::fill_n(MUX,nx*nz,0);\
 std::fill_n(MUZ,nx*nz,0);

/* Init global material */
void MATERIAL::init_for_full(PARAM & param)
{
  h=param.h;
  dt=param.dt;
  nx=param.nx;
  nz=param.nz;
  usetable=false;
  ALLOC_MAT;
}

void MATERIAL::get_model(PARAM & param)
{
  char * modelname=param.modelname;

  int fd, ix, iz, k;
  float *den, *vp, *vs, fac, mu, lam;
  double vpmin, vpmax, vsmin, vsmax, stab, coefsum, pts_per_wave;

  /* temporairily use state space */
  vp =new float[nx*nz];
  vs =new float[nx*nz];
  den =new float[nx*nz];

  if(extension(modelname,"mdl"))
  {
	/* Model is specified by object-based model */
	load_model(modelname);
	get_elas_model(vp,vs,den,nx,nz,h,0.0,0.0);

	if(param.output_material==1){
	  std::ofstream fid;
	  char filename[256];

	  sprintf(filename,"%s_VS",modelname);
	  fid.open(filename,std::ios::binary);
	  fid.write((char*)&nx,sizeof(int));
	  fid.write((char*)&nz,sizeof(int));
	  fid.write((char*)vs,nx*nz*sizeof(float));
	  fid.close();

	  sprintf(filename,"%s_DEN",modelname);
	  fid.open(filename,std::ios::binary);
	  fid.write((char*)&nx,sizeof(int));
	  fid.write((char*)&nz,sizeof(int));
	  fid.write((char*)den,nx*nz*sizeof(float));
	  fid.close();
	}
  }
  else
  {
	/* Model is specified by raster grids */
	read_model( vs,nx*nz,modelname,"vs");
	read_model(den,nx*nz,modelname,"den");
  }

  for(iz=0; iz<nz; iz++){
	for(ix=0; ix<nx; ix++)
	{
	  k= iz*nx + ix;
	  BV(ix,iz)= 1.0/den[k];

	  double m1,m2;
	  m1= vs[k]*vs[k]*den[k];
	  if(k+1<nx*nz){
		m2= vs[k+1]*vs[k+1]*den[k+1];
	  }else{
		m2=m1;
	  }
	  MUX(ix,iz)=float(2.0*m1*m2/(std::max(m1+m2,1E-10)));

	  m1= vs[k]*vs[k]*den[k];
	  if(k-nx>=0){
		m2= vs[k-nx]*vs[k-nx]*den[k-nx];
	  }else{
		m2= m1;
	  }
	  MUZ(ix,iz)=float(2.0*m1*m2/(std::max(m1+m2,1E-10)));
	}
  }

  fac= dt/h;
  for(iz=0; iz<nz; iz++){
	for(ix=0; ix<nx; ix++){
	  MUX(ix,iz) *= fac;
	  MUZ(ix,iz) *= fac;
	  BV(ix,iz) *= fac;
	}
  }

  /* check stability condition */
  vpmax= vsmax= -99999.0;
  vpmin= vsmin=  99999.0;
  for(k=0; k<nx*nz; k++)
  {
	if(vs[k] > vsmax) vsmax= vs[k];
	/* extra test to discount water */
	if(vs[k] < vsmin && vs[k] > 1.0) vsmin= vs[k];
  }

  coefsum= 1.0;
  int maxord=8;
  if(maxord == 4) coefsum= C1+C2;
  if(maxord == 6) coefsum= D1+D2+D3;
  if(maxord == 8) coefsum= E1+E2+E3+E4;
  stab= vsmax * dt * coefsum * sqrt(2.0)/h;
  fprintf(stdout,"model stability vmax= %8.4f stab= %8.4f (should be < 1)\n",
	  vsmax, stab);
  if(stab>1){
	exit(-1);
  }

  delete[]vp;
  delete[]vs;
  delete[]den;
}
