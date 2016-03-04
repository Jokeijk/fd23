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
  BU=new double[nx*nz];\
  BW=new double[nx*nz];\
  MU=new double[nx*nz];\
 MUA=new double[nx*nz];\
 LAM=new double[nx*nz];\
 std::fill_n(BU,nx*nz,0);\
 std::fill_n(BW,nx*nz,0);\
 std::fill_n(MU,nx*nz,0);\
 std::fill_n(MUA,nx*nz,0);\
 std::fill_n(LAM,nx*nz,0);

/* Init global material */
void MATERIAL::init_for_full(PARAM & param)
{
  h=param.h;
  dt=param.dt;
  nx=param.nx;
  nz=param.nz;
  srcden=param.srcden;
  srcsvel=param.srcsvel;
  srcpvel=param.srcpvel;
  usetable=false;
  ALLOC_MAT;
}

/* Init box material from global material, set homogeous box in global material  */
void MATERIAL::init_for_box_from_full(int bw,int start_x,int start_z,MATERIAL &mat)
{
  nx=bw;
  nz=bw;
  usetable=false;

  ALLOC_MAT;

  // Set box material to srcsvel srcpvel srcden if not equal
  srcden  = mat.srcden;
  srcsvel = mat.srcsvel;
  srcpvel = mat.srcpvel;

  dt = mat.dt;
   h = mat.h;

  double mu= srcsvel*srcsvel*srcden;
  double lam= (srcpvel*srcpvel*srcden-2.0*mu);
  double b= 1.0/srcden;

  mu=mu*dt/h;
  lam=lam*dt/h;
  b=b*dt/h;

  // check if source region is homogeous, if otherwise, set that to homogeous source
  double eps=1E-4;
  for(int j=0;j<nz;j++){
	for(int i=0;i<nx;i++){
	  if(
		  abs(b-mat.BU[ (i+start_x)+(j+start_z)*mat.nx ])>b*eps ||
		  abs(b-mat.BW[ (i+start_x)+(j+start_z)*mat.nx ])>b*eps ||
		  abs(mu-mat.MU[ (i+start_x)+(j+start_z)*mat.nx ])>mu*eps ||
		  abs(mu-mat.MUA[ (i+start_x)+(j+start_z)*mat.nx ])>mu*eps ||
		  abs(lam-mat.LAM[ (i+start_x)+(j+start_z)*mat.nx ])>lam*eps 
		)
	  {
		double tmu=mat.MU[ (i+start_x)+(j+start_z)*mat.nx ]*h/dt;
		double tlam=mat.LAM[ (i+start_x)+(j+start_z)*mat.nx ]*h/dt;
		double tbu=mat.BU[ (i+start_x)+(j+start_z)*mat.nx ]*h/dt;
		double tden=1.0/tbu;
		double tvs=sqrt(tmu*tbu);
		double tvp=sqrt((tlam+2.0*tmu)*tbu);
		fprintf(stderr,"No homogenous region, %d,%d,%f,%f,%f,%f,%f,%f\n",i+start_x,j+start_z,
			srcsvel, tvs,           
			srcpvel, tvp,           
			srcden,tden 
		   	);
		mat.BU[ (i+start_x)+(j+start_z)*mat.nx ] =b;
		mat.BW[ (i+start_x)+(j+start_z)*mat.nx ] =b;
		mat.MU[ (i+start_x)+(j+start_z)*mat.nx ] =mu;
		mat.MUA[ (i+start_x)+(j+start_z)*mat.nx ]=mu;
		mat.LAM[ (i+start_x)+(j+start_z)*mat.nx ]=lam;
	  }
	  BU[i+j*nx]=mat.BU[ (i+start_x)+(j+start_z)*mat.nx ];
	  BW[i+j*nx]=mat.BW[ (i+start_x)+(j+start_z)*mat.nx ];
	  MU[i+j*nx]=mat.MU[ (i+start_x)+(j+start_z)*mat.nx ];
	  MUA[i+j*nx]=mat.MUA[ (i+start_x)+(j+start_z)*mat.nx ];
	  LAM[i+j*nx]=mat.LAM[ (i+start_x)+(j+start_z)*mat.nx ];
	}
  }
}

void MATERIAL::get_model(PARAM & param)
{
  char * modelname=param.modelname;

  int fd, ix, iz, k;
  double *den, *vp, *vs, fac, mu, lam;
  double test, sqrt3;
  double vpmin, vpmax, vsmin, vsmax, stab, coefsum, pts_per_wave;

  /* temporairily use state space */
  vp =new double[nx*nz];
  vs =new double[nx*nz];
  den =new double[nx*nz];

  if(extension(modelname,"mdl"))
  {
	/* Model is specified by object-based model */
	load_model(modelname);
	get_elas_model(vp,vs,den,nx,nz,h,0.0,0.0);

	if(param.output_material==1){
	  std::ofstream fid;
	  char filename[256];
	  sprintf(filename,"%s_VP",modelname);
	  fid.open(filename,std::ios::binary);
	  fid.write((char*)&nx,sizeof(int));
	  fid.write((char*)&nz,sizeof(int));
	  fid.write((char*)vp,nx*nz*sizeof(double));
	  fid.close();

	  sprintf(filename,"%s_VS",modelname);
	  fid.open(filename,std::ios::binary);
	  fid.write((char*)&nx,sizeof(int));
	  fid.write((char*)&nz,sizeof(int));
	  fid.write((char*)vs,nx*nz*sizeof(double));
	  fid.close();

	  sprintf(filename,"%s_DEN",modelname);
	  fid.open(filename,std::ios::binary);
	  fid.write((char*)&nx,sizeof(int));
	  fid.write((char*)&nz,sizeof(int));
	  fid.write((char*)den,nx*nz*sizeof(double));
	  fid.close();
	}
  }
  else
  {
	/* Model is specified by raster grids */
	read_model( vp,nx*nz,modelname,"vp");
	read_model( vs,nx*nz,modelname,"vs");
	read_model(den,nx*nz,modelname,"den");
  }
  sqrt3= sqrt(3.0);
  fac= dt/h;

  for(iz=0; iz<nz; iz++)
	for(ix=0; ix<nx; ix++)
	{
	  k= iz*nx + ix;
	  mu= vs[k]*vs[k]*den[k];
	  lam= vp[k]*vp[k]*den[k] - 2.0*mu;
	  LAM(ix,iz)= lam;
	  MU(ix,iz)= mu;
	}

  int wst_width=5;
  //debug,make left and right and bottom boudary stable for water
  for(iz=0; iz<nz; iz++){
	for(ix=0; ix<wst_width; ix++)
	{
	  k= iz*nx + ix;
	  if(vs[k]/vp[k]<0.01){
		mu= vp[k]*vp[k]/3.0*den[k];
		lam= vp[k]*vp[k]*den[k] - 2.0*mu;
		LAM(ix,iz)= lam;
		MU(ix,iz)= mu;
	  }
	}
	for(ix=nx-wst_width; ix<nx; ix++)
	{
	  k= iz*nx + ix;
	  if(vs[k]/vp[k]<0.01){
		mu= vp[k]*vp[k]/3.0*den[k];
		lam= vp[k]*vp[k]*den[k] - 2.0*mu;
		LAM(ix,iz)= lam;
		MU(ix,iz)= mu;
	  }
	}
  }
  for(iz=nz-wst_width; iz<nz; iz++){
	for(ix=0; ix<nx; ix++)
	{
	  k= iz*nx + ix;
	  if(vs[k]/vp[k]<0.01){
		mu= vp[k]*vp[k]/3.0*den[k];
		lam= vp[k]*vp[k]*den[k] - 2.0*mu;
		LAM(ix,iz)= lam;
		MU(ix,iz)= mu;
	  }
	}
  }
  //end

  for(iz=1; iz<nz; iz++)
	for(ix=0; ix<nx; ix++)
	{
	  k= iz*nx + ix;
	  if(k+1<nx*nz){
		BU(ix,iz)= 2.0/(den[k] + den[k+1]);
	  }else{
		BU(ix,iz)= 2.0/(den[k] + den[k  ]);
	  }
	  if(k+nx<nx*nz){
		BW(ix,iz)= 2.0/(den[k] + den[k-nx]);
	  }else{
		BW(ix,iz)= 2.0/(den[k] + den[k+ 0]);
	  }
	  if(ix+1<nx && iz+1<nz){
		double m1,m2,m3,m4;
		m1=MU(ix,iz);
		m2=MU(ix+1,iz);
		m3=MU(ix,iz-1);
		m4=MU(ix+1,iz-1);
		MUA(ix,iz)=double(4.0*m1*m2*m3*m4/(std::max(m2*m3*m4+m1*m3*m4+m1*m2*m4+m2*m3*m4,1E-10)));
		/*
		MUA(ix,iz)= (MU(ix,iz)+MU(ix+1,iz)+MU(ix,iz-1)+MU(ix+1,iz-1))/4.0;
		test= MU(ix,iz)*MU(ix+1,iz)*MU(ix,iz-1)*MU(ix+1,iz-1);
		if(fabs(test) < 1.0e-4) MUA(ix,iz)= 0.0;
		*/
	  }else{
		MUA(ix,iz)=MU(ix,iz);
	  }
	}
  /*
  for(iz=1; iz<nz; iz++)
	for(ix=0; ix<nx; ix++)
	{
	  k= iz*nx + ix;
	  if(ix+1<nx){
		BU(ix,iz)= 2.0/(den[k] + den[k+1]);
	  }else{
		BU(ix,iz)= 2.0/(den[k] + den[k  ]);
	  }

	  BW(ix,iz)= 2.0/(den[k] + den[k-nx]);

	  if(ix+1<nx){
		MUA(ix,iz)= (MU(ix,iz)+MU(ix+1,iz)+MU(ix,iz-1)+MU(ix+1,iz-1))/4.0;
		test= MU(ix,iz)*MU(ix+1,iz)*MU(ix,iz-1)*MU(ix+1,iz-1);
		if(fabs(test) < 1.0e-4) MUA(ix,iz)= 0.0;
	  }else{
		MUA(ix,iz)=MU(ix,iz);
	  }
	}
  */

  /* on the top row we implement a free-surface as described in
	 Mittet, R., Free-Surface Boundary Conditions for Elastic
	 Staggered-Grid modeling schemes, Geophysics, 67, 5, p1616. */
  iz=0;
  //debug,make left and right boudary stable for water
  int wst_width2=1;
  for(ix= wst_width2; ix<nx-wst_width2; ix++)
  //for(ix=0; ix<nx; ix++)
  {
	if(ix+1<nx){
	  BU(ix,iz) = 4.0/(den[ix] + den[ix+1]);
	}else{
	  BU(ix,iz) = 4.0/(den[ix] + den[ix+0]);
	}
	LAM(ix,iz)= 0.0;
	MU(ix,iz)= 0.5 * MU(ix,iz);
  }

  for(iz=0; iz<nz; iz++)
	for(ix=0; ix<nx; ix++){
	  MU(ix,iz) *= fac;
	  LAM(ix,iz) *= fac;
	  BU(ix,iz) *= fac;
	  BW(ix,iz) *= fac;
	  MUA(ix,iz) *= fac;
	}

  /* check stability condition */
  vpmax= vsmax= -99999.0;
  vpmin= vsmin=  99999.0;
  for(k=0; k<nx*nz; k++)
  {
	if(vp[k] > vpmax) vpmax= vp[k];
	if(vp[k] < vpmin) vpmin= vp[k];
	if(vs[k] > vsmax) vsmax= vs[k];
	/* extra test to discount water */
	if(vs[k] < vsmin && vs[k] > 1.0) vsmin= vs[k];
  }

  coefsum= 1.0;
  int maxord=8;
  if(maxord == 4) coefsum= C1+C2;
  if(maxord == 6) coefsum= D1+D2+D3;
  if(maxord == 8) coefsum= E1+E2+E3+E4;
  stab= vpmax * dt * coefsum * sqrt(2.0)/h;
  fprintf(stdout,"model stability vmax= %8.4f stab= %8.4f (should be < 1)\n",
	  vpmax, stab);
  if(stab>1){
	exit(-1);
  }

  delete[]vp;
  delete[]vs;
  delete[]den;
}
