#include"fd.h"
#include"../param/const.h"
#include"../param/param.h"
#include"../material/material.h"
#include"../field/field.h"
#include"../cpml/cpml.h"

#include<algorithm>
using namespace std;
#define topabsorb false

#define  UDIFF \
  ux=g_coef[ord][0]*( U(ix,iz)   - U(ix-1,iz) )\
    -g_coef[ord][1]*( U(ix+1,iz) - U(ix-2,iz) )\
    +g_coef[ord][2]*( U(ix+2,iz) - U(ix-3,iz) )\
    -g_coef[ord][3]*( U(ix+3,iz) - U(ix-4,iz) );\
  uz=g_coef[ord][0]*( U(ix,iz)   - U(ix,iz-1) )\
    -g_coef[ord][1]*( U(ix,iz+1) - U(ix,iz-2) )\
    +g_coef[ord][2]*( U(ix,iz+2) - U(ix,iz-3) )\
    -g_coef[ord][3]*( U(ix,iz+3) - U(ix,iz-4) );

#define WDIFF \
  wx=g_coef[ord][0]*( W(ix+1,iz)   - W(ix,iz) ) \
    -g_coef[ord][1]*( W(ix+2,iz) - W(ix-1,iz) ) \
    +g_coef[ord][2]*( W(ix+3,iz) - W(ix-2,iz) ) \
    -g_coef[ord][3]*( W(ix+4,iz) - W(ix-3,iz) );\
  wz=g_coef[ord][0]*( W(ix,iz+1) - W(ix,iz) )\
    -g_coef[ord][1]*( W(ix,iz+2) - W(ix,iz-1) )\
    +g_coef[ord][2]*( W(ix,iz+3) - W(ix,iz-2) )\
    -g_coef[ord][3]*( W(ix,iz+4) - W(ix,iz-3) );

#define  DXTXX \
  dxTxx=g_coef[ord][0]*( Txx(ix+1,iz)   - Txx(ix,iz) )\
       -g_coef[ord][1]*( Txx(ix+2,iz) - Txx(ix-1,iz) )\
       +g_coef[ord][2]*( Txx(ix+3,iz) - Txx(ix-2,iz) )\
       -g_coef[ord][3]*( Txx(ix+4,iz) - Txx(ix-3,iz) );


#define  DZTZZ \
  dzTzz=g_coef[ord][0]*( Tzz(ix,iz)   - Tzz(ix,iz-1) )\
       -g_coef[ord][1]*( Tzz(ix,iz+1) - Tzz(ix,iz-2) )\
       +g_coef[ord][2]*( Tzz(ix,iz+2) - Tzz(ix,iz-3) )\
       -g_coef[ord][3]*( Tzz(ix,iz+3) - Tzz(ix,iz-4) );

#define DXDZTXZ \
  dxTxz=g_coef[ord][0]*( Txz(ix,iz)   - Txz(ix-1,iz) ) \
       -g_coef[ord][1]*( Txz(ix+1,iz) - Txz(ix-2,iz) ) \
       +g_coef[ord][2]*( Txz(ix+2,iz) - Txz(ix-3,iz) ) \
       -g_coef[ord][3]*( Txz(ix+3,iz) - Txz(ix-4,iz) );\
  dzTxz=g_coef[ord][0]*( Txz(ix,iz+1) - Txz(ix,iz) )\
       -g_coef[ord][1]*( Txz(ix,iz+2) - Txz(ix,iz-1) )\
       +g_coef[ord][2]*( Txz(ix,iz+3) - Txz(ix,iz-2) )\
       -g_coef[ord][3]*( Txz(ix,iz+4) - Txz(ix,iz-3) );


#define U(ix,iz) U[ix+(iz)*nx]
#define W(ix,iz) W[ix+(iz)*nx]
#define Txx(ix,iz) Txx[ix+(iz)*nx]
#define Txz(ix,iz) Txz[ix+(iz)*nx]
#define Tzz(ix,iz) Tzz[ix+(iz)*nx]


void step_stress_r(FIELD & fld,MATERIAL & mat)
{
  double ux,uz,wx,wz;

  int nx=fld.nx;
  int nz=fld.nz;

  double *U=fld.U;
  double *W=fld.W;
  double *Txx=fld.Txx;
  double *Txz=fld.Txz;
  double *Tzz=fld.Tzz;

  int ord=ORD8;
  for(int iz=4;iz<nz-4;iz++){
	for(int ix=4;ix<nx-4;ix++){
	  int in_idx=ix+iz*nx;
	  double LAM_ix_iz= mat.LAM[in_idx];
	  double MU_ix_iz = mat.MU [in_idx];
	  double MUA_ix_iz= mat.MUA[in_idx];
	  double gam = LAM_ix_iz + 2.0*MU_ix_iz;
	  UDIFF;
	  WDIFF;
	  Txx[in_idx] += gam*ux + LAM_ix_iz*wz ;
	  Tzz[in_idx] += LAM_ix_iz*ux + gam*wz ;
	  Txz[in_idx] += MUA_ix_iz*( wx +uz);
	}
  }
}

void step_stress(FIELD & fld,MATERIAL & mat,CPML & cpml)
{
  int ord,idx;
  double ux,uz,wx,wz;

  int nx=fld.nx;
  int nz=fld.nz;

  int npml=cpml.npml;

  double * __restrict__ U=fld.U;
  double * __restrict__ W=fld.W;
  double * __restrict__ Txx=fld.Txx;
  double * __restrict__ Txz=fld.Txz;
  double * __restrict__ Tzz=fld.Tzz;
  double * __restrict__ LAM=mat.LAM;
  double * __restrict__ MU=mat.MU;
  double * __restrict__ MUA=mat.MUA;

  double * __restrict__ psi_U_x=cpml.psi.U_x;
  double * __restrict__ psi_W_x=cpml.psi.W_x;
  double * __restrict__ psi_U_z=cpml.psi.U_z;
  double * __restrict__ psi_W_z=cpml.psi.W_z;

  double * __restrict__ b_U_x=cpml.b.U_x;
  double * __restrict__ b_W_x=cpml.b.W_x;
  double * __restrict__ b_U_z=cpml.b.U_z;
  double * __restrict__ b_W_z=cpml.b.W_z;

  double * __restrict__ c_U_x=cpml.c.U_x;
  double * __restrict__ c_W_x=cpml.c.W_x;
  double * __restrict__ c_U_z=cpml.c.U_z;
  double * __restrict__ c_W_z=cpml.c.W_z;

  double * __restrict__ k_U_x=cpml.k.U_x;
  double * __restrict__ k_W_x=cpml.k.W_x;
  double * __restrict__ k_U_z=cpml.k.U_z;
  double * __restrict__ k_W_z=cpml.k.W_z;



  ord=ORD8;
  for(int iz=0;iz<nz;iz++){
	for(int ix=0;ix<nx;ix++){
	  int in_idx=ix+iz*nx;
	  int ord=ix;
	  ord=min(ord,nx-1-ix);
	  ord=min(ord,nz-1-iz);
	  ord=min(ord,ORD8);


	  double LAM_ix_iz= LAM[in_idx];
	  double MU_ix_iz = MU [in_idx];
	  double MUA_ix_iz= MUA[in_idx];
	  double gam = LAM_ix_iz + 2.0*MU_ix_iz;

	  UDIFF;
	  WDIFF;
	  if(ix<npml){
		//left pml
		idx=ix+iz*2*npml;
		psi_U_x[idx] = b_U_x[ix] * psi_U_x[idx] + c_U_x[ix] * ux;
		Txx[in_idx] += gam*( ux * k_U_x[ix] + psi_U_x[idx] );
		Tzz[in_idx] += LAM_ix_iz*( ux * k_U_x[ix] + psi_U_x[idx] );

		psi_W_x[idx] = b_W_x[ix] * psi_W_x[idx] + c_W_x[ix] * wx;
		Txz[in_idx] += MUA_ix_iz*( wx * k_W_x[ix] + psi_W_x[idx] );

	  }else if(ix>=nx-npml){
		//right pml
		idx=npml+nx-1-ix+iz*2*npml;
		psi_U_x[idx] = b_U_x[ix] * psi_U_x[idx] + c_U_x[ix] * ux;
		Txx[in_idx] += gam*( ux * k_U_x[ix] + psi_U_x[idx] );
		Tzz[in_idx] += LAM_ix_iz*( ux * k_U_x[ix] + psi_U_x[idx] );

		psi_W_x[idx] = b_W_x[ix] * psi_W_x[idx] + c_W_x[ix] * wx;
		Txz[in_idx] += MUA_ix_iz*( wx * k_W_x[ix] + psi_W_x[idx] );
	  }else if(iz!=0){
		Txx[in_idx] += gam*( ux );
		Tzz[in_idx] += LAM_ix_iz*( ux );
		Txz[in_idx] += MUA_ix_iz*( wx );
	  }
	  if(iz<npml && topabsorb){
		//top pml
		idx=(iz*nx)+ix;
		psi_W_z[idx] = b_W_z[iz] * psi_W_z[idx] + c_W_z[iz] * wz;
		Txx[in_idx] += LAM_ix_iz*( wz * k_W_z[iz] + psi_W_z[idx] ) ;
		Tzz[in_idx] += gam      *( wz * k_W_z[iz] + psi_W_z[idx] ) ;

		psi_U_z[idx] = b_U_z[iz] * psi_U_z[idx] + c_U_z[iz] * uz;
		Txz[in_idx] += MUA_ix_iz*( uz * k_U_z[iz] + psi_U_z[idx] ) ;
	  }else if(iz>=nz-npml){
		// bottom
		idx=(npml+nz-1-iz)*nx+ix;
		psi_W_z[idx] = b_W_z[iz] * psi_W_z[idx] + c_W_z[iz] * wz;
		Txx[in_idx] += LAM_ix_iz*( wz * k_W_z[iz] + psi_W_z[idx] ) ;
		Tzz[in_idx] += gam      *( wz * k_W_z[iz] + psi_W_z[idx] ) ;

		psi_U_z[idx] = b_U_z[iz] * psi_U_z[idx] + c_U_z[iz] * uz;
		Txz[in_idx] += MUA_ix_iz*( uz * k_U_z[iz] + psi_U_z[idx] ) ;
	  }else if(iz==0){
		// lam=0,gam=2*mu, so gam*ux+lam*wz=2*mu*ux
		Txx[in_idx]+= 2.0*MU_ix_iz*ux;
		Tzz[in_idx]= 0.0;
		Txz[in_idx]= 0.0;
	  }else{
		Txx[in_idx] += LAM_ix_iz*( wz );
		Tzz[in_idx] += gam      *( wz );
		Txz[in_idx] += MUA_ix_iz*( uz );
	  }

	}
  }

}
void step_velocity_r(FIELD & fld,MATERIAL & mat)
{

  double dxTxx,dzTzz,dxTxz,dzTxz;
  int nx=fld.nx;
  int nz=fld.nz;

  double *U=fld.U;
  double *W=fld.W;
  double *Txx=fld.Txx;
  double *Txz=fld.Txz;
  double *Tzz=fld.Tzz;

  int ord=ORD8;
  for(int iz=4;iz<nz-4;iz++){
	for(int ix=4;ix<nx-4;ix++){
	  int in_idx=ix+iz*nx;
	  double BU_ix_iz= mat.BU[in_idx];
	  double BW_ix_iz= mat.BW[in_idx];
	  DXTXX;
	  DZTZZ;
	  DXDZTXZ;
	  U[in_idx] += BU_ix_iz*(dxTxx+dzTxz);
	  W[in_idx] += BW_ix_iz*(dxTxz+dzTzz);
	}
  }
}

void step_velocity(FIELD & fld,MATERIAL & mat,CPML & cpml)
{
  int ord,idx;
  double dxTxx,dzTzz,dxTxz,dzTxz;

  int nx=fld.nx;
  int nz=fld.nz;

  int npml=cpml.npml;

  double *U=fld.U;
  double *W=fld.W;
  double *Txx=fld.Txx;
  double *Txz=fld.Txz;
  double *Tzz=fld.Tzz;

#define psi cpml.psi
#define b cpml.b
#define c cpml.c
#define k cpml.k

  ord=ORD8;
  for(int iz=0;iz<nz;iz++){
	for(int ix=0;ix<nx;ix++){
	  int in_idx=ix+iz*nx;
	  double BU_ix_iz= mat.BU[in_idx];
	  double BW_ix_iz= mat.BW[in_idx];
	  ord=ix; // maybe top layer can also use 8th ,this reduce dispersion 
	  ord=min(ord,nx-1-ix);
	  ord=min(ord,nz-1-iz);
	  ord=min(ord,ORD8);
	  DXTXX;
	  DZTZZ;
	  DXDZTXZ;
	  if(ix<npml){
		//left pml
		idx=ix+iz*2*npml;
		psi.Txx_x[idx] = b.Txx_x[ix] * psi.Txx_x[idx] + c.Txx_x[ix] * dxTxx;
		fld.U[in_idx] += BU_ix_iz*( dxTxx * k.Txx_x[ix] + psi.Txx_x[idx] );

		psi.Txz_x[idx] = b.Txz_x[ix] * psi.Txz_x[idx] + c.Txz_x[ix] * dxTxz;
		fld.W[in_idx] += BW_ix_iz*( dxTxz * k.Txz_x[ix] + psi.Txz_x[idx] );
	  }else if(ix>=nx-npml){
		//right pml
		idx=npml+nx-1-ix+iz*2*npml;
		psi.Txx_x[idx] = b.Txx_x[ix] * psi.Txx_x[idx] + c.Txx_x[ix] * dxTxx;
		fld.U[in_idx] += BU_ix_iz*( dxTxx * k.Txx_x[ix] + psi.Txx_x[idx] );

		psi.Txz_x[idx] = b.Txz_x[ix] * psi.Txz_x[idx] + c.Txz_x[ix] * dxTxz;
		fld.W[in_idx] += BW_ix_iz*( dxTxz * k.Txz_x[ix] + psi.Txz_x[idx] );
	  }else if(iz!=0){
		fld.U[in_idx] += BU_ix_iz*(dxTxx);
		fld.W[in_idx] += BW_ix_iz*(dxTxz);
	  } 
	  if(iz<npml && topabsorb){
		//top pml
		idx=(iz*nx)+ix;
		psi.Txz_z[idx] = b.Txz_z[iz] * psi.Txz_z[idx] + c.Txz_z[iz] * dzTxz;
		psi.Tzz_z[idx] = b.Tzz_z[iz] * psi.Tzz_z[idx] + c.Tzz_z[iz] * dzTzz;

		fld.U[in_idx] += BU_ix_iz*(dzTxz * k.Txz_z[iz] + psi.Txz_z[idx] );
		fld.W[in_idx] += BW_ix_iz*(dzTzz * k.Tzz_z[iz] + psi.Tzz_z[idx] );

	  }else if(iz>=nz-npml){
		// bottom
		idx=(npml+nz-1-iz)*nx+ix;
		psi.Txz_z[idx] = b.Txz_z[iz] * psi.Txz_z[idx] + c.Txz_z[iz] * dzTxz;
		psi.Tzz_z[idx] = b.Tzz_z[iz] * psi.Tzz_z[idx] + c.Tzz_z[iz] * dzTzz;

		fld.U[in_idx] += BU_ix_iz*(dzTxz * k.Txz_z[iz] + psi.Txz_z[idx] );
		fld.W[in_idx] += BW_ix_iz*(dzTzz * k.Tzz_z[iz] + psi.Tzz_z[idx] );
	  }else if(iz==0){
		fld.U[in_idx] += BU_ix_iz*(dxTxx+dzTxz);
		fld.W[in_idx]=0.0;
	  }else{
		fld.U[in_idx] += BU_ix_iz*(dzTxz);
		fld.W[in_idx] += BW_ix_iz*(dzTzz);
	  }
	}
  }
}
