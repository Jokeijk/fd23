#include"fd.h"
#include"../param/const.h"
#include"../param/param.h"
#include"../material/material.h"
#include"../field/field.h"
#include"../cpml/cpml.h"

#include<algorithm>
using namespace std;
#define topabsorb false

#define  TxDIFF \
  Tx_x=g_coef[ord][0]*( Tx(ix,iz)   - Tx(ix-1,iz) )\
    -g_coef[ord][1]*( Tx(ix+1,iz) - Tx(ix-2,iz) )\
    +g_coef[ord][2]*( Tx(ix+2,iz) - Tx(ix-3,iz) )\
    -g_coef[ord][3]*( Tx(ix+3,iz) - Tx(ix-4,iz) );

#define TzDIFF \
  Tz_z=g_coef[ord][0]*( Tz(ix,iz+1) - Tz(ix,iz) )\
    -g_coef[ord][1]*( Tz(ix,iz+2) - Tz(ix,iz-1) )\
    +g_coef[ord][2]*( Tz(ix,iz+3) - Tz(ix,iz-2) )\
    -g_coef[ord][3]*( Tz(ix,iz+4) - Tz(ix,iz-3) );

#define  DxV \
  dxV=g_coef[ord][0]*( V(ix+1,iz)   - V(ix,iz) )\
       -g_coef[ord][1]*( V(ix+2,iz) - V(ix-1,iz) )\
       +g_coef[ord][2]*( V(ix+3,iz) - V(ix-2,iz) )\
       -g_coef[ord][3]*( V(ix+4,iz) - V(ix-3,iz) );


#define  DzV \
  dzV=g_coef[ord][0]*( V(ix,iz)   - V(ix,iz-1) )\
       -g_coef[ord][1]*( V(ix,iz+1) - V(ix,iz-2) )\
       +g_coef[ord][2]*( V(ix,iz+2) - V(ix,iz-3) )\
       -g_coef[ord][3]*( V(ix,iz+3) - V(ix,iz-4) );


#define Tx(ix,iz) Tx[ix+(iz)*nx]
#define Tz(ix,iz) Tz[ix+(iz)*nx]
#define V(ix,iz) V[ix+(iz)*nx]


void step_stress(FIELD & fld,MATERIAL & mat,CPML & cpml)
{
  int ord,idx;
  float dxV,dzV;

  int nx=fld.nx;
  int nz=fld.nz;

  int npml=cpml.npml;

  float * __restrict__ Tx=fld.Tx;
  float * __restrict__ Tz=fld.Tz;
  float * __restrict__ V=fld.V;

  float * __restrict__ MUX=mat.MUX;
  float * __restrict__ MUZ=mat.MUZ;

  float * __restrict__ psi_V_x=cpml.psi.V_x;
  float * __restrict__ psi_V_z=cpml.psi.V_z;

  float * __restrict__ b_V_x=cpml.b.V_x;
  float * __restrict__ b_V_z=cpml.b.V_z;

  float * __restrict__ c_V_x=cpml.c.V_x;
  float * __restrict__ c_V_z=cpml.c.V_z;

  float * __restrict__ k_V_x=cpml.k.V_x;
  float * __restrict__ k_V_z=cpml.k.V_z;


  ord=ORD8;
  for(int iz=0;iz<nz;iz++){
	for(int ix=0;ix<nx;ix++){
	  int in_idx=ix+iz*nx;
	  int ord=ix;
	  ord=min(ord,nx-1-ix);
	  ord=min(ord,nz-1-iz);
	  ord=min(ord,ORD8);


	  float MUX_ix_iz= MUX[in_idx];
	  float MUZ_ix_iz= MUZ[in_idx];

	  DxV;
	  DzV;
	  if(ix<npml){
		//left pml
		idx=ix+iz*2*npml;
		psi_V_x[idx] = b_V_x[ix] * psi_V_x[idx] + c_V_x[ix] * dxV;
		Tx[in_idx] += MUX_ix_iz*( dxV * k_V_x[ix] + psi_V_x[idx] );
	  }else if(ix>=nx-npml){
		//right pml
		idx=npml+nx-1-ix+iz*2*npml;
		psi_V_x[idx] = b_V_x[ix] * psi_V_x[idx] + c_V_x[ix] * dxV;
		Tx[in_idx] += MUX_ix_iz*( dxV * k_V_x[ix] + psi_V_x[idx] );
	  }else{
		Tx[in_idx] += MUX_ix_iz*( dxV );
	  }

	  if(iz<npml && topabsorb){
		//top pml
		idx=(iz*nx)+ix;
		psi_V_z[idx] = b_V_z[iz] * psi_V_z[idx] + c_V_z[iz] * dzV;
		Tz[in_idx] += MUZ_ix_iz*( dzV * k_V_z[iz] + psi_V_z[idx] ) ;
	  }else if(iz>=nz-npml){
		// bottom
		idx=(npml+nz-1-iz)*nx+ix;
		psi_V_z[idx] = b_V_z[iz] * psi_V_z[idx] + c_V_z[iz] * dzV;
		Tz[in_idx] += MUZ_ix_iz*( dzV * k_V_z[iz] + psi_V_z[idx] ) ;
	  }else{
		Tz[in_idx] += MUZ_ix_iz*( dzV );
	  }

	}
  }

}

void step_velocity(FIELD & fld,MATERIAL & mat,CPML & cpml)
{
  int ord,idx;
  float Tx_x,Tz_z;

  int nx=fld.nx;
  int nz=fld.nz;

  int npml=cpml.npml;

  float * __restrict__ Tx=fld.Tx;
  float * __restrict__ Tz=fld.Tz;
  float * __restrict__ V=fld.V;

  float * __restrict__ BV=mat.BV;

  float * __restrict__ psi_Tx_x=cpml.psi.Tx_x;
  float * __restrict__ psi_Tz_z=cpml.psi.Tz_z;

  float * __restrict__ b_Tx_x=cpml.b.Tx_x;
  float * __restrict__ b_Tz_z=cpml.b.Tz_z;

  float * __restrict__ c_Tx_x=cpml.c.Tx_x;
  float * __restrict__ c_Tz_z=cpml.c.Tz_z;

  float * __restrict__ k_Tx_x=cpml.k.Tx_x;
  float * __restrict__ k_Tz_z=cpml.k.Tz_z;

  for(int i=0;i<nx;i++){
	/* stress imaging */
	Tz[i]=-Tz[i+nx];  
  }

  ord=ORD8;
  for(int iz=0;iz<nz;iz++){
	for(int ix=0;ix<nx;ix++){
	  int in_idx=ix+iz*nx;
	  float BV_ix_iz= mat.BV[in_idx];
	  ord=ix; // maybe top layer can also use 8th ,this reduce dispersion 
	  ord=min(ord,nx-1-ix);
	  ord=min(ord,nz-1-iz);
	  ord=min(ord,ORD8);
	  TxDIFF;
	  TzDIFF;
	  if(ix<npml){
		//left pml
		idx=ix+iz*2*npml;
		psi_Tx_x[idx] = b_Tx_x[ix] * psi_Tx_x[idx] + c_Tx_x[ix] * Tx_x;
		V[in_idx] += BV_ix_iz*( Tx_x * k_Tx_x[ix] + psi_Tx_x[idx] );
	  }else if(ix>=nx-npml){
		//right pml
		idx=npml+nx-1-ix+iz*2*npml;
		psi_Tx_x[idx] = b_Tx_x[ix] * psi_Tx_x[idx] + c_Tx_x[ix] * Tx_x;
		V[in_idx] += BV_ix_iz*( Tx_x * k_Tx_x[ix] + psi_Tx_x[idx] );
	  }else{
		V[in_idx] += BV_ix_iz*( Tx_x );
	  } 
	  if(iz<npml && topabsorb){
		//top pml
		idx=(iz*nx)+ix;
		psi_Tz_z[idx] = b_Tz_z[iz] * psi_Tz_z[idx] + c_Tz_z[iz] * Tz_z;
		V[in_idx] += BV_ix_iz*(Tz_z * k_Tz_z[iz] + psi_Tz_z[idx] );
	  }else if(iz>=nz-npml){
		// bottom
		idx=(npml+nz-1-iz)*nx+ix;
		psi_Tz_z[idx] = b_Tz_z[iz] * psi_Tz_z[idx] + c_Tz_z[iz] * Tz_z;
		V[in_idx] += BV_ix_iz*(Tz_z * k_Tz_z[iz] + psi_Tz_z[idx] );
	  }else{
		V[in_idx] += BV_ix_iz*(Tz_z);
	  }
	}
  }
}
