#include"box.h"
#include<iostream>
#include<cstdlib>
#include"../param/const.h"
#include"../param/param.h"
#include"../material/material.h"
#include"../field/field.h"

using namespace std;
void BOX::init_box(PARAM & param,MATERIAL &mat)
{
  int nx=param.nx;
  int xs=param.xs;
  int zs=param.zs;
  char * sourcefilename=param.sourcefile;
  dt= param.dt;
  h =param.h;

  srcden=param.srcden;
  srcsvel=param.srcsvel;
  srcpvel=param.srcpvel;

  mu=param.srcden*param.srcsvel*param.srcsvel;
  gam=param.srcden*param.srcpvel*param.srcpvel;
  lam=gam-2.0*mu;

  bw=param.bw;
  boxdevice=param.boxdevice;

  int s[8];
  int idxend=(bw+1)/2;        // 11
  int idxstart=idxend+1-bw; //-10
  s[0]=idxstart;
  s[1]=s[0]+3;
  s[2]=s[0]+4;
  s[3]=s[0]+7;
  s[7]=idxend;
  s[6]=s[7]-3;
  s[5]=s[7]-4;
  s[4]=s[7]-7;

  //for the box maxord,which index start form 0, source is zs_b
  int xs_b=-idxstart; //10, 
  int zs_b=-idxstart; //10

  start_x=xs+idxstart;
  start_z=zs+idxstart;

  /*
   * bw         : box leading diminsion
   * xs_b,zs_b: source position in box matrix
   * xs,zs    : source position in total field 
   * set idx accoring to s0---s7,xs,zs
   */
  int *tmp1=new int[bw*bw];
  int *tmp2=new int[bw*bw];
  int i,j,k;

  //sub idx s0-s1, s6-s7
  k=-1;
  for (j=s[0];j<=s[7];j++)
	for (i=s[0];i<=s[7];i++)
	{
	  if( ! (i>s[1] && i<s[6] && j>s[1] && j<s[6]) )
	  {
		k=k+1;
		tmp1[k]=(j+zs_b)*bw+i+xs_b;
		tmp2[k]=(j+zs)*nx+i+xs;
	  }
	}

  n_sub=k+1;
  r_sub_idx=new int[n_sub];
  a_sub_idx=new int[n_sub];
  for(i=0;i<n_sub;i++) r_sub_idx[i]=tmp1[i];
  for(i=0;i<n_sub;i++) a_sub_idx[i]=tmp2[i];

  //add idx, s2-s3 s4-s5
  k=-1;
  for (j=s[2];j<=s[5];j++)
	for (i=s[2];i<=s[5];i++)
	{
	  if( ! (i>s[3] && i<s[4] && j>s[3] && j<s[4]) )
	  {
		k=k+1;
		tmp1[k]=(j+zs_b)*bw+i+xs_b;
		tmp2[k]=(j+zs)*nx+i+xs;
	  }
	}
  n_add=k+1;
  r_add_idx=new int[n_add];
  a_add_idx=new int[n_add];

  for(i=0;i<n_add;i++) r_add_idx[i]=tmp1[i];
  for(i=0;i<n_add;i++) a_add_idx[i]=tmp2[i];

  delete[] tmp1;
  delete[] tmp2;


  u1=new double[(bw+8)*(bw+8)];
  u2=new double[(bw+8)*(bw+8)];
  w1=new double[(bw+8)*(bw+8)];
  w2=new double[(bw+8)*(bw+8)];
  
  sourcefile.open(sourcefilename,std::ios::binary);
  sourcefile.read((char*)u1,sizeof(double)*(bw+8)*(bw+8));
  sourcefile.read((char*)w1,sizeof(double)*(bw+8)*(bw+8));
  uold=u1;
  wold=w1;
  u1_read_new=false;


  p_sfield=new FIELD();
  p_sfield->init_for_box(param);
  p_bmat=new MATERIAL();
  p_bmat->init_for_box_from_full(bw,start_x,start_z,mat);

}
void BOX::setboxsource()
{
  int i,j;
  if(u1_read_new){
	sourcefile.read((char*)u1,sizeof(double)*(bw+8)*(bw+8));
	sourcefile.read((char*)w1,sizeof(double)*(bw+8)*(bw+8));
	unew=u1;
	wnew=w1;
	uold=u2;
	wold=w2;
	u1_read_new=false;
  }else{
	sourcefile.read((char*)u2,sizeof(double)*(bw+8)*(bw+8));
	sourcefile.read((char*)w2,sizeof(double)*(bw+8)*(bw+8));
	unew=u2;
	wnew=w2;
	uold=u1;
	wold=w1;
	u1_read_new=true;
  }

#define Txx(i,j) p_sfield->Txx[i+(j)*bw]
#define Txz(i,j) p_sfield->Txz[i+(j)*bw]
#define Tzz(i,j) p_sfield->Tzz[i+(j)*bw]
#define Udt(i,j) p_sfield->U  [i+(j)*bw]
#define Wdt(i,j) p_sfield->W  [i+(j)*bw]
#define  un(i,j) unew[i+(j)*(bw+8)]
#define  uo(i,j) uold[i+(j)*(bw+8)]
#define  wn(i,j) wnew[i+(j)*(bw+8)]
#define  wo(i,j) wold[i+(j)*(bw+8)]


  for(int ty=4;ty<bw+4;ty++){
	for(int tx=4;tx<bw+4;tx++){

	  int ord=min(4,max( abs(ty-3-(bw+1)/2)-1 ,abs(tx-3-(bw+1)/2)-1));
	  ord=max(ord,0);

	  Udt(tx-4,ty-4) =(un(tx,ty)-uo(tx,ty))/dt;
	  Wdt(tx-4,ty-4) =(wn(tx,ty)-wo(tx,ty))/dt;

	  double ux=g_coef[ord][0]*( un(tx  ,ty) - un(tx-1,ty) )
		-g_coef[ord][1]*( un(tx+1,ty) - un(tx-2,ty) )
		+g_coef[ord][2]*( un(tx+2,ty) - un(tx-3,ty) )
		-g_coef[ord][3]*( un(tx+3,ty) - un(tx-4,ty) );
	  double uz=g_coef[ord][0]*( un(tx  ,ty) - un(tx,ty-1) )
		-g_coef[ord][1]*( un(tx,ty+1) - un(tx,ty-2) )
		+g_coef[ord][2]*( un(tx,ty+2) - un(tx,ty-3) )
		-g_coef[ord][3]*( un(tx,ty+3) - un(tx,ty-4) );

	  double wx=g_coef[ord][0]*( wn(tx+1,ty) - wn(tx,ty  ) ) 
		-g_coef[ord][1]*( wn(tx+2,ty) - wn(tx-1,ty) ) 
		+g_coef[ord][2]*( wn(tx+3,ty) - wn(tx-2,ty) ) 
		-g_coef[ord][3]*( wn(tx+4,ty) - wn(tx-3,ty) );
	  double wz=g_coef[ord][0]*( wn(tx,ty+1) - wn(tx,ty) )
		-g_coef[ord][1]*( wn(tx,ty+2) - wn(tx,ty-1) )
		+g_coef[ord][2]*( wn(tx,ty+3) - wn(tx,ty-2) )
		-g_coef[ord][3]*( wn(tx,ty+4) - wn(tx,ty-3) );

	  // Source region is homogenous
	  /*
	  Txx(tx-4,ty-4) =( gam*ux + lam*wz )/h;
	  Tzz(tx-4,ty-4) =( lam*ux + gam*wz )/h;
	  Txz(tx-4,ty-4) =( mu*(wx +uz)     )/h;
	  */
	  
	  // Test using inhomogenous source
	  double LAM=p_bmat->LAM[ (tx-4) + (ty-4)*bw ];
	  double MU =p_bmat->MU [ (tx-4) + (ty-4)*bw ];
	  double MUA=p_bmat->MUA[ (tx-4) + (ty-4)*bw ];
	  double GAM=LAM+2.0*MU;
	  Txx(tx-4,ty-4) =( GAM*ux + LAM*wz )/dt;
	  Tzz(tx-4,ty-4) =( LAM*ux + GAM*wz )/dt;
	  Txz(tx-4,ty-4) =( MUA*(wx +uz)    )/dt;

	}                                      
  }
}
void BOX::subsource(double *pa,double* pr,double* ps)
{
  for(int i=0;i<n_sub;i++){
	pr[r_sub_idx[i]]=pa[a_sub_idx[i]]-ps[r_sub_idx[i]];
  }
}

void BOX::addsource(double *pa,double* pr,double* ps)
{
  for(int i=0;i<n_add;i++){
	pa[a_add_idx[i]]=pr[r_add_idx[i]]+ps[r_add_idx[i]];
  }
}
