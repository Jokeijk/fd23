#include"../param/const.h"
#include"../param/param.h"
#include"waveform.h"
#include<cmath>
#include<complex>
#include<fstream>
using namespace std;

void compute_QW_complex(double r,double z,double pvel,double svel,double den,double dt,int nt,double *Q,double *W,double *T1,double *T2,char comp,double A1,double A2,double A3)
{

  double R=sqrt(r*r+z*z);

  double G=sqrt(2.0)/(4.0*PI*PI*den);

  int n1_a=floor(R/pvel/dt);
  int n1_b=floor(R/svel/dt);

  double a=(n1_a+1)-R/pvel/dt+0.5;
  double b=(n1_b+1)-R/svel/dt+0.5;

  double first_a=sqrt(a*(2.0*R/pvel/dt+a));
  double first_b=sqrt(b*(2.0*R/svel/dt+b));

  for(int i=0;i<=n1_a;i++) T1[i]=0;
  for(int i=0;i<=n1_b;i++) T2[i]=0;


  T1[n1_a+1]=first_a;
  T2[n1_b+1]=first_b;

  for(int i=n1_a+2;i<nt;i++) {
	double t=i*dt;
	double ta=(1.0/(t*pvel)*R)*(1.0/(t*pvel)*R);
	T1[i]=1./(sqrt(1.0-ta));
  }

  for(int i=n1_b+2;i<nt;i++) {
	double t=i*dt;
	double tb=(1.0/(t*svel)*R)*(1.0/(t*svel)*R);
	T2[i]=1./(sqrt(1.0-tb));
  }

  double R2=r*r+z*z;
  double iR2=1.0/R2;
  double ipvel2=1./(pvel*pvel);
  double isvel2=1./(svel*svel);
  double R2ipvel2=R2*ipvel2;
  double R2isvel2=R2*isvel2;
  double absziR2=abs(z)*iR2;
  double riR2=r*iR2;
  double eps=z>0?1:-1;

  if(comp=='W'){
	for(int i=0;i<n1_a;i++) W[i]=0.0;

	for(int i=n1_a;i<nt;i++){
	  double t=i*dt;

	  double ta;
	  double tb;
	  if(t*t>R2ipvel2){ 
		ta=sqrt(t*t-R2ipvel2);
	  }else{
		ta=0.0;
	  };
	  if(t*t>R2isvel2){ 
		tb=sqrt(t*t-R2isvel2);
	  }else{
		tb=0.0;
	  };

	  complex<double> pa( r*iR2*t , ta*absziR2 );
	  complex<double> na( absziR2*t , -ta*riR2 );
	  complex<double> pb( r*iR2*t , tb*absziR2 );
	  complex<double> nb( absziR2*t , -tb*riR2 );
	  complex<double> dpdt_na(0.0,T1[i]/t);
	  complex<double> dpdt_nb(0.0,T2[i]/t);

	  complex<double> c1=-pa*pa;
	  complex<double> c2=2.0*eps*pa*na;
	  complex<double> c3=pa*pa-2.0*na*na;

	  complex<double> sv1=-eps*pb*nb;
	  complex<double> sv2=nb*nb-pb*pb;
	  complex<double> sv3=3.0*eps*pb*nb;

	  complex<double> wa=-eps*na;
	  complex<double> wb=pb;

	  complex<double> result=
	   (wa*c1*dpdt_na+wb*sv1*dpdt_nb)*A1
	  +(wa*c2*dpdt_na+wb*sv2*dpdt_nb)*A2
	  +(wa*c3*dpdt_na+wb*sv3*dpdt_nb)*A3;


	  W[i]=G*imag(result);
	}
  }

  if(comp=='Q'){
	for(int i=0;i<n1_a;i++) Q[i]=0.0;

	for(int i=n1_a;i<nt;i++){
	  double t=i*dt;

	  double ta;
	  double tb;
	  if(t*t>R2ipvel2){ 
		ta=sqrt(t*t-R2ipvel2);
	  }else{
		ta=0.0;
	  };
	  if(t*t>R2isvel2){ 
		tb=sqrt(t*t-R2isvel2);
	  }else{
		tb=0.0;
	  };


	  complex<double> pa( r*iR2*t , ta*absziR2 );
	  complex<double> na( absziR2*t , -ta*riR2 );
	  complex<double> pb( r*iR2*t , tb*absziR2 );
	  complex<double> nb( absziR2*t , -tb*riR2 );
	  complex<double> dpdt_na(0.0,T1[i]/t);
	  complex<double> dpdt_nb(0.0,T2[i]/t);

	  complex<double> c1=-pa*pa;
	  complex<double> c2=2.0*eps*pa*na;
	  complex<double> c3=pa*pa-2.0*na*na;

	  complex<double> sv1=-eps*pb*nb;
	  complex<double> sv2=nb*nb-pb*pb;
	  complex<double> sv3=3.0*eps*pb*nb;

	  complex<double> qa=-pa;
	  complex<double> qb=-eps*nb;

	  complex<double> result=
	   (qa*c1*dpdt_na+qb*sv1*dpdt_nb)*A1
	  +(qa*c2*dpdt_na+qb*sv2*dpdt_nb)*A2
	  +(qa*c3*dpdt_na+qb*sv3*dpdt_nb)*A3;


	  Q[i]=G*imag(result);
	}
  }
}

void compute_pQpW_complex(double r,double z,double pvel,double svel,double den,double dt,int nt,double *Q,double *W,double *T1,double *T2,char comp,double A1,double A2,double A3)
{

  double R=sqrt(r*r+z*z);

  double G=sqrt(2.0)/(4.0*PI*PI*den);

  int n1_a=floor(R/pvel/dt);
  int n1_b=floor(R/svel/dt);

  double a=(n1_a+1)-R/pvel/dt+0.5;
  double b=(n1_b+1)-R/svel/dt+0.5;

  double first_a=sqrt(a*(2.0*R/pvel/dt+a));
  double first_b=sqrt(b*(2.0*R/svel/dt+b));

  for(int i=0;i<=n1_a;i++) T1[i]=0;
  for(int i=0;i<=n1_b;i++) T2[i]=0;


  T1[n1_a+1]=first_a;
  T2[n1_b+1]=first_b;

  for(int i=n1_a+2;i<nt;i++) {
	double t=i*dt;
	double ta=(1.0/(t*pvel)*R)*(1.0/(t*pvel)*R);
	T1[i]=1./(sqrt(1.0-ta));
  }

  for(int i=n1_b+2;i<nt;i++) {
	double t=i*dt;
	double tb=(1.0/(t*svel)*R)*(1.0/(t*svel)*R);
	T2[i]=1./(sqrt(1.0-tb));
  }

  double R2=r*r+z*z;
  double iR2=1.0/R2;
  double ipvel2=1./(pvel*pvel);
  double isvel2=1./(svel*svel);
  double R2ipvel2=R2*ipvel2;
  double R2isvel2=R2*isvel2;
  double absziR2=abs(z)*iR2;
  double riR2=r*iR2;
  double eps=z>0?1:-1;

  if(comp=='W'){
	for(int i=0;i<n1_a;i++) W[i]=0.0;

	for(int i=n1_a;i<nt;i++){
	  double t=i*dt;

	  double ta;
	  double tb;
	  if(t*t>R2ipvel2){ 
		ta=sqrt(t*t-R2ipvel2);
	  }else{
		ta=0.0;
	  };
	  if(t*t>R2isvel2){ 
		tb=sqrt(t*t-R2isvel2);
	  }else{
		tb=0.0;
	  };

	  complex<double> pa( r*iR2*t , ta*absziR2 );
	  complex<double> na( absziR2*t , -ta*riR2 );
	  complex<double> pb( r*iR2*t , tb*absziR2 );
	  complex<double> nb( absziR2*t , -tb*riR2 );
	  complex<double> dpdt_na(0.0,T1[i]/t);
	  complex<double> dpdt_nb(0.0,T2[i]/t);

	  complex<double> c1=-pa*pa;
	  complex<double> c2=2.0*eps*pa*na;
	  complex<double> c3=pa*pa-2.0*na*na;

	  complex<double> sv1=-eps*pb*nb;
	  complex<double> sv2=nb*nb-pb*pb;
	  complex<double> sv3=3.0*eps*pb*nb;

	  complex<double> wa=-eps*na;
	  complex<double> wb=pb;

	  complex<double> result=
	   (wa*c1*dpdt_na*pa+wb*sv1*dpdt_nb*pb)*A1
	  +(wa*c2*dpdt_na*pa+wb*sv2*dpdt_nb*pb)*A2
	  +(wa*c3*dpdt_na*pa+wb*sv3*dpdt_nb*pb)*A3;


	  W[i]=G*imag(result);
	}
  }

  if(comp=='Q'){
	for(int i=0;i<n1_a;i++) Q[i]=0.0;

	for(int i=n1_a;i<nt;i++){
	  double t=i*dt;

	  double ta;
	  double tb;
	  if(t*t>R2ipvel2){ 
		ta=sqrt(t*t-R2ipvel2);
	  }else{
		ta=0.0;
	  };
	  if(t*t>R2isvel2){ 
		tb=sqrt(t*t-R2isvel2);
	  }else{
		tb=0.0;
	  };


	  complex<double> pa( r*iR2*t , ta*absziR2 );
	  complex<double> na( absziR2*t , -ta*riR2 );
	  complex<double> pb( r*iR2*t , tb*absziR2 );
	  complex<double> nb( absziR2*t , -tb*riR2 );
	  complex<double> dpdt_na(0.0,T1[i]/t);
	  complex<double> dpdt_nb(0.0,T2[i]/t);

	  complex<double> c1=-pa*pa;
	  complex<double> c2=2.0*eps*pa*na;
	  complex<double> c3=pa*pa-2.0*na*na;

	  complex<double> sv1=-eps*pb*nb;
	  complex<double> sv2=nb*nb-pb*pb;
	  complex<double> sv3=3.0*eps*pb*nb;

	  complex<double> qa=-pa;
	  complex<double> qb=-eps*nb;

	  complex<double> result=
	   (qa*c1*dpdt_na*pa+qb*sv1*dpdt_nb*pb)*A1
	  +(qa*c2*dpdt_na*pa+qb*sv2*dpdt_nb*pb)*A2
	  +(qa*c3*dpdt_na*pa+qb*sv3*dpdt_nb*pb)*A3;


	  Q[i]=G*imag(result);
	}
  }
}

void cu_con(int N1,double *f1,int N2,double *f2);
int main(int ac,char **av) 
{
  double src_shift_x= 0.25;
  double src_shift_z=-0.25;
  PARAM param;
  param.read_param(ac,av,"source");

  int bw=param.bw;
  double pvel=param.srcpvel;
  double svel=param.srcsvel;
  double den=param.srcden;
  double dt=param.dt;
  double h=param.h;
  int nt=param.nt;

  double delta=param.dip*PI/180.0;
  double lamda=param.rake*PI/180.0;
  double theta=(param.azimuth-param.strike)*PI/180.0;

  double A1=sin(2.0*theta)*cos(lamda)*sin(delta)+1.0/2.0*cos(2.0*theta)*sin(lamda)*sin(2.0*delta);
  double A2=cos(theta)*cos(lamda)*cos(delta)-sin(theta)*sin(lamda)*cos(2.0*delta);
  double A3=1.0/2.0*sin(lamda)*sin(2.0*delta);
  
  double *src;
  int lsrc;
  

  if(param.sourcetime[0]=='T' || param.sourcetime[0]=='t'){
	double trap1=param.trap1;
	double trap2=param.trap2;
	double trap3=param.trap3;
	trapozoid(dt,trap1,trap2,trap3,lsrc,src);
  }else if(param.sourcetime[0]=='G' || param.sourcetime[0]=='g'){
	int alpha=param.alpha;
	gaussian(dt,alpha ,lsrc,src);
  }

  int space=(bw+8)*(bw+8);
  int idxend=(bw+1)/2+4;  // index,or position of txx,tzz
  int idxstart=idxend+1-(bw+8);
  double *Q=new double[nt*space];
  double *W=new double[nt*space];
  double *T1=new double[nt];
  double *T2=new double[nt];
  int k=-1;
  for(int isz=idxstart;isz<=idxend;isz++){
	  for(int isx=idxstart;isx<=idxend;isx++){
		// source is src_shift_x,src_shift_z,
		// u is at (isx+0.5,isz)
		k=k+1;
		double x=(isx+0.5)-(src_shift_x);
	  	double z=isz-src_shift_z;
		if(param.stype==13){
		  compute_pQpW_complex( x*h, z*h, pvel, svel, den, dt, nt, Q+k*nt, W+k*nt, T1, T2,'Q', A1, A2, A3);
		}else{
		  compute_QW_complex( x*h, z*h, pvel, svel, den, dt, nt, Q+k*nt, W+k*nt, T1, T2,'Q', A1, A2, A3);
		}
        cu_con(lsrc,src,nt,Q+k*nt);

		// w is at (isx,isz-0.5)
		x=isx-src_shift_x;
	  	z=(isz-0.5)-src_shift_z;

		if(param.stype==13){
		  compute_pQpW_complex( x*h, z*h, pvel, svel, den, dt, nt, Q+k*nt, W+k*nt, T1, T2,'W', A1, A2, A3);
		}else{
		  compute_QW_complex( x*h, z*h, pvel, svel, den, dt, nt, Q+k*nt, W+k*nt, T1, T2,'W', A1, A2, A3);
		}
        cu_con(lsrc,src,nt,W+k*nt);
	  }
  }

  std::ofstream fid;
  fid.open(param.sourcefile,std::ios::binary);

  double *buf=new double[2*space];
  for(int it=0;it<nt;it++){
	for(int i=0;i<space;i++){
	  buf[i]=double(Q[it+i*nt]);
	  buf[i+space]=double(W[it+i*nt]);
	}
	fid.write((char*)buf,sizeof(double)*2*space);
  }
  fid.close();
  printf("finish!");
}
