#include "../param/param.h"
#include "waveform.h"
#include<cmath>
#include<fstream>

void compute_QW(double r,double z,double pvel,double svel,double den,double dt,int nt,double *Q,double *W,double *T1,double *T2,char comp,double A1,double A2,double A3)
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


  double R6=R*R*R*R*R*R;
  double R4=R*R*R*R;

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

  if(comp=='W'){
	for(int i=0;i<n1_a;i++) W[i]=0.0;

	for(int i=n1_a;i<=std::min(nt-1,2*n1_b);i++){
	  double t=i*dt;

	  double Phia=G*t*t/R6*T1[i];
	  double Phib=G*t*t/R6*T2[i];

	  double PhiaTa=G/R4/(pvel*pvel)*T1[i];
	  double PhibTb=G/R4/(svel*svel)*T2[i];

	  double W1= -z*(  Phia*(z*z-3.0f*r*r)+PhiaTa*(2.0f*r*r-z*z) 
		  +Phib*(3.0f*r*r-z*z)+PhibTb*(z*z-2.0f*r*r) ) ;

	  double W2=  r*(  Phia*(2.0f*r*r-6.0f*z*z)+PhiaTa*(4.0f*z*z-2.0f*r*r)
		  +Phib*(6.0f*z*z-2*r*r   )+PhibTb*(r*r-5.0f*z*z)  );

	  double W3= -z*(  Phia*(9.0f*r*r-3.0f*z*z)+PhiaTa*(z*z-8.0f*r*r)
		  +Phib*(3.0f*z*z-9.0f*r*r)+PhibTb*(6.0f*r*r-3.0f*z*z) );

	  W[i]=A1*W1+A2*W2+A3*W3;
	}

	for(int i=std::min(nt-1,2*n1_b)+1;i<nt;i++){
	  double t=i*dt;
	  double ta=(1.0/(t*pvel)*R)*(1.0/(t*pvel)*R);
	  double tb=(1.0/(t*svel)*R)*(1.0/(t*svel)*R);
	  double sa=sqrt(1.0-ta);
	  double sb=sqrt(1.0-tb);

	  double Phia_Phib=G/R4*( 1.0/(pvel*pvel)-1.0/(svel*svel) )/
		( (sa+sb)*sa*sb );

	  double PhiaTa=G/R4/(pvel*pvel)*T1[i];
	  double PhibTb=G/R4/(svel*svel)*T2[i];

	  double W1= -z*(  Phia_Phib*(z*z-3.0f*r*r)
		  +PhiaTa*(2.0f*r*r-z*z) 
		  +PhibTb*(z*z-2.0f*r*r) ) ;

	  double W2=  r*(  Phia_Phib*(2.0f*r*r-6.0f*z*z)
		  +PhiaTa*(4.0f*z*z-2.0f*r*r)
		  +PhibTb*(r*r-5.0f*z*z)  );

	  double W3= -z*(  Phia_Phib*(9.0f*r*r-3.0f*z*z)
		  +PhiaTa*(z*z-8.0f*r*r)
		  +PhibTb*(6.0f*r*r-3.0f*z*z) );

	  W[i]=A1*W1+A2*W2+A3*W3;
	}

  }
  if(comp=='Q'){
	for(int i=0;i<n1_a;i++) Q[i]=0.0;

	for(int i=n1_a;i<=std::min(nt-1,2*n1_b);i++){
	  double t=i*dt;

	  double Phia=G*t*t/R6*T1[i];
	  double Phib=G*t*t/R6*T2[i];

	  double PhiaTa=G/R4/(pvel*pvel)*T1[i];
	  double PhibTb=G/R4/(svel*svel)*T2[i];

	  double Q1= r* (  Phia*(r*r-3*z*z) +3*PhiaTa*z*z
		  +Phib*(3*z*z-r*r)+PhibTb*(r*r-2*z*z)       );

	  double Q2= -z*(  Phia*(6*r*r-2*z*z)+PhiaTa*(2*z*z-4*r*r)
		  +Phib*(2*z*z-6*r*r)+PhibTb*(5*r*r-z*z)           );

	  double Q3=  r*(  Phia*(9*z*z-3*r*r)+PhiaTa*(2*r*r-7*z*z)
		  +Phib*(3*r*r-9*z*z)+PhibTb*(6*z*z-3*r*r)         );

	  Q[i]=A1*Q1+A2*Q2+A3*Q3;
	}

	for(int i=std::min(nt-1,2*n1_b)+1;i<nt;i++){
	  double t=i*dt;
	  double ta=(1.0/(t*pvel)*R)*(1.0/(t*pvel)*R);
	  double tb=(1.0/(t*svel)*R)*(1.0/(t*svel)*R);
	  double sa=sqrt(1.0-ta);
	  double sb=sqrt(1.0-tb);

	  double Phia_Phib=G/R4*( 1.0/(pvel*pvel)-1.0/(svel*svel) )/
		( (sa+sb)*sa*sb );

	  double PhiaTa=G/R4/(pvel*pvel)*T1[i];
	  double PhibTb=G/R4/(svel*svel)*T2[i];

	  double Q1= r* (  Phia_Phib*(r*r-3*z*z) 
		  +3*PhiaTa*z*z
		  +PhibTb*(r*r-2*z*z)       );

	  double Q2= -z*(  Phia_Phib*(6*r*r-2*z*z)
		  +PhiaTa*(2*z*z-4*r*r)
		  +PhibTb*(5*r*r-z*z)           );

	  double Q3=  r*(  Phia_Phib*(9*z*z-3*r*r)
		  +PhiaTa*(2*r*r-7*z*z)
		  +PhibTb*(6*z*z-3*r*r)         );

	  Q[i]=A1*Q1+A2*Q2+A3*Q3;

	}

  }

}

void compute_QW_mathematica(double r,double z,double pvel,double svel,double den,double dt,int nt,double *Q,double *W,double *T1,double *T2,char comp,double A1,double A2,double A3)
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

  double r2=r*r;
  double r4=r2*r2;
  double z2=z*z;
  double z4=z2*z2;
  double R2=r2+z2;
  double iR6=1./(R*R*R*R*R*R);
  double ipvel2=1./(pvel*pvel);
  double isvel2=1./(svel*svel);

  if(comp=='W'){
	for(int i=0;i<n1_a;i++) W[i]=0.0;

	for(int i=n1_a;i<nt;i++){
	  double t=i*dt;
	  double t2=t*t;
	  double Ta=T1[i];
	  double Tb=T2[i];

	  double W1=iR6*t2*z*(Ta*(3*r2 - z2) + Tb*(-3*r2 + z2)) + iR6*z*(isvel2*Tb*(2*r4 + r2*z2 - z4) + ipvel2*Ta*(-2*r4 - r2*z2 + z4));
	  double W2=iR6*r*t2*(2*Ta*(r2 - 3*z2) - 2*Tb*(r2 - 3*z2)) + iR6*r*(isvel2*Tb*(r4 - 4*r2*z2 - 5*z4) + ipvel2*Ta*(-2*r4 + 2*r2*z2 + 4*z4));
	  double W3=iR6*t2*z*(-3*Ta*(3*r2 - z2) + 3*Tb*(3*r2 - z2)) + iR6*z*(ipvel2*Ta*(8*r4 + 7*r2*z2 - z4) + isvel2*Tb*(-6*r4 - 3*r2*z2 + 3*z4));


	  W[i]=G*(A1*W1+A2*W2+A3*W3);
	}
  }
  if(comp=='Q'){
	for(int i=0;i<n1_a;i++) Q[i]=0.0;

	for(int i=n1_a;i<nt;i++){
	  double t=i*dt;
	  double t2=t*t;
	  double Ta=T1[i];
	  double Tb=T2[i];

	  double Q1=iR6*r*t2*(Ta*(r2 - 3*z2) - Tb*(r2 - 3*z2)) + iR6*r*(3*ipvel2*R2*Ta*z2 + isvel2*Tb*(r4 - r2*z2 - 2*z4));
	  double Q2=iR6*t2*z*(-2*Ta*(3*r2 - z2) + 2*Tb*(3*r2 - z2)) + iR6*z*(ipvel2*Ta*(4*r4 + 2*r2*z2 - 2*z4) + isvel2*Tb*(-5*r4 - 4*r2*z2 + z4));
	  double Q3=iR6*r*t2*(-3*Ta*(r2 - 3*z2) + 3*Tb*(r2 - 3*z2)) + iR6*r*(ipvel2*Ta*(2*r4 - 5*r2*z2 - 7*z4) + isvel2*Tb*(-3*r4 + 3*r2*z2 + 6*z4));


	  Q[i]=G*(A1*Q1+A2*Q2+A3*Q3);
	}
  }

}

void compute_pQpW_mathematica(double r,double z,double pvel,double svel,double den,double dt,int nt,double *Q,double *W,double *T1,double *T2,char comp,double A1,double A2,double A3)
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

  double r2=r*r;
  double r4=r2*r2;
  double r6=r2*r2*r2;

  double z2=z*z;
  double z4=z2*z2;
  double z6=z2*z2*z2;
  double z8=z4*z4;

  double R2=r2+z2;
  double R4=R2*R2;

  double iR8=1./(R2*R2*R2*R2);

  double ipvel2=1./(pvel*pvel);
  double isvel2=1./(svel*svel);
  double ipvel4=ipvel2*ipvel2;
  double isvel4=isvel2*isvel2;

  if(comp=='W'){
	for(int i=0;i<n1_a;i++) W[i]=0.0;

	for(int i=n1_a;i<nt;i++){
	  double t=i*dt;
	  double t3=t*t*t;
	  double Ta=T1[i];
	  double Tb=T2[i];

	  double W1=iR8*r*t3*z*(4*Ta*(r2 - z2) - 4*Tb*(r2 - z2)) + (iR8*r*z*(-(ipvel4*R4*Ta*z2) + isvel4*R4*Tb*z2))/t + iR8*r*t*z*(-(ipvel2*Ta*(3*r4 - 2*r2*z2 - 5*z4)) + isvel2*Tb*(3*r4 - 2*r2*z2 - 5*z4));
	  double W2=(iR8*(-2*ipvel4*r2*R4*Ta*z2 + isvel4*R4*Tb*(r2 - z2)*z2))/t + iR8*t3*(2*Ta*(r4 - 6*r2*z2 + z4) - 2*Tb*(r4 - 6*r2*z2 + z4)) + iR8*t*(-2*ipvel2*Ta*(r6 - 5*r4*z2 - 5*r2*z4 + z6) + isvel2*Tb*(r6 - 11*r4*z2 - 9*r2*z4 + 3*z6));
	  double W3=iR8*r*t3*z*(-12*Ta*(r2 - z2) + 12*Tb*(r2 - z2)) + iR8*r*t*z*(ipvel2*Ta*(13*r4 + 2*r2*z2 - 11*z4) - 3*isvel2*Tb*(3*r4 - 2*r2*z2 - 5*z4)) + (iR8*r*z*(-3*isvel4*R4*Tb*z2 + ipvel4*Ta*(-2*r6 - 3*r4*z2 + z6)))/t;


	  W[i]=G*(A1*W1+A2*W2+A3*W3);
	}
  }
  if(comp=='Q'){
	for(int i=0;i<n1_a;i++) Q[i]=0.0;

	for(int i=n1_a;i<nt;i++){
	  double t=i*dt;
	  double t3=t*t*t;
	  double Ta=T1[i];
	  double Tb=T2[i];

	  double Q1=(iR8*(isvel4*r2*R4*Tb*z2 + ipvel4*R4*Ta*z4))/t + iR8*t3*(Ta*(r4 - 6*r2*z2 + z4) - Tb*(r4 - 6*r2*z2 + z4)) + iR8*t*(2*ipvel2*Ta*(3*r4*z2 + 2*r2*z4 - z6) + isvel2*Tb*(r6 - 5*r4*z2 - 5*r2*z4 + z6));
	  double Q2=iR8*r*t3*z*(8*Tb*(r - z)*(r + z) - 8*Ta*(r2 - z2)) + (iR8*r*z*(isvel4*R4*Tb*(r - z)*(r + z) + 2*ipvel4*R4*Ta*z2))/t + iR8*r*t*z*(2*ipvel2*Ta*(3*r4 - 2*r2*z2 - 5*z4) + 8*isvel2*Tb*(-r4 + z4));
	  double Q3=iR8*t3*(-3*Ta*(r4 - 6*r2*z2 + z4) + 3*Tb*(r4 - 6*r2*z2 + z4)) + iR8*t*(-3*isvel2*Tb*(r6 - 5*r4*z2 - 5*r2*z4 + z6) + 2*ipvel2*Ta*(r6 - 8*r4*z2 - 7*r2*z4 + 2*z6)) + (iR8*(-3*isvel4*r2*R4*Tb*z2 + ipvel4*Ta*(2*r6*z2 + 3*r4*z4 - z8)))/t;


	  Q[i]=G*(A1*Q1+A2*Q2+A3*Q3);
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
        compute_pQpW_mathematica( x*h, z*h, pvel, svel, den, dt, nt, Q+k*nt, W+k*nt, T1, T2,'Q', A1, A2, A3);
        cu_con(lsrc,src,nt,Q+k*nt);

		/*
		printf("lsrc and nt is:%d,%d",lsrc,nt);
		std::ofstream tempfid;
		tempfid.open("temp",std::ios::binary);
		tempfid.write((char*)Q,sizeof(double)*nt);
		tempfid.write((char*)src,sizeof(double)*lsrc);
		tempfid.write((char*)Q,sizeof(double)*nt);
		tempfid.close();
		return 0;
		*/

		// w is at (isx,isz-0.5)
		x=isx-src_shift_x;
	  	z=(isz-0.5)-src_shift_z;
        compute_pQpW_mathematica( x*h, z*h, pvel, svel, den, dt, nt, Q+k*nt, W+k*nt, T1, T2,'W', A1, A2, A3);
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
