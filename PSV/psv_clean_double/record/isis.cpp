#include<pthread.h>
#include<fstream>
#include"isis.h"
#include"../param/param.h"
#include"record.h"
#include<cstring>

void RECORD::toisis(PARAM &param,double shiftx,double shiftz)
{

  struct traceinfo isis;
  char filename[256];
  double arrivetime=0;

  std::ifstream fid_in;
  std::ofstream fid_out;
  std::ifstream timefid;


  int ntrecord=param.nt/param.itrecord;

  double *buf=new double[nrec*ntrecord];
  double *tmp=new double[ntrecord];


  fid_in.open(recordfile,std::ios::binary);
  fid_in.read((char*)buf,sizeof(double)*nrec*ntrecord);
  fid_in.close();

  sprintf(filename,"%s.db.isis",recordfile);
  fid_out.open(filename,std::ios::binary);


  for(int ir=0; ir<nrec; ir++){

	int ixr= rx[ir];
	int izr= rz[ir];

	isis.xr= (ixr+shiftx) * param.h;
	isis.yr= 0.0;
	isis.zr= (izr+shiftz) * param.h;
	isis.xs= (double)(param.xs+param.src_shift_x) * param.h;
	isis.ys= 0.0;
	isis.zs= (double)(param.zs+param.src_shift_z) * param.h;
	isis.nt= ntrecord;
	isis.samplerate= param.dt * (double)(param.itrecord);
	isis.t0= 0.0;
	isis.status= 0;
	isis.gain=param.dt/param.h; //useful for p=k/omega
	isis.dtstatic= 0.0;
	isis.ampstatic= 1.0;
	isis.srcflagnum= param.xs;
	isis.recflagnum= ixr;
	isis.user[1]=param.bw;
	isis.user[2]=123456789;
	fid_out.write((char*)&isis,TRSIZE);
	for(int k=0;k<ntrecord;k++) tmp[k]=buf[k*nrec+ir];
	fid_out.write((char*)tmp,sizeof(double)*ntrecord);
  }
  fid_out.close();
}

