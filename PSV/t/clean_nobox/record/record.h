/* Author Dunzhu Li  dli@caltech.edu 
*/
#ifndef _RECORD_H_
#define _RECORD_H_
#include<fstream>
#include"../param/const.h"
class PARAM;
class FIELD;
class MATERIAL;

class RECORD{
  public:

	float * d_recbuf[MAXCUDAGPU];   // in GPU
	float * g_d_recbuf[MAXCUDAGPU]; // in global
	float * g_rec;                  // continous in global

	int nrec;
	int *rx;
	int *rz;

	int *d_rx[MAXCUDAGPU];
	int *d_rz[MAXCUDAGPU];

	int *ind_rec[MAXCUDAGPU];
	int lrec[MAXCUDAGPU];

	char recordfile[256];

	std::ofstream fid;


	void init_line(int nrec, int ix0, int iz0, int idx, int idz,const char* output); //line record
	void init_list(const char * sta_file,const char* output);                        //read station from file
	void cu_init_record(int deviceid,PARAM &param);
	void cu_record(int deviceid,int nx, float *U); 
	void cu_flush(int deviceid);

	void close();
	void toisis(PARAM & param,float shiftx,float shiftz);

};

void snapshot(FIELD & fld,int it);
void snapshot(MATERIAL & fld);
int cpy_backhost(int deviceid,PARAM & param,FIELD & fld,FIELD & d_fld);
#endif
