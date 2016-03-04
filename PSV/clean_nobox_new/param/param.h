#ifndef _PARAM_H_
#define _PARAM_H_

#include "const.h"

class PARAM{
  public:
	//flags
	int restart;
	int nt0;
	int savestate;

	int usepunch;
	int usetable;

	//Plan
	int ngpu;
	int gpuid[MAXCUDAGPU];
	int a_nz1[MAXCUDAGPU];
	int a_nz2[MAXCUDAGPU];


	int snapbox[6];

	//dimension
	int	nx;
	int	nz;
	int	nt;
	float	h;
	float	dt;

	//earthquake
	float strike;
	float dip;
	float rake;
	float azimuth;

	char  srctype;   // 1 or p   
	char  srctime[STRLEN];
	float trap1;
	float trap2;
	float trap3;
	float alpha;

	//source
	int xs;
	int zs;


	//for pml
	int   npml;
	float pml_r;
	float pml_v;
	float pml_fc;
	float pml_dt;

	//for record
	int	itrecord;
	int	nrec;
	int	ixrec0;
	int	izrec0_u;
	int	izrec0_w;
	int	idxrec;
	int	idzrec;
	char output[STRLEN];
	char U_sta_file[STRLEN];
	char W_sta_file[STRLEN];
	bool use_sta_file;

	//info
	int	itprint;
	int ntsnap;

	//model
	char modelname[STRLEN];

	//slide
	float before;
	float after;
	char timefile[STRLEN];

	void read_param(int ac, char**av, const char* param_type);
	void make_plan();
};
#endif
