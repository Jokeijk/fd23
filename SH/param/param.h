#ifndef _PARAM_H_
#define _PARAM_H_

#include "const.h"

class PARAM{
  public:
	//flags
	int usepunch;

	// output material
	int output_material;


	//Plan
	int ngpu;
	int a_nz1[MAXCUDAGPU];
	int a_nz2[MAXCUDAGPU];
	int boxdevice;

	//dimension
	int	nx;
	int	nz;
	int	nt;
	float	h;
	float	dt;

	//earthquake
	int stype;
	float strike;
	float dip;
	float rake;
	float azimuth;

	//source time
	char sourcetime[10];
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
	int	izrec0_v;
	int	idxrec;
	int	idzrec;
	char output[STRLEN];
	char V_sta_file[STRLEN];
	bool use_sta_file;

	//info
	int	itprint;
	int ntsnap;
	char ntsnapfile[STRLEN];

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
