#ifndef _PARAM_H_
#define _PARAM_H_

#include "const.h"

class PARAM{
  public:
	//flags
	int usebox;
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
	double	h;
	double	dt;

	//earthquake
	double strike;
	double dip;
	double rake;
	double azimuth;

	// source box

	int   stype;
	double srcden;
	double srcsvel;
	double srcpvel;
	char  sourcefile[STRLEN];
	int   bw;

	double src_shift_x;
	double src_shift_z;

	//source time
	char sourcetime[10];
	double trap1;
	double trap2;
	double trap3;
	double alpha;

	//source
	int xs;
	int zs;


	//for pml
	int   npml;
	double pml_r;
	double pml_v;
	double pml_fc;
	double pml_dt;

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
	double before;
	double after;
	char timefile[STRLEN];

	void read_param(int ac, char**av, const char* param_type);
	void make_plan();
};
#endif
