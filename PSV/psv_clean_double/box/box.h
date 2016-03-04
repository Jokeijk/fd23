#ifndef _BOX_H_
#define _BOX_H_
#include<fstream>

class PARAM;
class MATERIAL;
class FIELD;

class BOX{
  public:
	int * r_add_idx;
	int * r_sub_idx;
	int * a_add_idx;
	int * a_sub_idx;
	int  n_add;
	int  n_sub;

	FIELD * p_sfield;
	MATERIAL * p_bmat;

	int bw;

    void init_box(PARAM & param,MATERIAL & mat);
    void subsource(double *a,double* r,double* s);
    void addsource(double *a,double* r,double* s);


	void init_gpu(int deviceid,BOX &box,PARAM &param);
    void cu_subsource(int deviceid,double *a,double* r,double* s);
    void cu_addsource(int deviceid,double *a,double* r,double* s);

	void setboxsource();
	void cu_copy_source(int deviceid,BOX &box);

	int start_x; //box start index in global 
	int start_z;

	int boxdevice;

  private:
	std::ifstream sourcefile;
	double *u1,*u2,*w1,*w2;
	double *uold,*unew,*wold,*wnew;
	bool u1_read_new;
	double dt,h,mu,gam,lam,srcsvel,srcpvel,srcden;

}; 
#endif
