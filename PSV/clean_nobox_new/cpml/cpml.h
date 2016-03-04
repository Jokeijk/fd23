#ifndef _CPML_H_
#define _CPML_H_
class PARAM;
class PMLVARIABLE{
  public:
	float* Txx_x;
	float* Txz_x;
	float*   U_x;
	float*   W_x;

	float* Tzz_z;
	float* Txz_z;
	float*   U_z;
	float*   W_z;
};
class CPML {
  public:
	CPML(PARAM & param);
	CPML(int deviceid,CPML &cpml);
	PMLVARIABLE psi,b,c,k;
	int npml;
	int nx;
	int nz;
	float pml_dt;
	float pml_r;
	float pml_v;
	float pml_fc;

	void cu_load_restart(int deviceid,CPML &cpml);
	void cu_save_state(int deviceid,CPML &cpml);
};
#endif
