#ifndef _CPML_H_
#define _CPML_H_
class PARAM;
class PMLVARIABLE{
  public:
	double* Txx_x;
	double* Txz_x;
	double*   U_x;
	double*   W_x;

	double* Tzz_z;
	double* Txz_z;
	double*   U_z;
	double*   W_z;
};
class CPML {
  public:
	CPML(PARAM & param);
	CPML(int deviceid,CPML &cpml);
	PMLVARIABLE psi,b,c,k;
	int npml;
	int nx;
	int nz;
	double pml_dt;
	double pml_r;
	double pml_v;
	double pml_fc;
};
#endif
