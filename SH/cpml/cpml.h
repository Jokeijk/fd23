#ifndef _CPML_H_
#define _CPML_H_
class PARAM;
class PMLVARIABLE{
  public:
	float* Tx_x;
	float*  V_x;

	float* Tz_z;
	float*  V_z;
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
};
#endif
