#ifndef _MODEL_H_
#define _MODEL_H_

int load_model(char*);
/*x0,z0 is the first point position of the grid, the other grid point position is compute from x=x0+ix*h*/
void get_elas_model( float *vp,	float *vs,float *den,int nx,int nz,	double h,double x0,	double z0);

#endif
