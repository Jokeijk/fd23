#ifndef _WAVEFORM_H_
#define _WAVEFORM_H_
void trapozoid(float dt,float trap1,float trap2,float trap3,int &lsrc,float * & src);
void gaussian(float dt,int alpha ,int &lsrc,float * & src);

void trapozoid(double dt,double trap1,double trap2,double trap3,int &lsrc,double * & src);
void gaussian(double dt,int alpha ,int &lsrc,double * & src);
#endif
