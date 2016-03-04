#ifndef _FD_H_
#define _FD_H_

class FIELD;
class MATERIAL;
class CPML;
class PARAM;
class SLIDE;

void step_stress(FIELD & fld,MATERIAL & mat,	CPML & cpml);
void step_velocity(FIELD & fld,MATERIAL & mat,	CPML & cpml);

void cu_step_stress(int deviceid,PARAM &param,FIELD & fld,MATERIAL & mat,	CPML & cpml,SLIDE&slide);
void cu_step_velocity(int deviceid,PARAM &param,FIELD & fld,MATERIAL & mat,	CPML & cpml,SLIDE& slide);
#endif
