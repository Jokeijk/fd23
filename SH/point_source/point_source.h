#ifndef _POINT_SOURCE_H_
#define _POINT_SOURCE_H_
class PARAM;
class FIELD;
/* Set Point Source */
void cu_point_source(int deviceid,int it,int lsrc,float *src, PARAM &param, FIELD & fld);
void cu_point_source_p(int deviceid,int it,int lsrc,float *src, PARAM &param, FIELD & fld);
void cu_point_force(int deviceid,int it,int lsrc,float *src, PARAM &param, FIELD & fld);
/*cu_point_force was added by Kangchen*/
#endif
