/* Author Dunzhu Li  dli@caltech.edu 
*/
#ifndef _THREAD_H_
#define _THREAD_H_
#include"../param/const.h"
class PARAM;
class FIELD;
#include<pthread.h>
class THREAD{
  public:
  int exchange_stress  (int deviceid,PARAM &param,FIELD &d_fld,FIELD &g_fld);
  int exchange_velocity(int deviceId,PARAM &param,FIELD &d_fld,FIELD &g_fld);
  void wait();
  THREAD(int ngpu);
  pthread_barrier_t barr;
  pthread_t thr[MAXCUDAGPU];
  void selectdevice(int deviceid);
};
#endif
