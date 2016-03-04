#ifndef _FIELD_H_
#define _FIELD_H_
class PARAM;
class FIELD{
  public:
  float * Tx;
  float * Tz;
  float * V;
  int nx;
  int nz;

  void init_for_full(PARAM & param);
  void init_gpu_full(int deviceid,PARAM & param);
};
#endif
