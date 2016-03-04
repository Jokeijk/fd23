#ifndef _FIELD_H_
#define _FIELD_H_
class PARAM;
class FIELD{
  public:
  float * Txx;
  float * Txz;
  float * Tzz;
  float * U;
  float * W;
  int nx;
  int nz;

  void init_for_full(PARAM & param);
  void init_gpu_full(int deviceid,PARAM & param);

  void cu_load_restart(int deviceid,PARAM &param,FIELD & g_fld);
  void cu_save_state(int deviceid,PARAM &param,FIELD & g_fld);
};
#endif
