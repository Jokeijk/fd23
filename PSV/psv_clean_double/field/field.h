#ifndef _FIELD_H_
#define _FIELD_H_
class PARAM;
class FIELD{
  public:
  double * Txx;
  double * Txz;
  double * Tzz;
  double * U;
  double * W;
  int nx;
  int nz;

  void init_for_box(PARAM &param);
  void init_for_full(PARAM & param);

  void init_gpu_full(int deviceid,PARAM & param);
  void init_gpu_box (int deviceid,PARAM & param);
};
#endif
