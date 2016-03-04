#ifndef _MATERIAL_H_
#define _MATERIAL_H_
class PARAM;

class MATERIAL{
  public:
  float *  BV;
  float * MUZ;
  float * MUX;

  //table, will automatically determine if use table will help
  bool usetable; 
  int * index;
  int num_mat;
  float * tbl_BV;
  float * tbl_MUZ;
  float * tbl_MUX;

  float h;
  float dt;

  int nx;
  int nz;

  void init_for_full(PARAM & param);

  void init_gpu_full(int deviceid,PARAM &param,MATERIAL & mat);

  void mktable();
  void get_model(PARAM & param);
};

#define  BV(i,j)  BV[i+(j)*nx]
#define MUZ(i,j) MUZ[i+(j)*nx]
#define MUX(i,j) MUX[i+(j)*nx]

#endif
