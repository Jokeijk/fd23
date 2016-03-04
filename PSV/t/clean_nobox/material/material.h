/* Author Dunzhu Li  dli@caltech.edu 
*/
#ifndef _MATERIAL_H_
#define _MATERIAL_H_
class PARAM;

class MATERIAL{
  public:
  float *  BU;
  float *  BW;
  float *  MU;
  float * MUA;
  float * LAM;

  //table, will automatically determine if use table will help
  bool usetable; 
  int * index;
  int num_mat;
  float * tbl_BU;
  float * tbl_BW;
  float * tbl_MU;
  float * tbl_MUA;
  float * tbl_LAM;

  float h;
  float dt;

  int nx;
  int nz;

  void init_for_full(PARAM & param);

  void init_gpu_full(int deviceid,PARAM &param,MATERIAL & mat);

  void mktable();
  void get_model(PARAM & param);
};

#define  BU(i,j)  BU[i+(j)*nx]
#define  BW(i,j)  BW[i+(j)*nx]
#define  MU(i,j)  MU[i+(j)*nx]
#define MUA(i,j) MUA[i+(j)*nx]
#define LAM(i,j) LAM[i+(j)*nx]

#endif
