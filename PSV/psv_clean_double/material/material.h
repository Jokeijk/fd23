#ifndef _MATERIAL_H_
#define _MATERIAL_H_
class PARAM;

class MATERIAL{
  public:
  double *  BU;
  double *  BW;
  double *  MU;
  double * MUA;
  double * LAM;

  //table, will automatically determine if use table will help
  bool usetable; 
  int * index;
  int num_mat;
  double * tbl_BU;
  double * tbl_BW;
  double * tbl_MU;
  double * tbl_MUA;
  double * tbl_LAM;

  double h;
  double dt;
  double srcden;
  double srcpvel;
  double srcsvel;

  int nx;
  int nz;

  void init_for_full(PARAM & param);
  void init_for_box_from_full(int bw,int start_x,int start_z,MATERIAL &mat);

  void init_gpu_full(int deviceid,PARAM &param,MATERIAL & mat);
  void init_gpu_box(int deviceid,MATERIAL & b_mat);

  void mktable();
  void get_model(PARAM & param);
};

#define  BU(i,j)  BU[i+(j)*nx]
#define  BW(i,j)  BW[i+(j)*nx]
#define  MU(i,j)  MU[i+(j)*nx]
#define MUA(i,j) MUA[i+(j)*nx]
#define LAM(i,j) LAM[i+(j)*nx]

#endif
