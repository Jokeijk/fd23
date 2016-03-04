#include"field.h"
#include"../param/param.h"
#include<algorithm>

#define ALLOC_FIELD  \
  Txx=new float[nx*(nz+8)];\
  Txz=new float[nx*(nz+8)];\
  Tzz=new float[nx*(nz+8)];\
  U=new float[nx*(nz+8)];\
  W=new float[nx*(nz+8)];\
  std::fill_n(Txx,nx*(nz+8),0);\
  std::fill_n(Txz,nx*(nz+8),0);\
  std::fill_n(Tzz,nx*(nz+8),0);\
  std::fill_n(U,nx*(nz+8),0);\
  std::fill_n(W,nx*(nz+8),0);\
  Txx+=nx*4;\
  Txz+=nx*4;\
  Tzz+=nx*4;\
  U+=nx*4;\
  W+=nx*4;

void FIELD::init_for_full(PARAM & param){
  nx=param.nx;
  nz=param.nz;

  ALLOC_FIELD;
}
