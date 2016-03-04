#include"field.h"
#include"../param/param.h"
#include<algorithm>

#define ALLOC_FIELD  \
  V=new float[nx*(nz+8)];\
  Tx=new float[nx*(nz+8)];\
  Tz=new float[nx*(nz+8)];\
  std::fill_n(V,nx*(nz+8),0);\
  std::fill_n(Tx,nx*(nz+8),0);\
  std::fill_n(Tz,nx*(nz+8),0);\
  V+=nx*4;\
  Tx+=nx*4;\
  Tz+=nx*4;

void FIELD::init_for_full(PARAM & param){
  nx=param.nx;
  nz=param.nz;

  ALLOC_FIELD;
}
