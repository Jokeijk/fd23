#include"../material/material.h"
#include"../param/param.h"

int main(int ac, char **av)
{
  PARAM  param;
  param.read_param(ac,av,"Output_material");

  MATERIAL mat;
  mat.init_for_full(param);

  mat.get_model(param);
}
