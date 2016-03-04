#include"util.h"
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<fstream>
void zap(double *x, int n)
{
  int i;
  for(i=0; i<n; i++) x[i]= 0.0;
}

void read_model(double *field, int size, char *model,const char *ext)
{
  char name[128];

  sprintf(name,"%s.%s",model, ext);

  std::ifstream fd;
  fd.open(name,std::ios::binary);

  float *tmp = new float[size];
  fd.read((char*)tmp,sizeof(float)*size);
  int i;
  for(i=0;i<size;i++) field[i]=(double)tmp[i];
  
  //fd.read((char*)field,sizeof(double)*size);
  fd.close();
}

int extension(char *name, const char *ext)
{
  /* Return 1 if name has final extension = 'ext', or 0 if not */
  char *p;
  for(p= name; *p != '\0'; p++);
  for( ;p >= name; p--) if(*p == '.') break;
  if(*p != '.') return(0);
  p++;
  if(strcmp(p,ext) == 0) return(1);
  return(0);
}
