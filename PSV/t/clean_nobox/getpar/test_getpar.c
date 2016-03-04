#include"getpar.h"

int main(int argc, char*argv[])
{
  int a[3];
  setpar(argc,argv);  
  mstpar("a","d",&a);  
  endpar();
  fprintf(stdout,"a=%d\n",a[1]);
}
