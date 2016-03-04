#include"getpar.h"

int main(int argc, char*argv[])
{
  int a;
  setpar(argc,argv);                                                                       mstpar("a","d",&a);  
  getpar("a","d",&a);
  endpar();
  fprintf(stdout,"a=%d\n",a);
}
