#ifndef _GETPAR_H_
#define _GETPAR_H_

#include <stdio.h>

int getpar(char *name,char *type,const void *val);
int mstpar(char *name,char *type,const void *val);
int setpar(int ac,char **av);
void endpar();

#endif /* _GETPAR_H_ */

