#ifndef _GETPAR_H_
#define _GETPAR_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

int getpar(char *name,char *type,const void *val);
int mstpar(char *name,char *type,const void *val);
int setpar(int ac,char **av);
void endpar();

#ifdef __cplusplus
}
#endif

#endif /* _GETPAR_H_ */

