#ifndef _UTIL_H_
#define _UTIL_H_
int extension(char *name, const char *ext);
void read_model(double *field, int size, char *model, const char *ext);
void zap(double *x,int n);
#endif
