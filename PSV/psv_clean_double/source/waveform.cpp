#include"waveform.h"
#include<cmath>
#include<cstdlib>
/*
 * int s(t) dt = 1
 * We compute is s(t)dt , not s(t)
 * */
void trapozoid(double dt,double trap1,double trap2,double trap3,int &lsrc,double * & src)
{
  int t1 = (int) (trap1 / dt);
  int t2 = (int) (trap2 / dt);
  int t3 = (int) (trap3 / dt);
  int i;

  t2 += t1;
  lsrc = t2 + t3;
  src = new double[lsrc + 1];

  double amp = (trap2 + 0.5 * (trap1 + trap3)) / dt;

  double *pxx = src;
  for (i = 0 ; i < t1; i++) *pxx++ = i / ( t1 * amp);
  for (      ; i < t2; i++) *pxx++ = 1.0 / amp;
  for (i = t3; i > 0 ; i--) *pxx++ = i / ( t3 * amp);

}
void gaussian(double dt,int alpha ,int &lsrc,double * & src)
{
  /* gaussian time function */
  /* alpha<0, normalized one-sided gaussian */
  /* alpha>0, two-sided gaussian */
  /*     fprintf(stderr, "trap => 0.0\n"); */
  int aint = abs((int)(alpha * 3.0));
  lsrc = 2 * aint + 1;
  src = new double[lsrc + 1];
  if (alpha >= 0.0)
    for (int i = -aint; i <= aint; i++)
      src[i+aint] = (double)i * exp( double(-i) * double(i) / (alpha * alpha));
  else {
    double sum = 0.0;
    for (int i = -aint; i <= aint; i++) {
      src[i+aint] = exp( double (-i) * double(i) / (alpha * alpha));
      sum += src[i+aint];
    }
    for (int i = 0; i <= lsrc; i++) src[i] /= sum;
  }
}

//void trapozoid(double dt,double trap1,double trap2,double trap3,int &lsrc,double * & src)
//{
//  int t1 = (int) (trap1 / dt);
//  int t2 = (int) (trap2 / dt);
//  int t3 = (int) (trap3 / dt);
//  int i;
//
//  t2 += t1;
//  lsrc = t2 + t3;
//  src = new double[lsrc + 1];
//
//  double amp = (trap2 + 0.5 * (trap1 + trap3)) / dt;
//
//  double *pxx = src;
//  for (i = 0 ; i < t1; i++) *pxx++ = i / ( t1 * amp);
//  for (      ; i < t2; i++) *pxx++ = 1.0 / amp;
//  for (i = t3; i > 0 ; i--) *pxx++ = i / ( t3 * amp);
//
//}
//void gaussian(double dt,int alpha ,int &lsrc,double * & src)
//{
//  /* gaussian time function */
//  /* alpha<0, normalized one-sided gaussian */
//  /* alpha>0, two-sided gaussian */
//  /*     fprintf(stderr, "trap => 0.0\n"); */
//  int aint = abs((int)(alpha * 3.0));
//  lsrc = 2 * aint + 1;
//  src = new double[lsrc + 1];
//  if (alpha >= 0.0)
//    for (int i = -aint; i <= aint; i++)
//      src[i+aint] = (double)i * exp( double(-i) * double(i) / (alpha * alpha));
//  else {
//    double sum = 0.0;
//    for (int i = -aint; i <= aint; i++) {
//      src[i+aint] = exp( double (-i) * double(i) / (alpha * alpha));
//      sum += src[i+aint];
//    }
//    for (int i = 0; i <= lsrc; i++) src[i] /= sum;
//  }
//}
