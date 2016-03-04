/* Author Dunzhu Li  dli@caltech.edu 
*/
#ifndef _CONST_H_
#define _CONST_H_

#define C1       9./8. 
#define C2       1./24. 
#define D1      75./64. 
#define D2      25./384. 
#define D3       3./640. 
#define E1    1225./1024. 
#define E2     245./3072. 
#define E3      49./5120. 
#define E4       5./7168.
const float g_coef[5][4]=
{ 1.0, 0.0, 0.0, 0.0,
  1.0, 0.0, 0.0, 0.0,
   C1,  C2, 0.0, 0.0,
   D1,  D2,  D3, 0.0,
   E1,  E2,  E3, E4 }; 
#define ORD8 4
#define BDIMX 32
#define BDIMY 16
#define radius 4
#define MAXCUDAGPU 3

#define STRLEN 256

#define PI 3.14159265358979
#endif
