#include"material.h"
#include<cstdio>
#include<algorithm>

inline int inside_table(MATERIAL & mat,float BV,float MUX,float MUZ)
{
  int check_remember_max=10000; // only check the last 10000 material to save time
  int start_check=mat.num_mat-check_remember_max>0?mat.num_mat-check_remember_max:0;

  for(int i=start_check;i<mat.num_mat;i++){
	if (BV==mat.tbl_BV[i] && 
		MUX==mat.tbl_MUX[i] && MUZ==mat.tbl_MUZ[i])
	  return i;
  }
  return -1;
}
inline int insert_table(MATERIAL & mat,float BV,float MUX,float MUZ)
{
  int n=mat.num_mat;
  mat.tbl_BV[n]=BV;
  mat.tbl_MUX[n]=MUX;
  mat.tbl_MUZ[n]=MUZ;
  mat.num_mat+=1;
}
inline void resize(int n,float * & a,float* temp)
{
  std::copy(a,a+n,temp);
  delete [] a;
  a=new float[n];
  std::copy(temp,temp+n,a);
}
inline void resize(int n,int * & a,int* temp)
{
  std::copy(a,a+n,temp);
  delete [] a;
  a=new int[n];
  std::copy(temp,temp+n,a);
}
void MATERIAL::mktable()
{
  tbl_BV=new float[nx*nz];
  tbl_MUX=new float[nx*nz];
  tbl_MUZ=new float[nx*nz];
  index=new int[nx*nz];

  int TABLE_MAX=(nx*nz)*0.01;


  num_mat=0;
  for (int iz=0;iz<nz;iz++){
	for(int ix=0;ix<nx;ix++){
	  int ind=iz*nx+ix;
	  index[ind]=inside_table(*this,BV[ind],MUX[ind],MUZ[ind]);
	  if(index[ind]<0){
		index[ind]=num_mat;
		insert_table(*this,BV[ind],MUX[ind],MUZ[ind]);
	  }
	  if(num_mat == TABLE_MAX )
	  {
		printf("Too many materials, abort make table\n");
		return;
	  }
	}
//	printf("%5.2f percent done",(float)iz/(float)nz*100);
	printf("%5.2f%% done %06d",(float)iz/(float)nz*100,num_mat);
	printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
  }

  //resize array
  float *temp=new float[num_mat];
  resize(num_mat,tbl_BV,temp);
  resize(num_mat,tbl_MUX,temp);
  resize(num_mat,tbl_MUZ,temp);
  delete [] temp;
  usetable=true;
  printf("\n");

   //check table
 for(int i=0;i<nx*nz;i++){
   if(  BV[i]!=  tbl_BV[index[i]]
     || MUX[i]!=  tbl_MUX[index[i]]
     || MUZ[i]!= tbl_MUZ[index[i]]
       ){
     printf("error %07d: %16.9e %16.9e",i,BV[i],tbl_BV[index[i]]);
     printf("error %07d: %16.9e %16.9e",i,MUX[i],tbl_MUX[index[i]]);
     printf("error %07d: %16.9e %16.9e",i,MUZ[i],tbl_MUZ[index[i]]);
   }
 }
 printf("check passed\n");

  printf("---------material table-------\n");
  for(int i=0;i<num_mat;i++){
	printf("%05d %16.9e %16.9e %16.9e\n",i,tbl_BV[i],tbl_MUX[i],tbl_MUZ[i]);
  }
}
