#include"material.h"
#include<cstdio>
#include<algorithm>

inline int inside_table(MATERIAL & mat,double BU,double BW,double MU,double MUA,double LAM)
{
  int check_remember_max=10000; // only check the last 10000 material to save time
  int start_check=mat.num_mat-check_remember_max>0?mat.num_mat-check_remember_max:0;

  for(int i=start_check;i<mat.num_mat;i++){
	if (BU==mat.tbl_BU[i] && BW==mat.tbl_BW[i] && 
		MU==mat.tbl_MU[i] && MUA==mat.tbl_MUA[i] && 
		LAM==mat.tbl_LAM[i])
	  return i;
  }
  return -1;
}
inline int insert_table(MATERIAL & mat,double BU,double BW,double MU,double MUA,double LAM)
{
  int n=mat.num_mat;
  mat.tbl_BU[n]=BU;
  mat.tbl_BW[n]=BW;
  mat.tbl_MU[n]=MU;
  mat.tbl_MUA[n]=MUA;
  mat.tbl_LAM[n]=LAM;
  mat.num_mat+=1;
}
inline void resize(int n,double * & a,double* temp)
{
  std::copy(a,a+n,temp);
  delete [] a;
  a=new double[n];
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
  tbl_BU=new double[nx*nz];
  tbl_BW=new double[nx*nz];
  tbl_MU=new double[nx*nz];
  tbl_MUA=new double[nx*nz];
  tbl_LAM=new double[nx*nz];
  index=new int[nx*nz];

  int TABLE_MAX=(nx*nz)*0.01;


  num_mat=0;
  for (int iz=0;iz<nz;iz++){
	for(int ix=0;ix<nx;ix++){
	  int ind=iz*nx+ix;
	  index[ind]=inside_table(*this,BU[ind],BW[ind],MU[ind],MUA[ind],LAM[ind]);
	  if(index[ind]<0){
		index[ind]=num_mat;
		insert_table(*this,BU[ind],BW[ind],MU[ind],MUA[ind],LAM[ind]);
	  }
	  if(num_mat == TABLE_MAX )
	  {
		printf("Too many materials, abort make table\n");
		return;
	  }
	}
//	printf("%5.2f percent done",(double)iz/(double)nz*100);
	printf("%5.2f%% done %06d",(double)iz/(double)nz*100,num_mat);
	printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
  }

  //resize array
  double *temp=new double[num_mat];
  resize(num_mat,tbl_BU,temp);
  resize(num_mat,tbl_BW,temp);
  resize(num_mat,tbl_MU,temp);
  resize(num_mat,tbl_MUA,temp);
  resize(num_mat,tbl_LAM,temp);
  delete [] temp;
  usetable=true;
  printf("\n");

   //check table
 for(int i=0;i<nx*nz;i++){
   if(  BU[i]!=  tbl_BU[index[i]]
     || BW[i]!=  tbl_BW[index[i]]
     || MU[i]!=  tbl_MU[index[i]]
     || MUA[i]!= tbl_MUA[index[i]]
     || LAM[i]!= tbl_LAM[index[i]]
       ){
     printf("error %07d: %16.9e %16.9e",i,BU[i],tbl_BU[index[i]]);
     printf("error %07d: %16.9e %16.9e",i,BW[i],tbl_BW[index[i]]);
     printf("error %07d: %16.9e %16.9e",i,MU[i],tbl_MU[index[i]]);
     printf("error %07d: %16.9e %16.9e",i,MUA[i],tbl_MUA[index[i]]);
     printf("error %07d: %16.9e %16.9e",i,LAM[i],tbl_LAM[index[i]]);
   }
 }
 printf("check passed\n");

  printf("---------material table-------\n");
  for(int i=0;i<num_mat;i++){
	printf("%05d %16.9e %16.9e %16.9e %16.9e %16.9e\n",i,tbl_BU[i],tbl_BW[i],tbl_MU[i],tbl_MUA[i],tbl_LAM[i]);
  }
}
