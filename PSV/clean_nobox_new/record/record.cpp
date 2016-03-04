#include<fstream>
#include"record.h"
#include"../material/material.h"
#include"../field/field.h"
using namespace std;
void snapshot(MATERIAL & fld){
  int nx=fld.nx;
  int nz=fld.nz;
  ofstream snap;
  char snapfile[256];

  sprintf(snapfile,"BU.snap");
  snap.open(snapfile,ios::binary);
//  snap.write((char*)&nx,4);
//  snap.write((char*)&nz,4);
  snap.write((char*)fld.BU,4*nx*nz);
  snap.close();

  sprintf(snapfile,"MU.snap");
  snap.open(snapfile,ios::binary);
//  snap.write((char*)&nx,4);
//  snap.write((char*)&nz,4);
  snap.write((char*)fld.BU,4*nx*nz);
  snap.close();
}

void snapshot(FIELD & fld,int it,int*snapbox){
  int nx1=snapbox[0];
  int nx2=snapbox[1];
  int nz1=snapbox[2];
  int nz2=snapbox[3];

  int nx=fld.nx;
  int nz=fld.nz;

  ofstream snap;
  char snapfile[256];

  sprintf(snapfile,"./OUTPUT/U%05d",it);
  snap.open(snapfile,ios::binary);
  for(int i=nz1;i<=nz2;i++){
	snap.write((char*)(fld.U+i*nx+nx1),4*(nx2-nx1+1));
  }
  snap.close();

  sprintf(snapfile,"./OUTPUT/W%05d",it);
  snap.open(snapfile,ios::binary);
  for(int i=nz1;i<=nz2;i++){
	snap.write((char*)(fld.W+i*nx+nx1),4*(nx2-nx1+1));
  }
  snap.close();
}
void snapshot(FIELD & fld,int it){
  int nx=fld.nx;
  int nz=fld.nz;
  ofstream snap;
  char snapfile[256];

  sprintf(snapfile,"./OUTPUT/U%05d",it);
  snap.open(snapfile,ios::binary);
//  snap.write((char*)&nx,4);
//  snap.write((char*)&nz,4);
  snap.write((char*)fld.U,4*nx*nz);
  snap.close();

  sprintf(snapfile,"./OUTPUT/W%05d",it);
  snap.open(snapfile,ios::binary);
//  snap.write((char*)&nx,4);
//  snap.write((char*)&nz,4);
  snap.write((char*)fld.W,4*nx*nz);
  snap.close();
/*

  sprintf(snapfile,"Txx%05d",it);
  snap.open(snapfile,ios::binary);
  snap.write((char*)&nx,4);
  snap.write((char*)&nz,4);
  snap.write((char*)fld.Txx,4*nx*nz);
  snap.close();

  sprintf(snapfile,"Txz%05d",it);
  snap.open(snapfile,ios::binary);
  snap.write((char*)&nx,4);
  snap.write((char*)&nz,4);
  snap.write((char*)fld.Txz,4*nx*nz);
  snap.close();

  sprintf(snapfile,"Tzz%05d",it);
  snap.open(snapfile,ios::binary);
  snap.write((char*)&nx,4);
  snap.write((char*)&nz,4);
  snap.write((char*)fld.Tzz,4*nx*nz);
  snap.close();
*/
}
