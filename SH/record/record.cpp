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

  sprintf(snapfile,"BV.snap");
  snap.open(snapfile,ios::binary);
//  snap.write((char*)&nx,4);
//  snap.write((char*)&nz,4);
  snap.write((char*)fld.BV,4*nx*nz);
  snap.close();

  sprintf(snapfile,"MUX.snap");
  snap.open(snapfile,ios::binary);
//  snap.write((char*)&nx,4);
//  snap.write((char*)&nz,4);
  snap.write((char*)fld.MUX,4*nx*nz);
  snap.close();

  sprintf(snapfile,"MUZ.snap");
  snap.open(snapfile,ios::binary);
//  snap.write((char*)&nx,4);
//  snap.write((char*)&nz,4);
  snap.write((char*)fld.MUZ,4*nx*nz);
  snap.close();
}
void snapshot(FIELD & fld,int it){
  int nx=fld.nx;
  int nz=fld.nz;
  ofstream snap;
  char snapfile[256];

  sprintf(snapfile,"V%05d",it);
  snap.open(snapfile,ios::binary);
//  snap.write((char*)&nx,4);
//  snap.write((char*)&nz,4);
  snap.write((char*)fld.V,4*nx*nz);
  snap.close();
}
