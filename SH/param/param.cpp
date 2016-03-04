#include "param.h"
#include "../getpar/getpar.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>

void PARAM::read_param(int ac, char **av, const char* param_level)
{
  
  setpar(ac,av);

  output_material=0;
  getpar("output_material" ,"d",&output_material);

  mstpar("usepunch","d",&usepunch);
  mstpar("ngpu","d",&ngpu);

  mstpar("nx","d",&nx);
  mstpar("nz","d",&nz);
  mstpar("nt","d",&nt);
  mstpar("h","f",&h);
  mstpar("dt","f",&dt);


  mstpar("sourcetime","s",sourcetime);
  if(sourcetime[0]=='T' || sourcetime[0]=='t'){
	mstpar("trap1","f",&trap1);
	mstpar("trap2","f",&trap2);
	mstpar("trap3","f",&trap3);
  };
  if(sourcetime[0]=='G' || sourcetime[0]=='g'){
	mstpar("alpha","f",&alpha);
  }

  stype=3;
  getpar("stype" ,"d",&stype);
  mstpar("strike" ,"f",&strike);
  mstpar("dip"    ,"f",&dip);
  mstpar("rake"   ,"f",&rake);
  mstpar("azimuth","f",&azimuth);

  mstpar("xs","d",&xs);
  mstpar("zs","d",&zs);

  itprint=100;
  getpar("itprint","d",&itprint);



  mstpar("npml","d",&npml);
  mstpar("pml_dt","f",&pml_dt);
  mstpar("pml_v","f",&pml_v);
  mstpar("pml_r","f",&pml_r);
  mstpar("pml_fc","f",&pml_fc);


  mstpar("itrecord","d",&itrecord);
  mstpar("output"  ,"s",output);

  V_sta_file[0]=0;
  use_sta_file=true;
  getpar("V_sta_file","s",V_sta_file);

  if(V_sta_file[0]==0){
	use_sta_file=false;
    mstpar("nrec",    "d",&nrec);
	mstpar("ixrec0",  "d",&ixrec0);
	mstpar("izrec0_v","d",&izrec0_v);
	mstpar("idxrec"  ,"d",&idxrec);
	mstpar("idzrec"  ,"d",&idzrec);
  }

  mstpar("model","s",modelname);

  mstpar("ntsnap","d",&ntsnap);
  ntsnapfile[0]=0;
  getpar("ntsnapfile","s",ntsnapfile);


  if(usepunch==1){
	mstpar("before","f",&before);
	mstpar("after","f",&after);
	mstpar("timefile","s",timefile);
  }
  endpar();
}




void PARAM::make_plan()
{
  boxdevice=-1;
  if(ngpu<=0) return;

  int dn,i;
  if( nx%BDIMX !=0 || nz%BDIMY !=0 ){
	fprintf(stderr,"nx ny must be multiple of %d,%d\n",BDIMX,BDIMY);
	exit(-1);
  }

  dn=nz/BDIMY/ngpu;

  for(i=0;i<ngpu;i++){
	a_nz1[i]=dn*i*BDIMY;
	a_nz2[i]=dn*(i+1)*BDIMY-1;
  }
  a_nz2[ngpu-1]=nz-1;

  for(i=0;i<ngpu;i++){
	printf("nz1,nz2:%d,%d\n",a_nz1[i],a_nz2[i]);
  }

  for(i=0;i<ngpu;i++){
	if(a_nz2[i]-zs>=0){
	  boxdevice=i;
	  printf("boxdevice= %d\n",boxdevice);
	  break;
	}
  }
}

