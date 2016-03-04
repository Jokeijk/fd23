#include "param.h"
#include "../getpar/getpar.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>

void PARAM::read_param(int ac, char **av, const char* param_level)
{
  output_material=0;
  if(strcmp(param_level,"Output_material")==0){
	printf("Output VP VS DEN\n");
	output_material=1;
	setpar(ac,av);
	mstpar("nx","d",&nx);
	mstpar("nz","d",&nz);
	mstpar("h","f",&h);
	mstpar("model","s",modelname);
	endpar();
	return;
  }

  if(strcmp(param_level,"Source")==0){
	setpar(ac,av);

	mstpar("h","f",&h);
	mstpar("dt","f",&dt);

	mstpar("srcden","f",&srcden);
	mstpar("srcsvel","f",&srcsvel);
	mstpar("srcpvel","f",&srcpvel);
	mstpar("sourcefile","s",sourcefile);
	mstpar("bw","d",&bw);

	src_shift_x=0.25;
	src_shift_z=-0.25;

	getpar("src_shift_x","f",&src_shift_x);
	getpar("src_shift_z","f",&src_shift_z);
	printf("src_shift_x=%f,src_shift_z=%f\n",src_shift_x,src_shift_z);


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
	getpar("stype","d",&stype);

	mstpar("strike" ,"f",&strike);
	mstpar("dip"    ,"f",&dip);
	mstpar("rake"   ,"f",&rake);
	mstpar("azimuth","f",&azimuth);

	endpar();
	return;
  }

  setpar(ac,av);
  mstpar("usebox","d",&usebox);
  mstpar("usepunch","d",&usepunch);
  mstpar("ngpu","d",&ngpu);

  mstpar("nx","d",&nx);
  mstpar("nz","d",&nz);
  mstpar("nt","d",&nt);
  mstpar("h","f",&h);
  mstpar("dt","f",&dt);


  stype=3;
  getpar("stype","d",&stype);

  if(usebox==1){
	mstpar("srcden","f",&srcden);
	mstpar("srcsvel","f",&srcsvel);
	mstpar("srcpvel","f",&srcpvel);
	mstpar("sourcefile","s",sourcefile);
	mstpar("bw","d",&bw);
  }

  src_shift_x=0.25;
  src_shift_z=-0.25;

  getpar("src_shift_x","f",&src_shift_x);
  getpar("src_shift_z","f",&src_shift_z);
  printf("src_shift_x=%f,src_shift_z=%f\n",src_shift_x,src_shift_z);


  mstpar("sourcetime","s",sourcetime);
  if(sourcetime[0]=='T' || sourcetime[0]=='t'){
	mstpar("trap1","f",&trap1);
	mstpar("trap2","f",&trap2);
	mstpar("trap3","f",&trap3);
  };
  if(sourcetime[0]=='G' || sourcetime[0]=='g'){
	mstpar("alpha","f",&alpha);
  }

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

  U_sta_file[0]=0;
  W_sta_file[0]=0;
  use_sta_file=true;
  getpar("U_sta_file","s",U_sta_file);
  getpar("W_sta_file","s",W_sta_file);
  if(U_sta_file[0]==0 && W_sta_file[0]==0){
	use_sta_file=false;
    mstpar("nrec",    "d",&nrec);
	mstpar("ixrec0",  "d",&ixrec0);
	mstpar("izrec0_u","d",&izrec0_u);
	mstpar("izrec0_w","d",&izrec0_w);
	mstpar("idxrec"  ,"d",&idxrec);
	mstpar("idzrec"  ,"d",&idzrec);
  }

  mstpar("model","s",modelname);

  mstpar("ntsnap","d",&ntsnap);



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

  if((zs<bw/2+4 || nz-zs<bw/2+4) && usebox==1){
	fprintf(stderr,"zs %d is too near to the boundary\n",zs);
	exit(-1);
  }

  for(i=0;i<ngpu;i++){
	//adjust to let the box not in the bounday between compute domain
	if(a_nz2[i]-zs<(bw/2+1) && a_nz2[i]-zs> -(bw/2+1))
	{
	  a_nz2[i]  +=((bw/2+1)/BDIMY+1)*BDIMY;
	  a_nz1[i+1]+=((bw/2+1)/BDIMY+1)*BDIMY;
	}
  }

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

