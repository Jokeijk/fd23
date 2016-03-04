#include "param.h"
#include "../getpar/getpar.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>

void usage()
{

  char *sdoc[]= {
" restart=0          #  if it's restart run               ",
"     nt0=0          #  the first new time step, will load time step nt0-1   ",
" nt=                #                           ",
" savestate=0        #  save the last step (nt0+nt-1)  ",
" gpuid=1            #  gpuids [0,1,2], default use only one gpu(id is 1), use gpuid 0 will slow the desktop display",
" nx=                #  x direction, will be multiple of 32",
" nz=                #  z direction, will be multiple of 16",
" xs=                #  source position, the souce located on Txx,Tzz grid",
" zs=                #  source position",
" h=                 #  grid size ",
" dt=                #  time step size",
" srctype=1          #  1 for 2D line, p for p",
" srctime=           #  t for triangle, g for gaussian",
"  trap1=            #  trap1,trap2,trap3 for triangle",
"  trap2=            #   ",
"  trap3=            #   ",
"  alpha=            #  <0 for gaussian, bell length is about 6*alpha ",
"                    #      ",
" strike=            #   ",
" dip=               #   ",
" rake=              #     ",
" azimuth=           #   ",
"                    #      ",
" itprint=100        #     ",

" npml=32            #   pml need to be toned if you want good result",
" pml_dt=0.005       #   pml dt, similar to dt",
" pml_v=30.0         #   pml v, similar to vp",
" pml_r=1E-11        #   pml r, target reflection coefficient",
" pml_fc=2           #   pml frequency",

" itrecord=1         #   record every itrecord step",
" output=            #   tag name, output_U.isis",
" U_sta_file=        #   if specified, read receiver (z,x) pairs from file",
" W_sta_file=        #   if specified, read receiver (z,x) pairs from file",
"                    #       ",
" nrec=              #   for line record",
" ixrec0=            #   first x0",
" izrec0_u=          #   z0 for u",
" izrec0_w=          #   z0 for w",
" idxrec=            #   line record dx",
" idzrec=            #   line record dz",
"                    #      ",
" model=             #   model tag, read model.vp model.vs model.den",
" usetable=1         #   =1 make table, maybe slow",
" ntsnap=999999999   #   snapshot every ntstep",
" snapbox=0,nx-1,0,nz-1,0,nt-1   #   snapbox ",
" usepunch=0         #   use punch",
"  before=           #    ",
"  after=            #    ",
"  timefile=         #    ",
NULL};

  /* based on cwp */
  char **p = sdoc;
  FILE *fp;
  char *pager;
  char cmd[32];

  if ((pager=getenv("PAGER")) != (char *)NULL)
	sprintf(cmd,"%s 1>&2", pager);
  else 
	sprintf(cmd,"more 1>&2");


  fflush(stdout);
  fp = (FILE *) popen(cmd, "w");
  while(*p) (void)fprintf(fp, "%s\n", *p++);
  pclose(fp);

  exit(-1);
}


void PARAM::read_param(int ac, char **av, const char* param_level)
{

  if(ac==1) usage();


  setpar(ac,av);

  savestate=0;
  restart=0;
  nt0=0;

  getpar("restart","d",&restart);
  if(restart !=0){
	mstpar("nt0","d",&nt0);
  }

  getpar("savestate","d",&savestate);
  
  char strtemp[STRLEN];
  strtemp[0]='\0';
  getpar("gpuid","s",strtemp);
  if(strtemp[0]!='\0'){
	ngpu=sscanf(strtemp,"%d,%d,%d",&gpuid[0],&gpuid[1],&gpuid[2]);
  }else{
	ngpu=1;
	gpuid[0]=1;
  }
   
  
  mstpar("nx","d",&nx);
  mstpar("nz","d",&nz);
  mstpar("nt","d",&nt);
  mstpar("h","f",&h);
  mstpar("dt","f",&dt);


  srctype='1';
  strtemp[0]='\0';
  getpar("srctype","s",strtemp);
  if(strtemp[0]!='\0') srctype=strtemp[0];

  mstpar("srctime","s",srctime);
  
  if(srctime[0]=='T' || srctime[0]=='t'){
	mstpar("trap1","f",&trap1);
	mstpar("trap2","f",&trap2);
	mstpar("trap3","f",&trap3);
  };
  if(srctime[0]=='G' || srctime[0]=='g'){
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

  itrecord=1;
  getpar("itrecord","d",&itrecord);

  mstpar("output"  ,"s",output);

  U_sta_file[0]='\0';
  W_sta_file[0]='\0';
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

  usetable=0;
  getpar("usetable","d",&usetable);

  mstpar("ntsnap","d",&ntsnap);

  strtemp[0]='\0';
  getpar("snapbox","s",strtemp);
  if(strtemp[0]!='\0'){
	if(6!=sscanf(strtemp,"%d,%d,%d,%d,%d,%d",&snapbox[0],&snapbox[1],&snapbox[2],&snapbox[3],&snapbox[4],&snapbox[5])){
	  fprintf(stderr,"snapbox should be six integer");
	  exit(-1);
	}
  }else{
	snapbox[0]=0;
	snapbox[1]=nx-1;
	snapbox[2]=0;
	snapbox[3]=nz-1;
	snapbox[4]=0;
	snapbox[5]=nt;
  }


  usepunch=0;
  getpar("usepunch","d",&usepunch);
  if(usepunch==1){
	mstpar("before","f",&before);
	mstpar("after","f",&after);
	mstpar("timefile","s",timefile);
  }
  endpar();
}




void PARAM::make_plan()
{
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

}


