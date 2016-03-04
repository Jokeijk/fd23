#include"field/field.h"
#include"material/material.h"
#include"cpml/cpml.h"
#include"fd/fd.h"
#include"record/record.h"
#include"record/isis.h"
#include"thread/thread.h"
#include"slide/slide.h"
#include"timer.h"
#include"point_source/point_source.h"
#include"param/const.h"
#include"point_source/waveform.h"
#include<iostream>
#include<iomanip>




using namespace std;

struct ARGUMENT{
  int deviceid;

  THREAD * thread;
  RECORD * record_v;
  PARAM * param;

  CPML *cpml;
  CPML **d_cpml;

  FIELD * g_fld;
  MATERIAL *mat;

  FIELD ** d_fld;
  MATERIAL ** d_mat;
};

void * singlegpu(void * argument)
{
  struct ARGUMENT * arg=(struct ARGUMENT *)argument;

  int deviceid       = arg->deviceid;

  THREAD * thread    = arg->thread;
  RECORD * record_v  = arg->record_v;
  PARAM * param      = arg->param;

  CPML *cpml         = arg->cpml;
  CPML **d_cpml      = arg->d_cpml;

  FIELD *g_fld       = arg->g_fld;
  MATERIAL *mat      = arg->mat;


  FIELD ** d_fld     = arg->d_fld;
  MATERIAL ** d_mat  = arg->d_mat;


  thread->selectdevice(deviceid);
  d_fld[deviceid]=new FIELD();
  d_fld[deviceid]->init_gpu_full(deviceid, *param);

  d_mat[deviceid]=new MATERIAL();
  d_mat[deviceid]->init_gpu_full(deviceid, *param, *mat);

  d_cpml[deviceid]=new CPML(deviceid,*cpml);


  record_v->cu_init_record(deviceid,*param);

  SLIDE slide(*param);

  TIME time;

  /*Set up source time*/
  float *src;
  int lsrc;
  if(param->sourcetime[0]=='T' || param->sourcetime[0]=='t'){
	trapozoid(param->dt,param->trap1,param->trap2,param->trap3,lsrc,src);
  }else if(param->sourcetime[0]=='G' || param->sourcetime[0]=='g'){
	gaussian(param->dt,param->alpha,lsrc,src);
  }

  float *psrc=new float[lsrc];
  float sum=0;
  for(int i=0;i<lsrc;i++){
	sum+=src[i];
	psrc[i]=sum;
  }



  for(int it=0; it < param->nt; it++){

	slide.current_step=it;

	cu_step_stress(deviceid,*param,*d_fld[deviceid] ,*d_mat[deviceid] ,*d_cpml[deviceid],slide);
	thread->exchange_stress(deviceid,*param,*d_fld[deviceid],*g_fld);

	switch (param->stype)
	{
	  case 3:
		cu_point_source(deviceid,it,lsrc,src,*param,*d_fld[deviceid]);
		break;
	  case 13:
		cu_point_source_p(deviceid,it,lsrc,psrc,*param,*d_fld[deviceid]);
		break;
	  default:
		cu_point_source(deviceid,it,lsrc,src,*param,*d_fld[deviceid]);
        //        cu_point_force(deviceid,it,lsrc,src,*param,*d_fld[deviceid]);

	}
	cu_step_velocity(deviceid,*param,*d_fld[deviceid], *d_mat[deviceid],*d_cpml[deviceid],slide);
	if(it%param->itrecord==0){
	  //record this gpu portion, and wait for syncronize thread to flush
	  record_v->cu_record(deviceid,param->nx,d_fld[deviceid]->V);
	}
	thread->exchange_velocity(deviceid,*param,*d_fld[deviceid],*g_fld);
	if(it%param->itrecord==0){
	  //flush now,since there is also a barrier at exchange stress,so this is safe
	  record_v->cu_flush(deviceid);
	}


	if( ( it%param->ntsnap==0 || needsnap(*param,it) ) && it!=0){
	  cout<<it<<endl;
	  if(param->ngpu>0){
		cpy_backhost(deviceid,*param,*g_fld,*d_fld[deviceid]);
		thread->wait();
	  }
	  if(deviceid==0){
		snapshot(*g_fld,it);
	  }
	}
	if(deviceid==0 && it%param->itprint==0){
	  cout<<setw(10)<<it<<" TOTAL "
		<<setw(10)<<time.elips()/1000000.0<<" sec" 
		<<setw(10)<<time.elips()/(it+1)/1000.0<<" msec per step\n";
	  cout.flush();
	}
  }
}

int main(int ac, char **av)
{
  PARAM    param;
  param.read_param(ac,av,"full");
  param.make_plan();

  MATERIAL mat;
  mat.init_for_full(param);
  mat.get_model(param);

  mat.mktable();


  FIELD    g_fld;
  g_fld.init_for_full( param );

  CPML    cpml(param);

  FIELD * d_fld[MAXCUDAGPU];
  MATERIAL * d_mat[MAXCUDAGPU];
  CPML * d_cpml[MAXCUDAGPU];

  struct ARGUMENT arg[MAXCUDAGPU];

  THREAD thread(param.ngpu);


  char filename[256];
  RECORD record_v;

  sprintf(filename,"%s_V",param.output);
  if(param.use_sta_file){
	record_v.init_list(param.V_sta_file,filename);
  }else{
	record_v.init_line(param.nrec,param.ixrec0,param.izrec0_v,param.idxrec,param.idzrec,filename);
  }



  for(int i=0;i<param.ngpu;i++){
	arg[i].deviceid    =  i;
	arg[i].thread      = &thread;
	arg[i].record_v    = &record_v;
	arg[i].param       = &param;

	arg[i].cpml        = &cpml;
	arg[i].d_cpml      = d_cpml;

	arg[i].g_fld       = &g_fld;
	arg[i].mat         = &mat;

	arg[i].d_fld       = d_fld;
	arg[i].d_mat       = d_mat;

  }

  if(param.ngpu<=0){
	fprintf(stderr,"To be implement on non-GPU device!\n");
	return -1;
  }else{
	for(int i=0;i<param.ngpu;i++){
	  if(pthread_create(&thread.thr[i],NULL,&singlegpu,(void*)&arg[i]))
	  {
		fprintf(stderr,"failure to create thread\n");
	  }
	}
	for(int i=0;i<param.ngpu;i++){
	  if(pthread_join(thread.thr[i],NULL))
	  {
		fprintf(stderr,"failure to join thread\n");
	  }
	}
  }
  record_v.close();
  record_v.toisis(param,0.0,0.0);
  cout<<"finish\n";
}
