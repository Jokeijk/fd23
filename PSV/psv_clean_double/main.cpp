#include"field/field.h"
#include"material/material.h"
#include"cpml/cpml.h"
#include"box/box.h"
#include"fd/fd.h"
#include"record/record.h"
#include"record/isis.h"
#include"thread/thread.h"
#include"slide/slide.h"
#include"timer.h"
#include"point_source/point_source.h"
#include"param/const.h"
#include"source/waveform.h"
#include<iostream>
#include<iomanip>




using namespace std;

struct ARGUMENT{
  int deviceid;

  THREAD * thread;
  RECORD * record_u;
  RECORD * record_w;
  PARAM * param;

  CPML *cpml;
  CPML **d_cpml;

  FIELD * g_fld;
  MATERIAL *mat;

  FIELD ** d_fld;
  MATERIAL ** d_mat;

  BOX * box;
  BOX * d_box;

  FIELD * d_r_fld;
};

void * singlegpu(void * argument)
{
  struct ARGUMENT * arg=(struct ARGUMENT *)argument;

  int deviceid       = arg->deviceid;

  THREAD * thread    = arg->thread;
  RECORD * record_u  = arg->record_u;
  RECORD * record_w  = arg->record_w;
  PARAM * param      = arg->param;

  CPML *cpml         = arg->cpml;
  CPML **d_cpml      = arg->d_cpml;

  FIELD *g_fld       = arg->g_fld;
  MATERIAL *mat      = arg->mat;


  FIELD ** d_fld     = arg->d_fld;
  MATERIAL ** d_mat  = arg->d_mat;

  BOX * box          = arg->box;
  BOX * d_box        = arg->d_box;

  FIELD * d_r_fld    = arg->d_r_fld;


  thread->selectdevice(deviceid);
  d_fld[deviceid]=new FIELD();
  d_fld[deviceid]->init_gpu_full(deviceid, *param);

  d_mat[deviceid]=new MATERIAL();
  d_mat[deviceid]->init_gpu_full(deviceid, *param, *mat);

  d_cpml[deviceid]=new CPML(deviceid,*cpml);


  d_box=new BOX();
  if(param->usebox==1){
	d_box->init_gpu(deviceid,*box,*param);
  }

  d_r_fld=new FIELD();
  d_r_fld->init_gpu_box(deviceid,*param);

  record_u->cu_init_record(deviceid,*param);
  record_w->cu_init_record(deviceid,*param);

  SLIDE slide(*param);

  TIME time;

  /*Set up source time*/
  double *src;
  int lsrc;
  if(param->sourcetime[0]=='T' || param->sourcetime[0]=='t'){
	trapozoid(param->dt,param->trap1,param->trap2,param->trap3,lsrc,src);
  }else if(param->sourcetime[0]=='G' || param->sourcetime[0]=='g'){
	gaussian(param->dt,param->alpha,lsrc,src);
  }

  double *psrc=new double[lsrc];
  double sum=0;
  for(int i=0;i<lsrc;i++){
	sum+=src[i];
	psrc[i]=sum;
  }






  for(int it=0; it < param->nt; it++){

	slide.current_step=it;

	/*-------------------------------------------------------------------------*/
	/* quick add nobox situation*/
	if(param->usebox==0){
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
	  }


	  cu_step_velocity(deviceid,*param,*d_fld[deviceid], *d_mat[deviceid],*d_cpml[deviceid],slide);
	  if(it%param->itrecord==0){
		//record this gpu portion, and wait for syncronize thread to flush
		record_u->cu_record(deviceid,param->nx,d_fld[deviceid]->U);
		record_w->cu_record(deviceid,param->nx,d_fld[deviceid]->W);
	  }
	  thread->exchange_velocity(deviceid,*param,*d_fld[deviceid],*g_fld);
	  if(it%param->itrecord==0){
		//flush now,since there is also a barrier at exchange stress,so this is safe
		record_u->cu_flush(deviceid);
		record_w->cu_flush(deviceid);
	  }
	}else{
	  /*-------------------------------------------------------------------------*/
	  d_box->cu_copy_source(deviceid,*box);

	  if(param->boxdevice==deviceid){
		d_box->cu_addsource(deviceid,d_fld[deviceid]->U , d_r_fld->U , d_box->p_sfield->U);
		d_box->cu_addsource(deviceid,d_fld[deviceid]->W , d_r_fld->W , d_box->p_sfield->W);
	  }

	  cu_step_stress(deviceid,*param,*d_fld[deviceid] ,*d_mat[deviceid] ,*d_cpml[deviceid],slide);
	  thread->exchange_stress(deviceid,*param,*d_fld[deviceid],*g_fld);

	  if(param->boxdevice==deviceid){
		d_box->cu_subsource(deviceid,d_fld[deviceid]->U , d_r_fld->U , d_box->p_sfield->U);
		d_box->cu_subsource(deviceid,d_fld[deviceid]->W , d_r_fld->W , d_box->p_sfield->W);
	  }

	  cu_step_stress_r(deviceid,*param,*d_r_fld,*d_box->p_bmat);

	  if(param->boxdevice==deviceid){
		d_box->cu_addsource(deviceid,d_fld[deviceid]->Txx , d_r_fld->Txx , d_box->p_sfield->Txx);
		d_box->cu_addsource(deviceid,d_fld[deviceid]->Txz , d_r_fld->Txz , d_box->p_sfield->Txz);
		d_box->cu_addsource(deviceid,d_fld[deviceid]->Tzz , d_r_fld->Tzz , d_box->p_sfield->Tzz);
	  }

	  cu_step_velocity(deviceid,*param,*d_fld[deviceid], *d_mat[deviceid],*d_cpml[deviceid],slide);


	  if(it%param->itrecord==0){
		//record this gpu portion, and wait for syncronize thread to flush
		record_u->cu_record(deviceid,param->nx,d_fld[deviceid]->U);
		record_w->cu_record(deviceid,param->nx,d_fld[deviceid]->W);
	  }
	  thread->exchange_velocity(deviceid,*param,*d_fld[deviceid],*g_fld);
	  if(it%param->itrecord==0){
		//flush now,since there is also a barrier at exchange stress,so this is safe
		record_u->cu_flush(deviceid);
		record_w->cu_flush(deviceid);
	  }

	  if(param->boxdevice==deviceid){
		d_box->cu_subsource(deviceid,d_fld[deviceid]->Txx , d_r_fld->Txx , d_box->p_sfield->Txx);
		d_box->cu_subsource(deviceid,d_fld[deviceid]->Txz , d_r_fld->Txz , d_box->p_sfield->Txz);
		d_box->cu_subsource(deviceid,d_fld[deviceid]->Tzz , d_r_fld->Tzz , d_box->p_sfield->Tzz);
	  }

	  cu_step_velocity_r(deviceid,*param,*d_r_fld,*d_box->p_bmat);
	}



	if(it%param->ntsnap==0 && it!=0){
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

  BOX  box;
  if(param.usebox==1){
	box.init_box(param,mat);
  }

//  mat.mktable();


  FIELD    g_fld;
  g_fld.init_for_full( param );

  CPML    cpml(param);

  FIELD    r_fld;
  r_fld.init_for_box(param);


  FIELD * d_fld[MAXCUDAGPU];
  MATERIAL * d_mat[MAXCUDAGPU];
  CPML * d_cpml[MAXCUDAGPU];

  BOX * d_box;
  FIELD * d_r_fld;

  struct ARGUMENT arg[MAXCUDAGPU];

  THREAD thread(param.ngpu);


  char filename[256];
  RECORD record_u, record_w;

  if(param.use_sta_file){
	sprintf(filename,"%s_U",param.output);
	record_u.init_list(param.U_sta_file,filename);

	sprintf(filename,"%s_W",param.output);
	record_w.init_list(param.W_sta_file,filename);
  }else{
	sprintf(filename,"%s_U",param.output);
	record_u.init_line(param.nrec,param.ixrec0,param.izrec0_u,param.idxrec,param.idzrec,filename);
	sprintf(filename,"%s_W",param.output);
	record_w.init_line(param.nrec,param.ixrec0,param.izrec0_w,param.idxrec,param.idzrec,filename);
  }



  for(int i=0;i<param.ngpu;i++){
	arg[i].deviceid    =  i;
	arg[i].thread      = &thread;
	arg[i].record_u    = &record_u;
	arg[i].record_w    = &record_w;
	arg[i].param       = &param;

	arg[i].cpml        = &cpml;
	arg[i].d_cpml      = d_cpml;

	arg[i].g_fld       = &g_fld;
	arg[i].mat         = &mat;

	arg[i].box         = &box;
	arg[i].d_box       = d_box;

	arg[i].d_fld       = d_fld;
	arg[i].d_mat       = d_mat;

	arg[i].d_r_fld     = d_r_fld;
  }

  if(param.ngpu<=0){
	for(int it=0;it<param.nt;it++){
	  box.setboxsource();
	  box.addsource(g_fld.U , r_fld.U , box.p_sfield->U);
	  box.addsource(g_fld.W , r_fld.W , box.p_sfield->W);

	  step_stress(g_fld , mat , cpml);

	  box.subsource(g_fld.U , r_fld.U , box.p_sfield->U);
	  box.subsource(g_fld.W , r_fld.W , box.p_sfield->W);

	  step_stress_r(r_fld , *box.p_bmat);

	  box.addsource(g_fld.Txx , r_fld.Txx , box.p_sfield->Txx);
	  box.addsource(g_fld.Txz , r_fld.Txz , box.p_sfield->Txz);
	  box.addsource(g_fld.Tzz , r_fld.Tzz , box.p_sfield->Tzz);

	  step_velocity(g_fld, mat,	cpml);

	  box.subsource(g_fld.Txx , r_fld.Txx , box.p_sfield->Txx);
	  box.subsource(g_fld.Txz , r_fld.Txz , box.p_sfield->Txz);
	  box.subsource(g_fld.Tzz , r_fld.Tzz , box.p_sfield->Tzz);

	  step_velocity_r(r_fld , *box.p_bmat);

	  if(it%param.ntsnap==0 && it!=0){
		cout<<it<<endl;
		snapshot(g_fld,it);
	  }
	}
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
  record_u.close();
  record_w.close();
  record_u.toisis(param,0.5, 0.0);
  record_w.toisis(param,0.0,-0.5);
  cout<<"finish\n";
}
