/* Author Dunzhu Li  dli@caltech.edu 
*/
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
#include"source/waveform.h"
#include<iostream>
#include<iomanip>




using namespace std;

struct ARGUMENT{
  int deviceid;
  int truedevice;

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

};

void * singlegpu(void * argument)
{
  struct ARGUMENT * arg=(struct ARGUMENT *)argument;

  int deviceid       = arg->deviceid;
  int truedevice     = arg->truedevice;

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




  thread->selectdevice(truedevice);
  d_fld[deviceid]=new FIELD();
  d_fld[deviceid]->init_gpu_full(deviceid, *param);

  d_mat[deviceid]=new MATERIAL();
  d_mat[deviceid]->init_gpu_full(deviceid, *param, *mat);

  d_cpml[deviceid]=new CPML(deviceid,*cpml);




  record_u->cu_init_record(deviceid,*param);
  record_w->cu_init_record(deviceid,*param);

  SLIDE slide(*param);

  TIME time;

  /*Set up source time*/
  float *src;
  int lsrc;
  if(param->srctime[0]=='T' || param->srctime[0]=='t'){
	trapozoid(param->dt,param->trap1,param->trap2,param->trap3,lsrc,src);
  }else if(param->srctime[0]=='G' || param->srctime[0]=='g'){
	gaussian(param->dt,param->alpha,lsrc,src);
  }

  float *psrc=new float[lsrc];
  float sum=0;
  for(int i=0;i<lsrc;i++){
	sum+=src[i];
	psrc[i]=sum;
  }


  /* restart */
  if(param->restart ==1 ){
	d_fld[deviceid]->cu_load_restart(deviceid,*param,*g_fld);
	d_cpml[deviceid]->cu_load_restart(deviceid,*cpml);
	thread->exchange_velocity(deviceid,*param,*d_fld[deviceid],*g_fld);
	thread->exchange_stress(deviceid,*param,*d_fld[deviceid],*g_fld);
  }
  /* end restart */





  for(int it=0+param->nt0; it < param->nt+param->nt0; it++){

	slide.current_step=it;

	/*-------------------------------------------------------------------------*/
	/* quick add nobox situation*/
	cu_step_stress(deviceid,*param,*d_fld[deviceid] ,*d_mat[deviceid] ,*d_cpml[deviceid],slide);
	thread->exchange_stress(deviceid,*param,*d_fld[deviceid],*g_fld);

	switch (param->srctype)
	{
	  case '1':
		cu_point_source(deviceid,it,lsrc,src,*param,*d_fld[deviceid]);
		break;
	  case 'p':
	  case 'P':
		cu_point_source_p(deviceid,it,lsrc,psrc,*param,*d_fld[deviceid]);
		break;
	  default:
		fprintf(stderr,"Wrong srctype! Specify 1 or p.");
		exit(-1);
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
		<<setw(10)<<time.elips()/(it-param->nt0+1)/1000.0<<" msec per step\n";
	  cout.flush();
	}
  }

  /* restart */
  if(param->savestate ==1 ){
	thread->exchange_velocity(deviceid,*param,*d_fld[deviceid],*g_fld);
	thread->exchange_stress(deviceid,*param,*d_fld[deviceid],*g_fld);
	d_fld[deviceid]->cu_save_state(deviceid,*param,*g_fld);
	d_cpml[deviceid]->cu_save_state(deviceid,*cpml);
  }
  /* end restart */

  return 0;

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

  /* g_cpml[0] save one cpml, which corresonding to d_cpml[0] */
  CPML * g_cpml[MAXCUDAGPU];
  for (int i=0;i<MAXCUDAGPU;i++){
	g_cpml[i] = new CPML(param);
  }
  


  FIELD * d_fld[MAXCUDAGPU];
  MATERIAL * d_mat[MAXCUDAGPU];
  CPML * d_cpml[MAXCUDAGPU];

  struct ARGUMENT arg[MAXCUDAGPU];

  THREAD thread(param.ngpu);

  char filename[STRLEN];
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
	arg[i].deviceid    = i;
    arg[i].truedevice  = param.gpuid[i];
	arg[i].thread      = &thread;
	arg[i].record_u    = &record_u;
	arg[i].record_w    = &record_w;
	arg[i].param       = &param;

	arg[i].cpml        = g_cpml[i];
	arg[i].d_cpml      = d_cpml;

	arg[i].g_fld       = &g_fld;
	arg[i].mat         = &mat;

	arg[i].d_fld       = d_fld;
	arg[i].d_mat       = d_mat;

  }


  if(param.restart ==1) {

	 ifstream statefid;
	 char statefile[256];
	 sprintf(statefile,"state%06d",param.nt0-1);
	 statefid.open(statefile,ios::binary);
	 size_t tsize=sizeof(float)*param.nx*param.nz;
	 statefid.read((char*)g_fld.U,  tsize);
	 statefid.read((char*)g_fld.W,  tsize);
	 statefid.read((char*)g_fld.Txx,tsize);
	 statefid.read((char*)g_fld.Txz,tsize);
	 statefid.read((char*)g_fld.Tzz,tsize);
	 for(int i=0;i<param.ngpu;i++){
	   tsize=sizeof(float)*param.npml*2*param.nz;
	   statefid.read((char*)g_cpml[i]->psi.Txx_x,tsize);
	   statefid.read((char*)g_cpml[i]->psi.Txz_x,tsize);
	   statefid.read((char*)g_cpml[i]->psi.U_x,  tsize);
	   statefid.read((char*)g_cpml[i]->psi.W_x,  tsize);
	   tsize=sizeof(float)*param.npml*2*param.nx;                    
	   statefid.read((char*)g_cpml[i]->psi.Tzz_z,tsize);
	   statefid.read((char*)g_cpml[i]->psi.Txz_z,tsize);
	   statefid.read((char*)g_cpml[i]->psi.U_z,  tsize);
	   statefid.read((char*)g_cpml[i]->psi.W_z,  tsize);
	 }
	 statefid.close();
  }

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

  if(param.savestate ==1) {

	 ofstream statefid;
	 char statefile[256];
	 sprintf(statefile,"state%06d",param.nt+param.nt0-1);
	 statefid.open(statefile,ios::binary);
	 size_t tsize=sizeof(float)*param.nx*param.nz;
	 statefid.write((char*)g_fld.U,  tsize);
	 statefid.write((char*)g_fld.W,  tsize);
	 statefid.write((char*)g_fld.Txx,tsize);
	 statefid.write((char*)g_fld.Txz,tsize);
	 statefid.write((char*)g_fld.Tzz,tsize);
	 for(int i=0;i<param.ngpu;i++){
	   tsize=sizeof(float)*param.npml*2*param.nz;
	   statefid.write((char*)g_cpml[i]->psi.Txx_x,tsize);
	   statefid.write((char*)g_cpml[i]->psi.Txz_x,tsize);
	   statefid.write((char*)g_cpml[i]->psi.U_x,  tsize);
	   statefid.write((char*)g_cpml[i]->psi.W_x,  tsize);
	   tsize=sizeof(float)*param.npml*2*param.nx;                     
	   statefid.write((char*)g_cpml[i]->psi.Tzz_z,tsize);
	   statefid.write((char*)g_cpml[i]->psi.Txz_z,tsize);
	   statefid.write((char*)g_cpml[i]->psi.U_z,  tsize);
	   statefid.write((char*)g_cpml[i]->psi.W_z,  tsize);
	 }
	 statefid.close();
  }

  record_u.close();
  record_w.close();

  /* source is the the Mxx,Mzz position, 
   * but u and w is recorded at the staggered position*/
  record_u.toisis(param,0.5, 0.0);
  record_w.toisis(param,0.0,-0.5);

  cout<<"finish\n";
}
