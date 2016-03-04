/* Author Dunzhu Li  dli@caltech.edu 
*/
#ifndef _TIMER_H_
#define _TIMER_H_

#include<sys/time.h>
class TIME{
  public:
	struct timeval starttime;
	struct timeval nowtime;
	TIME(){
	  gettimeofday(&starttime,0);
	}
	double elips(){
	  gettimeofday(&nowtime,0);
	  double time0=starttime.tv_sec*1000000+starttime.tv_usec;
	  double time1=nowtime.tv_sec*1000000+nowtime.tv_usec;
	  return time1-time0;
	}

};

#endif
