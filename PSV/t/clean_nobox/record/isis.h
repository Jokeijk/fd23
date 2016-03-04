/* Author Dunzhu Li  dli@caltech.edu 
*/
#ifndef _ISIS_H_
#define _ISIS_H_

#define NUSER	4
#define ushort unsigned short
#include"../param/param.h"
struct traceinfo	/* the on disk or tape trace header */
   {
	float	xs,ys,zs;	/* source location */
	float	xr,yr,zr;	/* receiver location */
	float	t0;		/* time at start of data */
	int	nt;		/* number of samples */
	float	samplerate;	/* sample rate i.e. 4 mils = 0.004 */
	int	srcflagnum;	/* source survey flag number */
	int	recflagnum;	/* receiver survey flag number */
	int	status;		/* status bits for trace */
	float	gain;		/* gain of recorder */
	float   dtstatic;       /* trace time static */
	float   ampstatic;      /* trace amplitude static */
	ushort	recorder;	/* recorder # */
	ushort	chan;		/* recording channel */
	ushort	scomp;		/* source comp. */
	ushort	rcomp;		/* receiver comp. */
	int	user[NUSER];	/* user space */
   };
#define TRSIZE	sizeof(struct traceinfo)

void output_isis(char * output,float * rec,int nrec,int ntrecord,int ixrec0,int izrec0,float h,float dt,int xs,int zs);
#endif
