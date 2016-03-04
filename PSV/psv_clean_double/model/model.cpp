/* Make a raster grid from a polygon description of a model
 *  *  R. Clayton, Aug, 2007
 *   */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<algorithm>
#include"model.h"

/* definitions for model specifications */
struct vert
   {
   	double x, z;
   };
#define VERT struct vert

struct poly
   {
   	int	n;
	double	vp;
	double	vs;
	double	den;
	VERT	*vert;
   };

#define NDEFS	1000
struct md_def
   {
   	char name[16];
	double	xp;
	double	zp;
   };

int getargs(char *line, char **args, int nargmax, char *buf);
int inside(double  xp,double  yp,VERT *v,int n);
int my_getline(char *line, int max, FILE *fd);
void unmy_getline(char *line, int n);
int lookup(char *name, struct md_def *defs,int ndefs, double *xp, double *zp);
void uppercase(char *str);

double
amax1(double a1, double a2)
{
  return (a1 > a2 ? a1 : a2);
}

  double
amin1(double a1, double a2)
{
  return (a1 < a2 ? a1 : a2);
}

int getargs(char *line, char **args, int nargmax, char *buf)
{
  char *p, *pp;
  int narg;

  p = line;
  if(buf != NULL) /* make a copy of line */
  {
	for(p=buf; *line != '\0'; ) *p++ = *line++;
	p= buf;
  }
  /* map- newline into null */
  for(pp=p; *pp != '\0'; pp++)
	if(*pp == '\n')
	{
	  *pp= '\0';
	  break;
	}
  narg= 0;
  while( *p != '\0' &&  narg < nargmax)
  {
	while(*p == ' ' || *p == '\t') p++;
	if(*p == '\0') break;
	if(*p == '"')	/* quoted token */
	{
	  p++;
	  if(*p == '\0') break;
	  args[narg]= p;
	  while(*p != '"' && *p != '\0') p++;
	}
	else
	{
	  args[narg]= p;
	  while(*p != ' ' && *p != '\t' && *p != '\0') p++;
	}
	narg++;
	if(*p == '\0') break;
	*p++ = '\0';
  }
  return(narg);
}

struct poly md_p[4000];
int md_np       =0;

struct md_def  defs[NDEFS];
int ndef        =0;

VERT vert[10000];

#define  VP(ix,iz)  vp[ix + nx*(iz)]
#define  VS(ix,iz)  vs[ix + nx*(iz)]
#define DEN(ix,iz) d[ix + nx*(iz)]

struct poly *
get_model_poly()
{
  return(md_p);
}

int get_model_npoly()
{
  return(md_np);
}

struct md_def *
get_model_defs()
{
  return(defs);
}

int get_model_ndefs()
{
  return(ndef);
}

void get_acou_model( double *vp,
	double *d,
	int nx,
	int nz,
	double h,
	double x0,
	double z0)
{
  double xp, zp, vel, den, fraction;
  int ix, iz, i;

  if(md_np <= 0)
  {
	fprintf(stderr,"must load model first\n");
	exit(-1);
  }
  /* paint the background */
  vel= md_p[0].vp;
  den= md_p[0].den;
  for (ix=0 ; ix<nx ; ix++)
	for (iz=0 ; iz<nz ; iz++)
	{
	  VP(ix,iz) = vel;
	  DEN(ix,iz)= den;
	}

  for (ix=0 ; ix<nx ; ix++)
  {
	for (iz=0 ; iz<nz ; iz++)
	{
	  xp = x0+h*((double) ix);
	  zp = z0+h*((double) iz);
	  for (i=0 ; i<md_np ; i++)
	  {
		if(inside(xp,zp,md_p[i].vert,md_p[i].n))
		{
		  vel = md_p[i].vp;
		  den = md_p[i].den;
		  VP(ix,iz)= vel;
		  DEN(ix,iz)= den;
		}
	  }
	}
  }
}
void get_elas_model( double *vp,
	double *vs,
	double *d,
	int nx,
	int nz,
	double h,
	double x0,
	double z0)
{
  double xp, zp, velp, vels, den, fraction;
  int ix, iz, i, nuncover;

  if(md_np <= 0)
  {
	fprintf(stderr,"must load model first\n");
	exit(-1);
  }
  /* paint the background */
  velp= md_p[0].vp;
  vels= md_p[0].vs;
  den= md_p[0].den;
  for (ix=0 ; ix<nx ; ix++)
	for (iz=0 ; iz<nz ; iz++)
	{
	  VP(ix,iz) = velp;
	  VS(ix,iz) = vels;
	  DEN(ix,iz)= den;
	}


  //When polygon are intercept,assume the later one is true polygon
  for(int i=1;i<md_np;i++){
	double x1=1E28;
	double x2=-1E28;
	double z1=1E28;
	double z2=-1E28;
	for(int j=0;j<md_p[i].n;j++){
	  x1=std::min(x1,md_p[i].vert[j].x);
	  x2=std::max(x2,md_p[i].vert[j].x);
	  z1=std::min(z1,md_p[i].vert[j].z);
	  z2=std::max(z2,md_p[i].vert[j].z);
	}
	int ix1=x1/h-1; ix1=ix1<0?0:ix1; ix1=ix1>nx-1?nx-1:ix1;
	int ix2=x2/h+1; ix2=ix2<0?0:ix2; ix2=ix2>nx-1?nx-1:ix2;
	int iz1=z1/h-1; iz1=iz1<0?0:iz1; iz1=iz1>nz-1?nz-1:iz1;
	int iz2=z2/h+1; iz2=iz2<0?0:iz2; iz2=iz2>nz-1?nz-1:iz2;
	for(int iz=iz1;iz<=iz2;iz++){
	  for(int ix=ix1;ix<=ix2;ix++){
		xp = x0+h*((double) ix);
		zp = z0+h*((double) iz);
		if(inside(xp,zp,md_p[i].vert,md_p[i].n))
		{
		  velp = md_p[i].vp;
		  vels = md_p[i].vs;
		  den = md_p[i].den;
		  VP(ix,iz)= velp;
		  VS(ix,iz)= vels;
		  DEN(ix,iz)= den;
		}
	  }
	}
  }
}

#define MAXARGS 100
int load_model(char *file)
{
  FILE *fd;
  int i, n, linenum, k, nvert;
  VERT *pvert;
  double xp, zp, vp, vs, den;
  char line[128], key[16], sym[16], name[16];
  char *args[MAXARGS];
  int nargs;

  if( (fd= fopen(file,"r")) == NULL )
  {
	fprintf(stderr,"cannot open file = %s\n",file);
	exit(-1);
  }

  linenum= 0;
  ndef= 0;
  pvert= vert;
  while( (n= my_getline(line,128,fd)) >= 0)
  {

	linenum++;
	if(n == 0) continue;
	if(line[0] == '#') continue;
	if( (nargs= getargs(line,args,MAXARGS,NULL)) < 1)
	{
	  fprintf(stderr,"syntax error at line %d in file %s\n",
		  linenum, file);
	  fprintf(stderr,"no keyword\n");
	  exit(-1);
	}
	uppercase(args[0]);
	if(strcmp("DEF",args[0]) == 0)
	{
	  if(nargs !=4)
	  {
		fprintf(stderr,"syntax error at line %d in file %s\n",linenum,file);
		fprintf(stderr,"insufficient args in def command\n");
		exit(-1);
	  }
	  strcpy(defs[ndef].name,args[1]);
	  defs[ndef].xp= atof(args[2]);
	  defs[ndef].zp= atof(args[3]);
	  ndef++;
	  continue;
	}
	if(strcmp("BACK",args[0]) == 0)
	{
	  if(nargs !=4)
	  {
		fprintf(stderr,"syntax error at line %d in file %s\n",linenum,file);
		fprintf(stderr,"insufficient args in back command\n");
		exit(-1);
	  }
	  vp = atof(args[1]);
	  vs = atof(args[2]);
	  den= atof(args[3]);
	  if(vs < 0.0) vs= -vp * vs;
	  if(den < 0.0) den= -vp * den;
	  md_p[0].vp = vp;
	  md_p[0].vs = vs;
	  md_p[0].den = den;
	  md_p[0].vert= pvert;
	  /*
	   *                         md_p[0].vert[0].x= xmin; md_p[0].vert[0].z= zmin;
	   *                                                 md_p[0].vert[1].x= xmax; md_p[0].vert[1].z= zmin;
	   *                                                                         md_p[0].vert[2].x= xmax; md_p[0].vert[2].z= zmax;
	   *                                                                                                 md_p[0].vert[3].x= xmin; md_p[0].vert[3].z= zmax;
	   *                                                                                                                         */
	  md_p[0].n= 4;
	  md_np=1;
	  pvert += 4;
	  continue;
	}

	if(strcmp("POLY",args[0]) == 0)
	{
	  double z1, z2;
	  if(nargs !=4)
	  {
		fprintf(stderr,"syntax error at line %d in file %s\n",linenum,file);
		fprintf(stderr,"insufficient args in HLAYER command\n");
		exit(-1);
	  }
	  vp = atof(args[1]);
	  vs = atof(args[2]);
	  den= atof(args[3]);
	  if(vs < 0.0) vs= -vp * vs;
	  if(den < 0.0) den= -vp * den;
	  md_p[md_np].vp= vp;
	  md_p[md_np].vs= vs;
	  md_p[md_np].den = den;
	  md_p[md_np].vert= pvert;
	  nvert= 0;
	  while( (n= my_getline(line,128,fd)) >= 0)
	  {
		if(line[0]!= '\t')
		{
		  unmy_getline(line,n);
		  break;
		}
		nargs= getargs(line,args,MAXARGS,NULL);
		for(i=0; i<nargs; i++)
		{
		  if(isalpha(args[0][0]))
			lookup(args[i],defs,ndef,&xp,&zp);
		  else
		  {
			xp= atof(args[i]);
			zp= atof(args[i+1]);
			i++;
		  }

		  pvert[nvert].x= xp;
		  pvert[nvert].z= zp;
		  nvert++;
		}
	  }
	  md_p[md_np].n= nvert;
	  pvert += nvert;
	  md_np++;
	  continue;
	}
	fprintf(stderr,"bad model command at line %d in file %s\n",
		linenum,file);
	exit(-1);
  }
  fclose(fd);
}

int lookup(char *name, struct md_def *defs,int ndefs, double *xp, double *zp)
{
  int i;
  for(i=0; i<ndefs; i++)
  {
	if(strcmp(name,defs[i].name) == 0)
	{
	  *xp= defs[i].xp;
	  *zp= defs[i].zp;
	  return 0;
	}
  }
  fprintf(stderr,"undefined symbol= %s\n",name);
  exit(-1);
}


void uppercase(char *str)
{
  int i;
  while(*str != '\0')
  {
	if(*str >= 'a' && *str <= 'z') *str += 'A' - 'a';
	str++;
  }
}

int do_unget    =0;
int n_unget     =0;
void unmy_getline(char *line, int n)
{
  do_unget= 1;
  n_unget= n;
}
int my_getline(char *line, int max, FILE *fd)
{
  char c;
  int n;

  if(do_unget)
  {
	n= n_unget;
	n_unget= 0;
	do_unget= 0;
	return(n);
  }
  n=0;
  while( (c=getc(fd)) != EOF )
  {
	if(c == '\n')
	{
	  *line= '\0';
	  return(n);
	}
	*line++ = c;
	n++;
	if(n >= max)
	{
	  fprintf(stderr,"line too long\n");
	  exit(-1);
	}
  }

  return(EOF);
}

int inside(double xp, double yp, VERT *v, int n)
{
  double x1, x2, y1, y2, xi;
  int nleft, j1;

  /*
   * fprintf(stdout,"vert n= %d xp=%4.1f yp=%4.1f ",n,xp,yp);
   * for(j1=0; j1<n; j1++) fprintf(stdout,"(%4.1f,%4.1f) ",v[j1].x,v[j1].z);
   * fprintf(stdout,"\n");
   * */
  nleft = 0;
  x2 = v[n - 1].x;
  y2 = v[n - 1].z;
  for (j1=0; j1<n; j1++)
  {
	x1 = x2;
	y1 = y2;
	x2 = v[j1].x;
	y2 = v[j1].z;
	if (amin1(y1, y2) >= yp) continue;
	if (amax1(y1, y2) < yp) continue;
	if (y1 == y2) goto hop27;
	xi = x1 + (yp - y1) * (x2 - x1) / (y2 - y1);
	if (xi == xp) return 1;
	if (xi > xp) nleft = nleft + 1;
	continue;
hop27:
	if (xp > amax1(x1, x2)) continue;
	if (xp >= amin1(x1, x2)) return 1;
	nleft = nleft + 1;
  }
  return (nleft % 2);
}

