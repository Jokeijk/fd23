#include"../param/param.h"
#include<fstream>
#include<iostream>
#include<iterator>
#include<cstdio>
#include<cstdlib>

bool needsnap(PARAM & param,int it)
{
    if(param.ntsnapfile[0]==0) return false;

	std::ifstream file(param.ntsnapfile);

	if(!file.is_open())
	{
	  fprintf(stderr,"%s Not exist!",param.ntsnapfile);
	  exit(-1);
	}

	std::istream_iterator<int> eos;
	std::istream_iterator<int> iit(file);

	while (iit!=eos) { 
	  int tmp=*iit;
	  iit++;
	  if( it==tmp ) return true;
	}

	file.close();
	return false;
}
