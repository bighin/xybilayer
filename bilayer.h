#ifndef __BILAYER_H__
#define __BILAYER_H__

#include "xy.h"

#define LOWER_LAYER	(0)
#define UPPER_LAYER	(1)

struct bilayer_t
{
	struct spin2d_t *layers[2];
	int lx,ly;
	
	double J[2],K;
};

struct bilayer_t *bilayer_init(int x,int y,double Jup,double Jdown,double K);
void bilayer_fini(struct bilayer_t *b);

double pcc_bilayer(int x,int y,double beta,double Jup,double Jdown,double K);

#endif
