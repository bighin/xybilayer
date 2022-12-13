/*
	Probability changing cluster for the two-dimensional XY model.

	For an introduction to the method, please refer to: Y. Tomita and Y. Okabe, Phys. Rev. B 65, 184405 (2002).
*/

#include <assert.h>

#include "utils.h"
#include "xy.h"

double pcc(int x,int y,double beta,double J)
{
	struct spin2d_t *cfgt;
	double starting_delta,target_delta,delta,chi;
	double average;
	int c;

#define THERMALIZATION	(20000)
#define POST_THERMALIZATION	(10000)
#define UPDATE_INTERVAL	(10)

	assert(THERMALIZATION>UPDATE_INTERVAL);

	cfgt=spin2d_init(x,y);
	spin2d_random_configuration(cfgt);

	starting_delta=0.1f;
	target_delta=0.000680f;

	/*
		\delta (the \beta variation after each step) is adjusted each UPDATE_INTERVAL steps.

		The adjustment is performed such that in the beginning \delta is starting_delta,
		whereas at the end is target_delta.
	*/

	delta=starting_delta;
	chi=pow(target_delta/starting_delta,-1.0f/(THERMALIZATION/UPDATE_INTERVAL));

	for(c=0;c<THERMALIZATION;c++)
	{
		short percolating;

		percolating=swendsen_wang_step(cfgt,beta,J);

		beta+=((percolating==1)?(-delta):(+delta));

		if((c>0)&&((c%UPDATE_INTERVAL)==0))
			delta/=chi;
	}

	average=0.0f;
	for(c=0;c<POST_THERMALIZATION;c++)
	{
		short percolating;

		percolating=swendsen_wang_step(cfgt,beta,J);

		beta+=((percolating==1)?(-delta):(+delta));

		average+=beta;
	}

	spin2d_fini(cfgt);
	
	return average/((double)(POST_THERMALIZATION));
}

/*
	Probability changing cluster for the two-dimensional Ising model
*/

double pcc_ising(int x,int y,double beta,double J)
{
	struct ising2d_t *cfgt;
	struct bond2d_t *bonds;
	double starting_delta,target_delta,delta,chi;
	double average;
	int c,d;

#define ISING_THERMALIZATION		(20000)
#define ISING_POST_THERMALIZATION	(10000)
#define ISING_UPDATE_INTERVAL		(10)

	assert(ISING_THERMALIZATION>ISING_UPDATE_INTERVAL);

	cfgt=ising2d_init(x,y);

	for(c=0;c<x;c++)
		for(d=0;d<y;d++)
			ising2d_set_spin(cfgt,c,d,-1+2*(gen_random_int()%2));
			
	bonds=bond2d_init(x,y);
	
	for(c=0;c<x;c++)
	{
		for(d=0;d<y;d++)
		{
			bond2d_set_value(bonds,c,d,DIR_X,J);
			bond2d_set_value(bonds,c,d,DIR_Y,J);
		}
	}

	starting_delta=0.1f;
	target_delta=0.000680f;

	delta=starting_delta;
	chi=pow(target_delta/starting_delta,-1.0f/(ISING_THERMALIZATION/ISING_UPDATE_INTERVAL));

	for(c=0;c<ISING_THERMALIZATION;c++)
	{
		short percolating;

		percolating=swendsen_wang_ising_step(cfgt,bonds,beta);

		beta+=((percolating==1)?(-delta):(+delta));

		if(beta<=delta)
			beta=delta;

		if((c>0)&&((c%ISING_UPDATE_INTERVAL)==0))
		{
			delta/=chi;
		}
	}

	average=0.0f;
	for(c=0;c<ISING_POST_THERMALIZATION;c++)
	{
		short percolating;

		percolating=swendsen_wang_ising_step(cfgt,bonds,beta);

		beta+=((percolating==1)?(-delta):(+delta));

		average+=beta;
	}

	ising2d_fini(cfgt);
	
	return average/((double)(ISING_POST_THERMALIZATION));
}
