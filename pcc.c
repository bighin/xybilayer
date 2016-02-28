/*
	Probability changing cluster for the 2D XY model

	Reference: PhysRevB.65.184405
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

#define TOTAL_SWEEPS	(20000)
#define THERMALIZATION	(10000)
#define UPDATE_INTERVAL	(10)

	assert(TOTAL_SWEEPS>UPDATE_INTERVAL);

	cfgt=spin2d_init(x,y);
	spin2d_random_configuration(cfgt);

	starting_delta=0.1f;
	target_delta=0.000680f;

	/*
		delta (la variazione di beta dopo ogni passo) viene aggiustato ogni
		UPDATE_INTERVAL passi. La variazione chi viene calcolata in modo tale
		che inizialmente delta valga starting_delta e alla fine valga target_delta;
	*/

	delta=starting_delta;
	chi=pow(target_delta/starting_delta,-1.0f/(TOTAL_SWEEPS/UPDATE_INTERVAL));

	for(c=0;c<TOTAL_SWEEPS;c++)
	{
		short percolating;

		percolating=swendsen_wang_step(cfgt,beta,J);

		beta+=((percolating==1)?(-delta):(+delta));

		if((c>0)&&((c%UPDATE_INTERVAL)==0))
			delta/=chi;
	}

	average=0.0f;
	for(c=0;c<THERMALIZATION;c++)
	{
		short percolating;

		percolating=swendsen_wang_step(cfgt,beta,J);

		beta+=((percolating==1)?(-delta):(+delta));

		average+=beta;
	}

	spin2d_fini(cfgt);
	
	return average/((double)(THERMALIZATION));
}

/*
	Probability changing cluster for the 2D Ising model
*/

double pcc_ising(int x,int y,double beta,double J)
{
	struct ising2d_t *cfgt;
	struct bond2d_t *bonds;
	double starting_delta,target_delta,delta,chi;
	double average;
	int c,d;

#define ISING_TOTAL_SWEEPS	(20000)
#define ISING_THERMALIZATION	(10000)
#define ISING_UPDATE_INTERVAL	(10)

	assert(ISING_TOTAL_SWEEPS>ISING_UPDATE_INTERVAL);

	cfgt=ising2d_init(x,y);

	for(c=0;c<x;c++)
		for(d=0;d<y;d++)
			ising2d_set_spin(cfgt,x,y,-1+2*(gen_random_int()%2));
			
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
	chi=pow(target_delta/starting_delta,-1.0f/(ISING_TOTAL_SWEEPS/ISING_UPDATE_INTERVAL));

	for(c=0;c<ISING_TOTAL_SWEEPS;c++)
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
	for(c=0;c<ISING_THERMALIZATION;c++)
	{
		short percolating;

		percolating=swendsen_wang_ising_step(cfgt,bonds,beta);

		beta+=((percolating==1)?(-delta):(+delta));

		average+=beta;
	}

	ising2d_fini(cfgt);
	
	return average/((double)(ISING_THERMALIZATION));
}