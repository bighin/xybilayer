/*
	Two-dimensional XY bilayer
*/

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "bilayer.h"
#include "xy.h"
#include "common.h"

struct bilayer_t *bilayer_init(int x,int y,double Jup,double Jdown,double K)
{
	struct bilayer_t *ret;

	if(!(ret=malloc(sizeof(struct bilayer_t))))
		return NULL;

	ret->layers[UPPER_LAYER]=spin2d_init(x,y);
	ret->layers[LOWER_LAYER]=spin2d_init(x,y);

	if((!(ret->layers[LOWER_LAYER]))||(!(ret->layers[UPPER_LAYER])))
	{
		if(ret->layers[LOWER_LAYER])
			spin2d_fini(ret->layers[LOWER_LAYER]);

		if(ret->layers[UPPER_LAYER])
			spin2d_fini(ret->layers[UPPER_LAYER]);

		return NULL;
	}

	ret->lx=x;
	ret->ly=y;

	ret->J[LOWER_LAYER]=Jdown;
	ret->J[UPPER_LAYER]=Jup;
	ret->K=K;

	return ret;
}

void bilayer_fini(struct bilayer_t *b)
{
	if(b)
	{
		if(b->layers[LOWER_LAYER])
			spin2d_fini(b->layers[LOWER_LAYER]);

		if(b->layers[UPPER_LAYER])
			spin2d_fini(b->layers[UPPER_LAYER]);

		free(b);
	}
}

void wolff_embedded_cluster(struct spin2d_t *cfgt,struct vec2d_t **parallel,struct vec2d_t **orthogonal,
                            double alpha,struct ising2d_t **epsilon1,struct ising2d_t **epsilon2,
			    struct bond2d_t **Jij1,struct bond2d_t **Jij2,double J)
{
	int x,y;

	/*
		Wolff's embedded cluster

		Separo i gradi di libertà di spin paralleli e ortogonali al vettore unitario
		e costruisco un modello di Ising a partire dai gradi di spin paralleli a r.
	*/

	*epsilon1=ising2d_init(cfgt->lx,cfgt->ly);
	*epsilon2=ising2d_init(cfgt->lx,cfgt->ly);

	*parallel=vec2d_init(cfgt->lx,cfgt->ly);
	*orthogonal=vec2d_init(cfgt->lx,cfgt->ly);

	spin2d_projection(cfgt,alpha,*parallel,*orthogonal,*epsilon1,*epsilon2);

	/*
		Dapprima costruisco le variabil Jij, definite su ciascun bond
	*/

	*Jij1=bond2d_init(cfgt->lx,cfgt->ly);
	*Jij2=bond2d_init(cfgt->lx,cfgt->ly);

	for(x=0;x<cfgt->lx;x++)
	{
		for(y=0;y<cfgt->ly;y++)
		{
			double si[2],sj[2],r[2],q[2],z1,z2;

			r[0]=cos(alpha);
			r[1]=sin(alpha);

			q[0]=cos(alpha+M_PI/2.0f);
			q[1]=sin(alpha+M_PI/2.0f);

			si[0]=cos(spin2d_get_spin(cfgt,x,y));
			si[1]=sin(spin2d_get_spin(cfgt,x,y));

			if((x+1)<cfgt->lx)
			{
				sj[0]=cos(spin2d_get_spin(cfgt,x+1,y));
				sj[1]=sin(spin2d_get_spin(cfgt,x+1,y));

				z1=J*fabs(si[0]*r[0]+si[1]*r[1])*fabs(sj[0]*r[0]+sj[1]*r[1]);
				z2=J*fabs(si[0]*q[0]+si[1]*q[1])*fabs(sj[0]*q[0]+sj[1]*q[1]);

				bond2d_set_value(*Jij1,x,y,DIR_X,z1);
				bond2d_set_value(*Jij2,x,y,DIR_X,z2);
			}

			if((y+1)<cfgt->ly)
			{
				sj[0]=cos(spin2d_get_spin(cfgt,x,y+1));
				sj[1]=sin(spin2d_get_spin(cfgt,x,y+1));

				z1=J*fabs(si[0]*r[0]+si[1]*r[1])*fabs(sj[0]*r[0]+sj[1]*r[1]);
				z2=J*fabs(si[0]*q[0]+si[1]*q[1])*fabs(sj[0]*q[0]+sj[1]*q[1]);

				bond2d_set_value(*Jij1,x,y,DIR_Y,z1);
				bond2d_set_value(*Jij2,x,y,DIR_Y,z2);
			}
		}
	}
}

short swendsen_wang_ising_bilayer_step(struct ising2d_t *epsilon[2],struct bond2d_t *Jij[2],
                                       struct vbond2d_t *vbonds,double beta)
{
	struct ibond2d_t *bonds[2];
	struct ivbond2d_t *ivbonds;
	struct bclusters_t *bclusters;

	int c,x,y;
	int nr_clusters;
	short percolating;

	/*
		A partire da Jij ora costruisco in modo stocastico i bond tra i vari siti.
	*/

	for(c=0;c<=1;c++)
	{
		bonds[c]=ibond2d_init(epsilon[c]->lx,epsilon[c]->ly);

		for(x=0;x<epsilon[c]->lx;x++)
		{
			for(y=0;y<epsilon[c]->ly;y++)
			{
				if((x+1)<epsilon[c]->lx)
				{
					ibond2d_set_value(bonds[c],x,y,DIR_X,0);

					/*
						Se due spin adiacenti lungo x sono concordi allora
						attivo il bond con probabilità p.
					*/

					if(ising2d_get_spin(epsilon[c],x,y)==ising2d_get_spin(epsilon[c],x+1,y))
					{
						double p=1.0f-exp(-2.0f*beta*bond2d_get_value(Jij[c],x,y,DIR_X));

						if(gen_random_number()<p)
							ibond2d_set_value(bonds[c],x,y,DIR_X,1);
					}
				}

				if((y+1)<epsilon[c]->ly)
				{
					/*
						Stessa cosa lungo y.
					*/

					ibond2d_set_value(bonds[c],x,y,DIR_Y,0);

					if(ising2d_get_spin(epsilon[c],x,y)==ising2d_get_spin(epsilon[c],x,y+1))
					{
						double p=1.0f-exp(-2.0f*beta*bond2d_get_value(Jij[c],x,y,DIR_Y));

						if(gen_random_number()<p)
							ibond2d_set_value(bonds[c],x,y,DIR_Y,1);
					}
				}
			}
		}
	}

	/*
		Creo anche i bonds verticali
	*/

	ivbonds=ivbond2d_init(epsilon[LOWER_LAYER]->lx,epsilon[LOWER_LAYER]->ly);

	for(x=0;x<epsilon[LOWER_LAYER]->lx;x++)
	{
		for(y=0;y<epsilon[LOWER_LAYER]->ly;y++)
		{
			ivbond2d_set_value(ivbonds,x,y,0);

			if(ising2d_get_spin(epsilon[LOWER_LAYER],x,y)==ising2d_get_spin(epsilon[UPPER_LAYER],x,y))
			{
				double p=1.0f-exp(-2.0f*beta*vbond2d_get_value(vbonds,x,y));

				if(gen_random_number()<p)
					ivbond2d_set_value(ivbonds,x,y,1);
			}
		}
	}

	/*
		Identifico i clusters
	*/


	bclusters=bclusters_init(epsilon[LOWER_LAYER]->lx,epsilon[LOWER_LAYER]->ly);
	nr_clusters=ising2d_identify_bclusters(bonds,ivbonds,bclusters);

	/*
		Valuto le dimensioni di ciascun cluster e stabilisco se c'è stata percolazione
	*/

	percolating=0;
	for(c=1;c<nr_clusters;c++)
	{
		int xdim,ydim;
		short upper,lower,this_percolating=0;

		bcluster_dimensions(epsilon,bclusters,c,&xdim,&ydim);

#define TOPOLOGICAL
#ifdef TOPOLOGICAL
		if(xdim==epsilon[0]->lx)
		{
			for(y=0;y<epsilon[0]->ly;y++)
			{
				int a,b;

				a=bclusters_get_value(bclusters,0,y,LOWER_LAYER);
				b=bclusters_get_value(bclusters,epsilon[0]->lx-1,y,LOWER_LAYER);

				if((a==b)&&(a==c))
					this_percolating=1;
			}

			for(y=0;y<epsilon[0]->ly;y++)
			{
				int a,b;

				a=bclusters_get_value(bclusters,0,y,UPPER_LAYER);
				b=bclusters_get_value(bclusters,epsilon[0]->lx-1,y,UPPER_LAYER);

				if((a==b)&&(a==c))
					this_percolating=1;
			}
		}

		if(ydim==epsilon[0]->ly)
		{
			for(x=0;x<epsilon[0]->lx;x++)
			{
				int a,b;

				a=bclusters_get_value(bclusters,x,0,LOWER_LAYER);
				b=bclusters_get_value(bclusters,x,epsilon[0]->ly-1,LOWER_LAYER);

				if((a==b)&&(a==c))
					this_percolating=1;
			}

			for(x=0;x<epsilon[0]->lx;x++)
			{
				int a,b;

				a=bclusters_get_value(bclusters,x,0,UPPER_LAYER);
				b=bclusters_get_value(bclusters,x,epsilon[0]->ly-1,UPPER_LAYER);

				if((a==b)&&(a==c))
					this_percolating=1;
			}
		}
#else
		if((xdim==epsilon[0]->lx)||(ydim==epsilon[0]->ly))
				this_percolating=1;
#endif

		upper=lower=0;
		for(x=0;x<epsilon[0]->lx;x++)
		{
			for(y=0;y<epsilon[0]->ly;y++)
			{
				if(bclusters_get_value(bclusters,x,y,UPPER_LAYER)==c)
					upper=1;

				if(bclusters_get_value(bclusters,x,y,LOWER_LAYER)==c)
					lower=1;
			}
		}

		/*
			FIXME: i criteri di percolazione forse vanno rivisti
		*/

		//if((upper==1)&&(this_percolating==1))
		//	percolating=1;

		if(this_percolating==1)
			percolating=1;
	}

	/*
		Flippo i clusters appena creati
	*/

	for(c=1;c<nr_clusters;c++)
	{
		if(gen_random_number()<0.5f)
			reset_bcluster(epsilon,bclusters,c,1);
		else
			reset_bcluster(epsilon,bclusters,c,-1);
	}

	ibond2d_fini(bonds[LOWER_LAYER]);
	ibond2d_fini(bonds[UPPER_LAYER]);

	bclusters_fini(bclusters);
	ivbond2d_fini(ivbonds);

	return percolating;
}

short swendsen_wang_step_bilayer(struct bilayer_t *cfgt,double beta)
{
	struct vec2d_t *parallel[2],*orthogonal[2];
	struct ising2d_t *epsilon1[2],*epsilon2[2];
	struct bond2d_t *Jij1[2],*Jij2[2];
	struct vbond2d_t *vbonds1,*vbonds2;

	double alpha;
	int l,x,y;
	short percolating;

	/*
		Scelgo un vettore unitario casuale r, definito dall'angolo alpha
	*/

	alpha=gen_random_number()*2.0f*M_PI;

	/*
		Applico l'algoritmo embedded cluster su ciascun layer
	*/

	for(l=0;l<=1;l++)
	{
		wolff_embedded_cluster(cfgt->layers[l],&parallel[l],&orthogonal[l],alpha,
	                               &epsilon1[l],&epsilon2[l],&Jij1[l],&Jij2[l],cfgt->J[l]);
	}

	/*
		Mancano da costruire i bond verticali
	*/

	vbonds1=vbond2d_init(cfgt->lx,cfgt->ly);
	vbonds2=vbond2d_init(cfgt->lx,cfgt->ly);

	for(x=0;x<cfgt->lx;x++)
	{
		for(y=0;y<cfgt->ly;y++)
		{
			double si[2],sj[2],r[2],q[2],z1,z2;

			r[0]=cos(alpha);
			r[1]=sin(alpha);

			q[0]=cos(alpha+M_PI/2.0f);
			q[1]=sin(alpha+M_PI/2.0f);

			si[0]=cos(spin2d_get_spin(cfgt->layers[LOWER_LAYER],x,y));
			si[1]=sin(spin2d_get_spin(cfgt->layers[LOWER_LAYER],x,y));

			sj[0]=cos(spin2d_get_spin(cfgt->layers[UPPER_LAYER],x,y));
			sj[1]=sin(spin2d_get_spin(cfgt->layers[UPPER_LAYER],x,y));

			z1=cfgt->K*fabs(si[0]*r[0]+si[1]*r[1])*fabs(sj[0]*r[0]+sj[1]*r[1]);
			z2=cfgt->K*fabs(si[0]*q[0]+si[1]*q[1])*fabs(sj[0]*q[0]+sj[1]*q[1]);

			vbond2d_set_value(vbonds1,x,y,z1);
			vbond2d_set_value(vbonds2,x,y,z2);
		}
	}

	/*
		Avendo creato due modelli di Ising effettivi uso l'algoritmo di Swendsen-Wang su quelli.
	*/

	percolating=swendsen_wang_ising_bilayer_step(epsilon1,Jij1,vbonds1,beta);
	percolating+=swendsen_wang_ising_bilayer_step(epsilon2,Jij2,vbonds2,beta);

	/*
		Ricostruisco gli spin iniziali.
	*/

	for(l=0;l<=1;l++)
		spin2d_reconstruct(cfgt->layers[l],alpha,parallel[l],orthogonal[l],epsilon1[l],epsilon2[l]);

	/*
		Cleaning up and returning the percolation status.
	*/

	for(l=0;l<=1;l++)
	{
		vec2d_fini(parallel[l]);
		vec2d_fini(orthogonal[l]);

		ising2d_fini(epsilon1[l]);
		ising2d_fini(epsilon2[l]);
		bond2d_fini(Jij1[l]);
		bond2d_fini(Jij2[l]);
	}

	vbond2d_fini(vbonds1);
	vbond2d_fini(vbonds2);

	return (percolating>0)?(1):(0);

}

/*
	Probability changing cluster for a 2D XY bilayer
*/

double pcc_bilayer(int x,int y,double beta,double Jup,double Jdown,double K)
{
	struct bilayer_t *cfgt;
	double starting_delta,target_delta,delta,chi;
	double average;
	int c;

#define TOTAL_SWEEPS	(20000)
#define THERMALIZATION	(10000)
#define UPDATE_INTERVAL	(10)

	assert(TOTAL_SWEEPS>UPDATE_INTERVAL);

	cfgt=bilayer_init(x,y,Jup,Jdown,K);
	spin2d_random_configuration(cfgt->layers[LOWER_LAYER]);
	spin2d_random_configuration(cfgt->layers[UPPER_LAYER]);

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

		percolating=swendsen_wang_step_bilayer(cfgt,beta);

		beta+=((percolating==1)?(-delta):(+delta));

		if((c>0)&&((c%UPDATE_INTERVAL)==0))
			delta/=chi;
	}

	average=0.0f;
	for(c=0;c<THERMALIZATION;c++)
	{
		short percolating;

		percolating=swendsen_wang_step_bilayer(cfgt,beta);

		beta+=((percolating==1)?(-delta):(+delta));

		average+=beta;
	}

	bilayer_fini(cfgt);

	return average/((double)(THERMALIZATION));
}
