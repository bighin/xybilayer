/*
	Monte Carlo simulation of a two-dimensional XY bilayer model
*/

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "bilayer.h"
#include "xy.h"
#include "common.h"
#include "stat.h"

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

		After choosing a reference vector, the degrees of freedom of the system are split
		as parallel and orthogonal to the vector, leading to two effective Ising models.
	*/

	*epsilon1=ising2d_init(cfgt->lx,cfgt->ly);
	*epsilon2=ising2d_init(cfgt->lx,cfgt->ly);

	*parallel=vec2d_init(cfgt->lx,cfgt->ly);
	*orthogonal=vec2d_init(cfgt->lx,cfgt->ly);

	spin2d_projection(cfgt,alpha,*parallel,*orthogonal,*epsilon1,*epsilon2);

	/*
		We first build the Jij variables, defined on each bond
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
		Now the Jij's are used to stochastically activate bonds between sites.
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
						If two spins are first neighbours along the x direction and have the same sign,
						the bond is activated with probability p.
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
					ibond2d_set_value(bonds[c],x,y,DIR_Y,0);

					/*
						Same along the y direction.
					*/

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
		Here vertical bonds are created, as well.
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
		Let us identify the clusters
	*/


	bclusters=bclusters_init(epsilon[LOWER_LAYER]->lx,epsilon[LOWER_LAYER]->ly);
	nr_clusters=ising2d_identify_bclusters(bonds,ivbonds,bclusters);

	/*
		We evaluate the dimensions of each cluster, establishing if percolation has happened.
	*/

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

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

		if(this_percolating==1)
			percolating=1;
	}

#pragma GCC diagnostic pop

	/*
		The clusters just created are flipped
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
		We choose a random two-dimensional unitary vector, defined by the angle \alpha
	*/

	alpha=gen_random_number()*2.0f*M_PI;

	/*
		The embedded cluster algorithm is applied to each layer.
	*/

	for(l=0;l<=1;l++)
	{
		wolff_embedded_cluster(cfgt->layers[l],&parallel[l],&orthogonal[l],alpha,
	                               &epsilon1[l],&epsilon2[l],&Jij1[l],&Jij2[l],cfgt->J[l]);
	}

	/*
		We also need to build the vertical bonds.
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
		After having obtained two effective Ising models, we apply the Swendsen-Wang on each one.
	*/

	percolating=swendsen_wang_ising_bilayer_step(epsilon1,Jij1,vbonds1,beta);
	percolating+=swendsen_wang_ising_bilayer_step(epsilon2,Jij2,vbonds2,beta);

	/*
		The initial spins are reconstructed.
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
	Probability changing cluster for a two-dimensional XY bilayer
*/

double pcc_bilayer(int x,int y,double beta,double Jup,double Jdown,double K)
{
	struct bilayer_t *cfgt;
	double starting_delta,target_delta,delta,chi;
	double average;
	int c;

#define PCC_THERMALIZATION		(20000)
#define PCC_POST_THERMALIZATION		(10000)
#define PCC_UPDATE_INTERVAL		(10)

	assert(PCC_THERMALIZATION>PCC_UPDATE_INTERVAL);

	cfgt=bilayer_init(x,y,Jup,Jdown,K);
	spin2d_random_configuration(cfgt->layers[LOWER_LAYER]);
	spin2d_random_configuration(cfgt->layers[UPPER_LAYER]);

	starting_delta=0.1f;
	target_delta=0.000680f;

	/*
		\delta (the \beta variation after each step) is adjusted each UPDATE_INTERVAL steps.

		The adjustment is performed such that in the beginning \delta is starting_delta,
		whereas at the end is target_delta.
	*/

	delta=starting_delta;
	chi=pow(target_delta/starting_delta,-1.0f/(PCC_THERMALIZATION/PCC_UPDATE_INTERVAL));

	for(c=0;c<PCC_THERMALIZATION;c++)
	{
		short percolating;

		percolating=swendsen_wang_step_bilayer(cfgt,beta);

		beta+=((percolating==1)?(-delta):(+delta));

		if((c>0)&&((c%PCC_UPDATE_INTERVAL)==0))
			delta/=chi;
	}

	average=0.0f;
	for(c=0;c<PCC_POST_THERMALIZATION;c++)
	{
		short percolating;

		percolating=swendsen_wang_step_bilayer(cfgt,beta);

		beta+=((percolating==1)?(-delta):(+delta));

		average+=beta;
	}

	bilayer_fini(cfgt);

	return average/((double)(PCC_POST_THERMALIZATION));
}

/*
	c(k), i.e. the in-plane correlator, mathematically:

	\sum \exp(i psi_i - i psi_j)

	where k is the distance between i and j, while the psi are
	the XY variables defined either on the upper or on the lower layer.

	The sum extends over all pairs of sites (i, j) at a distance k,
	taking into account periodic boundary conditions.
*/

int ck(struct bilayer_t *cfgt,int k,short layer,double *res)
{
	int i,j;
	double rc,ic;
	
	rc=ic=0.0f;

	for(i=0;i<cfgt->lx;i++)
	{
		for(j=0;j<cfgt->ly;j++)
		{
			double psii,psij,phase;
			double lrc,lic;
			int ip,jp;

			lrc=lic=0;
		
			psii=spin2d_get_spin(cfgt->layers[layer],i,j);
			
			ip=(i+k)%cfgt->lx;
			jp=j;
			
			psij=spin2d_get_spin(cfgt->layers[layer],ip,jp);

			phase=psii-psij;
			lrc+=cos(phase);
			lic+=sin(phase);

			ip=i;
			jp=(j+k)%cfgt->ly;

			psij=spin2d_get_spin(cfgt->layers[layer],ip,jp);

			phase=psii-psij;
			lrc+=cos(phase);
			lic+=sin(phase);

			ip=(i+cfgt->lx-k)%cfgt->lx;
			jp=j;

			psij=spin2d_get_spin(cfgt->layers[layer],ip,jp);

			phase=psii-psij;
			lrc+=cos(phase);
			lic+=sin(phase);

			ip=i;
			jp=(j+cfgt->ly-k)%cfgt->ly;

			psij=spin2d_get_spin(cfgt->layers[layer],ip,jp);

			phase=psii-psij;
			lrc+=cos(phase);
			lic+=sin(phase);		
			
			lrc/=4.0f;
			lic/=4.0f;
			
			rc+=lrc;
			ic+=lic;
		}
	}

	rc/=cfgt->lx*cfgt->ly;
	ic/=cfgt->lx*cfgt->ly;
	
	res[0]=rc;
	res[1]=ic;

	return 0;
}


/*
	z(k), i.e. the 'anomalous' correlator proposed by Andrea:

	\sum \exp(i phi_i + i psi_i - i phi_j - i psi_j)

	where k is the distance between i and j, while phi e psi are XY variables
	defined, respectively, on the lower and upper layer.
	
	The sum extends over all pairs of sites (i, j) at a distance k,
	taking into account periodic boundary conditions.
*/

int zk(struct bilayer_t *cfgt,int k,double *res)
{
	int i,j;
	double rz,iz;
	
	rz=iz=0.0f;

	for(i=0;i<cfgt->lx;i++)
	{
		for(j=0;j<cfgt->ly;j++)
		{
			double psii,psij,phii,phij,phase;
			double lrz,liz;
			int ip,jp;

			lrz=liz=0;
		
			psii=spin2d_get_spin(cfgt->layers[LOWER_LAYER],i,j);
			phii=spin2d_get_spin(cfgt->layers[UPPER_LAYER],i,j);
			
			ip=(i+k)%cfgt->lx;
			jp=j;
			
			psij=spin2d_get_spin(cfgt->layers[LOWER_LAYER],ip,jp);
			phij=spin2d_get_spin(cfgt->layers[UPPER_LAYER],ip,jp);

			phase=phii+psii-phij-psij;
			lrz+=cos(phase);
			liz+=sin(phase);

			ip=i;
			jp=(j+k)%cfgt->ly;

			psij=spin2d_get_spin(cfgt->layers[LOWER_LAYER],ip,jp);
			phij=spin2d_get_spin(cfgt->layers[UPPER_LAYER],ip,jp);

			phase=phii+psii-phij-psij;
			lrz+=cos(phase);
			liz+=sin(phase);

			ip=(i+cfgt->lx-k)%cfgt->lx;
			jp=j;

			psij=spin2d_get_spin(cfgt->layers[LOWER_LAYER],ip,jp);
			phij=spin2d_get_spin(cfgt->layers[UPPER_LAYER],ip,jp);

			phase=phii+psii-phij-psij;
			lrz+=cos(phase);
			liz+=sin(phase);

			ip=i;
			jp=(j+cfgt->ly-k)%cfgt->ly;

			psij=spin2d_get_spin(cfgt->layers[LOWER_LAYER],ip,jp);
			phij=spin2d_get_spin(cfgt->layers[UPPER_LAYER],ip,jp);

			phase=phii+psii-phij-psij;
			lrz+=cos(phase);
			liz+=sin(phase);		
			
			lrz/=4.0f;
			liz/=4.0f;
			
			rz+=lrz;
			iz+=liz;
		}
	}

	rz/=cfgt->lx*cfgt->ly;
	iz/=cfgt->lx*cfgt->ly;
	
	res[0]=rz;
	res[1]=iz;

	return 0;
}

int get_total_channels(int maxk)
{
	return maxk*NR_CORRELATORS+NR_SCALARS;
}

int get_channel_nr(int maxk,char *desc,int k)
{
	if(strstr(desc,"zk"))
	{
		if(strstr(desc,"real"))
			return 2*k;

		if(strstr(desc,"imag"))
			return 2*k+1;

		return -1;
	}

	if(strstr(desc,"ck"))
	{
		int baseline;
	
		if(strstr(desc,"upper"))
			baseline=2*maxk;
		else if(strstr(desc,"lower"))
			baseline=4*maxk;
		else
			return -1;

		if(strstr(desc,"real"))
			return baseline+2*k;

		if(strstr(desc,"imag"))
			return baseline+2*k+1;

		return -1;
	}
	
	return -1;
}

/*
	Swendsen-Wang algorithm for a bilayer.

	After thermalization the correlation functions c(k) and z(k) are given as a output.
*/

struct sampling_ctx_t *sw_bilayer(int x,int y,double beta,double Jup,double Jdown,double K,int maxk)
{
	struct bilayer_t *cfgt;
	struct sampling_ctx_t *sctx;
	int c,k;

	assert(maxk>0);
	assert(maxk<=x);
	assert(maxk<=y);

	/*
		Initialization of various struct's
	*/

	cfgt=bilayer_init(x,y,Jup,Jdown,K);
	spin2d_random_configuration(cfgt->layers[LOWER_LAYER]);
	spin2d_random_configuration(cfgt->layers[UPPER_LAYER]);

	/*
		There are 6*maxk channels in the sampling context, as we are going to take
		samples for three correlators, real and imaginary parts, and every correlator
		will be sampled from k=0 to k=(maxk-1)
	
		Notare che se cambio questo valore devo cambiare anche la funzione channel_nr()
	*/

	sctx=sampling_ctx_init(get_total_channels(maxk));

	/*
		Note that the SW_THERMALIZATION and SW_POST_THERMALIZATION values
		are just tentative and need to be adjusted, in particular as a
		function of the lattice dimensions.
	*/

#define SW_THERMALIZATION		(1000)
#define SW_POST_THERMALIZATION		(500)

	for(c=0;c<SW_THERMALIZATION;c++)
	{
		swendsen_wang_step_bilayer(cfgt,beta);
	}

	for(c=0;c<SW_POST_THERMALIZATION;c++)
	{
		double z[2],c[2][2];

		swendsen_wang_step_bilayer(cfgt,beta);

		for(k=0;k<maxk;k++)
		{
			z[0]=z[1]=0.0f;
			c[UPPER_LAYER][0]=c[UPPER_LAYER][1]=0.0f;
			c[LOWER_LAYER][0]=c[LOWER_LAYER][1]=0.0f;

			zk(cfgt,k,z);
			ck(cfgt,k,UPPER_LAYER,c[UPPER_LAYER]);
			ck(cfgt,k,LOWER_LAYER,c[LOWER_LAYER]);
		
			sampling_ctx_add_entry_to_channel(sctx,get_channel_nr(maxk,"zk real",k),z[0]);
			sampling_ctx_add_entry_to_channel(sctx,get_channel_nr(maxk,"zk imag",k),z[1]);

			sampling_ctx_add_entry_to_channel(sctx,get_channel_nr(maxk,"ck upper real",k),c[UPPER_LAYER][0]);
			sampling_ctx_add_entry_to_channel(sctx,get_channel_nr(maxk,"ck upper imag",k),c[UPPER_LAYER][1]);

			sampling_ctx_add_entry_to_channel(sctx,get_channel_nr(maxk,"ck lower real",k),c[LOWER_LAYER][0]);
			sampling_ctx_add_entry_to_channel(sctx,get_channel_nr(maxk,"ck lower imag",k),c[LOWER_LAYER][1]);
		}
	}

	bilayer_fini(cfgt);

	return sctx;
}
