/*
	xy.c
	
	Swendsen-Wang algorithm and probability changing cluster
	for a 2D XY model.
	
	Giacomo Bighin <bighin@gmail.com>

	Copyright (c) 2016, Giacomo Bighin
	All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions are met:

	1. Redistributions of source code must retain the above copyright
	   notice, this list of conditions and the following disclaimer.
	2. Redistributions in binary form must reproduce the above copyright
	   notice, this list of conditions and the following disclaimer in the
	   documentation and/or other materials provided with the distribution.
	3. All advertising materials mentioning features or use of this software
	   must display the following acknowledgement:
	   This product includes software developed by the <organization>.
	4. Neither the name of the <organization> nor the
	   names of its contributors may be used to endorse or promote products
	   derived from this software without specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER ''AS IS'' AND ANY
	EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
	WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
	DISCLAIMED. IN NO EVENT THE COPYRIGHT HOLDER SHALL BE LIABLE FOR ANY
	DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
	(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
	LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
	ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "xy.h"

/*
	Observables:

	- Energy on a single site
	- Total energy
	- Magnetization
	- Magnetization per spin
	- Energy per spin
*/

double lsite_energy(struct spin2d_t *cfgt,int x,int y)
{
	double theta1,theta2,xy=0.0f;

	theta1=spin2d_get_spin(cfgt,x,y);
	theta2=spin2d_get_spin(cfgt,(x+1)%cfgt->lx,y);

	xy-=cos(theta1-theta2);

	theta1=spin2d_get_spin(cfgt,x,y);
	theta2=spin2d_get_spin(cfgt,x,(y+1)%cfgt->ly);

	xy-=cos(theta1-theta2);

	return xy;
}

double configuration_energy(struct spin2d_t *cfgt)
{
	int i,j;
	double ret=0.0f;
	
	for(i=0;i<cfgt->lx;i++)
		for(j=0;j<cfgt->ly;j++)
			ret+=lsite_energy(cfgt,i,j);

	return ret;
}

double magnetization(struct spin2d_t *cfgt)
{
	int i,j;
	double a,b;

	a=b=0.0f;
	for(i=0;i<cfgt->lx;i++)
	{
		for(j=0;j<cfgt->ly;j++)
		{
			a+=cos(spin2d_get_spin(cfgt,i,j));
			b+=sin(spin2d_get_spin(cfgt,i,j));
		}
	}

	return sqrt(pow(a,2.0f)+pow(b,2.0f));
}

double magnetization_per_spin(struct spin2d_t *cfgt)
{
	double spins=cfgt->lx*cfgt->ly;
	
	return magnetization(cfgt)/spins;
}

double energy_per_spin(struct spin2d_t *cfgt)
{
	double spins=cfgt->lx*cfgt->ly;
	
	return configuration_energy(cfgt)/spins;
}

/*
	Wolff's embedded cluster.

	Given a two-dimensional spin-configuration in cfgt, the following function separates
	the orthogonal and parallel parts with respect to the unit vector defined by the angle
	alpha.

	Two corresponding Ising models are created.
*/

void spin2d_projection(struct spin2d_t *cfgt,double alpha,struct vec2d_t *par,struct vec2d_t *ortho,
                       struct ising2d_t *epsilon1,struct ising2d_t *epsilon2)
{
	int x,y;
	double r[2],q[2];
	
	r[0]=cos(alpha);
	r[1]=sin(alpha);

	q[0]=cos(alpha+M_PI/2.0f);
	q[1]=sin(alpha+M_PI/2.0f);

	for(x=0;x<cfgt->lx;x++)
	{
		for(y=0;y<cfgt->ly;y++)
		{
			double s[2],m,p[2],o[2];
			
			s[0]=cos(spin2d_get_spin(cfgt,x,y));
			s[1]=sin(spin2d_get_spin(cfgt,x,y));
		
			m=fabs(s[0]*r[0]+s[1]*r[1]);
		
			p[0]=m*r[0];
			p[1]=m*r[1];

			m=fabs(s[0]*q[0]+s[1]*q[1]);

			o[0]=m*q[0];
			o[1]=m*q[1];

			ising2d_set_spin(epsilon1,x,y,((s[0]*r[0]+s[1]*r[1])>0.0f)?(+1):(-1));
			ising2d_set_spin(epsilon2,x,y,((s[0]*q[0]+s[1]*q[1])>0.0f)?(+1):(-1));

			vec2d_set_vector(par,x,y,p[0],p[1]);
			vec2d_set_vector(ortho,x,y,o[0],o[1]);
		}
	}
}

/*
	Wolff's embedded cluster.

	Given the parallel and orthogonal parts we reconstruct the initial spin configuration.
*/

void spin2d_reconstruct(struct spin2d_t *cfgt,double alpha,struct vec2d_t *par,struct vec2d_t *ortho,
                        struct ising2d_t *epsilon1,struct ising2d_t *epsilon2)
{
	int x,y;

	for(x=0;x<cfgt->lx;x++)
	{
		for(y=0;y<cfgt->ly;y++)
		{
			double s[2],p[2],o[2];

			p[0]=vec2d_get_vx(par,x,y);
			p[1]=vec2d_get_vy(par,x,y);

			o[0]=vec2d_get_vx(ortho,x,y);
			o[1]=vec2d_get_vy(ortho,x,y);

			if(ising2d_get_spin(epsilon1,x,y)==-1)
			{
				p[0]=-p[0];
				p[1]=-p[1];
			}

			if(ising2d_get_spin(epsilon2,x,y)==-1)
			{
				o[0]=-o[0];
				o[1]=-o[1];
			}

			s[0]=o[0]+p[0];
			s[1]=o[1]+p[1];

			spin2d_set_spin(cfgt,x,y,atan2(s[1],s[0]));
		}
	}
}

/*
	A Swendsen-Wang step for the Ising model
*/

short swendsen_wang_ising_step(struct ising2d_t *epsilon,struct bond2d_t *Jij,double beta)
{
	struct ibond2d_t *bonds;
	struct clusters_t *clusters;

	int c,x,y;
	int nr_clusters;
	short percolating;

	/*
		A partire da Jij ora costruisco in modo stocastico i bond tra i vari siti.
	*/

	bonds=ibond2d_init(epsilon->lx,epsilon->ly);

	for(x=0;x<epsilon->lx-1;x++)
	{
		for(y=0;y<epsilon->ly-1;y++)
		{
			ibond2d_set_value(bonds,x,y,DIR_X,0);
			ibond2d_set_value(bonds,x,y,DIR_Y,0);

			/*
				Se due spin adiacenti lungo x sono concordi allora
				attivo il bond con probabilità p.
			*/

			if(ising2d_get_spin(epsilon,x,y)==ising2d_get_spin(epsilon,x+1,y))
			{
				double p=1.0f-exp(-2.0f*beta*bond2d_get_value(Jij,x,y,DIR_X));

				if(gen_random_number()<p)
					ibond2d_set_value(bonds,x,y,DIR_X,1);
			}

			/*
				Stessa cosa lungo y.
			*/

			if(ising2d_get_spin(epsilon,x,y)==ising2d_get_spin(epsilon,x,y+1))
			{
				double p=1.0f-exp(-2.0f*beta*bond2d_get_value(Jij,x,y,DIR_Y));

				if(gen_random_number()<p)
					ibond2d_set_value(bonds,x,y,DIR_Y,1);
			}
		}
	}

	/*
		Identifico i clusters
	*/

	clusters=clusters_init(epsilon->lx,epsilon->ly);
	nr_clusters=ising2d_identify_clusters(bonds,clusters);

	/*
		Valuto le dimensioni di ciascun cluster e stabilisco se c'è stata percolazione
	*/

	percolating=0;
	for(c=1;c<nr_clusters;c++)
	{
		int xdim,ydim;

		cluster_dimensions(epsilon,clusters,c,&xdim,&ydim);

		/*
			If a cluster has the right dimensions (extension rule)
			we check if the cluster is connected using the periodic
			boundary conditions (topological rule)
		*/

#define TOPOLOGICAL
#ifdef TOPOLOGICAL
		if(xdim==epsilon->lx)
		{	
			for(y=0;y<epsilon->ly;y++)
				if(clusters_get_value(clusters,0,y)==clusters_get_value(clusters,epsilon->lx-1,y))
					percolating=1;
		}

		if(ydim==epsilon->ly)
		{	
			for(x=0;x<epsilon->lx;x++)
				if(clusters_get_value(clusters,x,0)==clusters_get_value(clusters,x,epsilon->ly-1))
					percolating=1;
		}
#else
		if((xdim==epsilon->lx)||(ydim==epsilon->ly))
			percolating=1;
#endif
	}

	/*
		Flippo i clusters appena creati
	*/
		
	for(c=1;c<nr_clusters;c++)
	{
		if(gen_random_number()<0.5f)
			reset_cluster(epsilon,clusters,c,1);
		else
			reset_cluster(epsilon,clusters,c,-1);
	}
	
	ibond2d_fini(bonds);
	clusters_fini(clusters);
	
	return percolating;
}

/*
	A Swendsen-Wang step for the 2D XY model.

	See also "Monte Carlo Simulations of Spin Systems", by Wolfhard Janke
*/

short swendsen_wang_step(struct spin2d_t *cfgt,double beta,double J)
{
	struct vec2d_t *parallel,*orthogonal;
	struct ising2d_t *epsilon1,*epsilon2;
	struct bond2d_t *Jij1,*Jij2;
	
	double alpha;
	int x,y;
	short percolating;

	/*
		Scelgo un vettore unitario casuale r, definito dall'angolo alpha
	*/

	alpha=gen_random_number()*2.0f*M_PI;

	/*
		Wolff's embedded cluster

		Separo i gradi di libertà di spin paralleli e ortogonali al vettore unitario
		e costruisco un modello di Ising a partire dai gradi di spin paralleli a r.
	*/

	epsilon1=ising2d_init(cfgt->lx,cfgt->ly);
	epsilon2=ising2d_init(cfgt->lx,cfgt->ly);

	parallel=vec2d_init(cfgt->lx,cfgt->ly);
	orthogonal=vec2d_init(cfgt->lx,cfgt->ly);

	spin2d_projection(cfgt,alpha,parallel,orthogonal,epsilon1,epsilon2);

	/*
		Dapprima costruisco le variabil Jij, definite su ciascun bond
	*/

	Jij1=bond2d_init(cfgt->lx,cfgt->ly);
	Jij2=bond2d_init(cfgt->lx,cfgt->ly);

	for(x=0;x<cfgt->lx-1;x++)
	{
		for(y=0;y<cfgt->ly-1;y++)
		{
			double si[2],sj[2],r[2],q[2],z1,z2;

			r[0]=cos(alpha);
			r[1]=sin(alpha);

			q[0]=cos(alpha+M_PI/2.0f);
			q[1]=sin(alpha+M_PI/2.0f);

			si[0]=cos(spin2d_get_spin(cfgt,x,y));
			si[1]=sin(spin2d_get_spin(cfgt,x,y));

			sj[0]=cos(spin2d_get_spin(cfgt,x+1,y));
			sj[1]=sin(spin2d_get_spin(cfgt,x+1,y));

			z1=J*fabs(si[0]*r[0]+si[1]*r[1])*fabs(sj[0]*r[0]+sj[1]*r[1]);
			z2=J*fabs(si[0]*q[0]+si[1]*q[1])*fabs(sj[0]*q[0]+sj[1]*q[1]);

			bond2d_set_value(Jij1,x,y,DIR_X,z1);
			bond2d_set_value(Jij2,x,y,DIR_X,z2);

			sj[0]=cos(spin2d_get_spin(cfgt,x,y+1));
			sj[1]=sin(spin2d_get_spin(cfgt,x,y+1));

			z1=J*fabs(si[0]*r[0]+si[1]*r[1])*fabs(sj[0]*r[0]+sj[1]*r[1]);
			z2=J*fabs(si[0]*q[0]+si[1]*q[1])*fabs(sj[0]*q[0]+sj[1]*q[1]);

			bond2d_set_value(Jij1,x,y,DIR_Y,z1);
			bond2d_set_value(Jij2,x,y,DIR_Y,z2);

		}
	}

	/*
		Avendo creato due modelli di Ising effettivi uso l'algoritmo di Swendsen-Wang su quelli.
	*/

	percolating=swendsen_wang_ising_step(epsilon1,Jij1,beta);
	percolating+=swendsen_wang_ising_step(epsilon2,Jij2,beta);

	/*
		Ricostruisco gli spin iniziali.
	*/

	spin2d_reconstruct(cfgt,alpha,parallel,orthogonal,epsilon1,epsilon2);

	/*
		Cleaning up and return the percolation status.
	*/

	vec2d_fini(parallel);
	vec2d_fini(orthogonal);

	ising2d_fini(epsilon1);
	ising2d_fini(epsilon2);
	bond2d_fini(Jij1);
	bond2d_fini(Jij2);
	
	return (percolating>0)?(1):(0);
}
