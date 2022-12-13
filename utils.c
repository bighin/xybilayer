#include <stdio.h>

#include "utils.h"
#include "common.h"

void init_prng(void)
{
	srand48(getpid());
}

double gen_random_number(void)
{
	return drand48();
}

int gen_random_int(void)
{
	return lrand48();
}

/*
	This struct represents a two-dimensional vector lattice
*/

struct vec2d_t *vec2d_init(int x,int y)
{
	struct vec2d_t *ret;
	
	if(!(ret=malloc(sizeof(struct vec2d_t))))
		return NULL;
	
	ret->vx=malloc(sizeof(double)*x*y);
	ret->vy=malloc(sizeof(double)*x*y);

	if((!ret->vx)||(!ret->vy))
	{
		if(ret->vx)
			free(ret->vx);

		if(ret->vy)
			free(ret->vy);

		if(ret)
			free(ret);
		
		return NULL;
	}
	
	ret->lx=x;
	ret->ly=y;
	
	return ret;
}

void vec2d_fini(struct vec2d_t *v)
{
	if(v)
	{
		if(v->vx)
			free(v->vx);

		if(v->vy)
			free(v->vy);
		
		free(v);
	}
}

double vec2d_get_vx(struct vec2d_t *v,int x,int y)
{
	return v->vx[MAKE_INDEX(v,x,y)];
}

double vec2d_get_vy(struct vec2d_t *v,int x,int y)
{
	return v->vy[MAKE_INDEX(v,x,y)];
}

void vec2d_set_vector(struct vec2d_t *v,int x,int y,double vx,double vy)
{
	v->vx[MAKE_INDEX(v,x,y)]=vx;
	v->vy[MAKE_INDEX(v,x,y)]=vy;
}

/*
	This struct represents a two-dimensional spin configuration
*/

struct spin2d_t *spin2d_init(int x,int y)
{
	struct spin2d_t *ret;
	
	if(!(ret=malloc(sizeof(struct spin2d_t))))
		return NULL;
	
	ret->theta=malloc(sizeof(double)*x*y);

	if(!ret->theta)
	{
		if(ret)
			free(ret);
		
		return NULL;
	}
	
	ret->lx=x;
	ret->ly=y;
	
	return ret;
}

void spin2d_fini(struct spin2d_t *s)
{
	if(s)
	{
		if(s->theta)
			free(s->theta);
		
		free(s);
	}
}

double spin2d_get_spin(struct spin2d_t *s,int x,int y)
{
	return s->theta[MAKE_INDEX(s,x,y)];
}

void spin2d_set_spin(struct spin2d_t *s,int x,int y,double alpha)
{
	s->theta[MAKE_INDEX(s,x,y)]=alpha;
}

void spin2d_clear_configuration(struct spin2d_t *cfgt)
{
	int i,j;
	
	for(i=0;i<cfgt->lx;i++)
		for(j=0;j<cfgt->ly;j++)
			spin2d_set_spin(cfgt,i,j,0.0f);
}

void spin2d_random_configuration(struct spin2d_t *cfgt)
{
	int i,j;
	
	for(i=0;i<cfgt->lx;i++)
		for(j=0;j<cfgt->ly;j++)
			spin2d_set_spin(cfgt,i,j,gen_random_number()*2.0f*M_PI);
}

/*
	This struct represents a two-dimensional Ising-like configuration.
*/

struct ising2d_t *ising2d_init(int x,int y)
{
	struct ising2d_t *ret;
	
	if(!(ret=malloc(sizeof(struct ising2d_t))))
		return NULL;
	
	ret->spins=malloc(sizeof(int)*x*y);

	if(!ret->spins)
	{
		if(ret)
			free(ret);
		
		return NULL;
	}
	
	ret->lx=x;
	ret->ly=y;
	
	return ret;
}

void ising2d_fini(struct ising2d_t *s)
{
	if(s)
	{
		if(s->spins)
			free(s->spins);
		
		free(s);
	}
}

int ising2d_get_spin(struct ising2d_t *s,int x,int y)
{
	return s->spins[MAKE_INDEX(s,x,y)];
}

void ising2d_set_spin(struct ising2d_t *s,int x,int y,int spin)
{
	s->spins[MAKE_INDEX(s,x,y)]=spin;
}

/*
	A floating point quantity defined on each bond in a two-dimensional lattice
*/

struct bond2d_t *bond2d_init(int x,int y)
{
	struct bond2d_t *ret;
	
	if(!(ret=malloc(sizeof(struct bond2d_t))))
		return NULL;
	
	ret->vals[0]=malloc(sizeof(double)*x*y);
	ret->vals[1]=malloc(sizeof(double)*x*y);

	if((!ret->vals[0])||(!ret->vals[1]))
	{
		if(ret->vals[0])
			free(ret->vals[0]);

		if(ret->vals[1])
			free(ret->vals[1]);

		if(ret)
			free(ret);
		
		return NULL;
	}
	
	ret->lx=x;
	ret->ly=y;
	
	return ret;
}

void bond2d_fini(struct bond2d_t *b)
{
	if(b)
	{
		if(b->vals[0])
			free(b->vals[0]);

		if(b->vals[1])
			free(b->vals[1]);
		
		free(b);
	}
}

double bond2d_get_value(struct bond2d_t *b,int x,int y,short direction)
{
	return b->vals[direction][MAKE_INDEX(b,x,y)];
}

void bond2d_set_value(struct bond2d_t *b,int x,int y,short direction,double value)
{
	b->vals[direction][MAKE_INDEX(b,x,y)]=value;
}

/*
	An integer quantity defined on each bond in a two-dimensional lattice
*/

struct ibond2d_t *ibond2d_init(int x,int y)
{
	struct ibond2d_t *ret;
	
	if(!(ret=malloc(sizeof(struct ibond2d_t))))
		return NULL;
	
	ret->vals[0]=malloc(sizeof(int)*x*y);
	ret->vals[1]=malloc(sizeof(int)*x*y);

	if((!ret->vals[0])||(!ret->vals[1]))
	{
		if(ret->vals[0])
			free(ret->vals[0]);

		if(ret->vals[1])
			free(ret->vals[1]);

		if(ret)
			free(ret);
		
		return NULL;
	}
	
	ret->lx=x;
	ret->ly=y;
	
	return ret;
}

void ibond2d_fini(struct ibond2d_t *b)
{
	if(b)
	{
		if(b->vals[0])
			free(b->vals[0]);

		if(b->vals[1])
			free(b->vals[1]);
		
		free(b);
	}
}

int ibond2d_get_value(struct ibond2d_t *b,int x,int y,short direction)
{
	return b->vals[direction][MAKE_INDEX(b,x,y)];
}

void ibond2d_set_value(struct ibond2d_t *b,int x,int y,short direction,int value)
{
	b->vals[direction][MAKE_INDEX(b,x,y)]=value;
}

/*
	vbond2d_t and ivbond2d_t represents a lattice of vertical bond
	variables (defined between the upper and lower layers)
	taking double and integer values, respectively.
*/

struct vbond2d_t *vbond2d_init(int x,int y)
{
	struct vbond2d_t *ret;
	
	if(!(ret=malloc(sizeof(struct vbond2d_t))))
		return NULL;
	
	ret->vals=malloc(sizeof(double)*x*y);

	if(!ret->vals)
	{
		if(ret)
			free(ret);
		
		return NULL;
	}
	
	ret->lx=x;
	ret->ly=y;
	
	return ret;
}

void vbond2d_fini(struct vbond2d_t *vb)
{
	if(vb)
	{
		if(vb->vals)
			free(vb->vals);
		
		free(vb);
	}
}

double vbond2d_get_value(struct vbond2d_t *vb,int x,int y)
{
	return vb->vals[MAKE_INDEX(vb,x,y)];
}

void vbond2d_set_value(struct vbond2d_t *vb,int x,int y,double val)
{
	vb->vals[MAKE_INDEX(vb,x,y)]=val;
}

struct ivbond2d_t *ivbond2d_init(int x,int y)
{
	struct ivbond2d_t *ret;
	
	if(!(ret=malloc(sizeof(struct ivbond2d_t))))
		return NULL;
	
	ret->vals=malloc(sizeof(int)*x*y);

	if(!ret->vals)
	{
		if(ret)
			free(ret);
		
		return NULL;
	}
	
	ret->lx=x;
	ret->ly=y;
	
	return ret;
}

void ivbond2d_fini(struct ivbond2d_t *vb)
{
	if(vb)
	{
		if(vb->vals)
			free(vb->vals);
		
		free(vb);
	}
}

int ivbond2d_get_value(struct ivbond2d_t *vb,int x,int y)
{
	return vb->vals[MAKE_INDEX(vb,x,y)];
}

void ivbond2d_set_value(struct ivbond2d_t *vb,int x,int y,int val)
{
	vb->vals[MAKE_INDEX(vb,x,y)]=val;
}
