#ifndef __UTILS_H__
#define __UTILS_H__

#include <math.h>
#include <stdlib.h>
#include <unistd.h>

#ifndef MIN
#define MIN(x,y)	(((x)<(y))?(x):(y))
#endif

#ifndef MAX
#define MAX(x,y)	(((x)>(y))?(x):(y))
#endif

void init_prng(void);
double gen_random_number(void);
int gen_random_int(void);

#define MAKE_INDEX(ctx,x,y)	((x)+ctx->lx*(y))

struct vec2d_t
{
	double *vx,*vy;
	int lx,ly;
};

struct vec2d_t *vec2d_init(int x,int y);
void vec2d_fini(struct vec2d_t *v);
double vec2d_get_vx(struct vec2d_t *v,int x,int y);
double vec2d_get_vy(struct vec2d_t *v,int x,int y);
void vec2d_set_vector(struct vec2d_t *v,int x,int y,double vx,double vy);

struct spin2d_t
{
	double *theta;
	int lx,ly;
};

struct spin2d_t *spin2d_init(int x,int y);
void spin2d_fini(struct spin2d_t *s);
double spin2d_get_spin(struct spin2d_t *s,int x,int y);
void spin2d_set_spin(struct spin2d_t *s,int x,int y,double alpha);
void spin2d_clear_configuration(struct spin2d_t *cfgt);
void spin2d_random_configuration(struct spin2d_t *cfgt);

struct ising2d_t
{
	int *spins;
	int lx,ly;
};

struct ising2d_t *ising2d_init(int x,int y);
void ising2d_fini(struct ising2d_t *s);
int ising2d_get_spin(struct ising2d_t *s,int x,int y);
void ising2d_set_spin(struct ising2d_t *s,int x,int y,int spin);

#define DIR_X	(0)
#define DIR_Y	(1)

struct bond2d_t
{
	double *vals[2];
	int lx,ly;
};

struct bond2d_t *bond2d_init(int x,int y);
void bond2d_fini(struct bond2d_t *b);
double bond2d_get_value(struct bond2d_t *b,int x,int y,short direction);
void bond2d_set_value(struct bond2d_t *b,int x,int y,short direction,double value);

struct ibond2d_t
{
	int *vals[2];
	int lx,ly;
};

struct ibond2d_t *ibond2d_init(int x,int y);
void ibond2d_fini(struct ibond2d_t *b);
int ibond2d_get_value(struct ibond2d_t *b,int x,int y,short direction);
void ibond2d_set_value(struct ibond2d_t *b,int x,int y,short direction,int value);

struct vbond2d_t
{
	double *vals;
	int lx,ly;
};

struct vbond2d_t *vbond2d_init(int x,int y);
void vbond2d_fini(struct vbond2d_t *vb);
double vbond2d_get_value(struct vbond2d_t *vb,int x,int y);
void vbond2d_set_value(struct vbond2d_t *vb,int x,int y,double val);

struct ivbond2d_t
{
	int *vals;
	int lx,ly;
};

struct ivbond2d_t *ivbond2d_init(int x,int y);
void ivbond2d_fini(struct ivbond2d_t *vb);
int ivbond2d_get_value(struct ivbond2d_t *vb,int x,int y);
void ivbond2d_set_value(struct ivbond2d_t *vb,int x,int y,int val);

#endif
