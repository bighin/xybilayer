#ifndef __BILAYER_H__
#define __BILAYER_H__

#include "xy.h"
#include "stat.h"

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

void wolff_embedded_cluster(struct spin2d_t *cfgt,struct vec2d_t **parallel,struct vec2d_t **orthogonal,
                            double alpha,struct ising2d_t **epsilon1,struct ising2d_t **epsilon2,
			    struct bond2d_t **Jij1,struct bond2d_t **Jij2,double J);

short swendsen_wang_ising_bilayer_step(struct ising2d_t *epsilon[2],struct bond2d_t *Jij[2],
                                       struct vbond2d_t *vbonds,double beta);

short swendsen_wang_step_bilayer(struct bilayer_t *cfgt,double beta);
double pcc_bilayer(int x,int y,double beta,double Jup,double Jdown,double K);

/*
	Number of correlators (complex entries are counted twice)
	and number of scalar samples.
*/

#define NR_CORRELATORS	(6)
#define NR_SCALARS	(0)

int get_total_channels(int maxk);
int get_channel_nr(int maxk,char *desc,int k);
struct sampling_ctx_t *sw_bilayer(int x,int y,double beta,double Jup,double Jdown,double K,int maxk);

#endif
