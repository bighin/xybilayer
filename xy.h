#ifndef __XY_H__
#define __XY_H__

#include "utils.h"
#include "clusters.h"

double lsite_energy(struct spin2d_t *cfgt,int x,int y);
double configuration_energy(struct spin2d_t *cfgt);
double magnetization(struct spin2d_t *cfgt);
double magnetization_per_spin(struct spin2d_t *cfgt);
double energy_per_spin(struct spin2d_t *cfgt);

void spin2d_projection(struct spin2d_t *cfgt,double alpha,struct vec2d_t *par,struct vec2d_t *ortho,
                       struct ising2d_t *epsilon1,struct ising2d_t *epsilon2);

void spin2d_reconstruct(struct spin2d_t *cfgt,double alpha,struct vec2d_t *par,struct vec2d_t *ortho,
                        struct ising2d_t *epsilon1,struct ising2d_t *epsilon2);

short swendsen_wang_ising_step(struct ising2d_t *epsilon,struct bond2d_t *Jij,double beta);
short swendsen_wang_step(struct spin2d_t *cfgt,double beta,double J);

#endif
