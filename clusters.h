#ifndef __CLUSTERS_H__
#define __CLUSTERS_H__

#include "utils.h"
#include "bilayer.h"

struct clusters_t
{
	int *vals;
	int lx,ly;
};

struct clusters_t *clusters_init(int x,int y);
void clusters_fini(struct clusters_t *clusters);
int clusters_get_value(struct clusters_t *clusters,int x,int y);
void clusters_set_value(struct clusters_t *clusters,int x,int y,int value);
void ising2d_grow_clusters(struct ibond2d_t *bonds,struct clusters_t *clusters,int x,int y,int nr_clusters);
int ising2d_identify_clusters(struct ibond2d_t *bonds,struct clusters_t *clusters);
void cluster_dimensions(struct ising2d_t *epsilon,struct clusters_t *clusters,int target,int *xdim,int *ydim);
void reset_cluster(struct ising2d_t *epsilon,struct clusters_t *clusters,int target,int value);

/*
	Clusters on a bilayer!
*/

struct bclusters_t
{
	int *vals[2];
	int lx,ly;
};

struct bclusters_t *bclusters_init(int x,int y);
void bclusters_fini(struct bclusters_t *bc);
int bclusters_get_value(struct bclusters_t *bclusters,int x,int y,int layer);
void bclusters_set_value(struct bclusters_t *bclusters,int x,int y,int layer,int value);
void ising2d_grow_bclusters(struct ibond2d_t *bonds[2],struct ivbond2d_t *ivbonds,struct bclusters_t *bclusters,int x,int y,int l,int nr_clusters);
int ising2d_identify_bclusters(struct ibond2d_t *bonds[2],struct ivbond2d_t *ivbonds,struct bclusters_t *bclusters);
void bcluster_dimensions(struct ising2d_t *epsilon[2],struct bclusters_t *bclusters,int target,int *xdim,int *ydim);
void reset_bcluster(struct ising2d_t *epsilon[2],struct bclusters_t *bclusters,int target,int value);

#endif
