#include "utils.h"
#include "clusters.h"

/*
	A struct for identifying and classify clusters on a Ising model, given a bond configuration,
	as in the Swendsen-Wang approach.
*/

struct clusters_t *clusters_init(int x,int y)
{
	struct clusters_t *ret;
	
	if(!(ret=malloc(sizeof(struct clusters_t))))
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

void clusters_fini(struct clusters_t *clusters)
{
	if(clusters)
	{
		if(clusters->vals)
			free(clusters->vals);
		
		free(clusters);
	}
}

int clusters_get_value(struct clusters_t *clusters,int x,int y)
{
	return clusters->vals[MAKE_INDEX(clusters,x,y)];
}

void clusters_set_value(struct clusters_t *clusters,int x,int y,int value)
{
	clusters->vals[MAKE_INDEX(clusters,x,y)]=value;
}

void ising2d_grow_clusters(struct ibond2d_t *bonds,struct clusters_t *clusters,int x,int y,int nr_clusters)
{
	if(clusters_get_value(clusters,x,y)!=0)
		return;

	clusters_set_value(clusters,x,y,nr_clusters);

	if((ibond2d_get_value(bonds,x,y,DIR_X)==1)&&((x+1)<clusters->lx))
		ising2d_grow_clusters(bonds,clusters,x+1,y,nr_clusters);

	if((ibond2d_get_value(bonds,x,y,DIR_Y)==1)&&((y+1)<clusters->ly))
		ising2d_grow_clusters(bonds,clusters,x,y+1,nr_clusters);

	if(x>0)
		if(ibond2d_get_value(bonds,x-1,y,DIR_X)==1)
			ising2d_grow_clusters(bonds,clusters,x-1,y,nr_clusters);

	if(y>0)
		if(ibond2d_get_value(bonds,x,y-1,DIR_Y)==1)
			ising2d_grow_clusters(bonds,clusters,x,y-1,nr_clusters);
}

int ising2d_identify_clusters(struct ibond2d_t *bonds,struct clusters_t *clusters)
{
	int x,y;
	int nr_clusters=1;

	for(x=0;x<clusters->lx;x++)
		for(y=0;y<clusters->ly;y++)
			clusters_set_value(clusters,x,y,0);

	for(x=0;x<clusters->lx;x++)
	{
		for(y=0;y<clusters->ly;y++)
		{
			if(clusters_get_value(clusters,x,y)==0)
			{
				ising2d_grow_clusters(bonds,clusters,x,y,nr_clusters);

				nr_clusters++;
			}
		}
	}
	
	return nr_clusters;
}

void cluster_dimensions(struct ising2d_t *epsilon,struct clusters_t *clusters,int target,int *xdim,int *ydim)
{
	int maxx,maxy,minx,miny;
	int x,y;
	
	maxx=maxy=0;
	minx=epsilon->lx;
	miny=epsilon->ly;

	for(x=0;x<epsilon->lx;x++)
	{
		for(y=0;y<epsilon->ly;y++)
		{
			if(clusters->vals[MAKE_INDEX(clusters,x,y)]==target)
			{
				if(x>maxx)
					maxx=x;

				if(y>maxy)
					maxy=y;

				if(x<minx)
					minx=x;

				if(y<miny)
					miny=y;
			}
		}
	}

	*xdim=maxx-minx+1;
	*ydim=maxy-miny+1;

	if(*xdim<0)
		*xdim=0;
	
	if(*ydim<0)
		*ydim=0;
}

void reset_cluster(struct ising2d_t *epsilon,struct clusters_t *clusters,int target,int value)
{
	int x,y;
		
	for(x=0;x<epsilon->lx;x++)
		for(y=0;y<epsilon->ly;y++)
			if(clusters_get_value(clusters,x,y)==target)
				ising2d_set_spin(epsilon,x,y,value);
}

/*
	Clusters on a bilayer!
*/

struct bclusters_t *bclusters_init(int x,int y)
{
	struct bclusters_t *ret;
	
	if(!(ret=malloc(sizeof(struct bclusters_t))))
		return NULL;
	
	ret->vals[LOWER_LAYER]=malloc(sizeof(int)*x*y);
	ret->vals[UPPER_LAYER]=malloc(sizeof(int)*x*y);

	if((!ret->vals[UPPER_LAYER])||(!ret->vals[LOWER_LAYER]))
	{
		if(ret->vals[LOWER_LAYER])
			free(ret->vals[LOWER_LAYER]);

		if(ret->vals[UPPER_LAYER])
			free(ret->vals[UPPER_LAYER]);
	
		if(ret)
			free(ret);
		
		return NULL;
	}
	
	ret->lx=x;
	ret->ly=y;
	
	return ret;
}

void bclusters_fini(struct bclusters_t *bc)
{
	if(bc)
	{
		if(bc->vals[LOWER_LAYER])
			free(bc->vals[LOWER_LAYER]);

		if(bc->vals[UPPER_LAYER])
			free(bc->vals[UPPER_LAYER]);
	
		free(bc);
	}
}

int bclusters_get_value(struct bclusters_t *bclusters,int x,int y,int layer)
{
	return bclusters->vals[layer][MAKE_INDEX(bclusters,x,y)];
}

void bclusters_set_value(struct bclusters_t *bclusters,int x,int y,int layer,int value)
{
	bclusters->vals[layer][MAKE_INDEX(bclusters,x,y)]=value;
}

void ising2d_grow_bclusters(struct ibond2d_t *bonds[2],struct ivbond2d_t *ivbonds,struct bclusters_t *bclusters,int x,int y,int l,int nr_clusters)
{
	if(bclusters_get_value(bclusters,x,y,l)!=0)
		return;

	bclusters_set_value(bclusters,x,y,l,nr_clusters);

	if((ibond2d_get_value(bonds[l],x,y,DIR_X)==1)&&((x+1)<bclusters->lx))
		ising2d_grow_bclusters(bonds,ivbonds,bclusters,x+1,y,l,nr_clusters);

	if((ibond2d_get_value(bonds[l],x,y,DIR_Y)==1)&&((y+1)<bclusters->ly))
		ising2d_grow_bclusters(bonds,ivbonds,bclusters,x,y+1,l,nr_clusters);

	if(x>0)
		if(ibond2d_get_value(bonds[l],x-1,y,DIR_X)==1)
			ising2d_grow_bclusters(bonds,ivbonds,bclusters,x-1,y,l,nr_clusters);

	if(y>0)
		if(ibond2d_get_value(bonds[l],x,y-1,DIR_Y)==1)
			ising2d_grow_bclusters(bonds,ivbonds,bclusters,x,y-1,l,nr_clusters);

	if(ivbond2d_get_value(ivbonds,x,y)==1)
	{
		int otherlayer=l^1;
	
		ising2d_grow_bclusters(bonds,ivbonds,bclusters,x,y,otherlayer,nr_clusters);
	}
}

int ising2d_identify_bclusters(struct ibond2d_t *bonds[2],struct ivbond2d_t *ivbonds,struct bclusters_t *bclusters)
{
	int x,y,l;
	int nr_clusters=1;

	for(x=0;x<bclusters->lx;x++)
		for(y=0;y<bclusters->ly;y++)
			for(l=0;l<=1;l++)
				bclusters_set_value(bclusters,x,y,l,0);

	for(x=0;x<bclusters->lx;x++)
	{
		for(y=0;y<bclusters->ly;y++)
		{
			for(l=0;l<=1;l++)
			{
				if(bclusters_get_value(bclusters,x,y,l)==0)
				{
					ising2d_grow_bclusters(bonds,ivbonds,bclusters,x,y,l,nr_clusters);

					nr_clusters++;
				}
			}
		}
	}
	
	return nr_clusters;	
}

void bcluster_dimensions(struct ising2d_t *epsilon[2],struct bclusters_t *bclusters,int target,int *xdim,int *ydim)
{
	int maxx,maxy,minx,miny;
	int x,y,l;
	
	maxx=maxy=0;
	minx=epsilon[LOWER_LAYER]->lx;
	miny=epsilon[LOWER_LAYER]->ly;

	for(x=0;x<epsilon[LOWER_LAYER]->lx;x++)
	{
		for(y=0;y<epsilon[LOWER_LAYER]->ly;y++)
		{
			for(l=0;l<=1;l++)
			{
				if(bclusters->vals[l][MAKE_INDEX(bclusters,x,y)]==target)
				{
					if(x>maxx)
						maxx=x;

					if(y>maxy)
						maxy=y;

					if(x<minx)
						minx=x;

					if(y<miny)
						miny=y;
				}
			}
		}
	}

	*xdim=maxx-minx+1;
	*ydim=maxy-miny+1;

	if(*xdim<0)
		*xdim=0;
	
	if(*ydim<0)
		*ydim=0;
}

void reset_bcluster(struct ising2d_t *epsilon[2],struct bclusters_t *bclusters,int target,int value)
{
	int x,y,l;

	for(x=0;x<epsilon[LOWER_LAYER]->lx;x++)
		for(y=0;y<epsilon[LOWER_LAYER]->ly;y++)
			for(l=0;l<=1;l++)
				if(bclusters_get_value(bclusters,x,y,l)==target)
					ising2d_set_spin(epsilon[l],x,y,value);

}
