/*
	This struct helps in collecting and (very simply) analyzing Monte Carlo samples
*/

#include <stdlib.h>
#include <math.h>

#include "stat.h"

struct samples_t *samples_init(void)
{
	struct samples_t *ret;
	
	if(!(ret=malloc(sizeof(struct samples_t))))
		return NULL;

	ret->nAlloced=1024;
	ret->next=0;
	ret->data=malloc(sizeof(double)*ret->nAlloced);
	
	if(!(ret->data))
	{
		if(ret)
			free(ret);
		
		return NULL;
	}
	
	return ret;
}

void samples_fini(struct samples_t *smpls)
{
	if(smpls)
	{
		if(smpls->data)
			free(smpls->data);
	
		free(smpls);
	}
}

void samples_add_entry(struct samples_t *smpls,double x)
{
	if(smpls->next==smpls->nAlloced)
	{
		smpls->nAlloced+=1024;
		smpls->data=realloc(smpls->data,sizeof(double)*smpls->nAlloced);
	}

	smpls->data[smpls->next]=x;
	smpls->next++;
}

double samples_get_average(struct samples_t *smpls)
{
	int c;
	double total=0.0f;
		
	for(c=0;c<smpls->next;c++)
		total+=smpls->data[c];
	
	return total/((double)(smpls->next));
}

double samples_get_variance(struct samples_t *smpls)
{
	int c;
	double total,average,variance,n;
	
	total=variance=0.0f;
	n=smpls->next;

	for(c=0;c<smpls->next;c++)
		total+=smpls->data[c];

	average=total/n;

	for(c=0;c<smpls->next;c++)
		variance+=pow(smpls->data[c]-average,2.0f);

	variance*=(1.0f/(n-1));

	return variance;
}

/*
	Linear regression
*/

struct linreg_ctx_t *init_linreg_ctx(void)
{
	struct linreg_ctx_t *lct=malloc(sizeof(struct linreg_ctx_t));
	
	if(!lct)
		return NULL;
	
	lct->n=0;
	lct->sx=lct->sy=lct->sxx=lct->syy=lct->sxy=0.0f;

	return lct;
}

void fini_linreg_ctx(struct linreg_ctx_t *lct)
{
	if(lct)
		free(lct);
}

void linreg_add_entry(struct linreg_ctx_t *lct,double x,double y)
{
	lct->sx+=x;
	lct->sy+=y;
	lct->sxx+=x*x;
	lct->syy+=y*y;
	lct->sxy+=x*y;
	lct->n++;
}

double slope(struct linreg_ctx_t *lct)
{
	double a,b,xmean,ymean;
	
	xmean=lct->sx/((double)(lct->n));
	ymean=lct->sy/((double)(lct->n));

	a=lct->sxy-ymean*lct->sx-xmean*lct->sy+lct->n*xmean*ymean;
	b=lct->sxx+lct->n*xmean*xmean-2.0f*xmean*lct->sx;

	return a/b;
}

double intercept(struct linreg_ctx_t *lct)
{
	double xmean,ymean;
	
	xmean=lct->sx/((double)(lct->n));
	ymean=lct->sy/((double)(lct->n));	

	return ymean-slope(lct)*xmean;
}

double slope_error(struct linreg_ctx_t *lct)
{
	double a,b,ret=0.0f;
	
	a=intercept(lct);
	b=slope(lct);

	ret+=lct->syy+lct->n*a*a+b*b*lct->sxx;
	ret+=2.0f*(-lct->sy*a+lct->sxy*b-a*b*lct->sx);

	return sqrt(ret/((double)(lct->n)));
}
