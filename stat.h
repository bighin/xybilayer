#ifndef __STAT_H__
#define __STAT_H__

struct samples_t
{
	int nAlloced,next;

	double *data;
};

struct samples_t *samples_init(void);
void samples_fini(struct samples_t *smpls);
void samples_add_entry(struct samples_t *smpls,double x);
double samples_get_average(struct samples_t *smpls);
double samples_get_variance(struct samples_t *smpls);

struct linreg_ctx_t
{
	int n;
	double sx,sy,sxx,syy,sxy;
};

struct linreg_ctx_t *init_linreg_ctx(void);
void fini_linreg_ctx(struct linreg_ctx_t *lct);
void linreg_add_entry(struct linreg_ctx_t *lct,double x,double y);
double slope(struct linreg_ctx_t *lct);
double intercept(struct linreg_ctx_t *lct);
double slope_error(struct linreg_ctx_t *lct);

#endif
