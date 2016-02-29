#include <stdio.h>
#include <math.h>

#include "xy.h"
#include "pcc.h"
#include "stat.h"
#include "bilayer.h"
#include "libprogressbar/progressbar.h"

int main_xy(void)
{
	struct samples_t *tc;
	int c,d;

	init_prng();

	tc=samples_init();

#define NR_DIMENSIONS	(5)
#define AVGSAMPLES	(100)

	for(d=0;d<NR_DIMENSIONS;d++)
	{
		int dimensions[NR_DIMENSIONS]={8,16,24,32,48};
		progressbar *progress;
		char description[128];
		
		snprintf(description,128,"D=%d",dimensions[d]);

		progress=progressbar_new(description,AVGSAMPLES);

#pragma omp parallel for

		for(c=0;c<AVGSAMPLES;c++)
		{
			double localtc;

			localtc=1.0f/pcc(dimensions[d],dimensions[d],1.5f,1.0f);

#pragma omp critical

			{
				samples_add_entry(tc,localtc);
				progressbar_inc(progress);
			}
		}

		progressbar_finish(progress);

		fprintf(stdout,"%d %f +- %f\n",dimensions[d],samples_get_average(tc),sqrt(samples_get_variance(tc)));
		fflush(stdout);
	}

	samples_fini(tc);

	return 0;
}

int main_ising(void)
{
	int c;
	
	init_prng();

	for(c=0;c<10;c++)
		printf("%f\n",1.0f/pcc_ising(16,16,1.5f,1.0f));

	printf("\n");

	for(c=0;c<10;c++)
		printf("%f\n",1.0f/pcc_ising(24,24,1.5f,1.0f));

	printf("\n");

	for(c=0;c<10;c++)
		printf("%f\n",1.0f/pcc_ising(32,32,1.5f,1.0f));

	printf("\n");

	for(c=0;c<10;c++)
		printf("%f\n",1.0f/pcc_ising(48,48,1.5f,1.0f));

	return 0;
}

int main_bilayer(void)
{
	struct samples_t *tc;
	int c,d;

	init_prng();

	tc=samples_init();

#define NR_DIMENSIONS_BILAYER	(5)
#define AVGSAMPLES_BILAYER	(20)

	for(d=0;d<NR_DIMENSIONS_BILAYER;d++)
	{
		int dimensions[NR_DIMENSIONS_BILAYER]={8,16,24,32,48};
		progressbar *progress;
		char description[128];
		
		snprintf(description,128,"D=%d",dimensions[d]);

		progress=progressbar_new(description,AVGSAMPLES_BILAYER);

#define PARALLEL

#ifdef PARALLELA
#pragma omp parallel for
#endif

		for(c=0;c<AVGSAMPLES_BILAYER;c++)
		{
			double localtc;

			localtc=1.0f/pcc_bilayer(dimensions[d],dimensions[d],1.5f,1.0f,1.0f,0.0f);

#ifdef PARALLEL
#pragma omp critical
#endif
			{
				samples_add_entry(tc,localtc);
				progressbar_inc(progress);
			}
		}

		progressbar_finish(progress);

		fprintf(stdout,"%d %f +- %f\n",dimensions[d],samples_get_average(tc),sqrt(samples_get_variance(tc)));
		fflush(stdout);
	}

	samples_fini(tc);

	return 0;
}

int main(void)
{
	return main_ising();
}
