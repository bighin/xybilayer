#include <stdio.h>
#include <math.h>

#include "xy.h"
#include "pcc.h"
#include "stat.h"
#include "bilayer.h"

#define VERBOSE

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

#pragma omp parallel for

		for(c=0;c<AVGSAMPLES;c++)
		{
			double localtc;

			localtc=1.0f/pcc(dimensions[d],dimensions[d],1.5f,1.0f);

#pragma omp critical

			{
				samples_add_entry(tc,localtc);

#ifdef VERBOSE
				fprintf(stderr,".");
#endif
			}
		}

#ifdef VERBOSE
				fprintf(stderr,"\n");
#endif

		printf("%d %f +- %f\n",dimensions[d],samples_get_average(tc),sqrt(samples_get_variance(tc)));
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

int main(void)
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

#pragma omp parallel for

		for(c=0;c<AVGSAMPLES;c++)
		{
			double localtc;

			localtc=1.0f/pcc_bilayer(dimensions[d],dimensions[d],1.5f,1.0f,1.0f,0.0f);

#pragma omp critical

			{
				samples_add_entry(tc,localtc);

#ifdef VERBOSE
				fprintf(stderr,"%f\n",localtc);
#endif
			}
		}

#ifdef VERBOSE
				fprintf(stderr,"\n");
#endif

		printf("%d %f +- %f\n",dimensions[d],samples_get_average(tc),sqrt(samples_get_variance(tc)));
	}

	samples_fini(tc);

	return 0;
}
