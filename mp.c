/*
	main_bilayer_phase_diagram(), massively parallel version

	FIXME: codice del tutto sperimentale, mai provato seriamente né debuggato.
*/

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <pthread.h>

#include "thr_pool/thr_pool.h"

#include "xy.h"
#include "pcc.h"
#include "stat.h"
#include "bilayer.h"

pthread_mutex_t write_mutex;

#define AVGSAMPLES_BILAYER_K_MP		(32)
#define NR_DIMENSIONS_BILAYER_K_MP	(6)

struct thread_args_t
{
	FILE *out;

	int millik;
	const int *dimensions;
	int nthreads;
};

void *phasediagram_mp_output_line(void *args)
{
	int d;
	struct thread_args_t *targs=args;
	
	for(d=0;d<NR_DIMENSIONS_BILAYER_K_MP;d++)
	{
		int c;
		double J,K;

		struct samples_t *tc;		

		J=1.0f;
		K=targs->millik*0.001f;

		tc=samples_init();

#pragma omp parallel num_threads(targs->nthreads)

		for(c=0;c<AVGSAMPLES_BILAYER_K_MP;c++)
		{
			double localtc;

			localtc=1.0f/pcc_bilayer(targs->dimensions[d],targs->dimensions[d],1.5f,J,J,K);

#pragma omp critical
			{
				samples_add_entry(tc,localtc);
			}
		}

		pthread_mutex_lock(&write_mutex);

		fprintf(targs->out,"%f %f ",K/samples_get_average(tc),J/samples_get_average(tc));
		fflush(targs->out);

		pthread_mutex_unlock(&write_mutex);

		samples_fini(tc);
	}

	pthread_mutex_lock(&write_mutex);

	fprintf(targs->out,"\n");
	fflush(targs->out);

	pthread_mutex_unlock(&write_mutex);

	if(targs)
		free(targs);

	return NULL;
}

int main_bilayer_phase_diagram_mp(int argc,char *argv[])
{
	int millik,deltamillik,d;
	int nworkers=16,nthreads=8;
	FILE *out;

	int dimensions[NR_DIMENSIONS_BILAYER_K_MP]={8,12,16,20,24,32};

	if(argc!=2)
	{
		printf("Usage: %s <logfile>\n",argv[0]);
		return 0;
	}

	if(!(out=fopen(argv[1],"w+")))
	{
		printf("Couldn't open %s for writing!\n",argv[1]);
		return 0;
	}

	thr_pool_t *thp=thr_pool_create(nworkers,nworkers,16,NULL);

	printf("Using %d worker threads in the pool, and %d inner parallel threads,\n",nworkers,nthreads);
	printf("for a total of %d concurrent threads.\n",nworkers*nthreads);

        pthread_mutex_init(&write_mutex,NULL);

	fprintf(out,"# Phase diagram: pairs of J/T and K/T are given on each line.\n");
	fprintf(out,"# Every pair is calculated for a different lattice size, they are:\n");
	fprintf(out,"# ");

	for(d=0;d<(NR_DIMENSIONS_BILAYER_K_MP-1);d++)
		fprintf(out,"2x%dx%d, ",dimensions[d],dimensions[d]);

	fprintf(out,"2x%dx%d\n",dimensions[NR_DIMENSIONS_BILAYER_K_MP-1],dimensions[NR_DIMENSIONS_BILAYER_K_MP-1]);
	fflush(out);

	deltamillik=125;
	for(millik=0;millik<=20000;millik+=deltamillik)
	{
		struct thread_args_t *ta=malloc(sizeof(struct thread_args_t));

		ta->out=out;
		ta->millik=millik;
		ta->dimensions=dimensions;
		ta->nthreads=nthreads;

		thr_pool_queue(thp,phasediagram_mp_output_line,ta);

		deltamillik=(millik<2000)?(125):(1000);
	}

	thr_pool_wait(thp);
	
	if(out)
		fclose(out);

	pthread_mutex_destroy(&write_mutex);

	thr_pool_destroy(thp);

	return 0;
}
