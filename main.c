#include <stdio.h>
#include <math.h>
#include <assert.h>

#undef PARALLEL

#include "xy.h"
#include "pcc.h"
#include "stat.h"
#include "bilayer.h"
#include "libprogressbar/progressbar.h"

int main_xy(int argc,char *argv[])
{
	int c,d;
	FILE *out;

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

	init_prng();

#define NR_DIMENSIONS	(5)
#define AVGSAMPLES	(100)

	for(d=0;d<NR_DIMENSIONS;d++)
	{
		int dimensions[NR_DIMENSIONS]={8,16,24,32,48};
		struct samples_t *tc;
		progressbar *progress;

		char description[128];

		snprintf(description,128,"D=%d",dimensions[d]);

		tc=samples_init();
		progress=progressbar_new(description,AVGSAMPLES);

#ifdef PARALLEL
#pragma omp parallel for
#endif
		for(c=0;c<AVGSAMPLES;c++)
		{
			double localtc;

			localtc=1.0f/pcc(dimensions[d],dimensions[d],1.5f,1.0f);

#ifdef PARALLEL
#pragma omp critical
#endif
			{
				samples_add_entry(tc,localtc);
				progressbar_inc(progress);
			}
		}

		progressbar_finish(progress);

		fprintf(out,"%d %f ",dimensions[d],samples_get_average(tc));
		fflush(out);

		fprintf(stdout,"%d %f +- %f\n",dimensions[d],samples_get_average(tc),sqrt(samples_get_variance(tc)));
		fflush(stdout);

		samples_fini(tc);
	}

	if(out)
		fclose(out);

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

int main_bilayer_crosscheck(void)
{
	int c,d;

	init_prng();

#define NR_DIMENSIONS_BILAYER	(5)
#define AVGSAMPLES_BILAYER	(20)

	for(d=0;d<NR_DIMENSIONS_BILAYER;d++)
	{
		int dimensions[NR_DIMENSIONS_BILAYER]={8,16,24,32,48};
		struct samples_t *tc;
		progressbar *progress;

		char description[128];

		snprintf(description,128,"D=%d",dimensions[d]);

		tc=samples_init();
		progress=progressbar_new(description,AVGSAMPLES);

#ifdef PARALLEL
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

		samples_fini(tc);
	}

	return 0;
}

int main_bilayer_phase_diagram(int argc,char *argv[])
{
	int c,millik,deltamillik,d;
	FILE *out;

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

	init_prng();

#define AVGSAMPLES_BILAYER_K	(32)
#define NR_DIMENSIONS_BILAYER_K	(4)

	deltamillik=250;
	for(millik=0;millik<=20000;millik+=deltamillik)
	{
		for(d=0;d<NR_DIMENSIONS_BILAYER_K;d++)
		{
			int dimensions[NR_DIMENSIONS_BILAYER]={8,12,16,20};

			progressbar *progress;
			struct samples_t *tc;

			char description[128];
		
			double J,K;

			J=1.0f;
			K=millik*0.001f;

			deltamillik=(millik<2000)?(250):(1000);
		
			snprintf(description,128,"(2x%dx%d lattice, K=%f)",dimensions[d],dimensions[d],K);

			progress=progressbar_new(description,AVGSAMPLES_BILAYER);
			tc=samples_init();

#ifdef PARALLEL
#pragma omp parallel for
#endif

			for(c=0;c<AVGSAMPLES_BILAYER_K;c++)
			{
				double localtc;

				localtc=1.0f/pcc_bilayer(dimensions[d],dimensions[d],1.5f,J,J,K);

#ifdef PARALLEL
#pragma omp critical
#endif
				{
					samples_add_entry(tc,localtc);
					progressbar_inc(progress);
				}
			}

			progressbar_finish(progress);

			fprintf(stdout,"beta*K=%f beta*J=%f\n",K/samples_get_average(tc),J/samples_get_average(tc));
			fflush(stdout);

			fprintf(out,"%f %f ",K/samples_get_average(tc),J/samples_get_average(tc));
			fflush(out);

			samples_fini(tc);
		}

		fprintf(out,"\n");
		fflush(out);
	}
	
	if(out)
		fclose(out);

	return 0;
}

int main_bilayer_correlator(int argc,char *argv[])
{
	int x,y,k,maxk,c,runs;
	double beta,Jup,Jdown,K;

	double *zks,*localzks,*localszks;
	FILE *out;

	progressbar *progress;
	char description[1024];

	if(argc!=7)
	{
		printf("Usage: %s <logfile> <lattice-size> <runs> <beta> <J> <K>\n",argv[0]);
		return 0;
	}

	if(!(out=fopen(argv[1],"w+")))
	{
		printf("Couldn't open %s for writing!\n",argv[1]);
		return 0;
	}

	/*
		Inizializziamo i parametri della simulazione leggendo dalla riga di comando,
		e controlliamo che non siano totalmente assurdi.
	*/

	x=atoi(argv[2]);
	y=atoi(argv[2]);

	assert(x>0);
	assert(x<=1024);

	maxk=MIN(x,y)/2;
	runs=atoi(argv[3]);
	
	assert(runs>0);
	assert(runs<=1024*4);

	beta=atof(argv[4]);

	Jup=atof(argv[5]);
	Jdown=atof(argv[5]);
	K=atof(argv[6]);

	assert(Jup>0.0f);
	assert(K>0.0f);

	/*
		Inizializziamo gli array che conterranno il correlatore
	*/

	zks=malloc(sizeof(double)*maxk);
	localzks=malloc(sizeof(double)*maxk);
	localszks=malloc(sizeof(double)*maxk);

	for(k=0;k<maxk;k++)
		zks[k]=0.0f;

	/*
		Inizializziamo la progressbar e siamo pronti!
	*/

	snprintf(description,1024,"(2x%dx%d lattice, beta=%f, J=%f, K=%f)",x,y,beta,Jup,K);
	description[1023]='\0';

	progress=progressbar_new(description,AVGSAMPLES_BILAYER);

	for(c=0;c<runs;c++)
	{
		sw_bilayer(x,y,beta,Jup,Jdown,K,localzks,localszks,maxk);

		for(k=0;k<maxk;k++)
			zks[k]+=localzks[k];
		
		progressbar_inc(progress);
	}

	progressbar_finish(progress);

	for(k=0;k<maxk;k++)
		zks[k]/=runs;

	for(k=0;k<maxk;k++)
		printf("%d %f\n",k,zks[k]);

	if(zks)
		free(zks);

	if(localzks)
		free(localzks);

	if(localszks)
		free(localszks);

	if(out)
		fclose(out);

	return 0;
}

int main(int argc,char *argv[])
{
	return main_bilayer_correlator(argc,argv);

	//return main_bilayer_phasediagram(argc,argv);
	//return main_xy(argc,argv);
}
