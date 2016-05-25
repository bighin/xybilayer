#include <stdio.h>
#include <math.h>
#include <assert.h>

#define PARALLEL

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

#define AVGSAMPLES_BILAYER_K	(32)
#define NR_DIMENSIONS_BILAYER_K	(6)

	int dimensions[NR_DIMENSIONS_BILAYER_K]={8,12,16,20,24,32};

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

	fprintf(out,"# Phase diagram: pairs of J/T and K/T are given on each line.\n");
	fprintf(out,"# Every pair is calculated for a different lattice size, they are:\n");
	fprintf(out,"# ");

	for(d=0;d<(NR_DIMENSIONS_BILAYER_K-1);d++)
		fprintf(out,"2x%dx%d, ",dimensions[d],dimensions[d]);

	fprintf(out,"2x%dx%d\n",dimensions[NR_DIMENSIONS_BILAYER_K-1],dimensions[NR_DIMENSIONS_BILAYER_K-1]);
	fflush(out);

	deltamillik=125;
	for(millik=0;millik<=20000;millik+=deltamillik)
	{
		for(d=0;d<NR_DIMENSIONS_BILAYER_K;d++)
		{

			progressbar *progress;
			struct samples_t *tc;

			char description[128];
		
			double J,K;

			J=1.0f;
			K=millik*0.001f;

			deltamillik=(millik<2000)?(125):(1000);
		
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

int do_correlators(FILE *out,int x,int y,int runs,double beta,double Jup,double Jdown,double K)
{
	double *wmeans,*variances;
	int c,maxk,total_channels;

	progressbar *progress;
	char description[1024];

	assert(x>0);
	assert(y>0);

	maxk=MIN(x,y)/2;
	total_channels=get_total_channels(maxk);

	/*
		Inizializziamo l'array dove salveremo le medie pesate dei risultati e le
		relative deviazioni standard.
	*/

	wmeans=malloc(sizeof(double)*get_total_channels(maxk));
	variances=malloc(sizeof(double)*get_total_channels(maxk));

	for(c=0;c<total_channels;c++)
	{
		wmeans[c]=0.0f;
		variances[c]=0.0f;
	}

	/*
		Inizializziamo la progressbar, scriviamo l'header sul file e siamo pronti!
	*/

	snprintf(description,1024,"(2x%dx%d lattice, beta=%f, J=%f, K=%f)",x,y,beta,Jup,K);
	description[1023]='\0';

	progress=progressbar_new(description,runs);

	fprintf(out,"# 2x%dx%d lattice, runs=%d, beta=%f, J=%f, K=%f\n",x,y,runs,beta,Jup,K);
	fprintf(out,"# k z(k) c_up(k) c_lo(k) sigma(z(k)) sigma(c_up(k)) sigma(c_lo(k))\n");

#ifdef PARALLEL
#pragma omp parallel for
#endif

	for(c=0;c<runs;c++)
	{
		struct sampling_ctx_t *sctx;
		
		sctx=sw_bilayer(x,y,beta,Jup,Jdown,K,maxk);

#ifdef PARALLEL
#pragma omp critical
#endif

		{
			double *average=malloc(sizeof(double)*get_total_channels(maxk));
			double *variance=malloc(sizeof(double)*get_total_channels(maxk));
			int d;
			
			sampling_ctx_to_tuple(sctx,average,variance);

			for(d=0;d<total_channels;d++)
			{
				double weight=1.0f/variance[d];
			
				wmeans[d]+=average[d]*weight;
				variances[d]+=weight;
			}

			progressbar_inc(progress);

			if(average)
				free(average);

			if(variance)
				free(variance);
		}
		
		sampling_ctx_fini(sctx);
	}

	/*
		Completiamo il calcolo della media pesata e stampiamo i risultati
	*/

	for(c=0;c<total_channels;c++)
	{
		wmeans[c]/=variances[c];
		variances[c]=1.0f/variances[c];
	}

	for(c=0;c<maxk;c++)
	{
		fprintf(out,"%d ",c);
		fprintf(out,"%f ",wmeans[get_channel_nr(maxk,"zk real",c)]);
		fprintf(out,"%f ",wmeans[get_channel_nr(maxk,"ck upper real",c)]);
		fprintf(out,"%f ",wmeans[get_channel_nr(maxk,"ck lower real",c)]);
		fprintf(out,"%f ",sqrt(variances[get_channel_nr(maxk,"zk real",c)]));
		fprintf(out,"%f ",sqrt(variances[get_channel_nr(maxk,"ck upper real",c)]));
		fprintf(out,"%f\n",sqrt(variances[get_channel_nr(maxk,"ck lower real",c)]));
	}

	progressbar_finish(progress);

	if(wmeans)
		free(wmeans);

	if(variances)
		free(variances);
	
	return 0;
}

int main_bilayer_correlators(int argc,char *argv[])
{
	int x,y,runs;
	double beta,Jup,Jdown,K;
	FILE *out;

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
		Let's go!
	*/
	
	do_correlators(out,x,y,runs,beta,Jup,Jdown,K);

	if(out)
		fclose(out);

	return 0;
}

#undef DRYRUN

int main_bilayer_correlators_phase_diagram(int argc,char *argv[])
{
	int centiJ,centiK,x,y,runs;
	double beta,J,K;
	
	char *basename="correlators";
	
	x=y=32;
	runs=24;
	beta=1.0f;

	centiJ=35;
	centiK=1000;

	for(centiK=1000;centiK>=50;centiK-=50)
	{
		FILE *out;
		char fname[1024];
		
		snprintf(fname,1024,"%s.%d.%d.dat",basename,centiJ,centiK);
		fname[1023]='\0';

#ifdef DRYRUN
		printf("%s\n",fname);
		continue;
#endif

		J=0.01f*centiJ;
		K=0.01f*centiK;

		if(!(out=fopen(fname,"w+")))
		{
			printf("Couldn't open %s for writing!\n",argv[1]);
			return 0;
		}

		do_correlators(out,x,y,runs,beta,J,J,K);

		if(out)
			fclose(out);
	}

	centiJ=60;
	centiK=50;

	for(centiJ=60;centiJ<=260;centiJ+=25)
	{
		FILE *out;
		char fname[1024];
				
		snprintf(fname,1024,"%s.%d.%d.dat",basename,centiJ,centiK);
		fname[1023]='\0';

		J=0.01f*centiJ;
		K=0.01f*centiK;

#ifdef DRYRUN
		printf("%s\n",fname);
		continue;
#endif

		if(!(out=fopen(fname,"w+")))
		{
			printf("Couldn't open %s for writing!\n",argv[1]);
			return 0;
		}

		do_correlators(out,x,y,runs,beta,J,J,K);

		if(out)
			fclose(out);
	}

	centiJ=260;
	centiK=100;

	for(centiK=100;centiK<=1000;centiK+=50)
	{
		FILE *out;
		char fname[1024];
		
		snprintf(fname,1024,"%s.%d.%d.dat",basename,centiJ,centiK);
		fname[1023]='\0';

#ifdef DRYRUN
		printf("%s\n",fname);
		continue;
#endif

		J=0.01f*centiJ;
		K=0.01f*centiK;
		
		if(!(out=fopen(fname,"w+")))
		{
			printf("Couldn't open %s for writing!\n",argv[1]);
			return 0;
		}

		do_correlators(out,x,y,runs,beta,J,J,K);

		if(out)
			fclose(out);
	}

	centiJ=235;
	centiK=1000;

	for(centiJ=235;centiJ>=60;centiJ-=25)
	{
		FILE *out;
		char fname[1024];
		
		snprintf(fname,1024,"%s.%d.%d.dat",basename,centiJ,centiK);
		fname[1023]='\0';

#ifdef DRYRUN
		printf("%s\n",fname);
		continue;
#endif

		J=0.01f*centiJ;
		K=0.01f*centiK;
		
		if(!(out=fopen(fname,"w+")))
		{
			printf("Couldn't open %s for writing!\n",argv[1]);
			return 0;
		}

		do_correlators(out,x,y,runs,beta,J,J,K);

		if(out)
			fclose(out);
	}

	return 0;
}

int main_bilayer_correlators_phase_diagram(int argc,char *argv[])
{
	int centiJ,centiK,x,y,runs;
	double beta,J,K;

	char *basename="cx48";

	x=y=48;
	runs=24;
	beta=1.0f;

	for(centiK=25;centiK<=1000;centiK+=25)
	{
		for(centiJ=60;centiJ<=260;centiJ+=25)
		{
			FILE *out;
			char fname[1024];

			snprintf(fname,1024,"%s.%d.%d.dat",basename,centiJ,centiK);
			fname[1023]='\0';

#ifdef DRYRUN
			printf("%s\n",fname);
			continue;
#endif
			J=0.01f*centiJ;
			K=0.01f*centiK;

			if(!(out=fopen(fname,"w+")))
			{
				printf("Couldn't open %s for writing!\n",argv[1]);
				return 0;
			}

			do_correlators(out,x,y,runs,beta,J,J,K);

			if(out)
				fclose(out);
		}
	}

	return 0;
}

int main(int argc,char *argv[])
{
	init_prng();

	//main_bilayer_phase_diagram_mp(argc,argv);
	return main_bilayer_correlators(argc,argv);
	//return main_bilayer_phase_diagram(argc,argv);
	//return main_xy(argc,argv);
}
