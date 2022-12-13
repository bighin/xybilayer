#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

#define PARALLEL

#include "xy.h"
#include "pcc.h"
#include "stat.h"
#include "bilayer.h"
#include "libprogressbar/progressbar.h"

int do_correlators(FILE *out,int x,int y,int runs,double beta,double Jup,double Jdown,double K)
{
	double *wmeans,*variances;
	int c,maxk,total_channels;

	assert(x>0);
	assert(y>0);

	maxk=MIN(x,y)/2;
	total_channels=get_total_channels(maxk);

	/*
		We initialize an array where the means and variances of sampled
		quantities will be saved.
	*/

	wmeans=malloc(sizeof(double)*get_total_channels(maxk));
	variances=malloc(sizeof(double)*get_total_channels(maxk));

	for(c=0;c<total_channels;c++)
	{
		wmeans[c]=0.0f;
		variances[c]=0.0f;
	}

	/*
		We write a nice header on the output file and we are good to go!
	*/

	fprintf(out,"# 2x%dx%d lattice, runs=%d, beta=%f, J=%f, K=%f\n",x,y,runs,beta,Jup,K);
	fprintf(out,"# k z(k) c_up(k) c_lo(k) sigma(z(k)) sigma(c_up(k)) sigma(c_lo(k))\n");

	/*
		Do all the Monte Carlo runs, saving and the results for calculating averages.
	*/

	for(c=0;c<runs;c++)
	{
		struct sampling_ctx_t *sctx;
		double *average=malloc(sizeof(double)*get_total_channels(maxk));
		double *variance=malloc(sizeof(double)*get_total_channels(maxk));
		int d;

		sctx=sw_bilayer(x,y,beta,Jup,Jdown,K,maxk);

		sampling_ctx_to_tuple(sctx,average,variance);

		for(d=0;d<total_channels;d++)
		{
			double weight=1.0f/variance[d];
			
			wmeans[d]+=average[d]*weight;
			variances[d]+=weight;
		}

		if(average)
			free(average);

		if(variance)
			free(variance);
		
		sampling_ctx_fini(sctx);
	}

	/*
		Finally we calculate a weighted average and print the results.
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

	if(wmeans)
		free(wmeans);

	if(variances)
		free(variances);
	
	return 0;
}

int main_bilayer_correlators_phase_diagram(int argc,char *argv[])
{
	int centiJ,centiK,x,y,runs;
	double beta,J,K;

	progressbar *progress;
	char description[1024];

	char *basename="cx20";

	x=y=20;
	runs=10;
	beta=1.0f;

	snprintf(description,1024,"(2x%dx%d lattice, beta=%.3f, J=%.3f, K=%.3f, prefix=%s, nthreads=%d)",x,y,beta,J,K,basename,omp_get_max_threads());
	description[1023]='\0';

	progress=progressbar_new(description,(250/5)*(200/5));

#pragma omp parallel for

	for(centiK=0;centiK<=250;centiK+=5)
	{

#pragma omp parallel for

		for(centiJ=60;centiJ<=260;centiJ+=5)
		{
			FILE *out;
			char fname[1024];

			snprintf(fname,1024,"%s.%d.%d.dat",basename,centiJ,centiK);
			fname[1023]='\0';

			J=0.01f*centiJ;
			K=0.01f*centiK;

			if(!(out=fopen(fname,"w+")))
			{
				printf("Couldn't open %s for writing!\n",argv[1]);
				continue;
			}

			do_correlators(out,x,y,runs,beta,J,J,K);

			if(out)
				fclose(out);

#pragma omp critical

			{
				progressbar_inc(progress);
			}
		}
	}

	progressbar_finish(progress);

	return 0;
}

int main(int argc,char *argv[])
{
	init_prng();

	return main_bilayer_correlators_phase_diagram(argc,argv);
}
