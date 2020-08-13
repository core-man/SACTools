/*
 *  Stack SAC files to a single one.
 *
 *  Author: Jiayuan Yao @ NTU
 *
 *  Revisions:
 *      2015-01-24  Jiayuan Yao  Initial Coding
 *      2020-08-13  Jiayuan Yao  Better function naming and coding style
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <malloc.h>
#include <math.h>
#include "sacio.h"
#include "sac.h"

#define MAX_FNAME  256

typedef struct initial_vals
{
	char str[MAX_FNAME];
    double tt;
    double evla;
    double evlo;
    double stla;
    double stlo;
    double cela;
    double celo;
	double Gc0;
} IVAL;

void usage(void);
int get_file_line(char *fname);
int read_file(char *fname, int line_num, IVAL *IV);
int read_sac_data(IVAL *IV, int len, float **data, int *npts, float *b, float *depmin, float *depmax, float *delta, float *tt, int tmark);
int sac_norm(IVAL *IV, int len, float **data, int *npts, float *b, float *delta, float ts, float tw_norm);
int sac_stack(int len, float **data, int *npts, float *b, float *depmin, float *depmax, float *delta, IVAL *IV, int num, float *ccf, double slowness, float ts_stack);
int GCinit(double lat1, double lon1,double lat2, double lon2, double *x1, double *yp1, double *z1, double *x2, double *y2, double *z2, double *GCarc);
float *hilbert(float *in, int npts, float *out);

void usage() {
    fprintf(stderr, "Stack SAC files in a specified time window.        \n");
    fprintf(stderr, "                                                   \n");
    fprintf(stderr, "Usage:                                             \n");
    fprintf(stderr, "  sacstack [-tdatalist] [-wtime-window] [-ewf-dir] \n");
    fprintf(stderr, "           [-eenv-dir] [-llist-dir] [-btest-dir]   \n");
    fprintf(stderr, "           [-nnorm] [-mtmark]                      \n");
    fprintf(stderr, "                                                   \n");
    fprintf(stderr, "Options:                                           \n");
    fprintf(stderr, "  -t: data list                                    \n");
    fprintf(stderr, "  -w: time window in second                        \n");
    fprintf(stderr, "  -o: output waveform directory                    \n");
    fprintf(stderr, "  -e: output envelope directory                    \n");
    fprintf(stderr, "  -l: list directory                               \n");
    fprintf(stderr, "  -b: hilbert directory                            \n");
    fprintf(stderr, "  -n: normalization (1: yes; 0: no)                \n");
    fprintf(stderr, "  -m: tmark                                        \n");
    fprintf(stderr, "  -h: show usage                                   \n");
}

int main(int argc, char *argv[])
{
    int i,c;

    SACHEAD hdr;                        /* SAC headers */
	int tmark, npts_stack;
	int *npts;
	float *b, *delta, *tt, *depmin, *depmax;

	float ts_stack=-5.0, tw_stack; 	    /* stacking time window */
	float ts_norm=-10.0, tw_norm=20.0;   /* normalization time window */
	double slowness;                    /* slowness */

    int norm;                           /* normalization flag */
	float **data, **data_env;           /* waveforms and envelopes */
	float *ccf, *env;                   /* stacked waveform and envelope */

    int len=0;                                  /* file number */
	char datalist[MAX_FNAME], fname[MAX_FNAME]; /* data list and output file name */
    char wf_dir[MAX_FNAME], env_dir[MAX_FNAME]; /* output stacking direcotry */
	IVAL *IV;                                   /* input vars */

	//float delay;                              /* vars used to test codes */
    //FILE *fp;
    char list_dir[MAX_FNAME];

    //char fname_test[MAX_FNAME];       /* vars used to test hilbert tranform */
	//int n_scratch;
    //float *scratch;
    char test_dir[MAX_FNAME];
	float *hilbert, *env_test;

	double lat1,lat2,lon1,lon2,lat,lon;     /* vars used to calculate gcarc */
	double x1,yp1,z1,x2,y2,z2,GCarc,GCarc1;


    while ((c = getopt(argc, argv, "t:w:o:e:lbn:m:h")) != -1) {
        switch(c) {
			case 't':
	    	    sscanf(optarg, "%s", datalist);
			    break;
			case 'w':
		        sscanf(optarg, "%f", &tw_stack);
			    break;
			case 'o':
	    	    sscanf(optarg, "%s", wf_dir);
			    break;
			case 'e':
	    	    sscanf(optarg, "%s", env_dir);
			    break;
			case 'l':
	    	    sscanf(optarg, "%s", list_dir);
			    break;
			case 'b':
	    	    sscanf(optarg, "%s", test_dir);
			    break;
			case 'n':
	    	    sscanf(optarg, "%d", &norm);
	    	    break;
			case 'm':
	    	    sscanf(optarg, "%d", &tmark);
	    	    break;
            case 'h':
                usage();
                return -1;
        	default:
                usage();
                return -1;
        }
    }

    if (argc-optind != 0) {
        usage();
        return -1;
    }


    /*********** Read Data Part ***********/
    /* get file line number */
    len = get_file_line(datalist);

    /* set memory for struct vars */
    if ((IV = malloc(sizeof(IVAL)*len)) == NULL) {
        fprintf(stderr, "malloc memory error for IV\n");
        return -1;
    }

    /* read data list into strcut vars */
    len = read_file(datalist, len, IV);

    /* set memory for SAC headers and data */
	if ((npts = (int*)malloc(sizeof(int)*len)) == NULL) {
        fprintf(stderr, "Error in allocating memory for npts\n");
        return -1;
    }
	if ((b = (float*)malloc(sizeof(float)*len)) == NULL) {
        fprintf(stderr, "Error in allocating memory for b\n");
        return -1;
    }
	if ((depmin = (float*)malloc(sizeof(float)*len)) == NULL) {
        fprintf(stderr, "Error in allocating memory for depmin\n");
        return -1;
    }
	if ((depmax = (float*)malloc(sizeof(float)*len)) == NULL) {
        fprintf(stderr, "Error in allocating memory for depmax\n");
        return -1;
    }
	if ((delta = (float*)malloc(sizeof(float)*len)) == NULL) {
        fprintf(stderr, "Error in allocating memory for delta\n");
        return -1;
    }
	if ((tt = (float*)malloc(sizeof(float)*len)) == NULL) {     /* seems not used at all */
        fprintf(stderr, "Error in allocating memory for tt\n");
        return -1;
    }

	if ((data = (float**)malloc(sizeof(float*)*len)) == NULL) {
        fprintf(stderr, "Error in allocating memory for data\n");
        return -1;
    }

    /* read sac headers and data */
	read_sac_data(IV, len, data, npts, b, depmin, depmax, delta, tt, tmark);

    /* set memory and calculate data envelope */
	if ((data_env = (float**)malloc(sizeof(float*)*len)) == NULL) {
        fprintf(stderr, "Error in allocating memory for envelope\n");
        return -1;
    }
	for (i=0; i<len; i++) {
 		if ((data_env[i] = malloc(npts[i]*sizeof(*data_env[i]))) == NULL) {
    		fprintf(stderr, "Error in allocating memory for envelope %d\n", i);
        	return -1;
    	}
		envelope(npts[i], data[i], data_env[i]);
	}


    /*********** Data Processing Part ***********/
    /* normalize data and envelope with max amplitude
     * from tt + ts_norm with a time window of tw_norm sec */
    if (norm == 1) {
    	sac_norm(IV, len, data,     npts, b, delta, ts_norm, tw_norm);
	    sac_norm(IV, len, data_env, npts, b, delta, ts_norm, tw_norm);
    }

    /* time window used in stacking waveforms */
	npts_stack = tw_stack / delta[0];

    /* set memory for stacked data and envelope */
	if ((ccf = malloc(npts_stack*sizeof(*ccf))) == NULL) {
        fprintf(stderr, "Error in allocating memory for ccf\n");
        fprintf(stderr, "%d %d %f %f\n", npts_stack, npts[0], tw_stack, delta[0]);
		return -1;
	}
	if ((env = malloc(npts_stack*sizeof(*env))) == NULL) {
    	fprintf(stderr, "Error in allocating memory for envelope\n");
        return -1;
	}

	/* delta=delta[0], npts=npts_stack, b=IV[0].tt+ts_stack */
    hdr = new_sac_head(delta[0], npts_stack, IV[0].tt+ts_stack);
	//hdr.t0 = IV[0].tt;

    /* set memory for hilbert transfrom */
	if ((hilbert = malloc(npts_stack*sizeof(*hilbert))) == NULL) {
    	fprintf(stderr, "Error in allocating memory for hilbert\n");
        return -1;
	}
	if ((env_test = malloc(npts_stack*sizeof(*env_test))) == NULL) {
    	fprintf(stderr, "Error in allocating memory for envelope\n");
        return -1;
	}
    /*
	n_scratch = 5*1024;
  	scratch = (float *)calloc(n_scratch, sizeof(float));
   	if (scratch == NULL) {
        fprintf(stderr, "Error in allocating memory for hilbert transform\n");
        fprintf(stderr, "npts %d\n", n_scratch);
        return -1;
    }
    */

    /* calculate great-cricle distances */
	for (i=0; i<len; i++) {
	    lat = IV[i].cela; lat1 = IV[i].evla; lat2 = IV[i].stla;
        lon = IV[i].celo; lon1 = IV[i].evlo; lon2 = IV[i].stlo;
        GCinit(lat1,lon1,lat2,lon2,&x1,&yp1,&z1,&x2,&y2,&z2,&GCarc);
        GCinit(lat1,lon1,lat,lon,&x1,&yp1,&z1,&x2,&y2,&z2,&GCarc1);
        IV[i].Gc0 = GCarc - GCarc1;
    }


    /*********** Data Stacking Part ***********/
    /* slowness loop */
    //for (slowness = 3; slowness >= -3; slowness -= 0.1) {
	for (slowness = 0; slowness >= 0; slowness -= 0.1) {
        /* stack waveform wiggles */
		for (i=0; i<npts_stack; i++)
			ccf[i] = 0.0;
		sac_stack(len, data, npts, b, depmin, depmax, delta, IV, npts_stack, ccf, slowness, ts_stack);
    	snprintf(fname, MAX_FNAME, "%s/%.4f", wf_dir, slowness);
		write_sac(fname, hdr, ccf);

        /* stack envelope */
		for (i=0; i<npts_stack; i++)
			ccf[i] = 0.0;
		sac_stack(len, data_env, npts, b, depmin, depmax, delta, IV, npts_stack, ccf, slowness, ts_stack);
		snprintf(fname, MAX_FNAME, "%s/%.4f", env_dir, slowness);
		write_sac(fname, hdr, ccf);

        /* calcualte envelope for stacked waveform
        envelope(npts_stack, ccf, env);
		snprintf(fname, MAX_FNAME, "%s/%.4f", env_dir, slowness);
		write_sac(fname, hdr, env);
        */

        /*** test if delayed waveforms are aligned ***/
        /*
        snprintf(fname_test, MAX_FNAME, "%s/%.4f", list_dir, slowness);
		if ((fp = fopen(fname_test, "w")) == NULL) {
            fprintf(stderr, "open %s failed\n", fname_test);
			exit(-1);
        }
        fprintf(fp, "%s %f\n", fname, IV[0].tt);
		for (i=0; i<len; i++) {
			delay = IV[i].Gc0*slowness + IV[i].tt;
        	fprintf(fp, "%s %f\n", IV[i].str, delay);
		}
     	fclose(fp);
        */

        /*** test hilbert transfer ***/
        /* initialization for hilbert tranfrom */
        /*
		for (i=0; i<npts_stack; i++) {
			hilbert[i]  = 0.0;
			env_test[i] = 0.0;
		}

        //firtrn("HILBERT", ccf, npts_stack, scratch, hilbert);
        hilbert = hilbert(ccf, npts_stack);

		for (i=0; i<npts_stack; i++) {
			env_test[i] = sqrt(ccf[i]*ccf[i] + hilbert[i]*hilbert[i]);
		}
		snprintf(fname_test, MAX_FNAME, "%s/%.4f", test_dir, slowness);
        //write_sac(fname_test, hdr, hilbert);
		write_sac(fname_test, hdr, env_test);
        */
	} /* slowness loop ends*/


    /* free memory */
	free(IV);

	free(npts);
    free(b);
	free(depmin);
    free(depmax);
	free(delta);
    free(tt);

	for (i=0; i<len; i++) {
         free(data[i]);
         free(data_env[i]);
	}
	free(data);
    free(data_env);

    free(ccf);
	free(env);

    //free(scratch);
    free(hilbert);
    free(env_test);


    return(1);
}




/* **************************************************** */
/* **************************************************** */

/*
 * get file line number
 */
int get_file_line(char *fname) {
    int line_num=0;
    char c;
    FILE *fp;

    while ((fp = fopen(fname, "r")) == NULL) {
        fprintf(stderr, "Can not open file %s in read_file\n", fname);
        exit(-1);
    }

    for (c = getc(fp); c != EOF; c = getc(fp))
        if (c == '\n') line_num++;
    //fprintf(stdout, "line number: %d\n", line_num);

    fclose(fp);

    return line_num;
}



/*
 * read file into variables
 */
int read_file(char *fname, int line_num, IVAL *IV) {
    int i;
    FILE *fp;

    while ((fp = fopen(fname, "r")) == NULL) {
        fprintf(stderr, "Can not open file %s in read_file\n", fname);
        exit(-1);
    }

    for (i=0; i<line_num; i++) {
	    fscanf(fp, "%s %lf %lf %lf %lf %lf %lf %lf\n", IV[i].str,
               &IV[i].stlo,&IV[i].stla,
               &IV[i].evlo,&IV[i].evla,
               &IV[i].celo,&IV[i].cela,
               &IV[i].tt);
    }

    fclose(fp);

    if (i != line_num) {
        fprintf(stderr, "Error in reading file %s\n", fname);
        exit(-1);
    }

    return i;
}



/*
 * read SAC headers and data
 */
int read_sac_data(IVAL *IV,
                  int len,
                  float **data,
                  int *npts,
                  float *b,
                  float *depmin,
                  float *depmax,
                  float *delta,
                  float *tt,
                  int tmark)
{
	int i;
	SACHEAD hdr;

	for (i=0; i<len; i++) {
        data[i] = read_sac(IV[i].str, &hdr);

		npts[i]   = hdr.npts;
		b[i]      = hdr.b;
		depmin[i] = hdr.depmin;
		depmax[i] = hdr.depmax;
		delta[i]  = hdr.delta;
        tt[i]     = *((float *) &hdr + TMARK + tmark);
        //fprintf(stdout, "npts: %d\n", npts[i]);
        //fprintf(stdout, "tt: %f\n", tt[i]);
    }

    return(1);
}


/*
 * normalize data
 */
int sac_norm(IVAL *IV, int len, float **data, int *npts, float *b,
			 float *delta, float ts, float tw_norm)
{
	int i,j,k;
    int ns,nwin;
	float ts_rel;
	double sum;

    /* time window for normalization */
    //tw_norm = 15;

	for (i=0; i<len; i++) {
		ts_rel = (IV[i].tt + ts) - b[i];
		if (ts_rel <= 0) ts_rel = 0;

		ns   = ts_rel / delta[i];
		nwin = tw_norm / delta[i];

        /* find maximum */
		sum = 0.0;
		for (j=ns,k=0; j<npts[i] && k<nwin; j++,k++)
			sum = (fabs(data[i][j]) > sum) ? fabs(data[i][j]) : sum;

		/* normalize data */
		for (j = 0; j < npts[i]; j++)
			data[i][j] /= sum;
	}

	return 1;
}


/*
 * stack SAC waveforms
 */
int sac_stack(int len, float **data, int *npts, float *b, float *depmin,
              float *depmax, float *delta, IVAL *IV, int num, float *ccf,
              double slowness, float ts_stack)
{
	int i,j,k;
    //int nd;
	int ns=0, sta_num=0;
	float ts_rel;
    //float nth=4;
	double delay;

    for (i=0; i<len; i++) {
		if (isnan(depmin[i]) || isinf(depmax[i])) continue;
		if (delta[i] != delta[0]) continue; /* non-equale sampling rates */

        /* use tt as reference time
         * skip if begin time is too early. Possible bias!!! TODO */
		delay  = IV[i].Gc0 * slowness;
		ts_rel = (IV[i].tt + ts_stack) + delay - b[i];
		if (ts_rel <= 0) continue;

		ns = ts_rel / delta[i];
		for (j=0, k=ns; j<num && k<npts[i]; j++, k++) {
            /* linear stack */
			ccf[j] = ccf[j] + data[i][k];

			/* nth root stack */
			/*if (data[i][k] >= 0) {
				ccf[j] = ccf[j] + pow(fabs(data[i][k]), 1/nth);
			} else {
				ccf[j] = ccf[j] - pow(fabs(data[i][k]), 1/nth);
			}*/
		}

		sta_num++;
	}

	if (sta_num < 1) return(-1);    /* no data used in stacking */

    /* calculate stacked waveform */
	for (j = 0; j < num; j++) {
		ccf[j] = ccf[j] / sta_num;
        /* nth root stack */
        /*if (ccf[j] >= 0) {
			ccf[j] = pow(fabs(ccf[j]), nth );
		} else {
			ccf[j] = 0 - pow(fabs(ccf[j]), nth );
		}*/
	}

	return sta_num;
}



/*
 * calculate great-circle distance
 */
int GCinit(double lat1, double lon1,             /*IN: (lat,lon) in degrees */
           double lat2, double lon2,
           double *x1, double *yp1, double *z1,  /* OUT: xyz coordinates */
           double *x2, double *y2, double *z2,   /* and  distance in radius */
           double *GCarc)
{
	double the1,phe1,the2,phe2;
    double R2D = 57.2957795130823208768;
    double D2R = 0.017453292519943295769237;

	the1=(90.0-lat1)*D2R;		/* convert to radius */
  	phe1=lon1*D2R;
  	the2=(90.0-lat2)*D2R;
  	phe2=lon2*D2R;

  	*x1=sin(the1)*cos(phe1);
  	*yp1=sin(the1)*sin(phe1);
  	*z1=cos(the1);
  	*x2=sin(the2)*cos(phe2);
  	*y2=sin(the2)*sin(phe2);
  	*z2=cos(the2);
  	*GCarc=acos((*x1)*(*x2)+(*yp1)*(*y2)+(*z1)*(*z2));
    /*
  	if (fabs(*GCarc-M_PI) <= 1.e-16) {
        fprintf(stderr, " The great circle is not determined!\n");
        return(-1);
  	}
  	if (*GCarc <= 1.e-16) {
        fprintf(stderr,"Two same points. Program exits!\n");
        return(-1);
  	}*/
    *GCarc *= R2D;

    return(1);
}



/*
 * call firtrn to calcualte hilbert tranfrom
 */
float *hilbert(float *in, int npts, float *out) {
    float *buffer;
	//int n_scratch = 5*1024;
	int n_scratch = 50000;

    buffer = (float *)malloc(sizeof(float)*n_scratch);
    firtrn("HILBERT", in, npts, buffer, out);

    free(buffer);
    return out;
}


