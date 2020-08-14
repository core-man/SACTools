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

#define MAX_FNAME  256

typedef struct initial_vals
{
	char str[MAX_FNAME];
    double tt;
} IVAL;

void usage(void);
int get_file_line(char *fname);
int read_file(char *fname, int line_num, IVAL *IV);
int read_sac_data(IVAL *IV, int len, float **data, int *npts, float *b, float *depmin, float *depmax, float *delta, float *tt, int tmark);
int sac_norm(IVAL *IV, int len, float **data, int *npts, float *b, float *delta, float ts, float tw_norm);
int sac_stack(int len, float **data, int *npts, float *b, float *depmin, float *depmax, float *delta, IVAL *IV, int num, float *data_stack, float ts_stack);

void usage() {
    fprintf(stderr, "Stack SAC files in a specified time window.        \n");
    fprintf(stderr, "                                                   \n");
    fprintf(stderr, "Usage:                                             \n");
    fprintf(stderr, "  sacstack [-Ddatalist] [-Ttmark/ts/tw] [-Ooutifle]\n");
    fprintf(stderr, "           [-Nnorm/ts/tw]                          \n");
    fprintf(stderr, "                                                   \n");
    fprintf(stderr, "Options:                                           \n");
    fprintf(stderr, "  -D: data list                                    \n");
    fprintf(stderr, "  -T: tmark/begin time (sec)/time window (sec)     \n");
    fprintf(stderr, "  -O: output file                                  \n");
    fprintf(stderr, "  -N: normalization(1: yes; 0: no)/begin time (sec)/time window (sec)\n");
    fprintf(stderr, "  -h: show usage                                   \n");
}

int main(int argc, char *argv[])
{
    int c;
    int error;
    int i;

    int len=0;                                    /* file number */
	char datalist[MAX_FNAME], outfile[MAX_FNAME]; /* data list and output file */
	IVAL *IV;                                     /* input vars */

    SACHEAD hdr;                        /* SAC headers */
	int tmark, npts_stack;
	int *npts;
	float *b, *delta, *tt, *depmin, *depmax;

	float ts_stack=-5.0, tw_stack=30.0; /* stacking time window */
    int norm;                           /* normalization flag */
	float ts_norm=-10.0, tw_norm=20.0;  /* normalization time window */
	float **data, *data_stack;          /* waveforms and stacked waveform */


    error = 0;
    while ((c = getopt(argc, argv, "D:T:O:N:h")) != -1) {
        switch(c) {
			case 'D':
	    	    sscanf(optarg, "%s", datalist);
			    break;
			case 'T':
		        if (sscanf(optarg, "%d/%f/%f", &tmark, &ts_stack, &tw_stack) != 3)
                    error++;
			    break;
			case 'O':
	    	    sscanf(optarg, "%s", outfile);
			    break;
			case 'N':
		        if (sscanf(optarg, "%d/%f/%f", &norm, &ts_norm, &tw_norm) != 3)
                    error++;
	    	    break;
            case 'h':
                usage();
                return -1;
        	default:
                usage();
                return -1;
        }
    }

    if (argc-optind != 0 || error) {
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
	if ((data = (float**)malloc(sizeof(float*)*len)) == NULL) {
        fprintf(stderr, "Error in allocating memory for data\n");
        return -1;
    }
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
    /* tt is not used now */
	if ((tt = (float*)malloc(sizeof(float)*len)) == NULL) {
        fprintf(stderr, "Error in allocating memory for tt\n");
        return -1;
    }

    /* read sac headers and data */
	read_sac_data(IV, len, data, npts, b, depmin, depmax, delta, tt, tmark);


    /*********** Data Processing Part ***********/
    /* normalize data and envelope with max amplitude
     * from tt + ts_norm with a time window of tw_norm sec */
    if (norm == 1) sac_norm(IV, len, data, npts, b, delta, ts_norm, tw_norm);


    /*********** Data Stacking Part ***********/
    /* time window used in stacking waveforms */
    /* using delta in the first file!!! */
	npts_stack = tw_stack / delta[0];

    /* set memory for stacked data and envelope */
	if ((data_stack = malloc(npts_stack*sizeof(*data_stack))) == NULL) {
        fprintf(stderr, "Error in allocating memory for data_stack\n");
		return -1;
	}

    /* stack waveforms */
	for (i=0; i<npts_stack; i++) data_stack[i] = 0.0;
	sac_stack(len, data, npts, b, depmin, depmax, delta, IV, npts_stack, data_stack, ts_stack);

    /* output stacked waveform */
	/* delta=delta[0]
     * npts=npts_stack
     * b=IV[0].tt+ts_stack */
    hdr = new_sac_head(delta[0], npts_stack, IV[0].tt+ts_stack);
	write_sac(outfile, hdr, data_stack);


    /*** free memory ***/
	free(IV);

	free(npts);   free(b);
	free(depmin); free(depmax);
	free(delta);  free(tt);

	for (i=0; i<len; i++) free(data[i]);
	free(data);
    free(data_stack);

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

    for (i=0; i<line_num; i++)
	    fscanf(fp, "%s %lf\n", IV[i].str, &IV[i].tt);

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
              float *depmax, float *delta, IVAL *IV, int num, float *data_stack,
              float ts_stack)
{
	int i,j,k;
    //int nd;
	int ns=0, sta_num=0;
	float ts_rel;
    //float nth=4;

    for (i=0; i<len; i++) {
		if (isnan(depmin[i]) || isinf(depmax[i])) continue;
		if (delta[i] != delta[0]) continue; /* non-equale sampling rates */

        /* use tt as reference time
         * skip if begin time is too early. Possible bias!!! TODO */
		ts_rel = (IV[i].tt + ts_stack) - b[i];
		if (ts_rel <= 0) continue;

		ns = ts_rel / delta[i];
		for (j=0, k=ns; j<num && k<npts[i]; j++, k++) {
            /* linear stack */
			data_stack[j] = data_stack[j] + data[i][k];

			/* nth root stack */
			/*if (data[i][k] >= 0) {
				data_stack[j] = data_stack[j] + pow(fabs(data[i][k]), 1/nth);
			} else {
				data_stack[j] = data_stack[j] - pow(fabs(data[i][k]), 1/nth);
			}*/
		}

		sta_num++;
	}

	if (sta_num < 1) return(-1);    /* no data used in stacking */

    /* calculate stacked waveform */
	for (j = 0; j < num; j++) {
		data_stack[j] = data_stack[j] / sta_num;
        /* nth root stack */
        /*if (data_stack[j] >= 0) {
			data_stack[j] = pow(fabs(data_stack[j]), nth );
		} else {
			data_stack[j] = 0 - pow(fabs(data_stack[j]), nth );
		}*/
	}

	return sta_num;
}


