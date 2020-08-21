/*
 *  Get amplitude at the given time.
 *
 *  Author: Jiayuan Yao @ NTU
 *
 *  Revisions:
 *      2019-01-20  Jiayuan Yao  Initial Coding
 *      2020-08-17  Jiayuan Yao  Better function naming and coding style
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <float.h>
#include <math.h>
#include "sacio.h"
#include "const.h"

void usage(void);

void usage() {
    fprintf(stderr, "Get amplitude at the given time.                   \n");
    fprintf(stderr, "                                                   \n");
    fprintf(stderr, "Usage:                                             \n");
    fprintf(stderr, "  sacstack -Tt[/tmark/ts/tw]                       \n");
    fprintf(stderr, "                                                   \n");
    fprintf(stderr, "Options:                                           \n");
    fprintf(stderr, "  -T: arrival time [/tmark/begin time (sec)/time window (sec)     \n");
    fprintf(stderr, "  -h: show usage                                   \n");
}

int main(int argc, char *argv[])
{
    int c;
    int error;
    int cut = 0;   /* cut a time window or not */
    int tmark;
    float t0, t1, t;
    int i;

    error = 0;
    while ((c = getopt(argc, argv, "D:T:h")) != -1) {
        switch(c) {
            case 'T':
                if (sscanf(optarg, "%f/%d/%f/%f", &t, &tmark, &t0, &t1) == 4) {
                    cut = 1;
                } else if (sscanf(optarg, "%f", &t) == 1) {

                } else {
                    error++;
                }
                break;
            case 'h':
                usage();
                return -1;
            default:
                usage();
                return -1;
        }
    }

    if (argc-optind < 1 || error) {
        usage();
        return -1;
    }

    for (i=optind; i<argc; i++) {  /* loop over files */
        float *data;
        SACHEAD hd;

        if (cut) data = read_sac_pdw(argv[i], &hd, tmark, t0, t1);
        else     data = read_sac(argv[i], &hd);

        if (t < hd.b) {
            fprintf(stderr, "%s: t is smaller than begin time.\n", argv[i]);
        } else if (t > hd.e) {
            fprintf(stderr, "%s: t is larger than end time.\n", argv[i]);
        } else {
            printf("%s %f\n", argv[i], data[(int)((t - hd.b)/hd.delta)]);
        }
    }

    return 0;
}


