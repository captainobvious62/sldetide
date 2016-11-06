/*
 * Copyright (c) 2012 Institute of Geological & Nuclear Sciences Ltd.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 */

/* system includes */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <math.h>

/* libmseed library includes */
#include <libmseed.h>
#include <libtidal.h>

#define PROGRAM "msdetide" /* program name */

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "xxx"
#endif

/*
 * msdetide: process raw miniseed tidal records to convert from recorded counts into sea level heights, with and without the tidal components
 *
 */

/* program variables */
static char *program_name = PROGRAM;
static char *program_version = PROGRAM " (" PACKAGE_VERSION ") (c) GNS 2012 (m.chadwick@gns.cri.nz)";
static char *program_usage = PROGRAM " [-hv][-A <alpha>][-B <beta>][-O <orient>][-L <latitude>][-Z <zone>][-T <label/amp/lag> ...][<files> ... ]";
static char *program_prefix = "[" PROGRAM "] ";

static int verbose = 0; /* program verbosity */

static char *orient = "T";
static double alpha = 0.0;
static double beta = 10.0;

static double zone = 0.0;
static double latitude = 0.0;
static int num_tides = 0;

static tidal_t tides[LIBTIDAL_MAX_CONSTITUENTS];

static void log_print(char *message) {
  if (verbose)
    fprintf(stderr, "%s", message);
}

static void err_print(char *message) {
  fprintf(stderr, "error: %s", message);
}

static void record_handler (char *record, int reclen, void *extra) {
  if (fwrite(record, reclen, 1, stdout) != 1) {
    ms_log (2, "error writing mseed record to stdout\n");
  }
}

static int detide_record(MSRecord *msr) {
    int n;

    MSTraceGroup *mstg;

    double epoch = 0.0;
    double height = 0.0;
    double value = 0.0;

    long precords = 0;
    long psamples = 0;

    /* need at least a sample */
    if (msr->samplecnt < 1)
        return 0;

    /* we need an integer value */
    if (msr->sampletype != 'i')
        return 0;

    /* problem with rate */
    if (msr->samprate == 0.0)
        return 0;

    /* convert to detided values ... */
    if (orient != NULL)
        msr->channel[2] = (*orient);

    epoch = (double) MS_HPTIME2EPOCH((msr->starttime));
    for (n = 0; n < msr->numsamples; n++) {
        libtidal_height(num_tides, tides, epoch + (double) n / msr->samprate, latitude, zone, &height);
        ((int *) msr->datasamples)[n] = (int) rint((double) ((int *) msr->datasamples)[n] - (alpha + beta * height));
    }

    if ((mstg = mst_initgroup(NULL)) == NULL)
        return -1;

    mst_addmsrtogroup(mstg, msr, 0, -1, -1);
    mst_printtracelist(mstg, 1, 1, 1);
    precords = mst_packgroup (mstg, record_handler, NULL, 512, DE_STEIM2, 1, &psamples, 1, 0, NULL);

    mst_freegroup (&mstg);

    return psamples;
}

int main(int argc, char **argv) {

  MSRecord *msr = NULL;

  int rc;
  int option_index = 0;
  struct option long_options[] = {
    {"help", 0, 0, 'h'},
    {"verbose", 0, 0, 'v'},
    {"alpha", 1, 0, 'A'},
    {"beta", 1, 0, 'B'},
    {"orient", 1, 0, 'O'},
    {"latitude", 1, 0, 'L'},
    {"zone", 1, 0, 'Z'},
    {"tide", 1, 0, 'T'},
    {0, 0, 0, 0}
  };

  /* adjust output logging ... -> syslog maybe? */
  ms_loginit (log_print, program_prefix, err_print, program_prefix);

  while ((rc = getopt_long(argc, argv, "hvA:B:O:T:L:Z:", long_options, &option_index)) != EOF) {
    switch(rc) {
    case '?':
      (void) fprintf(stderr, "usage: %s\n", program_usage);
      exit(-1); /*NOTREACHED*/
    case 'h':
      (void) fprintf(stderr, "\n[%s] miniseed tidal correction\n\n", program_name);
      (void) fprintf(stderr, "usage:\n\t%s\n", program_usage);
      (void) fprintf(stderr, "version:\n\t%s\n", program_version);
      (void) fprintf(stderr, "options:\n");
      (void) fprintf(stderr, "\t-h --help\tcommand line help (this)\n");
      (void) fprintf(stderr, "\t-v --verbose\trun program in verbose mode\n");
      (void) fprintf(stderr, "\t-A --alpha\tadd offset to calculated tidal heights [%g]\n", alpha);
      (void) fprintf(stderr, "\t-B --beta\tscale calculated tidal heights [%g]\n", beta);
      (void) fprintf(stderr, "\t-O --orient\talternative orientation code [%s]\n", orient);
      (void) fprintf(stderr, "\t-L --latitude\tprovide reference latitude [%g]\n", latitude);
      (void) fprintf(stderr, "\t-Z --zone\tprovide reference time zone offet [%g]\n", zone);
      (void) fprintf(stderr, "\t-T --tide\tprovide tidal constants [<label>/<amplitude>/<lag>]\n");
      exit(0); /*NOTREACHED*/
    case 'v':
      verbose++;
      break;
    case 'A':
      alpha = atof(optarg);
      break;
    case 'B':
      beta = atof(optarg);
      break;
    case 'L':
      latitude = atof(optarg);
      break;
    case 'Z':
      zone = atof(optarg);
      break;
    case 'O':
      orient = optarg;
      break;
    case 'T':
            if (num_tides < LIBTIDAL_MAX_CONSTITUENTS) {
                strncpy(tides[num_tides].name, strtok(strdup(optarg), "/"), LIBTIDAL_CHARLEN - 1);
                tides[num_tides].amplitude = atof(strtok(NULL, "/"));
                tides[num_tides].lag = atof(strtok(NULL, "/")) / 360.0;
                num_tides++;
            }
      break;
    }
  }

  /* report the program version */
  if (verbose)
    ms_log (0, "%s\n", program_version);

  if (verbose) {
      ms_log (0, "tidal [%s] zone=%g latitude=%g alpha=%g beta=%g\n", orient, zone, latitude, alpha, beta);
        for (rc = 0; rc < num_tides; rc++) {
      ms_log (0, "\t[%s] %g (%6.3f)\n", tides[rc].name, tides[rc].amplitude, tides[rc].lag);
        }
    }

    do {
        if (verbose)
      ms_log (0, "process miniseed data from %s\n", (optind < argc) ? argv[optind] : "<stdin>");

    while ((rc = ms_readmsr (&msr, (optind < argc) ? argv[optind] : "-", 0, NULL, NULL, 1, 1, (verbose > 1) ? 1 : 0)) == MS_NOERROR) {
      if (verbose > 1)
        msr_print(msr, (verbose > 2) ? 1 : 0);
      if (detide_record(msr) < 0)
                break;
    }
    if (rc != MS_ENDOFFILE )
        ms_log (2, "error reading stdin: %s\n", ms_errorstr(rc));

    /* Cleanup memory and close file */
    ms_readmsr (&msr, NULL, 0, NULL, NULL, 0, 0, (verbose > 1) ? 1 : 0);
    } while((++optind) < argc);

  /* closing down */
  if (verbose)
    ms_log (0, "terminated\n");

  /* done */
  return(0);
}
