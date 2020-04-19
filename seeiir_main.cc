/*
 * seeiir_main.cc
 *
 * Stochastic SEIR with two E and to I states.  Population separated
 * in families.  Simulated in continuous time (Gillespie algorithm).
 *
 * This is the driver for SEEIIR simulation (main and parameter
 * reading).  Use the #defines below (SEEIIR_IMPLEMENTATION) to choose
 * an implementation.  Implementations are included from files
 * seeir_i?.cc.
 *
 * In this model the states are S->E1->E2->I1->I2->R.  Exposed
 * individuals are those that will developthe illness but are still
 * not contagious (maybe not symptomatic).  The two E and I states are
 * completely identical, the serve to produce a non-exponential decay
 * from E to I and from I to R.  Individuals are grouped in families,
 * and a different infection probability is assigned to in-familiy and
 * out-of family contacts.  The transition rates (per individual) are
 *
 *  W(S->E1) = beta_out * (I1+I2) / (N-1)  + beta_in * (I1[f] + I2[f])
 *  W(E1->E2) = 2 * sigma
 *  W(E2->I1) = 2 * sigma
 *  W(I1->I2) = 2 * gamma
 *  W(I2->R) = 2 * gamma
 *
 * where N is the number of individuals, I1 and I2 the total number of
 * individuals in those states, and I1[f] the number of individuals in
 * state I1 belonging to the same family as the individual whose rate
 * is being computed.
 *
 * Command-line argumets are parameter file, seed, maximum time,
 * number of runs.  If number of runs is greater than 1, output is the
 * average and variance over all runs.  Ouput is to stdout.
 *
 * The parameter file specifies the number of families and the
 * distribution of family sizes (see seeiir_par.dat for the format),
 * the rate constants beta_in, beta_out, sigma, and gamma, and a file
 * to read imported infections.
 *
 * The population is initialized to all S, so infected individuals
 * must be introduced in order to start the epidemic.  The external
 * infections are read from a two-column file (see
 * imported_infections.dat) giving time and number of imported cases.
 * The second column is interpreted as giving cumulative total cases,
 * not new cases.  At the specifed time a number of S individuals are
 * forced to state E1 so that the total count of imported cases
 * matches the number given in the file.  Both times and cases must be
 * monotonically increasing.
 *
 *
 *
 * This file is part of COVIDm.
 *
 * COVIDm is copyright (C) 2020 by the authors (see file AUTHORS)
 * 
 * COVIDm is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License (GPL) as
 * published by the Free Software Foundation. You can use either
 * version 3, or (at your option) any later version.
 * 
 * COVIDm is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * For details see the file LICENSE.
 *
 */

// #undef SEEIIR_IMPLEMENTATION_1
// #undef SEEIIR_IMPLEMENTATION_2
// #undef SEEIIR_IMPLEMENTATION_3

#include <iostream>
#include <cstdio>
#include <list>
#include <queue>
#include <cmath>
#include <limits>
#include <algorithm>

#include <stdlib.h>
#include <string.h>

#include "qdrandom.hh"
#include "popstate.hh"
#include "bsearch.hh"


///////////////////////////////////////////////////////////////////////////////
//
// simulation options and parameters

struct opt {
  int    last_arg_read;
  
  char   *ifile;
  int    Nruns;
  int    steps;
  long   seed;

  int         Nfamilies;    // Total number of families
  int         Mmax;         // Maximum family size
  double      *PM;          // Family size distribution
  std::string eifile;       // File to read imported infected cases
  std::string betafile;     // File to read beta_out vs time

  // imported infections
  struct ei {double time; int I;} ;
  typedef std::queue<ei>                imported_infections_t;
  imported_infections_t                 imported_infections;
  // beta vs time
  struct eb {double time; double beta;} ;
  typedef std::queue<eb>                beta_vs_time_t;
  beta_vs_time_t                        beta_vs_time;

  // Epidemic parameters
  double beta_in,beta_out,sigma,gamma;

  opt() : last_arg_read(0) {}
  ~opt() {delete[] PM;}

} options;

///////////////////////////////////////////////////////////////////////////////
//
// Read parameters from command-line and file, and compute
// derived parameters

#include "read_arg.hh"

static int nargs=4;

void show_usage(char *prog)
{
  std::cerr << "usage: " << prog << " parameterfile seed steps Nruns\n\n"
    ;
  exit(1);
}

char *readbuf(FILE *f)
{
  static char buf[1000];
  char *s;
  do
    s=fgets(buf,1000,f);
  while (*buf=='#');
  return buf;
}

void read_imported_infections();
void read_beta_vs_time();

void read_parameters(int argc,char *argv[])
{
  if (argc!=nargs+1) show_usage(argv[0]);
  read_arg(argv,options.ifile);
  read_arg(argv,options.seed);
  read_arg(argv,options.steps);
  read_arg(argv,options.Nruns);

  FILE *f=fopen(options.ifile,"r");
  if (f==0) throw std::runtime_error(strerror(errno));
  char *buf=readbuf(f);
  sscanf(buf,"%d %d",&options.Nfamilies,&options.Mmax);
  options.PM=new double[options.Mmax+1];
  options.PM[0]=0;
  for (int M=1; M<=options.Mmax; ++M) {
    buf=readbuf(f);
    sscanf(buf,"%lg",options.PM+M);
  }
  buf=readbuf(f);
  sscanf(buf,"%lg %lg %lg %lg",&options.beta_in,&options.beta_out,&options.sigma,&options.gamma);
  buf=readbuf(f);
  options.eifile=buf;
  options.eifile.erase(options.eifile.end()-1);   // remove trailing newline
  read_imported_infections();
  if (options.beta_out<0) {                      // beta_out comes from beta vs time file
    buf=readbuf(f);
    options.betafile=buf;
    options.betafile.erase(options.betafile.end()-1);   // remove trailing newline
    read_beta_vs_time();
  }
  fclose(f);

  printf("##### Parameters\n");
  printf("# beta_in = %g\n",options.beta_in);
  printf("# beta_out = %g\n",options.beta_out);
  printf("# sigma = %g\n",options.sigma);
  printf("# gamma = %g\n",options.gamma);
  printf("# Nfamilies = %d\n",options.Nfamilies);
  printf("# Mmax      = %d\n",options.Mmax);
  for (int i=1; i<=options.Mmax; ++i)
    printf("# P[%d]    = %g\n",i,options.PM[i]);
  printf("# Nruns = %d\n",options.Nruns);
  printf("# Imported infections:\n");
  printf("# Time   Cases\n");
  opt::ei r;
  opt::imported_infections_t ii=options.imported_infections;
  while (!ii.empty()) {
    r=ii.front();
    printf("# %g %d\n",r.time,r.I);
    ii.pop();
  }
  if (options.beta_out<0) {
  printf("#\n# Beta_out:\n");
  printf("# Time   Beta_out\n");
    opt::beta_vs_time_t ib=options.beta_vs_time;
    opt::eb r;
    while (!ib.empty()) {
      r=ib.front();
      printf("# %g %g\n",r.time,r.beta);
      ib.pop();
    }
  }
  printf("#\n# Nruns = %d\n",options.Nruns);
}

void read_imported_infections()
{
  FILE *f=fopen(options.eifile.c_str(),"r");
  if (f==0) throw std::runtime_error(strerror(errno));

  opt::ei ei;
  while (ungetc(fgetc(f),f)!=EOF) {
    char *buf=readbuf(f);
    if (sscanf(buf,"%lg %d",&ei.time,&ei.I)!=2) {
      std::cerr  << "couldn't read record: " << buf << "\n";
      throw std::runtime_error(strerror(errno));}
    options.imported_infections.push(ei);
  }
  
  fclose(f);
}

void read_beta_vs_time()
{
  FILE *f=fopen(options.betafile.c_str(),"r");
  if (f==0) throw std::runtime_error(strerror(errno));

  opt::eb bi;
  while (ungetc(fgetc(f),f)!=EOF) {
    char *buf=readbuf(f);
    if (sscanf(buf,"%lg %lg",&bi.time,&bi.beta)!=2) {
      std::cerr  << "couldn't read record: " << buf << "\n";
      throw std::runtime_error(strerror(errno));}
    options.beta_vs_time.push(bi);
  }
  fclose(f);
}

///////////////////////////////////////////////////////////////////////////////
//
// merge_events()

struct event {
  double time;
  enum   {infection,beta_change} kind;
  int    I;
  double beta;
} ;

typedef std::queue<event> event_queue_t;
event_queue_t             event_queue;

//
// build a time ordered queue of beta and imported infection changes
// closes with dummy event at infinite time
void merge_events()
{
  while (!event_queue.empty()) event_queue.pop();
  
  opt::imported_infections_t ii=options.imported_infections;
  opt::beta_vs_time_t        bt=options.beta_vs_time;

  while (!ii.empty() || !bt.empty() ) {

    while (!ii.empty() && (bt.empty() || ii.front().time<=bt.front().time) ) {
      event_queue.push({ii.front().time,event::infection,ii.front().I,-1});
      ii.pop();
    }

    while (!bt.empty() && (ii.empty() || bt.front().time<=ii.front().time) ) {
      event_queue.push({bt.front().time,event::beta_change,-1000,bt.front().beta});
      bt.pop();
    }

  }

  // dummy event
  event_queue.push({std::numeric_limits<double>::max(),event::infection,0,0});
}

#ifdef SEEIIR_IMPLEMENTATION_1
#include "seeiir_i1.cc"
#endif
#ifdef SEEIIR_IMPLEMENTATION_2
#include "seeiir_i2.cc"
#endif
#ifdef SEEIIR_IMPLEMENTATION_3
#include "seeiir_i3.cc"
#endif

///////////////////////////////////////////////////////////////////////////////
//
// main

int main(int argc,char *argv[])
{
  read_parameters(argc,argv);
  Random_number_generator RNG(options.seed);
  Uniform_real ran(0,1);

  // Prepare global state (for output) and population
  SEEIIRstate *state;
  state = options.Nruns>1 ?
          new SEEIIRstate_av : new SEEIIRstate;
  std::cout << state->header() << '\n';

  SEIRPopulation pop(options.Nfamilies,options.beta_in,options.beta_out,
		     options.sigma,options.gamma,options.Mmax,options.PM);

  merge_events();

  // Do runs and print results
  for (int n=0; n<options.Nruns; ++n) {
    // std::cout << "# N = " << pop.gstate.N << '\n';
    run(pop,state);
    pop.set_all_S();
  }

  if (options.Nruns>1)
    std::cout << *state;

  delete state;
}
