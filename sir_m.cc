/*
 * sir_m.cc
 *
 * Continuous time (Gillespie) Monte Carlo simulation of simple
 * mean-field stochastic SIR.  This differs from sir.cc in that it can
 * compute mean and variance over many runs of the same system.  Also
 * implementation details differ, with this file using the Popstate
 * classes for printing and averaging.
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

#include <iostream>
#include <cstdio>
#include <math.h>

#include <stdlib.h>
#include <string.h>

#include "qdrandom.hh"
#include "popstate.hh"
#include "gillespie_sampler.hh"

///////////////////////////////////////////////////////////////////////////////
//
// simulation options and parameters

struct opt {
  int    last_arg_read;
  
  char   *ifile;
  int    Nruns;
  int    N;
  int    steps;
  long   seed;

  double R0;
  double inf_time;   // infection time (days)
  double S0,I0;      // initial susceptible / infected fraction  

  double beta;
  double gamma;

  opt() : last_arg_read(0) {}

} options;

static int nargs=5;

///////////////////////////////////////////////////////////////////////////////
//
// Read parameters from command-line and file, and compute
// derived parameters

#include "read_arg.hh"

void show_usage(char *prog)
{
  std::cerr << "usage: " << prog << " parameterfile seed N steps Nruns\n\n";
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

void read_parameters(int argc,char *argv[])
{
  if (argc!=nargs+1) show_usage(argv[0]);

  read_arg(argv,options.ifile);
  read_arg(argv,options.seed);
  read_arg(argv,options.N);
  read_arg(argv,options.steps);
  read_arg(argv,options.Nruns);

  FILE *f=fopen(options.ifile,"r");
  if (f==0) throw std::runtime_error(strerror(errno));
  char *buf=readbuf(f);
  sscanf(buf,"%lg %lg",&options.R0,&options.inf_time);
  buf=readbuf(f);
  sscanf(buf,"%lg %lg",&options.S0,&options.I0);
  fclose(f);

  // Compute derived parameters
  options.gamma=1./options.inf_time;
  options.beta=options.R0*options.gamma;

  printf("##### Parameters\n");
  printf("# R0 = %g\n",options.R0);
  printf("# inf_time = %g\n",options.inf_time);
  printf("# S0 = %g\n",options.S0);
  printf("# I0 = %g\n",options.I0);
  printf("# beta = %g\n",options.beta);
  printf("# gamma = %g\n",options.gamma);
}

///////////////////////////////////////////////////////////////////////////////
//
// Simulation

/*
 * Version with Gillespie (kinetic MC) algorithm
 *
 */

template <bool multipleruns>
void run(SIRstate *state)
{
  Uniform_real ran(0,1.);

  int S=options.S0*options.N;
  int I=options.I0*options.N;
  int R=options.N-S-I;

  Gillespie_sampler<SIRstate,SIRistate,multipleruns> gsamp(*state,0,options.steps,1);

  SIRistate gstate;
  gstate.S=(double) S/options.N;
  gstate.I=(double) I/options.N;
  gstate.R=(double) R/options.N;
  gsamp.push_data(gstate);

  double time=0;
  double last=0;
  Exponential_distribution rexp;
  while (time<options.steps) {

    // compute transition probabilities
    double pinf=options.beta*I*S/options.N;
    double prec=options.gamma*I;
    double pany=pinf+prec;

    // advance time
    double deltat=rexp(1./pany);
    time+=deltat;
    gsamp.push_time(time);

    // choose the transition and apply it
    if (ran()<pinf/pany) {  // chose infeccion
      --S;
      ++I;
    } else { // chose recovery
      --I;
      ++R;
    }

    gstate.S=(double) S/options.N;
    gstate.I=(double) I/options.N;
    gstate.R=(double) R/options.N;
    gsamp.push_data(gstate);

  }
}


///////////////////////////////////////////////////////////////////////////////
//
// main

int main(int argc,char *argv[])
{
  read_parameters(argc,argv);
  Random_number_generator RNG(options.seed);

  SIRstate_av state;
  std::cout << state.header() << '\n';

  // Do runs and print results
  if (options.Nruns>1)
    for (int n=0; n<options.Nruns; ++n)
      run<true>(&state);
  else
    run<false>(&state);

  std::cout << state;
}
