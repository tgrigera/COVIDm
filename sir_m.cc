/*
 * sir_m.cc
 *
 * Monte Carlo simulation of simple mean-field stochastic SIR.  This
 * differs from sir.cc in that it can compute mean and variance over
 * many runs of the same system.  Also implementation details differ,
 * with this file using the Popstate classes for printing and
 * averaging.
 *
 * Choose discrete or continuous time (Gillespie) from the command line
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

#undef LATTICE_MC

///////////////////////////////////////////////////////////////////////////////
//
// simulation options and parameters

struct opt {
  int    last_arg_read;
  
  char   *ifile;
  bool   gillespie;
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

static int nargs=6;

///////////////////////////////////////////////////////////////////////////////
//
// Read parameters from command-line and file, and compute
// derived parameters

#include "read_arg.hh"

void show_usage(char *prog)
{
  std::cerr << "usage: " << prog << " [G or M] parameterfile seed N steps Nruns\n\n"
	    << "The first argument is G for Gillespie algorithm or M for discrete-time Monte Carlo\n";

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

  char *A;
  read_arg(argv,A);
  options.gillespie= (*A=='G');
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
 * Version with discrete-time MC
 *
 */

void run_MC(SIRstate *state)
{
  Uniform_real ran(0,1.);

  int S=options.S0*options.N;
  int I=options.I0*options.N;
  int R=options.N-S-I;
  double fS=(double) S/options.N;
  double fI=(double) I/options.N;
  double fR=(double) R/options.N;

  state->push(0,fS,fI,fR);

#ifdef LATTICE_MC
  short *lat=new short[options.N];
  int n=0;
  while (n<S) lat[n++]=0;
  while (n<S+I) lat[n++]=1;
  while (n<S+I+R) lat[n++]=2;

  for (int t=1; t<=options.steps; ++t) {
    for (n=0; n<options.N; n++)
      switch(lat[n]) {
      case 0:         // susceptible
  	if (ran()<options.beta*fI)
  	  { lat[n]=1; --S; ++I; }
  	break;
      case 1:         // infected
  	if (ran()<options.gamma)
  	  { lat[n]=2; --I; ++R; }
  	break;
      }

    fS=(double) S/options.N;
    fI=(double) I/options.N;
    fR=(double) R/options.N;
    state->push(t,fS,fI,fR);
  }
  delete[] lat;
#else
  for (int t=1; t<=options.steps; ++t) {
    for (int i=0; i<S; ++i)
      if (ran()<options.beta*I/options.N) { --S; ++I; }
    for (int i=0; i<I; ++i)
      if (ran()<options.gamma) { --I; ++R; }
    state->push(t,(double) S/options.N,(double) I/options.N,(double) R/options.N);
  }
#endif
}


/*
 * Version with Gillespie (kinetic MC) algorithm
 *
 */

void run_gillespie(SIRstate *state)
{
  Uniform_real ran(0,1.);

  int S=options.S0*options.N;
  int I=options.I0*options.N;
  int R=options.N-S-I;

  state->push(0,(double) S/options.N,(double) I/options.N,(double) R/options.N);

  double time=0;
  double last=0;
  Exponential_distribution rexp;
  while (time<options.steps) {

    // compute transition probabilities
    double pinf=options.beta*I*S/options.N;
    double prec=options.gamma*I;
    double pany=pinf+prec;
    if (pany==0.) break;

    // advance time
    double deltat=rexp(1./pany);
    time+=deltat;

    // choose the transition and apply it
    if (ran()<pinf/pany) {  // chose infeccion
      --S;
      ++I;
    } else { // chose recovery
      --I;
      ++R;
    }
    
    if (time>last+1.) {
      last=time;
      state->push(time,(double) S/options.N,(double) I/options.N,(double) R/options.N);
    }

  }
}


///////////////////////////////////////////////////////////////////////////////
//
// main

int main(int argc,char *argv[])
{
  read_parameters(argc,argv);
  Random_number_generator RNG(options.seed);

  SIRstate *state;
  state = options.Nruns>1 ?
    new SIRstate_av : new SIRstate;

  // Choose algorithm
  void (*runf)(SIRstate*);
  runf= options.gillespie ? run_gillespie : run_MC;

  std::cout << (options.gillespie ?
    "#\n# ***** Using Gillespie algorithm *****\n" :
		"#\n# ***** Using Monte Carlo algorithm with fixed time steps *****\n");
  std::cout << state->header() << '\n';

  // Do runs and print results
  for (int n=0; n<options.Nruns; ++n)
    (*runf)(state);

  if (options.Nruns>1)
    std::cout << *state;

  delete state;
}



