/*
 * sir_f.cc
 *
 * Monte Carlo simulation of stochastic SIR with population separated
 * in families.
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

#define LATTICE_MC

///////////////////////////////////////////////////////////////////////////////
//
// simulation options and parameters

struct opt {
  int    last_arg_read;
  
  char   *ifile;
  int    Nruns;
  int    steps;
  long   seed;

  int    Nfamilies; // Total number of families
  int    M;
  double beta_in,beta_out,gamma;
  double S0,I0;      // initial susceptible / infected fraction  

  opt() : last_arg_read(0) {}

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
  sscanf(buf,"%d %d",&options.Nfamilies,&options.M);
  buf=readbuf(f);
  sscanf(buf,"%lg %lg %lg",&options.beta_in,&options.beta_out,&options.gamma);
  buf=readbuf(f);
  sscanf(buf,"%lg %lg",&options.S0,&options.I0);
  fclose(f);

  printf("##### Parameters\n");
  printf("# beta_in = %g\n",options.beta_in);
  printf("# beta_out = %g\n",options.beta_out);
  printf("# gamma = %g\n",options.gamma);
  printf("# Nfamilies = %d\n",options.Nfamilies);
  printf("# M         = %d\n",options.M);
  printf("# S0 = %g\n",options.S0);
  printf("# I0 = %g\n",options.I0);
  printf("# Nruns = %d\n",options.Nruns);
}

///////////////////////////////////////////////////////////////////////////////
//
// Simulation

struct Gstate {
  int    N;         // Total population
  int    S,I,R;     // Total in state S,I,R
} gstate;

struct Fstate {
  int  M;
  int  S,I,R;
} ;

std::vector <Fstate> families;

void create_families()
{
  Uniform_real ran(0,1.);
  Fstate fam;

  double pI=options.I0;       // pI = probability of infection at t=0
  double pR=1-options.S0;     // pR = pI + probability of recovered at t=0

  gstate.N=gstate.S=gstate.I=gstate.R=0;

  for (int f=0; f<options.Nfamilies; ++f) {
    fam.S=fam.M=options.M;
    fam.I=fam.R=0;
    for (int i=0; i<fam.S; ++i) {
      double r=ran();
      if (r<pI) {fam.S--; fam.I++;}
      else if (r<pR) {fam.S--; fam.R++;}
    }
    gstate.N+=fam.M;
    gstate.S+=fam.S;
    gstate.I+=fam.I;
    gstate.R+=fam.R;
    families.push_back(fam);
  }

  printf("# N     = %d\n",gstate.N);
}

/* 
 * Version with discrete-time MC
 *
 */

#ifdef LATTICE_MC
void run(SIRstate *state)
{
  Uniform_real ran(0,1.);

  create_families();
  double fS=(double) gstate.S/gstate.N;
  double fI=(double) gstate.I/gstate.N;
  double fR=(double) gstate.R/gstate.N;

  state->push(0,fS,fI,fR);

  for (int t=1; t<=options.steps; ++t) {

    double pinf_ext=options.beta_out*fI;

    for (int f=1; f<families.size(); ++f) {
      int Sfam=families[f].S;
      int Ifam=families[f].I;
      for (int i=0; i<Sfam; ++i) {
	double pinf=pinf_ext + options.beta_in*Ifam/families[f].M;
	if (ran()<pinf) {
	  families[f].S--;
	  families[f].I++;
	  gstate.S--;
	  gstate.I++;
	}
      }
      for (int i=0; i<Ifam; ++i)
	if (ran()<options.gamma) {
	  families[f].I--;
	  families[f].R++;
	  gstate.I--;
	  gstate.R++;
	}
    }

     fS=(double) gstate.S/gstate.N;
     fI=(double) gstate.I/gstate.N;
     fR=(double) gstate.R/gstate.N;
     state->push(t,fS,fI,fR);
  }
}

#else

/*
 * Version with Gillespie (kinetic MC) algorithm
 *
 */

void run(SIRstate *state)
{
  std::cerr << "Unimplemented\n";
}

#endif /* LATTICE_MC */

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

#ifdef LATTICE_MC
  std::cout << "#\n# ***** Using Monte Carlo algorithm with fixed time steps *****\n";
#else
  std::cout << "#\n# ***** Using Gillespie algorithm *****\n";
#endif

  std::cout << state->header() << '\n';

  // Do runs and print results
  for (int n=0; n<options.Nruns; ++n)
    run(state);

  if (options.Nruns>1)
    std::cout << *state;

  delete state;
}
