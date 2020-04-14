/*
 * seeiir_main.cc
 *
 * This is the driver for SEEIIR simulation (main and parameter
 * reading).  Use the #defines below (SEEIIR_IMPLEMENTATION) to choose
 * an imlementation.
 *
 * Stochastic SEIR with two E and to I states.  Population separated
 * in families.  Simulated in continuous time (Gillespie algorithm).
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
// #undef  SEEIIR_IMPLEMENTATION_2

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

  // imported infections
  struct ei {double time; int I;} ;
  typedef std::queue<ei>                imported_infections_t;
  imported_infections_t                 imported_infections;

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
  fclose(f);

  read_imported_infections();

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
  printf("# Nruns = %d\n",options.Nruns);
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

#ifdef SEEIIR_IMPLEMENTATION_1
#include "seeiir_i1.cc"
#endif
#ifdef SEEIIR_IMPLEMENTATION_2
#include "seeiir_i2.cc"
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
