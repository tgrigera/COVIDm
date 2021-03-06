/*
 * sir_f.cc
 *
 * Stochastic SIR with population separated in families, simulated in
 * continuous time (Gillespie algorithm).
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
#include "bsearch.hh"
#include "gillespie_sampler.hh"


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
  int    Mmax;
  double *PM;
  double beta_in,beta_out,gamma;
  double I0;      // initial infected fraction  

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
  sscanf(buf,"%lg %lg %lg",&options.beta_in,&options.beta_out,&options.gamma);
  buf=readbuf(f);
  sscanf(buf,"%lg",&options.I0);
  fclose(f);

  printf("##### Parameters\n");
  printf("# beta_in = %g\n",options.beta_in);
  printf("# beta_out = %g\n",options.beta_out);
  printf("# gamma = %g\n",options.gamma);
  printf("# Nfamilies = %d\n",options.Nfamilies);
  printf("# Mmax      = %d\n",options.Mmax);
  for (int i=1; i<=options.Mmax; ++i)
    printf("# P[%d]    = %g\n",i,options.PM[i]);
  printf("# I0 = %g\n",options.I0);
  printf("# Nruns = %d\n",options.Nruns);
}

///////////////////////////////////////////////////////////////////////////////
//
// Simulation


/*
 * Global and family state
 *
 */
struct Gstate {
  int    N;         // Total population
  int    S,I,R;     // Total in state S,I,R
} gstate;

struct Family {
  int  M;
  int  S,I,R;
} ;

std::ostream& operator<<(std::ostream& o,const Family &f)
{
  std::cout << "M = " << f.M  << " S,I,R = " << f.S << ' ' << f.I << ' ' << f.R;
  return o;
}

/*
 * class Population holds per family information, computes rates and
 * performs individual state switchs (function event)
 *
 */
class Population {
public:
  Population(int NFamilies,double beta_in,double beta_out,double gamma,
	     int Mmmax,double PM[]);
  ~Population() { delete Mdist;}
  
  void rebuild_families();
  void set_all_S();
  void compute_rates();      // compute all rates from scratch
  void event(int f,double r);  // perform an event in family f, update family's rates
  void infect(int f);          // Force infection in family f (unless no susceptibles)
  
  double beta_in,beta_out,gamma;
  int    NFamilies,Mmax;
  Discrete_distribution *Mdist;

  Gstate              gstate;
  std::vector<Family> families;

  std::vector<double> cumrate;
  double              total_rate;
} ;
		 
Population::Population(int NFamilies,double beta_in,double beta_out,double gamma,int Mmax,
		       double *P) :
  beta_in(beta_in),
  beta_out(beta_out),
  gamma(gamma),
  NFamilies(NFamilies),
  Mmax(Mmax),
  Mdist(0)
{
  Mdist = new Discrete_distribution(Mmax+1,P);
  rebuild_families();
}

void Population::rebuild_families()
{
  gstate.N=gstate.S=gstate.I=gstate.R=0;

  Family fam;
  families.clear();

  for (int f=0; f<NFamilies; ++f) {
    fam.S=fam.M=(*Mdist)();
    fam.I=fam.R=0;
    gstate.N+=fam.M;
    gstate.S+=fam.S;
    families.push_back(fam);
  }
  cumrate.resize(families.size()+1,0.);
}

inline void Population::infect(int f)  // Infect someone in family f
{
  if (families[f].S==0) return;
  families[f].S--;
  families[f].I++;
  gstate.S--;
  gstate.I++;
}

void Population::set_all_S()
{
  gstate.S=gstate.N;
  gstate.I=gstate.R=0;

  for (auto &f: families) {
    f.S=f.M;
    f.I=f.R=0;
  }
}

/*
 * compute_rates() computes the cumulative rates for all families,
 * event performs the inteded state change within the given family,
 * given a random number that is compared to the individual rates
 *
 */
void Population::compute_rates()
{
  int N1=gstate.N-1;
  cumrate[0]=0;
  for (int f=0; f<families.size(); ++f) {
    double famrate = families[f].S * (beta_out*gstate.I/N1 + beta_in*families[f].I ) +
      families[f].I*gamma;
    cumrate[f+1]=cumrate[f]+famrate;
  }
  total_rate=cumrate[families.size()];
}

void Population::event(int f,double r)   // Do infection or recovery on family f
                                         // and update the family rate
{
  double rr = families[f].I*gamma;

  if (r<rr) {                   // recovery
    families[f].I--;
    gstate.I--;
    families[f].R++;
    gstate.R++;
  } else {                      // infection
    families[f].S--;
    gstate.S--;
    families[f].I++;
    gstate.I++;
  }
}

void run(Population &pop,SIRstate *state)
{
  SIRistate istate;
  Gillespie_sampler<SIRstate,SIRistate> gsamp(*state,0.,options.steps,1.);
  
  istate.S=(double) pop.gstate.S/pop.gstate.N;
  istate.I=(double) pop.gstate.I/pop.gstate.N;
  istate.R=(double) pop.gstate.R/pop.gstate.N;
  gsamp.push_data(istate);

  double time=0;
  Exponential_distribution rexp;
  Uniform_real ran(0,1.);

  while (time<options.steps) {

    // compute transition probabilities
    pop.compute_rates();
    double mutot=pop.total_rate;
    // advance time
    double deltat=rexp(1./mutot);
    time+=deltat;
    gsamp.push_time(time);

    // choose the transition and apply it
    double r=ran()*mutot;
    int f=bsearch(r,pop.cumrate);      // choose family
    pop.event(f,r-pop.cumrate[f]);     // choose and apply event
    
    istate.S=(double) pop.gstate.S/pop.gstate.N;
    istate.I=(double) pop.gstate.I/pop.gstate.N;
    istate.R=(double) pop.gstate.R/pop.gstate.N;
    gsamp.push_data(istate);
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// main

int main(int argc,char *argv[])
{
  read_parameters(argc,argv);
  Random_number_generator RNG(options.seed);
  Uniform_real ran(0,1);

  // Prepare global state (for output) and population
  SIRstate *state;
  state = options.Nruns>1 ?
    new SIRstate_av : new SIRstate;
  Population pop(options.Nfamilies,options.beta_in,options.beta_out,options.gamma,
		 options.Mmax,options.PM);

  std::cout << "# N = " << pop.gstate.N << '\n';
  std::cout << state->header() << '\n';

  // Do runs and print results
  for (int n=0; n<options.Nruns; ++n) {
    // seed infected
    pop.set_all_S();
    for (int f=0; f<pop.families.size(); ++f) {
      for (int i=0; i<pop.families[f].M; ++i)
	if (ran()<options.I0) pop.infect(f);
    }

    // run and accumulate averages
    run(pop,state);
  }

  if (options.Nruns>1)
    std::cout << *state;

  delete state;
}
