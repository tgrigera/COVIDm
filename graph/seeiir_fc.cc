/*
 * seeiir_fc.cc -- SEEIIR model on the fully-connected graph with random weights
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

#include <stdlib.h>
#include <string.h>

#include "emodel.hh"
#include "seir_collector.hh"

///////////////////////////////////////////////////////////////////////////////
//
// simulation options and parameters

struct opt {
  int    last_arg_read;
  
  char   *ifile;
  int    Nruns;
  int    steps;
  double deltat;
  long   seed;

  int    Nnodes;
  enum {exp} beta_distribution;
  double exp_mu;
  

  // Forced transitions
  typedef std::vector<Forced_transition> forced_transition_t;
  forced_transition_t                    forced_transitions;
  std::string eifile;             // File to read imported infected cases

  // rates vs time
  typedef std::vector<Rate_constant_change<MWFCGraph>> rates_vs_time_t;
  rates_vs_time_t                           rates_vs_time;

  opt() : last_arg_read(0), deltat(1.) {}

} options;

static int nargs=5;

///////////////////////////////////////////////////////////////////////////////
//
// Read parameters from command-line and file, and compute
// derived parameters

#include "../read_arg.hh"

void show_usage(char *prog)
{
  std::cerr << "usage: " << prog << " parameterfile seed steps Nruns delta_t\n\n";
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
void read_rates_vs_time(FILE*);

void read_parameters(int argc,char *argv[])
{
  if (argc!=nargs+1) show_usage(argv[0]);

  read_arg(argv,options.ifile);
  read_arg(argv,options.seed);
  read_arg(argv,options.steps);
  read_arg(argv,options.Nruns);
  read_arg(argv,options.deltat);

  FILE *f=fopen(options.ifile,"r");
  if (f==0) throw std::runtime_error(strerror(errno));
  char *buf;
  printf("##### Parameters\n");

  buf=readbuf(f);
  sscanf(buf,"%d",&options.Nnodes);
  printf("# N = %d\n",options.Nnodes);

  buf=readbuf(f);
  if (strncmp(buf,"exp",3)==0) {
    sscanf(buf+3,"%lg",&options.exp_mu);
    printf("# Drawing betas from exponential distribution with mu = %g\n",options.exp_mu);
  } else throw std::runtime_error("Invalid random distribution specified in parameter file");
  

  printf("#\n# Nruns = %d\n",options.Nruns);

  // Imported infections
  buf=readbuf(f);
  options.eifile=buf;
  options.eifile.erase(options.eifile.end()-1);   // remove trailing newline
  read_imported_infections();

  printf("# Infections and recoveries:\n");
  printf("# Time   Infections   Recoveries\n");
  int II=0;
  int RR=0;
  for (auto iir: options.forced_transitions)
    printf("# %g %d %d\n",iir.time,II+=iir.new_infected,RR+=iir.new_recovered);

  read_rates_vs_time(f);
  fclose(f);
  printf("#\n# Rate constants:\n");
  printf("# time beta_0 sigma_1 sigma_2 gamma_1 gamma_2\n");
  for (auto r:options.rates_vs_time)
    printf("# %g %g %g %g %g %g\n",r.time,r.beta,r.sigma1,r.sigma2,r.gamma1,r.gamma2);
}

void read_imported_infections()
{
  FILE *f=fopen(options.eifile.c_str(),"r");
  if (f==0) {
    std::cerr << "Error opening file (" << options.eifile << ")\n";
    throw std::runtime_error(strerror(errno));
  }

  double etime;
  int   eI,eIold=0;
  int   eR,eRold=0;
  while (ungetc(fgetc(f),f)!=EOF) {
    char *buf=readbuf(f);
    if (sscanf(buf,"%lg %d %d",&etime,&eI,&eR)!=3) {
      std::cerr  << "couldn't read record: " << buf << "\n";
      throw std::runtime_error(strerror(errno));}
    options.forced_transitions.push_back(Forced_transition(etime,eI-eIold,eR-eRold));
    eIold=eI;
    eRold=eR;
  }
  
  fclose(f);
}

void read_rates_vs_time(FILE *f)
{
  double time,b,s1,s2,g1,g2;
  while (ungetc(fgetc(f),f)!=EOF) {
    char *buf=readbuf(f);
    if (sscanf(buf,"%lg %lg %lg %lg %lg %lg",&time,&b,&s1,&s2,&g1,&g2)!=6) {
      std::cerr  << "couldn't read record: " << buf << "\n";
      throw std::runtime_error(strerror(errno));}
    options.rates_vs_time.push_back(Rate_constant_change<MWFCGraph>(time,b,s1,s2,g1,g2));
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// merge_events()

event_queue_t event_queue;

/*
 * build a time ordered queue of beta and imported infection changes
 * closes with dummy event at infinite time
 *
 */

void merge_events()
{
  while (!event_queue.empty()) event_queue.pop();
  
  auto ii_begin=options.forced_transitions.begin();
  auto ii_end=options.forced_transitions.end();
  auto ir_begin=options.rates_vs_time.begin();
  auto ir_end=options.rates_vs_time.end();
  
  auto ii=ii_begin;
  auto ir=ir_begin;
  while ( ii!=ii_end  || ir!=ir_end ) {

    while (ii!=ii_end && (ir==ir_end || ii->time <= ir->time ) ) {
      event_queue.push(&(*ii));
      ++ii;
    }

    while (ir!=ir_end && (ii==ii_end || ir->time <= ii->time) ) {
      event_queue.push(&(*ir));
      ++ir;
    }

  }
}

double beta_from_exp_dist()
{
  static Exponential_distribution edr(options.exp_mu);
  return edr();
}

int main(int argc,char* argv[])
{
  read_parameters(argc,argv);
  Random_number_generator RNG(options.seed);

  std::cerr << "# Building graph...\n";
  MWFCGraph* egraph = MWFCGraph::create(options.Nnodes);

  std::cerr << "#      ...setting weights\n";
  switch (options.beta_distribution) {
  case opt::exp:
    egraph->set_weights_random_multiplicative(beta_from_exp_dist,options.exp_mu);
    break;
  }

  std::cerr << "# Additional setup...\n";
  SEEIIR_model<MWFCGraph> *SEEIIR = new SEEIIR_model<MWFCGraph>(*egraph);
  SEEIIRcollector<MWFCGraph> *collector =
    options.Nruns > 1 ?
    new SEEIIRcollector_av<MWFCGraph>(*SEEIIR,options.deltat) :
    new SEEIIRcollector<MWFCGraph>(*SEEIIR);
  
  std::cerr << "# Starting run\n";
  std::cout << collector->header() << '\n';
  for (int n=0; n<options.Nruns; ++n) {
    merge_events();
    Sampler *sampler =  new Gillespie_sampler(0,options.steps,options.deltat,collector);
    run(SEEIIR,sampler,event_queue,options.steps);
    delete sampler;
  }
  if (options.Nruns>1) std::cout << *collector;
  
  delete collector;
  delete SEEIIR;
  delete egraph;
}
