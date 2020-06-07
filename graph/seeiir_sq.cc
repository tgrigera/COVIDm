/*
 * seeiir_sq.cc -- SEEIIR model on the square lattice
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
  long   seed;

  int    Lx,Ly;

  // imported infections
  typedef std::vector<Imported_infection> imported_infections_t;
  imported_infections_t                   imported_infections;
  std::string eifile;             // File to read imported infected cases

  // rates vs time
  typedef std::vector<Rate_constant_change<SQGraph>> rates_vs_time_t;
  rates_vs_time_t                           rates_vs_time;

  opt() : last_arg_read(0) {}

} options;

static int nargs=4;

///////////////////////////////////////////////////////////////////////////////
//
// Read parameters from command-line and file, and compute
// derived parameters

#include "../read_arg.hh"

void show_usage(char *prog)
{
  std::cerr << "usage: " << prog << " parameterfile seed steps Nruns\n\n";
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

  FILE *f=fopen(options.ifile,"r");
  if (f==0) throw std::runtime_error(strerror(errno));
  char *buf;
  printf("##### Parameters\n");

  buf=readbuf(f);
  sscanf(buf,"%d %d",&options.Lx,&options.Ly);
  printf("# Lx = %d\n",options.Lx);
  printf("# Ly = %d\n",options.Ly);

  printf("#\n# Nruns = %d\n",options.Nruns);

  // Imported infections
  buf=readbuf(f);
  options.eifile=buf;
  options.eifile.erase(options.eifile.end()-1);   // remove trailing newline
  read_imported_infections();

  printf("# Imported infections:\n");
  printf("# Time   Cases\n");
  int II=0;
  for (auto iir: options.imported_infections)
    printf("# %g %d\n",iir.time,II+=iir.new_cases);

  read_rates_vs_time(f);
  fclose(f);
  printf("#\n# Rate constatst:\n");
  printf("# time beta sigma_1 sigma_2 gamma_1 gamma_2\n");
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
  while (ungetc(fgetc(f),f)!=EOF) {
    char *buf=readbuf(f);
    if (sscanf(buf,"%lg %d",&etime,&eI)!=2) {
      std::cerr  << "couldn't read record: " << buf << "\n";
      throw std::runtime_error(strerror(errno));}
    options.imported_infections.push_back(Imported_infection(etime,eI-eIold));
    eIold=eI;
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
    options.rates_vs_time.push_back(Rate_constant_change<SQGraph>(time,b,s1,s2,g1,g2));
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
  
  auto ii_begin=options.imported_infections.begin();
  auto ii_end=options.imported_infections.end();
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

int main(int argc,char* argv[])
{
  read_parameters(argc,argv);
  Random_number_generator RNG(options.seed);

  SQGraph* egraph = SQGraph::create(options.Lx,options.Ly);
  SEEIIR_model<SQGraph> SEEIIR(*egraph);
  SEEIIRcollector<SQGraph> *collector =
    options.Nruns > 1 ?
    new SEEIIRcollector_av<SQGraph>(SEEIIR,1.) :
    new SEEIIRcollector<SQGraph>(SEEIIR);

  std::cout << collector->header() << '\n';
  for (int n=0; n<options.Nruns; ++n) {
    merge_events();
    Sampler *sampler =  new Gillespie_sampler(0,options.steps,1.,collector);
    run(&SEEIIR,sampler,event_queue,options.steps);
    delete sampler;
  }
  if (options.Nruns>1) std::cout << *collector;
}
