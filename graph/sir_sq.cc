/*
 * sir_sq.cc -- SIR model on the square lattice
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
  int    I0;

  double beta;
  double gamma;

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
  sscanf(buf,"%lg %lg",&options.beta,&options.gamma);
  buf=readbuf(f);
  sscanf(buf,"%d %d %d",&options.Lx,&options.Ly,&options.I0);
  fclose(f);

  printf("##### Parameters\n");
  printf("# beta = %g\n",options.beta);
  printf("# gamma = %g\n",options.gamma);
  printf("# Lx, Ly = %d, %d\n",options.Lx,options.Ly);
  printf("# I0 = %d\n",options.I0);
  printf("# Number of runs = %d\n",options.Nruns);
}


int main(int argc,char* argv[])
{
  read_parameters(argc,argv);
  Random_number_generator RNG(options.seed);

#ifdef DEBUG_FCGRAPH
  FCGraph* egraph = FCGraph::create(options.Lx*options.Ly);
#else
  SQGraph* egraph = SQGraph::create(options.Lx,options.Ly);
#endif

  { // open scope so that SIR_model object is destroyed before calling
    // delete on the graph
    
#ifdef DEBUG_FCGRAPH
    SIR_model<FCGraph> SIR(*egraph);
    SIRcollector_av<FCGraph> collector(SIR);
#else
    SIR_model<SQGraph> SIR(*egraph);
    SIRcollector_av<SQGraph> collector(SIR);
#endif
    SIR.set_beta(options.beta);
    SIR.set_gamma(options.gamma);

    event_queue_t events;
    Imported_infection iinf(1,options.I0);
    events.push(&iinf);

    std::cout << collector.header() << '\n';
    for (int n=0; n<options.Nruns; ++n) {
      Sampler *sampler = options.Nruns>1 ?
	static_cast<Sampler*>( new Gillespie_sampler(0,options.steps,1.,&collector)  ) :
	new Passthrough_sampler(0,options.steps,&collector) ;
      run(&SIR,sampler,events,options.steps);
      delete sampler;
    }
    std::cout << collector;
  }
  
  delete egraph;
}
