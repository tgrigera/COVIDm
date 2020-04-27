/*
 * seeiir_h.cc
 *
 * Hierarchical stochastic SEIR with two E and to I states.
 * Population grouped in families, nieghbourhoods, towns, etc..
 * Simulated in continuous time (Gillespie algorithm).
 *
 * This implementation uses LEMON graphs to build the hierarchy tree.
 *
 * Hierachy tree has NL levels, with level 0 being the leaves
 * (individuals).  Leaves are not actually stored, only the level 1
 * nodes (families).  Each node holds a cumulative count of
 * individuals and their epidemiological state in the subtree to which
 * the node is root.
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
#include <list>
#include <queue>
#include <cmath>
#include <limits>
#include <algorithm>
#include <type_traits>

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <lemon/list_graph.h>

#include "qdrandom.hh"
#include "popstate.hh"
#include "bsearch.hh"

///////////////////////////////////////////////////////////////////////////////
//
// reading and setting options and parameters

/*
 * rates_t
 *
 */
struct rates_t {
  double time;
  std::vector<double> beta;
  double sigma1,sigma2,gamma1,gamma2;

  rates_t(int levels) :
    time(0.), sigma1(0), sigma2(0), gamma1(0), gamma2(0)
  {beta.resize(levels+1,0.);}
} ;

std::ostream& operator<<(std::ostream& o,const rates_t &r)
{
  o << "time = " << r.time << " beta_1 ... beta_" << r.beta.size()-1 << " = ";
  for (int i=1; i<r.beta.size(); ++i)  o << r.beta[i] << " ";
  o << "sigma_1 sigma_2 " << r.sigma1 << " " << r.sigma2
    << " gamma_1 gamma_2 " << r.gamma1 << " " << r.gamma2 << "\n";
  return o; 
}


/*
 * options
 *
 */
struct opt {
  int    last_arg_read;
  
  char   *ifile;
  int    Nruns;
  int    steps;
  long   seed;

  int         levels;       // tree depth (not counting individuals)
  int         Nfamilies;    // Total number of families
  int         Mmax;         // Maximum family size
  double      *PM;          // Family size distribution
  std::string eifile;       // File to read imported infected cases

  // imported infections
  struct ei {double time; int I;} ;
  typedef std::vector<ei>               imported_infections_t;
  imported_infections_t                 imported_infections;
  // rates vs time
  typedef std::vector<rates_t>          rates_vs_time_t;
  rates_vs_time_t                       rates_vs_time;

  opt() : last_arg_read(0) {}
  ~opt() {delete[] PM;}

} options;

/*
 * Read parameters from command-line and file, and compute
 * derived parameters
 *
 */
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
  static char buf[5000];
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

  options.levels=2;

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
  options.eifile=buf;
  options.eifile.erase(options.eifile.end()-1);   // remove trailing newline
  read_rates_vs_time(f);
  fclose(f);
  read_imported_infections();

  printf("##### Parameters\n");
  printf("# Nfamilies = %d\n",options.Nfamilies);
  printf("# Mmax      = %d\n",options.Mmax);
  for (int i=1; i<=options.Mmax; ++i)
    printf("# P[%d]    = %g\n",i,options.PM[i]);
  printf("# Nruns = %d\n",options.Nruns);
  printf("# Imported infections:\n");
  printf("# Time   Cases\n");
  for (auto iir: options.imported_infections)
    printf("# %g %d\n",iir.time,iir.I);

  printf("#\n# Rate constatst:\n");
  printf("# time ");
  for (int i=1; i<=options.levels; ++i) printf("beta_%d ",i);
  printf("sigma_1 sigma_2 gamma_1 gamma_2\n");
  for (auto r:options.rates_vs_time) {
    printf("# %g ",r.time);
    for (int i=1; i<=options.levels; ++i) printf("%g ",r.beta[i]);
    printf("%g %g %g %g\n",r.sigma1,r.sigma2,r.gamma1,r.gamma2);
  }

  printf("#\n# Nruns = %d\n",options.Nruns);
}

void read_imported_infections()
{
  FILE *f=fopen(options.eifile.c_str(),"r");
  if (f==0) {
    std::cerr << "Error opening file (" << options.eifile << ")\n";
    throw std::runtime_error(strerror(errno));
  }

  opt::ei ei;
  while (ungetc(fgetc(f),f)!=EOF) {
    char *buf=readbuf(f);
    if (sscanf(buf,"%lg %d",&ei.time,&ei.I)!=2) {
      std::cerr  << "couldn't read record: " << buf << "\n";
      throw std::runtime_error(strerror(errno));}
    options.imported_infections.push_back(ei);
  }
  
  fclose(f);
}

void read_rates_vs_time(FILE *f)
{
  rates_t rates(options.levels);
  std::string fmt;
  int ncread=0;
  
  while (ungetc(fgetc(f),f)!=EOF) {
    char *buf=readbuf(f);
    if (sscanf(buf,"%lg %n",&rates.time,&ncread)!=1) {
      std::cerr  << "couldn't read record: " << buf << "\n";
      throw std::runtime_error(strerror(errno));}
    for (int i=1; i<=options.levels; ++i) {
      buf+=ncread;
      if (sscanf(buf,"%lg %n",&(rates.beta[i]),&ncread)!=1 ) {
	std::cerr  << "couldn't read record: " << buf << "\n";
	throw std::runtime_error(strerror(errno));}
    }
    buf+=ncread;
    if (sscanf(buf,"%lg %lg %lg %lg",&(rates.sigma1),&(rates.sigma2),&(rates.gamma1),
	       &(rates.gamma2) )!=4) {
	std::cerr  << "couldn't read record: " << buf << "\n";
	throw std::runtime_error(strerror(errno));}
    options.rates_vs_time.push_back(rates);
  }
}


///////////////////////////////////////////////////////////////////////////////
//
// merge_events()

struct event {
  double  time;
  enum    {infection,rate_change} kind;
  long    enumber;
} ;

typedef std::queue<event> event_queue_t;
event_queue_t             event_queue;

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
      event_queue.push({ii->time,event::infection,ii-ii_begin});
      ++ii;
    }

    while (ir!=ir_end && (ii==ii_end || ir->time <= ii->time) ) {
      event_queue.push({ir->time,event::rate_change,ir-ir_begin});
      ++ir;
    }

  }

  // dummy event
  event_queue.push({std::numeric_limits<double>::max(),event::infection,0});
}


///////////////////////////////////////////////////////////////////////////////
//
// SEIRPopulation: holds population state, computes rates, performs
// transitions

typedef lemon::ListDigraph            graph_t;
typedef graph_t::Node                 node_t;
typedef std::vector<node_t>           node_list_t;

struct global_data {
  int              infections_imported;
  std::vector<int> infections_level;

  global_data(int levels) :
    infections_imported(0)
  {infections_level.resize(levels+1,0);}
} ;

struct node_data {
  int     level;
  int     M;   // number of direct descendents 
  int     N;   // number of cumulative descendents
  int     S,E1,E2,I1,I2,R;
  size_t  first_S_in_list;
  size_t  level_nodes_in_list;
  size_t  infected_nodes_in_list;

  node_data() :
    level(0), N(0), M(0),
    S(0), E1(0), E2(0), I1(0), I2(0), R(0),
    first_S_in_list(-1), level_nodes_in_list(-1), infected_nodes_in_list(-1)
  {}
} ;

std::ostream& operator<<(std::ostream& o,const node_data &nd)
{
  std::cout << "level = " << nd.level << " N = " << nd.N << " M = " << nd.M
    	    << "\n       S, E1, E2, I1, I2, R "
	    << nd.S << ' ' << nd.E1 << ' ' << nd.E2 << ' '
	    << nd.I1 << ' ' << nd.I2 << ' ' << nd.R << '\n';
  return o; 
}

struct epidemiological_event {
  enum {SE1,E1E2,E2I1,I1I2,I2R} type;
  lemon::ListDigraph::Node      node;
} ;

/*
 * class SEIRPopulation
 *
 */
class SEIRPopulation {
public:

  SEIRPopulation(int levels,int (*noffspring)(int));
  void rebuild_hierarchy();
  void set_all_S();
  void set_rate_parameters(rates_t& r) {rates=r;}
  void add_imported(int I);
  void recompute_counts();
  void check_structures();
  void compute_rates();
  void apply_event(int evn);

  int                 levels;
  std::vector<double> cumrate;
  std::vector<epidemiological_event>  events;
  double              total_rate;

  graph_t                     tree;
  node_t                      root;
  graph_t::NodeMap<node_data> treemap;
  global_data                 gdata;
  rates_t                     rates;

private:
  int (*noffspring)(int);
  Uniform_integer                        ran;
  std::vector<node_list_t>               level_nodes;
  std::vector<node_t>                    listS,listE1,listE2,listI1,listI2;
  std::vector<node_t>                    infected_nodes;
  
  node_t build_tree(int level);

  struct readS {
    static int& field(node_data &nd) {return nd.S;}
  } ;

  struct readE1 {
    static int& field(node_data &nd) {return nd.E1;}
  } ;

  struct readE2 {
    static int& field(node_data &nd) {return nd.E2;}
  } ;

  struct readI1 {
    static int& field(node_data &nd) {return nd.I1;}
  } ;
  
  struct readI2 {
    static int& field(node_data &nd) {return nd.I2;}
  } ;

  struct readR {
    static int& field(node_data &nd) {return nd.R;}
  } ;
  
  template <typename readF1,typename readF2>
  void update_counts(node_t l0node);
  void update_after_erase_susceptible(node_t node);
} ;

SEIRPopulation::SEIRPopulation(int levels,int (*noffspring)(int) ) :
  levels(levels),
  noffspring(noffspring),
  rates(levels),
  gdata(levels),
  treemap(tree)
{
  rebuild_hierarchy();
}

void SEIRPopulation::rebuild_hierarchy()
{
  tree.clear();
  level_nodes.clear();
  level_nodes.resize(levels+1);
  root=build_tree(levels);
  set_all_S();
}

// Recursively build a tree starting at given level
node_t SEIRPopulation::build_tree(int level)
{
  node_t subtree=tree.addNode();
  treemap[subtree].level=level;
  treemap[subtree].level_nodes_in_list=level_nodes[level].size();
  level_nodes[level].push_back(subtree);
  int M=noffspring(level);
  treemap[subtree].M=M;

  if (level>1) {  // we don't store the leaves
    for (int i=0; i<M; ++i) {
      node_t n;
      n=build_tree(level-1);
      tree.addArc(subtree,n);
    }
  }

  return subtree;
}

void SEIRPopulation::set_all_S()
{
  for (graph_t::NodeIt node(tree); node!=lemon::INVALID; ++node) {
    node_data& noded=treemap[node];
    noded.N = noded.level==1 ? noded.M : 0;
    noded.S=noded.N;
    noded.E1=noded.E2=noded.I1=noded.I2=noded.R=0;
  }
  recompute_counts();
  gdata.infections_imported=0;
  gdata.infections_level.resize(levels+1,0);
}

void SEIRPopulation::recompute_counts()
{
  listS.clear();
  listE1.clear();
  listE2.clear();
  listI1.clear();
  listI2.clear();
  infected_nodes.clear();

  for (auto &node: level_nodes[1]) {
    node_data& noded=treemap[node];
    noded.N=noded.M;
    noded.first_S_in_list=listS.size();
    for (int i=0; i<noded.S; ++i) listS.push_back(node);
    for (int i=0; i<noded.E1; ++i) listE1.push_back(node);
    for (int i=0; i<noded.E2; ++i) listE2.push_back(node);
    for (int i=0; i<noded.I1; ++i) listI1.push_back(node);
    for (int i=0; i<noded.I2; ++i) listI2.push_back(node);
    if (noded.I1+noded.I2>0) infected_nodes.push_back(node);
  }
    
  for (int level=2; level<=levels; ++level) {
    int Nprev=0;
    for (auto &node: level_nodes[level]) {
      node_data& noded=treemap[node];
      noded.N=noded.S=noded.E1=noded.E2=noded.I1=noded.I2=noded.R=0;
      noded.first_S_in_list=Nprev;
      for (graph_t::OutArcIt arc(tree,node); arc!=lemon::INVALID; ++arc) {
	node_t son=tree.target(arc);
	node_data& sond=treemap[son];
	noded.N+=sond.N;
	noded.S+=sond.S;
	noded.E1+=sond.E1;
	noded.E2+=sond.E2;
	noded.I1+=sond.I1;
	noded.I2+=sond.I2;
	noded.R+=sond.R;
      }
      if (noded.I1+noded.I2>0) infected_nodes.push_back(node);
      Nprev=noded.N;
    }
  }
}

void SEIRPopulation::check_structures()
{
  node_data& rootd=treemap[root];
  
  std::cerr << rootd;
  // for (lemon::ListDigraph::NodeIt n(tree); n != lemon::INVALID; ++n)
  //   std::cout << treemap[n];

  assert(listS.size()==rootd.S);
  for (node_t &node: listS) {
    assert(treemap[node].level==1);
    assert(treemap[node].S>0);
  }
  assert(listE1.size()==rootd.E1);
  for (node_t &node: listE1) {
    assert(treemap[node].level==1);
    assert(treemap[node].E1>0);
  }
  assert(listE2.size()==rootd.E2);
  for (node_t &node: listE2) {
    assert(treemap[node].level==1);
    assert(treemap[node].E2>0);
  }
  assert(listI1.size()==rootd.I1);
  for (node_t &node: listI1) {
    assert(treemap[node].level==1);
    assert(treemap[node].I1>0);
  }
  assert(listI2.size()==rootd.I2);
  for (node_t &node: listI2) {
    assert(treemap[node].level==1);
    assert(treemap[node].I2>0);
  }
  for (node_t &node: infected_nodes) {
    assert(treemap[node].I1+treemap[node].I2>0);
  }
  
}

/*
 * Rates and events
 *
 */
void SEIRPopulation::compute_rates()
{
  cumrate.clear();
  cumrate.push_back(0.);
  events.clear();

  epidemiological_event ev = {epidemiological_event::SE1,root};
  double cr=0;

  // compute infection rate at all levels
  //  for (int in=0; in<infected_nodes.size(); ++in) {
  for (node_t node: infected_nodes) {
    ev.node=node;
    node_data& noded=treemap[node];
    double norm = noded.level >1 ? 1./(noded.N-1) : 1;
    cr += noded.S * rates.beta[noded.level] * (noded.I1 + noded.I2) * norm;
    cumrate.push_back(cr);
    events.push_back(ev);
  }
  // for (int lev=levels; lev>0; --lev) {
  //   for (node_t &node: level_nodes[lev]) {
  //     ev.node=node;
  //     node_data& noded=treemap[node];
  //     double norm = lev>1 ? 1./(noded.N-1) : 1;
  //     cr += noded.S * rates.beta[lev] * (noded.I1 + noded.I2) * norm;
  //     cumrate.push_back(cr);
  //     events.push_back(ev);
  //   }
  // }

  // the other events are only global
  ev.node=root;
  // E1->E2
  ev.type=epidemiological_event::E1E2;
  cr += treemap[root].E1 * rates.sigma1;
  cumrate.push_back(cr);
  events.push_back(ev);
  // E2->I1
  ev.type=epidemiological_event::E2I1;
  cr += treemap[root].E2 * rates.sigma2;
  cumrate.push_back(cr);
  events.push_back(ev);
  // I1->I2
  ev.type=epidemiological_event::I1I2;
  cr += treemap[root].I1 * rates.gamma1;
  cumrate.push_back(cr);
  events.push_back(ev);
  // I2->R
  ev.type=epidemiological_event::I2R;
  cr += treemap[root].I2 * rates.gamma2;
  cumrate.push_back(cr);
  events.push_back(ev);

  total_rate=cumrate.back();
}

void SEIRPopulation::apply_event(int evn)
{
  epidemiological_event     &ev=events[evn];
  node_data &noded=treemap[ev.node];
  node_t    l1node;
  int       noden;
  auto      listi=listE1.begin();

  switch(ev.type) {
  case epidemiological_event::SE1:
    noden=ran(noded.S);                // Choose a susceptible at random within the level
    listi = listS.begin() + noded.first_S_in_list + noden;
    l1node=*listi;                     // update lists
    listE1.push_back(l1node);
    listS.erase(listi);
    update_counts<readS,readE1>(l1node);
    update_after_erase_susceptible(l1node);
    break;
    
  case epidemiological_event::E1E2:    // from here on, noded.node must be root
    noden=ran(noded.E1);              // randomly choose a level 1 node (family)
    listi = listE1.begin()+noden;
    l1node=*listi;
    listE2.push_back(l1node);         // update lists
    listE1.erase(listi);
    update_counts<readE1,readE2>(l1node);
    break;

  case epidemiological_event::E2I1:
    noden=ran(noded.E2);              // randomly choose a level 0 node
    listi = listE2.begin()+noden;
    l1node=*listi;
    listI1.push_back(l1node);         // update lists
    listE2.erase(listi);
    update_counts<readE2,readI1>(l1node);
    break;

  case epidemiological_event::I1I2:
    noden=ran(noded.I1);              // randomly choose a level 0 node
    listi = listI1.begin()+noden;
    l1node=*listi;
    listI2.push_back(l1node);         // update lists
    listI1.erase(listi);
    update_counts<readI1,readI2>(l1node);
    break;
    
  case epidemiological_event::I2R:
    noden=ran(noded.I2);              // randomly choose a level 0 node
    listi = listI2.begin()+noden;
    l1node=*listi;
    listI2.erase(listi);              // update lists
    update_counts<readI2,readR>(l1node);
    break;
    
  }
}

template <typename readF1,typename readF2>
void SEIRPopulation::update_counts(node_t cnode)
{
  lemon::ListDigraph::InArcIt arc(tree,cnode);
  do {
    node_data& cnoded=treemap[cnode];
    ( readF1::field(cnoded) )--;
    ( readF2::field(cnoded) )++;

    if (std::is_same<readF2,SEIRPopulation::readI1>::value) {
      if (cnoded.I1+cnoded.I2 == 1) {
	cnoded.infected_nodes_in_list=infected_nodes.size();
	infected_nodes.push_back(cnode);
      }
    }
    if (std::is_same<readF1,SEIRPopulation::readI2>::value) {
      if (cnoded.I1+cnoded.I2 == 0) {
	auto it = infected_nodes.begin() + cnoded.infected_nodes_in_list;
	cnoded.infected_nodes_in_list=-1;
	for (auto it2=it+1; it2!=infected_nodes.end(); ++it2)
	  treemap[*it2].infected_nodes_in_list--;
	infected_nodes.erase(it);
      }
    }

    arc=lemon::ListDigraph::InArcIt(tree,cnode);
  } while ( arc != lemon::INVALID && (cnode=tree.source(arc)) != lemon::INVALID ) ;
}

void SEIRPopulation::update_after_erase_susceptible(node_t node)
{
  node_data &noded=treemap[node];
  assert(noded.level==1);
  auto arc = graph_t::InArcIt(tree,node);

  do {
    node_data &noded=treemap[node];
    int level=noded.level;
    for (int in=noded.level_nodes_in_list+1; in<level_nodes[level].size(); ++in) {
      node_data &niterd=treemap[level_nodes[level][in]];
      niterd.first_S_in_list--;
    }
    arc=graph_t::InArcIt(tree,node);
  }  while (arc!=lemon::INVALID &&  (node=tree.source(arc)) != lemon::INVALID );
}

void SEIRPopulation::add_imported(int I)
{
  I-=gdata.infections_imported;  // This is the number of new cases
  if (I<0)
    {std::cerr << "Error in imported infections file: external infections must be monotonically increasing\n"; exit(1);}

  node_data& rootd=treemap[root];
  if (I>rootd.S)
    {std::cerr << "Cannot add imported, too few suscetibles\n"; exit(1);}

  // Randomly choose and infect I individuals
  for (int infn=0; infn<I; ++infn) {
    int noden=ran(rootd.S);
    // find in family and infect in state I1
    auto listi= listS.begin() + noden;
    node_t l1node=*listi;
    listI1.push_back(l1node);
    listS.erase(listi);
    update_counts<readS,readI1>(l1node);
    update_after_erase_susceptible(l1node);
  }
  gdata.infections_imported+=I;
}



///////////////////////////////////////////////////////////////////////////////
//
// simulation driver: uses a given SEIRpopulation object to drive the
// dynamics (Gillespie).  Output through a SEEIIRstate object

void run(SEIRPopulation &pop,SEEIIRstate *state)
{
  Exponential_distribution rexp;
  Uniform_real ran(0,1.);
  double deltat,time=0;
  double last=-10;

  event_queue_t events=event_queue;


  // while (!events.empty()) {
  //   std::cerr  << "time type number " << events.front().time << " " << events.front().kind
  // 	       << " " << events.front().enumber << '\n';
  //   events.pop();
  // }

  // exit(1);


  SEEIIRistate gstate;

  while (time<options.steps) {

    // compute transition probabilities
    pop.compute_rates();
    double mutot=pop.total_rate;

    // advance time
    deltat=rexp(1./mutot);
    time+=deltat;

    if (time>=events.front().time) {              // imported infections or beta change

      if (events.size()==1) break;
      time=events.front().time;

      switch (events.front().kind) {
      case event::infection:
  	pop.add_imported(options.imported_infections[events.front().enumber].I);
  	break;
      case event::rate_change:
  	pop.set_rate_parameters(options.rates_vs_time[events.front().enumber]);
  	break;
      }
      events.pop();

    } else {

      // choose the transition and apply it
      double r=ran()*mutot;
      int e=bsearch(r,pop.cumrate);
      pop.apply_event(e);
	
    }

    if (time>=last+1.) {  // print or accumulate average
      last=time;
      node_data &rootd=pop.treemap[pop.root];
      gstate.N=rootd.N;
      gstate.S=rootd.S;
      gstate.E1=rootd.E1;
      gstate.E2=rootd.E2;
      gstate.I1=rootd.I1;
      gstate.I2=rootd.I2;
      gstate.R=rootd.R;
      gstate.inf_imported=pop.gdata.infections_imported;
      gstate.inf_close=pop.gdata.infections_level[1];
      gstate.inf_community=pop.gdata.infections_level[2];
      gstate.beta_out=pop.rates.beta[2];
      state->push(time,gstate);
    }

  }
} 



///////////////////////////////////////////////////////////////////////////////
//
// main and noffspring
//
// noffspring provides the number of descendants at each tree level

Discrete_distribution *Mdist;

int noffspring(int level)
{
  switch(level) {
  case 2: return options.Nfamilies;
  case 1: return (*Mdist)();
  default: return 0;
  }
}

int main(int argc,char *argv[])
{
  read_parameters(argc,argv);
  Random_number_generator RNG(options.seed);
  Mdist = new Discrete_distribution(options.Mmax+1,options.PM);
  merge_events();

  // Prepare global state (for output) and population
  SEEIIRstate *state;
  state = options.Nruns>1 ?
          new SEEIIRstate_av : new SEEIIRstate;
  std::cout << state->header() << '\n';

  SEIRPopulation pop(options.levels,noffspring);

  // pop.check_structures();
  // return 1;

  // Do runs and print results
  for (int n=0; n<options.Nruns; ++n) {
    // std::cout << "# N = " << pop.gstate.N << '\n';
    run(pop,state);
    pop.set_all_S();
  }

  if (options.Nruns>1)
    std::cout << *state;

  delete state;
  delete Mdist;
}
