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

#define NDEBUG

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
#include "gillespie_sampler.hh"
#include "avevar.hh"

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
    time(0.),
    beta(levels+1,0.),
    sigma1(0), sigma2(0), gamma1(0), gamma2(0)
  {}
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

  int               levels;       // tree depth (not counting individuals)
  std::vector<int>  M;            // number of descendants at each level (negative=fluctuating)
  std::vector<std::vector<double>> PM; // Distribution of descendants
  std::string eifile;             // File to read imported infected cases

  // imported infections
  struct ei {double time; int I; int R;} ;
  typedef std::vector<ei>               imported_infections_t;
  imported_infections_t                 imported_infections;
  // rates vs time
  typedef std::vector<rates_t>          rates_vs_time_t;
  rates_vs_time_t                       rates_vs_time;

  char *dfile;        // detailed level info file
  int  detail_level;  // print detail info down to level, negative means don't print

  opt() : last_arg_read(0), detail_level(-1) {}

} options;

/*
 * Read parameters from command-line and file, and compute
 * derived parameters
 *
 */
#include "read_arg.hh"

static int nargs=6;

void show_usage(char *prog)
{
  std::cerr << "usage: " << prog << " parameterfile seed steps Nruns\n\n"
	    << "    or " << prog << " parameterfile seed steps Nruns detail_level detail_file\n\n"
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
  if (argc!=nargs+1 && argc!=nargs-1) show_usage(argv[0]);
  read_arg(argv,options.ifile);
  read_arg(argv,options.seed);
  read_arg(argv,options.steps);
  read_arg(argv,options.Nruns);
  if (argc==nargs+1) {
    read_arg(argv,options.detail_level);
    read_arg(argv,options.dfile);
  }

  FILE *f=fopen(options.ifile,"r");
  if (f==0) throw std::runtime_error(strerror(errno));

  char *buf=readbuf(f);
  sscanf(buf,"%d",&options.levels);
  printf("##### Parameters\n");
  printf("# Nlevels = %d\n",options.levels);

  options.M.resize(options.levels+1);
  options.PM.resize(options.levels+1);
  for (int lev=options.levels; lev>0; --lev) {
    buf=readbuf(f);
    sscanf(buf,"%d",&(options.M[lev]));
    if (options.M[lev]<0) {
      options.PM[lev].resize(-options.M[lev]+1);
      options.PM[lev][0]=0.;
      for (int M=1; M<=-options.M[lev]; ++M) {
	buf=readbuf(f);
	sscanf(buf,"%lg",&(options.PM[lev][M]));
      }
    }
  }

  for (int lev=options.levels; lev>0; --lev) {
    printf("# Number of descendants at level %d = ",lev);
    if (options.M[lev]>0) printf("%d\n",options.M[lev]);
    else {
      printf(" 1 to %d, with weights: \n",-options.M[lev]);
      for (int i=1; i<options.PM[lev].size(); ++i)
	printf("#       %d:   %g\n",i,options.PM[lev][i]);
    }
  }

  printf("#\n# Nruns = %d\n",options.Nruns);
  if (options.detail_level>0)
    printf("# Writing detail down to level %d to file %s\n",options.detail_level,options.dfile);

  buf=readbuf(f);
  options.eifile=buf;
  options.eifile.erase(options.eifile.end()-1);   // remove trailing newline
  read_imported_infections();

  printf("# Imported infections:\n");
  printf("# Time   Imported_inf  Forced_R\n");
  for (auto iir: options.imported_infections)
    printf("# %g %d %d\n",iir.time,iir.I,iir.R);

  read_rates_vs_time(f);
  fclose(f);

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
  int nread;
  
  FILE *f=fopen(options.eifile.c_str(),"r");
  if (f==0) {
    std::cerr << "Error opening file (" << options.eifile << ")\n";
    throw std::runtime_error(strerror(errno));
  }

  opt::ei ei;
  while (ungetc(fgetc(f),f)!=EOF) {
    char *buf=readbuf(f);
    if ( (nread=sscanf(buf,"%lg  %d  %d",&ei.time,&ei.I,&ei.R))!=3) {
      std::cerr  << "couldn't read record: " << buf << "\n";
      std::cerr << "shit\n";
      std::cerr << "correctly read " << nread << " values\n";
      std::cerr << ei.time << ' ' << ei.I << ' ' << ei.R << '\n';
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
  int              forcibly_recovered;
  int              Eacc;
  std::vector<int> infections_level;

  global_data(int levels) :
    infections_imported(0),
    forcibly_recovered(0),
    infections_level(levels+1,0)
  {}
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
    	    << "         S, E1, E2, I1, I2, R "
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
  void force_infection_recover(int I,int R);
  void add_imported(int I);
  void force_recover(int R);
  void unrecover(int S);
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
  std::vector<node_t>                    listS,listE1,listE2,listI1,listI2,listR;
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
  void count_infection_kind(node_t node);

  friend class SEEIIR_observer;

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
  gdata.forcibly_recovered=0;
  gdata.infections_level.resize(levels+1,0);
  gdata.Eacc=0;
  std::fill(gdata.infections_level.begin(),gdata.infections_level.end(),0.);
}

// This recomputes all cumulative counts and rebuilds lists
// except listR
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
    if (noded.I1+noded.I2>0) {
      noded.infected_nodes_in_list=infected_nodes.size();
      infected_nodes.push_back(node);
    }
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
      if (noded.I1+noded.I2>0) {
	noded.infected_nodes_in_list=infected_nodes.size();
	infected_nodes.push_back(node);
      }
      Nprev+=noded.S;
    }
  }
}

void SEIRPopulation::check_structures()
{
  node_data& rootd=treemap[root];

  std::cerr << "Checking tree\n";
  
  // for (int l=levels; l>0; --l) {
  //   std::cout << "***** Level " << l << '\n';
  //   for (auto nn: level_nodes[l]) {
  //     std::cout << treemap[nn];
  //     std::cout << "      firstS " << treemap[nn].first_S_in_list << '\n';
  //   }
  // }

  assert(listS.size()==rootd.S);
  for (int is=0; is<listS.size(); ++is) {
    node_t node=listS[is];
    assert(treemap[node].level==1);
    assert(treemap[node].S>0);
    assert(treemap[node].first_S_in_list<=is);
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

  assert(listR.size()==gdata.forcibly_recovered);

  for (int l=levels; l>0; --l) {
    for (auto nn: level_nodes[l]) {
      auto nnd=treemap[nn];
      if (nnd.S>0) {
	// std::cerr << nnd;
	// std::cerr << "asserting,  " << treemap[listS[nnd.first_S_in_list]].first_S_in_list
	// 	  <<  ">=" << nnd.first_S_in_list << '\n';
	assert(treemap[listS[nnd.first_S_in_list]].first_S_in_list>=nnd.first_S_in_list);
      }
    }
  }

  for (int l=levels; l>0; --l) {
    for (auto nn: level_nodes[l]) {
      auto nnd=treemap[nn];
      if (nnd.I1 + nnd.I2>0) {
	assert(infected_nodes[nnd.infectod_nodes_in_list]==nn);
      }
    }
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
    gdata.Eacc++;
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
    count_infection_kind(l1node);
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

// new infection in node, count kind
void SEIRPopulation::count_infection_kind(node_t node)
{
  auto arc = graph_t::InArcIt(tree,node);
  do {
    node_data &noded=treemap[node];
    int level=noded.level;
    if (noded.I1+noded.I2+noded.R>1 || level==levels) {         // counts have already been updated, so there must be at least one infected
      gdata.infections_level[level]++;
      break;
    }
    arc=graph_t::InArcIt(tree,node);
  }  while (arc!=lemon::INVALID &&  (node=tree.source(arc)) != lemon::INVALID );
}


void SEIRPopulation::force_infection_recover(int I,int R)
{
  if (I!=gdata.infections_imported)
    add_imported(I);
  if (R>gdata.forcibly_recovered)
    force_recover(R-gdata.forcibly_recovered);
  else if (R<gdata.forcibly_recovered)
    unrecover(gdata.forcibly_recovered-R);
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

#ifdef FORCE_RECOVER_WHOLE_FAMILIES

void SEIRPopulation::force_recover(int R)  // Move aprox R individuals from S to R, recovering whole families
{
  node_data& rootd=treemap[root];
  if (R>rootd.S)
    {std::cerr << "Cannot recover, too few suscetibles\n"; exit(1);}

  // Randomly choose and recover aprox R individuals
  int infn=0;
  while (infn<R) {
    int noden=ran(rootd.S);
    // find in family and recover
    auto listi= listS.begin() + noden;
    node_t l1node=*listi;
    node_data& nd=treemap[l1node];
    if (nd.S+nd.R<nd.N) continue;  // Look for a family with only S and R
    int rec=nd.S;
    for (int i=0; i<rec; ++i) {
      listR.push_back(l1node);           // listR tracks only the focibly recovered
      listi = listS.begin() + nd.first_S_in_list;
      listS.erase(listi);
      update_counts<readS,readR>(l1node);
      update_after_erase_susceptible(l1node);
    }
    infn+=rec;
  }
  gdata.forcibly_recovered+=infn;
}

void SEIRPopulation::unrecover(int S)  // Make aprox S of the forcibly recovered susceptible again
{
  if (S>listR.size()) 
    {std::cerr << "Requested too many unrecovers\n"; exit(1);}

  int isus=0;
  while (isus<S) {
    int noden=ran(listR.size());
    // find in family 
    auto listi= listR.begin() + noden;
    node_t l1node=*listi;

    // Now make S all of the forced recoveries of the same family
    node_data& nd=treemap[l1node];
    auto efirst = std::remove_if(listR.begin(),listR.end(),[l1node](node_t x) {return x==l1node;} ); // moves all instances of l1node in listR to the end
    int nrec=listR.end()-efirst;
    nd.S+=nrec;
    nd.R-=nrec;
    isus+=nrec;
    listR.erase(efirst, listR.end() );  // erase all instances of l1node in listR
  }
  recompute_counts();
  gdata.forcibly_recovered-=isus;
}

#else /* FORCE_RECOVER_WHOLE_FAMILIES */

void SEIRPopulation::force_recover(int R)  // Move R individuals from S to R
{
  node_data& rootd=treemap[root];
  if (R>rootd.S)
    {std::cerr << "Cannot recover, too few suscetibles\n"; exit(1);}

  // Randomly choose and recover R individuals
  for (int infn=0; infn<R; ++infn) {
    int noden=ran(rootd.S);
    // find in family and recover
    auto listi= listS.begin() + noden;
    node_t l1node=*listi;
    listR.push_back(l1node);           // listR tracks only the focibly recovered, so that the can be turned susceptible afterwards
    listS.erase(listi);
    update_counts<readS,readR>(l1node);
    update_after_erase_susceptible(l1node);
  }
  gdata.forcibly_recovered+=R;
}

void SEIRPopulation::unrecover(int S)  // Make S of the forcibly recovered susceptible again
{
  if (S>listR.size()) 
    {std::cerr << "Requested too many unrecovers\n"; exit(1);}

  for (int isus=0; isus<S; ++isus) {
    int noden=ran(listR.size());
    // find in family and recover
    auto listi= listR.begin() + noden;
    node_t l1node=*listi;
    node_data& nd=treemap[l1node];
    nd.R--;
    nd.S++;
    listR.erase(listi);
  }
  recompute_counts();
  gdata.forcibly_recovered-=S;
}

#endif /* FORCE_RECOVER_WHOLE_FAMILIES */


///////////////////////////////////////////////////////////////////////////////
//
// commputing and printing aggregate data at the different hierarchical levels

class SEEIIR_observer {
public:
  SEEIIR_observer(SEEIIRstate *state,SEIRPopulation& pop,int dlevel,char *dfile);
  ~SEEIIR_observer();

  void push(double time,SEIRPopulation& pop);
  
private:
  SEEIIRstate  *state;
  SEEIIRistate gstate;
  int  dlevel;
  char *file;
  FILE *f;
} ;

SEEIIR_observer::SEEIIR_observer(SEEIIRstate *state,SEIRPopulation& pop,int dlevel,char *dfile) :
    state(state), dlevel(dlevel), file(dfile)
{
  if (dlevel<0) return;
  f=fopen(file,"w");

  int nc=1+2*(pop.levels-1);
  for (int l=pop.levels; l>=dlevel; --l)
    nc+=pop.level_nodes[l].size();

  fprintf(f,"#     ( 1)|");
  for (int i=2; i<=nc; ++i)
    fprintf(f," |     (%2d)|",i);
  fprintf(f,"\n");
  
  fprintf(f,"#           |------ Level %2d -----| ",pop.levels-1);
  for (int i=pop.levels-2; i>0; --i)
    fprintf(f,"|------ Level %2d -----| ",i);
  for (int l=pop.levels; l>=dlevel; --l) {
    int width=11*pop.level_nodes[l].size() + pop.level_nodes[l].size() - 1;
    std::string fill1((width-10)/2,'-');
    std::string fill2(width-10-fill1.size(),'-');
    fprintf(f,"|%sLevel %2d%s| ",fill1.c_str(),l,fill2.c_str());
  }
  fprintf(f,"\n");

  fprintf(f,"#      time ");
  for (int i=0; i<pop.levels-1; ++i)
    fprintf(f,"        ave         var ");
  for (int l=pop.levels; l>=dlevel; --l)
    for (int n=0; n<pop.level_nodes[l].size(); ++n)
      fprintf(f,"   Node %3d ",n);

  fprintf(f,"\n");
}

SEEIIR_observer::~SEEIIR_observer()
{
  if (dlevel>0) fclose(f);
}

void SEEIIR_observer::push(double time,SEIRPopulation& pop)
{
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
  gstate.inf_community=pop.gdata.infections_level[pop.levels];
  gstate.beta_out=pop.rates.beta[2];
  gstate.Eacc=pop.gdata.Eacc;
  gstate.tinf= 1./pop.rates.gamma1 + 1./pop.rates.gamma2;
  state->push(time,gstate);

  if (dlevel<0) return;

  fprintf(f,"%11.6g ",time);
  
  AveVar<false> av;
  for (int l=pop.levels-1; l>0; --l) {
    av.clear();
    for (node_t node: pop.level_nodes[l]) {
      node_data &noded=pop.treemap[node];
      av.push(noded.I1+noded.I2);
    }
    fprintf(f,"%11.6g %11.6g ",av.ave(),av.var());
  }

  for (int l=pop.levels; l>=dlevel; --l) {
    for (node_t node: pop.level_nodes[l]) {
      node_data &noded=pop.treemap[node];
      fprintf(f,"%11d ",noded.I1+noded.I2);
    }
  }

  fprintf(f,"\n");
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
  SEEIIR_observer observer(state,pop,options.detail_level,options.dfile);
  Gillespie_sampler<SEEIIR_observer,SEIRPopulation> gsamp(observer,0.,options.steps,1.);
  gsamp.push_data(pop);

  while (time<=options.steps) {

    // compute transition probabilities
    pop.compute_rates();
    double mutot=pop.total_rate;
    // advance time
    deltat=rexp(1./mutot);
    time+=deltat;

    if (time>=events.front().time) {              // imported infections or beta change

      time=events.front().time;
      gsamp.push_time(time);
      if (events.size()==1) break;

      switch (events.front().kind) {
      case event::infection:
  	pop.force_infection_recover(options.imported_infections[events.front().enumber].I,options.imported_infections[events.front().enumber].R);
  	break;
      case event::rate_change:
  	pop.set_rate_parameters(options.rates_vs_time[events.front().enumber]);
  	break;
      }
      events.pop();

    } else {

      gsamp.push_time(time);
      // choose the transition and apply it
      double r=ran()*mutot;
      int e=bsearch(r,pop.cumrate);
      pop.apply_event(e);
	
    }

    gsamp.push_data(pop);

  }
} 

///////////////////////////////////////////////////////////////////////////////
//
// main and noffspring
//
// noffspring provides the number of descendants at each tree level

std::vector<Discrete_distribution*> Mdist;

void prepare_noffspring()
{
  Mdist.resize(options.levels);
  for (int l=options.levels; l>0; --l) {
    if (options.M[l]<0) 
      Mdist[l]=new Discrete_distribution(options.PM[l].size(),&(options.PM[l][0]));
  }
}

int noffspring(int level)
{
  if (options.M[level]>0) return options.M[level];
  return (*(Mdist[level]))();
}

int main(int argc,char *argv[])
{
  read_parameters(argc,argv);
  Random_number_generator RNG(options.seed);
  merge_events();

  // Prepare global state (for output) and population
  SEEIIRstate *state;
  state = options.Nruns>1 ?
          new SEEIIRstate_av : new SEEIIRstate;
  std::cout << state->header() << '\n';

  prepare_noffspring();
  SEIRPopulation pop(options.levels,noffspring);

  // pop.check_structures();
  // return 1;

  // Do runs and print results
  for (int n=0; n<options.Nruns; ++n) {
    // std::cout << "# N = " << pop.gstate.N << '\n';
    run(pop,state);
    // pop.check_structures();
    pop.set_all_S();
  }

  if (options.Nruns>1)
    std::cout << *state;

  delete state;
}
