/*
 * seeiir_h.cc
 *
 * Hierarchical stochastic SEIR with two E and to I states.
 * Population grouped in families, nieghbourhoods, towns, etc..
 * Simulated in continuous time (Gillespie algorithm).
 *
 * This is implementation uses structs and pointers to build the
 * hierarchy tree.
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

  int               levels;       // tree depth (not counting individuals)
  std::vector<int>  M;            // number of descendants at each level (negative=fluctuating)
  std::vector<std::vector<double>> PM; // Distribution of descendants
  std::string eifile;             // File to read imported infected cases

  // imported infections
  struct ei {double time; int I;} ;
  typedef std::vector<ei>               imported_infections_t;
  imported_infections_t                 imported_infections;
  // rates vs time
  typedef std::vector<rates_t>          rates_vs_time_t;
  rates_vs_time_t                       rates_vs_time;

  opt() : last_arg_read(0) {}

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

  buf=readbuf(f);
  options.eifile=buf;
  options.eifile.erase(options.eifile.end()-1);   // remove trailing newline
  read_imported_infections();

  printf("# Imported infections:\n");
  printf("# Time   Cases\n");
  for (auto iir: options.imported_infections)
    printf("# %g %d\n",iir.time,iir.I);

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

struct node_data;
typedef node_data*                    node_t;
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

  node_t  parent;

  node_data() :
    level(0), N(0), M(0),
    S(0), E1(0), E2(0), I1(0), I2(0), R(0),
    first_S_in_list(-1), level_nodes_in_list(-1), infected_nodes_in_list(-1),
    parent(0)
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
  node_t                        node;
} ;

/*
 * class SEIRPopulation
 *
 */
class SEIRPopulation {
public:

  SEIRPopulation(int levels,int (*noffspring)(int));
  ~SEIRPopulation();
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

  node_t                      root;
  global_data                 gdata;
  rates_t                     rates;

private:
  int (*noffspring)(int);
  Uniform_integer                        ran;
  std::vector<node_list_t>               level_nodes;
  std::vector<node_t>                    listS,listE1,listE2,listI1,listI2;
  std::vector<node_t>                    infected_nodes;
  
  node_t build_tree(int level);
  void   tree_clear(node_t tree);

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
} ;

SEIRPopulation::SEIRPopulation(int levels,int (*noffspring)(int) ) :
  levels(levels),
  noffspring(noffspring),
  rates(levels),
  gdata(levels),
  root(0)
{
  rebuild_hierarchy();
}

inline SEIRPopulation::~SEIRPopulation()
{
  tree_clear(root);
}

void SEIRPopulation::tree_clear(node_t tree)
{
  if (tree==0) return;
  for (int lev=levels; lev>0; --lev)
    for (auto &node: level_nodes[lev] ) {
      delete node;
    }
}
  
void SEIRPopulation::rebuild_hierarchy()
{
  tree_clear(root);
  level_nodes.clear();
  level_nodes.resize(levels+1);
  root=build_tree(levels);
  set_all_S();
}

// Recursively build a tree starting at given level
node_t SEIRPopulation::build_tree(int level)
{
  node_t subtree=new node_data;
  subtree->level=level;
  subtree->level_nodes_in_list=level_nodes[level].size();
  level_nodes[level].push_back(subtree);
  int M=noffspring(level);
  subtree->M=M;

  if (level>1) {  // we don't store the leaves
    for (int i=0; i<M; ++i) {
      node_t n=build_tree(level-1);
      n->parent=subtree;
    }
  }

  return subtree;
}

void SEIRPopulation::set_all_S()
{
  for (int lev=levels; lev>0; --lev)
    for (auto &node: level_nodes[lev] ) {
      node->N = node->level==1 ? node->M : 0;
      node->S=node->N;
      node->E1=node->E2=node->I1=node->I2=node->R=0;
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

  for (int level=2; level<=levels; ++level) {
    for (auto &node: level_nodes[level]) {
      node->N=node->S=node->E1=node->E2=node->I1=node->I2=node->R=0;
    }
  }

  for (auto node: level_nodes[1]) {
    node->N=node->M;
    node->first_S_in_list=listS.size();
    for (int i=0; i<node->S; ++i) listS.push_back(node);
    for (int i=0; i<node->E1; ++i) listE1.push_back(node);
    for (int i=0; i<node->E2; ++i) listE2.push_back(node);
    for (int i=0; i<node->I1; ++i) listI1.push_back(node);
    for (int i=0; i<node->I2; ++i) listI2.push_back(node);
    if (node->I1+node->I2>0) infected_nodes.push_back(node);

    node_t pnode=node;
    while (pnode->parent!=0) {
      pnode->parent->N += node->N;
      pnode->parent->S += node->S;
      pnode->parent->E1 += node->E1;
      pnode->parent->E2 += node->E2;
      pnode->parent->I1 += node->I1;
      pnode->parent->I2 += node->I2;
      pnode->parent->R += node->R;
      pnode=pnode->parent;
    }

  }
    
  for (int level=2; level<=levels; ++level) {
    int Nprev=0;
    for (auto &node: level_nodes[level]) {
      node->first_S_in_list=Nprev;
      if (node->I1+node->I2>0) infected_nodes.push_back(node);
      Nprev+=node->S;
    }
  }
}

#ifdef NDEBUG
inline void SEIRPopulation::check_structures() {}
#else

void SEIRPopulation::check_structures()
{
  node_data& rootd=*root;

  // for (int l=levels; l>0; --l) {
  //   std::cout << "***** Level " << l << '\n';
  //   for (auto nn: level_nodes[l]) {
  //     std::cout << *nn;
  //     std::cout << "      firstS " << nn->first_S_in_list << '\n';
  //   }
  // }

  assert(listS.size()==rootd.S);
  for (int is=0; is<listS.size(); ++is) {
    node_t node=listS[is];
    assert(node->level==1);
    assert(node->S>0);
    assert(node->first_S_in_list<=is);
  }
  assert(listE1.size()==rootd.E1);
  for (node_t &node: listE1) {
    assert(node->level==1);
    assert(node->E1>0);
  }
  assert(listE2.size()==rootd.E2);
  for (node_t &node: listE2) {
    assert(node->level==1);
    assert(node->E2>0);
  }
  assert(listI1.size()==rootd.I1);
  for (node_t &node: listI1) {
    assert(node->level==1);
    assert(node->I1>0);
  }
  assert(listI2.size()==rootd.I2);
  for (node_t &node: listI2) {
    assert(node->level==1);
    assert(node->I2>0);
  }
  for (node_t &node: infected_nodes) {
    assert(node->I1+node->I2>0);
  }

  for (int l=levels; l>0; --l) {
    for (auto nn: level_nodes[l]) {
      if (nn->S>0) {
  	// std::cerr << nnd;
  	// std::cerr << "asserting,  " << treemap[listS[nnd.first_S_in_list]].first_S_in_list
  	// 	  <<  ">=" << nnd.first_S_in_list << '\n';
  	assert(listS[nn->first_S_in_list]->first_S_in_list>=nn->first_S_in_list);
      }
    }
  }
}

#endif // NDEBUG


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
  for (node_t node: infected_nodes) {
    ev.node=node;
    double norm = node->level >1 ? 1./(node->N-1) : 1;
    cr += node->S * rates.beta[node->level] * (node->I1 + node->I2) * norm;
    cumrate.push_back(cr);
    events.push_back(ev);
  }

  // the other events are only global
  ev.node=root;
  // E1->E2
  ev.type=epidemiological_event::E1E2;
  cr += root->E1 * rates.sigma1;
  cumrate.push_back(cr);
  events.push_back(ev);
  // E2->I1
  ev.type=epidemiological_event::E2I1;
  cr += root->E2 * rates.sigma2;
  cumrate.push_back(cr);
  events.push_back(ev);
  // I1->I2
  ev.type=epidemiological_event::I1I2;
  cr += root->I1 * rates.gamma1;
  cumrate.push_back(cr);
  events.push_back(ev);
  // I2->R
  ev.type=epidemiological_event::I2R;
  cr += root->I2 * rates.gamma2;
  cumrate.push_back(cr);
  events.push_back(ev);

  total_rate=cumrate.back();
}

void SEIRPopulation::apply_event(int evn)
{
  epidemiological_event     &ev=events[evn];
  node_data &noded=*(ev.node);
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
  do {
    ( readF1::field(*cnode) )--;
    ( readF2::field(*cnode) )++;

    if (std::is_same<readF2,SEIRPopulation::readI1>::value) {
      if (cnode->I1+cnode->I2 == 1) {  // first infection, add to infected nodes list
	cnode->infected_nodes_in_list=infected_nodes.size();
	infected_nodes.push_back(cnode);
      }
    }
    if (std::is_same<readF1,SEIRPopulation::readI2>::value) {
      if (cnode->I1+cnode->I2 == 0) {  // no more infected, remove from list
	auto it = infected_nodes.begin() + cnode->infected_nodes_in_list;
	cnode->infected_nodes_in_list=-1;
	for (auto it2=it+1; it2!=infected_nodes.end(); ++it2)
	  (*it2)->infected_nodes_in_list--;
	infected_nodes.erase(it);
      }
    }

    cnode=cnode->parent;
  } while ( cnode!=0 );
}

void SEIRPopulation::update_after_erase_susceptible(node_t node)
{
  node_data &noded=*node;
  assert(noded.level==1);

  do {
    int level=node->level;
    for (int in=node->level_nodes_in_list+1; in<level_nodes[level].size(); ++in) {
      node_data &niterd=*(level_nodes[level][in]);
      niterd.first_S_in_list--;
    }
    node=node->parent;
  } while (node!=0);
}

// new infection in node, count kind
void SEIRPopulation::count_infection_kind(node_t node)
{
  do {
    if (node->I1+node->I2>0) {
      gdata.infections_level[node->level]++;
      break;
    }
    node=node->parent;
  } while (node!=0);
}

void SEIRPopulation::add_imported(int I)
{
  I-=gdata.infections_imported;  // This is the number of new cases
  if (I<0)
    {std::cerr << "Error in imported infections file: external infections must be monotonically increasing\n"; exit(1);}

  //node_data& rootd=treemap[root];
  node_data& rootd=*root;
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
      node_data &rootd=*(pop.root);
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
      state->push(time,gstate);
    }

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
    pop.set_all_S();
  }

  if (options.Nruns>1)
    std::cout << *state;

  delete state;
}
