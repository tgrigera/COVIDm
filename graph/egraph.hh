/*
 * egraph.hh -- classes for epidemiological models on graphs (header)
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

#ifndef EGRAPH_HH
#define EGRAPH_HH

//
// Print source context (for debugging)
#define I_AM_HERE \
  std::cerr << "(DD) Reached " << __FILE__ << ":" << __LINE__ << " (in function " \
  <<	 __func__ << ")\n";

#include <lemon/full_graph.h>
#include <lemon/grid_graph.h>
#include <lemon/maps.h>

#include "../qdrandom.hh"

template <typename Graph>
class Graph_base {
public:
  typedef          Graph graph_t;
  typedef typename Graph::Node  node_t;

  graph_t                 &graph;
  typename graph_t::template ArcMap<double> arc_weight;
  int                     node_count;

  node_t node(int id);
  node_t random_node();
  int    id(typename graph_t::Node);

protected:
  Graph_base(Graph* graphp,double default_arc_weight=1.);
  ~Graph_base() {delete graphp;}
  lemon::IdMap<graph_t,typename graph_t::Node>     idmap;
  lemon::RangeIdMap<graph_t,typename graph_t::Node>  rangemap;

  graph_t         *graphp;
  Uniform_integer ran;
} ;

template <typename Graph>
Graph_base<Graph>::Graph_base(Graph* graphp,double default_arc_weight) :
  graphp(graphp),
  graph(*graphp),
  arc_weight(graph,default_arc_weight),
  idmap(graph),
  rangemap(graph)
{
  node_count=lemon::countNodes(graph);
}

template <typename Graph>
inline typename Graph_base<Graph>::node_t Graph_base<Graph>::node(int id)
{
  return idmap(id);
}

template <typename Graph>
inline int Graph_base<Graph>::id(Graph_base<Graph>::node_t node)
{
  return idmap[node];
}

template <typename Graph>
inline typename Graph_base<Graph>::node_t Graph_base<Graph>::random_node()
{
  return rangemap(ran(rangemap.size()));
}


///////////////////////////////////////////////////////////////////////////////
//
// Fully-connected graph

class FCGraph : public Graph_base<lemon::FullGraph>  { 
public:

  static FCGraph* create(int N);

private:
  FCGraph(graph_t* g,double def_arc_w);
} ;

inline FCGraph::FCGraph(graph_t* g,double aw) :
  Graph_base<lemon::FullGraph>(g,aw)
{}

inline FCGraph* FCGraph::create(int N)
{
  graph_t *g = new graph_t(N);
  return new FCGraph(g,1./N);
}
  

// class SQGraph {
// public:
//   typedef lemon::GridGraph graph_t;
//   typedef lemon::GridGraph::Node  node_t;

//   SQGraph(int Lx,int Ly);
//   node_t node(int id);
//   node_t random_node();
//   int   id(graph_t::Node);
  
//   graph_t                  graph;
//   graph_t::ArcMap<double>  arc_weight;

//   int node_count;

// private:
//   lemon::IdMap<graph_t,graph_t::Node>     idmap;
//   lemon::RangeIdMap<graph_t,graph_t::Node>  rangemap;

//   Uniform_integer ran;
// } ;

// inline SQGraph::SQGraph(int Lx,int Ly) :
//   graph(Lx,Ly),
//   arc_weights(graph,1.),
//   idmap(graph),
//   rangemap(graph)
// {
//   node_count=lemon::countNodes(graph);
// }

// }


#endif /* EGRAPH_HH */
