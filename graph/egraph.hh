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
  typedef typename graph_t::Arc arc_t;

  graph_t &graph;
  double  arc_weight(arc_t arc);
  node_t  node(int id);
  node_t  random_node();
  int     id(typename graph_t::Node);

  int     node_count;

protected:
  Graph_base(Graph* graphp,double default_arc_weight=1.);
  ~Graph_base() {delete graphp;}

  graph_t                                           *graphp;
  lemon::IdMap<graph_t,typename graph_t::Node>      idmap;
  lemon::RangeIdMap<graph_t,typename graph_t::Node> rangemap;
  Uniform_integer                                   ran;
  double                                            default_arc_weight;
  typename graph_t::template ArcMap<double>         arcmap;

} ;

template <typename Graph>
Graph_base<Graph>::Graph_base(Graph* graphp,double default_arc_weight) :
  graph(*graphp),
  graphp(graphp),
  idmap(graph),
  rangemap(graph),
  default_arc_weight(default_arc_weight),
  arcmap(graph,default_arc_weight)
{
  node_count=lemon::countNodes(graph);
}

template <typename Graph>
inline double Graph_base<Graph>::arc_weight(Graph_base<Graph>::arc_t arc)
{
  return arcmap[arc];
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

///////////////////////////////////////////////////////////////////////////////
//
// Square Lattice

class SQGraph : public Graph_base<lemon::GridGraph> {
public:
  static SQGraph* create(int Lx,int Ly);

private:
  SQGraph(graph_t *g,double def_arc_w);

} ;

inline SQGraph::SQGraph(graph_t* g,double aw) :
  Graph_base<lemon::GridGraph>(g,aw)
{}

inline SQGraph* SQGraph::create(int Lx,int Ly)
{
  graph_t* g=new graph_t(Lx,Ly);
  return new SQGraph(g,1.);
}

#endif /* EGRAPH_HH */
