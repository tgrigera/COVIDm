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

#include <lemon/grid_graph.h>
#include <lemon/maps.h>

#include "../qdrandom.hh"

class SQGraph {
public:
  typedef lemon::GridGraph graph_t;
  typedef lemon::GridGraph::Node  node_t;

  SQGraph(int Lx,int Ly);
  node_t node(int id);
  node_t random_node();
  int   id(graph_t::Node);
  
  graph_t                  graph;
  graph_t::ArcMap<double>  arc_weights;

  int node_count;

private:
  lemon::IdMap<graph_t,graph_t::Node>     idmap;
  lemon::RangeIdMap<graph_t,graph_t::Node>  rangemap;

  Uniform_integer ran;
} ;

inline SQGraph::SQGraph(int Lx,int Ly) :
  graph(Lx,Ly),
  arc_weights(graph,1.),
  idmap(graph),
  rangemap(graph)
{
  node_count=lemon::countNodes(graph);
}

inline SQGraph::node_t SQGraph::node(int id)
{
  return idmap(id);
}

inline int SQGraph::id(SQGraph::node_t node)
{
  return idmap[node];
}

inline SQGraph::graph_t::Node SQGraph::random_node()
{
  return rangemap(ran(rangemap.size()));
}


#endif /* EGRAPH_HH */
