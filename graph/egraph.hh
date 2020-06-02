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
#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/maps.h>
#include <lemon/adaptors.h>

#include "../qdrandom.hh"

//typedef lemon::ListDigraph digraph_t;
typedef lemon::SmartDigraph digraph_t;

///////////////////////////////////////////////////////////////////////////////
//
// Graph_base
//
// Base class for all epidemiological graphs
//
// Provides two graphs to the model-implementing classes:
//
//   - igraph: A graph which includes only nodes representing individuals,
//             and (weighted) links (bidirectional arcs) among them representing
//             their interaction.
//
//   - hgraph: A hierarchical graph which includes all nodes (individual and
//             aggreate) but with only (directed) arcs representing the aggregation
//             hierarchy.  No arcs between individual nodes.
//
//   - fgraph: The full graph from which the other two are obtained through filtering

class Graph_base {
public:
  typedef digraph_t                     fgraph_t;
  typedef lemon::FilterArcs<digraph_t>  hgraph_t;
  typedef lemon::FilterNodes<digraph_t> igraph_t;
  typedef typename fgraph_t::Node       node_t;
  typedef typename fgraph_t::Arc        arc_t;

  igraph_t &igraph;               
  hgraph_t &hgraph;
  double   arc_weight(arc_t arc);
  node_t   inode(int id);
  node_t   random_inode();
  int      id(node_t);

  template <typename Fun>
  void     for_each_anode(node_t,Fun fun);

  node_t   hroot;
  int      inode_count;

protected:
  Graph_base(fgraph_t* fgraphp,igraph_t* igraphp,hgraph_t* hgraphp,
	     node_t hroot,double default_arc_weight=1.);
  ~Graph_base();

  fgraph_t                                           *fgraphp;
  igraph_t                                           *igraphp;
  hgraph_t                                           *hgraphp;
  fgraph_t                                           &fgraph;
  lemon::IdMap<igraph_t,igraph_t::Node>     idmap;
  lemon::RangeIdMap<igraph_t,igraph_t::Node> rangemap;
  Uniform_integer                                    ran;
  double                                             default_arc_weight;
  typename igraph_t::template ArcMap<double>         arcmap;

} ;

inline Graph_base::Graph_base(fgraph_t* fgraphp,igraph_t* igraphp,hgraph_t* hgraphp,
			      node_t hroot,double default_arc_weight) :
  fgraphp(fgraphp), igraphp(igraphp), hgraphp(hgraphp),
  fgraph(*fgraphp), igraph(*igraphp), hgraph(*hgraphp),
  hroot(hroot),
  idmap(igraph),
  rangemap(igraph),
  default_arc_weight(default_arc_weight),
  arcmap(igraph,default_arc_weight)
{
  inode_count=lemon::countNodes(igraph);
}

inline Graph_base::~Graph_base()
{
  delete fgraphp;
  delete igraphp;
  delete hgraphp;
}

inline double Graph_base::arc_weight(Graph_base::arc_t arc)
{
  return arcmap[arc];
}

inline Graph_base::node_t Graph_base::inode(int id)
{
  return idmap(id);
}

inline int Graph_base::id(Graph_base::node_t node)
{
  return idmap[node];
}

inline Graph_base::node_t Graph_base::random_inode()
{
  return rangemap(ran(rangemap.size()));
}

template <typename Fun>
void Graph_base::for_each_anode(node_t inode,Fun fun)
{
  hgraph_t::InArcIt arc(hgraph,inode);
  while ( (arc = hgraph_t::InArcIt(hgraph,inode) ) != lemon::INVALID ) {
    inode = hgraph.source(arc);
    fun(inode);
  }    
}

///////////////////////////////////////////////////////////////////////////////
//
// Fully-connected graph

class FCGraph : public Graph_base  { 
public:
  static FCGraph* create(int N);
  ~FCGraph();

private:
  FCGraph(fgraph_t* fg,igraph_t* ig,hgraph_t *hg,node_t hroot,double def_arc_w,
	  fgraph_t::ArcMap<bool> *hmap, fgraph_t::NodeMap<bool> *imap);

  fgraph_t::ArcMap<bool>  *hmap;   // hierarchical graph has same nodes, less arcs
  fgraph_t::NodeMap<bool> *imap;  // individual graph has only the individual nodes
  
} ;

inline FCGraph::FCGraph(fgraph_t* fg,igraph_t* ig,hgraph_t *hg,
			node_t hroot, double def_arc_w,
			fgraph_t::ArcMap<bool> *hmap,
			fgraph_t::NodeMap<bool> *imap) :
  Graph_base(fg,ig,hg,hroot,def_arc_w),
  hmap(hmap), imap(imap)
{}

inline FCGraph::~FCGraph()
{
  delete hmap;
  delete imap;
}

///////////////////////////////////////////////////////////////////////////////
//
// Square Lattice

class SQGraph : public Graph_base {
public:
  static SQGraph* create(int Lx,int Ly);
  ~SQGraph();

private:
  SQGraph(fgraph_t* fg,igraph_t* ig,hgraph_t *hg,node_t hroot,double def_arc_w,
	  fgraph_t::ArcMap<bool> *hmap, fgraph_t::NodeMap<bool> *imap);

  fgraph_t::ArcMap<bool>  *hmap;   // hierarchical graph has same nodes, less arcs
  fgraph_t::NodeMap<bool> *imap;  // individual graph has only the individual nodes
} ;

inline SQGraph::SQGraph(fgraph_t* fg,igraph_t* ig,hgraph_t *hg,node_t hroot,double def_arc_w,
			fgraph_t::ArcMap<bool> *hmap, fgraph_t::NodeMap<bool> *imap) : 
  Graph_base(fg,ig,hg,hroot,def_arc_w),
  hmap(hmap), imap(imap)
{}

#endif /* EGRAPH_HH */
