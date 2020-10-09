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
// Base class for most epidemiological graphs
//
// NOTE that the class is written assuming that the graph is static (i.e. links and
// nodes are not added or deleted during the simulation).
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
// There is also a (protected) graph
//
//   - fgraph: The full graph from which the other two are obtained
//             through filtering.  Implementations based on this class will build this as an
//             actual LEMON graph, but others might not; the full graph is a conceptual
//             construction.
//
// The only safe way to obtain the hierarchical parent of an individual node is by calling
// for_each_anode() and related methods.
//
// Node identities:
//
//   nodes can be selected with:
//       -  a node object (type hnode_t or inode_t)
//       -  a node id (lemon's idmap, obtained through id methods), this does not change
//       -  an integer (lemon's rangemap, obtained through rangeid methods), this can change
//          if nodes are added or deleted
//
//  We provide methods id() and inode() to convert from node to ids (for the individual graph)
//  more easily.  Unfortunately the RangeIdMap does not seem to work correctly for the individual
//  graph, which is an adapted graph (lemon::FilterNodes).  We thus provide custom-implemented
//  random_inode() method that returns a randomly chose individual node.

class Graph_base {
public:
  typedef lemon::FilterArcs<digraph_t>  hgraph_t;
  typedef lemon::FilterNodes<digraph_t> igraph_t;
  typedef typename hgraph_t::Node       hnode_t;
  typedef typename igraph_t::Arc        iarc_t;
  typedef typename igraph_t::Node       inode_t;

  igraph_t &igraph;               
  hgraph_t &hgraph;
  virtual  double arc_weight(iarc_t arc);
  inode_t  random_inode();
  int      id(inode_t);            // returns Lemon's unique id 
  inode_t  inode(int id);          // returns a node referenced by id (obtained through id() method)

  template <typename Fun>
  void     for_each_anode(inode_t,Fun fun);

  hnode_t   hroot;
  int       inode_count;

protected:
  typedef digraph_t                     fgraph_t;
  typedef typename fgraph_t::Node       node_t;

  Graph_base(fgraph_t* fgraphp,igraph_t* igraphp,hgraph_t* hgraphp,
	     node_t hroot,double default_arc_weight=1.);
  virtual ~Graph_base();

  fgraph_t                                   *fgraphp;
  igraph_t                                   *igraphp;
  hgraph_t                                   *hgraphp;
  fgraph_t                                   &fgraph;

  Uniform_integer                              ran;
  double                                       default_arc_weight;
  std::vector<inode_t>                         rangemap;  // must build a custom one because Lemon's does not work with graph adaptors

} ;

inline Graph_base::Graph_base(fgraph_t* fgraphp,igraph_t* igraphp,hgraph_t* hgraphp,
			      node_t hroot,double default_arc_weight) :
  fgraphp(fgraphp), igraphp(igraphp), hgraphp(hgraphp),
  fgraph(*fgraphp), igraph(*igraphp), hgraph(*hgraphp),
  hroot(hroot),
  default_arc_weight(default_arc_weight)
{
  inode_count=lemon::countNodes(igraph);
  rangemap.reserve(inode_count);
  for (igraph_t::NodeIt n(igraph); n!=lemon::INVALID; ++n)
    rangemap.push_back(n);
}

inline int Graph_base::id(Graph_base::inode_t node)
{
  return igraph.id(node);
}

inline Graph_base::inode_t Graph_base::inode(int id)
{
  return igraph.nodeFromId(id);
}

inline Graph_base::inode_t Graph_base::random_inode()
{
  return rangemap[ran(rangemap.size())];
}

// Applies a function to all aggregate (hierarchichal) nodes starting at inode.
// assumes each node has only one parent (incoming arc) of higher hierarchy
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
// Fully-connected graph and weighted fully-connected graph

// class FCGraph : public Graph_base  { 
// public:
//   static FCGraph* create(int N);
//   virtual ~FCGraph();

// protected:
//   struct FCGraph_ctor_data {
//     fgraph_t* fg;
//     igraph_t* ig;
//     hgraph_t *hg;
//     node_t hroot;
//     double def_arc_w;
//     fgraph_t::ArcMap<bool> *hmap;
//     fgraph_t::NodeMap<bool> *imap;
//   } ;

//   FCGraph(FCGraph_ctor_data *cdata);
//   static FCGraph_ctor_data* create_ctor_data(int N);

// private:
//   fgraph_t::ArcMap<bool>  *hmap;   // hierarchical graph has same nodes, less arcs
//   fgraph_t::NodeMap<bool> *imap;  // individual graph has only the individual nodes
  
// } ;

// inline FCGraph::FCGraph(FCGraph_ctor_data *cdata) :
//   Graph_base(cdata->fg,cdata->ig,cdata->hg,cdata->hroot,cdata->def_arc_w),
//   hmap(cdata->hmap), imap(cdata->imap)
// {
//   delete cdata;
// }

// inline FCGraph* FCGraph::create(int N)
// {
//   return new FCGraph(create_ctor_data(N));
// }


class FCGraph { 
public:
  typedef digraph_t                     hgraph_t;
  typedef typename hgraph_t::Node       hnode_t;
  typedef lemon::FullGraph              igraph_t;
  typedef typename igraph_t::Arc        iarc_t;
  typedef typename igraph_t::Node       inode_t;

  igraph_t &igraph;               
  hgraph_t &hgraph;
  virtual  double arc_weight(iarc_t arc);
  virtual  double arc_weight(inode_t i,inode_t j);
  inode_t  random_inode();
  int      id(inode_t);            // returns Lemon's unique id 
  inode_t  inode(int id);          // returns a node referenced by id (obtained through id() method)

  template <typename Fun>
  void     for_each_anode(inode_t,Fun fun);

  hnode_t   hroot;
  int       inode_count;

  static FCGraph* create(int N);
  virtual ~FCGraph();

protected:
  struct FCGraph_ctor_data {
    igraph_t* ig;
    hgraph_t *hg;
    hnode_t   hroot;
    double    def_arc_w;
  } ;

  igraph_t *igraphp;
  hgraph_t *hgraphp;
  double   default_arc_weight;
  Uniform_integer                              ran;
  lemon::RangeIdMap<igraph_t,igraph_t::Node>   rangemap;
  //std::vector<inode_t>                         rangemap;
  
  FCGraph(FCGraph_ctor_data *cdata);
  static FCGraph_ctor_data* create_ctor_data(int N);

private:
  FCGraph_ctor_data *cdata;
} ;

inline FCGraph::FCGraph(FCGraph_ctor_data *cdata) :
  igraph(*(cdata->ig)),
  hgraph(*(cdata->hg)),
  igraphp(cdata->ig),
  hgraphp(cdata->hg),
  hroot(cdata->hroot),
  default_arc_weight(cdata->def_arc_w),
  rangemap(*(cdata->ig))
{
  delete cdata;
  inode_count=lemon::countNodes(igraph);
  // rangemap.reserve(inode_count);
  // for (igraph_t::NodeIt n(igraph); n!=lemon::INVALID; ++n)
  //   rangemap.push_back(n);
}

inline FCGraph* FCGraph::create(int N)
{
  return new FCGraph(create_ctor_data(N));
}

inline int FCGraph::id(FCGraph::inode_t node)
{
  return igraph.id(node);
}

inline FCGraph::inode_t FCGraph::inode(int id)
{
  return igraph.nodeFromId(id);
}

inline FCGraph::inode_t FCGraph::random_inode()
{
  return rangemap(ran(rangemap.size()));
  //  return rangemap[ran(rangemap.size())];
}

// Applies a function to all aggregate (hierarchichal) nodes starting at inode.
// assumes each node has only one parent (incoming arc) of higher hierarchy
template <typename Fun>
void FCGraph::for_each_anode(inode_t inode,Fun fun)
{
  fun(hroot);
}

//
// FC graph with multiplicative weights
//

class MWFCGraph : public FCGraph {
public:
  static MWFCGraph* create(int N);
  double arc_weight(iarc_t arc);
  double arc_weight(inode_t i,inode_t j);
  void set_weights_random_multiplicative(double (*betadist)(),double scale);

protected:
  MWFCGraph(FCGraph_ctor_data *cdata);

private:
  igraph_t::NodeMap<double> wfactor;
} ;

inline MWFCGraph::MWFCGraph(FCGraph::FCGraph_ctor_data *cdata) :
  FCGraph(cdata),
  wfactor(igraph)
{}

inline MWFCGraph* MWFCGraph::create(int N)
{
  return new MWFCGraph(FCGraph::create_ctor_data(N));
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


///////////////////////////////////////////////////////////////////////////////
//
// Hierarchical fully-connected  graph

// class HFCGraph : public Graph_base {
// public:
//   static HFCGraph* create(int Lx,int Ly);
//   ~HFCGraph();

// private:
//   HFCGraph(fgraph_t* fg,igraph_t* ig,hgraph_t *hg,node_t hroot,double def_arc_w,
// 	   fgraph_t::ArcMap<bool> *hmap, fgraph_t::NodeMap<bool> *imap);

//   fgraph_t::ArcMap<bool>  *hmap;   // hierarchical graph has same nodes, less arcs
//   fgraph_t::NodeMap<bool> *imap;  // individual graph has only the individual nodes
// } ;

// inline HFCGraph::HFCGraph(fgraph_t* fg,igraph_t* ig,hgraph_t *hg,node_t hroot,double def_arc_w,
// 			  fgraph_t::ArcMap<bool> *hmap, fgraph_t::NodeMap<bool> *imap) : 
//   Graph_base(fg,ig,hg,hroot,def_arc_w),
//   hmap(hmap), imap(imap)
// {}


#endif /* EGRAPH_HH */
