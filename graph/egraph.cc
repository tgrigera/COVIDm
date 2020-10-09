/*
 * egraph.cc -- classes for epidemiological models on graphs (header)
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

#include <math.h>

#include "egraph.hh"

///////////////////////////////////////////////////////////////////////////////
//
// Graph_base

Graph_base::~Graph_base()
{
  delete igraphp;
  delete hgraphp;
  delete fgraphp;
}

double Graph_base::arc_weight(Graph_base::iarc_t arc)
{
  return default_arc_weight;
}

///////////////////////////////////////////////////////////////////////////////
//
// Fully-connected graph

// FCGraph::FCGraph_ctor_data* FCGraph::create_ctor_data(int N)
// {
//   FCGraph_ctor_data *cdata=new FCGraph_ctor_data;
  
//   // The graph is built copying the structure of a lemon FullGraph
//   lemon::FullGraph fcg(N);
//   cdata->fg = new fgraph_t;
//   lemon::digraphCopy(fcg,*(cdata->fg)).run();

//   // Now add the hierarchical (aggregate) node, and the filtered graphs
//   // for the individual and hierarchical vision of the structure
//   cdata->hmap =   // hierarchical graph has same nodes, less arcs
//     new fgraph_t::ArcMap<bool>(*(cdata->fg),false);
//   cdata->imap =   // individual graph has only the individual nodes
//     new fgraph_t::NodeMap<bool>(*(cdata->fg),true);

//   cdata->hroot=cdata->fg->addNode();
//   (*(cdata->imap))[cdata->hroot]=false;
//   for (fgraph_t::NodeIt inode(*(cdata->fg)); inode!=lemon::INVALID; ++inode) {
//     if (inode==cdata->hroot) continue;
//     auto arc=cdata->fg->addArc(cdata->hroot,inode);
//     (*(cdata->hmap))[arc]=true;
//   }

//   cdata->ig = new igraph_t(*(cdata->fg),*(cdata->imap));
//   cdata->hg = new hgraph_t(*(cdata->fg),*(cdata->hmap));

//   return cdata;
// }

FCGraph::FCGraph_ctor_data* FCGraph::create_ctor_data(int N)
{
  FCGraph_ctor_data *cdata=new FCGraph_ctor_data;
  
  cdata->ig = new igraph_t(N);
  cdata->hg = new hgraph_t;
  cdata->hroot=cdata->hg->addNode();
  cdata->def_arc_w=1;

  return cdata;
}

FCGraph::~FCGraph()
{
  delete igraphp;
  delete hgraphp;
}

double FCGraph::arc_weight(iarc_t arc)
{
  return default_arc_weight;
}

double FCGraph::arc_weight(inode_t i,inode_t j)
{
  return default_arc_weight;
}

// FCGraph::~FCGraph()
// {
//   delete hmap;
//   delete imap;
// }

///////////////////////////////////////////////////////////////////////////////
//
// Weighted FC Graphs

// Multiplicative weights for the FC grah
//
// Assigns symmetrical weights to each arc according to J_ij = beta_i * beta_j / A
// betas are taken from user-provided random distribution betadist() and
// A= inode_count * beta_scale
//

void MWFCGraph::set_weights_random_multiplicative(double (*betadist)(),double beta_scale)
{
  double wnorm=sqrt(beta_scale*inode_count);
  // obtain random betas
  for (igraph_t::NodeIt node(igraph); node!=lemon::INVALID; ++node) {
    wfactor[node]=betadist()/wnorm;
  }
}

double MWFCGraph::arc_weight(iarc_t arc)
{
  inode_t i=igraph.source(arc);
  inode_t j=igraph.target(arc);
  return wfactor[i]*wfactor[j];
}

double MWFCGraph::arc_weight(inode_t i,inode_t j)
{
  return wfactor[i]*wfactor[j];
}

///////////////////////////////////////////////////////////////////////////////
//
// Square Lattice

SQGraph* SQGraph::create(int Lx,int Ly)
{
  // The graph is built copying the structure of a lemon GridGraph
  lemon::GridGraph grg(Lx,Ly);
  fgraph_t *fgr = new fgraph_t;
  lemon::digraphCopy(grg,*fgr).run();

  // Now add the hierarchical (aggregate) node, and the filtered graphs
  // for the individual and hierarchical vision of the structure
  fgraph_t::ArcMap<bool> *hmap =   // hierarchical graph has same nodes, less arcs
    new fgraph_t::ArcMap<bool>(*fgr,false);
  fgraph_t::NodeMap<bool> *imap =   // individual graph has only the individual nodes
    new fgraph_t::NodeMap<bool>(*fgr,true);
  
  fgraph_t::Node hnode=fgr->addNode();
  (*imap)[hnode]=false;
  for (fgraph_t::NodeIt inode(*fgr); inode!=lemon::INVALID; ++inode) {
    if (inode==hnode) continue;
    auto arc=fgr->addArc(hnode,inode);
    (*hmap)[arc]=true;
  }
  igraph_t *igr = new igraph_t(*fgr,*imap);
  hgraph_t *hgr = new hgraph_t(*fgr,*hmap);
 
  return new SQGraph(fgr,igr,hgr,hnode,1.,hmap,imap);
}

SQGraph::~SQGraph()
{
  delete hmap;
  delete imap;
}

///////////////////////////////////////////////////////////////////////////////
//
// Hierarchical fully-connected  graph

// HFCGraph* HFCGraph::create( Lx,int Ly)
// {
//   // Build the graph by going down the hierarchy.  The arcs added at
//   // this point belong only to the hierarchical graph.  Only the
//   // lowest-level (leaf) nodes belong to the individual graph.  Links
//   // among individuals will be added later

//   fgraph_t *fgr = new fgraph_t;
//   fgraph_t::ArcMap<bool> *hmap =   // hierarchical graph has same nodes, less arcs
//     new fgraph_t::ArcMap<bool>(*fgr,false);
//   fgraph_t::NodeMap<bool> *imap =   // individual graph has only the individual nodes
//     new fgraph_t::NodeMap<bool>(*fgr,true);
//   fgr=build_tree(levels,hmap,imap);

//   // Now add the links among individuals
  
//   for (fgraph_t::NodeIt inode(*fgr); inode!=lemon::INVALID; ++inode)gi {
//     if (inode==hnode) continue;
//     auto arc=fgr->addArc(hnode,inode);
//     (*hmap)[arc]=true;
//   }
//   igraph_t *igr = new igraph_t(*fgr,*imap);
//   hgraph_t *hgr = new hgraph_t(*fgr,*hmap);
 
//   return new SQGraph(fgr,igr,hgr,hnode,1.,hmap,imap);
// }


// // Recursively build a tree starting at given level
// node_t HFCGraph::build_tree(int level)
// {
//   node_t subtree=fgr->addNode();
//   (*imap)[subtree]= level==0;

//   if (level>0) {
//     int M=noffspring(level);
//     for (int i=0; i<M; ++i) {
//       node_t n;
//       n=build_tree(level-1);
//       auto arc=fgr->addArc(subtree,n);
//       (*hmap)[arc]=true;
//     }
//   }

//   return subtree;
// }
