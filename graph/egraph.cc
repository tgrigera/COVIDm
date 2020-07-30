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

#include "egraph.hh"

///////////////////////////////////////////////////////////////////////////////
//
// Fully-connected graph

FCGraph* FCGraph::create(int N)
{
  // The graph is built copying the structure of a lemon FullGraph
  lemon::FullGraph fcg(N);
  fgraph_t *fgr = new fgraph_t;
  lemon::digraphCopy(fcg,*fgr).run();

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

  return new FCGraph(fgr,igr,hgr,hnode,1./N,hmap,imap);
}

//
// Random multiplicative weights for the FC grah
//
// Assigns symmetrical weights to each arc according to J_ij = beta_i * beta_j / A
// betas are taken from user-provided random distribution betadist() and
// A= inode_count * beta_scale
//
void FCGraph::set_weights_random_multiplicative(double (*betadist)(),double beta_scale)
{
  // obtain random betas
  igraph_t::NodeMap<double> beta(igraph);
  for (igraph_t::NodeIt node(igraph); node!=lemon::INVALID; ++node) {
    beta[node]=betadist();
  }

  // loop over arcs and set weight
  double aN=beta_scale*inode_count;
  for (igraph_t::ArcIt arc(igraph); arc!=lemon::INVALID; ++arc) {
    inode_t i=igraph.source(arc);
    inode_t j=igraph.target(arc);
    arcmap[arc] = beta[i]*beta[j]/aN;
  }
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
