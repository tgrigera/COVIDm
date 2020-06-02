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

