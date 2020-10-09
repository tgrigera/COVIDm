/*
 * seirmodel.cc -- classes for the SEIR model on graphs
 *
 * Here we provide some template specializaitons for the seirmodel class.
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

#include "seirmodel.hh"

template<>
void SEEIIR_model<MWFCGraph>::compute_rates(typename MWFCGraph::igraph_t::Node node)
{
  auto &noded=inodemap[node];
  double w,rate;

  switch(noded.state) {
  case SEEIIR_node::S:
    w=0;
    for (typename MWFCGraph::igraph_t::NodeIt snode(egraph.igraph); snode!=lemon::INVALID; ++snode) {
      if (snode==node) continue;
      auto tnode=inodemap[snode];
      if (tnode.state==SEEIIR_node::I1 || tnode.state==SEEIIR_node::I2) {
	w+=egraph.arc_weight(node,snode);
      }
    }
    rate=beta*w;
    break;
  case SEEIIR_node::E1:
    rate=sigma1;
    break;
  case SEEIIR_node::E2:
    rate=sigma2;
    break;
  case SEEIIR_node::I1:
    rate=gamma1;
    break;
  case SEEIIR_node::I2:
    rate=gamma2;
    break;
  case SEEIIR_node::R:
    rate=0;
    break;
  }

  transitions[noded.itransition].rate=rate;
}

template<>
void SEEIIR_model<MWFCGraph>::apply_transition(int itran)
{
  static int n_rate_updates=0;

  auto node=egraph.inode(transitions[itran].nodeid);
  auto &noded=inodemap[node];
  enum {none, new_infection, new_recovery} need_recomp=none;
  
  switch(noded.state) {
  case SEEIIR_node::S:
    noded.state=SEEIIR_node::E1;
    egraph.for_each_anode(node,
			  [this](MWFCGraph::hnode_t hnode) ->void
			  {aggregate_data* anode=this->anodemap[hnode]; anode->NS--; anode->NE1++; anode->Eacc++;} );
    break;

  case SEEIIR_node::E1:
    noded.state=SEEIIR_node::E2;
    egraph.for_each_anode(node,
			  [this](MWFCGraph::hnode_t hnode) ->void
			  {aggregate_data* anode=this->anodemap[hnode]; anode->NE1--; anode->NE2++;} );
    break;

  case SEEIIR_node::E2:
    noded.state=SEEIIR_node::I1;
    egraph.for_each_anode(node,
			  [this](MWFCGraph::hnode_t hnode) ->void
			  {aggregate_data* anode=this->anodemap[hnode]; anode->inf_accum++; anode->inf_close++; anode->NE2--; anode->NI1++;} );
    need_recomp=new_infection;
    break;

  case SEEIIR_node::I1:
    noded.state=SEEIIR_node::I2;
    egraph.for_each_anode(node,
			  [this](MWFCGraph::hnode_t hnode) ->void
			  {aggregate_data* anode=this->anodemap[hnode]; anode->NI1--; anode->NI2++;} );
    break;

  case SEEIIR_node::I2:
    noded.state=SEEIIR_node::R;
    egraph.for_each_anode(node,
			  [this](MWFCGraph::hnode_t hnode) ->void
			  {aggregate_data* anode=this->anodemap[hnode]; anode->NI2--; anode->NR++;} );
    need_recomp=new_recovery;
    break;
    
  case SEEIIR_node::R:
    std::cerr << "Internal error: should not be R\n";
    I_AM_HERE;
    exit(10);
    break;
  }

  compute_rates(node);
  if (need_recomp==none) return;
  n_rate_updates++;
  if (n_rate_updates>egraph.inode_count/10) {
    for (typename MWFCGraph::igraph_t::InArcIt arc(egraph.igraph,node); arc!=lemon::INVALID; ++arc)
      compute_rates(egraph.igraph.source(arc));
    n_rate_updates=0;
    return;
  }

  double rsign=1;
  if (need_recomp==new_recovery) rsign=-1;

  for (MWFCGraph::igraph_t::InArcIt arc(egraph.igraph,node); arc!=lemon::INVALID; ++arc) {
    auto node2=egraph.igraph.source(arc);
    auto &node2d=inodemap[node2];
    if (node2d.state==SEEIIR_node::S) {
      transitions[node2d.itransition].rate += rsign*beta*egraph.arc_weight(arc);
      if (transitions[node2d.itransition].rate<0)            // The correction could bring the rate to less than 0 due to numerical error
	transitions[node2d.itransition].rate=0;              // we correct this because negative rates cause errors in the binary search
    }
  }

}

// This had to be specialized to avoid using arcs with the fully connected graph
// because lemon cannot correctly handle arcs in graphs larger than 47000 nodes or so

template<>
void SEEIIR_model<MWFCGraph>::add_imported(Forced_transition* ii)
{
  typename MWFCGraph::igraph_t::Node node;

  if (ii->new_infected > anodemap[hroot]->NS)
    throw std::runtime_error("Too many imported infections");
  for (int i=0; i<ii->new_infected; ++i) {
    do node=egraph.random_inode(); while(inodemap[node].state!=SEEIIR_node::S);
    inodemap[node].state=SEEIIR_node::I1;
    egraph.for_each_anode(node,
		   [this](typename MWFCGraph::hnode_t hnode)
		   {aggregate_data* anode=this->anodemap[hnode];
		     anode->NS--; anode->NI1++; anode->inf_imported++; anode->inf_accum++; }  );
    compute_rates(node);
    for (typename MWFCGraph::igraph_t::NodeIt snode(egraph.igraph); snode!=lemon::INVALID; ++snode) {
      if (snode==node) continue;
      compute_rates(snode);
    }
  }

  if (ii->new_recovered > anodemap[hroot]->NS)
    throw std::runtime_error("Too many imported infections");
  for (int i=0; i<ii->new_recovered; ++i) {
    do node=egraph.random_inode(); while(inodemap[node].state!=SEEIIR_node::S);
    inodemap[node].state=SEEIIR_node::R;
    egraph.for_each_anode(node,
		   [this](typename MWFCGraph::hnode_t hnode)
		   {aggregate_data* anode=this->anodemap[hnode];
		     anode->NS--; anode->NR++; }  );
    compute_rates(node);
  }

}
