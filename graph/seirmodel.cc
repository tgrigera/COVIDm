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
