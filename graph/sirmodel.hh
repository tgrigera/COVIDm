/*
 * sirmodel.hh -- classes for the SIR model on graphs
 *
 * The classes here implement different epidemiological models.  They
 * are coded generically and can be used over any graph provided by
 * egraph.hh.
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

#ifndef SIRMODEL_HH
#define SIRMODEL_HH

#include <assert.h>

#include "emodel.hh"
#include "egraph.hh"

///////////////////////////////////////////////////////////////////////////////
//
// SIR model

template<typename EGraph>
class SIR_model : public Epidemiological_model_graph_base<EGraph> {
public:
  SIR_model(EGraph&);
  void set_all_susceptible();
  void apply_transition(int);
  void compute_rates(typename EGraph::igraph_t::Node);
  void add_imported_infections(Imported_infection*);
  void set_beta(double b) {beta=b;}
  void set_gamma(double g) {gamma=g;}

  struct SIR_node {
    enum {S,I,R} state;
    int  itransition;
  } ;

  struct aggregate_node {
    int  NS,NI,NR;
  } ;

  typename EGraph::node_t   hroot;
  typename EGraph::hgraph_t::template NodeMap<aggregate_node*> anodemap;
  

private:
  double beta,gamma;

  typename EGraph::igraph_t::template NodeMap<SIR_node> inodemap;
  using Epidemiological_model_graph_base<EGraph>::egraph;
  using Epidemiological_model_graph_base<EGraph>::transitions;

  void recompute_counts(typename EGraph::node_t);
} ;

template<typename EGraph>
SIR_model<EGraph>::SIR_model(EGraph& egraph) :
  Epidemiological_model_graph_base<EGraph>(egraph),
  hroot(egraph.hroot),
  anodemap(egraph.hgraph,0),
  inodemap(egraph.igraph)
{
  transitions.clear();
  for (typename EGraph::igraph_t::NodeIt inode(egraph.igraph); inode!=lemon::INVALID; ++inode) {
    inodemap[inode].state=SIR_node::S;
    inodemap[inode].itransition=transitions.size();
    Epidemiological_model::transition tr={egraph.id(inode),0,0};
    transitions.push_back(tr);
  }
}

template<typename EGraph>
void SIR_model<EGraph>::set_all_susceptible()
{
  for (typename EGraph::igraph_t::NodeIt inode(egraph.igraph); inode!=lemon::INVALID; ++inode) {
    inodemap[inode].state=SIR_node::S;
  }
  recompute_counts(egraph.hroot);

  assert(anodemap[hroot]->NS==egraph.inode_count);
}

template<typename EGraph>
void SIR_model<EGraph>::recompute_counts(typename EGraph::node_t lroot)
{
  aggregate_node* anode=anodemap[lroot];
  if (anode==0) {
    anode=new aggregate_node;
    anodemap[lroot]=anode;
  }
  anode->NS=anode->NI=anode->NR=0;
  for (typename EGraph::hgraph_t::OutArcIt arc(egraph.hgraph,lroot); arc!=lemon::INVALID; ++arc) {
    auto lnode=egraph.hgraph.target(arc);
    if (countOutArcs(egraph.hgraph,lnode)>0) {
      recompute_counts(lnode);
      anode->NS+=anodemap[lnode]->NS;
      anode->NI+=anodemap[lnode]->NI;
      anode->NR+=anodemap[lnode]->NR;
    } else {
      switch(inodemap[lnode].state) {
      case SIR_node::S:
	anode->NS++;
	break;
      case SIR_node::I:
	anode->NI++;
	break;
      case SIR_node::R:
	anode->NR++;
	break;
      }
    }
  }
}

template<typename EGraph>
void SIR_model<EGraph>::compute_rates(typename EGraph::igraph_t::Node node)
{
  auto &noded=inodemap[node];
  double w,rate;

  switch(noded.state) {
  case SIR_node::S:
    w=0;
    for (typename EGraph::igraph_t::OutArcIt arc(egraph.igraph,node); arc!=lemon::INVALID; ++arc) {
      auto tnode=inodemap[egraph.igraph.target(arc)];
      if (tnode.state==SIR_node::I) w+=egraph.arc_weight(arc);
    }
    rate=beta*w;
    break;
  case SIR_node::I:
    rate=gamma;
    break;
  case SIR_node::R:
    rate=0;
    break;
  }

  transitions[noded.itransition].rate=rate;
}

template<typename EGraph>
void SIR_model<EGraph>::apply_transition(int itran)
{
  auto node=egraph.inode(transitions[itran].nodeid);
  auto &noded=inodemap[node];
  switch(noded.state) {
  case SIR_node::S:
    noded.state=SIR_node::I;
    egraph.for_each_anode(node,
			  [this](typename EGraph::node_t hnode) ->void
			  {aggregate_node* anode=this->anodemap[hnode]; anode->NS--; anode->NI++;}    );
    break;

  case SIR_node::I:
    noded.state=SIR_node::R;
    egraph.for_each_anode(node,
		   [this](typename EGraph::node_t hnode)
		   {aggregate_node* anode=this->anodemap[hnode]; anode->NI--; anode->NR++;}    );
    break;

  case SIR_node::R:
    std::cerr << "Internal error: should not be R\n";
    I_AM_HERE;
    exit(10);
    break;
  }

  compute_rates(node);
  for (typename EGraph::igraph_t::InArcIt arc(egraph.igraph,node); arc!=lemon::INVALID; ++arc)
    compute_rates(egraph.igraph.source(arc));
}

template<typename EGraph>
void SIR_model<EGraph>::add_imported_infections(Imported_infection* ii)
{
  typename EGraph::igraph_t::Node node;

  for (int i=0; i<ii->new_cases; ++i) {
    do node=egraph.random_inode(); while(inodemap[node].state!=SIR_node::S);
    inodemap[node].state=SIR_node::I;
    egraph.for_each_anode(node,
		   [this](typename EGraph::node_t hnode)
		   {aggregate_node* anode=this->anodemap[hnode]; anode->NS--; anode->NI++;}    );
    compute_rates(node);
    for (typename EGraph::igraph_t::InArcIt arc(egraph.igraph,node); arc!=lemon::INVALID; ++arc)
      compute_rates(egraph.igraph.source(arc));
  }
}

#endif /* SIRMODEL_HH */
