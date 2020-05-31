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

#include "emodel.hh"
#include "egraph.hh"

///////////////////////////////////////////////////////////////////////////////
//
// SIR model

template<typename IGraph>
class SIR_model : public Epidemiological_model_graph_base<IGraph> {
public:
  SIR_model(IGraph&);
  void set_all_susceptible();
  void apply_transition(int);
  void compute_rates(typename IGraph::graph_t::Node);
  void add_imported_infections(Imported_infection*);
  void set_beta(double b) {beta=b;}
  void set_gamma(double g) {gamma=g;}

  struct SIR_node {
    enum {S,I,R} state;
    int  itransition;
  } ;

  int NS,NI,NR;

private:
  double beta,gamma;

  typename IGraph::graph_t::template NodeMap<SIR_node> inodemap;
  using Epidemiological_model_graph_base<IGraph>::igraph;
  using Epidemiological_model_graph_base<IGraph>::transitions;
} ;

template<typename IGraph>
SIR_model<IGraph>::SIR_model(IGraph& igraph) :
  Epidemiological_model_graph_base<IGraph>(igraph),
  inodemap(igraph.graph),
  NS(igraph.node_count),
  NI(0),
  NR(0)
{
  transitions.clear();
  for (typename IGraph::graph_t::NodeIt inode(igraph.graph); inode!=lemon::INVALID; ++inode) {
    inodemap[inode].state=SIR_node::S;
    inodemap[inode].itransition=transitions.size();
    Epidemiological_model::transition tr={igraph.id(inode),0,0};
    transitions.push_back(tr);
  }
}

template<typename IGraph>
void SIR_model<IGraph>::set_all_susceptible()
{
  for (typename IGraph::graph_t::NodeIt inode(igraph.graph); inode!=lemon::INVALID; ++inode)
    inodemap[inode].state=SIR_node::S;
  NS=igraph.node_count;
  NI=NR=0;
}

template<typename IGraph>
void SIR_model<IGraph>::compute_rates(typename IGraph::graph_t::Node node)
{
  auto &noded=inodemap[node];
  double w,rate;

  switch(noded.state) {
  case SIR_node::S:
    w=0;
    for (typename IGraph::graph_t::OutArcIt arc(igraph.graph,node); arc!=lemon::INVALID; ++arc) {
      auto tnode=inodemap[igraph.graph.target(arc)];
      if (tnode.state==SIR_node::I) w+=igraph.arc_weight(arc);
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

template<typename IGraph>
void SIR_model<IGraph>::apply_transition(int itran)
{
  auto node=igraph.node(transitions[itran].nodeid);
  auto &noded=inodemap[node];
  switch(noded.state) {
  case SIR_node::S:
    noded.state=SIR_node::I;
    NS--;
    NI++;
    break;

  case SIR_node::I:
    noded.state=SIR_node::R;
    NI--;
    NR++;
    break;

  case SIR_node::R:
    std::cerr << "Internal error: should not be R\n";
    I_AM_HERE;
    exit(10);
    break;
  }

  compute_rates(node);
  for (typename IGraph::graph_t::InArcIt arc(igraph.graph,node); arc!=lemon::INVALID; ++arc)
    compute_rates(igraph.graph.source(arc));
}

template<typename IGraph>
void SIR_model<IGraph>::add_imported_infections(Imported_infection* ii)
{
  typename IGraph::graph_t::Node node;

  for (int i=0; i<ii->new_cases; ++i) {
    do node=igraph.random_node(); while(inodemap[node].state!=SIR_node::S);
    inodemap[node].state=SIR_node::I;
    NS--;
    NI++;
    compute_rates(node);
    for (typename IGraph::graph_t::InArcIt arc(igraph.graph,node); arc!=lemon::INVALID; ++arc)
      compute_rates(igraph.graph.source(arc));
  }
}

#endif /* SIRMODEL_HH */
