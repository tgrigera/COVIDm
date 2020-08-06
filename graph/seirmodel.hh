/*
 * seirmodel.hh -- classes for the SEIR model on graphs
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

#ifndef SEIRMODEL_HH
#define SEIRMODEL_HH

///////////////////////////////////////////////////////////////////////////////
//
// SEEIIR model

template<typename EGraph>
class SEEIIR_model : public Epidemiological_model_graph_base<EGraph> {
public:
  SEEIIR_model(EGraph&);
  ~SEEIIR_model();
  void set_all_susceptible();
  void apply_transition(int);
  void compute_rates(typename EGraph::igraph_t::Node);
  void add_imported_infections(Imported_infection*);
  void set_rate_constants(double beta,double sigma1,double sigma2,double gamma1,
			  double gamma2);

  double tinf() { return 1./gamma1 + 1./gamma2;}

  struct SEEIIR_node {
    enum {S,E1,E2,I1,I2,R} state;
    int  itransition;
  } ;
  struct aggregate_data {
    int Ntot;
    int NS,NE1,NE2,NI1,NI2,NR;
    int inf_accum;
    int inf_imported,inf_close,inf_community;
  } ;

  typename EGraph::node_t   hroot;
  typename EGraph::hgraph_t::template NodeMap<aggregate_data*> anodemap;

private:
  double beta,sigma1,sigma2,gamma1,gamma2;

  typename EGraph::igraph_t::template NodeMap<SEEIIR_node> inodemap;
  using Epidemiological_model_graph_base<EGraph>::egraph;
  using Epidemiological_model_graph_base<EGraph>::transitions;

  void recompute_counts(typename EGraph::node_t);
} ;

template<typename EGraph>
SEEIIR_model<EGraph>::SEEIIR_model(EGraph& egraph) :
  Epidemiological_model_graph_base<EGraph>(egraph),
  hroot(egraph.hroot),
  anodemap(egraph.hgraph,0),
  inodemap(egraph.igraph)
{
  transitions.clear();
  for (typename EGraph::igraph_t::NodeIt inode(egraph.igraph); inode!=lemon::INVALID; ++inode) {
    inodemap[inode].state=SEEIIR_node::S;
    inodemap[inode].itransition=transitions.size();
    Epidemiological_model::transition tr={egraph.id(inode),0,0};
    transitions.push_back(tr);
  }
}

template<typename EGraph>
SEEIIR_model<EGraph>::~SEEIIR_model()
{
  for (typename EGraph::hgraph_t::NodeIt anode(egraph.hgraph); anode!=lemon::INVALID; ++anode)
    delete anodemap[anode];
}

template<typename EGraph>
void SEEIIR_model<EGraph>::set_all_susceptible()
{
  for (typename EGraph::igraph_t::NodeIt inode(egraph.igraph); inode!=lemon::INVALID; ++inode) {
    inodemap[inode].state=SEEIIR_node::S;
  }
  recompute_counts(egraph.hroot);

  assert(anodemap[hroot]->NS==egraph.inode_count);
}

template<typename EGraph>
inline void SEEIIR_model<EGraph>::set_rate_constants(double beta_,double sigma1_,
					      double sigma2_,double gamma1_,double gamma2_)
{
  beta=beta_;
  sigma1=sigma1_;
  sigma2=sigma2_;
  gamma1=gamma1_;
  gamma2=gamma2_;
}

template<typename EGraph>
void SEEIIR_model<EGraph>::recompute_counts(typename EGraph::node_t lroot)
{
  aggregate_data* anode=anodemap[lroot];
  if (anode==0) {
    anode=new aggregate_data;
    anodemap[lroot]=anode;
  }
  anode->NS=anode->NE1=anode->NE2=anode->NI1=anode->NI2=anode->NR=0;
  anode->inf_accum=anode->inf_imported=anode->inf_close=anode->inf_community=0;
  anode->Ntot=0;
  for (typename EGraph::hgraph_t::OutArcIt arc(egraph.hgraph,lroot); arc!=lemon::INVALID; ++arc) {
    auto lnode=egraph.hgraph.target(arc);
    if (countOutArcs(egraph.hgraph,lnode)>0) {
      recompute_counts(lnode);
      anode->Ntot+=anodemap[lnode]->Ntot;
      anode->NS+=anodemap[lnode]->NS;
      anode->NE1+=anodemap[lnode]->NE1;
      anode->NE2+=anodemap[lnode]->NE2;
      anode->NI1+=anodemap[lnode]->NI1;
      anode->NI2+=anodemap[lnode]->NI2;
      anode->NR+=anodemap[lnode]->NR;
    } else {
      anode->Ntot++;
      switch(inodemap[lnode].state) {
      case SEEIIR_node::S:
	anode->NS++;
	break;
      case SEEIIR_node::E1:
	anode->NE1++;
	break;
      case SEEIIR_node::E2:
	anode->NE2++;
	break;
      case SEEIIR_node::I1:
	anode->NI1++;
	break;
      case SEEIIR_node::I2:
	anode->NI2++;
	break;
      case SEEIIR_node::R:
	anode->NR++;
	break;
      }
    }
  }
}

template<typename EGraph>
void SEEIIR_model<EGraph>::compute_rates(typename EGraph::igraph_t::Node node)
{
  auto &noded=inodemap[node];
  double w,rate;

  switch(noded.state) {
  case SEEIIR_node::S:
    w=0;
    for (typename EGraph::igraph_t::OutArcIt arc(egraph.igraph,node); arc!=lemon::INVALID; ++arc) {
      auto tnode=inodemap[egraph.igraph.target(arc)];
      if (tnode.state==SEEIIR_node::I1 || tnode.state==SEEIIR_node::I2) {
	w+=egraph.arc_weight(arc);
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

template<typename EGraph>
void SEEIIR_model<EGraph>::apply_transition(int itran)
{
  auto node=egraph.inode(transitions[itran].nodeid);
  auto &noded=inodemap[node];
  bool need_recomp=false;
  
  switch(noded.state) {
  case SEEIIR_node::S:
    noded.state=SEEIIR_node::E1;
    egraph.for_each_anode(node,
			  [this](typename EGraph::node_t hnode) ->void
			  {aggregate_data* anode=this->anodemap[hnode]; anode->NS--; anode->NE1++;} );
    break;

  case SEEIIR_node::E1:
    noded.state=SEEIIR_node::E2;
    egraph.for_each_anode(node,
			  [this](typename EGraph::node_t hnode) ->void
			  {aggregate_data* anode=this->anodemap[hnode]; anode->NE1--; anode->NE2++;} );
    break;

  case SEEIIR_node::E2:
    noded.state=SEEIIR_node::I1;
    egraph.for_each_anode(node,
			  [this](typename EGraph::node_t hnode) ->void
			  {aggregate_data* anode=this->anodemap[hnode]; anode->inf_accum++; anode->inf_close++; anode->NE2--; anode->NI1++;} );
    need_recomp=true;
    break;

  case SEEIIR_node::I1:
    noded.state=SEEIIR_node::I2;
    egraph.for_each_anode(node,
			  [this](typename EGraph::node_t hnode) ->void
			  {aggregate_data* anode=this->anodemap[hnode]; anode->NI1--; anode->NI2++;} );
    break;

  case SEEIIR_node::I2:
    noded.state=SEEIIR_node::R;
    egraph.for_each_anode(node,
			  [this](typename EGraph::node_t hnode) ->void
			  {aggregate_data* anode=this->anodemap[hnode]; anode->NI2--; anode->NR++;} );
    need_recomp=true;
    break;
    
  case SEEIIR_node::R:
    std::cerr << "Internal error: should not be R\n";
    I_AM_HERE;
    exit(10);
    break;
  }

  compute_rates(node);
  if (need_recomp)
    for (typename EGraph::igraph_t::InArcIt arc(egraph.igraph,node); arc!=lemon::INVALID; ++arc)
      compute_rates(egraph.igraph.source(arc));
}

template<typename EGraph>
void SEEIIR_model<EGraph>::add_imported_infections(Imported_infection* ii)
{
  typename EGraph::igraph_t::Node node;

  if (ii->new_cases > anodemap[hroot]->NS)
    throw std::runtime_error("Too many imported infections");
  for (int i=0; i<ii->new_cases; ++i) {
    do node=egraph.random_inode(); while(inodemap[node].state!=SEEIIR_node::S);
    inodemap[node].state=SEEIIR_node::I1;
    egraph.for_each_anode(node,
		   [this](typename EGraph::node_t hnode)
		   {aggregate_data* anode=this->anodemap[hnode];
		     anode->NS--; anode->NI1++; anode->inf_imported++; anode->inf_accum++; }  );
    compute_rates(node);
    for (typename EGraph::igraph_t::InArcIt arc(egraph.igraph,node); arc!=lemon::INVALID; ++arc)
      compute_rates(egraph.igraph.source(arc));
  }


}


///////////////////////////////////////////////////////////////////////////////
//
// event: change rate constants

template<typename EGraph>
class Rate_constant_change : public Event {
public:
  Rate_constant_change(double time,double beta,double sigma1,double sigma2,
		       double gamma1,double gamma2) :
    Event(time), beta(beta), sigma1(sigma1), sigma2(sigma2),
    gamma1(gamma1), gamma2(gamma2) {}
  virtual void apply(Epidemiological_model *);

  double beta,sigma1,sigma2,gamma1,gamma2;
} ;

template<typename EGraph>
inline void Rate_constant_change<EGraph>::apply(Epidemiological_model *em)
{
  (dynamic_cast<SEEIIR_model<EGraph>*>(em))->set_rate_constants(beta,sigma1,sigma2,gamma1,gamma2);
  (dynamic_cast<SEEIIR_model<EGraph>*>(em))->compute_all_rates();
}

#endif /* SEIRMODEL_HH */
